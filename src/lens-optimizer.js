import { getUpdateStrategies } from "./update-strategies.js";
import { polar } from "./util.js";

/**
 * Each lens pixel has:
 *   - baseIndex (lens.x)
 *   - dispersion (lens.y)
 *   - chi2       (lens.z)
 *   - chi3       (lens.w)
 *
 * This class aggregates wave data from fundamental/shg fields,
 * updates a radial-sliced lens representation, and returns
 * a 4-channel Float32Array for the GPU.
 */
export class LensOptimizer {
  constructor(config) {
    this.config = config;

    // lensModes[z][s] => { baseIndex, dispersion, chi2, chi3 }
    this.lensModes = Array(config.fresnelZones)
      .fill()
      .map(() =>
        Array(config.numSectors)
          .fill()
          .map(() => ({
            baseIndex: 1.0, // e.g. n=1.0 "vacuum" or baseline
            dispersion: 0.0, // no dispersion by default
            chi2: 0.0, // no chi2 by default
            chi3: 0.0, // no chi3 by default
          })),
      );

    // We track "metrics" from fundamental and SHG
    this.fundMetrics = Array(config.fresnelZones)
      .fill()
      .map(() =>
        Array(config.numSectors)
          .fill()
          .map(() => ({
            phaseGradMag: 0,
            kerrStrength: 0,
          })),
      );

    this.shgMetrics = Array(config.fresnelZones)
      .fill()
      .map(() =>
        Array(config.numSectors)
          .fill()
          .map(() => ({
            phaseMatch: 0,
            chiRatio: 0,
          })),
      );

    // Initialize ADAM optimizer state for 4 parameters each:
    this.adamState = adam.createState([config.fresnelZones, config.numSectors]);

    this.updateStrategies = getUpdateStrategies(config);

    // Store a list of progress metrics over time
    this.progressHistory = [];
  }

  /**
   * updateLens(...) accumulates fundamental and SHG data in
   * radial slices, then calls a chosen updateStrategy to adjust
   * lensModes. We use ADAM to do gradient-based updates on
   * (baseIndex, dispersion, chi2, chi3).
   */
  updateLens(fundamentalData, shgData, width, height) {
    // 1) Accumulate complex sums per zone/sector
    const fundSum = Array(this.config.fresnelZones)
      .fill()
      .map(() =>
        Array(this.config.numSectors)
          .fill()
          .map(() => ({ real: 0, imag: 0 })),
      );
    const shgSum = Array(this.config.fresnelZones)
      .fill()
      .map(() =>
        Array(this.config.numSectors)
          .fill()
          .map(() => ({ real: 0, imag: 0 })),
      );
    const modeCount = Array(this.config.fresnelZones)
      .fill()
      .map(() => Array(this.config.numSectors).fill(0));

    // Clear metric accumulators (kept for structure, not used for .zw anymore)
    this.fundMetrics.forEach((zone) =>
      zone.forEach((sector) => {
        sector.phaseGradMag = 0;
        sector.kerrStrength = 0;
      }),
    );
    this.shgMetrics.forEach((zone) =>
      zone.forEach((sector) => {
        sector.phaseMatch = 0;
        sector.chiRatio = 0;
      }),
    );

    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    // (A) Bucket complex sums
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const dx = x - centerX;
        const dy = y - centerY;
        const r = Math.hypot(dx, dy);
        if (r <= this.config.lensRadius) {
          const coords = polar.getZoneAndSector(
            x,
            y,
            centerX,
            centerY,
            this.config.lensRadius,
            this.config.fresnelZones,
            this.config.numSectors,
          );
          if (!coords) continue;

          const idx = (y * width + x) * 4;

          // Only take real/imag; ignore .zw (ADE polarization in the sim)
          const fRe = fundamentalData[idx + 0];
          const fIm = fundamentalData[idx + 1];
          const sRe = shgData[idx + 0];
          const sIm = shgData[idx + 1];

          fundSum[coords.zone][coords.sector].real += fRe;
          fundSum[coords.zone][coords.sector].imag += fIm;
          shgSum[coords.zone][coords.sector].real += sRe;
          shgSum[coords.zone][coords.sector].imag += sIm;

          modeCount[coords.zone][coords.sector]++;
        }
      }
    }

    // (B) Strategy + ADAM
    this.adamState.t += 1;
    const { bc1, bc2 } = adam.getBiasCorrection(this.adamState);

    const chosenStrategy = this.config.updateStrategy;
    const updateStrategy = this.updateStrategies[chosenStrategy];
    if (!updateStrategy) {
      alert("Couldn't find update strategy: " + chosenStrategy);
      return this.lensModes;
    }

    let totalLoss = 0;
    let countUpdates = 0;

    for (let z = 0; z < this.config.fresnelZones; z++) {
      for (let s = 0; s < this.config.numSectors; s++) {
        const count = modeCount[z][s];
        if (count <= 0) continue;

        // Average complex field per bucket
        const fRe = fundSum[z][s].real / count;
        const fIm = fundSum[z][s].imag / count;
        const sRe = shgSum[z][s].real / count;
        const sIm = shgSum[z][s].imag / count;

        // Derived metrics (proxies), now that we don't read .zw
        const fundAmp = Math.hypot(fRe, fIm);
        const shgAmp = Math.hypot(sRe, sIm);
        const fundPhase = Math.atan2(fIm, fRe);
        const shgPhase = Math.atan2(sIm, sRe);

        const fund = {
          real: fRe,
          imag: fIm,
          // Use intensity as a Kerr proxy; leave phaseGradMag at 0 unless you add a gradient buffer
          kerrStrength: fundAmp * fundAmp,
          phaseGradMag: 0.0,
        };

        const shg = {
          real: sRe,
          imag: sIm,
          // Coherence proxy: +1 when shg ≈ e^{i2φ1}
          phaseMatch: Math.cos(shgPhase - 2 * fundPhase),
          // Rough χ ratio proxy: |E2| / |E1|^2
          chiRatio: shgAmp / (fundAmp * fundAmp + 1e-6),
        };

        const lensParams = this.lensModes[z][s];

        const { gradients, loss } = updateStrategy(
          fund,
          shg,
          count,
          z,
          s,
          lensParams,
          this.config,
        );

        // ---- ADAM moments ----
        this.adamState.m_baseIndex[z][s] =
          this.adamState.beta1 * this.adamState.m_baseIndex[z][s] +
          (1 - this.adamState.beta1) * gradients.baseIndex;
        this.adamState.v_baseIndex[z][s] =
          this.adamState.beta2 * this.adamState.v_baseIndex[z][s] +
          (1 - this.adamState.beta2) *
            gradients.baseIndex *
            gradients.baseIndex;

        this.adamState.m_dispersion[z][s] =
          this.adamState.beta1 * this.adamState.m_dispersion[z][s] +
          (1 - this.adamState.beta1) * gradients.dispersion;
        this.adamState.v_dispersion[z][s] =
          this.adamState.beta2 * this.adamState.v_dispersion[z][s] +
          (1 - this.adamState.beta2) *
            gradients.dispersion *
            gradients.dispersion;

        this.adamState.m_chi2[z][s] =
          this.adamState.beta1 * this.adamState.m_chi2[z][s] +
          (1 - this.adamState.beta1) * gradients.chi2;
        this.adamState.v_chi2[z][s] =
          this.adamState.beta2 * this.adamState.v_chi2[z][s] +
          (1 - this.adamState.beta2) * gradients.chi2 * gradients.chi2;

        this.adamState.m_chi3[z][s] =
          this.adamState.beta1 * this.adamState.m_chi3[z][s] +
          (1 - this.adamState.beta1) * gradients.chi3;
        this.adamState.v_chi3[z][s] =
          this.adamState.beta2 * this.adamState.v_chi3[z][s] +
          (1 - this.adamState.beta2) * gradients.chi3 * gradients.chi3;

        // ---- Bias-corrected step ----
        const mBase_hat = this.adamState.m_baseIndex[z][s] * bc1;
        const vBase_hat = this.adamState.v_baseIndex[z][s] * bc2;

        const mDisp_hat = this.adamState.m_dispersion[z][s] * bc1;
        const vDisp_hat = this.adamState.v_dispersion[z][s] * bc2;

        const mChi2_hat = this.adamState.m_chi2[z][s] * bc1;
        const vChi2_hat = this.adamState.v_chi2[z][s] * bc2;

        const mChi3_hat = this.adamState.m_chi3[z][s] * bc1;
        const vChi3_hat = this.adamState.v_chi3[z][s] * bc2;

        const lr = this.config.learningRate / (1 + 0.01 * this.adamState.t);

        const lp = this.lensModes[z][s];
        lp.baseIndex +=
          (lr * mBase_hat) / (Math.sqrt(vBase_hat) + this.adamState.epsilon);
        lp.dispersion +=
          (lr * mDisp_hat) / (Math.sqrt(vDisp_hat) + this.adamState.epsilon);
        lp.chi2 +=
          (lr * mChi2_hat) / (Math.sqrt(vChi2_hat) + this.adamState.epsilon);
        lp.chi3 +=
          (lr * mChi3_hat) / (Math.sqrt(vChi3_hat) + this.adamState.epsilon);

        totalLoss += loss;
        countUpdates++;
      }
    }

    // (C) Progress metric
    if (countUpdates > 0) {
      this.progressHistory.push({
        iteration: this.adamState.t,
        metric: totalLoss / countUpdates,
      });
    }
  }

  /**
   * getLensData(...) returns a Float32Array (width x height x 4)
   * with the final lens parameters for each pixel.
   *   channel 0 => baseIndex
   *   channel 1 => dispersion
   *   channel 2 => chi2
   *   channel 3 => chi3
   */
  getLensData(width, height) {
    const data = new Float32Array(width * height * 4);
    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const idx = (y * width + x) * 4;
        const coords = polar.getZoneAndSector(
          x,
          y,
          centerX,
          centerY,
          this.config.lensRadius,
          this.config.fresnelZones,
          this.config.numSectors,
        );

        if (coords) {
          const mode = this.lensModes[coords.zone][coords.sector];

          // Map your learned parameters to the shader's packing:
          const n_o = mode.baseIndex; // ordinary index
          const n_e = mode.baseIndex; // isotropic fallback (n_e == n_o)
          const axisAngle = 0.0; // no optic-axis model (radians)
          const dispersion = mode.dispersion;

          data[idx + 0] = n_o; // R: n_o
          data[idx + 1] = n_e; // G: n_e
          data[idx + 2] = axisAngle; // B: optic axis angle
          data[idx + 3] = dispersion; // A: dispersion coefficient
        } else {
          // Outside the lens: vacuum-ish, no dispersion
          data[idx + 0] = 1.0; // n_o
          data[idx + 1] = 1.0; // n_e
          data[idx + 2] = 0.0; // axis angle
          data[idx + 3] = 0.0; // dispersion
        }
      }
    }

    return data;
  }

  setBaseRMask(baseRMask) {
    // Expect a Float32Array of length width*height from createGainMaskData(...)
    this.baseRMask = baseRMask;
  }

  getGainMaskData(width, height) {
    if (!this.baseRMask) {
      throw new Error(
        "LensOptimizer.getGainMaskData: baseRMask not set. Call setBaseRMask(...) after getInitialFieldState.",
      );
    }

    // ---- Safety knobs (good defaults) ----
    const chi2Scale = this.config.chi2Scale ?? 0.02; // keep weights small
    const chi3Scale = this.config.chi3Scale ?? 0.02;
    const MAX_W = this.config.maxChiWeight ?? 0.25; // hard cap (post-scale)
    const EMA_ALPHA = this.config.gainMaskEma ?? 0.2; // smoothing (0..1)
    const SLEW = this.config.gainMaskSlew ?? 0.02; // max step per update

    const softplus = (x) =>
      Math.log1p(Math.exp(Math.max(-40, Math.min(40, x))));

    const data = new Float32Array(width * height * 4);
    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    // Allocate previous GB (for smoothing) on first run
    if (!this._prevG || this._prevG.length !== width * height) {
      this._prevG = new Float32Array(width * height);
      this._prevB = new Float32Array(width * height);
      // seed with tiny values so first frame isn’t a jump
      for (let i = 0; i < this._prevG.length; i++) {
        this._prevG[i] = 0.0001;
        this._prevB[i] = 0.0001;
      }
    }

    const clamp01 = (v) => Math.max(0, Math.min(1, v));
    const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));

    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const idx4 = (y * width + x) * 4;
        const idx1 = y * width + x;

        const coords = polar.getZoneAndSector(
          x,
          y,
          centerX,
          centerY,
          this.config.lensRadius,
          this.config.fresnelZones,
          this.config.numSectors,
        );

        const R = this.baseRMask[idx1] || 0.0;
        let G = 0.0,
          B = 0.0;

        if (coords) {
          const mode = this.lensModes[coords.zone][coords.sector];

          // raw (nonnegative) weights from optimizer params
          let rawG = chi2Scale * softplus(mode.chi2);
          let rawB = chi3Scale * softplus(mode.chi3);

          // hard cap
          rawG = Math.min(rawG, MAX_W);
          rawB = Math.min(rawB, MAX_W);

          // slew-limit per update to avoid instant spikes
          const prevG = this._prevG[idx1];
          const prevB = this._prevB[idx1];
          const stepG = clamp(rawG - prevG, -SLEW, SLEW);
          const stepB = clamp(rawB - prevB, -SLEW, SLEW);
          const slG = prevG + stepG;
          const slB = prevB + stepB;

          // EMA smoothing
          G = EMA_ALPHA * slG + (1 - EMA_ALPHA) * prevG;
          B = EMA_ALPHA * slB + (1 - EMA_ALPHA) * prevB;

          // final clamp to [0, MAX_W]
          G = clamp(G, 0, MAX_W);
          B = clamp(B, 0, MAX_W);

          // store for next frame
          this._prevG[idx1] = G;
          this._prevB[idx1] = B;
        }

        data[idx4 + 0] = clamp01(R); // R: linear gain mask (unchanged)
        data[idx4 + 1] = G; // G: χ(2) weight (safe)
        data[idx4 + 2] = B; // B: χ(3) weight (safe)
        data[idx4 + 3] = 1.0; // A
      }
    }

    return data;
  }

  getProgress() {
    // Return the array of recorded progress steps
    return this.progressHistory;
  }

  resetProgress() {
    this.progressHistory = [];
  }
}

// =========================================================
// ADAM optimizer for the four lens parameters
// =========================================================
const adam = {
  createState(dimensions) {
    const [Z, S] = dimensions; // #fresnelZones, #numSectors
    // For each parameter, we store m_..., v_...
    //   baseIndex, dispersion, chi2, chi3
    return {
      beta1: 0.9,
      beta2: 0.999,
      epsilon: 1e-8,
      // baseIndex:
      m_baseIndex: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      v_baseIndex: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      // dispersion:
      m_dispersion: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      v_dispersion: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      // chi2:
      m_chi2: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      v_chi2: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      // chi3:
      m_chi3: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),
      v_chi3: Array(Z)
        .fill()
        .map(() => Array(S).fill(0)),

      t: 0,
    };
  },

  getBiasCorrection(state) {
    const bc1 = 1.0 / (1.0 - Math.pow(state.beta1, state.t));
    const bc2 = 1.0 / (1.0 - Math.pow(state.beta2, state.t));
    return { bc1, bc2 };
  },
};
