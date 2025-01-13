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
    // 1) Allocate accumulators for fundamental + SHG wave sums
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

    // 2) Reset the “extra” metrics we store in this.fundMetrics / this.shgMetrics
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

    // --------------------------------------------
    // (A) ACCUMULATE WAVE DATA INTO ZONE/SECTOR BUCKETS
    // --------------------------------------------
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        // Check if (x,y) is inside the lens radius
        const dx = x - centerX;
        const dy = y - centerY;
        const r = Math.sqrt(dx * dx + dy * dy);
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

          if (coords) {
            const idx = (y * width + x) * 4;

            // fundamentalData, shgData assumed to be RGBA or [real, imag, extra1, extra2]
            const fundReal = fundamentalData[idx + 0];
            const fundImag = fundamentalData[idx + 1];
            const shgReal = shgData[idx + 0];
            const shgImag = shgData[idx + 1];

            // Accumulate amplitude
            fundSum[coords.zone][coords.sector].real += fundReal;
            fundSum[coords.zone][coords.sector].imag += fundImag;
            shgSum[coords.zone][coords.sector].real += shgReal;
            shgSum[coords.zone][coords.sector].imag += shgImag;

            // Accumulate “metrics” (phaseGradMag, kerrStrength, phaseMatch, chiRatio)
            this.fundMetrics[coords.zone][coords.sector].phaseGradMag +=
              fundamentalData[idx + 2];
            this.fundMetrics[coords.zone][coords.sector].kerrStrength +=
              fundamentalData[idx + 3];

            this.shgMetrics[coords.zone][coords.sector].phaseMatch +=
              shgData[idx + 2];
            this.shgMetrics[coords.zone][coords.sector].chiRatio +=
              shgData[idx + 3];

            modeCount[coords.zone][coords.sector]++;
          }
        }
      }
    }

    // --------------------------------------------
    // (B) UPDATE LENS MODES VIA CHOSEN “STRATEGY” + ADAM
    // --------------------------------------------
    // increment ADAM iteration
    this.adamState.t += 1;
    const { bc1, bc2 } = adam.getBiasCorrection(this.adamState); // e.g. (1 - β1^t), etc.

    // console.log(bc1, bc2);

    const chosenStrategy = this.config.updateStrategy;
    const updateStrategy = this.updateStrategies[chosenStrategy];
    if (!updateStrategy) {
      alert("Couldn't find update strategy: " + chosenStrategy);
      return this.lensModes; // fallback
    }

    let totalLoss = 0;
    let countUpdates = 0;

    // Traverse each zone/sector
    for (let z = 0; z < this.config.fresnelZones; z++) {
      for (let s = 0; s < this.config.numSectors; s++) {
        const count = modeCount[z][s];
        if (count > 0) {
          // Average wave data
          const fund = {
            real: fundSum[z][s].real / count,
            imag: fundSum[z][s].imag / count,
            phaseGradMag: this.fundMetrics[z][s].phaseGradMag / count,
            kerrStrength: this.fundMetrics[z][s].kerrStrength / count,
          };
          const shg = {
            real: shgSum[z][s].real / count,
            imag: shgSum[z][s].imag / count,
            phaseMatch: this.shgMetrics[z][s].phaseMatch / count,
            chiRatio: this.shgMetrics[z][s].chiRatio / count,
          };

          // Current lens parameters in this zone/sector
          // (Should already exist as an object in lensModes[z][s])
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

          // 1) Accumulate ADAM moments
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

          // 2) Compute bias-corrected
          const mBase_hat = this.adamState.m_baseIndex[z][s] * bc1;
          const vBase_hat = this.adamState.v_baseIndex[z][s] * bc2;

          const mDisp_hat = this.adamState.m_dispersion[z][s] * bc1;
          const vDisp_hat = this.adamState.v_dispersion[z][s] * bc2;

          const mChi2_hat = this.adamState.m_chi2[z][s] * bc1;
          const vChi2_hat = this.adamState.v_chi2[z][s] * bc2;

          const mChi3_hat = this.adamState.m_chi3[z][s] * bc1;
          const vChi3_hat = this.adamState.v_chi3[z][s] * bc2;

          // 3) Final param update
          const lr = this.config.learningRate;
          lensParams.baseIndex +=
            (lr * mBase_hat) / (Math.sqrt(vBase_hat) + this.adamState.epsilon);
          lensParams.dispersion +=
            (lr * mDisp_hat) / (Math.sqrt(vDisp_hat) + this.adamState.epsilon);
          lensParams.chi2 +=
            (lr * mChi2_hat) / (Math.sqrt(vChi2_hat) + this.adamState.epsilon);
          lensParams.chi3 +=
            (lr * mChi3_hat) / (Math.sqrt(vChi3_hat) + this.adamState.epsilon);

          totalLoss += loss;
          countUpdates++;
        }
      }
    }

    // (C) Record progress metric if any updates happened
    if (countUpdates > 0) {
      const avgLoss = totalLoss / countUpdates;
      this.progressHistory.push({
        iteration: this.adamState.t,
        metric: avgLoss,
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
          // Each pixel inherits the zone/sector's lens parameters
          const mode = this.lensModes[coords.zone][coords.sector];
          data[idx + 0] = mode.baseIndex;
          data[idx + 1] = mode.dispersion;
          data[idx + 2] = mode.chi2;
          data[idx + 3] = mode.chi3;
        } else {
          // Outside lens radius => e.g. vacuum or zero
          data[idx + 0] = 1.0; // baseline
          data[idx + 1] = 0.0; // no dispersion
          data[idx + 2] = 0.0; // no chi2
          data[idx + 3] = 0.0; // no chi3
        }
      }
    }

    return data;
  }

  getProgress() {
    // Return the array of recorded progress steps
    return this.progressHistory;
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
