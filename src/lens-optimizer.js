import { getUpdateStrategies } from "../update-strategies.js";
import { polar } from "./util.js";

export class LensOptimizer {
  constructor(config) {
    this.config = config;

    // Initialize lens modes
    this.lensModes = Array(config.fresnelZones)
      .fill()
      .map(() =>
        Array(config.numSectors)
          .fill()
          .map(() => ({ real: 1.0, imag: 0.0 })),
      );

    // Initialize Adam optimizer state
    this.adamState = adam.createState([config.fresnelZones, config.numSectors]);

    this.updateStrategies = getUpdateStrategies(config);

    // Store progress metrics over time
    this.progressHistory = [];
  }

  updateLens(fundamentalData, shgData, width, height) {
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

    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    // Accumulate
    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
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
            const fundReal = fundamentalData[idx];
            const fundImag = fundamentalData[idx + 1];
            const shgReal = shgData[idx];
            const shgImag = shgData[idx + 1];

            fundSum[coords.zone][coords.sector].real += fundReal;
            fundSum[coords.zone][coords.sector].imag += fundImag;
            shgSum[coords.zone][coords.sector].real += shgReal;
            shgSum[coords.zone][coords.sector].imag += shgImag;
            modeCount[coords.zone][coords.sector]++;
          }
        }
      }
    }

    // Update lens
    this.adamState.t += 1;
    const { bc1, bc2 } = adam.getBiasCorrection(this.adamState);
    const chosenStrategy = this.config.updateStrategy;
    const updateStrategy = this.updateStrategies[chosenStrategy];

    if (!updateStrategy)
      alert("couldn't find update strategy: " + this.config.updateStrategy);

    let totalUpdateMagnitude = 0;
    let countUpdates = 0;

    for (let z = 0; z < this.config.fresnelZones; z++) {
      for (let s = 0; s < this.config.numSectors; s++) {
        if (modeCount[z][s] > 0) {
          const fund = {
            real: fundSum[z][s].real / modeCount[z][s],
            imag: fundSum[z][s].imag / modeCount[z][s],
          };
          const shg = {
            real: shgSum[z][s].real / modeCount[z][s],
            imag: shgSum[z][s].imag / modeCount[z][s],
          };

          let update;
          if (updateStrategy.length > 3) {
            // hierarchicalZones or others needing zone info
            update = updateStrategy(
              fund,
              shg,
              modeCount[z][s],
              z,
              s,
              this.lensModes,
            );
          } else {
            update = updateStrategy(fund, shg, modeCount[z][s]);
          }

          // Update momentum estimates
          this.adamState.m_real[z][s] =
            this.adamState.beta1 * this.adamState.m_real[z][s] +
            (1 - this.adamState.beta1) * update.real;
          this.adamState.v_real[z][s] =
            this.adamState.beta2 * this.adamState.v_real[z][s] +
            (1 - this.adamState.beta2) * update.real * update.real;

          this.adamState.m_imag[z][s] =
            this.adamState.beta1 * this.adamState.m_imag[z][s] +
            (1 - this.adamState.beta1) * update.imag;
          this.adamState.v_imag[z][s] =
            this.adamState.beta2 * this.adamState.v_imag[z][s] +
            (1 - this.adamState.beta2) * update.imag * update.imag;

          // Compute bias-corrected estimates
          const m_real_hat = this.adamState.m_real[z][s] * bc1;
          const v_real_hat = this.adamState.v_real[z][s] * bc2;
          const m_imag_hat = this.adamState.m_imag[z][s] * bc1;
          const v_imag_hat = this.adamState.v_imag[z][s] * bc2;

          // Compute updates to modes
          const deltaReal =
            (this.config.learningRate * m_real_hat) /
            (Math.sqrt(v_real_hat) + this.adamState.epsilon);
          const deltaImag =
            (this.config.learningRate * m_imag_hat) /
            (Math.sqrt(v_imag_hat) + this.adamState.epsilon);

          this.lensModes[z][s].real += deltaReal;
          this.lensModes[z][s].imag += deltaImag;

          // Accumulate update magnitude for progress metric
          const updateMag = Math.sqrt(
            deltaReal * deltaReal + deltaImag * deltaImag,
          );
          totalUpdateMagnitude += updateMag;
          countUpdates++;
        }
      }
    }

    // Record progress metric: average update magnitude per iteration
    if (countUpdates > 0) {
      const avgUpdateMag = totalUpdateMagnitude / countUpdates;
      this.progressHistory.push({
        iteration: this.adamState.t,
        metric: avgUpdateMag,
      });
    }

    return this.lensModes;
  }

  getLensData(width, height) {
    const data = new Float32Array(width * height * 4);
    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    for (let y = 0; y < height; y++) {
      for (let x = 0; x < width; x++) {
        const coords = polar.getZoneAndSector(
          x,
          y,
          centerX,
          centerY,
          this.config.lensRadius,
          this.config.fresnelZones,
          this.config.numSectors,
        );

        const idx = (y * width + x) * 4;
        if (coords) {
          const mode = this.lensModes[coords.zone][coords.sector];
          data[idx] = mode.real;
          data[idx + 1] = mode.imag;
        } else {
          data[idx] = 1.0;
          data[idx + 1] = 0.0;
        }
        data[idx + 2] = 0.0;
        data[idx + 3] = 1.0;
      }
    }

    return data;
  }

  getProgress() {
    // Return a copy or a reference to the current progress history
    return this.progressHistory;
  }
}

const adam = {
  createState(dimensions) {
    return {
      beta1: 0.9,
      beta2: 0.999,
      epsilon: 1e-8,
      m_real: Array(dimensions[0])
        .fill()
        .map(() => Array(dimensions[1]).fill(0)),
      m_imag: Array(dimensions[0])
        .fill()
        .map(() => Array(dimensions[1]).fill(0)),
      v_real: Array(dimensions[0])
        .fill()
        .map(() => Array(dimensions[1]).fill(0)),
      v_imag: Array(dimensions[0])
        .fill()
        .map(() => Array(dimensions[1]).fill(0)),
      t: 0,
    };
  },

  getBiasCorrection(state) {
    const bc1 = 1.0 / (1.0 - Math.pow(state.beta1, state.t));
    const bc2 = 1.0 / (1.0 - Math.pow(state.beta2, state.t));
    return { bc1, bc2 };
  },
};
