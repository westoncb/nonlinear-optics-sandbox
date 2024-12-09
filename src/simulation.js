import { polar, adam, stats } from "./util.js";

export class Simulation {
  constructor(config) {
    this.config = config;

    // Initialize lens modes
    this.lensModes = Array(config.fresnelZones)
      .fill()
      .map(() =>
        Array(config.numSectors)
          .fill()
          .map(() => ({
            real: 1.0,
            imag: 0.0,
          })),
      );

    // Initialize Adam optimizer state
    this.adamState = adam.createState([config.fresnelZones, config.numSectors]);

    // Measurement accumulators
    this.intensityRatios = Array(config.fresnelZones)
      .fill()
      .map(() => Array(config.numSectors).fill(0));
    this.phaseDiffs = Array(config.fresnelZones)
      .fill()
      .map(() => Array(config.numSectors).fill(0));
    this.modeCount = Array(config.fresnelZones)
      .fill()
      .map(() => Array(config.numSectors).fill(0));
  }

  updateLens(fundamentalData, shgData, width, height) {
    // Reset accumulators
    for (let z = 0; z < this.config.fresnelZones; z++) {
      for (let s = 0; s < this.config.numSectors; s++) {
        this.intensityRatios[z][s] = 0;
        this.phaseDiffs[z][s] = 0;
        this.modeCount[z][s] = 0;
      }
    }

    const centerX = Math.floor(width / 2);
    const centerY = Math.floor(height / 2);

    // Gather measurements
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

        if (coords) {
          const idx = (y * width + x) * 4;

          const fundReal = fundamentalData[idx];
          const fundImag = fundamentalData[idx + 1];
          const fundIntensity = fundReal * fundReal + fundImag * fundImag;

          const shgReal = shgData[idx];
          const shgImag = shgData[idx + 1];
          const shgIntensity = shgReal * shgReal + shgImag * shgImag;

          const fundPhase = Math.atan2(fundImag, fundReal);
          const shgPhase = Math.atan2(shgImag, shgReal);
          let phaseDiff = shgPhase - 2 * fundPhase;
          while (phaseDiff > Math.PI) phaseDiff -= 2 * Math.PI;
          while (phaseDiff < -Math.PI) phaseDiff += 2 * Math.PI;

          const intensityRatio =
            fundIntensity > 1e-10
              ? Math.min(shgIntensity / fundIntensity, 10)
              : 0;

          this.intensityRatios[coords.zone][coords.sector] += intensityRatio;
          this.phaseDiffs[coords.zone][coords.sector] += phaseDiff;
          this.modeCount[coords.zone][coords.sector]++;
        }
      }
    }

    // Update lens modes using Adam optimizer
    this.adamState.t += 1;
    const { bc1, bc2 } = adam.getBiasCorrection(this.adamState);

    for (let z = 0; z < this.config.fresnelZones; z++) {
      for (let s = 0; s < this.config.numSectors; s++) {
        if (this.modeCount[z][s] > 0) {
          const avgIntensityRatio =
            this.intensityRatios[z][s] / this.modeCount[z][s];
          const avgPhaseDiff = this.phaseDiffs[z][s] / this.modeCount[z][s];

          const magnitude = Math.tanh(avgIntensityRatio);
          const updateReal = magnitude * Math.cos(avgPhaseDiff);
          const updateImag = magnitude * Math.sin(avgPhaseDiff);

          // Update momentum estimates
          this.adamState.m_real[z][s] =
            this.adamState.beta1 * this.adamState.m_real[z][s] +
            (1 - this.adamState.beta1) * updateReal;
          this.adamState.v_real[z][s] =
            this.adamState.beta2 * this.adamState.v_real[z][s] +
            (1 - this.adamState.beta2) * updateReal * updateReal;

          this.adamState.m_imag[z][s] =
            this.adamState.beta1 * this.adamState.m_imag[z][s] +
            (1 - this.adamState.beta1) * updateImag;
          this.adamState.v_imag[z][s] =
            this.adamState.beta2 * this.adamState.v_imag[z][s] +
            (1 - this.adamState.beta2) * updateImag * updateImag;

          // Compute bias-corrected estimates and update modes
          const m_real_hat = this.adamState.m_real[z][s] * bc1;
          const v_real_hat = this.adamState.v_real[z][s] * bc2;
          const m_imag_hat = this.adamState.m_imag[z][s] * bc1;
          const v_imag_hat = this.adamState.v_imag[z][s] * bc2;

          this.lensModes[z][s].real +=
            (this.config.learningRate * m_real_hat) /
            (Math.sqrt(v_real_hat) + this.adamState.epsilon);
          this.lensModes[z][s].imag +=
            (this.config.learningRate * m_imag_hat) /
            (Math.sqrt(v_imag_hat) + this.adamState.epsilon);
        }
      }
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
}
