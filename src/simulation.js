import { polar, adam } from "./util.js";

export class Simulation {
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

    // Update strategies
    this.updateStrategies = {
      original: (fund, shg, modeCount) => {
        const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
        const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

        const fundPhase = Math.atan2(fund.imag, fund.real);
        const shgPhase = Math.atan2(shg.imag, shg.real);
        let phaseDiff = shgPhase - 2 * fundPhase;
        while (phaseDiff > Math.PI) phaseDiff -= 2 * Math.PI;
        while (phaseDiff < -Math.PI) phaseDiff += 2 * Math.PI;

        const intensityRatio =
          fundIntensity > 1e-10
            ? Math.min(shgIntensity / fundIntensity, 10)
            : 0;

        const hillPower = this.config.hillPower || 1.0;
        const magnitude =
          Math.pow(intensityRatio, hillPower) /
          (1 + Math.pow(intensityRatio, hillPower));

        return {
          real: magnitude * Math.cos(phaseDiff),
          imag: magnitude * Math.sin(phaseDiff),
        };
      },

      phaseMatching: (fund, shg, modeCount) => {
        const fundPhase = Math.atan2(fund.imag, fund.real);
        const shgPhase = Math.atan2(shg.imag, shg.real);
        const phaseMismatch = shgPhase - 2 * fundPhase;

        const fundAmp = Math.sqrt(
          fund.real * fund.real + fund.imag * fund.imag,
        );
        const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);
        const expectedShgAmp = fundAmp * fundAmp;

        const coupling = {
          real: (shgAmp - expectedShgAmp) * Math.cos(phaseMismatch),
          imag: (shgAmp - expectedShgAmp) * Math.sin(phaseMismatch),
        };

        const hillPower = this.config.hillPower || 1.0;
        const couplingMagnitude = Math.sqrt(
          coupling.real * coupling.real + coupling.imag * coupling.imag,
        );
        const magnitude =
          Math.pow(couplingMagnitude, hillPower) /
          (1 + Math.pow(couplingMagnitude, hillPower));

        const couplingPhase = Math.atan2(coupling.imag, coupling.real);

        return {
          real: magnitude * Math.cos(couplingPhase),
          imag: magnitude * Math.sin(couplingPhase),
        };
      },

      multiCoupling: (fund, shg, modeCount) => {
        const fundSquared = {
          real: fund.real * fund.real - fund.imag * fund.imag,
          imag: 2 * fund.real * fund.imag,
        };

        const directCoupling = {
          real: shg.real - fundSquared.real,
          imag: shg.imag - fundSquared.imag,
        };

        const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;
        const fundSquaredIntensity =
          fundSquared.real * fundSquared.real +
          fundSquared.imag * fundSquared.imag;
        const intensityCoupling = {
          real: shgIntensity - fundSquaredIntensity,
          imag: 0,
        };

        const coupling = {
          real: 0.7 * directCoupling.real + 0.3 * intensityCoupling.real,
          imag: 0.7 * directCoupling.imag + 0.3 * intensityCoupling.imag,
        };

        const hillPower = this.config.hillPower || 1.0;
        const couplingMagnitude = Math.sqrt(
          coupling.real * coupling.real + coupling.imag * coupling.imag,
        );
        const magnitude =
          Math.pow(couplingMagnitude, hillPower) /
          (1 + Math.pow(couplingMagnitude, hillPower));

        const couplingPhase = Math.atan2(coupling.imag, coupling.real);

        return {
          real: magnitude * Math.cos(couplingPhase),
          imag: magnitude * Math.sin(couplingPhase),
        };
      },

      hierarchicalZones: (
        fund,
        shg,
        modeCount,
        zoneIndex,
        sectorIndex,
        lensModes,
      ) => {
        const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
        const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

        const fundPhase = Math.atan2(fund.imag, fund.real);
        const shgPhase = Math.atan2(shg.imag, shg.real);
        const phaseMismatch = shgPhase - 2 * fundPhase;

        let innerZoneInfluence = { real: 0, imag: 0 };
        let outerZoneInfluence = { real: 0, imag: 0 };

        if (zoneIndex > 0) {
          const innerMode = lensModes[zoneIndex - 1][sectorIndex];
          innerZoneInfluence = {
            real: 0.3 * innerMode.real,
            imag: 0.3 * innerMode.imag,
          };
        }

        if (zoneIndex < lensModes.length - 1) {
          const outerMode = lensModes[zoneIndex + 1][sectorIndex];
          outerZoneInfluence = {
            real: 0.2 * outerMode.real,
            imag: 0.2 * outerMode.imag,
          };
        }

        const hillPower = this.config.hillPower || 1.0;
        const localCoupling = shgIntensity / (fundIntensity + 1e-10);
        const hierarchicalWeight =
          Math.pow(localCoupling, hillPower) /
          (1 + Math.pow(localCoupling, hillPower));

        return {
          real:
            hierarchicalWeight *
            (Math.cos(phaseMismatch) +
              innerZoneInfluence.real +
              outerZoneInfluence.real),
          imag:
            hierarchicalWeight *
            (Math.sin(phaseMismatch) +
              innerZoneInfluence.imag +
              outerZoneInfluence.imag),
        };
      },

      o1FirstOptimal: function (
        fund,
        shg,
        modeCount,
        zoneIndex,
        sectorIndex,
        lensModes,
      ) {
        // Extract configuration and parameters
        const hillPower = config.hillPower || 1.0;

        // Compute fundamental amplitude and SHG amplitude
        const fundAmp = Math.sqrt(
          fund.real * fund.real + fund.imag * fund.imag,
        );
        const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

        // Compute the desired SHG field = (fund)^2
        // SHG_desired = fund² = (fund.real + i fund.imag)²
        // = (fund.real² - fund.imag²) + i(2 * fund.real * fund.imag)
        const desiredShg = {
          real: fund.real * fund.real - fund.imag * fund.imag,
          imag: 2 * fund.real * fund.imag,
        };

        // Error is the difference between actual SHG and desired SHG
        const error = {
          real: shg.real - desiredShg.real,
          imag: shg.imag - desiredShg.imag,
        };

        const errorMag = Math.sqrt(
          error.real * error.real + error.imag * error.imag,
        );

        // Define a "coupling" factor as a ratio of SHG intensity to fundamental intensity
        const fundIntensity = fundAmp * fundAmp;
        const shgIntensity = shgAmp * shgAmp;
        const coupling = shgIntensity / (fundIntensity + 1e-10);

        // Apply a logistic-like weighting based on coupling intensity
        const couplingWeight =
          Math.pow(coupling, hillPower) / (1 + Math.pow(coupling, hillPower));

        // Similarly weight the error to avoid over-amplifying large deviations
        const errorWeight =
          Math.pow(errorMag, hillPower) / (1 + Math.pow(errorMag, hillPower));

        // Direction of update: move against the error vector (gradient-like step)
        // Normalize the error direction and scale by combined weights
        const errorAngle = Math.atan2(error.imag, error.real);
        const updateReal = -errorWeight * couplingWeight * Math.cos(errorAngle);
        const updateImag = -errorWeight * couplingWeight * Math.sin(errorAngle);

        return { real: updateReal, imag: updateImag };
      },

      o1SecondOptimal: function (
        fund,
        shg,
        modeCount,
        zoneIndex,
        sectorIndex,
        lensModes,
      ) {
        const hillPower = config.hillPower || 1.0;

        // Fundamental and SHG amplitudes
        const fundAmp = Math.sqrt(
          fund.real * fund.real + fund.imag * fund.imag,
        );
        const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

        // Desired SHG amplitude is fundAmp²
        const desiredShgAmp = fundAmp * fundAmp;

        // Compute amplitude mismatch and phase mismatch
        const ampMismatch = shgAmp - desiredShgAmp;
        const fundPhase = Math.atan2(fund.imag, fund.real);
        const shgPhase = Math.atan2(shg.imag, shg.real);
        let phaseMismatch = shgPhase - 2.0 * fundPhase;
        // Normalize phase into (-π, π)
        while (phaseMismatch > Math.PI) phaseMismatch -= 2.0 * Math.PI;
        while (phaseMismatch < -Math.PI) phaseMismatch += 2.0 * Math.PI;

        // Form a complex error signal: error = SHG - fund²
        const desiredShg = {
          real: fundAmp * fundAmp * Math.cos(2 * fundPhase),
          imag: fundAmp * fundAmp * Math.sin(2 * fundPhase),
        };
        const error = {
          real: shg.real - desiredShg.real,
          imag: shg.imag - desiredShg.imag,
        };
        const errorMag = Math.sqrt(
          error.real * error.real + error.imag * error.imag,
        );

        // Compute coupling factor as before
        const fundIntensity = fundAmp * fundAmp;
        const shgIntensity = shgAmp * shgAmp;
        const coupling = shgIntensity / (fundIntensity + 1e-10);

        // Weighting function to stabilize updates
        const weight =
          Math.pow(coupling, hillPower) / (1 + Math.pow(coupling, hillPower));

        // Another weighting based on error magnitude
        // Ensures that we don't overcorrect on large errors or underreact to small ones
        const errorWeight =
          Math.pow(errorMag, hillPower) / (1 + Math.pow(errorMag, hillPower));

        // The baseline direction: move opposite to the error vector
        // This pushes SHG towards desiredShg.
        const errorAngle = Math.atan2(error.imag, error.real);

        // To incorporate PDE insight (phase matching focus), we can slightly bias
        // the correction towards reducing phaseMismatch first:
        // If we let the update direction be errorAngle shifted by a fraction of -phaseMismatch,
        // it prioritizes aligning phase.
        //
        // For example, rotate the error vector by -0.5 * phaseMismatch to emphasize phase correction.
        // (This is a heuristic choice; tuning the factor may help.)
        const correctionPhase = errorAngle - 0.5 * phaseMismatch;
        const updateMagnitude = errorWeight * weight;

        const update = {
          real: -updateMagnitude * Math.cos(correctionPhase),
          imag: -updateMagnitude * Math.sin(correctionPhase),
        };

        // This approach:
        // - Uses the direct error between actual and desired SHG fields.
        // - Weighs updates by nonlinear coupling intensity and error magnitude.
        // - Adjusts update direction to prioritize reducing phase mismatch.
        //
        // Over iterations, this should guide the lens modes toward configurations
        // that yield well phase-matched and correctly amplified SHG fields,
        // given how the PDE utilizes the lens in modifying the wave propagation.

        return update;
      },
    };

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

    // Update lens
    this.adamState.t += 1;
    const { bc1, bc2 } = adam.getBiasCorrection(this.adamState);
    const chosenStrategy = this.config.updateStrategy || "original";
    const updateStrategy =
      this.updateStrategies[chosenStrategy] || this.updateStrategies.original;

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
