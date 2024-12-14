export const getUpdateStrategies = (config) => {
  return {
    original: (fund, shg, modeCount) => {
      const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
      const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);
      let phaseDiff = shgPhase - 2 * fundPhase;
      while (phaseDiff > Math.PI) phaseDiff -= 2 * Math.PI;
      while (phaseDiff < -Math.PI) phaseDiff += 2 * Math.PI;

      const intensityRatio =
        fundIntensity > 1e-10 ? Math.min(shgIntensity / fundIntensity, 10) : 0;

      const hillPower = config.hillPower || 1.0;
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

      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);
      const expectedShgAmp = fundAmp * fundAmp;

      const coupling = {
        real: (shgAmp - expectedShgAmp) * Math.cos(phaseMismatch),
        imag: (shgAmp - expectedShgAmp) * Math.sin(phaseMismatch),
      };

      const hillPower = config.hillPower || 1.0;
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

      const hillPower = config.hillPower || 1.0;
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

      const hillPower = config.hillPower || 1.0;
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
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
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
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
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

    constructiveInterference: function (
      fund,
      shg,
      modeCount,
      zoneIndex,
      sectorIndex,
      lensModes,
    ) {
      // Basic field measures
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);

      // Phase matching term
      let phaseMismatch = shgPhase - 2.0 * fundPhase;
      while (phaseMismatch > Math.PI) phaseMismatch -= 2.0 * Math.PI;
      while (phaseMismatch < -Math.PI) phaseMismatch += 2.0 * Math.PI;

      // Compute interference efficiency
      // This measures how well the phases align for energy transfer
      const interferenceQuality = Math.cos(phaseMismatch);

      // Compute relative field strength ratio
      // We want this near 1 for optimal energy transfer
      const fieldRatio = (shgAmp + 1e-10) / (fundAmp * fundAmp + 1e-10);
      const balanceFactor = Math.exp(-Math.abs(1 - fieldRatio));

      // Zone-dependent phase sensitivity
      // Outer zones need more precise phase control
      const zoneWeight = Math.pow(1 + zoneIndex, 0.5) / config.fresnelZones;

      // Combine factors to determine update magnitude
      const updateMagnitude =
        zoneWeight * balanceFactor * Math.abs(interferenceQuality);

      // Update direction based on both phase mismatch and amplitude balance
      const updatePhase =
        phaseMismatch * (1 - balanceFactor) +
        Math.PI * (1 - interferenceQuality);

      return {
        real: updateMagnitude * Math.cos(updatePhase),
        imag: updateMagnitude * Math.sin(updatePhase),
      };
    },

    angularMomentum: (
      fund,
      shg,
      modeCount,
      zoneIndex,
      sectorIndex,
      lensModes,
    ) => {
      const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
      const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

      // Basic phase matching terms
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);
      let phaseMismatch = shgPhase - 2 * fundPhase;
      while (phaseMismatch > Math.PI) phaseMismatch -= 2 * Math.PI;
      while (phaseMismatch < -Math.PI) phaseMismatch += 2 * Math.PI;

      // Get neighboring sectors
      const numSectors = lensModes[0].length;
      const prevSector = (sectorIndex - 1 + numSectors) % numSectors;
      const nextSector = (sectorIndex + 1) % numSectors;

      // Get modes of neighboring sectors
      const prevMode = lensModes[zoneIndex][prevSector];
      const nextMode = lensModes[zoneIndex][nextSector];

      // Calculate phase gradients between sectors
      let prevPhaseDiff =
        Math.atan2(prevMode.imag, prevMode.real) -
        Math.atan2(
          lensModes[zoneIndex][sectorIndex].imag,
          lensModes[zoneIndex][sectorIndex].real,
        );
      let nextPhaseDiff =
        Math.atan2(nextMode.imag, nextMode.real) -
        Math.atan2(
          lensModes[zoneIndex][sectorIndex].imag,
          lensModes[zoneIndex][sectorIndex].real,
        );

      // Normalize phase differences
      while (prevPhaseDiff > Math.PI) prevPhaseDiff -= 2 * Math.PI;
      while (prevPhaseDiff < -Math.PI) prevPhaseDiff += 2 * Math.PI;
      while (nextPhaseDiff > Math.PI) nextPhaseDiff -= 2 * Math.PI;
      while (nextPhaseDiff < -Math.PI) nextPhaseDiff += 2 * Math.PI;

      // Add rotational bias (positive for clockwise, negative for counter-clockwise)
      const rotationalBias = 0.2;
      const targetPhaseDiff = rotationalBias;

      // Calculate angular momentum terms
      const angularCorrection = {
        prev: prevPhaseDiff - targetPhaseDiff,
        next: nextPhaseDiff + targetPhaseDiff,
      };

      // Weight for radial position - stronger effect at larger radii
      const radialWeight = (zoneIndex + 1) / lensModes.length;

      // Coupling strength based on intensities
      const hillPower = config.hillPower || 1.0;
      const couplingStrength = shgIntensity / (fundIntensity + 1e-10);
      const weight =
        Math.pow(couplingStrength, hillPower) /
        (1 + Math.pow(couplingStrength, hillPower));

      // Combine phase matching and angular momentum terms
      const baseUpdate = {
        real: weight * Math.cos(phaseMismatch),
        imag: weight * Math.sin(phaseMismatch),
      };

      // Add weighted angular momentum correction
      const angularWeight = 0.95 * radialWeight * weight;
      const angularUpdate = {
        real:
          -angularWeight *
          (Math.cos(angularCorrection.prev) + Math.cos(angularCorrection.next)),
        imag:
          -angularWeight *
          (Math.sin(angularCorrection.prev) + Math.sin(angularCorrection.next)),
      };

      return {
        real: baseUpdate.real + angularUpdate.real,
        imag: baseUpdate.imag + angularUpdate.imag,
      };
    },

    aestheticChaos: (
      fund,
      shg,
      modeCount,
      zoneIndex,
      sectorIndex,
      lensModes,
    ) => {
      const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
      const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

      // Get phases
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);

      // Calculate various competing terms

      // 1. Anti-uniformity term: encourage intensity oscillations
      const targetIntensity =
        Math.sin(zoneIndex * 0.7) * Math.cos(sectorIndex * 0.15) * 0.5 + 0.5;
      const intensityDiff =
        shgIntensity / (fundIntensity + 1e-10) - targetIntensity;

      // 2. Phase pattern term: create spiral-like patterns
      const targetPhase = (zoneIndex * 0.2 + sectorIndex * 0.1) % (2 * Math.PI);
      let phaseDiff = shgPhase - targetPhase;
      while (phaseDiff > Math.PI) phaseDiff -= 2 * Math.PI;
      while (phaseDiff < -Math.PI) phaseDiff += 2 * Math.PI;

      // 3. Neighbor coupling with instability
      const numSectors = lensModes[0].length;
      const prevSector = (sectorIndex - 1 + numSectors) % numSectors;
      const nextSector = (sectorIndex + 1) % numSectors;

      const prevMode = lensModes[zoneIndex][prevSector];
      const nextMode = lensModes[zoneIndex][nextSector];

      // Add slight instability by amplifying differences
      const neighborDiff = {
        real: (nextMode.real - prevMode.real) * Math.sin(zoneIndex * 0.4),
        imag: (nextMode.imag - prevMode.imag) * Math.cos(sectorIndex * 0.3),
      };

      // 4. Zone-based oscillation
      const zoneOscillation =
        Math.sin(zoneIndex * 0.8) * Math.cos(sectorIndex * 0.2);

      // Combine all terms with time-varying weights
      const t = Date.now() * 0.001; // Slow time variation
      const w1 = Math.sin(t * 0.1) * 0.5 + 0.5;
      const w2 = Math.cos(t * 0.15) * 0.5 + 0.5;

      // Final update combines everything
      const update = {
        real:
          -0.3 * intensityDiff * Math.cos(phaseDiff) +
          0.4 * neighborDiff.real * zoneOscillation * w1 +
          0.3 * Math.cos(targetPhase) * w2,
        imag:
          -0.3 * intensityDiff * Math.sin(phaseDiff) +
          0.4 * neighborDiff.imag * zoneOscillation * w1 +
          0.3 * Math.sin(targetPhase) * w2,
      };

      // Apply a "chaos factor" that increases with radius
      const chaosFactor = 0.5 + (zoneIndex / lensModes.length) * 0.5;

      return {
        real: update.real * chaosFactor,
        imag: update.imag * chaosFactor,
      };
    },

    conservativeUpdate: (
      fund,
      shg,
      modeCount,
      zoneIndex,
      sectorIndex,
      lensModes,
    ) => {
      // Calculate field intensities
      const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
      const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;
      const totalIntensity = fundIntensity + shgIntensity;

      // Calculate local phase relationships
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);
      let phaseDiff = shgPhase - 2 * fundPhase;
      while (phaseDiff > Math.PI) phaseDiff -= 2 * Math.PI;
      while (phaseDiff < -Math.PI) phaseDiff += 2 * Math.PI;

      // Get neighbor information
      const numSectors = lensModes[0].length;
      const prevSector = (sectorIndex - 1 + numSectors) % numSectors;
      const nextSector = (sectorIndex + 1) % numSectors;
      const prevMode = lensModes[zoneIndex][prevSector];
      const nextMode = lensModes[zoneIndex][nextSector];

      // Calculate average mode strength of neighbors
      const neighborAvg = {
        real: (prevMode.real + nextMode.real) * 0.5,
        imag: (prevMode.imag + nextMode.imag) * 0.5,
      };

      // Current mode
      const currentMode = lensModes[zoneIndex][sectorIndex];

      // Calculate deviation from neighbor average
      const deviation = {
        real: currentMode.real - neighborAvg.real,
        imag: currentMode.imag - neighborAvg.imag,
      };

      // Inverse scaling with total intensity (saturates response at high intensities)
      const intensityScaling = 1.0 / (1.0 + totalIntensity * 0.1);

      // Competitive inhibition: stronger fields reduce update magnitude
      const neighborIntensity = Math.sqrt(
        neighborAvg.real * neighborAvg.real +
          neighborAvg.imag * neighborAvg.imag,
      );
      const competitiveInhibition = 1.0 / (1.0 + neighborIntensity * 0.2);

      // Combine phase matching with neighbor coupling
      const update = {
        real:
          (-0.3 * deviation.real * competitiveInhibition +
            0.7 * Math.cos(phaseDiff)) *
          intensityScaling,
        imag:
          (-0.3 * deviation.imag * competitiveInhibition +
            0.7 * Math.sin(phaseDiff)) *
          intensityScaling,
      };

      // Add slight randomization that decreases with intensity
      const randomFactor = intensityScaling * 0.1;
      update.real += (Math.random() - 0.5) * randomFactor;
      update.imag += (Math.random() - 0.5) * randomFactor;

      return update;
    },

    freeEnergyMin: (
      fund,
      shg,
      modeCount,
      zoneIndex,
      sectorIndex,
      lensModes,
    ) => {
      const numZones = lensModes.length;
      const numSectors = lensModes[0].length;

      // Calculate local energies
      const fundIntensity = fund.real * fund.real + fund.imag * fund.imag;
      const shgIntensity = shg.real * shg.real + shg.imag * shg.imag;

      // Get neighboring sectors for local entropy calculation
      const prevSector = (sectorIndex - 1 + numSectors) % numSectors;
      const nextSector = (sectorIndex + 1) % numSectors;
      const neighborModes = [
        lensModes[zoneIndex][prevSector],
        lensModes[zoneIndex][nextSector],
      ];

      // Add zones above and below if they exist
      if (zoneIndex > 0) {
        neighborModes.push(lensModes[zoneIndex - 1][sectorIndex]);
      }
      if (zoneIndex < numZones - 1) {
        neighborModes.push(lensModes[zoneIndex + 1][sectorIndex]);
      }

      // Calculate local "entropy" based on field variation from neighbors
      let localEntropy = 0;
      let neighborAvg = { real: 0, imag: 0 };

      neighborModes.forEach((mode) => {
        neighborAvg.real += mode.real / neighborModes.length;
        neighborAvg.imag += mode.imag / neighborModes.length;
      });

      const currentMode = lensModes[zoneIndex][sectorIndex];
      const deviation = Math.sqrt(
        Math.pow(currentMode.real - neighborAvg.real, 2) +
          Math.pow(currentMode.imag - neighborAvg.imag, 2),
      );

      // Higher deviation = higher entropy
      localEntropy = Math.log(1 + deviation);

      // Include phase mismatch in energy term
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);
      let phaseMismatch = shgPhase - 2 * fundPhase;
      while (phaseMismatch > Math.PI) phaseMismatch -= 2 * Math.PI;
      while (phaseMismatch < -Math.PI) phaseMismatch += 2 * Math.PI;

      // Energy term includes both intensity and phase mismatch
      const localEnergy =
        (fundIntensity + shgIntensity) *
        (1 + Math.abs(phaseMismatch) / Math.PI);

      // Define our effective temperature (could be made dynamic)
      const T = 0.3 * (1 + zoneIndex / numZones); // Temperature increases with radius

      // Free energy gradient: dF = dE - T*dS
      // We want to move down the gradient
      const freeEnergyGradient = {
        real: -neighborAvg.real * (localEnergy - T * localEntropy),
        imag: -neighborAvg.imag * (localEnergy - T * localEntropy),
      };

      // Add a small stochastic term that scales with temperature
      const noise = T * 0.1;
      freeEnergyGradient.real += (Math.random() - 0.5) * noise;
      freeEnergyGradient.imag += (Math.random() - 0.5) * noise;

      return {
        real: freeEnergyGradient.real,
        imag: freeEnergyGradient.imag,
      };
    },
  };
};
