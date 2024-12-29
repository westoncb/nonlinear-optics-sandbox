const computeAmplitudeGradients = (real, imag) => {
  const amp = Math.sqrt(real * real + imag * imag);
  if (amp < 1e-10) return { real: 0, imag: 0 };
  return {
    real: real / amp,
    imag: imag / amp,
  };
};

const computePhaseGradients = (real, imag) => {
  const r2i2 = real * real + imag * imag;
  if (r2i2 < 1e-10) return { real: 0, imag: 0 };
  // ∂(atan2(i,r))/∂r = -i/(r² + i²)
  // ∂(atan2(i,r))/∂i = r/(r² + i²)
  return {
    real: -imag / r2i2,
    imag: real / r2i2,
  };
};

export const getUpdateStrategies = (config) => {
  return {
    phaseMatch: (fund, shg, count, zone, sector, lensParams, config) => {
      // 1. Calculate basic field measures
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

      // 2. Consider spatial position in lens
      const r = zone / config.fresnelZones;
      const theta = (2 * Math.PI * sector) / config.numSectors;

      // 3. Phase matching quality measures
      const conversionEfficiency = shgAmp / (fundAmp + 1e-6);
      const phaseQuality = shg.phaseMatch;

      // 4. Position-dependent objectives
      const zoneWeight = {
        chi2: Math.exp(-2 * r * r),
        chi3: Math.sin(Math.PI * r),
        base: r,
      };

      // 5. Calculate composite loss
      const loss =
        1.0 * fundAmp +
        0.5 * fund.kerrStrength -
        2.0 * conversionEfficiency -
        1.0 * phaseQuality +
        0.1 * Math.sin(4 * theta);

      // 6. Convert the original updates to gradients
      const baseScale = 0.015;
      return {
        gradients: {
          baseIndex: baseScale * loss * zoneWeight.base,
          dispersion: 0.2 * baseScale * loss * (1 - r),
          chi2: 2.0 * baseScale * loss * zoneWeight.chi2,
          chi3: baseScale * loss * zoneWeight.chi3,
        },
        loss,
      };
    },

    phaseMatchAsymmetric: (
      fund,
      shg,
      count,
      zone,
      sector,
      lensParams,
      config,
    ) => {
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

      const r = zone / config.fresnelZones;
      const theta = (2 * Math.PI * sector) / config.numSectors;

      // Create alternating zones that favor either fundamental or SHG
      const zonePreference = Math.sin(3 * Math.PI * r); // oscillates between -1 and 1

      const zoneWeight = {
        chi2: Math.exp(-2 * r * r) * (1 + 0.5 * zonePreference),
        chi3: Math.sin(Math.PI * r) * (1 - 0.5 * zonePreference),
        base: r * (1 + 0.3 * Math.cos(6 * theta)),
      };

      // Asymmetric loss function
      const loss =
        (zonePreference > 0 ? 2.0 : 0.5) * fundAmp +
        0.5 * fund.kerrStrength -
        (zonePreference < 0 ? 2.0 : 0.5) * shgAmp -
        1.0 * shg.phaseMatch +
        0.1 * Math.sin(4 * theta);

      const baseScale = 0.015;
      return {
        gradients: {
          baseIndex: baseScale * loss * zoneWeight.base,
          dispersion: 0.2 * baseScale * loss * (1 - r),
          chi2: 2.0 * baseScale * loss * zoneWeight.chi2,
          chi3: baseScale * loss * zoneWeight.chi3,
        },
        loss,
      };
    },

    SHGTest: (fund, shg, count, zone, sector, lensParams, config) => {
      // Compute both amplitude and phase gradients
      const fundAmpGrads = computeAmplitudeGradients(fund.real, fund.imag);
      const shgAmpGrads = computeAmplitudeGradients(shg.real, shg.imag);
      const fundPhaseGrads = computePhaseGradients(fund.real, fund.imag);
      const shgPhaseGrads = computePhaseGradients(shg.real, shg.imag);

      // Spatial factors (encourage structure formation)
      const zoneFactor = Math.exp(-zone / config.fresnelZones); // stronger in center
      const radialPhase = (2 * Math.PI * zone) / config.fresnelZones;

      // Compute actual phase values for phase matching terms
      const fundPhase = Math.atan2(fund.imag, fund.real);
      const shgPhase = Math.atan2(shg.imag, shg.real);
      const phaseMismatch = shgPhase - 2 * fundPhase;

      // Basic amplitude gradients
      const fundAmpGrad =
        fundAmpGrads.real * fund.real + fundAmpGrads.imag * fund.imag;
      const shgAmpGrad =
        shgAmpGrads.real * shg.real + shgAmpGrads.imag * shg.imag;

      // Phase-sensitive terms
      const fundPhaseGrad =
        fundPhaseGrads.real * fund.real + fundPhaseGrads.imag * fund.imag;
      const shgPhaseGrad =
        shgPhaseGrads.real * shg.real + shgPhaseGrads.imag * shg.imag;

      // Weights with spatial dependence
      const w_fund = 1.0 * zoneFactor;
      const w_shg = 2.0 * (1 - zoneFactor); // stronger SHG in outer zones
      const w_kerr = 0.5 * Math.exp(-fund.kerrStrength); // adaptive Kerr weight
      const w_pm = 1.0 * (1 + 0.5 * Math.cos(radialPhase)); // phase-dependent

      // Enhanced gradient terms
      const fundTerm = w_fund * (fundAmpGrad + 0.2 * fundPhaseGrad);
      const shgTerm = -w_shg * (shgAmpGrad + 0.2 * shgPhaseGrad);
      const kerrTerm =
        w_kerr * fund.kerrStrength * (1 + Math.abs(fundPhaseGrad));
      const phaseTerm =
        -w_pm * (shg.phaseMatch + 0.1 * Math.cos(phaseMismatch));

      // Base gradient combines all terms
      const baseGrad = fundTerm + shgTerm + kerrTerm + phaseTerm;

      // Parameter-specific modifications based on position and phase
      const baseIndex_mod = 1.0 + 0.2 * Math.cos(radialPhase);
      const dispersion_mod = 1.0 + 0.3 * Math.sin(phaseMismatch);
      const chi2_mod = 1.0 + 0.4 * zoneFactor * Math.cos(2 * fundPhase);
      const chi3_mod = 1.0 + 0.2 * (1 - zoneFactor) * Math.sin(shgPhase);

      return {
        gradients: {
          baseIndex: baseGrad * baseIndex_mod + 0.1 * shg.phaseMatch,
          dispersion: 0.5 * baseGrad * dispersion_mod - 0.2 * phaseMismatch,
          chi2: 2.0 * baseGrad * chi2_mod - w_shg * shgAmpGrad,
          chi3:
            baseGrad * chi3_mod +
            0.5 * kerrTerm * (1 - Math.abs(fundPhaseGrad)),
        },
        loss: fundTerm + shgTerm + kerrTerm + phaseTerm,
      };
    },
  };
};
