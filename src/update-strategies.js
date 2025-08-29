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

    // Strategy focusing purely on coherent phase matching
    EnhancedPhaseMatch: (
      fund,
      shg,
      count,
      zone,
      sector,
      lensParams,
      config,
    ) => {
      const fundAmp = Math.hypot(fund.real, fund.imag);
      const shgAmp = Math.hypot(shg.real, shg.imag);
      const conversionEfficiency = shgAmp / (fundAmp + 1e-6);

      const r = zone / config.fresnelZones;
      const theta = (2 * Math.PI * sector) / config.numSectors;
      const radialPhase = Math.PI * r;

      const phaseMismatch = Math.tanh(2.0 * Math.abs(shg.phaseMatch));

      const loss =
        0.8 * fundAmp + // penalize fundamental buildup
        0.4 * fund.kerrStrength - // penalize excessive Kerr nonlinearities
        2.5 * conversionEfficiency + // strongly reward SHG efficiency
        1.5 * phaseMismatch + // strongly penalize phase mismatch
        0.05 * Math.sin(4 * theta); // angular symmetry-breaking

      const baseScale = 0.01;

      const zoneWeight = {
        chi2: Math.exp(-3 * r * r), // stronger focus at center
        chi3: Math.sin(radialPhase), // smoothly varying radial modulation
        base: 0.5 * (1 - Math.cos(radialPhase)), // smooth radial growth outward
      };

      return {
        gradients: {
          baseIndex: baseScale * loss * zoneWeight.base,
          dispersion: baseScale * loss * (1 - r),
          chi2: baseScale * loss * zoneWeight.chi2,
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

    PhaseLockedSHG: (fund, shg, count, zone, sector, lensParams, config) => {
      // --- Derived field metrics ---
      const eps = 1e-6;
      const fAmp = Math.hypot(fund.real, fund.imag);
      const sAmp = Math.hypot(shg.real, shg.imag);

      const fPhase = Math.atan2(fund.imag, fund.real);
      const sPhase = Math.atan2(shg.imag, shg.real);

      // Δφ = φ2 - 2φ1 wrapped to [-π, π]
      let dphi = sPhase - 2 * fPhase;
      if (dphi > Math.PI) dphi -= 2 * Math.PI;
      if (dphi < -Math.PI) dphi += 2 * Math.PI;

      // Coherence and efficiency proxies
      const pm = Math.cos(dphi); // ∈ [-1,1], 1 = perfect lock
      const pm_pos = 0.5 * (pm + 1.0); // ∈ [0,1]
      const eta = Math.max(0, Math.min(1, sAmp / (fAmp * fAmp + eps)));

      // Spatial factors
      const Z = Math.max(1, (config.fresnelZones || 1) - 1);
      const r = zone / Math.max(1, config.fresnelZones || 1); // [0,1] center→edge
      const theta =
        (2 * Math.PI * sector) / Math.max(1, config.numSectors || 1);

      // Weights (can be tweaked via config.strategyWeights?.PhaseLockedSHG)
      const w = Object.assign(
        {
          w_pm: 1.0,
          w_eta: 1.2,
          w_kerr: 0.1, // loss weights
          kb: 0.8,
          kd: 0.6,
          kc2: 1.1,
          kc3: 0.8, // per-parameter gains
          lambda_base: 0.02,
          lambda_disp: 0.02, // L2 priors to keep params sane
        },
        (config.strategyWeights && config.strategyWeights.PhaseLockedSHG) || {},
      );

      // Loss for logging/ADAM scaling (lower is better)
      const loss =
        w.w_pm * (1.0 - pm) + // phase misalignment
        w.w_eta * (1.0 - eta) + // low SHG efficiency
        w.w_kerr * (fund.kerrStrength || 0); // strong Kerr

      // Base step scale (optionally tunable via config.strategyScales?.PhaseLockedSHG)
      const baseScale =
        (config.strategyScales && config.strategyScales.PhaseLockedSHG) || 0.01;

      // --- Gradients (heuristic but directional) ---
      // Use sin(Δφ) to move mismatch toward zero; modulate spatially.
      const sinD = Math.sin(dphi);

      // Index: push more in outer zones; add mild azimuthal symmetry-breaking
      const g_base =
        w.kb * sinD * (0.5 + 0.5 * r) * (1.0 + 0.15 * Math.cos(2 * theta)) -
        // L2 prior toward baseIndex ≈ 1.0
        w.lambda_base * (lensParams.baseIndex - 1.0);

      // Dispersion: adjust more near center (acts like "fine tuner")
      const g_disp =
        w.kd * sinD * (1.0 - r) -
        // L2 prior toward dispersion ≈ 0
        w.lambda_disp * lensParams.dispersion;

      // χ(2): reward when phase is aligned and efficiency is lagging
      const g_chi2 = w.kc2 * pm_pos * (1.0 - eta) * Math.exp(-3 * r * r);

      // χ(3): damp if Kerr is high; allow a small positive nudge otherwise
      const kerr = fund.kerrStrength || 0.0;
      const g_chi3 = -w.kc3 * kerr + 0.05 * (1.0 - eta) * (1.0 - pm_pos);

      return {
        gradients: {
          baseIndex: baseScale * g_base,
          dispersion: baseScale * g_disp,
          chi2: baseScale * g_chi2,
          chi3: baseScale * g_chi3,
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
