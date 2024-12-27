export const getUpdateStrategies = (config) => {
  return {
    phaseMatch: (fund, shg, count, zone, sector, lensParams) => {
      // 1. Calculate basic field measures
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);
      const fundPhase = Math.atan2(fund.imag, fund.real);

      // 2. Consider spatial position in lens
      const r = zone / config.fresnelZones; // normalized radius
      const theta = (2 * Math.PI * sector) / config.numSectors;

      // 3. Phase matching quality measures
      const conversionEfficiency = shgAmp / (fundAmp + 1e-6);
      const phaseQuality = shg.phaseMatch;

      // 4. Position-dependent objectives
      // - Inner zone: prioritize χ(2) for conversion
      // - Middle zone: balance χ(2) and χ(3)
      // - Outer zone: focus on phase matching via base index
      const zoneWeight = {
        chi2: Math.exp(-2 * r * r), // strongest in center
        chi3: Math.sin(Math.PI * r), // peaks in middle
        base: r, // increases toward edge
      };

      // 5. Calculate composite loss
      const loss =
        1.0 * fundAmp + // reduce fundamental
        0.5 * fund.kerrStrength - // limit kerr effects
        2.0 * conversionEfficiency - // encourage conversion
        1.0 * phaseQuality + // improve phase matching
        0.1 * Math.sin(4 * theta); // subtle angular modulation

      // 6. Position-aware updates
      const baseScale = 0.015;
      const update = {
        baseIndex: -baseScale * loss * zoneWeight.base,
        dispersion: -0.2 * baseScale * loss * (1 - r), // stronger in center
        chi2: -2.0 * baseScale * loss * zoneWeight.chi2,
        chi3: -baseScale * loss * zoneWeight.chi3,
      };

      return { update, loss };
    },

    phaseMatchAsymmetric: (fund, shg, count, zone, sector, lensParams) => {
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
        (zonePreference > 0 ? 2.0 : 0.5) * fundAmp + // varying fundamental penalty
        0.5 * fund.kerrStrength -
        (zonePreference < 0 ? 2.0 : 0.5) * shgAmp - // varying SHG reward
        1.0 * shg.phaseMatch +
        0.1 * Math.sin(4 * theta);

      const baseScale = 0.015;
      const update = {
        baseIndex: -baseScale * loss * zoneWeight.base,
        dispersion: -0.2 * baseScale * loss * (1 - r),
        chi2: -2.0 * baseScale * loss * zoneWeight.chi2,
        chi3: -baseScale * loss * zoneWeight.chi3,
      };

      return { update, loss };
    },

    debug: (fund, shg, count, zone, sector, lensParams) => {
      // 1) Compute amplitude
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

      // 2) Weighted sum for "loss"
      const wKerr = 0.2;
      const wPM = 0.5;
      const loss =
        fundAmp + wKerr * fund.kerrStrength - (shgAmp + wPM * shg.phaseMatch);

      // 3) If the absolute loss is super small, let's force a nonzero update
      //    for debugging.  (Remove if you only want "legit" updates.)
      let forcedLoss = loss;
      if (Math.abs(forcedLoss) < 1e-6) {
        forcedLoss = 1.0; // Force a sign or value for debugging
        console.log(" -- Force: Setting forcedLoss=1.0 for debug");
      }

      // 4) Build the update
      const scale = 0.02;
      const update = {
        baseIndex: -scale * forcedLoss,
        dispersion: 0, // not adjusting dispersion in this example
        chi2: -1.5 * scale * forcedLoss,
        chi3: -scale * forcedLoss,
      };

      return { update, loss };
    },

    // A lens strategy that aims to maximize SHG amplitude,
    // reduce fundamental amplitude, and ensure good phase matching.
    // It also penalizes extremely high Kerr strength to avoid blowups.
    SHGTest: (fund, shg, count, zone, sector, lensParams) => {
      /*
          fund: { real, imag, phaseGradMag, kerrStrength }
          shg:  { real, imag, phaseMatch,   chiRatio }
          lensParams: { baseIndex, dispersion, chi2, chi3 }
        */

      // 1) Basic amplitude measures
      const fundAmp = Math.sqrt(fund.real * fund.real + fund.imag * fund.imag);
      const shgAmp = Math.sqrt(shg.real * shg.real + shg.imag * shg.imag);

      // 2) Weighted multi-objective loss
      //    We want to:
      //      - Minimize fundamental amplitude
      //      - Maximize SHG amplitude => negative sign in the formula
      //      - Maximize phaseMatch => negative sign
      //      - Minimize kerrStrength => positive sign
      //
      //    So we define something like:
      //
      //      loss = w_fund * fundAmp
      //           + w_kerr * fund.kerrStrength
      //           - w_shg * shgAmp
      //           - w_pm * shg.phaseMatch
      //
      //    Because we do gradient descent on this "loss,"
      //    the negative terms are effectively "rewards."
      const w_fund = 1.0; // penalty weight for leftover fundamental
      const w_shg = 2.0; // reward weight for big SHG amplitude
      const w_kerr = 0.5; // penalty weight for large kerr strength
      const w_pm = 1.0; // reward weight for good phase match

      const loss =
        w_fund * fundAmp +
        w_kerr * fund.kerrStrength -
        w_shg * shgAmp -
        w_pm * shg.phaseMatch;

      // 3) Very naive "gradient" from the loss
      //    In a real physics-based approach, each lens param's partial derivative
      //    might differ. But here we simply do:  updateParam = -α * loss
      //
      //    Some heuristics you might adopt:
      //      - big negative update for chi2 if we want to push more SHG,
      //        or big negative update for baseIndex if we suspect that
      //        adjusting refractive index helps phase matching, etc.
      //      - or treat dispersion differently if you want to manage GVD.
      //
      //    Below, we assume that adjusting chi2 is the main driver for SHG,
      //    so we weight it more strongly in the gradient.
      const alpha = 0.02; // base scale
      const update = {
        baseIndex: -alpha * loss,
        dispersion: -0.5 * alpha * loss, // maybe we tweak dispersion less
        chi2: -2.0 * alpha * loss, // put stronger emphasis on chi2
        chi3: -alpha * loss,
      };

      return { update, loss };
    },
  };
};
