// Shared vertex shader for all programs
export const vertexShaderSource = `#version 300 es
in vec2 position;
out vec2 uv;
void main() {
    uv = position * 0.5 + 0.5;
    gl_Position = vec4(position, 0.0, 1.0);
}`;

// Main simulation shader
export const getSimulationShaderSource = (config) => {
  return `#version 300 es
  precision highp float;
  precision highp sampler2D;

  /*----------------------------------------------------------
    0) Uniforms
  ----------------------------------------------------------*/
  uniform sampler2D u_current;
  uniform sampler2D u_previous;
  uniform sampler2D u_lens;
  uniform sampler2D u_fundamental;
  uniform sampler2D u_shg;

  // NEW: Additional uniforms for saturable gain
  uniform sampler2D u_gainMask;
  uniform float u_gain0;
  uniform float u_gainSat;
  uniform float u_linearLoss;

  // OPTIONAL: A mask for “neurons” that amplifies Kerr effects locally
  // uniform sampler2D u_neuronMask;

  // Grid, time, boundary
  uniform float u_dt;
  uniform float u_dx;
  uniform float u_damping;
  uniform float u_c;
  uniform float u_chi;
  uniform float u_chi_ratio;
  uniform float u_shg_Isat;
  uniform float u_kerr_Isat;
  uniform vec2  u_resolution;
  uniform float u_boundaryR0;
  uniform float u_boundaryAlpha;
  uniform float u_boundaryM;
  uniform int   u_updateTarget;
  uniform int   u_frameCount;
  uniform int   u_pulseInterval;

  // Phase matching & additional
  uniform float u_lambdaFund;
  uniform float u_lambdaSHG;
  uniform float u_temperature;
  uniform float u_phaseRef;
  uniform float u_crossKerrCoupling;
  uniform float u_conversionCoupling;

  const float PI = 3.14159265358979323846;

  /*----------------------------------------------------------
    1) Boundary logic
  ----------------------------------------------------------*/
  bool insideBoundary(vec2 pos) {
      vec2 center  = 0.5 * u_resolution;
      vec2 rel     = pos - center;
      float r      = length(rel);
      float theta  = atan(rel.y, rel.x);
      float radius = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));
      return (r <= radius);
  }

  /*
     Optionally, you may wish to implement a “soft” boundary or an absorbing
     layer. The below function returns reflectivity=1.0 for perfect reflection,
     but you could do something like:

         float reflectivity = 1.0 - smoothstep(R0, R1, r);

     if you want a tapered absorbing zone between [R0, R1].
  */
  /* Example uniform parameters:
     - u_useSoftBoundary: bool - enable the smooth radial fade
     - u_absorbRadiusInner, u_absorbRadiusOuter: float - radial region in which reflectivity goes from near 1.0 to near 0.0
     - u_usePositionDependentReflect: bool - enable an angular or radial reflectivity variation
     - u_reflectBase: float - baseline reflectivity
     - u_reflectAmplitude: float - amplitude of reflectivity modulation
     - u_reflectFreq: float - frequency (e.g. how many cycles around the boundary)
     - u_enableHardClipAtOuterBoundary: bool - if “true,” anything beyond outer boundary is forced to 0 reflectivity
  */

  uniform bool  u_useSoftBoundary;
  uniform float u_absorbRadiusInner;
  uniform float u_absorbRadiusOuter;

  uniform bool  u_usePositionDependentReflect;
  uniform float u_reflectBase;
  uniform float u_reflectAmplitude;
  uniform float u_reflectFreq;

  uniform bool  u_enableHardClipAtOuterBoundary;

  float getBoundaryReflectivity(vec2 pos)
  {
      // Define default values for the variables
      float u_useSoftBoundary = 1.0;       // 1.0 = true, 0.0 = false
      float u_absorbRadiusInner = 100.0;  // Inner radius for absorbing layer
      float u_absorbRadiusOuter = 200.0;  // Outer radius for absorbing layer
      float u_usePositionDependentReflect = 1.0; // 1.0 = true, 0.0 = false
      float u_reflectBase = 0.5;          // Base reflectivity
      float u_reflectAmplitude = 0.5;     // Amplitude of position-dependent reflectivity
      float u_reflectFreq = 6.28;         // Frequency for position-dependent reflectivity (e.g., radians)
      float u_enableHardClipAtOuterBoundary = 1.0; // 1.0 = true, 0.0 = false

      // 1) Basic radial distance from center
      vec2 center = 0.5 * u_resolution;
      float r     = length(pos - center);

      // 2) Start with a default reflectivity of 1.0
      float reflectVal = 1.0;

      // 3) If we want a soft boundary (absorbing layer),
      //    fade out reflectivity from r=in to r=out:
      if (u_useSoftBoundary > 0.5) {
          // The factor alpha = 0 in [0, absorbRadiusInner],
          // alpha=1 in [absorbRadiusOuter, ∞].
          // reflectVal = 1 at alpha=0, 0 at alpha=1.
          float alpha = clamp(
              (r - u_absorbRadiusInner) / (u_absorbRadiusOuter - u_absorbRadiusInner),
              0.0, 1.0
          );
          float fade = 1.0 - alpha;
          reflectVal *= fade;
      }

      // 4) If we want position-dependent reflectivity (angular or radial):
      //    For example, a sinusoidal function of the angle:
      if (u_usePositionDependentReflect > 0.5) {
          float angle = atan(pos.y - center.y, pos.x - center.x);
          // E.g. reflect pattern = reflectBase + reflectAmplitude * cos(reflectFreq * angle)
          // Adjust as you see fit, could also be radial-based or a custom function:
          float pattern = u_reflectBase
                        + u_reflectAmplitude * cos(u_reflectFreq * angle);
          // ensure we don’t go negative if amplitude is large
          pattern = max(0.0, pattern);
          reflectVal *= pattern;
      }

      // 5) If we do not want *any* reflection beyond outer boundary
      //    (e.g. a “hard clip”):
      if (u_enableHardClipAtOuterBoundary > 0.5
          && u_useSoftBoundary > 0.5
          && (r > u_absorbRadiusOuter)) {
          reflectVal = 0.0;
      }

      // 6) Return final reflectivity
      return reflectVal;
  }


  /*----------------------------------------------------------
    2) Phase derivative helpers with optional simple unwrap
  ----------------------------------------------------------*/
  float safeAtan2(float y, float x) {
      // Some hardware or older ES profiles can show quirks
      // on atan(y,x).  Usually fine in modern GLSL though.
      return atan(y, x);
  }

  /*
     A naive attempt at local phase unwrapping. We measure
     the difference, then keep it in (-π, π):
  */
  float wrappedPhaseDelta(float phase1, float phase2) {
      float diff = phase1 - phase2;
      // shift diff to (-π, π)
      diff = mod(diff + PI, 2.0 * PI) - PI;
      return diff;
  }

  /* Approximate partial derivatives in phase with simple wrap. */
  vec2 calculatePhaseGradients(vec2 pos, vec2 texel, vec4 center) {
      float centerPhase = safeAtan2(center.y, center.x);

      vec4 right = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 left  = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 up    = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 down  = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);

      float phaseR = safeAtan2(right.y, right.x);
      float phaseL = safeAtan2(left.y,  left.x);
      float phaseU = safeAtan2(up.y,    up.x);
      float phaseD = safeAtan2(down.y,  down.x);

      float dphase_dx = (wrappedPhaseDelta(phaseR, centerPhase)
                       - wrappedPhaseDelta(phaseL, centerPhase)) / (2.0 * u_dx);
      float dphase_dy = (wrappedPhaseDelta(phaseD, centerPhase)
                       - wrappedPhaseDelta(phaseU, centerPhase)) / (2.0 * u_dx);

      return vec2(dphase_dx, dphase_dy);
  }

  /*----------------------------------------------------------
    3) Laplacian stencils
  ----------------------------------------------------------*/
  vec2 ninePointLaplacian(vec2 pos, vec2 texel, vec4 center) {
      vec4 left   = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 right  = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 up     = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 down   = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);

      // Diagonals
      vec4 upleft    = texture(u_current, (pos + vec2(-1.0, -1.0)) * texel);
      vec4 upright   = texture(u_current, (pos + vec2( 1.0, -1.0)) * texel);
      vec4 downleft  = texture(u_current, (pos + vec2(-1.0,  1.0)) * texel);
      vec4 downright = texture(u_current, (pos + vec2( 1.0,  1.0)) * texel);

      // Weighted sum
      // The factor 2/3 for direct neighbors, 1/6 for diagonals, -4 for center
      // is one typical approach; you may want 1/(u_dx*u_dx) after summing.
      vec2 sumDirect    = (left + right + up + down).xy * (2.0 / 3.0);
      vec2 sumDiagonals = (upleft + upright + downleft + downright).xy * (1.0 / 6.0);
      vec2 centerTerm   = center.xy * -4.0;

      // Typically you'd multiply the entire sum by (1.0 / (u_dx*u_dx)) if
      // you want a physically consistent wave equation. Do that below.
      vec2 lap = sumDirect + sumDiagonals + centerTerm;
      return lap / (u_dx * u_dx);
  }

  vec2 fivePointLaplacian(vec2 pos, vec2 texel, vec4 center) {
      vec4 left  = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 right = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 up    = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 down  = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);

      vec2 lap = (left.xy + right.xy + up.xy + down.xy - 4.0 * center.xy);
      return lap / (u_dx * u_dx);
  }

  /*----------------------------------------------------------
    4) Temperature-dependent Sellmeier (same basic logic)
  ----------------------------------------------------------*/
  float calculateRefractiveIndex(float wavelength, float temperature, bool extraordinary) {
      float lam_um = wavelength * 1e6;

      // (Fictitious example coefficients)
      float A_o = 4.9048;
      float B_o = 0.11768;
      float C_o = 0.04750;
      float D_o = 1.5e-5;
      float T0  = 25.0;

      float A_e = 4.5820;
      float B_e = 0.09920;
      float C_e = 0.04438;
      float D_e = 2.0e-5;

      if (!extraordinary) {
          float n_o_sq = A_o + (B_o / (lam_um * lam_um - C_o)) - D_o*(temperature - T0);
          return sqrt(n_o_sq);
      } else {
          float n_e_sq = A_e + (B_e / (lam_um * lam_um - C_e)) - D_e*(temperature - T0);
          return sqrt(n_e_sq);
      }
  }

  /*----------------------------------------------------------
    5) Angle-dependent index
  ----------------------------------------------------------*/
  float calculateAngleDependentIndex(float theta, float ne, float no) {
      float ne2  = ne * ne;
      float no2  = no * no;
      float c2   = cos(theta); c2 *= c2;
      float s2   = sin(theta); s2 *= s2;
      return sqrt((ne2 * no2) / (ne2 * c2 + no2 * s2));
  }

  /*----------------------------------------------------------
    6) Phase mismatch term
  ----------------------------------------------------------*/
  vec4 calculatePhaseMismatchTerm(vec2 pos, vec2 texel, vec4 center) {
      float amp2 = dot(center.xy, center.xy);

      // For demonstration, use ordinary index for both:
      float nFund = calculateRefractiveIndex(u_lambdaFund, u_temperature, false);
      float nSHG  = calculateRefractiveIndex(u_lambdaSHG, u_temperature, false);

      float kFund = 2.0 * PI * nFund / u_lambdaFund;
      float kSHG  = 2.0 * PI * nSHG  / u_lambdaSHG;

      float deltaK = 2.0 * kFund - kSHG;

      // Incorporate a mild position dependence along y
      float phase    = deltaK * (pos.y * u_dx) + u_phaseRef;
      vec2 mismatch  = vec2(cos(phase), sin(phase));

      return vec4(mismatch, 0.0, 0.0);
  }

  /*----------------------------------------------------------
    7) Main update
  ----------------------------------------------------------*/
  out vec4 fragColor;

  void main() {
      vec2 pos   = gl_FragCoord.xy;
      vec2 texel = 1.0 / u_resolution;

      // Boundary condition check
      if (!insideBoundary(pos)) {
          vec4 current = texture(u_current, pos * texel);
          float boundaryReflectivity = getBoundaryReflectivity(pos);

          // Scale amplitude by reflectivity (and optionally invert normal comp).
          vec2 reflectedField = current.xy * boundaryReflectivity;

          // Could do e.g. “absorbing boundary” if we want:
          // reflectedField *= 0.95;  // damp a bit each frame
          // Or more advanced absorbing layers.

          fragColor = vec4(reflectedField, 0.0, 1.0);
          return;
      }

      // Fetch fields at center
      vec4 center   = texture(u_current, pos * texel);
      vec4 oldField = texture(u_previous, pos * texel);
      vec4 lens     = texture(u_lens, pos * texel);

      // Select which Laplacian to use
      vec2 laplacianRaw = ${
        config.use9PointStencil ? "nine" : "five"
      }PointLaplacian(pos, texel, center);

      // “Lens” or local wavefront factor:
      // multiply (laplacianRaw.x + i laplacianRaw.y) by (lens.x + i lens.y)
      vec2 lensedLaplacian = vec2(
          laplacianRaw.x * lens.x - laplacianRaw.y * lens.y,
          laplacianRaw.x * lens.y + laplacianRaw.y * lens.x
      );

      // === Saturable gain calculation ===
      float amp2        = dot(center.xy, center.xy);
      float localGain   = u_gain0 / (1.0 + amp2 / u_gainSat) - u_linearLoss;
      float gainMaskVal = texture(u_gainMask, pos * texel).r;
      localGain        *= gainMaskVal;  // zero or scaled outside pump region

      // === OPTIONAL “neuron” factor for localized Kerr amplification ===
      // float neuronMaskVal = texture(u_neuronMask, pos * texel).r;
      // float neuronFactor   = 1.0 + 2.0 * neuronMaskVal;
      float neuronFactor = 1.0;  // just identity if not used

      /*------------------------------------------------------------------
        Fundamental vs. SHG update branches
      ------------------------------------------------------------------*/
      if (u_updateTarget == 0) {
          // -------------- Fundamental --------------
          vec4 shg       = texture(u_shg, pos * texel);
          float shg_amp2 = dot(shg.xy, shg.xy);

          // Basic Kerr
          float kerrSat  = 1.0 / (1.0 + amp2 / u_kerr_Isat);
          float n2_eff   = u_chi * u_chi * u_chi_ratio;
          float kerrSelf = n2_eff * amp2 * kerrSat;
          float kerrCross= u_crossKerrCoupling * n2_eff * shg_amp2 * kerrSat;
          float totalNon = (kerrSelf + kerrCross) * neuronFactor; // “nonlinearity”

          // Effective speed
          float c2_local = u_c * u_c * (1.0 + totalNon);

          // SHG backreaction
          vec4 mismatchTerm   = calculatePhaseMismatchTerm(pos, texel, center);
          float shgAmp2_satur = shg_amp2;
          vec2 backReaction   = -u_conversionCoupling * u_chi * shgAmp2_satur
                                * mismatchTerm.xy;

          // Wave equation update (leapfrog)
          vec2 newField = (2.0 * center.xy - oldField.xy)
                        + c2_local * u_dt * u_dt * lensedLaplacian
                        + u_dt * u_dt * backReaction
                        + u_dt * u_dt * localGain * center.xy;
          newField *= u_damping;

          // Add pulsed injection
          if (u_pulseInterval > 0 && (u_frameCount % u_pulseInterval == 0)) {
              float w0 = 4.0;
              vec2 centerGrid = 0.5 * u_resolution;
              vec2 rel = pos - centerGrid;
              float r2 = dot(rel, rel);

              // amplitude ~ Gaussian
              float amplitude = 0.01 * exp(-r2/(w0*w0));
              float phase     = 0.0; // could modulate phase if desired

              vec2 pulseField = amplitude * vec2(cos(phase), sin(phase));
              // lens it
              vec2 lensedPulse = vec2(
                  pulseField.x * lens.x - pulseField.y * lens.y,
                  pulseField.x * lens.y + pulseField.y * lens.x
              );
              newField += lensedPulse;
          }

          // Phase gradient measure
          vec2 phaseGrad = calculatePhaseGradients(pos, texel, center);

          // Final color
          fragColor = vec4(newField,
                           length(phaseGrad),
                           totalNon);

      } else {
          // -------------- SHG --------------
          float shg_amp2 = dot(center.xy, center.xy);

          float kerrSat     = 1.0 / (1.0 + shg_amp2 / u_kerr_Isat);
          float n2_eff      = u_chi * u_chi * u_chi_ratio;
          float kerrSelf    = n2_eff * shg_amp2 * kerrSat;

          // Cross-Kerr from fundamental
          vec4 fund       = texture(u_fundamental, pos * texel);
          float fund_amp2 = dot(fund.xy, fund.xy);
          float kerrCross = u_crossKerrCoupling * n2_eff * fund_amp2 * kerrSat;

          float totalNon  = (kerrSelf + kerrCross) /* * neuronMaskVal, if used */;
          float c2_local  = u_c * u_c * (1.0 + totalNon);

          // SHG saturation vs. fundamental amplitude
          float shgSaturFactor = 1.0 / (1.0 + fund_amp2 / u_shg_Isat);

          // Mismatch
          vec4 mismatchTerm = calculatePhaseMismatchTerm(pos, texel, center);

          // Forward SHG gen
          vec2 sourceTerm = u_conversionCoupling * u_chi * fund_amp2
                          * shgSaturFactor * mismatchTerm.xy;

          // Back-conversion
          vec2 backConversionTerm = -u_conversionCoupling * u_chi
                                    * shg_amp2 * mismatchTerm.xy;

          // Wave eq
          vec2 newField = (2.0 * center.xy - oldField.xy)
                        + c2_local * u_dt * u_dt * lensedLaplacian
                        + u_dt * u_dt * (sourceTerm + backConversionTerm)
                        + u_dt * u_dt * localGain * center.xy;
          newField *= u_damping;

          // Example local measure of phase match quality
          float localPhaseMatch = dot(mismatchTerm.xy, sourceTerm)
                                / (length(mismatchTerm.xy) * length(sourceTerm) + 1e-10);

          // Example ratio or dimensionless measure of SHG efficiency
          float chiRatio = (u_chi * fund_amp2 * shgSaturFactor)
                         / (n2_eff * shg_amp2 * kerrSat + 1e-10);

          fragColor = vec4(
              newField,
              localPhaseMatch,
              totalNon
          );
      }
  }
`;
};

export const displayShaderSource = `#version 300 es
precision highp float;
precision highp sampler2D;

uniform sampler2D u_field;
uniform int u_displayMode;  // 0 = wave, 1 = shg, 2 = lens
uniform float u_lensDisplayMin;
uniform float u_lensDisplayMax;
uniform float u_lensRadius;
uniform vec2 u_resolution;

in vec2 uv;
out vec4 fragColor;

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main() {
  if (u_displayMode == 2) {  // Lens display
      vec2 center = vec2(0.5, 0.5);
      vec2 fromCenter = uv - center;

      // this works because we always scale the lens to the full extent of the canvas
      if (length(fromCenter) > .5) {
          fragColor = vec4(0.0, 0.0, 0.0, 1.0);
          return;
      }

      // Scale our sampling coordinates to only look within the lens radius
      vec2 scaledUV = center + fromCenter * (u_lensRadius / (u_resolution.x * 0.5));

      // Check if we're outside the canvas
      if (scaledUV.x < 0.0 || scaledUV.x > 1.0 ||
          scaledUV.y < 0.0 || scaledUV.y > 1.0) {
          fragColor = vec4(0.0, 0.0, 0.0, 1.0);
          return;
      }

      vec4 field = texture(u_field, scaledUV);
      float amp = length(field.xy);
      float phase = atan(field.y, field.x);

      float normAmp = clamp(
          (amp - u_lensDisplayMin) / (u_lensDisplayMax - u_lensDisplayMin),
          0.0, 1.0
      );

      float gamma = 0.4;
      float val = pow(normAmp, gamma);
      float hue = (phase + 3.14159) / 6.28318;
      float sat = 1.0;

      fragColor = vec4(hsv2rgb(vec3(hue, sat, val)), 1.0);
  } else {  // Wave or SHG display
        vec4 field = texture(u_field, uv);
        float amp = length(field.xy);
        float phase = atan(field.y, field.x);

        float hue = (phase + 3.14159) / 6.28318;
        float sat = 1.0;
        float val = amp / (1.0 + amp);

        fragColor = vec4(hsv2rgb(vec3(hue, sat, val)), 1.0);
    }
}`;
