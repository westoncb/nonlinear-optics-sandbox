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

  // Saturable gain
  uniform sampler2D u_gainMask;
  uniform float u_gain0;
  uniform float u_gainSat;
  uniform float u_linearLoss;

  // Grid, time, boundary
  uniform float u_dt;
  uniform float u_dx;
  uniform float u_damping;
  uniform float u_c;            // baseline speed in vacuum or reference medium

  uniform float u_chi;          // global scaling for chi(2) and chi(3)
  uniform float u_chi_ratio;
  uniform float u_chi2_ratio;
  uniform float u_shg_Isat;
  uniform float u_kerr_Isat;
  uniform float u_crossKerrCoupling;
  uniform float u_conversionCoupling;

  uniform vec2  u_resolution;
  uniform float u_boundaryR0;
  uniform float u_boundaryAlpha;
  uniform float u_boundaryM;
  uniform float u_boundaryReflectivity;
  uniform int   u_updateTarget;
  uniform int   u_frameCount;
  uniform int   u_pulseInterval;

  // Phase matching
  uniform float u_lambdaFund;
  uniform float u_lambdaSHG;
  uniform float u_temperature;
  uniform float u_phaseRef;

  const float PI = 3.14159265358979323846;

  out vec4 fragColor;

  /*----------------------------------------------------------
    2) Helpers
  ----------------------------------------------------------*/
  float wrappedPhaseDelta(float phase1, float phase2) {
      float diff = phase1 - phase2;
      // shift diff to (-π, π)
      diff = mod(diff + PI, 2.0 * PI) - PI;
      return diff;
  }

  /*
    Approximate partial derivatives in phase with simple wrap.
  */
  vec2 calculatePhaseGradients(vec2 pos, vec2 texel, vec4 center) {
      // Need 2 points on each side
      vec4 far_left   = texture(u_current, (pos + vec2(-2.0,  0.0)) * texel);
      vec4 left       = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 right      = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 far_right  = texture(u_current, (pos + vec2( 2.0,  0.0)) * texel);

      vec4 far_down   = texture(u_current, (pos + vec2( 0.0,  2.0)) * texel);
      vec4 down       = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);
      vec4 up         = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 far_up     = texture(u_current, (pos + vec2( 0.0, -2.0)) * texel);

      // 4th order coefficients: (-1/12, 8/12, -8/12, 1/12)
      float dphase_dx =
          (-wrappedPhaseDelta(atan(far_right.y, far_right.x), atan(far_left.y, far_left.x)) / 12.0 +
            wrappedPhaseDelta(atan(right.y, right.x), atan(left.y, left.x)) * (8.0/12.0)) / (u_dx);

      float dphase_dy =
          (-wrappedPhaseDelta(atan(far_down.y, far_down.x), atan(far_up.y, far_up.x)) / 12.0 +
            wrappedPhaseDelta(atan(down.y, down.x), atan(up.y, up.x)) * (8.0/12.0)) / (u_dx);

      return vec2(dphase_dx, dphase_dy);
  }

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

  vec2 boundaryNormal(vec2 pos) {
      // Get position relative to center
      vec2 center = 0.5 * u_resolution;
      vec2 rel = pos - center;
      float r = length(rel);

      // Calculate theta and relevant trig terms
      float theta = atan(rel.y, rel.x);
      float cosTheta = cos(theta);
      float sinTheta = sin(theta);
      float cosMTheta = cos(u_boundaryM * theta);
      float sinMTheta = sin(u_boundaryM * theta);

      // For a polar curve r = R(θ), the normal vector is:
      // N = (R(θ) - R'(θ))ȓ + R(θ)θ̂
      // where ȓ = (cos θ, sin θ) and θ̂ = (-sin θ, cos θ)

      // Base radius with modulation
      float R = u_boundaryR0 * (1.0 + u_boundaryAlpha * cosMTheta);

      // Derivative of R with respect to theta
      float Rprime = -u_boundaryR0 * u_boundaryAlpha * u_boundaryM * sinMTheta;

      // Construct normal vector components
      vec2 rHat = vec2(cosTheta, sinTheta);
      vec2 thetaHat = vec2(-sinTheta, cosTheta);

      vec2 normal = (R - Rprime) * rHat + R * thetaHat;

      // Normalize and return (pointing inward by convention)
      return -normalize(normal);
  }

  // A helper to do a simpler TE Fresnel reflectivity:
  float fresnelReflectivity_TE(float n_in, float n_out, float theta_i) {
      // Snell's law: n_in * sin(theta_i) = n_out * sin(theta_t)
      float sin_t = (n_in / n_out) * sin(theta_i);
      // clamp to avoid domain errors
      sin_t = clamp(sin_t, -1.0, 1.0);

      float theta_t = asin(sin_t);

      // TE reflectivity: R_TE = |(n_in cos(theta_i) - n_out cos(theta_t)) /
      //                          (n_in cos(theta_i) + n_out cos(theta_t))|^2
      float cos_i = cos(theta_i);
      float cos_t = cos(theta_t);

      float r_num = (n_in * cos_i) - (n_out * cos_t);
      float r_den = (n_in * cos_i) + (n_out * cos_t);
      float r     = r_num / r_den;
      return r * r;
  }

  float fresnelReflectivity_TM(float n_in, float n_out, float theta_i) {
      float sin_t = (n_in / n_out) * sin(theta_i);
      sin_t = clamp(sin_t, -1.0, 1.0);
      float theta_t = asin(sin_t);

      float cos_i = cos(theta_i);
      float cos_t = cos(theta_t);

      float r_num = (n_out * cos_i) - (n_in * cos_t);
      float r_den = (n_out * cos_i) + (n_in * cos_t);
      float r     = r_num / r_den;
      return r * r;
  }

  float computeAngleDependentReflectivity(vec2 pos, vec2 texel) {
      vec2 centerVal = texture(u_current, pos * texel).xy;

      // 1) If amplitude is too small, skip
      if (dot(centerVal, centerVal) < 1e-9) return 0.0;

      // 2) Get local wave-vector direction, etc.
      vec2 kvec = calculatePhaseGradients(pos, texel, vec4(centerVal, 0.0, 0.0));
      float cos_in = dot(normalize(kvec), boundaryNormal(pos)); // define a boundaryNormal(pos)
      float theta_i = acos(clamp(cos_in, -1.0, 1.0));

      // 3) Local index inside vs. outside
      float n_local = max(texture(u_lens, pos * texel).x, 1e-6);

      float R_TE = fresnelReflectivity_TE(n_local, 1.0, theta_i);
      float R_TM = fresnelReflectivity_TM(n_local, 1.0, theta_i);

      // 4) Weighted or averaged reflection
      float R_avg = 0.5 * (R_TE + R_TM);

      return R_avg;
  }

  /*----------------------------------------------------------
    3) 9-point Laplacian (fourth-order accurate)
  ----------------------------------------------------------*/
  vec2 ninePointLaplacian(vec2 pos, vec2 texel, vec4 center) {
      // Direct neighbors
      vec4 left   = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 right  = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 up     = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 down   = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);

      // Diagonals
      vec4 ul  = texture(u_current, (pos + vec2(-1.0, -1.0)) * texel);
      vec4 ur  = texture(u_current, (pos + vec2( 1.0, -1.0)) * texel);
      vec4 dl  = texture(u_current, (pos + vec2(-1.0,  1.0)) * texel);
      vec4 dr  = texture(u_current, (pos + vec2( 1.0,  1.0)) * texel);

      // Standard 9-point 4th-order weights:
      // ∇²f ≈ (1 / (6 h^2)) * [4*(N,S,E,W) + (NE,NW,SE,SW) - 20*fC ]
      // sum of weights = 0 for constant field.
      vec2 lap = (4.0 * (left.xy + right.xy + up.xy + down.xy)
                + (ul.xy + ur.xy + dl.xy + dr.xy)
                - 20.0 * center.xy) / (6.0 * (u_dx * u_dx));

      return lap;
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
    4) Refractive index from local geometry
  ----------------------------------------------------------*/

  float computeAngleDependentBaseIndex(vec2 pos, vec2 texel) {
      // 1) fetch and validate lens data
      vec4 lensVal = texture(u_lens, pos * texel);

      float n_o = max(lensVal.x, 1.0);  // ordinary index shouldn't be less than vacuum
      float n_e = max(lensVal.y, 1.0);  // extraordinary index shouldn't be less than vacuum
      float axisAngle = lensVal.z;      // orientation of optic axis at this pixel

      // Early exit if indices are effectively equal (isotropic case)
      if (abs(n_e - n_o) < 1e-6) {
          return n_o;
      }

      // 2) compute local wave-vector direction using the provided gradient function
      vec4 cVal = texture(u_current, pos * texel);
      vec2 gradPhi = calculatePhaseGradients(pos, texel, cVal);
      float gradLen = length(gradPhi);

      // if gradient is too small, return ordinary index
      const float MIN_GRAD_LENGTH = 1e-9;
      if (gradLen < MIN_GRAD_LENGTH) {
          return n_o;
      }

      // 3) compute direction of wave vector and angle with optic axis
      vec2 k_hat = gradPhi / gradLen;
      vec2 axisDir = vec2(cos(axisAngle), sin(axisAngle));
      float cosAngle = dot(k_hat, axisDir);
      cosAngle = clamp(cosAngle, -1.0, 1.0);

      // Optimization: avoid acos/sin calls by working with cosine directly
      float cosT = cosAngle;
      float sinT_sq = 1.0 - cosT * cosT;  // sin²(θ) = 1 - cos²(θ)

      // 4) uniaxial formula
      // n(θ) = (n_o n_e)/sqrt(n_e² cos²(θ) + n_o² sin²(θ))
      float denom = n_e * n_e * cosT * cosT + n_o * n_o * sinT_sq;
      denom = max(denom, 1e-12);  // Prevent division by zero

      return (n_o * n_e) / sqrt(denom);
  }

  float getWavelengthDependentIndex(float baseIndex, float dispersion, float wavelength) {
      // Simple dispersion model: n(λ) = n_base + dispersion * (1/λ² - 1/λ_ref²)
      float lambda_ref = u_lambdaFund;
      return baseIndex + dispersion * (1.0/(wavelength*wavelength) - 1.0/(lambda_ref*lambda_ref));
  }

  /*----------------------------------------------------------
    5) Phase mismatch
  ----------------------------------------------------------*/
  vec4 calculatePhaseMismatchTerm(vec2 pos, vec2 texel, vec4 center) {
      // 1. First get angle-dependent base indices
      float n_base_fund = computeAngleDependentBaseIndex(pos, texel);
      vec4 shgVal = texture(u_shg, pos * texel);
      float n_base_shg = computeAngleDependentBaseIndex(pos, texel);

      // 2. Get dispersion coefficient from lens texture
      vec4 lensVal = texture(u_lens, pos * texel);
      float dispersionCoeff = lensVal.y;

      // 3. Combine angle and wavelength dependence
      float n_fund = getWavelengthDependentIndex(n_base_fund, dispersionCoeff, u_lambdaFund);
      float n_shg = getWavelengthDependentIndex(n_base_shg, dispersionCoeff, u_lambdaSHG);

      // 4. Calculate k-vectors and phase mismatch as before
      float kFund = 2.0 * PI * n_fund / u_lambdaFund;
      float kSHG = 2.0 * PI * n_shg / u_lambdaSHG;
      float deltaK = 2.0 * kFund - kSHG;

      float phase = deltaK * (pos.y * u_dx) + u_phaseRef;
      return vec4(cos(phase), sin(phase), 0.0, 0.0);
  }

  /*----------------------------------------------------------
    6) Local wave speed, dispersion
  ----------------------------------------------------------*/
  // We'll interpret lens.x as the local base index n_base.
  float computeLocalSpeed(float n_local) {
      // Avoid division by zero
      n_local = max(n_local, 1e-6);
      // wave speed in this region
      return u_c / n_local;
  }

  // A simplistic approach to "dispersion" as a second time derivative adjustment
  vec2 applyDispersion(
      vec2 centerField, vec2 oldField,
      float dispersionCoeff, float dt
  ) {
      // Very ad-hoc: new += dispCoeff*(E(t)-2E(t-dt)) * dt^2
      vec2 secondTimeDerivativeApprox = centerField - 2.0 * oldField;
      return dispersionCoeff * secondTimeDerivativeApprox * (dt * dt);
  }

  // =======================================
  // Main
  // =======================================
  void main() {
      vec2 pos   = gl_FragCoord.xy;
      vec2 texel = 1.0 / u_resolution;

      // -----------------------------------
      // 1) Boundary logic with angle-dependent reflection
      // -----------------------------------
      if (!insideBoundary(pos)) {
          // If outside boundary, reflect wave with angle-dependent Fresnel logic
          vec4 current = texture(u_current, pos * texel);

          // For a real mirrored cavity boundary:
          float boundaryReflectivity = computeAngleDependentReflectivity(pos, texel);

          // Multiply field by reflection coefficient
          vec2 reflectedField = current.xy * boundaryReflectivity;

          // Typically we won't compute localPhaseMatch or totalNon outside
          // the boundary, so set them to zero (or some sentinel).
          fragColor = vec4(reflectedField, 0.0, 0.0);
          return;
      }

      // -----------------------------------
      // 2) Fetch fields
      // -----------------------------------
      vec4 center   = texture(u_current, pos * texel);
      vec4 oldField = texture(u_previous, pos * texel);
      vec4 lensVal  = texture(u_lens, pos * texel);

      // For convenience
      float amp2 = dot(center.xy, center.xy);

      // -----------------------------------
      // 3) Angle-dependent base index
      // -----------------------------------
      float n_baseAngle = computeAngleDependentBaseIndex(pos, texel);

      // -----------------------------------
      // 4) Kerr shift
      //     n_kerr = chi3_local * amp^2 / (1 + amp^2/I_sat)
      // -----------------------------------
      float chi3_local = u_chi * u_chi * u_chi_ratio * lensVal.w;  // lensVal.w used for local chi(3)
      float kerrSatFactor  = 1.0 / (1.0 + amp2 / u_kerr_Isat);
      float n_kerr         = chi3_local * amp2 * kerrSatFactor;


      // -----------------------------------
      // 5) Combine into total index => local speed
      // -----------------------------------
      // We'll add cross-Kerr separately. Here is just the base + self-Kerr:
      float n_eff = n_baseAngle + n_kerr;
      n_eff = max(n_eff, 1e-6);  // clamp to avoid div by zero

      float c_local = u_c / n_eff;
      float c2_local = c_local * c_local;

      float dispersionCoeff = lensVal.y;

      // -----------------------------------
      // 6) Saturable gain
      // localGain = gain0/(1+amp2/gainSat) - linearLoss
      // -----------------------------------
      float localGain = u_gain0 / (1.0 + amp2 / u_gainSat) - u_linearLoss;
      float gainMaskVal = texture(u_gainMask, pos * texel).r;
      localGain *= gainMaskVal;

      // We'll store placeholders for localPhaseMatch & totalNon for final output.
      float localPhaseMatch = 0.0;
      float totalNon        = 0.0;

      // -----------------------------------
      // 7) Switch: Fundamental or SHG update
      // -----------------------------------
      if (u_updateTarget == 0) {
          // === Fundamental update ===
          vec4 shg = texture(u_shg, pos * texel);
          float shg_amp2 = dot(shg.xy, shg.xy);
          float fund_amp2 = amp2;

          // Base refractive index
          float n_baseAngle = computeAngleDependentBaseIndex(pos, texel);
          float n_fund = getWavelengthDependentIndex(n_baseAngle, dispersionCoeff, u_lambdaFund);

          // Kerr effects
          float chi3_local = u_chi * u_chi * u_chi_ratio * lensVal.w;
          float fundKerrSatFactor = 1.0 / (1.0 + fund_amp2 / u_kerr_Isat);
          float kerrSelf = chi3_local * fund_amp2 * fundKerrSatFactor;
          float kerrCross = u_crossKerrCoupling * chi3_local * shg_amp2 * fundKerrSatFactor;
          totalNon = kerrSelf + kerrCross;

          // Total effective index and speed
          float n_eff_f = n_fund + totalNon;
          n_eff_f = max(n_eff_f, 1e-6);
          float c_eff_f = u_c / n_eff_f;
          float c2_eff_f = c_eff_f * c_eff_f;
          c2_eff_f = min(c2_eff_f, 1.);

          // Phase mismatch and conversion
          vec4 mismatchTerm = calculatePhaseMismatchTerm(pos, texel, center);
          float chi2_local = u_chi * u_chi2_ratio * lensVal.z;  // lensVal.z used for local chi(2)

          // Fundamental branch - should lose twice the energy it contributes to SHG
          vec2 upConversion = -2.0 * u_conversionCoupling * chi2_local * fund_amp2
                            * (1.0 / (1.0 + fund_amp2 / u_shg_Isat))
                            * mismatchTerm.xy;
          vec2 downConversion = u_conversionCoupling * chi2_local * shg_amp2
                             * mismatchTerm.xy;

          // Phase matching measure using both up/down conversion
          float upMatch = dot(mismatchTerm.xy, upConversion)
                       / (length(mismatchTerm.xy) * length(upConversion) + 1e-10);
          float downMatch = dot(mismatchTerm.xy, downConversion)
                         / (length(mismatchTerm.xy) * length(downConversion) + 1e-10);
          localPhaseMatch = sign(upMatch) * sqrt(abs(upMatch * downMatch));

          // Standard leapfrog
          vec2 lap = ${config.use9PointStencil ? "nine" : "five"}PointLaplacian(pos, texel, center);
          vec2 newField = (2.0 * center.xy - oldField.xy)
                        + c2_eff_f * (u_dt * u_dt) * lap
                        + u_dt * u_dt * downConversion
                        + u_dt * u_dt * localGain * center.xy;

          // dispersion
          newField += applyDispersion(center.xy, oldField.xy, dispersionCoeff, u_dt);

          // Optional pulsed injection
          if (u_pulseInterval > 0 && (u_frameCount % u_pulseInterval == 0)) {
              float w0 = 4.0;
              vec2 centerGrid = 0.5 * u_resolution;
              vec2 rel = pos - centerGrid;
              float r2 = dot(rel, rel);
              float amplitude = 0.01 * exp(-r2/(w0*w0));
              float phase = 0.0;
              vec2 pulseField = amplitude * vec2(cos(phase), sin(phase));
              newField += pulseField;
          }

          newField *= u_damping;
          fragColor = vec4(newField, localPhaseMatch, totalNon);

      } else {
          // === SHG update ===
          vec4 fund = texture(u_fundamental, pos * texel);
          float fund_amp2 = dot(fund.xy, fund.xy);
          float shg_amp2 = amp2;

          // Base refractive index combining angle
          float n_baseAngle = computeAngleDependentBaseIndex(pos, texel);
          float n_shg = getWavelengthDependentIndex(n_baseAngle, dispersionCoeff, u_lambdaSHG);

          // Kerr effects
          float chi3_local = u_chi * u_chi * u_chi_ratio * lensVal.w;
          float shgKerrSatFactor = 1.0 / (1.0 + shg_amp2 / u_kerr_Isat);
          float kerrSelf = chi3_local * shg_amp2 * shgKerrSatFactor;
          float kerrCross = u_crossKerrCoupling * chi3_local * fund_amp2 * shgKerrSatFactor;
          totalNon = kerrSelf + kerrCross;

          // Total effective index and speed
          float n_eff_shg = n_shg + totalNon;
          n_eff_shg = max(n_eff_shg, 1e-6);
          float c_eff_shg = u_c / n_eff_shg;
          float c2_eff_shg = c_eff_shg * c_eff_shg;
          c2_eff_shg = min(c2_eff_shg, 1.);

          // Phase mismatch and conversion
          vec4 mismatchTerm = calculatePhaseMismatchTerm(pos, texel, center);
          float chi2_local = u_chi * lensVal.z;  // lensVal.z used for local chi(2)

          // SHG branch - gets one photon for every two fundamental photons
          vec2 upConversion = u_conversionCoupling * chi2_local * fund_amp2
                            * (1.0 / (1.0 + fund_amp2 / u_shg_Isat))
                            * mismatchTerm.xy;
          vec2 downConversion = -u_conversionCoupling * chi2_local * shg_amp2
                             * mismatchTerm.xy;

          // Phase matching measure using both up/down conversion
          float upMatch = dot(mismatchTerm.xy, upConversion)
                       / (length(mismatchTerm.xy) * length(upConversion) + 1e-10);
          float downMatch = dot(mismatchTerm.xy, downConversion)
                         / (length(mismatchTerm.xy) * length(downConversion) + 1e-10);
          localPhaseMatch = sign(upMatch) * sqrt(abs(upMatch * downMatch));

          // Standard leapfrog
          vec2 lap = ${config.use9PointStencil ? "nine" : "five"}PointLaplacian(pos, texel, center);
          vec2 newField = (2.0 * center.xy - oldField.xy)
                        + c2_eff_shg * (u_dt * u_dt) * lap
                        + u_dt * u_dt * (upConversion + downConversion)
                        + u_dt * u_dt * localGain * center.xy;

          newField += applyDispersion(center.xy, oldField.xy, dispersionCoeff, u_dt);
          newField *= u_damping;

          fragColor = vec4(newField, localPhaseMatch, totalNon);
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
uniform vec2 u_baseIndexRange;
uniform vec2 u_chi2Range;
uniform vec2 u_chi3Range;
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

   if (length(fromCenter) > .5) {
       fragColor = vec4(0.0, 0.0, 0.0, 1.0);
       return;
   }

   vec2 scaledUV = center + fromCenter * (u_lensRadius / (u_resolution.x * 0.5));
   vec4 lens = texture(u_field, scaledUV);

   // Sample neighboring points
   float eps = 1.0 / u_resolution.x;
   vec4 dx = texture(u_field, scaledUV + vec2(eps, 0.0)) - lens;
   vec4 dy = texture(u_field, scaledUV + vec2(0.0, eps)) - lens;

   // Calculate relative spatial variations
   float baseGrad = length(vec2(dx.x, dy.x)) / max(abs(lens.x - 1.0), 1e-6);
   float dispGrad = length(vec2(dx.y, dy.y)) / max(lens.y, 1e-6);
   float chi2Grad = length(vec2(dx.z, dy.z)) / max(lens.z, 1e-6);
   float chi3Grad = length(vec2(dx.w, dy.w)) / max(lens.w, 1e-6);

   // Calculate color components
   float hueBias = 0.4;
   float hueScale = 3.0;   // Increased spread for more color variation
   float hue = hueBias + (atan(chi3Grad, chi2Grad) / (2.0 * 3.14159)) * hueScale;

   // Combine nonlinear and dispersion gradients for saturation
   float nonlinearStrength = sqrt(chi2Grad * chi2Grad + chi3Grad * chi3Grad);
   float combinedGrad = mix(nonlinearStrength, dispGrad, 0.3);  // 30% dispersion weight
   float sat = sqrt(combinedGrad);

   // Base value on index gradient but modulate with dispersion
   float val = 0.2 + 0.8 * (1.0 - exp(-baseGrad * 100.0));
   val = mix(val, val * (1.0 - 0.2 * dispGrad), 0.2);  // 20% dispersion influence

   // Enhance local contrast
   sat = pow(sat, 0.3);
   val = pow(val, 0.7);

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
