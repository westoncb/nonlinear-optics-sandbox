import { gaussianRandom, mulberry32, polar } from "./util";

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

  uniform float u_beamWidth;
  uniform float u_subsequentPulseAmplitude;

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
  uniform float u_boundaryTransitionWidth;
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

  float boundaryDampingFactor(vec2 pos) {
      // 1) Get position relative to center in polar coordinates
      vec2 center = 0.5 * u_resolution;
      vec2 rel = pos - center;
      float r = length(rel);
      float theta = atan(rel.y, rel.x);

      // 2) Calculate modulated boundary radius at this angle
      float boundaryRadius = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));

      // 3) If inside modulated boundary, no damping
      if (r <= boundaryRadius) {
          return 1.0;
      }

      // 4) Calculate distance past boundary
      float distPast = r - boundaryRadius;
      if (distPast >= u_boundaryTransitionWidth) {
          return 0.0;  // Far outside => complete damping
      }

      // 5) Smooth transition from 1.0 to 0.0
      // Optional: Could use smoother transition like cosine instead of linear
      float t = distPast / u_boundaryTransitionWidth;
      return 1.0 - t;  // Linear falloff
      // return 0.5 * (1.0 + cos(PI * t));  // Cosine falloff (smoother)
  }

  float fresnelReflectivity_TE(float n_in, float n_out, float theta_i) {
      // Ensure valid indices and angle
      n_in = max(n_in, 1e-6);
      n_out = max(n_out, 1e-6);
      theta_i = clamp(theta_i, 0.0, PI/2.0);  // Restrict to [0, π/2]

      // Check for total internal reflection
      float crit_angle = (n_out < n_in) ? asin(n_out/n_in) : PI;
      if (theta_i > crit_angle) {
          return 1.0;
      }

      // Fresnel calculation with safeguards
      float sin_t = (n_in / n_out) * sin(theta_i);
      sin_t = clamp(sin_t, -1.0, 1.0);

      float theta_t = asin(sin_t);
      float cos_i = cos(theta_i);
      float cos_t = cos(theta_t);

      // Add small epsilon to prevent division by zero
      const float eps = 1e-6;
      float r_num = (n_in * cos_i) - (n_out * cos_t);
      float r_den = (n_in * cos_i) + (n_out * cos_t) + eps;
      float r = r_num / r_den;
      return r * r;  // Ensure output is valid reflectivity
  }

  float fresnelReflectivity_TM(float n_in, float n_out, float theta_i) {
      // Ensure valid indices and angle
      n_in = max(n_in, 1e-6);
      n_out = max(n_out, 1e-6);
      theta_i = clamp(theta_i, 0.0, PI/2.0);  // Restrict to [0, π/2]

      // Check for total internal reflection
      float crit_angle = (n_out < n_in) ? asin(n_out/n_in) : PI;
      if (theta_i > crit_angle) {
          return 1.0;
      }

      // Fresnel calculation with safeguards
      float sin_t = (n_in / n_out) * sin(theta_i);
      sin_t = clamp(sin_t, -1.0, 1.0);
      float theta_t = asin(sin_t);

      float cos_i = cos(theta_i);
      float cos_t = cos(theta_t);

      // Add small epsilon to prevent division by zero
      const float eps = 1e-6;
      float r_num = (n_out * cos_i) - (n_in * cos_t);
      float r_den = (n_out * cos_i) + (n_in * cos_t) + eps;
      float r = r_num / r_den;
      return r * r;  // Ensure output is valid reflectivity
  }

  struct BoundaryProps {
      vec2  normal;      // Boundary normal vector (points inward)
      float theta_i;     // Angle of incidence
      float n_local;     // Local refractive index
      vec2  kvec;        // Local wave vector direction
      float amp2;        // Field intensity
      bool  validField;  // Whether field is strong enough for meaningful calculations
  };

  BoundaryProps computeBoundaryProperties(vec2 pos, vec2 texel) {
      BoundaryProps props;

      // 1) Get local field
      vec4 centerVal = texture(u_current, pos * texel);
      props.amp2 = dot(centerVal.xy, centerVal.xy);

      // Check if field is strong enough for meaningful calculations
      const float MIN_FIELD_STRENGTH = 1e-9;
      props.validField = (props.amp2 >= MIN_FIELD_STRENGTH);

      if (!props.validField) {
          // Set defaults for weak field case
          props.normal = vec2(0.0);
          props.theta_i = 0.0;
          props.n_local = 1.0;
          props.kvec = vec2(0.0);
          return props;
      }

      // 2) Get local wave-vector direction
      props.kvec = calculatePhaseGradients(pos, texel, centerVal);

      // 3) Get boundary normal at this position
      props.normal = boundaryNormal(pos);

      // 4) Calculate angle of incidence
      float cos_in = dot(normalize(props.kvec), props.normal);
      props.theta_i = acos(clamp(cos_in, -1.0, 1.0));

      // 5) Get local refractive index
      vec4 lensVal = texture(u_lens, pos * texel);
      props.n_local = max(lensVal.x, 1.0);  // Ensure not less than vacuum

      return props;
  }

  struct FresnelCoeffs {
      float R_TE;
      float R_TM;
      float R_avg;
  };

  FresnelCoeffs computeFresnelCoefficients(BoundaryProps props) {
      FresnelCoeffs coeffs;

      if (!props.validField) {
          coeffs.R_TE = 0.0;
          coeffs.R_TM = 0.0;
          coeffs.R_avg = 0.0;
          return coeffs;
      }

      // Compute both polarization states
      coeffs.R_TE = fresnelReflectivity_TE(props.n_local, 1.0, props.theta_i);
      coeffs.R_TM = fresnelReflectivity_TM(props.n_local, 1.0, props.theta_i);
      coeffs.R_avg = 0.5 * (coeffs.R_TE + coeffs.R_TM);

      return coeffs;
  }

  vec4 computeBoundaryField(vec2 pos, vec2 texel) {
      // 1) Get boundary properties
      BoundaryProps props = computeBoundaryProperties(pos, texel);

      if (!props.validField) {
          return vec4(0.0);
      }

      // 2) Compute Fresnel coefficients with validation
      FresnelCoeffs coeffs = computeFresnelCoefficients(props);

      // 3) Get current field
      vec4 current = texture(u_current, pos * texel);

      // 4) Compute reflected and transmitted parts with validation
      vec2 reflectPart = current.xy * coeffs.R_avg;
      vec2 transmitPart = current.xy * (1.0 - coeffs.R_avg);

      float damp = boundaryDampingFactor(pos);  // Smooth spatial transition

      // 5) Apply damping to transmitted part
      vec2 outField = reflectPart + transmitPart * damp;

      return vec4(outField, 0.0, 0.0);
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

  // Struct to hold common field properties
  struct FieldState {
      vec2 field;           // Current field value
      vec2 oldField;        // Previous field value
      float amp2;           // Field intensity
      float wavelength;     // Field wavelength
  };

  // Struct to hold material properties
  struct MaterialProps {
      float n_base;         // Base refractive index
      float n_eff;          // Effective index with nonlinear contributions
      float c2_eff;         // Effective speed squared
      float chi2_local;     // Local chi(2) coefficient
      float chi3_local;     // Local chi(3) coefficient
      float dispCoeff;      // Dispersion coefficient
  };

  // Struct to hold conversion terms
  struct ConversionTerms {
      vec2 upConversion;
      vec2 downConversion;
      float phaseMatch;
  };

  // Compute material properties for a given field
  MaterialProps computeMaterialProperties(vec2 pos, vec2 texel, FieldState state, float otherAmp2) {
      MaterialProps props;
      vec4 lensVal = texture(u_lens, pos * texel);

      // Base index including angle dependence
      props.n_base = computeAngleDependentBaseIndex(pos, texel);
      props.n_base = getWavelengthDependentIndex(
          props.n_base,
          lensVal.y,  // dispersion coefficient
          state.wavelength
      );

      // Nonlinear coefficients
      props.chi3_local = u_chi * u_chi * u_chi_ratio * lensVal.w;
      props.chi2_local = u_chi * u_chi2_ratio * lensVal.z;
      props.dispCoeff = lensVal.y;

      // Kerr effects
      float kerrSatFactor = 1.0 / (1.0 + state.amp2 / u_kerr_Isat);
      float kerrSelf = props.chi3_local * state.amp2 * kerrSatFactor;
      float kerrCross = u_crossKerrCoupling * props.chi3_local * otherAmp2 * kerrSatFactor;

      // Effective index and speed
      props.n_eff = props.n_base + kerrSelf + kerrCross;
      props.n_eff = max(props.n_eff, 1e-6);
      float c_eff = u_c / props.n_eff;
      props.c2_eff = min(c_eff * c_eff, 1.0);

      return props;
  }

  // Compute conversion terms between fundamental and SHG
  ConversionTerms computeConversionTerms(
      vec2 pos, vec2 texel,
      FieldState fund, FieldState shg,
      MaterialProps props, bool isFundamental
  ) {
      ConversionTerms terms;

      vec4 mismatchTerm = calculatePhaseMismatchTerm(pos, texel, vec4(fund.field, 0.0, 0.0));
      float saturationFactor = 1.0 / (1.0 + fund.amp2 / u_shg_Isat);

      if (isFundamental) {
          // Fundamental loses two photons in up-conversion
          terms.upConversion = -2.0 * u_conversionCoupling * props.chi2_local
                           * fund.amp2 * saturationFactor * mismatchTerm.xy;
          terms.downConversion = u_conversionCoupling * props.chi2_local
                             * shg.amp2 * mismatchTerm.xy;
      } else {
          // SHG gains one photon for every two fundamental photons
          terms.upConversion = u_conversionCoupling * props.chi2_local
                           * fund.amp2 * saturationFactor * mismatchTerm.xy;
          terms.downConversion = -u_conversionCoupling * props.chi2_local
                             * shg.amp2 * mismatchTerm.xy;
      }

      // Phase matching measure
      float upMatch = dot(mismatchTerm.xy, terms.upConversion)
                   / (length(mismatchTerm.xy) * length(terms.upConversion) + 1e-10);
      float downMatch = dot(mismatchTerm.xy, terms.downConversion)
                     / (length(mismatchTerm.xy) * length(terms.downConversion) + 1e-10);
      terms.phaseMatch = sign(upMatch) * sqrt(abs(upMatch * downMatch));

      return terms;
  }

  // Compute the interior field evolution
  vec4 computeInteriorField(vec2 pos, vec2 texel, bool isFundamental) {
      // 1) Get field states
      vec4 center = texture(u_current, pos * texel);
      vec4 oldField = texture(u_previous, pos * texel);

      FieldState primary, other;
      vec4 otherTex;

      if (isFundamental) {
          primary = FieldState(
              center.xy,
              oldField.xy,
              dot(center.xy, center.xy),
              u_lambdaFund
          );
          otherTex = texture(u_shg, pos * texel);
          other = FieldState(
              otherTex.xy,
              vec2(0.0),
              dot(otherTex.xy, otherTex.xy),
              u_lambdaSHG
          );
      } else {
          primary = FieldState(
              center.xy,
              oldField.xy,
              dot(center.xy, center.xy),
              u_lambdaSHG
          );
          otherTex = texture(u_fundamental, pos * texel);
          other = FieldState(
              otherTex.xy,
              vec2(0.0),
              dot(otherTex.xy, otherTex.xy),
              u_lambdaFund
          );
      }

      // 2) Compute material properties
      MaterialProps props = computeMaterialProperties(pos, texel, primary, other.amp2);

      // 3) Compute conversion terms
      ConversionTerms conv;
      if (isFundamental) {
          conv = computeConversionTerms(pos, texel, primary, other, props, true);
      } else {
          conv = computeConversionTerms(pos, texel, other, primary, props, false);
      }

      // 4) Compute gain
      float localGain = u_gain0 / (1.0 + primary.amp2 / u_gainSat) - u_linearLoss;
      localGain *= texture(u_gainMask, pos * texel).r;

      // 5) Evolution
      vec2 lap = ${config.use9PointStencil ? "nine" : "five"}PointLaplacian(pos, texel, center);
      vec2 newField = (2.0 * primary.field - primary.oldField)
                    + props.c2_eff * (u_dt * u_dt) * lap
                    + u_dt * u_dt * (conv.upConversion + conv.downConversion)
                    + u_dt * u_dt * localGain * primary.field;

      // 6) Add dispersion and damping
      newField += applyDispersion(primary.field, primary.oldField, props.dispCoeff, u_dt);

      // 7) Optional pulse injection for fundamental
      if (isFundamental && u_pulseInterval > 0 && (u_frameCount % u_pulseInterval == 0)) {
          float w0 = u_beamWidth;
          vec2 centerGrid = 0.5 * u_resolution;
          vec2 rel = pos - centerGrid;
          float r2 = dot(rel, rel);
          float amplitude = u_subsequentPulseAmplitude * exp(-r2/(w0*w0));
          newField += amplitude * vec2(1.0, 0.0);  // Real pulse
      }

      newField *= u_damping;

      return vec4(newField, conv.phaseMatch, props.n_eff - props.n_base);
  }

  // =======================================
  // Main
  // =======================================
  void main() {
      vec2 pos = gl_FragCoord.xy;
      vec2 texel = 1.0 / u_resolution;
      float damp = boundaryDampingFactor(pos);

      if (damp <= 0.0) {
          // Fully outside - no field
          fragColor = vec4(0.0);
          return;
      }

      // Always compute the wave equation update
      vec4 fieldUpdate = computeInteriorField(pos, texel, u_updateTarget == 0);

      // If we're in the transition region, apply partial reflection
      if (damp < 1.0) {
          BoundaryProps props = computeBoundaryProperties(pos, texel);

          if (props.validField) {
              FresnelCoeffs coeffs = computeFresnelCoefficients(props);

              // Get current field
              vec4 current = texture(u_current, pos * texel);

              // Scale reflection by how "outside" we are
              float outsideness = 1.0 - damp;

              // Reflected part (increases with outsideness)
              vec2 reflectedPart = current.xy * coeffs.R_avg * outsideness;

              // Transmitted part (from PDE update, decreases with outsideness)
              vec2 transmittedPart = fieldUpdate.xy * damp;

              // Combine them
              fieldUpdate = vec4(
                  reflectedPart + transmittedPart,
                  fieldUpdate.z,  // Keep phase match term
                  fieldUpdate.w   // Keep nonlinearity term
              );
          }
      }

      fragColor = fieldUpdate;
  }
`;
};

export const getInitialFieldState = (config) => {
  const fundamentalData = new Float32Array(
    config.gridSize * config.gridSize * 4,
  );
  const shgData = new Float32Array(config.gridSize * config.gridSize * 4);

  // Initialize PRNG for thermal noise
  const seed = Math.floor(Math.random() * 2 ** 32);
  const prng = mulberry32(seed);

  // Physical constants for thermal noise
  const kB = 1.380649e-23;
  const T = 295.15;
  const scalingFactor = config.initialNoiseScale;
  const sigma = Math.sqrt(kB * T) * scalingFactor;

  // Gaussian beam parameters
  const w0 = config.beamWidth;
  const z = 0;
  const lambda = 1;
  const k = (2 * Math.PI) / lambda;
  const zR = (Math.PI * w0 * w0) / lambda;

  const centerX = Math.floor(config.gridSize / 2);
  const centerY = Math.floor(config.gridSize / 2);

  // Initialize fundamental field

  // Gaussian beam with phase terms for fundamental
  for (let y = 0; y < config.gridSize; y++) {
    for (let x = 0; x < config.gridSize; x++) {
      const idx = (y * config.gridSize + x) * 4;
      const dx = x - centerX;
      const dy = y - centerY;
      const r2 = dx * dx + dy * dy;

      const w = w0 * Math.sqrt(1 + (z / zR) ** 2);
      const R = z === 0 ? Infinity : z * (1 + (zR / z) ** 2);
      const gouyPhase = Math.atan(z / zR);

      const phase =
        -k * z -
        (R === Infinity ? 0 : (k * r2) / (2 * R)) +
        gouyPhase +
        config.initialPulsePhaseShift;

      const amplitude =
        config.initialPulseAmplitude * (w0 / w) * Math.exp(-r2 / (w * w));

      // Add thermal noise to fundamental
      const noiseAmp = sigma * gaussianRandom(prng);
      const noisePhase = 2.0 * Math.PI * prng();

      fundamentalData[idx] =
        amplitude * Math.cos(phase) + noiseAmp * Math.cos(noisePhase);
      fundamentalData[idx + 1] =
        amplitude * Math.sin(phase) + noiseAmp * Math.sin(noisePhase);
    }
  }

  // Initialize SHG field with very small noise only
  const shgNoiseFactor = 0.01; // Much smaller noise for SHG
  for (let y = 0; y < config.gridSize; y++) {
    for (let x = 0; x < config.gridSize; x++) {
      const idx = (y * config.gridSize + x) * 4;
      const noiseAmp = sigma * shgNoiseFactor * gaussianRandom(prng);
      const noisePhase = 2.0 * Math.PI * prng();
      shgData[idx] = noiseAmp * Math.cos(noisePhase);
      shgData[idx + 1] = noiseAmp * Math.sin(noisePhase);
      // Initialize extra channels to zero
      shgData[idx + 2] = 0;
      shgData[idx + 3] = 0;
      fundamentalData[idx + 2] = 0;
      fundamentalData[idx + 3] = 0;
    }
  }

  // Initialize gain mask
  const gainMaskData = createGainMaskData(config);
  const rgbaData = new Float32Array(config.gridSize * config.gridSize * 4);
  for (let i = 0; i < gainMaskData.length; i++) {
    rgbaData[i * 4] = gainMaskData[i]; // R channel
    rgbaData[i * 4 + 1] = 0; // G channel
    rgbaData[i * 4 + 2] = 0; // B channel
    rgbaData[i * 4 + 3] = 1; // A channel
  }

  return { shgData, fundamentalData, gainMaskData: rgbaData };
};

function createGainMaskData(config) {
  const gridSize = config.gridSize;
  const gainMaskData = new Float32Array(gridSize * gridSize);
  const centerX = Math.floor(gridSize / 2);
  const centerY = Math.floor(gridSize / 2);

  for (let y = 0; y < gridSize; y++) {
    for (let x = 0; x < gridSize; x++) {
      const coords = polar.getZoneAndSector(
        x,
        y,
        centerX,
        centerY,
        config.lensRadius,
        config.fresnelZones,
        config.numSectors,
      );

      if (coords) {
        // Radial gradient (stronger at center)
        const normalizedRadius = coords.zone / config.fresnelZones;
        gainMaskData[y * gridSize + x] = 1.0 - normalizedRadius;
      } else {
        // Outside the lens disk - no gain
        gainMaskData[y * gridSize + x] = 0.0;
      }
    }
  }

  return gainMaskData;
}
