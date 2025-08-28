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

  // u_lens packing (consistent everywhere in this shader):
  // R: n_o        (ordinary index)
  // G: n_e        (extraordinary index)
  // B: axisAngle  (optic axis orientation, radians)
  // A: dispersion coefficient
  #define LENS_NO(l)      (l.r)
  #define LENS_NE(l)      (l.g)
  #define LENS_AXIS(l)    (l.b)
  #define LENS_DISP(l)    (l.a)

  const float PI = 3.14159265358979323846;

  out vec4 fragColor;

  /*----------------------------------------------------------
    2) Helpers
  ----------------------------------------------------------*/

  vec2 cmul(vec2 a, vec2 b) { return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
  vec2 cconj(vec2 a) { return vec2(a.x, -a.y); }

  float localIavg(vec2 pos, vec2 texel, sampler2D tex){
    float s=0.0;
    for(int j=-1;j<=1;++j){
      for(int i=-1;i<=1;++i){
        vec2 E = texture(tex,(pos+vec2(float(i),float(j)))*texel).xy;
        s += dot(E,E);
      }
    }
    return s/9.0;
  }

  // 2nd-order central diffs (reference + fallback)
  void grad2_E(vec2 pos, vec2 texel, sampler2D tex, out vec2 dEx, out vec2 dEy){
    vec2 E  = texture(tex, pos*texel).xy;
    vec2 Ex = 0.5*(texture(tex,(pos+vec2( 1.0, 0.0))*texel).xy
                - texture(tex,(pos+vec2(-1.0, 0.0))*texel).xy);
    vec2 Ey = 0.5*(texture(tex,(pos+vec2( 0.0, 1.0))*texel).xy
                - texture(tex,(pos+vec2( 0.0,-1.0))*texel).xy);
    dEx = Ex / u_dx;
    dEy = Ey / u_dx;
  }

  // 4th-order central diffs (±1, ±2 stencil)
  void grad4_E(vec2 pos, vec2 texel, sampler2D tex, out vec2 dEx, out vec2 dEy){
    vec2 Ex2  = texture(tex,(pos+vec2( 2.0, 0.0))*texel).xy;
    vec2 Ex1  = texture(tex,(pos+vec2( 1.0, 0.0))*texel).xy;
    vec2 Ex_1 = texture(tex,(pos+vec2(-1.0, 0.0))*texel).xy;
    vec2 Ex_2 = texture(tex,(pos+vec2(-2.0, 0.0))*texel).xy;

    vec2 Ey2  = texture(tex,(pos+vec2( 0.0, 2.0))*texel).xy;
    vec2 Ey1  = texture(tex,(pos+vec2( 0.0, 1.0))*texel).xy;
    vec2 Ey_1 = texture(tex,(pos+vec2( 0.0,-1.0))*texel).xy;
    vec2 Ey_2 = texture(tex,(pos+vec2( 0.0,-2.0))*texel).xy;

    // ( -f_{i+2} + 8 f_{i+1} - 8 f_{i-1} + f_{i-2} ) / (12 Δ)
    dEx = (-Ex2 + 8.0*Ex1 - 8.0*Ex_1 + Ex_2) * (1.0 / (12.0 * u_dx));
    dEy = (-Ey2 + 8.0*Ey1 - 8.0*Ey_1 + Ey_2) * (1.0 / (12.0 * u_dx));
  }

  // Auto-select (edge/low-SNR fallback)
  void grad_auto(vec2 pos, vec2 texel, sampler2D tex, vec2 E0, out vec2 dEx, out vec2 dEy){
    // Near rectangular texture edge? fall back to 2nd order
    bool nearEdge = (pos.x < 2.5) || (pos.y < 2.5)
                 || (pos.x > u_resolution.x - 2.5)
                 || (pos.y > u_resolution.y - 2.5);

    if (!nearEdge){
      vec2 d4x, d4y; grad4_E(pos, texel, tex, d4x, d4y);
      vec2 d2x, d2y; grad2_E(pos, texel, tex, d2x, d2y);

      // intensity-gated blend to keep noise tame when |E| is tiny
      float I  = dot(E0,E0);
      float Ia = localIavg(pos, texel, tex);
      float I0 = max(1e-12, 0.1 * Ia);     // threshold tracks local signal
      float w  = smoothstep(I0, 10.0*I0, I);
      dEx = mix(d2x, d4x, w);
      dEy = mix(d2y, d4y, w);

    } else {
      grad2_E(pos, texel, tex, dEx, dEy);
    }
  }

  // a(n) from lens (same isotropic base you use in varCoeffLaplacian)
  float a_from_lens(vec4 l){
    float n = max(LENS_NO(l), 1.0);
    return (u_c*u_c)/(n*n);
  }

  // relative spread of a over 8 neighbors (small => uniform)
  float localUniformity(vec2 pos, vec2 texel){
    vec4 lL = texture(u_lens,(pos+vec2(-1,0))*texel);
    vec4 lR = texture(u_lens,(pos+vec2( 1,0))*texel);
    vec4 lU = texture(u_lens,(pos+vec2( 0,-1))*texel);
    vec4 lD = texture(u_lens,(pos+vec2( 0, 1))*texel);
    vec4 lUL= texture(u_lens,(pos+vec2(-1,-1))*texel);
    vec4 lUR= texture(u_lens,(pos+vec2( 1,-1))*texel);
    vec4 lDL= texture(u_lens,(pos+vec2(-1, 1))*texel);
    vec4 lDR= texture(u_lens,(pos+vec2( 1, 1))*texel);

    float aL=a_from_lens(lL), aR=a_from_lens(lR), aU=a_from_lens(lU), aD=a_from_lens(lD);
    float aUL=a_from_lens(lUL), aUR=a_from_lens(lUR), aDL=a_from_lens(lDL), aDR=a_from_lens(lDR);

    float amin = min(min(min(min(aL,aR),min(aU,aD)),min(aUL,aUR)),min(aDL,aDR));
    float amax = max(max(max(max(aL,aR),max(aU,aD)),max(aUL,aUR)),max(aDL,aDR));
    float amean = 0.125*(aL+aR+aU+aD+aUL+aUR+aDL+aDR);
    float rel = (amax - amin) / max(amean, 1e-6);
    return rel; // smaller => more uniform
  }

  // constant-coeff 9-point Laplacian (isotropic), scaled by aC
  vec2 laplacian9_const(vec2 pos, vec2 texel, float aC){
    vec2 EC  = texture(u_current, pos*texel).xy;
    vec2 EXp = texture(u_current,(pos+vec2( 1,0))*texel).xy;
    vec2 EXm = texture(u_current,(pos+vec2(-1,0))*texel).xy;
    vec2 EYp = texture(u_current,(pos+vec2(0, 1))*texel).xy;
    vec2 EYm = texture(u_current,(pos+vec2(0,-1))*texel).xy;

    vec2 Epp = texture(u_current,(pos+vec2( 1, 1))*texel).xy;
    vec2 Epm = texture(u_current,(pos+vec2( 1,-1))*texel).xy;
    vec2 Emp = texture(u_current,(pos+vec2(-1, 1))*texel).xy;
    vec2 Emm = texture(u_current,(pos+vec2(-1,-1))*texel).xy;

    vec2 sumCross = EXp + EXm + EYp + EYm;
    vec2 sumDiag  = Epp + Epm + Emp + Emm;

    // (4*cross + diag - 20*center) / (6*dx^2)
    vec2 lap = (4.0*sumCross + sumDiag - 20.0*EC) * (1.0/(6.0*u_dx*u_dx));
    return aC * lap;
  }

  /*
    Approximate partial derivatives in phase with simple wrap.
  */
  vec2 calculatePhaseGradients(vec2 pos, vec2 texel, vec4 center) {
    vec2 E = center.xy;  // use caller-provided center
    vec2 dEx, dEy; grad_auto(pos, texel, u_current, E, dEx, dEy);
    float invI = 1.0 / (dot(E,E) + 1e-12);
    float dphidx = (E.x*dEx.y - E.y*dEx.x) * invI;
    float dphidy = (E.x*dEy.y - E.y*dEy.x) * invI;
    return vec2(dphidx, dphidy);
  }

  /*----------------------------------------------------------
    1) Boundary logic
  ----------------------------------------------------------*/

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
      props.n_local = max(LENS_NO(lensVal), 1.0);  // Ensure not less than vacuum

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

  /*----------------------------------------------------------
    4) Refractive index from local geometry
  ----------------------------------------------------------*/

  float computeAngleDependentBaseIndex(vec2 pos, vec2 texel) {
      // 1) fetch and validate lens data
      vec4 lensVal = texture(u_lens, pos * texel);

      float n_o = max(LENS_NO(lensVal), 1.0);  // ordinary index shouldn't be less than vacuum
      float n_e = max(LENS_NE(lensVal), 1.0);  // extraordinary index shouldn't be less than vacuum
      float axisAngle = LENS_AXIS(lensVal);    // orientation of optic axis at this pixel

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

  vec2 varCoeffLaplacian(vec2 pos, vec2 texel, float aC) {
    // field samples
    vec2 EC = texture(u_current, pos*texel).xy;
    vec2 EL = texture(u_current,(pos+vec2(-1.0,0.0))*texel).xy;
    vec2 ER = texture(u_current,(pos+vec2( 1.0,0.0))*texel).xy;
    vec2 EU = texture(u_current,(pos+vec2( 0.0,-1.0))*texel).xy;
    vec2 ED = texture(u_current,(pos+vec2( 0.0, 1.0))*texel).xy;

    // lens samples (cheap, isotropic base index for a(x))
    vec4 lC=texture(u_lens,pos*texel);
    vec4 lL=texture(u_lens,(pos+vec2(-1,0))*texel);
    vec4 lR=texture(u_lens,(pos+vec2( 1,0))*texel);
    vec4 lU=texture(u_lens,(pos+vec2( 0,-1))*texel);
    vec4 lD=texture(u_lens,(pos+vec2( 0, 1))*texel);

    float nC=max(LENS_NO(lC),1.0);
    float nL=max(LENS_NO(lL),1.0);
    float nR=max(LENS_NO(lR),1.0);
    float nU=max(LENS_NO(lU),1.0);
    float nD=max(LENS_NO(lD),1.0);

    float aL=(u_c*u_c)/(nL*nL), aR=(u_c*u_c)/(nR*nR);
    float aU=(u_c*u_c)/(nU*nU), aD=(u_c*u_c)/(nD*nD);

    // face averages
    float aW=0.5*(aC + aL), aE=0.5*(aC + aR);
    float aN=0.5*(aC + aU), aS=0.5*(aC + aD);

    vec2 flux = aE*(ER-EC) - aW*(EC-EL) + aS*(ED-EC) - aN*(EC-EU);
    return flux / (u_dx*u_dx);
  }

  float getWavelengthDependentIndex(float baseIndex, float dispersion, float wavelength) {
      // Simple dispersion model: n(λ) = n_base + dispersion * (1/λ² - 1/λ_ref²)
      float lambda_ref = u_lambdaFund;
      return baseIndex + dispersion * (1.0/(wavelength*wavelength) - 1.0/(lambda_ref*lambda_ref));
  }

  /*----------------------------------------------------------
    5) Phase mismatch
  ----------------------------------------------------------*/
  vec2 phaseGradientFrom(sampler2D tex, vec2 pos, vec2 texel) {
    vec2 E = texture(tex, pos*texel).xy;
    vec2 dEx, dEy; grad_auto(pos, texel, tex, E, dEx, dEy);
    float invI = 1.0/(dot(E,E)+1e-12);
    return vec2((E.x*dEx.y - E.y*dEx.x)*invI, (E.x*dEy.y - E.y*dEy.x)*invI);
  }

  float n_angle(vec4 lensVal, vec2 khat) {
    float n_o = max(LENS_NO(lensVal),1.0), n_e = max(LENS_NE(lensVal),1.0);
    if (abs(n_e-n_o)<1e-6) return n_o;
    vec2 axis = vec2(cos(LENS_AXIS(lensVal)), sin(LENS_AXIS(lensVal)));
    float c = clamp(dot(khat,axis), -1.0, 1.0);
    float s2 = 1.0 - c*c;
    float denom = n_e*n_e*c*c + n_o*n_o*s2;
    return (n_o*n_e)/sqrt(max(denom,1e-12));
  }

  vec4 calculatePhaseMismatchTerm(vec2 pos, vec2 texel, vec4 /*center*/) {
    vec4 lensVal = texture(u_lens, pos*texel);
    vec2 gF = phaseGradientFrom(u_fundamental, pos, texel);
    vec2 gS = phaseGradientFrom(u_shg,         pos, texel);
    float gFl = length(gF), gSl = length(gS);
    vec2 kF = (gFl > 1e-9) ? (gF / gFl) : vec2(0.0, 1.0);
    vec2 kS = (gSl > 1e-9) ? (gS / gSl) : vec2(0.0, 1.0);

    float nF_base = n_angle(lensVal, kF);
    float nS_base = n_angle(lensVal, kS);

    float nF = getWavelengthDependentIndex(nF_base, LENS_DISP(lensVal), u_lambdaFund);
    float nS = getWavelengthDependentIndex(nS_base, LENS_DISP(lensVal), u_lambdaSHG);

    float kFmag = 2.0*PI*nF / u_lambdaFund;
    float kSmag = 2.0*PI*nS / u_lambdaSHG;
    float deltaK = 2.0*kFmag - kSmag;

    // project onto *fundamental* local propagation direction
    vec2 r = (pos - 0.5 * u_resolution) * u_dx;   // meters-ish in-plane
    float zLocal = dot(kF, r);
    float phase  = deltaK * zLocal + u_phaseRef;

    return vec4(cos(phase), sin(phase), 0.0, 0.0);
  }


  /*----------------------------------------------------------
    6) Dispersion
  ----------------------------------------------------------*/
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
      vec4 gm = texture(u_gainMask, pos * texel);

      // Base index including angle dependence
      props.n_base = computeAngleDependentBaseIndex(pos, texel);
      props.n_base = getWavelengthDependentIndex(
          props.n_base,
          LENS_DISP(lensVal),  // dispersion coefficient
          state.wavelength
      );

      // Nonlinear coefficients — use gainMask G/B as chi2/chi3 spatial weights
      props.chi3_local = u_chi * u_chi * u_chi_ratio * gm.b;
      props.chi2_local = u_chi * u_chi2_ratio * gm.g;
      props.dispCoeff = LENS_DISP(lensVal);

      // Kerr effects
      float kerrSatFactor = 1.0 / (1.0 + state.amp2 / u_kerr_Isat);
      float kerrSelf = props.chi3_local * state.amp2 * kerrSatFactor;
      float kerrCross = u_crossKerrCoupling * props.chi3_local * otherAmp2 * kerrSatFactor;

      // n_eff only used to form the nonlinear phase kick later:
      props.n_eff = max(props.n_base + kerrSelf + kerrCross, 1e-6);

      // Linear operator uses base only:
      float c_lin = u_c / props.n_base;
      props.c2_eff = c_lin * c_lin;  // <-- no min(..., 1.0) clamp (see CFL below)

      return props;
  }

  // Compute conversion terms between fundamental and SHG
  ConversionTerms computeConversionTerms(
    vec2 pos, vec2 texel,
    FieldState fund, FieldState shg,
    MaterialProps props, bool isFundamental
  ){
    ConversionTerms T;
    vec4 mismatch = calculatePhaseMismatchTerm(pos, texel, vec4(fund.field,0,0));
    vec2 ph = mismatch.xy;                 // e^{+iΔkz}
    vec2 phConj = vec2(ph.x, -ph.y);       // e^{-iΔkz}

    float sat = 1.0 / (1.0 + dot(fund.field,fund.field)/u_shg_Isat);
    float k = u_conversionCoupling * props.chi2_local;

    if (isFundamental){
      // d²E1/dt² … + k * E2 * E1* * e^{-iΔkz}
      vec2 term = cmul(shg.field, cmul(cconj(fund.field), phConj));
      T.upConversion   = vec2(0.0);               // fundamental doesn't get +E1²
      T.downConversion = -k * term * sat;         // depletion by back-conversion
    } else {
      // d²E2/dt² … + k * E1² * e^{+iΔkz}
      vec2 term = cmul( cmul(fund.field, fund.field), ph );
      T.upConversion   =  k * term * sat;         // build-up of SHG
      T.downConversion = vec2(0.0);
    }

    // A simple, stable “phase-match quality” scalar for your debug channel
    vec2 e1sq = cmul(fund.field, fund.field);
    float L = length(e1sq);
    vec2 e1sq_hat = (L > 1e-12) ? (e1sq / L) : vec2(1.0, 0.0);
    T.phaseMatch = clamp(dot(ph, e1sq_hat), -1.0, 1.0);

    return T;
  }


  // Compute the interior field evolution
  vec4 computeInteriorField(vec2 pos, vec2 texel, bool isFundamental, float damp) {
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

      // Strang-split: use only the post-linear half-kick.
      // Do not apply a pre-linear rotation here so the Laplacian and neighbor samples
      // read the same (unmodified) u_current. χ(2) conversion terms were already
      // computed above from the un-kicked primary, so we leave them as-is.

      // 5) Evolution (linear propagation + χ(2) conversions + gain)
      // CFL safety clamp (limits c_eff * dt / dx to ~1/2 in 2D finite difference)
      float c_eff = sqrt(props.c2_eff);
      float s_cfl = min(1.0, (1.0 / sqrt(2.0)) * (u_dx / (c_eff * u_dt)));

      // ======= Laplacian blending =======
      const bool  USE_LAP9      = true;  // set false to disable
      const float LAP9_WEIGHT   = 0.6;   // 0..1 max blend when uniform
      const float UNIF_EPS      = 0.02;  // smaller = stricter "uniform a" test

      vec2 flux = varCoeffLaplacian(pos, texel, props.c2_eff); // conservative 5-pt (variable a)
      if (USE_LAP9){
        float rel = localUniformity(pos, texel);
        // 1 when uniform, 0 when not
        float wUniform = 1.0 - smoothstep(UNIF_EPS, 2.0*UNIF_EPS, rel);
        vec2 flux9 = laplacian9_const(pos, texel, props.c2_eff);
        float w = clamp(LAP9_WEIGHT * wUniform, 0.0, 1.0);
        w = clamp(w, 0.0, 0.95);  // avoid fully replacing conservative flux

        // slightly safer CFL as 9-pt share grows (it’s "stiffer" at Nyquist)
        s_cfl *= (1.0 - 0.2*w);

        flux = mix(flux, flux9, w);
      }

      vec2 linTerm = s_cfl * (u_dt * u_dt) * flux;


      vec2 newField = (2.0 * primary.field - primary.oldField)
                    + linTerm
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

      // Phase-preserving sponge (PML-lite): matched damping on time derivative
      // Compute once per pixel using the boundary "damp" (1 inside → 0 outside)
      float sigma = (1.0 - damp) * (1.0 - damp) * u_damping; // quadratic ramp
      vec2 dEdt   = (texture(u_current,pos*texel).xy - texture(u_previous,pos*texel).xy) / u_dt;
      // Add to update (with sign that damps energy)
      newField -= u_dt * sigma * dEdt;

      if (damp < 1.0) newField *= u_damping;  // keep interior energy

      // Post linear-step: recompute local nonlinear index using updated intensity,
      // then apply the second half of the nonlinear phase kick (Strang splitting).
      FieldState updatedPrimary = FieldState(newField, primary.oldField, dot(newField, newField), primary.wavelength);
      MaterialProps props2 = computeMaterialProperties(pos, texel, updatedPrimary, other.amp2);

      float phiNL2 = (2.0 * PI / updatedPrimary.wavelength) * (props2.n_eff - props2.n_base) * u_dt * 0.5;
      vec2 rot2 = vec2(cos(phiNL2), sin(phiNL2));
      newField = vec2(
        newField.x * rot2.x - newField.y * rot2.y,
        newField.x * rot2.y + newField.y * rot2.x
      );

      return vec4(newField, conv.phaseMatch, props2.n_eff - props2.n_base);
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
      vec4 fieldUpdate = computeInteriorField(pos, texel, u_updateTarget == 0, damp);

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
    rgbaData[i * 4 + 1] = 0.001; // chi2 weight
    rgbaData[i * 4 + 2] = 0.001; // chi3 weight
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
