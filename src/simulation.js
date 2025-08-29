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

  // ==== Dispersion (Lorentz ADE) defaults: {W0, Γ, α0, κ}
  const vec4 DISP_FUND = vec4(1.0, 0.5, 0.05, 0.05);
  const vec4 DISP_SHG  = vec4(1.0, 0.5, 0.05, 0.05);

  // ==== Walk-off controls
  const float WALK_TAN_CAP  = 0.35; // ~20°
  const float WALK_MAX_PIX  = 0.45;  // per-step lateral hop cap
  const float WALK_BLEND    = 0.5;  // how strongly we trust advected field

  // Max allowed nonlinear phase advance per *full* step (radians)
  const float PHI_NL_MAX = 1.5;

  float nlScaleForHalf(float phi_half){
    // scale such that |2*phi_half| <= PHI_NL_MAX
    float a = abs(phi_half) * 2.0;
    return (a > PHI_NL_MAX) ? (PHI_NL_MAX / a) : 1.0;
  }

  out vec4 fragColor;

  /* ==========================================================
     A) Math & Field Helpers (NO optics knowledge)
     ========================================================== */

  bool onInteriorInterface(vec2 pos) {
    vec2 center = 0.5 * u_resolution;
    vec2 rel = pos - center;
    float r = length(rel);
    float theta = atan(rel.y, rel.x);
    float R = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));
    // 0.5 = about one pixel thickness
    return (r >= R - 0.5) && (r <= R + 0.5);
  }


  // ---------- Complex & small vector helpers ----------
  vec2 cmul(vec2 a, vec2 b) { return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }
  vec2 cconj(vec2 a) { return vec2(a.x, -a.y); }

  vec2 vmin(vec2 a, vec2 b){ return vec2(min(a.x,b.x), min(a.y,b.y)); }
  vec2 vmax(vec2 a, vec2 b){ return vec2(max(a.x,b.x), max(a.y,b.y)); }

  vec2 sampleField(sampler2D tex, vec2 uv) {
    return texture(tex, uv).xy;
  }

  vec2 dirFrom(float ang){ return vec2(cos(ang), sin(ang)); }

  // Clamp a complex value’s magnitude to neighborhood max (phase-preserving)
  vec2 clampMagnitude(vec2 val, sampler2D tex, vec2 uv_dep, vec2 texel){
    float maxMag = 0.0;
    vec2 pts[5] = vec2[5](
      uv_dep,
      uv_dep + vec2( texel.x, 0.0),
      uv_dep + vec2(-texel.x, 0.0),
      uv_dep + vec2(0.0,  texel.y),
      uv_dep + vec2(0.0, -texel.y)
    );
    for (int i=0;i<5;++i){
      maxMag = max(maxMag, length(sampleField(tex, pts[i])));
    }
    float L = length(val);
    if (L > maxMag + 1e-12) return val * (maxMag / max(L, 1e-12));
    return val;
  }

  // ---------- Local statistics ----------
  float localIavg(vec2 pos, vec2 texel, sampler2D tex){
    float s = 0.0;
    for(int j=-1;j<=1;++j){
      for(int i=-1;i<=1;++i){
        vec2 E = texture(tex, (pos+vec2(float(i),float(j)))*texel).xy;
        s += dot(E,E);
      }
    }
    return s/9.0;
  }

  // ---------- Donor-cell limiter neighborhood (5-point cross) ----------
  void donorMinMax(sampler2D tex, vec2 uv_dep, vec2 texel, out vec2 vmin2, out vec2 vmax2){
    vec2 c  = sampleField(tex, uv_dep);
    vec2 px = sampleField(tex, uv_dep + vec2( texel.x, 0.0));
    vec2 nx = sampleField(tex, uv_dep + vec2(-texel.x, 0.0));
    vec2 py = sampleField(tex, uv_dep + vec2(0.0,  texel.y));
    vec2 ny = sampleField(tex, uv_dep + vec2(0.0, -texel.y));
    vmin2 = vmin(vmin(vmin(vmin(c,px),nx),py),ny);
    vmax2 = vmax(vmax(vmax(vmax(c,px),nx),py),ny);
  }

  // ---------- Gradients (auto 2nd/4th order with intensity gating) ----------
  void grad2_E(vec2 pos, vec2 texel, sampler2D tex, out vec2 dEx, out vec2 dEy){
    vec2 Ex_p = texture(tex,(pos+vec2( 1.0, 0.0))*texel).xy;
    vec2 Ex_m = texture(tex,(pos+vec2(-1.0, 0.0))*texel).xy;
    vec2 Ey_p = texture(tex,(pos+vec2( 0.0, 1.0))*texel).xy;
    vec2 Ey_m = texture(tex,(pos+vec2( 0.0,-1.0))*texel).xy;
    dEx = 0.5*(Ex_p - Ex_m) / u_dx;
    dEy = 0.5*(Ey_p - Ey_m) / u_dx;
  }

  void grad4_E(vec2 pos, vec2 texel, sampler2D tex, out vec2 dEx, out vec2 dEy){
    vec2 Ex2  = texture(tex,(pos+vec2( 2.0, 0.0))*texel).xy;
    vec2 Ex1  = texture(tex,(pos+vec2( 1.0, 0.0))*texel).xy;
    vec2 Ex_1 = texture(tex,(pos+vec2(-1.0, 0.0))*texel).xy;
    vec2 Ex_2 = texture(tex,(pos+vec2(-2.0, 0.0))*texel).xy;

    vec2 Ey2  = texture(tex,(pos+vec2( 0.0, 2.0))*texel).xy;
    vec2 Ey1  = texture(tex,(pos+vec2( 0.0, 1.0))*texel).xy;
    vec2 Ey_1 = texture(tex,(pos+vec2( 0.0,-1.0))*texel).xy;
    vec2 Ey_2 = texture(tex,(pos+vec2( 0.0,-2.0))*texel).xy;

    dEx = (-Ex2 + 8.0*Ex1 - 8.0*Ex_1 + Ex_2) * (1.0 / (12.0 * u_dx));
    dEy = (-Ey2 + 8.0*Ey1 - 8.0*Ey_1 + Ey_2) * (1.0 / (12.0 * u_dx));
  }

  void grad_auto(vec2 pos, vec2 texel, sampler2D tex, vec2 E0, out vec2 dEx, out vec2 dEy){
    bool nearEdge = (pos.x < 2.5) || (pos.y < 2.5)
                 || (pos.x > u_resolution.x - 2.5)
                 || (pos.y > u_resolution.y - 2.5);

    if (!nearEdge){
      vec2 d4x, d4y; grad4_E(pos, texel, tex, d4x, d4y);
      vec2 d2x, d2y; grad2_E(pos, texel, tex, d2x, d2y);

      float I  = dot(E0,E0);
      float Ia = localIavg(pos, texel, tex);
      float I0 = max(1e-12, 0.1 * Ia);
      float w  = smoothstep(I0, 10.0*I0, I);
      dEx = mix(d2x, d4x, w);
      dEy = mix(d2y, d4y, w);
    } else {
      grad2_E(pos, texel, tex, dEx, dEy);
    }
  }

  // Phase gradient ∇φ from complex field (u,v) with central differences
  vec2 phaseGradientFrom(sampler2D tex, vec2 pos, vec2 texel) {
    vec2 E   = texture(tex, pos*texel).xy;
    vec2 Ex  = texture(tex, (pos+vec2(1.0,0.0))*texel).xy;
    vec2 Ey  = texture(tex, (pos+vec2(0.0,1.0))*texel).xy;

    // dφ/dx ≈ Im{ ln(Ex * conj(E)) } / dx   (and same for y)
    vec2 rx  = cmul(Ex, cconj(E));
    vec2 ry  = cmul(Ey, cconj(E));
    float dphix = atan(rx.y, rx.x) / u_dx;
    float dphiy = atan(ry.y, ry.x) / u_dx;

    return vec2(dphix, dphiy);
  }

  // ---------- MacCormack/BFECC advection for complex fields ----------
  vec2 macCormackAdvect(sampler2D tex, vec2 pos, vec2 texel, vec2 offsetPixels){
    vec2 uv   = pos * texel;
    vec2 off  = offsetPixels * texel;

    // Backward trace (donor)
    vec2 uv_dep  = uv - off;
    vec2 phi_dep = sampleField(tex, uv_dep);

    // Forward trace from donor to estimate truncation error
    vec2 uv_fwd  = uv_dep + off;
    vec2 phi_fwd = sampleField(tex, uv_fwd);

    // BFECC correction
    vec2 phi0     = sampleField(tex, uv);
    vec2 phi_corr = phi_dep + 0.5 * (phi0 - phi_fwd);

    // Donor-cell limiter (component-wise)
    vec2 vmin2, vmax2;
    donorMinMax(tex, uv_dep, texel, vmin2, vmax2);
    phi_corr = clamp(phi_corr, vmin2, vmax2);

    float L0 = length(phi_dep);
    float Lc = length(phi_corr);
    float R_MAX = 2.5; // allow generous overshoot
    if (Lc > R_MAX * max(L0, 1e-12)) {
      phi_corr *= (R_MAX * L0) / Lc;
    }

    return phi_corr;
  }


  /* ==========================================================
     B) Optics Primitives (SINGLE SOURCE OF TRUTH)
     ========================================================== */

  // ---------- k-hat from phase ----------
  vec2 khatFrom(sampler2D tex, vec2 pos, vec2 texel){
    vec2 g = phaseGradientFrom(tex, pos, texel);
    float L = length(g);
    return (L > 1e-12) ? (g / L) : vec2(1.0, 0.0);
  }

  // ---------- Index vs. angle and wavelength ----------
  float n_angle(vec4 lensVal, vec2 khat) {
    float n_o = max(LENS_NO(lensVal), 1.0);
    float n_e = max(LENS_NE(lensVal), 1.0);
    if (abs(n_e - n_o) < 1e-6) return n_o;
    vec2 axis = vec2(cos(LENS_AXIS(lensVal)), sin(LENS_AXIS(lensVal)));
    float c = clamp(dot(khat, axis), -1.0, 1.0);
    float s2 = 1.0 - c*c;
    float denom = n_e*n_e*c*c + n_o*n_o*s2;
    return (n_o*n_e) / sqrt(max(denom, 1e-12));
  }

  float getWavelengthDependentIndex(float baseIndex, float dispersion, float wavelength) {
    float lambda_ref = u_lambdaFund;
    return baseIndex + dispersion * (1.0/(wavelength*wavelength) - 1.0/(lambda_ref*lambda_ref));
  }

  // ---------- ADE dispersion (Lorentz, per-pixel) ----------
  void lorentzADE_step(
    bool isFundamental, vec2 pos, vec2 texel,
    vec2 E_in,      // field used to drive P at this step
    vec2 P, vec2 P_prev,
    out vec2 P_new, out vec2 dE_disp
  ){
    vec4 Pm = (isFundamental ? DISP_FUND : DISP_SHG);
    float W0    = max(1e-6, Pm.x);
    float GAMMA = max(0.0,  Pm.y);
    float ALPHA0= Pm.z;
    float KAPPA = Pm.w;

    // stability clamps
    float wdt = W0 * u_dt; if (wdt > 1.8) W0 *= 1.8 / wdt;
    float gdt = GAMMA * u_dt; if (gdt > 1.0) GAMMA = 1.0 / u_dt;

    float alpha = ALPHA0 * LENS_DISP(texture(u_lens, pos*texel));
    float dt  = u_dt, dt2 = dt*dt;

    // (P^{n+1} - 2P^n + P^{n-1})/dt^2 + Γ (P^n - P^{n-1})/dt + W0^2 P^n = α E^n
    float A = 2.0 - (W0*W0)*dt2 - GAMMA*dt;
    float B = -1.0 + GAMMA*dt;
    vec2  C = (alpha * dt2) * E_in;

    P_new  = A*P + B*P_prev + C;
    vec2 P_tt = (P_new - 2.0*P + P_prev) / dt2;

    // wave equation uses: E^{n+1} += (-κ) dt^2 P_tt
    dE_disp = (-KAPPA) * dt2 * P_tt;
  }

  // ---------- Uniaxial walk-off primitives ----------
  float n_e_of_theta(float no, float ne, float cost, float sint){
    float D = ne*ne*cost*cost + no*no*sint*sint;
    return (no*ne) / sqrt(max(D, 1e-12));
  }

  float tanRho_uniaxial(float no, float ne, float cost, float sint){
    float D = ne*ne*cost*cost + no*no*sint*sint;
    return ((no*no - ne*ne) * sint * cost) / max(D, 1e-12);
  }

  vec2 walkoff_dir(vec2 khat, vec2 axis){
    vec2 t = axis - khat * dot(axis, khat);
    float L = length(t);
    if (L < 1e-12) return vec2(-khat.y, khat.x);
    return t / L;
  }

  // Offset in pixels for SHG walk-off during one time step
  vec2 computeWalkoffOffset_SHG(vec2 pos, vec2 texel){
    // local propagation from current band (SHG pass should bind u_current=u_shg)
    vec2 g = phaseGradientFrom(u_current, pos, texel);
    float gL = length(g);
    if (gL <= 1e-12) return vec2(0.0);
    vec2 khat = g / gL;

    vec4 lv   = texture(u_lens, pos * texel);
    float no  = max(LENS_NO(lv), 1.0);
    float ne  = max(LENS_NE(lv), 1.0);
    vec2 axis = vec2(cos(LENS_AXIS(lv)), sin(LENS_AXIS(lv)));

    float cost = clamp(dot(khat, axis), -1.0, 1.0);
    float sint = sqrt(max(0.0, 1.0 - cost*cost));

    float tanR    = abs(tanRho_uniaxial(no, ne, cost, sint));
    float n_base  = n_e_of_theta(no, ne, cost, sint); // extraordinary branch
    float stepLong = (u_c / n_base) * u_dt / u_dx;

    float tanEff = min(tanR, WALK_TAN_CAP);
    vec2 t_hat   = walkoff_dir(khat, axis);
    vec2 offset  = stepLong * tanEff * t_hat; // pixels

    // pixel cap
    float L = length(offset);
    if (L > WALK_MAX_PIX) offset *= (WALK_MAX_PIX / L);
    return offset;
  }

  // ---------- Phase mismatch per time-step (time-domain consistent) ----------
  vec4 calculatePhaseMismatchTerm(vec2 pos, vec2 texel, vec4 /*center*/) {
    vec2 kF = khatFrom(u_fundamental, pos, texel);
    vec2 kS = khatFrom(u_shg,         pos, texel);

    vec4 lensVal = texture(u_lens, pos*texel);

    // Angle-dependent effective indices at the *actual* directions
    float nF_base = n_angle(lensVal, kF);
    float nS_base = n_angle(lensVal, kS);

    float nF = getWavelengthDependentIndex(nF_base, LENS_DISP(lensVal), u_lambdaFund);
    float nS = getWavelengthDependentIndex(nS_base, LENS_DISP(lensVal), u_lambdaSHG);

    float kFmag = 2.0*PI*nF / u_lambdaFund;
    float kSmag = 2.0*PI*nS / u_lambdaSHG;

    // Vectorial mismatch projected along kF direction
    float cosTheta = clamp(dot(kF, kS), -1.0, 1.0);
    float deltaK   = 2.0*kFmag - kSmag * cosTheta;

    float dz    = (u_c / nF) * u_dt;
    float phase = deltaK * dz + u_phaseRef;
    return vec4(cos(phase), sin(phase), 0.0, 0.0);
  }


  // ---------- Material properties ----------
  struct FieldState {
    vec2  field;
    vec2  oldField;
    float amp2;
    float wavelength;
  };

  struct MaterialProps {
    float n_base;
    float n_eff;
    float c2_eff;
    float chi2_local;
    float chi3_local;
    float dispCoeff;
  };

  float angleDependentBaseIndex(vec2 pos, vec2 texel){
    vec4 lensVal = texture(u_lens, pos*texel);
    vec2 kh = khatFrom(u_current, pos, texel);
    return n_angle(lensVal, kh);
  }

  MaterialProps computeMaterialProperties(vec2 pos, vec2 texel, FieldState state, float otherAmp2) {
    MaterialProps props;

    vec4 lensVal = texture(u_lens, pos * texel);
    vec4 gm      = texture(u_gainMask, pos * texel);

    // Base index for NL phase rotation (angle-dependent, from current band's khat)
    props.n_base    = angleDependentBaseIndex(pos, texel);
    props.dispCoeff = LENS_DISP(lensVal);

    // Nonlinearity weights (mask channels clamped non-negative)
    float chi3w = max(gm.b, 0.0);
    float chi2w = gm.g;
    props.chi3_local = u_chi * u_chi * u_chi_ratio * chi3w;
    props.chi2_local = u_chi * u_chi2_ratio       * chi2w;

    // Nonlocal Kerr intensity (band-local blur)
    float u_nonlocalKerr = 0.0; // temp constant; wire to uniform later
    float I_local = state.amp2;
    float I_blur  = localIavg(pos, texel, u_current);
    float I_eff   = mix(I_local, I_blur, clamp(u_nonlocalKerr, 0.0, 1.0));

    // Kerr (self + cross) with saturation (sat keyed to primary intensity)
    float kerrSat   = 1.0 / (1.0 + I_eff / u_kerr_Isat);
    float kerrSelf  = props.chi3_local * I_eff * kerrSat;
    float kerrCross = u_crossKerrCoupling * props.chi3_local * otherAmp2 * kerrSat;

    props.n_eff = max(props.n_base + kerrSelf + kerrCross, 1e-6);

    // Linear wave-speed coefficient: use ORDINARY index only (self-consistent with varCoeffLaplacian)
    float n_o   = max(LENS_NO(lensVal), 1.0);
    float c_lin = u_c / n_o;
    props.c2_eff = c_lin * c_lin;

    return props;
  }


  // ---------- χ(2) conversion terms ----------
  struct ConversionTerms {
    vec2 upConversion;
    vec2 downConversion;
    float phaseMatch; // diagnostic scalar
  };

  ConversionTerms computeConversionTerms(
    vec2 pos, vec2 texel,
    FieldState fund, FieldState shg,
    MaterialProps props, bool isFundamental
  ){
    ConversionTerms T;

    // phase mismatch e^{+iΔkΔz}
    vec2 ph    = calculatePhaseMismatchTerm(pos, texel, vec4(fund.field, 0.0, 0.0)).xy;
    vec2 phConj= vec2(ph.x, -ph.y);

    float sat = 1.0 / (1.0 + dot(fund.field, fund.field) / u_shg_Isat);
    float k   = u_conversionCoupling * props.chi2_local;

    if (isFundamental){
      // fundamental depletion by back-conversion: ~ E2 * E1* * e^{-iΔkΔz}
      vec2 term = cmul(shg.field, cmul(cconj(fund.field), phConj));
      T.upConversion   = vec2(0.0);
      T.downConversion = -k * term * sat;
    } else {
      // SHG build-up: ~ E1^2 * e^{+iΔkΔz}
      vec2 term = cmul(cmul(fund.field, fund.field), ph);
      T.upConversion   =  k * term * sat;
      T.downConversion = vec2(0.0);
    }

    // diagnostic: alignment of e^{iΔkΔz} with E1^2
    vec2 e1sq = cmul(fund.field, fund.field);
    float L = length(e1sq);
    vec2 e1sq_hat = (L > 1e-12) ? (e1sq / L) : vec2(1.0, 0.0);
    T.phaseMatch = clamp(dot(ph, e1sq_hat), -1.0, 1.0);

    return T;
  }

  /* ==========================================================
     C) Linear Operator & Boundary Subsystem
     ========================================================== */

  /* ------------------ C1) Linear Operator (var-coeff Laplacian) ------------------ */

  // isotropic base wave-speed coefficient a(x) from lens (uses ordinary index)
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

  // conservative var-coeff Laplacian using face-averaged a(x)
  vec2 varCoeffLaplacian(vec2 pos, vec2 texel, float aC, float wavelength) {
    vec2 EC = texture(u_current, pos*texel).xy;
    vec2 EL = texture(u_current,(pos+vec2(-1.0,0.0))*texel).xy;
    vec2 ER = texture(u_current,(pos+vec2( 1.0,0.0))*texel).xy;
    vec2 EU = texture(u_current,(pos+vec2( 0.0,-1.0))*texel).xy;
    vec2 ED = texture(u_current,(pos+vec2( 0.0, 1.0))*texel).xy;

    vec4 lC=texture(u_lens,pos*texel);
    vec4 lL=texture(u_lens,(pos+vec2(-1,0))*texel);
    vec4 lR=texture(u_lens,(pos+vec2( 1,0))*texel);
    vec4 lU=texture(u_lens,(pos+vec2( 0,-1))*texel);
    vec4 lD=texture(u_lens,(pos+vec2( 0, 1))*texel);

    // wavelength-adjusted ordinary indices for neighbors (match center’s band)
    float nL_wl = getWavelengthDependentIndex(max(LENS_NO(lL),1.0), LENS_DISP(lL), wavelength);
    float nR_wl = getWavelengthDependentIndex(max(LENS_NO(lR),1.0), LENS_DISP(lR), wavelength);
    float nU_wl = getWavelengthDependentIndex(max(LENS_NO(lU),1.0), LENS_DISP(lU), wavelength);
    float nD_wl = getWavelengthDependentIndex(max(LENS_NO(lD),1.0), LENS_DISP(lD), wavelength);

    float aL=(u_c*u_c)/(nL_wl*nL_wl);
    float aR=(u_c*u_c)/(nR_wl*nR_wl);
    float aU=(u_c*u_c)/(nU_wl*nU_wl);
    float aD=(u_c*u_c)/(nD_wl*nD_wl);

    float aW=0.5*(aC + aL), aE=0.5*(aC + aR);
    float aN=0.5*(aC + aU), aS=0.5*(aC + aD);

    vec2 flux = aE*(ER-EC) - aW*(EC-EL) + aS*(ED-EC) - aN*(EC-EU);
    return flux / (u_dx*u_dx);
  }

  /* ------------------ C2) Boundary Subsystem ------------------ */

  // inward-pointing boundary normal (robust slightly off the curve)
  vec2 boundaryNormal(vec2 pos) {
    vec2 center = 0.5 * u_resolution;
    vec2 rel = pos - center;
    float r = max(length(rel), 1e-6);
    float theta = atan(rel.y, rel.x);

    float cosMTheta = cos(u_boundaryM * theta);
    float sinMTheta = sin(u_boundaryM * theta);

    float Rprime = -u_boundaryR0 * u_boundaryAlpha * u_boundaryM * sinMTheta;

    // unit polar basis at (r, theta)
    vec2 rHat     = rel / r;
    vec2 thetaHat = vec2(-rHat.y, rHat.x);

    // inward normal = -∇(r - R(θ)) normalized
    vec2 n = rHat - (Rprime / r) * thetaHat;
    return -normalize(n);
  }

  // ∂E/∂n at the boundary (component-wise, for complex E)
  vec2 normalDerivativeAt(vec2 pos, vec2 texel, vec2 E0){
    vec2 dEx, dEy;            // complex ∂E/∂x and ∂E/∂y
    grad_auto(pos, texel, u_current, E0, dEx, dEy);
    vec2 n_in = boundaryNormal(pos);        // inward-pointing normal
    return dEx * n_in.x + dEy * n_in.y;     // project gradient onto normal
  }

  // 1 inside, 0 far outside, linear ramp over transition width
  float boundaryDampingFactor(vec2 pos) {
    vec2 center = 0.5 * u_resolution;
    vec2 rel = pos - center;
    float r = length(rel);
    float theta = atan(rel.y, rel.x);

    float boundaryRadius = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));
    if (r <= boundaryRadius) return 1.0;

    float distPast = r - boundaryRadius;
    if (distPast >= u_boundaryTransitionWidth) return 0.0;

    float t = distPast / max(u_boundaryTransitionWidth, 1e-6);
    return 1.0 - t;
  }

  struct BoundaryProps {
    vec2  normal;      // inward normal
    float theta_i;     // incidence angle
    float n_local;     // local refractive index (inside)
    vec2  kvec;        // local phase gradient direction (not normalized magnitude)
    float amp2;        // |E|^2 at pixel
    bool  validField;  // is the field strong enough
  };

  BoundaryProps computeBoundaryProperties(vec2 pos, vec2 texel) {
    BoundaryProps props;

    vec4 centerVal = texture(u_current, pos * texel);
    props.amp2 = dot(centerVal.xy, centerVal.xy);

    const float MIN_FIELD_STRENGTH = 1e-9;
    props.validField = (props.amp2 >= MIN_FIELD_STRENGTH);

    if (!props.validField) {
      props.normal = vec2(0.0);
      props.theta_i = 0.0;
      props.n_local = 1.0;
      props.kvec = vec2(0.0);
      return props;
    }

    // local k-vector direction from current band
    props.kvec = phaseGradientFrom(u_current, pos, texel);

    // boundary normal
    props.normal = boundaryNormal(pos);

    // incidence angle (use absolute cosine to ignore sign)
    float kL = length(props.kvec);
    float cos_in = (kL > 1e-12) ? dot(props.kvec / kL, props.normal) : 0.0;
    props.theta_i = acos(clamp(abs(cos_in), 0.0, 1.0));


    // local inside refractive index (use ordinary component)
    vec4 lensVal = texture(u_lens, pos * texel);
    props.n_local = max(LENS_NO(lensVal), 1.0);

    return props;
  }

  // ---------- Fresnel: intensity reflectivity (TE/TM) ----------
  float fresnelReflectivity_TE(float n_in, float n_out, float theta_i) {
    n_in = max(n_in, 1e-6);
    n_out = max(n_out, 1e-6);
    theta_i = clamp(theta_i, 0.0, PI/2.0);

    float crit_angle = (n_out < n_in) ? asin(n_out/n_in) : PI;
    if (theta_i > crit_angle) {
      return 1.0; // TIR (intensity = 1)
    }

    float sin_t = (n_in / n_out) * sin(theta_i);
    sin_t = clamp(sin_t, -1.0, 1.0);
    float theta_t = asin(sin_t);

    float cos_i = cos(theta_i);
    float cos_t = cos(theta_t);

    const float eps = 1e-6;
    float r_num = (n_in * cos_i) - (n_out * cos_t);
    float r_den = (n_in * cos_i) + (n_out * cos_t) + eps;
    float r = r_num / r_den;
    return r * r;
  }

  float fresnelReflectivity_TM(float n_in, float n_out, float theta_i) {
    n_in = max(n_in, 1e-6);
    n_out = max(n_out, 1e-6);
    theta_i = clamp(theta_i, 0.0, PI/2.0);

    float crit_angle = (n_out < n_in) ? asin(n_out/n_in) : PI;
    if (theta_i > crit_angle) {
      return 1.0; // TIR
    }

    float sin_t = (n_in / n_out) * sin(theta_i);
    sin_t = clamp(sin_t, -1.0, 1.0);
    float theta_t = asin(sin_t);

    float cos_i = cos(theta_i);
    float cos_t = cos(theta_t);

    const float eps = 1e-6;
    float r_num = (n_out * cos_i) - (n_in * cos_t);
    float r_den = (n_out * cos_i) + (n_in * cos_t) + eps;
    float r = r_num / r_den;
    return r * r;
  }

  // Complex Fresnel with phase (scalar polarization compromise)
  vec2 fresnelR_complex_TE(float n1, float n2, float theta_i){
    // returns complex r = r_re + i r_im in vec2
    float s = sin(theta_i), c = cos(theta_i);
    float k = n1/n2;
    float st = k*s;
    float under = 1.0 - st*st;

    if (under >= 0.0) { // partial reflection, real cos(theta_t)
      float ct = sqrt(max(under, 0.0));
      float r = (n1*c - n2*ct) / (n1*c + n2*ct);
      return vec2(r, 0.0);
    } else { // TIR: cos(theta_t) = i*eta, |r|=1, phase != 0
      float eta = sqrt(-under);
      // r = (n1*c - i*n2*eta)/(n1*c + i*n2*eta) = e^{i*phi}
      float num_re =  n1*c;
      float num_im = -n2*eta;
      float den_re =  n1*c;
      float den_im =  n2*eta;
      // complex division (num/den)
      float den2   = den_re*den_re + den_im*den_im + 1e-12;
      float r_re = (num_re*den_re + num_im*den_im) / den2;
      float r_im = (num_im*den_re - num_re*den_im) / den2;
      // normalize to avoid drift
      float L = inversesqrt(max(r_re*r_re + r_im*r_im, 1e-12));
      return vec2(r_re, r_im) * L;
    }
  }

  vec2 fresnelR_complex_TM(float n1, float n2, float theta_i){
    float s = sin(theta_i), c = cos(theta_i);
    float k = n1/n2;
    float st = k*s;
    float under = 1.0 - st*st;

    if (under >= 0.0) {
      float ct = sqrt(max(under, 0.0));
      float r = (n2*c - n1*ct) / (n2*c + n1*ct);
      return vec2(r, 0.0);
    } else {
      float eta = sqrt(-under);
      float num_re =  n2*c;
      float num_im = -n1*eta;
      float den_re =  n2*c;
      float den_im =  n1*eta;
      float den2   = den_re*den_re + den_im*den_im + 1e-12;
      float r_re = (num_re*den_re + num_im*den_im) / den2;
      float r_im = (num_im*den_re - num_re*den_im) / den2;
      float L = inversesqrt(max(r_re*r_re + r_im*r_im, 1e-12));
      return vec2(r_re, r_im) * L;
    }
  }

  // Simple average of TE/TM complex r (keeps phase)
  vec2 complexReflectionAvg(BoundaryProps p){
    vec2 rTE = fresnelR_complex_TE(p.n_local, 1.0, p.theta_i);
    vec2 rTM = fresnelR_complex_TM(p.n_local, 1.0, p.theta_i);
    return 0.5*(rTE + rTM);
  }

  // Leapfrog-consistent characteristic reflection at the interior interface.
  // Uses inward normal n_in. Acts only on the **inner** side of the ring.
  vec2 applyBoundaryCharacteristic(vec2 predictedNextE, vec2 pos, vec2 texel, float damp){
    // Only correct exactly on the ring, and only from the interior side.
    if (!onInteriorInterface(pos) || damp < 0.999) return predictedNextE;

    // Local fields at t^n and t^{n-1}
    vec2 En   = texture(u_current,  pos*texel).xy;
    vec2 Enm1 = texture(u_previous, pos*texel).xy;

    // ∂E/∂t and ∂E/∂n at t^n
    vec2 dEdt  = (En - Enm1) / u_dt;
    vec2 dEndn = normalDerivativeAt(pos, texel, En);

    // Boundary optics
    BoundaryProps p = computeBoundaryProperties(pos, texel);
    if (!p.validField) return predictedNextE;

    float c_n = u_c / max(p.n_local, 1e-6);

    // Average complex Fresnel (guard magnitude ≤ 1)
    vec2 r = complexReflectionAvg(p);
    float rL = length(r);
    if (rL > 1.0) r *= (0.999 / rL);

    // Optional outward-only hysteresis to avoid grazing flicker
    bool outwardLikely = false;
    float kL = length(p.kvec);
    if (kL > 1e-6) {
      float cos_in = dot(p.kvec / kL, p.normal); // inward normal
      outwardLikely = (cos_in < -0.02);          // negative = outward
    }

    // Characteristics wrt inward normal:
    // W_out = E_t - c ∂_n E  (outward)
    // W_in  = E_t + c ∂_n E  (inward)
    vec2 W_out = dEdt - c_n * dEndn;
    vec2 W_in_new = outwardLikely ? cmul(r, W_out) : vec2(0.0);

    // Recombine to boundary-consistent **time-centered** derivative at t^n
    vec2 dEdt_bc = 0.5 * (W_out + W_in_new);

    // Leapfrog-consistent next value: (E^{n+1} - E^{n-1})/(2Δt) = E_t^n
    vec2 E_bc = Enm1 + 2.0 * u_dt * dEdt_bc;

    // Blend with interior predictor (set <1 if you still see edge chatter)
    const float CHAR_BLEND = 1.0; // e.g. 0.85 if needed
    return mix(predictedNextE, E_bc, CHAR_BLEND);
  }


  /* ==========================================================
     D) Band Kernels & Main
     ========================================================== */

  // Core per-pixel update for one band with Strang-split NL
  vec4 computeInteriorField(vec2 pos, vec2 texel, bool isFundamental, float damp) {
    // ---- Gather state for this band ----
    vec4 center   = texture(u_current,  pos * texel);   // this band, t^n (E.xy, P.zw)
    vec4 oldField = texture(u_previous, pos * texel);   // this band, t^{n-1}

    FieldState primary, other;
    vec4 otherTex;

    if (isFundamental) {
      primary = FieldState(center.xy, oldField.xy, dot(center.xy, center.xy), u_lambdaFund);
      otherTex = texture(u_shg, pos * texel);
      other   = FieldState(otherTex.xy, vec2(0.0), dot(otherTex.xy, otherTex.xy), u_lambdaSHG);
    } else {
      primary = FieldState(center.xy, oldField.xy, dot(center.xy, center.xy), u_lambdaSHG);
      otherTex = texture(u_fundamental, pos * texel);
      other   = FieldState(otherTex.xy, vec2(0.0), dot(otherTex.xy, otherTex.xy), u_lambdaFund);
    }

    // === Pre NL half-step: Kerr rotation + χ(2) half-kick ===
    MaterialProps props_pre = computeMaterialProperties(pos, texel, primary, other.amp2);

    ConversionTerms conv_pre;
    if (isFundamental) { conv_pre = computeConversionTerms(pos, texel, primary, other,  props_pre, true);
    } else             { conv_pre = computeConversionTerms(pos, texel, other,  primary, props_pre, false); }

    float phiNL_half = (2.0 * PI / primary.wavelength) * (props_pre.n_eff - props_pre.n_base) * (0.5 * u_dt);
    float sNL_pre = nlScaleForHalf(phiNL_half);
    float dtNL_pre = u_dt * sNL_pre;
    float phiNL_half_limited = phiNL_half * sNL_pre;

    vec2 rotH = vec2(cos(phiNL_half_limited), sin(phiNL_half_limited));
    vec2 field_nl = vec2(
      primary.field.x * rotH.x - primary.field.y * rotH.y,
      primary.field.x * rotH.y + primary.field.y * rotH.x
    );
    field_nl += 0.5 * (dtNL_pre * dtNL_pre) * (conv_pre.upConversion + conv_pre.downConversion);

    // === Amplitude (gain + nonlinear absorption) as multiplicative split (pre half) ===
    float Ipre = dot(field_nl, field_nl);
    float maskR_raw = texture(u_gainMask, pos * texel).r;
    // fade out gain as soon as the PML starts; 0.2 is a gentle taper
    float maskR = maskR_raw * smoothstep(0.2, 1.0, damp);
    float g_lin = (u_gain0 / (1.0 + Ipre / u_gainSat) - u_linearLoss) * maskR;

    // Physical absorption placeholders off (purely conservative unless gainMask says otherwise)
    float u_betaTPA  = 0.0;
    float u_gamma3PA = 0.0;
    float g_nl  = -(u_betaTPA) * Ipre - (u_gamma3PA) * Ipre * Ipre;

    float amp_pre = exp(0.5 * u_dt * (g_lin + g_nl));
    vec2  linBase  = field_nl * amp_pre;

    // === Linear full-step: conservative var-coeff Laplacian + centered ADE ===
    vec2  linOld = primary.oldField;

    // a(x) from wavelength-adjusted ordinary index at center
    vec4  lC      = texture(u_lens, pos*texel);
    float nC_o     = max(LENS_NO(lC), 1.0);
    float nC_wl    = getWavelengthDependentIndex(nC_o, LENS_DISP(lC), primary.wavelength);
    float aC_phys  = (u_c*u_c) / (nC_wl*nC_wl);

    // CFL safety scaling using the same dispersed speed used in aC_phys
    float c_eff = u_c / nC_wl;
    float s_cfl = min(1.0, (1.0 / sqrt(2.0)) * (u_dx / (c_eff * u_dt)));

    // spatial operator at t^n
    vec2 flux = varCoeffLaplacian(pos, texel, aC_phys, primary.wavelength);

    // Blend with isotropic 9-pt Laplacian when a(x) is locally uniform
    const float LAP9_WEIGHT = 0.6;
    const float UNIF_EPS    = 0.02;
    float relUniform = localUniformity(pos, texel);
    float wUniform = 1.0 - smoothstep(UNIF_EPS, 2.0*UNIF_EPS, relUniform);
    vec2  flux9 = laplacian9_const(pos, texel, aC_phys);
    float w = clamp(LAP9_WEIGHT * wUniform, 0.0, 0.95);
    s_cfl *= (1.0 - 0.2*w); // small safety since 9-pt is stiffer near Nyquist
    flux = mix(flux, flux9, w);

    // ADE dispersion forcing from E^n (linBase), included on RHS
    vec2 P_new, dEdisp;
    lorentzADE_step(isFundamental, pos, texel, linBase, center.zw, oldField.zw, P_new, dEdisp);

    vec2 newField = (2.0 * linBase - linOld) + s_cfl * (u_dt * u_dt) * (flux + dEdisp);

    // Optional pulse injection (fundamental only)
    if (isFundamental && u_pulseInterval > 0 && (u_frameCount % u_pulseInterval == 0)) {
      float w0 = u_beamWidth;
      vec2 centerGrid = 0.5 * u_resolution;
      vec2 r = pos - centerGrid;                // renamed to avoid shadowing
      float r2 = dot(r, r);
      float amplitude = u_subsequentPulseAmplitude * exp(-r2 / (w0 * w0));
      newField += amplitude * vec2(1.0, 0.0);   // real pulse
    }

    // SHG walk-off transport (no ADE here; ADE already in RHS)
    if (!isFundamental){
      vec2 shgHere = texture(u_shg, pos*texel).xy;                   // explicit SHG
      float I_here = dot(shgHere, shgHere);
      if (I_here > 1e-6) {
        vec2 offset = computeWalkoffOffset_SHG(pos, texel);
        if (offset.x != 0.0 || offset.y != 0.0){
          vec2 adv = macCormackAdvect(u_shg, pos, texel, offset);    // advect SHG
          float Ia = localIavg(pos, texel, u_shg);                   // SHG neighborhood
          float I0 = max(1e-12, 0.1 * Ia);
          float wadv  = smoothstep(I0, 10.0*I0, I_here);
          newField = mix(newField, adv, WALK_BLEND * wadv);
        }
      }
    }

    // PML-lite (matched damping on time derivative)
    float sigma = (1.0 - damp) * (1.0 - damp) * u_damping; // quadratic ramp
    vec2 dEdt   = (center.xy - oldField.xy) / u_dt;
    newField   -= u_dt * sigma * dEdt;
    // (No extra amplitude multiply by u_damping here)

    // === Post NL half-step: Kerr rotation + χ(2) half-kick ===
    FieldState updatedPrimary = FieldState(newField, primary.oldField, dot(newField,newField), primary.wavelength);
    MaterialProps props_post = computeMaterialProperties(pos, texel, updatedPrimary, other.amp2);

    float phiNL_half_post = (2.0 * PI / updatedPrimary.wavelength) * (props_post.n_eff - props_post.n_base) * (0.5 * u_dt);
    float sNL_post = nlScaleForHalf(phiNL_half_post);
    float dtNL_post = u_dt * sNL_post;
    float phiNL_half_post_limited = phiNL_half_post * sNL_post;

    ConversionTerms conv_post;
    if (isFundamental) { conv_post = computeConversionTerms(pos, texel, updatedPrimary, other,  props_post, true);
    } else             { conv_post = computeConversionTerms(pos, texel, other,          updatedPrimary, props_post, false); }

    vec2 rotH2 = vec2(cos(phiNL_half_post_limited), sin(phiNL_half_post_limited));
    newField = vec2(
      newField.x * rotH2.x - newField.y * rotH2.y,
      newField.x * rotH2.y + newField.y * rotH2.x
    );
    newField += 0.5 * (dtNL_post * dtNL_post) * (conv_post.upConversion + conv_post.downConversion);

    // === Amplitude split (post half) ===
    float Ipost = dot(newField, newField);
    float g_lin_p = (u_gain0 / (1.0 + Ipost / u_gainSat) - u_linearLoss) * maskR;
    float g_nl_p  = -(u_betaTPA) * Ipost - (u_gamma3PA) * Ipost * Ipost;
    float amp_post = exp(0.5 * u_dt * (g_lin_p + g_nl_p));
    newField *= amp_post;

    // writeback: .xy = field, .zw = ADE polarization (from centered update)
    return vec4(newField, P_new);
  }



  void main() {
    vec2 pos   = gl_FragCoord.xy;
    vec2 texel = 1.0 / u_resolution;
    float damp = boundaryDampingFactor(pos);

    if (damp <= 0.0) { fragColor = vec4(0.0); return; }

    vec4 fieldUpdate = computeInteriorField(pos, texel, u_updateTarget == 0, damp);

    // Reflect only on the **interior** interface pixels (damp≈1), not in the sponge.
    if (onInteriorInterface(pos) && damp >= 0.999) {
      fieldUpdate.xy = applyBoundaryCharacteristic(fieldUpdate.xy, pos, texel, damp);
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

  // --- Build the base 1-channel gain mask (this is "baseRMask") ---
  const baseRMask = createGainMaskData(config); // length = gridSize*gridSize

  // --- Expand to RGBA for the shader's u_gainMask texture (initial seed) ---
  const rgbaData = new Float32Array(config.gridSize * config.gridSize * 4);
  for (let i = 0; i < baseRMask.length; i++) {
    rgbaData[i * 4] = baseRMask[i]; // R channel (linear gain)
    rgbaData[i * 4 + 1] = 0.00001; // G: initial chi2 weight (tiny seed)
    rgbaData[i * 4 + 2] = 0.00001; // B: initial chi3 weight (tiny seed)
    rgbaData[i * 4 + 3] = 1.0; // A
  }

  // Return both: the RGBA texture data and the raw 1-channel base mask
  return {
    shgData,
    fundamentalData,
    gainMaskData: rgbaData, // unchanged key for your existing uploader
    baseRMask, // <-- pass this to optimizer.getGainMaskData(...)
  };
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
