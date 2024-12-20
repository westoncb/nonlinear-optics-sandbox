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

  uniform sampler2D u_current;
  uniform sampler2D u_previous;
  uniform sampler2D u_lens;
  uniform sampler2D u_fundamental;
  uniform float u_dt;
  uniform float u_dx;
  uniform float u_damping;
  uniform float u_c;
  uniform float u_chi;
  uniform float u_chi_ratio;
  uniform float u_shg_Isat;
  uniform float u_kerr_Isat;
  uniform vec2 u_resolution;
  uniform float u_boundaryR0;
  uniform float u_boundaryAlpha;
  uniform float u_boundaryM;
  uniform int u_updateTarget;
  uniform int u_frameCount;
  uniform int u_pulseInterval;

  const float PI = 3.14159265358979323846;

  out vec4 fragColor;

  bool insideBoundary(vec2 pos) {
      vec2 center = u_resolution * 0.5;
      vec2 rel = pos - center;
      float r = length(rel);
      float theta = atan(rel.y, rel.x);
      float boundaryRadius = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));
      return r <= boundaryRadius;
  }

  vec2 calculatePhaseGradients(vec2 pos, vec2 texel, vec4 center) {
      vec4 right = texture(u_current, (pos + vec2(1.0, 0.0)) * texel);
      vec4 up = texture(u_current, (pos + vec2(0.0, -1.0)) * texel);
      vec4 left = texture(u_current, (pos + vec2(-1.0, 0.0)) * texel);
      vec4 down = texture(u_current, (pos + vec2(0.0, 1.0)) * texel);

      float dphase_dx = (atan(right.y, right.x) - atan(left.y, left.x)) / (2.0 * u_dx);
      float dphase_dy = (atan(down.y, down.x) - atan(up.y, up.x)) / (2.0 * u_dx);

      return vec2(dphase_dx, dphase_dy);
  }

  vec2 ninePointLaplacian(vec2 pos, vec2 texel, vec4 center) {
      // Direct neighbors (sharing edges) - weight 2/3
      vec4 left   = texture(u_current, (pos + vec2(-1.0,  0.0)) * texel);
      vec4 right  = texture(u_current, (pos + vec2( 1.0,  0.0)) * texel);
      vec4 up     = texture(u_current, (pos + vec2( 0.0, -1.0)) * texel);
      vec4 down   = texture(u_current, (pos + vec2( 0.0,  1.0)) * texel);

      // Diagonal neighbors - weight 1/6
      vec4 upleft    = texture(u_current, (pos + vec2(-1.0, -1.0)) * texel);
      vec4 upright   = texture(u_current, (pos + vec2( 1.0, -1.0)) * texel);
      vec4 downleft  = texture(u_current, (pos + vec2(-1.0,  1.0)) * texel);
      vec4 downright = texture(u_current, (pos + vec2( 1.0,  1.0)) * texel);

      return (
          (left.xy + right.xy + up.xy + down.xy) * (2.0/3.0) +
          (upleft.xy + upright.xy + downleft.xy + downright.xy) * (1.0/6.0) +
          center.xy * -4.0
      );
  }

  vec2 fivePointLaplacian(vec2 pos, vec2 texel, vec4 center) {
      vec4 left = texture(u_current, (pos + vec2(-1.0, 0.0)) * texel);
      vec4 right = texture(u_current, (pos + vec2(1.0, 0.0)) * texel);
      vec4 up = texture(u_current, (pos + vec2(0.0, -1.0)) * texel);
      vec4 down = texture(u_current, (pos + vec2(0.0, 1.0)) * texel);

      return (left.xy + right.xy + up.xy + down.xy - 4.0 * center.xy) / (u_dx * u_dx);
  }

  float calculateRefractiveIndex(float wavelength, float temperature) {
      float B1 = 1.03961212;
      float B2 = 0.231792344;
      float B3 = 1.01046945;
      float C1 = 0.00600069867;
      float C2 = 0.0200179144;
      float C3 = 103.560653;

      float lambda2 = wavelength * wavelength;
      float n2 = 1.0 + (B1 * lambda2 / (lambda2 - C1)) +
                 (B2 * lambda2 / (lambda2 - C2)) +
                 (B3 * lambda2 / (lambda2 - C3));

      return sqrt(n2);
  }

  vec4 calculatePhaseMismatchTerm(vec2 pos, vec2 texel) {
        float wavelength_fund = 40.0;
        float wavelength_shg = wavelength_fund * 0.5;  // = 10.0
        const float temperature = 22.0;  // Room temperature in Celsius

        float n_fund = calculateRefractiveIndex(wavelength_fund, temperature);
        float n_shg = calculateRefractiveIndex(wavelength_shg, temperature);

        float k_fund_real = 2.0 * PI * n_fund / wavelength_fund;
        float k_fund_imag = 0.0;
        float k_shg_real = 2.0 * PI * n_shg / wavelength_shg;
        float k_shg_imag = 0.0;

        float delta_k_real = k_shg_real - 2.0 * k_fund_real;
        float delta_k_imag = k_shg_imag - 2.0 * k_fund_imag;

        float qpm_period = 2.0 * PI / sqrt(delta_k_real * delta_k_real + delta_k_imag * delta_k_imag);
        float qpm_modulation_real = cos(pos.x / qpm_period * delta_k_real);
        float qpm_modulation_imag = sin(pos.x / qpm_period * delta_k_imag);

        return vec4(qpm_modulation_real, qpm_modulation_imag, 0.0, 1.0);
    }

  void main() {
      vec2 pos = gl_FragCoord.xy;
      vec2 texel = 1.0 / u_resolution;

      if (!insideBoundary(pos)) {
          vec4 current = texture(u_current, pos * texel);
          vec2 reflectedField = current.xy * ${config.boundaryReflectivity};
          fragColor = vec4(reflectedField, 0.0, 1.0);
          return;
      }

      vec4 center = texture(u_current, pos * texel);
      vec4 old = texture(u_previous, pos * texel);
      vec4 lens = texture(u_lens, pos * texel);

      vec2 laplacianRaw = ${config.use9PointStencil ? "nine" : "five"}PointLaplacian(pos, texel, center);
      vec2 lensedLaplacian = vec2(
          laplacianRaw.x * lens.x - laplacianRaw.y * lens.y,
          laplacianRaw.x * lens.y + laplacianRaw.y * lens.x
      );

      if (u_updateTarget == 0) {
        // Fundamental field update
        float amp2 = dot(center.xy, center.xy);
        float kerrSaturationFactor = 1.0 / (1.0 + amp2/u_kerr_Isat);
        float n2_effective = u_chi * u_chi * u_chi_ratio;
        float kerrStrength = n2_effective * amp2 * kerrSaturationFactor;
        float local_c2 = u_c * u_c * (1.0 + kerrStrength);

        vec2 new_val = (2.0 * center.xy - old.xy + local_c2 * u_dt * u_dt * lensedLaplacian) * u_damping;

        // Add pulsed injection
        if (u_frameCount % u_pulseInterval == 0 && u_pulseInterval > 0) {
              vec4 lens = texture(u_lens, pos * texel);

              // Gaussian beam parameters
              float w0 = 3.0;
              vec2 center = u_resolution * 0.5;
              vec2 rel = pos - center;
              float r2 = dot(rel, rel);

              // Calculate pulse amplitude
              float amplitude = 1. * exp(-r2/(w0*w0));
              float phase = 0.0;

              // Create and lens the new field
              vec2 pulseField = amplitude * vec2(cos(phase), sin(phase));
              vec2 lensedPulse = vec2(
                  pulseField.x * lens.x - pulseField.y * lens.y,
                  pulseField.x * lens.y + pulseField.y * lens.x
              );

              new_val += lensedPulse;
        }

        // Calculate metrics for fundamental field
        vec2 phaseGrad = calculatePhaseGradients(pos, texel, center);

        fragColor = vec4(
            new_val,                // xy: updated field
            length(phaseGrad),      // z: phase gradient magnitude
            kerrStrength           // w: Kerr nonlinearity strength
        );
      } else {
          // SHG field update
          float shg_amp2 = dot(center.xy, center.xy);
          float kerrSaturationFactor = 1.0 / (1.0 + shg_amp2/u_kerr_Isat);
          float n2_effective = u_chi * u_chi * u_chi_ratio;
          float kerrStrength = n2_effective * shg_amp2 * kerrSaturationFactor;
          float local_c2 = u_c * u_c * (1.0 + kerrStrength);

          vec4 fund = texture(u_fundamental, pos * texel);
          float fund_amp2 = dot(fund.xy, fund.xy);
          float shgSaturationFactor = 1.0 / (1.0 + fund_amp2/u_shg_Isat);

          vec4 phaseMismatchTerm = calculatePhaseMismatchTerm(pos, texel);
          vec2 sourceTerm = vec2(
              u_chi * fund_amp2 * shgSaturationFactor * phaseMismatchTerm.x,
              u_chi * fund_amp2 * shgSaturationFactor * phaseMismatchTerm.y
          );

          // Calculate phase matching efficiency
          float localPhaseMatch = dot(phaseMismatchTerm.xy, sourceTerm) /
                                 (length(phaseMismatchTerm.xy) * length(sourceTerm) + 1e-10);

          // Calculate χ(2)/χ(3) ratio
          float chiRatio = (u_chi * fund_amp2 * shgSaturationFactor) /
                          (n2_effective * shg_amp2 * kerrSaturationFactor + 1e-10);

          vec2 new_val = (2.0 * center.xy - old.xy +
                        local_c2 * u_dt * u_dt * lensedLaplacian +
                        u_dt * u_dt * sourceTerm) * u_damping;

          fragColor = vec4(
              new_val,           // xy: updated field
              localPhaseMatch,   // z: phase matching efficiency
              chiRatio          // w: χ(2)/χ(3) ratio
          );
      }
  }`;
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
