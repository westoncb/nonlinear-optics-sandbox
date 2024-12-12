// Shared vertex shader for all programs
export const vertexShaderSource = `#version 300 es
in vec2 position;
out vec2 uv;
void main() {
    uv = position * 0.5 + 0.5;
    gl_Position = vec4(position, 0.0, 1.0);
}`;

// Main simulation shader
export const simulationShaderSource = `#version 300 es
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
uniform float u_n2;
uniform float u_Isat;
uniform float u_chi;
uniform vec2 u_resolution;
uniform float u_boundaryR0;
uniform float u_boundaryAlpha;
uniform float u_boundaryM;
uniform int u_updateTarget;

out vec4 fragColor;

bool insideBoundary(vec2 pos) {
    vec2 center = u_resolution * 0.5;
    vec2 rel = pos - center;
    float r = length(rel);
    float theta = atan(rel.y, rel.x);
    float boundaryRadius = u_boundaryR0 * (1.0 + u_boundaryAlpha * cos(u_boundaryM * theta));
    return r <= boundaryRadius;
}

void main() {
    vec2 pos = gl_FragCoord.xy;
    vec2 texel = 1.0 / u_resolution;

    if (!insideBoundary(pos)) {
        vec4 current = texture(u_current, pos * texel);
        fragColor = vec4(-current.xy, 0.0, 1.0);
        return;
    }

    vec4 center = texture(u_current, pos * texel);
    vec4 left = texture(u_current, (pos + vec2(-1.0, 0.0)) * texel);
    vec4 right = texture(u_current, (pos + vec2(1.0, 0.0)) * texel);
    vec4 up = texture(u_current, (pos + vec2(0.0, -1.0)) * texel);
    vec4 down = texture(u_current, (pos + vec2(0.0, 1.0)) * texel);
    vec4 old = texture(u_previous, pos * texel);
    vec4 lens = texture(u_lens, pos * texel);

    vec2 laplacianRaw = (left.xy + right.xy + up.xy + down.xy - 4.0 * center.xy) / (u_dx * u_dx);

    vec2 lensedLaplacian = vec2(
        laplacianRaw.x * lens.x - laplacianRaw.y * lens.y,
        laplacianRaw.x * lens.y + laplacianRaw.y * lens.x
    );

    if (u_updateTarget == 0) {
        // Fundamental field update
        float amp2 = dot(center.xy, center.xy);
        float local_c2 = u_c * u_c * (1.0 + (u_n2 * amp2) / (1.0 + amp2/u_Isat));

        vec2 new_val = (2.0 * center.xy - old.xy + local_c2 * u_dt * u_dt * lensedLaplacian) * u_damping;
        fragColor = vec4(new_val, 0.0, 1.0);
    } else {
        // SHG field update
        float shg_amp2 = dot(center.xy, center.xy);
        float local_c2 = u_c * u_c * (1.0 + (u_n2 * shg_amp2) / (1.0 + shg_amp2/u_Isat));

        vec4 fund = texture(u_fundamental, pos * texel);
        float fund_amp2 = dot(fund.xy, fund.xy);
        float saturationFactor = 1.0 / (1.0 + fund_amp2/u_Isat);
        vec2 sourceTerm = vec2(u_chi * fund_amp2 * saturationFactor, 0.0);

        vec2 new_val = (2.0 * center.xy - old.xy +
                      local_c2 * u_dt * u_dt * lensedLaplacian +
                      u_dt * u_dt * sourceTerm) * u_damping;
        fragColor = vec4(new_val, 0.0, 1.0);
    }
}`;

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
