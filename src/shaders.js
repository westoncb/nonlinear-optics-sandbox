// Shared vertex shader for all programs
export const vertexShaderSource = `#version 300 es
in vec2 position;
out vec2 uv;
void main() {
    uv = position * 0.5 + 0.5;
    gl_Position = vec4(position, 0.0, 1.0);
}`;

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
