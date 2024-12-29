// Color conversion
export function hsv2rgb(h, s, v) {
  const K = [1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0];
  const p = Math.abs(((h + K[0]) % 6) - K[3]);
  const q = Math.abs(((h + K[1]) % 6) - K[3]);
  const t = Math.abs(((h + K[2]) % 6) - K[3]);
  return [
    v * (1 - s * Math.max(0, Math.min(p, 1, 4 - p))),
    v * (1 - s * Math.max(0, Math.min(q, 1, 4 - q))),
    v * (1 - s * Math.max(0, Math.min(t, 1, 4 - t))),
  ].map((x) => Math.floor(x * 255));
}

// WebGL helpers
export const webgl = {
  createShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      const error = new Error(
        "Shader compile error: " + gl.getShaderInfoLog(shader),
      );
      gl.deleteShader(shader);
      throw error;
    }
    return shader;
  },

  createProgram(gl, vertexSource, fragmentSource) {
    const program = gl.createProgram();
    const vertexShader = this.createShader(gl, gl.VERTEX_SHADER, vertexSource);
    const fragmentShader = this.createShader(
      gl,
      gl.FRAGMENT_SHADER,
      fragmentSource,
    );
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
      throw new Error("Program link error: " + gl.getProgramInfoLog(program));
    }
    return program;
  },

  createTexture(gl, width, height) {
    const texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.RGBA32F,
      width,
      height,
      0,
      gl.RGBA,
      gl.FLOAT,
      null,
    );
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return texture;
  },

  createFramebuffer(gl, texture) {
    const fb = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(
      gl.FRAMEBUFFER,
      gl.COLOR_ATTACHMENT0,
      gl.TEXTURE_2D,
      texture,
      0,
    );
    const status = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (status !== gl.FRAMEBUFFER_COMPLETE) {
      throw new Error("Framebuffer is incomplete: " + status);
    }
    return fb;
  },
};

// Statistics helpers
export const stats = {
  calculateQuartiles(values) {
    if (values.length === 0) return { q25: 0, q75: 0 };
    const sorted = [...values].sort((a, b) => a - b);
    return {
      q25: sorted[Math.floor(values.length * 0.25)] || 0,
      q75: sorted[Math.floor(values.length * 0.75)] || 1,
    };
  },

  calculateDisplayRange(values) {
    const { q25, q75 } = this.calculateQuartiles(values);
    const iqr = q75 - q25;
    return {
      min: Math.max(0, q25 - 1.5 * iqr),
      max: q75 + 1.5 * iqr,
    };
  },
};

// Polar coordinate helpers
export const polar = {
  getZoneAndSector(x, y, centerX, centerY, radius, numZones, numSectors) {
    const dx = x - centerX;
    const dy = y - centerY;
    const r = Math.sqrt(dx * dx + dy * dy);

    if (r > radius) {
      return null; // This early return is correct
    } else {
      const zoneIndex = Math.floor((r / radius) * numZones);
      let theta = Math.atan2(dy, dx);
      let sectorIndex =
        Math.floor((theta + Math.PI) / ((2 * Math.PI) / numSectors)) %
        numSectors;

      return {
        zone: Math.min(numZones - 1, Math.max(0, zoneIndex)),
        sector: sectorIndex,
      };
    }
  },
};

export const mulberry32 = (seed) => {
  return function () {
    let t = (seed += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
};

// Gaussian random number generator using Box-Muller
export const gaussianRandom = (prng) => {
  const u1 = prng();
  const u2 = prng();
  const radius = Math.sqrt(-2.0 * Math.log(u1));
  const theta = 2.0 * Math.PI * u2;
  return radius * Math.cos(theta);
};

export const scaleDataForDisplay = (
  sourceData,
  sourceSize,
  targetBuffer,
  targetSize,
) => {
  const scale = sourceSize / targetSize;

  for (let y = 0; y < targetSize; y++) {
    for (let x = 0; x < targetSize; x++) {
      const sourceX = Math.floor(x * scale);
      const sourceY = Math.floor(y * scale);

      const sourceIdx = (sourceY * sourceSize + sourceX) * 4;
      const targetIdx = (y * targetSize + x) * 4;

      targetBuffer[targetIdx] = sourceData[sourceIdx];
      targetBuffer[targetIdx + 1] = sourceData[sourceIdx + 1];
      targetBuffer[targetIdx + 2] = sourceData[sourceIdx + 2];
      targetBuffer[targetIdx + 3] = sourceData[sourceIdx + 3];
    }
  }
};
