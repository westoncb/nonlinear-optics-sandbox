import { LensOptimizer } from "./lens-optimizer.js";
import { webgl, stats, scaleDataForDisplay } from "./util.js";
import {
  vertexShaderSource,
  getSimulationShaderSource,
  displayShaderSource,
} from "./shaders.js";

const DisplayMode = {
  FUNDAMENTAL: 0,
  SHG: 1,
  LENS: 2,
};

export class App {
  constructor(config) {
    this.config = config.data;
    this.lensOptimizer = new LensOptimizer(this.config);
    this.isRunning = false;

    console.log("CONFIG", this.config);

    // GL contexts
    // mainGL for simulation AND primary display:
    const mainCanvas = document.getElementById("primaryCanvas");
    this.mainGL = mainCanvas.getContext("webgl2");
    if (!this.mainGL)
      throw new Error("WebGL 2 not available for primaryCanvas");

    const preview1Canvas = document.getElementById("preview1Canvas");
    this.preview1GL = preview1Canvas.getContext("webgl2");
    if (!this.preview1GL)
      throw new Error("WebGL 2 not available for preview1Canvas");

    const preview2Canvas = document.getElementById("preview2Canvas");
    this.preview2GL = preview2Canvas.getContext("webgl2");
    if (!this.preview2GL)
      throw new Error("WebGL 2 not available for preview2Canvas");

    [this.mainGL, this.preview1GL, this.preview2GL].forEach((gl) => {
      const ext = gl.getExtension("EXT_color_buffer_float");
      if (!ext) throw new Error("EXT_color_buffer_float not available");
    });

    // Programs
    // Simulation in mainGL
    this.programs = {
      simulation: null,
      primaryDisplay: null, // display in mainGL
      preview1Display: null, // display in preview1GL
      preview2Display: null, // display in preview2GL
    };

    // Simulation textures and framebuffers
    this.fundamentalTextures = [];
    this.fundamentalFramebuffers = [];
    this.shgTextures = [];
    this.shgFramebuffers = [];
    this.lensTexture = null;

    // Display textures (reused each frame)
    // For primary (mainGL), we can render directly from simulation textures—no extra texture needed.
    this.preview1Texture = null;
    this.preview2Texture = null;

    // Uniform locations
    this.uniformLocations = {
      simulation: {},
      primaryDisplay: {},
      preview1Display: {},
      preview2Display: {},
    };

    // Animation state
    this.frameCount = 0;
    this.current = 0;
    this.previous = 2;
    this.next = 1;

    // Readback buffers
    this.readbackBuffer = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );
    this.shgDisplayBuffer = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );
    this.lensDisplayBuffer = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );

    // Lens stats
    this.lensStatistics = {
      displayMin: 0,
      displayMax: 1,
    };

    // Quad buffers
    this.quadBuffers = {
      main: null,
      preview1: null,
      preview2: null,
    };

    this.displayModes = {
      primary: DisplayMode.FUNDAMENTAL,
      preview1: DisplayMode.LENS,
      preview2: DisplayMode.SHG,
    };
  }

  async initialize() {
    try {
      this.initShaders();
      this.initBuffers();
      this.initTextures();
      console.log("Initialization complete");
      return true;
    } catch (error) {
      console.error("Initialization failed:", error);
      document.getElementById("statusDisplay").textContent =
        "Initialization failed: " + error.message;
      throw error;
    }
  }

  initShaders() {
    // Simulation
    this.programs.simulation = webgl.createProgram(
      this.mainGL,
      vertexShaderSource,
      getSimulationShaderSource(this.config),
    );

    // Display
    this.programs.primaryDisplay = webgl.createProgram(
      this.mainGL,
      vertexShaderSource,
      displayShaderSource,
    );

    this.programs.preview1Display = webgl.createProgram(
      this.preview1GL,
      vertexShaderSource,
      displayShaderSource,
    );

    this.programs.preview2Display = webgl.createProgram(
      this.preview2GL,
      vertexShaderSource,
      displayShaderSource,
    );

    this.cacheUniformLocations();
  }

  cacheUniformLocations() {
    // Simulation uniforms
    const simUniforms = [
      "u_dt",
      "u_dx",
      "u_damping",
      "u_c",
      "u_chi_ratio",
      "u_shg_Isat",
      "u_kerr_Isat",
      "u_chi",
      "u_resolution",
      "u_boundaryR0",
      "u_boundaryAlpha",
      "u_boundaryM",
      "u_updateTarget",
      "u_current",
      "u_previous",
      "u_lens",
      "u_fundamental",
    ];
    this.mainGL.useProgram(this.programs.simulation);
    simUniforms.forEach((name) => {
      const loc = this.mainGL.getUniformLocation(
        this.programs.simulation,
        name,
      );
      this.uniformLocations.simulation[name] = loc;
    });

    const displayUniforms = [
      "u_field",
      "u_displayMode",
      "u_lensDisplayMin",
      "u_lensDisplayMax",
      "u_lensRadius",
      "u_resolution",
    ];

    // Primary (mainGL)
    this.mainGL.useProgram(this.programs.primaryDisplay);
    displayUniforms.forEach((name) => {
      this.uniformLocations.primaryDisplay[name] =
        this.mainGL.getUniformLocation(this.programs.primaryDisplay, name);
    });

    // Preview1
    this.preview1GL.useProgram(this.programs.preview1Display);
    displayUniforms.forEach((name) => {
      this.uniformLocations.preview1Display[name] =
        this.preview1GL.getUniformLocation(this.programs.preview1Display, name);
    });

    // Preview2
    this.preview2GL.useProgram(this.programs.preview2Display);
    displayUniforms.forEach((name) => {
      this.uniformLocations.preview2Display[name] =
        this.preview2GL.getUniformLocation(this.programs.preview2Display, name);
    });
  }

  initBuffers() {
    this.quadBuffers.main = this.createQuadBuffer(this.mainGL);
    this.quadBuffers.preview1 = this.createQuadBuffer(this.preview1GL);
    this.quadBuffers.preview2 = this.createQuadBuffer(this.preview2GL);

    // Simulation attributes
    this.setupProgramAttributes(
      this.mainGL,
      this.programs.simulation,
      this.quadBuffers.main,
    );

    // Display attributes
    this.setupProgramAttributes(
      this.mainGL,
      this.programs.primaryDisplay,
      this.quadBuffers.main,
    );
    this.setupProgramAttributes(
      this.preview1GL,
      this.programs.preview1Display,
      this.quadBuffers.preview1,
    );
    this.setupProgramAttributes(
      this.preview2GL,
      this.programs.preview2Display,
      this.quadBuffers.preview2,
    );
  }

  createQuadBuffer(gl) {
    const quadBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quadBuffer);
    gl.bufferData(
      gl.ARRAY_BUFFER,
      new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]),
      gl.STATIC_DRAW,
    );
    return quadBuffer;
  }

  setupProgramAttributes(gl, program, quadBuffer) {
    gl.useProgram(program);
    gl.bindBuffer(gl.ARRAY_BUFFER, quadBuffer);
    const posLoc = gl.getAttribLocation(program, "position");
    if (posLoc >= 0) {
      gl.enableVertexAttribArray(posLoc);
      gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 0, 0);
    }
  }

  initTextures() {
    // Simulation textures
    this.fundamentalTextures = [
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
    ];
    this.fundamentalFramebuffers = this.fundamentalTextures.map((t) =>
      webgl.createFramebuffer(this.mainGL, t),
    );

    this.shgTextures = [
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
      webgl.createTexture(
        this.mainGL,
        this.config.gridSize,
        this.config.gridSize,
      ),
    ];
    this.shgFramebuffers = this.shgTextures.map((t) =>
      webgl.createFramebuffer(this.mainGL, t),
    );

    this.lensTexture = webgl.createTexture(
      this.mainGL,
      this.config.gridSize,
      this.config.gridSize,
    );

    // Display textures for previews
    // We'll allocate them once and reuse them by calling texSubImage2D each frame
    this.preview1Texture = webgl.createTexture(
      this.preview1GL,
      this.config.gridSize,
      this.config.gridSize,
    );
    this.preview2Texture = webgl.createTexture(
      this.preview2GL,
      this.config.gridSize,
      this.config.gridSize,
    );

    this.initializeFields();
  }

  initializeFields() {
    const fundamentalData = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );
    const zeroData = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );

    const centerX = Math.floor(this.config.gridSize / 2);
    const centerY = Math.floor(this.config.gridSize / 2);

    for (let j = -3; j <= 3; j++) {
      for (let i = -3; i <= 3; i++) {
        const x = centerX + i;
        const y = centerY + j;
        const idx = (y * this.config.gridSize + x) * 4;
        fundamentalData[idx] =
          Math.exp(-(i * i + j * j) / 4) * this.config.initialPulseAmplitude;
        fundamentalData[idx + 1] = this.config.initialPulsePhaseShift;
      }
    }

    // Upload initial fields
    this.fundamentalTextures.forEach((texture) => {
      this.mainGL.bindTexture(this.mainGL.TEXTURE_2D, texture);
      this.mainGL.texImage2D(
        this.mainGL.TEXTURE_2D,
        0,
        this.mainGL.RGBA32F,
        this.config.gridSize,
        this.config.gridSize,
        0,
        this.mainGL.RGBA,
        this.mainGL.FLOAT,
        fundamentalData,
      );
    });

    this.shgTextures.forEach((texture) => {
      this.mainGL.bindTexture(this.mainGL.TEXTURE_2D, texture);
      this.mainGL.texImage2D(
        this.mainGL.TEXTURE_2D,
        0,
        this.mainGL.RGBA32F,
        this.config.gridSize,
        this.config.gridSize,
        0,
        this.mainGL.RGBA,
        this.mainGL.FLOAT,
        zeroData,
      );
    });
  }

  readFieldData(texture) {
    const gl = this.mainGL;
    const fb = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(
      gl.FRAMEBUFFER,
      gl.COLOR_ATTACHMENT0,
      gl.TEXTURE_2D,
      texture,
      0,
    );
    gl.readPixels(
      0,
      0,
      this.config.gridSize,
      this.config.gridSize,
      gl.RGBA,
      gl.FLOAT,
      this.readbackBuffer,
    );
    gl.deleteFramebuffer(fb);
    return this.readbackBuffer;
  }

  updateLensTexture() {
    const lensData = this.lensOptimizer.getLensData(
      this.config.gridSize,
      this.config.gridSize,
    );
    this.mainGL.bindTexture(this.mainGL.TEXTURE_2D, this.lensTexture);
    this.mainGL.texImage2D(
      this.mainGL.TEXTURE_2D,
      0,
      this.mainGL.RGBA32F,
      this.config.gridSize,
      this.config.gridSize,
      0,
      this.mainGL.RGBA,
      this.mainGL.FLOAT,
      lensData,
    );
  }

  updateLens() {
    const fundamentalData = this.readFieldData(
      this.fundamentalTextures[this.current],
    );
    const shgData = this.readFieldData(this.shgTextures[this.current]);

    if (!this.config.disableAdaptation) {
      this.lensOptimizer.updateLens(
        fundamentalData,
        shgData,
        this.config.gridSize,
        this.config.gridSize,
      );
    }

    this.updateLensTexture();
    this.updateLensStatistics();
  }

  updateLensStatistics() {
    const lensData = this.readFieldData(this.lensTexture);
    let values = [];

    const centerX = Math.floor(this.config.gridSize / 2);
    const centerY = Math.floor(this.config.gridSize / 2);

    for (let y = 0; y < this.config.gridSize; y++) {
      for (let x = 0; x < this.config.gridSize; x++) {
        const dx = x - centerX;
        const dy = y - centerY;
        const r = Math.sqrt(dx * dx + dy * dy);
        if (r <= this.config.lensRadius) {
          const idx = (y * this.config.gridSize + x) * 4;
          const LR = lensData[idx];
          const LI = lensData[idx + 1];
          const amp = Math.sqrt(LR * LR + LI * LI);
          values.push(amp);
        }
      }
    }

    const range = stats.calculateDisplayRange(values);
    this.lensStatistics.displayMin = range.min;
    this.lensStatistics.displayMax = range.max;
  }

  render() {
    // Simulation step
    this.runSimulationStep();

    // Update lens occasionally
    if (this.frameCount % this.config.optimizationInterval === 0) {
      this.updateLens();
    }

    // Render primary (no CPU readback if FUNDAMENTAL)
    this.renderPrimary();

    // Render previews (may involve CPU readback)
    this.renderPreview(
      this.preview1GL,
      this.uniformLocations.preview1Display,
      this.displayModes.preview1,
      this.preview1Texture,
      this.preview1CanvasSize(),
    );
    this.renderPreview(
      this.preview2GL,
      this.uniformLocations.preview2Display,
      this.displayModes.preview2,
      this.preview2Texture,
      this.preview2CanvasSize(),
    );

    // Cycle buffers
    [this.previous, this.current, this.next] = [
      this.current,
      this.next,
      this.previous,
    ];
    this.frameCount++;
  }

  runSimulationStep() {
    const gl = this.mainGL;
    gl.useProgram(this.programs.simulation);

    // Fundamental update
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.fundamentalFramebuffers[this.next]);
    this.setSimulationUniforms(0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

    // SHG update
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.shgFramebuffers[this.next]);
    this.setSimulationUniforms(1);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  setSimulationUniforms(updateTarget) {
    const gl = this.mainGL;
    const uniforms = this.uniformLocations.simulation;

    gl.uniform1f(uniforms.u_dt, this.config.dt);
    gl.uniform1f(uniforms.u_dx, this.config.dx);
    gl.uniform1f(uniforms.u_damping, this.config.damping);
    gl.uniform1f(uniforms.u_c, this.config.c);
    gl.uniform1f(uniforms.u_chi_ratio, this.config.chi_ratio);
    gl.uniform1f(uniforms.u_shg_Isat, this.config.shg_Isat);
    gl.uniform1f(uniforms.u_kerr_Isat, this.config.kerr_Isat);
    gl.uniform1f(uniforms.u_chi, this.config.chi);
    gl.uniform2f(
      uniforms.u_resolution,
      this.config.gridSize,
      this.config.gridSize,
    );
    gl.uniform1f(uniforms.u_boundaryR0, this.config.boundaryR0);
    gl.uniform1f(uniforms.u_boundaryAlpha, this.config.boundaryAlpha);
    gl.uniform1f(uniforms.u_boundaryM, this.config.boundaryM);
    gl.uniform1i(uniforms.u_updateTarget, updateTarget);

    // current/previous textures
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(
      gl.TEXTURE_2D,
      updateTarget === 0
        ? this.fundamentalTextures[this.current]
        : this.shgTextures[this.current],
    );
    gl.uniform1i(uniforms.u_current, 0);

    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(
      gl.TEXTURE_2D,
      updateTarget === 0
        ? this.fundamentalTextures[this.previous]
        : this.shgTextures[this.previous],
    );
    gl.uniform1i(uniforms.u_previous, 1);

    // Lens
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, this.lensTexture);
    gl.uniform1i(uniforms.u_lens, 2);

    // Fundamental for SHG updates
    if (updateTarget === 1) {
      gl.activeTexture(gl.TEXTURE3);
      gl.bindTexture(gl.TEXTURE_2D, this.fundamentalTextures[this.next]);
      gl.uniform1i(uniforms.u_fundamental, 3);
    }
  }

  renderPrimary() {
    const gl = this.mainGL;
    gl.useProgram(this.programs.primaryDisplay);

    const uniforms = this.uniformLocations.primaryDisplay;
    const mode = this.displayModes.primary;
    gl.uniform1i(uniforms.u_displayMode, mode);

    if (mode === DisplayMode.FUNDAMENTAL) {
      // Directly use fundamentalTextures[this.current]
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.fundamentalTextures[this.current]);
      gl.uniform1i(uniforms.u_field, 0);

      // Set lens uniforms to default (not used)
      gl.uniform1f(uniforms.u_lensDisplayMin, 0.0);
      gl.uniform1f(uniforms.u_lensDisplayMax, 1.0);
      gl.uniform1f(uniforms.u_lensRadius, 0.0);

      // Full resolution
      gl.uniform2f(
        uniforms.u_resolution,
        this.config.gridSize,
        this.config.gridSize,
      );
    } else if (mode === DisplayMode.SHG) {
      // If we ever choose SHG for primary, we can similarily read directly from shgTextures[this.current]
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.shgTextures[this.current]);
      gl.uniform1i(uniforms.u_field, 0);

      gl.uniform1f(uniforms.u_lensDisplayMin, 0.0);
      gl.uniform1f(uniforms.u_lensDisplayMax, 1.0);
      gl.uniform1f(uniforms.u_lensRadius, 0.0);

      gl.uniform2f(
        uniforms.u_resolution,
        this.config.gridSize,
        this.config.gridSize,
      );
    } else if (mode === DisplayMode.LENS) {
      // For lens in primary, can directly use lensTexture
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.lensTexture);
      gl.uniform1i(uniforms.u_field, 0);

      gl.uniform1f(uniforms.u_lensDisplayMin, this.lensStatistics.displayMin);
      gl.uniform1f(uniforms.u_lensDisplayMax, this.lensStatistics.displayMax);
      gl.uniform1f(uniforms.u_lensRadius, this.config.lensRadius);
      gl.uniform2f(
        uniforms.u_resolution,
        this.config.gridSize,
        this.config.gridSize,
      );
    }

    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  preview1CanvasSize() {
    return this.config.gridSize;
  }

  preview2CanvasSize() {
    return this.config.gridSize;
  }

  renderPreview(gl, uniforms, mode, texture, previewSize) {
    gl.useProgram(
      this.programs[
        gl === this.preview1GL ? "preview1Display" : "preview2Display"
      ],
    );

    // For SHG and LENS modes (and even FUNDAMENTAL if chosen), we do CPU readback + scale
    // Then upload to the pre-created texture with texSubImage2D.

    let data = null;
    let lensSettings = { min: 0, max: 1, radius: 0 };
    const size = this.config.gridSize;

    if (mode === DisplayMode.FUNDAMENTAL) {
      // CPU readback for fundamental if displayed in preview:
      // This is less common, but we handle it similarly to SHG.
      this.mainGL.bindFramebuffer(
        this.mainGL.FRAMEBUFFER,
        this.fundamentalFramebuffers[this.current],
      );
      this.mainGL.readPixels(
        0,
        0,
        size,
        size,
        this.mainGL.RGBA,
        this.mainGL.FLOAT,
        this.readbackBuffer,
      );
      scaleDataForDisplay(
        this.readbackBuffer,
        size,
        this.shgDisplayBuffer,
        previewSize,
      );
      data = this.shgDisplayBuffer;
    } else if (mode === DisplayMode.SHG) {
      this.mainGL.bindFramebuffer(
        this.mainGL.FRAMEBUFFER,
        this.shgFramebuffers[this.current],
      );
      this.mainGL.readPixels(
        0,
        0,
        size,
        size,
        this.mainGL.RGBA,
        this.mainGL.FLOAT,
        this.readbackBuffer,
      );
      scaleDataForDisplay(
        this.readbackBuffer,
        size,
        this.shgDisplayBuffer,
        previewSize,
      );
      data = this.shgDisplayBuffer;
    } else if (mode === DisplayMode.LENS) {
      const lensData = this.readFieldData(this.lensTexture); // CPU readback from lens
      scaleDataForDisplay(lensData, size, this.lensDisplayBuffer, previewSize);
      data = this.lensDisplayBuffer;
      lensSettings = {
        min: this.lensStatistics.displayMin,
        max: this.lensStatistics.displayMax,
        radius: this.config.lensRadius,
      };
    }

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texSubImage2D(
      gl.TEXTURE_2D,
      0,
      0,
      0,
      previewSize,
      previewSize,
      gl.RGBA,
      gl.FLOAT,
      data,
    );

    gl.uniform1i(uniforms.u_field, 0);
    gl.uniform1i(uniforms.u_displayMode, mode);
    gl.uniform1f(uniforms.u_lensDisplayMin, lensSettings.min);
    gl.uniform1f(uniforms.u_lensDisplayMax, lensSettings.max);
    gl.uniform1f(uniforms.u_lensRadius, lensSettings.radius);
    gl.uniform2f(uniforms.u_resolution, previewSize, previewSize);

    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  setDisplayMode({ primary, preview1, preview2 }) {
    if (primary !== undefined) this.displayModes.primary = primary;
    if (preview1 !== undefined) this.displayModes.preview1 = preview1;
    if (preview2 !== undefined) this.displayModes.preview2 = preview2;
  }

  animate() {
    if (!this.isRunning) return;

    this.render();
    requestAnimationFrame(() => this.animate());
  }

  start() {
    this.isRunning = true;
    this.animate();
  }

  stop() {
    this.isRunning = false;
  }

  cleanup() {
    const glContexts = [
      {
        gl: this.mainGL,
        programs: [this.programs.simulation, this.programs.primaryDisplay],
        textures: [
          ...this.fundamentalTextures,
          ...this.shgTextures,
          this.lensTexture,
        ],
        buffers: [this.quadBuffers.main],
        framebuffers: [
          ...this.fundamentalFramebuffers,
          ...this.shgFramebuffers,
        ],
      },
      {
        gl: this.preview1GL,
        programs: [this.programs.preview1Display],
        textures: [this.preview1Texture],
        buffers: [this.quadBuffers.preview1],
      },
      {
        gl: this.preview2GL,
        programs: [this.programs.preview2Display],
        textures: [this.preview2Texture],
        buffers: [this.quadBuffers.preview2],
      },
    ];

    glContexts.forEach(
      ({ gl, programs, textures = [], buffers = [], framebuffers = [] }) => {
        if (!gl) return;

        // Delete programs
        programs.forEach((program) => {
          if (program) {
            // Detach and delete shaders first
            const attachedShaders = gl.getAttachedShaders(program);
            attachedShaders.forEach((shader) => {
              gl.detachShader(program, shader);
              gl.deleteShader(shader);
            });

            gl.deleteProgram(program);
          }
        });

        // Delete textures
        textures.forEach((texture) => {
          if (texture) gl.deleteTexture(texture);
        });

        // Delete buffers
        buffers.forEach((buffer) => {
          if (buffer) gl.deleteBuffer(buffer);
        });

        // Delete framebuffers
        framebuffers.forEach((framebuffer) => {
          if (framebuffer) gl.deleteFramebuffer(framebuffer);
        });

        // Unbind everything
        gl.bindBuffer(gl.ARRAY_BUFFER, null);
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
        gl.bindTexture(gl.TEXTURE_2D, null);
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.bindRenderbuffer(gl.RENDERBUFFER, null);
        gl.useProgram(null);
      },
    );

    // Clear references
    this.programs = {};
    this.quadBuffers = {};
    this.fundamentalTextures = [];
    this.fundamentalFramebuffers = [];
    this.shgTextures = [];
    this.shgFramebuffers = [];
    this.lensTexture = null;
    this.preview1Texture = null;
    this.preview2Texture = null;

    // Clear WebGL context references
    this.mainGL = null;
    this.preview1GL = null;
    this.preview2GL = null;
  }
}
