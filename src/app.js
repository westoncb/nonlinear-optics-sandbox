import { Config } from "./config.js";
import { Simulation } from "./simulation.js";
import { webgl, stats, scaleDataForDisplay } from "./util.js";
import {
  vertexShaderSource,
  simulationShaderSource,
  fieldDisplayShaderSource,
  lensDisplayShaderSource,
} from "./shaders.js";

export class App {
  constructor() {
    this.config = new Config();
    this.simulation = new Simulation(this.config);

    // WebGL contexts
    this.mainGL = null;
    this.lensGL = null;
    this.shgGL = null;

    // Programs
    this.programs = {
      simulation: null,
      fieldDisplay: null,
      lensDisplay: null,
      shgDisplay: null,
    };

    // Textures and framebuffers
    this.fundamentalTextures = [];
    this.fundamentalFramebuffers = [];
    this.shgTextures = [];
    this.shgFramebuffers = [];
    this.lensTexture = null;

    // Additional display textures
    this.shgDisplayTexture = null;
    this.lensDisplayTexture = null;

    // Uniform locations
    this.uniformLocations = {
      simulation: {},
      fieldDisplay: {},
      lensDisplay: {},
      shgDisplay: {},
    };

    // Animation state
    this.frameCount = 0;
    this.current = 0;
    this.previous = 2;
    this.next = 1;

    // Readback buffers
    this.readbackBuffer = null;
    this.shgDisplayBuffer = null;
    this.lensDisplayBuffer = null;

    // Lens statistics
    this.lensStatistics = {
      displayMin: 0,
      displayMax: 1,
    };

    // Quad buffers for each context
    this.quadBuffers = {
      main: null,
      lens: null,
      shg: null,
    };
  }

  async initialize() {
    try {
      await this.initWebGL();
      await this.initShaders();
      await this.initBuffers();
      await this.initTextures();
      this.initEvents();
      console.log("Initialization complete");
      return true;
    } catch (error) {
      console.error("Initialization failed:", error);
      document.getElementById("statusDisplay").textContent =
        "Initialization failed: " + error.message;
      throw error;
    }
  }

  async initWebGL() {
    // Initialize main canvas
    const mainCanvas = document.getElementById("waveCanvas");
    this.mainGL = mainCanvas.getContext("webgl2");
    if (!this.mainGL) throw new Error("WebGL 2 not available for main canvas");

    // Initialize lens canvas
    const lensCanvas = document.getElementById("lensCanvas");
    this.lensGL = lensCanvas.getContext("webgl2");
    if (!this.lensGL) throw new Error("WebGL 2 not available for lens canvas");

    // Initialize SHG canvas
    const shgCanvas = document.getElementById("shgCanvas");
    this.shgGL = shgCanvas.getContext("webgl2");
    if (!this.shgGL) throw new Error("WebGL 2 not available for SHG canvas");

    // Check for required extensions
    [this.mainGL, this.lensGL, this.shgGL].forEach((gl) => {
      const ext = gl.getExtension("EXT_color_buffer_float");
      if (!ext) throw new Error("EXT_color_buffer_float not available");
    });

    this.handleResize();
  }

  async initShaders() {
    // Create programs for each context
    this.programs.simulation = webgl.createProgram(
      this.mainGL,
      vertexShaderSource,
      simulationShaderSource,
    );

    this.programs.fieldDisplay = webgl.createProgram(
      this.mainGL,
      vertexShaderSource,
      fieldDisplayShaderSource,
    );

    this.programs.shgDisplay = webgl.createProgram(
      this.shgGL,
      vertexShaderSource,
      fieldDisplayShaderSource,
    );

    this.programs.lensDisplay = webgl.createProgram(
      this.lensGL,
      vertexShaderSource,
      lensDisplayShaderSource,
    );

    this.cacheUniformLocations();
  }

  cacheUniformLocations() {
    const simUniforms = [
      "u_dt",
      "u_dx",
      "u_damping",
      "u_c",
      "u_n2",
      "u_Isat",
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
      const location = this.mainGL.getUniformLocation(
        this.programs.simulation,
        name,
      );
      if (location === null) {
        console.warn(`Couldn't find simulation uniform: ${name}`);
      }
      this.uniformLocations.simulation[name] = location;
    });

    const fieldUniforms = ["u_field"];

    // Field display (main)
    this.mainGL.useProgram(this.programs.fieldDisplay);
    fieldUniforms.forEach((name) => {
      const location = this.mainGL.getUniformLocation(
        this.programs.fieldDisplay,
        name,
      );
      if (location === null) {
        console.warn(`Couldn't find field display uniform: ${name}`);
      }
      this.uniformLocations.fieldDisplay[name] = location;
    });

    // SHG display
    this.shgGL.useProgram(this.programs.shgDisplay);
    fieldUniforms.forEach((name) => {
      const location = this.shgGL.getUniformLocation(
        this.programs.shgDisplay,
        name,
      );
      if (location === null) {
        console.warn(`Couldn't find SHG display uniform: ${name}`);
      }
      this.uniformLocations.shgDisplay[name] = location;
    });

    const lensUniforms = [
      "u_field",
      "u_lensDisplayMin",
      "u_lensDisplayMax",
      "u_lensRadius",
      "u_resolution",
    ];

    this.lensGL.useProgram(this.programs.lensDisplay);
    lensUniforms.forEach((name) => {
      const location = this.lensGL.getUniformLocation(
        this.programs.lensDisplay,
        name,
      );
      if (location === null) {
        console.warn(`Couldn't find lens display uniform: ${name}`);
      }
      this.uniformLocations.lensDisplay[name] = location;
    });
  }

  async initBuffers() {
    // Create quad buffers for each context
    this.quadBuffers.main = this.createQuadBuffer(this.mainGL);
    this.quadBuffers.lens = this.createQuadBuffer(this.lensGL);
    this.quadBuffers.shg = this.createQuadBuffer(this.shgGL);

    // Set up attributes for each program now that we have quad buffers
    this.setupProgramAttributes(
      this.mainGL,
      this.programs.simulation,
      this.quadBuffers.main,
    );
    this.setupProgramAttributes(
      this.mainGL,
      this.programs.fieldDisplay,
      this.quadBuffers.main,
    );
    this.setupProgramAttributes(
      this.shgGL,
      this.programs.shgDisplay,
      this.quadBuffers.shg,
    );
    this.setupProgramAttributes(
      this.lensGL,
      this.programs.lensDisplay,
      this.quadBuffers.lens,
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

  async initTextures() {
    // Main context textures (for simulation)
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

    this.fundamentalFramebuffers = this.fundamentalTextures.map((texture) =>
      webgl.createFramebuffer(this.mainGL, texture),
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

    this.shgFramebuffers = this.shgTextures.map((texture) =>
      webgl.createFramebuffer(this.mainGL, texture),
    );

    this.lensTexture = webgl.createTexture(
      this.mainGL,
      this.config.gridSize,
      this.config.gridSize,
    );

    // Display textures
    this.shgDisplayTexture = webgl.createTexture(this.shgGL, 200, 200);
    this.lensDisplayTexture = webgl.createTexture(this.lensGL, 200, 200);

    // Initialize readback buffers
    this.readbackBuffer = new Float32Array(
      this.config.gridSize * this.config.gridSize * 4,
    );
    this.shgDisplayBuffer = new Float32Array(200 * 200 * 4);
    this.lensDisplayBuffer = new Float32Array(200 * 200 * 4);

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

    // Initialize fundamental field with gaussian pulse
    for (let j = -3; j <= 3; j++) {
      for (let i = -3; i <= 3; i++) {
        const x = centerX + i;
        const y = centerY + j;
        const idx = (y * this.config.gridSize + x) * 4;
        fundamentalData[idx] = Math.exp(-(i * i + j * j) / 4) * 15;
        fundamentalData[idx + 1] = 0.0001;
      }
    }

    // Upload initial data
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
    const lensData = this.simulation.getLensData(
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

    this.simulation.updateLens(
      fundamentalData,
      shgData,
      this.config.gridSize,
      this.config.gridSize,
    );

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
    // Update simulation
    this.mainGL.useProgram(this.programs.simulation);

    // Fundamental field update
    this.mainGL.bindFramebuffer(
      this.mainGL.FRAMEBUFFER,
      this.fundamentalFramebuffers[this.next],
    );
    this.setSimulationUniforms(0);
    this.mainGL.drawArrays(this.mainGL.TRIANGLE_STRIP, 0, 4);

    // SHG field update
    this.mainGL.bindFramebuffer(
      this.mainGL.FRAMEBUFFER,
      this.shgFramebuffers[this.next],
    );
    this.setSimulationUniforms(1);
    this.mainGL.drawArrays(this.mainGL.TRIANGLE_STRIP, 0, 4);

    // Render displays
    this.renderMainDisplay();
    this.renderSHGDisplay();
    this.renderLensDisplay();

    // Update lens if needed
    if (this.frameCount % this.config.optimizationInterval === 0) {
      this.updateLens();
    }

    // Cycle buffers
    [this.previous, this.current, this.next] = [
      this.current,
      this.next,
      this.previous,
    ];

    this.frameCount++;
  }

  renderMainDisplay() {
    const gl = this.mainGL;
    gl.useProgram(this.programs.fieldDisplay);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.fundamentalTextures[this.next]);
    gl.uniform1i(this.uniformLocations.fieldDisplay.u_field, 0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  renderSHGDisplay() {
    // Bind the framebuffer that contains the SHG data
    const gl = this.mainGL;
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.shgFramebuffers[this.next]);
    gl.readPixels(
      0,
      0,
      this.config.gridSize,
      this.config.gridSize,
      gl.RGBA,
      gl.FLOAT,
      this.readbackBuffer,
    );

    // Scale down the data for display
    scaleDataForDisplay(
      this.readbackBuffer,
      this.config.gridSize,
      this.shgDisplayBuffer,
      200,
    );

    // Update and render the SHG display texture
    const shgGL = this.shgGL;
    shgGL.useProgram(this.programs.shgDisplay);
    shgGL.activeTexture(shgGL.TEXTURE0);
    shgGL.bindTexture(shgGL.TEXTURE_2D, this.shgDisplayTexture);
    shgGL.texImage2D(
      shgGL.TEXTURE_2D,
      0,
      shgGL.RGBA32F,
      200,
      200,
      0,
      shgGL.RGBA,
      shgGL.FLOAT,
      this.shgDisplayBuffer,
    );
    shgGL.uniform1i(this.uniformLocations.shgDisplay.u_field, 0);
    shgGL.bindFramebuffer(shgGL.FRAMEBUFFER, null);
    shgGL.drawArrays(shgGL.TRIANGLE_STRIP, 0, 4);
  }

  renderLensDisplay() {
    // Scale the lens data for display (CPU -> lensDisplayBuffer)
    scaleDataForDisplay(
      this.simulation.getLensData(this.config.gridSize, this.config.gridSize),
      this.config.gridSize,
      this.lensDisplayBuffer,
      200,
    );

    const gl = this.lensGL;
    gl.useProgram(this.programs.lensDisplay);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.lensDisplayTexture);
    gl.texImage2D(
      gl.TEXTURE_2D,
      0,
      gl.RGBA32F,
      200,
      200,
      0,
      gl.RGBA,
      gl.FLOAT,
      this.lensDisplayBuffer,
    );
    gl.uniform1i(this.uniformLocations.lensDisplay.u_field, 0);
    this.setLensDisplayUniforms();
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
  }

  setSimulationUniforms(updateTarget) {
    const gl = this.mainGL;
    const uniforms = this.uniformLocations.simulation;

    gl.uniform1f(uniforms.u_dt, this.config.dt);
    gl.uniform1f(uniforms.u_dx, this.config.dx);
    gl.uniform1f(uniforms.u_damping, this.config.damping);
    gl.uniform1f(uniforms.u_c, this.config.c);
    gl.uniform1f(uniforms.u_n2, this.config.n2);
    gl.uniform1f(uniforms.u_Isat, this.config.Isat);
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

    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, this.lensTexture);
    gl.uniform1i(uniforms.u_lens, 2);

    if (updateTarget === 1) {
      gl.activeTexture(gl.TEXTURE3);
      gl.bindTexture(gl.TEXTURE_2D, this.fundamentalTextures[this.next]);
      gl.uniform1i(uniforms.u_fundamental, 3);
    }
  }

  setLensDisplayUniforms() {
    const gl = this.lensGL;
    const uniforms = this.uniformLocations.lensDisplay;

    gl.uniform1f(uniforms.u_lensDisplayMin, this.lensStatistics.displayMin);
    gl.uniform1f(uniforms.u_lensDisplayMax, this.lensStatistics.displayMax);
    gl.uniform1f(uniforms.u_lensRadius, this.config.lensRadius);
    gl.uniform2f(uniforms.u_resolution, 200, 200);
  }

  initEvents() {
    window.addEventListener("resize", () => this.handleResize());
  }

  handleResize() {
    const mainCanvas = document.getElementById("waveCanvas");
    const pixelRatio = window.devicePixelRatio || 1;
    mainCanvas.width = window.innerWidth * pixelRatio;
    mainCanvas.height = (window.innerHeight - 200) * pixelRatio;
    mainCanvas.style.width = window.innerWidth + "px";
    mainCanvas.style.height = window.innerHeight - 200 + "px";
    this.mainGL.viewport(0, 0, mainCanvas.width, mainCanvas.height);
  }

  animate() {
    this.render();
    requestAnimationFrame(() => this.animate());
  }

  start() {
    this.animate();
  }
}

// Start the application
const app = new App();
app.initialize().then(() => app.start());
