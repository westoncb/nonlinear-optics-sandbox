export const constants = {
  PREVIEW_SIZE: 200,
};

export class Config {
  constructor() {
    // Core simulation parameters
    this.dt = 0.35; // Slightly smaller timestep for stability
    this.dx = 1.0; // Keep spatial discretization simple
    this.gridSize = 512; // Larger domain for more complex pattern formation
    this.damping = 0.9999; // Slight damping to avoid unbounded growth but still allow some wave action
    this.c = 1.0; // Base wave speed

    // Nonlinear parameters
    this.n2 = 2.0; // Increased nonlinear index for stronger self-focusing
    this.Isat = 0.25; // Slightly lower saturation intensity to see nonlinear effects kick in earlier
    this.chi = 0.7; // Stronger second-harmonic generation coefficient to amplify nonlinear interactions

    // Boundary parameters
    this.boundaryAlpha = 0.4; // More pronounced boundary shaping
    this.boundaryM = 10; // Higher mode boundary modulation
    this.margin = 10; // A bit more margin for boundary conditions

    // Lens parameters
    this.lensRadius = 60; // Larger lens radius for more intricate focusing patterns
    this.fresnelZones = 32; // More zones for finer radial control
    this.numSectors = 1800; // Higher angular resolution to produce complex angular patterns
    this.hillPower = 1.5; // Slightly increased hill power for a more nuanced logistic weighting
    this.updateStrategy = "o1OptimalStrategy2"; // Using the PDE-informed strategy

    // Optimization parameters
    this.learningRate = 0.00004; // Slightly higher learning rate for more dynamic lens updates
    this.optimizationInterval = 1; // Keep frequent updates to lens modes

    // Initial pulse parameters
    this.initialPulseAmplitude = 30; // Strong initial field, but not too high
    this.initialPulsePhaseShift = 0.7; // Non-zero phase shift to create interesting initial conditions

    // Calculate derived values
    this.updateDerivedValues();

    // Initialize preset storage
    this.presets = new Map();
    this.loadPresets();
  }

  // Update any values that depend on other parameters
  updateDerivedValues() {
    this.boundaryR0 =
      (this.gridSize / 2 - this.margin) / (1 + this.boundaryAlpha);
  }

  // Save current configuration as a preset
  savePreset(name) {
    const preset = {
      dt: this.dt,
      dx: this.dx,
      gridSize: this.gridSize,
      damping: this.damping,
      c: this.c,
      n2: this.n2,
      Isat: this.Isat,
      chi: this.chi,
      boundaryAlpha: this.boundaryAlpha,
      boundaryM: this.boundaryM,
      margin: this.margin,
      lensRadius: this.lensRadius,
      fresnelZones: this.fresnelZones,
      numSectors: this.numSectors,
      learningRate: this.learningRate,
      optimizationInterval: this.optimizationInterval,
    };

    this.presets.set(name, preset);
    this.savePresets();
  }

  // Load a saved preset
  loadPreset(name) {
    const preset = this.presets.get(name);
    if (preset) {
      Object.assign(this, preset);
      this.updateDerivedValues();
      return true;
    }
    return false;
  }

  // Save presets to localStorage
  savePresets() {
    try {
      const presetsObj = Object.fromEntries(this.presets);
      localStorage.setItem("waveSimPresets", JSON.stringify(presetsObj));
    } catch (e) {
      console.warn("Failed to save presets:", e);
    }
  }

  // Load presets from localStorage
  loadPresets() {
    try {
      const saved = localStorage.getItem("waveSimPresets");
      if (saved) {
        const presetsObj = JSON.parse(saved);
        this.presets = new Map(Object.entries(presetsObj));
      }
    } catch (e) {
      console.warn("Failed to load presets:", e);
    }
  }

  // Update a single parameter
  updateParameter(name, value) {
    if (this.hasOwnProperty(name)) {
      this[name] = value;
      this.updateDerivedValues();
      return true;
    }
    return false;
  }

  // Get all parameter names
  getParameterNames() {
    return Object.keys(this).filter(
      (key) =>
        typeof this[key] !== "function" &&
        key !== "presets" &&
        key !== "boundaryR0", // Derived value
    );
  }

  // Get all preset names
  getPresetNames() {
    return Array.from(this.presets.keys());
  }

  // Export current configuration
  export() {
    const config = {};
    this.getParameterNames().forEach((name) => {
      config[name] = this[name];
    });
    return config;
  }

  // Import configuration
  import(config) {
    this.getParameterNames().forEach((name) => {
      if (config.hasOwnProperty(name)) {
        this[name] = config[name];
      }
    });
    this.updateDerivedValues();
  }

  // Create a copy of the current configuration
  clone() {
    const config = new Config();
    config.import(this.export());
    return config;
  }

  // Reset to default values
  reset() {
    Object.assign(this, new Config());
  }
}

// Example interesting presets
export const defaultPresets = {
  "High Energy": {
    dt: 0.4,
    n2: 1.5,
    Isat: 0.05,
    chi: 0.2,
    learningRate: 0.00002,
  },
  "Stable Solitons": {
    dt: 0.3,
    n2: 1.2,
    Isat: 0.15,
    chi: 0.08,
    learningRate: 0.000008,
  },
  "Chaos Mode": {
    dt: 0.45,
    n2: 1.8,
    Isat: 0.03,
    chi: 0.25,
    learningRate: 0.00003,
  },
};
