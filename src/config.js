export class Config {
  constructor() {
    // Core simulation parameters
    this.dt = 0.4;
    this.dx = 1;
    this.gridSize = 512;
    this.damping = 1.0;
    this.c = 1.0;

    // Nonlinear parameters
    this.n2 = 1.1;
    this.Isat = 0.1;
    this.chi = 0.1;

    // Boundary parameters
    this.boundaryAlpha = 0.1;
    this.boundaryM = 4;
    this.margin = 30;

    // Lens parameters
    this.lensRadius = 120;
    this.fresnelZones = 1000;
    this.numSectors = 2000;

    // Optimization parameters
    this.learningRate = 0.00001;
    this.optimizationInterval = 1;

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
