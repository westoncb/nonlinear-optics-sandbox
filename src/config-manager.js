export const constants = {
  GRID_SIZE: 512,
  PREVIEW_SIZE: 200,
};

// Default configuration values
const defaultConfigData = {
  // Core simulation parameters
  // Temporal step size for numerical integration of wave equation
  dt: 0.2,
  // Spatial grid spacing
  dx: 1,
  // Number of grid points in each dimension (square grid)
  gridSize: constants.GRID_SIZE,
  // Field amplitude preservation factor (1 = perfect preservation)
  damping: 1,
  // Speed of light in simulation units
  c: 1,
  // Toggle between 5-point and 9-point Laplacian stencil
  use9PointStencil: true,

  // Nonlinear optical parameters
  // Base strength of nonlinear optical effects
  chi: 0.8,
  // Relative strength of Kerr (χ(3)) nonlinearity
  chi_ratio: -12,
  // Relative strength of second-harmonic (χ(2)) nonlinearity
  chi2_ratio: 1.0,
  // Intensity at which SHG conversion begins to saturate
  shg_Isat: 0.02,
  // Intensity at which Kerr effect begins to saturate
  kerr_Isat: 0.03,

  // Coupling parameters
  // Strength of Kerr-mediated interaction between fundamental and SHG fields
  crossKerrCoupling: 0.001,
  // Efficiency of frequency conversion between fundamental and SHG
  conversionCoupling: 0.005,

  // Gain and loss parameters
  // Maximum small-signal gain coefficient
  gain0: 1,
  // Intensity at which gain medium begins to saturate
  gainSat: 0.25,
  // Rate of passive power loss in cavity
  linearLoss: 0.0,

  // Wave properties
  // Wavelength of fundamental field in grid units
  lambdaFund: 20.0,
  // Wavelength of second harmonic (must be lambdaFund/2 for phase matching)
  lambdaSHG: 10.0,
  // Reference phase for evaluating phase matching conditions
  phaseRef: Math.PI,

  // Cavity geometry
  // Amplitude of boundary shape modulation (0 = perfect circle)
  boundaryAlpha: 0.2,
  // Number of symmetric perturbations around cavity boundary
  boundaryM: 8,
  // Width of border region in grid points
  margin: 10,
  // Width of smooth transition region at boundary (in grid points)
  boundaryTransitionWidth: 20,

  // Adaptive lens parameters
  // Radius of lens region in grid points
  lensRadius: 64,
  // Number of phase wraps from center to edge of lens
  fresnelZones: 64,
  // Number of independently adjustable regions in lens
  numSectors: 200,

  // Optimization controls
  // Selected optimization algorithm
  updateStrategy: "SHGTest",
  // Step size for lens parameter updates
  learningRate: 1e-8,
  // Update lens every N simulation steps
  optimizationInterval: 1,
  // Toggle lens optimization
  disableAdaptation: false,

  // Number of steps between pulse injections (0 = single pulse)
  pulseInterval: 0,
  // Peak amplitude of all pulses after initial
  subsequentPulseAmplitude: 0.1,
  // Beam waist/width for a Gaussian pulse in grid units
  beamWidth: 4,

  // Initial pulse peak amplitude
  initialPulseAmplitude: 20,
  // Initial pulse phase offset
  initialPulsePhaseShift: 0,
  // Amplitude of random fluctuations in initial field
  initialNoiseScale: 1e-8,
};

export const referenceConfig = {
  data: defaultConfigData,
  metadata: {
    name: "Reference config",
    description: "default canonical config without derived values",
  },
};

export const updateDerivedValues = (fullConfig) => {
  // Create a new config object with gridSize enforced from constants
  const newData = {
    ...fullConfig.data,
    gridSize: constants.GRID_SIZE,
  };

  return {
    ...fullConfig,
    data: {
      ...newData,
      boundaryR0:
        (constants.GRID_SIZE / 2 - newData.margin) /
        (1 + newData.boundaryAlpha),
    },
  };
};

export const canonicalConfigs = [
  {
    data: defaultConfigData,
    metadata: {
      name: "Test Config",
      description: "typically what i use to test new ideas quickly",
    },
  },
].map(updateDerivedValues);

export class ConfigManager {
  constructor() {
    this.userConfigs = new Map();
    this.loadUserConfigs();
  }

  // Get all available configs (canonical + user)
  getAllConfigs() {
    const configs = new Map();

    // Add canonical configs with derived values
    canonicalConfigs.forEach((config) => {
      configs.set(config.metadata.name, {
        ...updateDerivedValues(config),
        isCanonical: true,
      });
    });

    // Add user configs (will override canonical if same name)
    this.userConfigs.forEach((config, name) => {
      configs.set(name, {
        ...updateDerivedValues(config),
        isCanonical: false,
      });
    });

    return configs;
  }

  // Get a single config with derived values
  getConfig(name) {
    // Check canonical configs first
    const canonicalConfig = canonicalConfigs.find(
      (c) => c.metadata.name === name,
    );
    if (canonicalConfig) {
      return {
        ...updateDerivedValues(canonicalConfig),
        isCanonical: true,
      };
    }

    // Then check user configs
    const userConfig = this.userConfigs.get(name);
    if (userConfig) {
      return {
        ...updateDerivedValues(userConfig),
        isCanonical: false,
      };
    }

    return null;
  }

  // Save a new user config - store raw values
  saveUserConfig(data, metadata) {
    const config = { data, metadata };
    this.userConfigs.set(metadata.name, config);
    this.persistUserConfigs();
    return this.getConfig(metadata.name); // Return with derived values
  }

  // Delete a user config
  deleteConfig(name) {
    // Check if it's a canonical config
    if (canonicalConfigs.some((c) => c.metadata.name === name)) {
      throw new Error(`Cannot delete canonical config "${name}"`);
    }

    // Check if the config exists
    if (!this.userConfigs.has(name)) {
      throw new Error(`User config "${name}" not found`);
    }

    this.userConfigs.delete(name);
    this.persistUserConfigs();
    return true;
  }

  // Load user configs from localStorage - store raw values
  loadUserConfigs() {
    try {
      const saved = localStorage.getItem("userConfigs");
      if (saved) {
        const configs = JSON.parse(saved);
        this.userConfigs = new Map(Object.entries(configs));
      }
    } catch (e) {
      console.warn("Failed to load user configs:", e);
      this.userConfigs = new Map();
    }
  }

  // Save user configs to localStorage - store raw values
  persistUserConfigs() {
    try {
      const configsObj = Object.fromEntries(this.userConfigs);
      localStorage.setItem("userConfigs", JSON.stringify(configsObj));
    } catch (e) {
      console.warn("Failed to save user configs:", e);
    }
  }

  // Helper to create a new config from existing one
  forkConfig(existingConfig, newName) {
    // Remove any derived values before forking
    const baseConfig = {
      data: { ...existingConfig.data },
      metadata: {
        name: newName,
        description: `Forked from ${existingConfig.metadata.name}`,
      },
    };

    return this.saveUserConfig(baseConfig.data, baseConfig.metadata);
  }

  // Create a new blank config
  createBlankConfig(name) {
    return this.saveUserConfig(
      { ...defaultConfigData },
      {
        name,
        description: "New configuration",
      },
    );
  }
}
