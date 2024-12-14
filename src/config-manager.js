export const constants = {
  PREVIEW_SIZE: 200,
  GRID_SIZE: 256, // because of the some react/canvas nuances this needs to be available before first render; putting it here as a workaround
};

// Default configuration values
const defaultConfigData = {
  // Core simulation parameters
  dt: 0.1,
  dx: 1.0,
  gridSize: constants.GRID_SIZE,
  damping: 0.9999,
  c: 1.0,
  use9PointStencil: true,

  // Nonlinear parameters
  n2: 2.5,
  Isat: 0.25,
  chi: 0.25,

  // Boundary parameters
  boundaryAlpha: 0.3,
  boundaryM: 8,
  margin: 10,
  boundaryReflectivity: "1.0",

  // Lens parameters
  lensRadius: 128,
  fresnelZones: 64,
  numSectors: 1000,
  hillPower: 2,
  updateStrategy: "hierarchicalZones",

  // Optimization parameters
  learningRate: 0.0000002,
  optimizationInterval: 1,

  // Initial pulse parameters
  initialPulseAmplitude: 90,
  initialPulsePhaseShift: 0.07,

  disableAdaptation: false,
};

const updateDerivedValues = (fullConfig) => {
  const config = fullConfig.data;

  return {
    ...fullConfig,
    data: {
      ...config,
      boundaryR0:
        (config.gridSize / 2 - config.margin) / (1 + config.boundaryAlpha),
    },
  };
};

// Canonical configs that ship with the code
export const canonicalConfigs = [
  {
    data: defaultConfigData,
    metadata: {
      name: "High Learning Rate",
      description: "Configuration with increased learning rate and hill power",
    },
  },
  // Add more canonical configs here
].map(updateDerivedValues);

export class ConfigManager {
  constructor() {
    this.userConfigs = new Map();
    this.loadUserConfigs();
  }

  // Get all available configs (canonical + user)
  getAllConfigs() {
    const configs = new Map();

    // Add canonical configs first
    canonicalConfigs.map(updateDerivedValues).forEach((config) => {
      configs.set(config.metadata.name, {
        ...config,
        isCanonical: true,
      });
    });

    // Add user configs (will override canonical if same name)
    this.userConfigs.forEach((config, name) => {
      configs.set(name, {
        ...config,
        isCanonical: false,
      });
    });

    return configs;
  }

  // Save a new user config
  saveUserConfig(data, metadata) {
    const config = { data, metadata };
    this.userConfigs.set(metadata.name, config);
    this.persistUserConfigs();
    return config;
  }

  // Load user configs from localStorage
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

  // Save user configs to localStorage
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
    return {
      data: { ...existingConfig.data },
      metadata: {
        name: newName,
        description: `Forked from ${existingConfig.metadata.name}`,
      },
    };
  }

  // Create a new blank config
  createBlankConfig(name) {
    return {
      data: { ...defaultConfigData },
      metadata: {
        name,
        description: "New configuration",
      },
    };
  }
}
