export const constants = {
  GRID_SIZE: 512,
  PREVIEW_SIZE: 200,
};

// Default configuration values
const defaultConfigData = {
  dt: 0.2,
  dx: 1,
  gridSize: constants.GRID_SIZE,
  damping: 1,
  c: 1,
  use9PointStencil: true,
  chi: 0.9, // Even stronger nonlinearity
  chi_ratio: -8, // Increased χ(3) effect
  shg_Isat: 0.06, // Lower saturation threshold
  kerr_Isat: 0.08, // Lower Kerr saturation
  boundaryAlpha: 0.1,
  boundaryM: 8, // Reduced for cleaner geometric patterns
  margin: 0,
  boundaryReflectivity: "1.",
  lensRadius: 64, // Larger interaction region
  fresnelZones: 12, // More radial resolution
  numSectors: 84, // Even finer angular resolution
  hillPower: 2,
  updateStrategy: "solitonStrategy",
  learningRate: 1e-5,
  optimizationInterval: 1,
  pulseInterval: 200,
  initialPulseAmplitude: 24,
  initialPulsePhaseShift: 0.0,
  initialNoiseScale: 0,
  initialBeamWidth: 3,
  useSimplePulse: false,
  disableAdaptation: false,
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
