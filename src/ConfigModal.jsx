import { useState, useEffect } from "react";
import { createPortal } from "react-dom";
import { canonicalConfigs } from "./config-manager";
import "./ConfigModal.css";

export const ConfigModal = ({
  isOpen,
  onClose,
  configManager,
  onSave,
  setConfig,
}) => {
  const [selectedConfig, setSelectedConfig] = useState(null);
  const [name, setName] = useState("");
  const [description, setDescription] = useState("");
  const [configData, setConfigData] = useState("");
  const [initialConfigData, setInitialConfigData] = useState("");
  const [isEditing, setIsEditing] = useState(false);
  const [error, setError] = useState("");

  const configs = configManager.getAllConfigs();

  useEffect(() => {
    if (isOpen && !selectedConfig) {
      const defaultConfig = {
        data: canonicalConfigs[0].data,
        metadata: { name: "", description: "" },
        isCanonical: false,
      };
      handleConfigSelect(defaultConfig);
    }

    const handleKeyDown = (event) => {
      if (event.key === "Escape") {
        handleCancel();
      } else if (event.key === "Enter" && !event.shiftKey) {
        event.preventDefault();
        handleRunOnly();
      }
    };

    if (isOpen) {
      document.body.style.overflow = "hidden";
      document.addEventListener("keydown", handleKeyDown);
    } else {
      document.body.style.overflow = "unset";
    }

    return () => {
      document.body.style.overflow = "unset";
      document.removeEventListener("keydown", handleKeyDown);
    };
  }, [isOpen]);

  const validateAndParseConfig = () => {
    try {
      const parsedData = JSON.parse(configData);
      return {
        success: true,
        data: parsedData,
        metadata: { name, description },
      };
    } catch (e) {
      setError("Invalid JSON configuration data");
      return { success: false };
    }
  };

  const handleConfigSelect = (config) => {
    setSelectedConfig(config);
    setName(config.metadata.name);
    setDescription(config.metadata.description);
    const formattedData = JSON.stringify(config.data, null, 2);
    setConfigData(formattedData);
    setInitialConfigData(formattedData);
    setIsEditing(!config.isCanonical);
    setError("");
  };

  const handleNew = () => {
    const newConfig = configManager.createBlankConfig("");
    handleConfigSelect({ ...newConfig, isCanonical: false });
  };

  const handleFork = () => {
    if (!selectedConfig) return;
    const forkedConfig = configManager.forkConfig(
      selectedConfig,
      `${selectedConfig.metadata.name} (Copy)`,
    );
    handleConfigSelect({ ...forkedConfig, isCanonical: false });
  };

  const handleCancel = () => {
    setConfigData(initialConfigData);
    onClose();
  };

  const handleRunOnly = () => {
    const result = validateAndParseConfig();
    if (result.success) {
      setConfig(result);
      onClose();
    }
  };

  const handleSave = () => {
    const result = validateAndParseConfig();
    if (result.success) {
      onSave(result.data, result.metadata);
      onClose();
    }
  };

  if (!isOpen) return null;

  return createPortal(
    <div className="modal-overlay" onClick={handleCancel}>
      <div className="modal-content" onClick={(e) => e.stopPropagation()}>
        <div className="modal-header">
          <h2 className="modal-title">Configuration Editor</h2>
        </div>

        <div className="modal-controls">
          <select
            value={selectedConfig?.metadata.name || ""}
            onChange={(e) => {
              const config = Array.from(configs.values()).find(
                (c) => c.metadata.name === e.target.value,
              );
              if (config) handleConfigSelect(config);
            }}
          >
            <option value="">Select a configuration</option>
            {Array.from(configs.entries()).map(([name, config]) => (
              <option key={name} value={name}>
                {name} {config.isCanonical ? "(Canonical)" : ""}
              </option>
            ))}
          </select>

          <button onClick={handleNew}>New</button>
          <button onClick={handleFork} disabled={!selectedConfig}>
            Fork
          </button>
        </div>

        <div className="modal-fields">
          <div className="field">
            <label htmlFor="name">Name:</label>
            <input
              id="name"
              value={name}
              onChange={(e) => setName(e.target.value)}
              disabled={!isEditing}
            />
          </div>

          <div className="field">
            <label htmlFor="description">Description:</label>
            <input
              id="description"
              value={description}
              onChange={(e) => setDescription(e.target.value)}
              disabled={!isEditing}
            />
          </div>
        </div>

        <div className="config-editor">
          <label htmlFor="config">Configuration Data</label>
          <textarea
            id="config"
            value={configData}
            onChange={(e) => setConfigData(e.target.value)}
            disabled={!isEditing}
            spellCheck="false"
            onKeyDown={(e) => {
              if (e.key === "Enter" && !e.shiftKey) {
                e.preventDefault();
                handleRunOnly();
              }
            }}
          />
          {error && <p className="error-message">{error}</p>}
        </div>

        <div className="modal-actions">
          <button onClick={handleCancel}>Cancel</button>
          <button
            onClick={handleRunOnly}
            disabled={!isEditing}
            className="run-only"
          >
            Run Only
          </button>
          <button
            onClick={handleSave}
            disabled={!isEditing}
            className="primary"
          >
            Save
          </button>
        </div>
      </div>
    </div>,
    document.body,
  );
};
