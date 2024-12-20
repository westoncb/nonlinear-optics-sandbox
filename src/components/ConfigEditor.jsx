import { useState, useEffect, useCallback, useRef } from "react";
import { createPortal } from "react-dom";
import { canonicalConfigs, updateDerivedValues } from "../config-manager";
import "./ConfigEditor.css";

// Detect OS for shortcut display
const isMac = navigator.platform.toLowerCase().includes("mac");
const modifierKey = isMac ? "⌘" : "Ctrl";

const formatConfigData = (data) => {
  // Convert to string with 2-space indentation
  const jsonString = JSON.stringify(data, null, 2);

  // Remove quotes around property names
  return jsonString.replace(/"([^"]+)":/g, "$1:");
};

const initialEditorState = {
  selectedConfigId: canonicalConfigs[0].metadata.name,
  form: {
    name: canonicalConfigs[0].metadata.name,
    description: canonicalConfigs[0].metadata.description,
    data: formatConfigData(canonicalConfigs[0].data),
  },
  validation: {
    message: { text: `Press ${modifierKey}+Enter to run`, isError: false },
    parsed: canonicalConfigs[0].data,
  },
  mode: "viewing",
  originalText: JSON.stringify(canonicalConfigs[0].data, null, 2),
};

export const ConfigEditor = ({
  isOpen,
  onClose,
  configManager,
  onSave,
  setConfig,
}) => {
  const [state, setState] = useState(initialEditorState);
  const configs = configManager.getAllConfigs();
  const textareaRef = useRef(null);

  const parseConfigData = (str) => {
    try {
      // Remove comments and clean up whitespace
      const sanitized = str
        .replace(/\/\/.*/g, "") // Remove single-line comments
        .replace(/\/\*[\s\S]*?\*\//g, "") // Remove multi-line comments
        .trim();

      // If it starts with a plain '{', wrap it in parentheses
      const wrapped = sanitized.startsWith("{") ? `(${sanitized})` : sanitized;

      // Use Function constructor to evaluate as JS object
      const data = new Function(`return ${wrapped}`)();
      return { success: true, data };
    } catch (e) {
      return {
        success: false,
        error: `Unable to parse configuration: ${e.message}`,
      };
    }
  };

  const validateConfigData = (configData) => {
    const validationErrors = validateConfig(configData);
    return Object.keys(validationErrors).length === 0;
  };

  const validateConfig = useCallback((data) => {
    const errors = {};
    const referenceConfig = canonicalConfigs[0];

    Object.keys(referenceConfig.data).forEach((key) => {
      if (!(key in data)) {
        errors[key] = `Missing required parameter: ${key}`;
      }
    });

    Object.keys(data).forEach((key) => {
      if (!(key in referenceConfig.data)) {
        errors[key] = `Unrecognized parameter: ${key}`;
      }
    });

    return errors;
  }, []);

  useEffect(() => {
    const handleKeyDown = (event) => {
      if (event.key === "Escape") {
        onClose();
      } else if (event.key === "Enter" && (event.metaKey || event.ctrlKey)) {
        event.preventDefault();
        handleRun();
      }
    };

    if (isOpen) {
      document.addEventListener("keydown", handleKeyDown);
      return () => document.removeEventListener("keydown", handleKeyDown);
    }
  }, [isOpen, state.selectedConfigId]);

  const updateForm = useCallback(
    (updates) => {
      setState((prev) => {
        const newForm = { ...prev.form, ...updates };
        const parseResult = parseConfigData(newForm.data);

        let message;
        if (!parseResult.success) {
          message = { text: parseResult.error, isError: true };
        } else {
          const validationErrors = validateConfig(
            updateDerivedValues(parseResult).data,
          );
          const errorMessages = Object.values(validationErrors);
          message =
            errorMessages.length > 0
              ? { text: errorMessages.join("; "), isError: true }
              : { text: `Press ${modifierKey}+Enter to run`, isError: false };
        }

        return {
          ...prev,
          form: newForm,
          validation: {
            message,
            parsed: parseResult.success ? parseResult.data : null,
          },
        };
      });
    },
    [validateConfig],
  );

  const handleConfigSelect = useCallback(
    (config) => {
      // First update the basic state that updateForm doesn't handle
      setState((prev) => ({
        ...prev,
        selectedConfigId: config.metadata.name,
        mode: config.isCanonical ? "viewing" : "editing",
        originalText: formatConfigData(config.data),
      }));

      // Then use updateForm to handle the form state and validation
      updateForm({
        name: config.metadata.name,
        description: config.metadata.description,
        data: formatConfigData(config.data),
      });
    },
    [updateForm],
  );

  const handleSave = () => {
    if (!state.validation.parsed) return;
    onSave(state.validation.parsed, {
      name: state.form.name,
      description: state.form.description,
    });
    setState((prev) => ({ ...prev, originalText: prev.form.data }));
  };

  const handleRun = () => {
    // Get the current value directly from the textarea
    const currentValue = textareaRef.current?.value || state.form.data;

    const parseResult = parseConfigData(currentValue);
    if (!parseResult.success) return;

    const configData = updateDerivedValues(parseResult).data;
    const validationErrors = validateConfig(configData);
    if (Object.keys(validationErrors).length > 0) {
      console.log("Validation Errors:", validationErrors, configData);
      return;
    }

    setConfig({
      data: configData,
      metadata: {
        name: state.form.name,
        description: state.form.description,
      },
    });
    onClose();
  };

  const handleDelete = () => {
    if (window.confirm("Are you sure you want to delete this configuration?")) {
      // Get list of configs before deletion
      const configList = Array.from(configs.entries());
      const currentIndex = configList.findIndex(
        ([name]) => name === state.selectedConfigId,
      );

      // Delete the config
      configManager.deleteConfig(state.selectedConfigId);

      // Get updated list post-deletion
      const updatedConfigs = configManager.getAllConfigs();
      const updatedConfigList = Array.from(updatedConfigs.entries());

      // If there are any configs left, select the next one (with wraparound)
      if (updatedConfigList.length > 0) {
        const nextIndex = currentIndex % updatedConfigList.length;
        const [nextConfigName, nextConfig] = updatedConfigList[nextIndex];
        handleConfigSelect(nextConfig);
      }
    }
  };

  if (!isOpen) return null;

  const hasErrors = state.validation.message.isError;
  const isCanonical = configs.get(state.selectedConfigId)?.isCanonical;
  const isDirty = state.form.data !== state.originalText;

  return createPortal(
    <div className="modal-overlay" onClick={onClose}>
      <div
        className="modal-content ${isDirty ? 'modified' : ''}"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="modal-header">
          <h2>Configuration Editor</h2>
        </div>

        <div className="modal-body">
          <div className="control-group">
            <select
              value={state.selectedConfigId}
              onChange={(e) => {
                const config = configs.get(e.target.value);
                if (config) handleConfigSelect(config);
              }}
            >
              {Array.from(configs.entries()).map(([name, config]) => {
                const isValid = validateConfigData(config.data);
                return (
                  <option
                    key={name}
                    value={name}
                    className={`config-option ${!isValid ? "config-option-error" : ""}`}
                  >
                    {name} {config.isCanonical ? "(Canonical)" : ""}
                    {!isValid ? " ⚠" : ""}{" "}
                  </option>
                );
              })}
            </select>

            <button
              onClick={() =>
                handleConfigSelect(configManager.createBlankConfig(""))
              }
            >
              New
            </button>
            <button
              onClick={() => {
                const config = configs.get(state.selectedConfigId);
                if (config) {
                  handleConfigSelect(
                    configManager.forkConfig(
                      config,
                      `${config.metadata.name} (Copy)`,
                    ),
                  );
                }
              }}
            >
              Fork
            </button>
          </div>

          <div className="form-fields">
            <div className="field inline">
              <label>Name:</label>
              <input
                value={state.form.name}
                onChange={(e) => updateForm({ name: e.target.value })}
                disabled={state.mode === "viewing"}
              />
            </div>

            <div className="field inline">
              <label>Description:</label>
              <input
                value={state.form.description}
                onChange={(e) => updateForm({ description: e.target.value })}
                disabled={state.mode === "viewing"}
              />
            </div>

            <div className="field">
              <label>Configuration Data:</label>
              <textarea
                ref={textareaRef}
                value={state.form.data}
                onChange={(e) => updateForm({ data: e.target.value })}
                onInput={(e) => {
                  if (e.target.value === state.originalText) {
                    setState((prev) => ({ ...prev, isDirty: false }));
                  }
                }}
                spellCheck={false}
              />
              <div
                className={`message ${state.validation.message.isError ? "error" : "info"}`}
              >
                {state.validation.message.text}
              </div>
            </div>
          </div>
        </div>

        <div className="modal-footer">
          <div className="left-buttons">
            <button
              onClick={handleSave}
              disabled={hasErrors || state.mode === "viewing" || !isDirty}
              className="save-button"
            >
              Save
            </button>

            <button
              onClick={handleDelete}
              disabled={isCanonical}
              className="delete-button"
            >
              Delete
            </button>
          </div>

          <div className="right-buttons">
            <button onClick={onClose} className="cancel-button">
              Cancel
            </button>
            <button
              onClick={handleRun}
              disabled={hasErrors}
              className="run-button"
            >
              Run
            </button>
          </div>
        </div>
      </div>
    </div>,
    document.body,
  );
};
