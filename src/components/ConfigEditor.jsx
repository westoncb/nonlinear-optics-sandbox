import { useState, useEffect, useCallback } from "react";
import { createPortal } from "react-dom";
import { canonicalConfigs } from "../config-manager";
import "./ConfigEditor.css";

// Detect OS for shortcut display
const isMac = navigator.platform.toLowerCase().includes("mac");
const modifierKey = isMac ? "⌘" : "Ctrl";

const initialEditorState = {
  selectedConfigId: canonicalConfigs[0].metadata.name,
  form: {
    name: canonicalConfigs[0].metadata.name,
    description: canonicalConfigs[0].metadata.description,
    data: JSON.stringify(canonicalConfigs[0].data, null, 2),
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

  const parseConfigData = (str) => {
    try {
      const sanitized = str
        .replace(/\/\/.*/g, "")
        .replace(/^\s*{/, "({")
        .replace(/}\s*$/, "})");
      const data = eval(sanitized);
      return { success: true, data };
    } catch (e) {
      return {
        success: false,
        error: `Unable to parse configuration: ${e.message}`,
      };
    }
  };

  const validateConfig = useCallback((data) => {
    const errors = {};
    const referenceConfig = canonicalConfigs[0].data;

    Object.keys(referenceConfig).forEach((key) => {
      if (!(key in data)) {
        errors[key] = `Missing required parameter: ${key}`;
      }
    });

    Object.keys(data).forEach((key) => {
      if (!(key in referenceConfig)) {
        errors[key] = `Unrecognized parameter: ${key}`;
      }
    });

    return errors;
  }, []);

  const handleConfigSelect = useCallback((config) => {
    const formattedData = JSON.stringify(config.data, null, 2);
    setState({
      selectedConfigId: config.metadata.name,
      form: {
        name: config.metadata.name,
        description: config.metadata.description,
        data: formattedData,
      },
      validation: {
        message: { text: `Press ${modifierKey}+Enter to run`, isError: false },
        parsed: config.data,
      },
      mode: config.isCanonical ? "viewing" : "editing",
      originalText: formattedData,
    });
  }, []);

  useEffect(() => {
    if (isOpen && !state.selectedConfigId) {
      handleConfigSelect(configs.get(canonicalConfigs[0].metadata.name));
    }

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
          const validationErrors = validateConfig(parseResult.data);
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

  const handleSave = () => {
    if (!state.validation.parsed) return;
    onSave(state.validation.parsed, {
      name: state.form.name,
      description: state.form.description,
    });
    setState((prev) => ({ ...prev, originalText: prev.form.data }));
  };

  const handleRun = () => {
    if (!state.validation.parsed) return;
    setConfig({
      data: state.validation.parsed,
      metadata: {
        name: state.form.name,
        description: state.form.description,
      },
    });
    onClose();
  };

  const handleDelete = () => {
    if (window.confirm("Are you sure you want to delete this configuration?")) {
      configManager.deleteConfig(state.selectedConfigId);
      onClose();
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
              {Array.from(configs.entries()).map(([name, config]) => (
                <option key={name} value={name}>
                  {name} {config.isCanonical ? "(Canonical)" : ""}
                </option>
              ))}
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
                value={state.form.data}
                onChange={(e) => updateForm({ data: e.target.value })}
                onInput={(e) => {
                  if (e.target.value === state.originalText) {
                    setState((prev) => ({ ...prev, isDirty: false }));
                  }
                }}
                disabled={state.mode === "viewing"}
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
