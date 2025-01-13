// @hmr-reset

import { useState, useEffect, useRef } from "react";
import { RenderSystem } from "../render-system";
import { ConfigManager, constants, canonicalConfigs } from "../config-manager";
import { ProgressGraph } from "./ProgressGraph";
import { ConfigEditor } from "./ConfigEditor";
import "./App.css";

export const App = () => {
  const mainContainerRef = useRef(null);
  const renderSystemRef = useRef(null);
  const [isConfigModalOpen, setIsConfigModalOpen] = useState(false);
  const [getProgress, setGetProgress] = useState(() => () => [
    {
      iteration: 0,
      metric: 0,
    },
  ]);

  const configManager = useRef(new ConfigManager()).current;
  const [config, setConfig] = useState(canonicalConfigs[0]);
  const [gridSize, setGridSize] = useState(config.data.gridSize);

  const [isFullscreen, setIsFullscreen] = useState(false);

  useEffect(() => {
    const handleEsc = (e) => {
      if (e.key === "Escape" && isFullscreen) {
        setIsFullscreen(false);
      }
    };
    window.addEventListener("keydown", handleEsc);
    return () => window.removeEventListener("keydown", handleEsc);
  }, [isFullscreen]);

  // Track which content is shown in each canvas
  const [canvasContent, setCanvasContent] = useState({
    primary: "Fundamental field",
    preview1: "Lens Profile",
    preview2: "2nd harmonics",
  });

  const handleConfigSave = (configData, metadata) => {
    const newConfig = configManager.saveUserConfig(configData, metadata);
    setConfig(newConfig);
  };

  // Initialize App
  useEffect(() => {
    const initApp = async () => {
      renderSystemRef.current = new RenderSystem(config);
      await renderSystemRef.current.initialize();
      setGridSize(renderSystemRef.current.config.gridSize);

      // Set initial sizes
      const size = calculatePrimarySize(mainContainerRef);
      if (size) {
        const primaryCanvas = document.getElementById("primaryCanvas");
        primaryCanvas.style.width = `${size}px`;
        primaryCanvas.style.height = `${size}px`;
      }

      // Set up progress function once app is initialized
      if (renderSystemRef.current?.lensOptimizer.getProgress) {
        setGetProgress(
          () => () => renderSystemRef.current?.lensOptimizer.getProgress(),
        );
      }

      renderSystemRef.current.start();
    };

    initApp();

    // Cleanup function
    return () => {
      if (renderSystemRef.current) {
        renderSystemRef.current.stop();
        renderSystemRef.current.cleanup();
        renderSystemRef.current = null;
      }
    };
  }, [config]);

  const calculatePrimarySize = (containerRef) => {
    if (!containerRef.current) return null;
    const containerHeight = containerRef.current.clientHeight;
    return Math.min(containerHeight, window.innerWidth);
  };

  const handlePreviewClick = (clickedPreviewNum) => {
    if (!renderSystemRef.current) return;

    const currentPrimary = renderSystemRef.current.displayModes.primary;
    const currentPreview1 = renderSystemRef.current.displayModes.preview1;
    const currentPreview2 = renderSystemRef.current.displayModes.preview2;

    // Update display modes and content labels
    if (clickedPreviewNum === 1) {
      renderSystemRef.current.setDisplayMode({
        primary: currentPreview1,
        preview1: currentPrimary,
      });

      setCanvasContent((prev) => ({
        ...prev,
        primary: prev.preview1,
        preview1: prev.primary,
      }));
    } else if (clickedPreviewNum === 2) {
      renderSystemRef.current.setDisplayMode({
        primary: currentPreview2,
        preview2: currentPrimary,
      });

      setCanvasContent((prev) => ({
        ...prev,
        primary: prev.preview2,
        preview2: prev.primary,
      }));
    }
  };

  // Handle resize
  useEffect(() => {
    if (!mainContainerRef.current) return;

    const observer = new ResizeObserver(() => {
      const size = calculatePrimarySize(mainContainerRef);
      if (size) {
        const primaryCanvas = document.getElementById("primaryCanvas");
        primaryCanvas.style.width = `${size}px`;
        primaryCanvas.style.height = `${size}px`;
      }
    });

    observer.observe(mainContainerRef.current);
    return () => observer.disconnect();
  }, []);

  return (
    <div className="app">
      <div className="top-bar">
        <div
          style={{
            display: "flex",
            flexDirection: "column",
            justifyContent: "center",
            flex: "1 1 auto", // grow, shrink, auto basis
            minWidth: 0, // allows shrinking
          }}
        >
          <div className="preview-label" style={{ marginBottom: "0.5rem" }}>
            Loss -{" "}
            {config.data.disableAdaptation ? "" : config.data.updateStrategy}
          </div>
          <div style={{ flex: "1 1 auto" }}>
            <ProgressGraph
              getProgress={getProgress}
              updateRate={50}
              maxPoints={40000}
            />
          </div>
          <hr
            style={{
              margin: "1rem",
              border: "1px solid #333",
            }}
          />
          <button
            className="config-button"
            onClick={() => setIsConfigModalOpen(true)}
            aria-label="Open configuration settings"
          >
            Configure
          </button>
        </div>
        <div className="spacer" />
        {/* Preview 1 */}
        <div className="labeled-container">
          <canvas
            id="preview1Canvas"
            width={gridSize}
            height={gridSize}
            style={{
              width: `${constants.PREVIEW_SIZE}px`,
              height: `${constants.PREVIEW_SIZE}px`,
            }}
            className="preview-canvas"
            onClick={() => handlePreviewClick(1)}
          />
          <span className="preview-label">{canvasContent.preview1}</span>
        </div>
        {/* Preview 2 */}
        <div className="labeled-container">
          <canvas
            id="preview2Canvas"
            width={gridSize}
            height={gridSize}
            style={{
              width: `${constants.PREVIEW_SIZE}px`,
              height: `${constants.PREVIEW_SIZE}px`,
            }}
            className="preview-canvas"
            onClick={() => handlePreviewClick(2)}
          />
          <span className="preview-label">{canvasContent.preview2}</span>
        </div>
      </div>

      <div
        className="main-container"
        ref={mainContainerRef}
        onClick={() => !isFullscreen && setIsFullscreen(true)}
      >
        <span className="main-label">{canvasContent.primary}</span>
        <div
          className={`canvas-wrapper ${isFullscreen ? "fullscreen-overlay" : ""}`}
        >
          {isFullscreen && (
            <div className="fullscreen-header">
              <button
                onClick={() => setIsFullscreen(false)}
                aria-label="Close fullscreen view"
              >
                ×
              </button>
            </div>
          )}
          <canvas
            id="primaryCanvas"
            width={gridSize}
            height={gridSize}
            className={`primary-canvas ${isFullscreen ? "fullscreen-canvas" : ""}`}
            onClick={(e) => {
              e.stopPropagation();
              !isFullscreen && setIsFullscreen(true);
            }}
            style={{ cursor: isFullscreen ? "default" : "pointer" }}
          />
        </div>
        <canvas
          id="simulationCanvas"
          width={gridSize}
          height={gridSize}
          style={{
            position: "absolute",
            top: 0,
            left: 0,
            visibility: "hidden",
          }}
          className="simulation-canvas"
        />
      </div>

      <ConfigEditor
        isOpen={isConfigModalOpen}
        onClose={() => setIsConfigModalOpen(false)}
        configManager={configManager}
        onSave={handleConfigSave}
        setConfig={setConfig}
      />
    </div>
  );
};
