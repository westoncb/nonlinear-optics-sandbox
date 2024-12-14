import { useState, useEffect, useRef } from "react";
import { App } from "./app";
import { ConfigManager, constants, canonicalConfigs } from "./config-manager";
import { ProgressGraph } from "./ProgressGraph";
import { ConfigModal } from "./ConfigModal";
import "./WaveSimulator.css";

export const WaveSimulator = () => {
  const mainContainerRef = useRef(null);
  const appRef = useRef(null);
  const [isConfigModalOpen, setIsConfigModalOpen] = useState(false);
  const [gridSize, setGridSize] = useState(constants.GRID_SIZE);
  const [getProgress, setGetProgress] = useState(() => () => [
    {
      iteration: 0,
      metric: 0,
    },
  ]);
  const [config, setConfig] = useState(canonicalConfigs[0]);
  const configManager = useRef(new ConfigManager()).current;

  // Track which content is shown in each canvas
  const [canvasContent, setCanvasContent] = useState({
    primary: "Fundamental field",
    preview1: "Lens Profile",
    preview2: "2nd harmonics",
  });

  const handleConfigSave = (configData, metadata) => {
    const newConfig = configManager.saveConfig(configData, metadata);
    setConfig(newConfig);
  };

  // Initialize App
  useEffect(() => {
    const initApp = async () => {
      appRef.current = new App(config);
      await appRef.current.initialize();
      setGridSize(appRef.current.config.gridSize);

      // Set initial sizes
      const size = calculatePrimarySize(mainContainerRef);
      if (size) {
        const primaryCanvas = document.getElementById("primaryCanvas");
        primaryCanvas.style.width = `${size}px`;
        primaryCanvas.style.height = `${size}px`;
      }

      // Set up progress function once app is initialized
      if (appRef.current?.simulation.getProgress) {
        setGetProgress(() => () => appRef.current.simulation.getProgress());
      }

      appRef.current.start();
    };

    initApp();

    // Cleanup function
    return () => {
      if (appRef.current) {
        appRef.current.stop();
        appRef.current.cleanup();
        appRef.current = null;
      }
    };
  }, [config]);

  const calculatePrimarySize = (containerRef) => {
    if (!containerRef.current) return null;
    const containerHeight = containerRef.current.clientHeight;
    return Math.min(containerHeight - 40, window.innerWidth);
  };

  const handlePreviewClick = (clickedPreviewNum) => {
    if (!appRef.current) return;

    const currentPrimary = appRef.current.displayModes.primary;
    const currentPreview1 = appRef.current.displayModes.preview1;
    const currentPreview2 = appRef.current.displayModes.preview2;

    // Update display modes and content labels
    if (clickedPreviewNum === 1) {
      appRef.current.setDisplayMode({
        primary: currentPreview1,
        preview1: currentPrimary,
      });

      setCanvasContent((prev) => ({
        ...prev,
        primary: prev.preview1,
        preview1: prev.primary,
      }));
    } else if (clickedPreviewNum === 2) {
      appRef.current.setDisplayMode({
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
    <div className="simulator">
      <div className="top-bar">
        <div
          style={{
            display: "flex",
            flexDirection: "column",
            justifyContent: "center",
            flexGrow: 1,
          }}
        >
          <div className="preview-label" style={{ marginBottom: "0.5rem" }}>
            convergence
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

      <div className="main-container" ref={mainContainerRef}>
        <span className="main-label">{canvasContent.primary}</span>
        <div className="canvas-wrapper">
          <canvas
            id="primaryCanvas"
            width={gridSize}
            height={gridSize}
            className="primary-canvas"
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

      <ConfigModal
        isOpen={isConfigModalOpen}
        onClose={() => setIsConfigModalOpen(false)}
        configManager={configManager}
        onSave={handleConfigSave}
        setConfig={setConfig}
      />
    </div>
  );
};
