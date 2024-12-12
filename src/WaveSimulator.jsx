import { useState, useEffect, useRef } from "react";
import { App } from "./app";
import { constants } from "./config";
import { ProgressGraph } from "./ProgressGraph";
import "./WaveSimulator.css";

export const WaveSimulator = () => {
  const mainContainerRef = useRef(null);
  const appRef = useRef(null);
  const [gridSize, setGridSize] = useState(512); // Default until App initializes
  const [progressData, setProgressData] = useState([]);

  // Initialize App
  useEffect(() => {
    const initApp = async () => {
      appRef.current = new App();
      await appRef.current.initialize();
      setGridSize(appRef.current.config.gridSize);

      // Set initial sizes
      const size = calculatePrimarySize(mainContainerRef);
      if (size) {
        const primaryCanvas = document.getElementById("primaryCanvas");
        primaryCanvas.style.width = `${size}px`;
        primaryCanvas.style.height = `${size}px`;
      }

      appRef.current.start();
    };

    initApp();
  }, []);

  useEffect(() => {
    const interval = setInterval(() => {
      // Poll simulation progress
      const p = appRef.current?.simulation.getProgress();
      setProgressData([...p]);
    }, 50);
    return () => clearInterval(interval);
  }, [appRef.current?.simulation]);

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

    if (clickedPreviewNum === 1) {
      appRef.current.setDisplayMode({
        primary: currentPreview1,
        preview1: currentPrimary,
      });
    } else if (clickedPreviewNum === 2) {
      appRef.current.setDisplayMode({
        primary: currentPreview2,
        preview2: currentPrimary,
      });
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
        <div style={{ flex: "1 1 auto" }}>
          <ProgressGraph data={progressData} />
        </div>
        <div className="spacer" />
        {/* Preview 1 */}
        <div className="preview-container">
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
          <span className="preview-label">Lens Profile</span>
        </div>

        {/* Preview 2 */}
        <div className="preview-container">
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
          <span className="preview-label">SHG Output</span>
        </div>
      </div>

      {/* <div className="status-bar">
        <div className="status-content">
          Status: Running | Frame: 1000 | Time: 1.2ms
        </div>
      </div> */}

      <div className="main-container" ref={mainContainerRef}>
        <span className="main-label">Wave Propagation</span>
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
    </div>
  );
};
