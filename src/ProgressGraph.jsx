import React, { useEffect, useRef, useState } from "react";
import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";

export const ProgressGraph = ({
  getProgress,
  updateRate = 100,
  maxPoints = 20000,
}) => {
  const containerRef = useRef(null);
  const plotRef = useRef(null);
  const timeoutRef = useRef(null);
  const [data, setData] = useState([]);

  // Store last values to avoid duplicate updates
  const lastValuesRef = useRef({ iteration: 0, metric: 0 });

  // Initialize plot
  useEffect(() => {
    if (!containerRef.current) return;

    const rect = containerRef.current.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) return;

    const opts = {
      width: rect.width,
      height: rect.height,
      padding: [20, 8, 8, 8],
      cursor: {
        show: false,
      },
      series: [
        {
          show: false,
        },
        {
          stroke: "#4af",
          width: 2,
          points: {
            show: false,
          },
          paths: uPlot.paths.linear(),
        },
      ],
      scales: {
        x: {
          time: false,
          auto: true,
        },
        y: {
          auto: true,
        },
      },
      axes: [{ show: false }, { show: false }],
      grid: {
        show: true,
        stroke: "#444",
        width: 1,
      },
      legend: {
        show: false,
      },
      focus: {
        show: false,
      },
    };

    // Initial empty data
    plotRef.current = new uPlot(opts, [[], []], containerRef.current);

    // Handle resize
    const resizeObserver = new ResizeObserver((entries) => {
      if (!plotRef.current) return;

      const entry = entries[0];
      if (entry.contentRect.width > 0 && entry.contentRect.height > 0) {
        plotRef.current.setSize({
          width: entry.contentRect.width,
          height: entry.contentRect.height,
        });
      }
    });

    resizeObserver.observe(containerRef.current);

    return () => {
      resizeObserver.disconnect();
      if (plotRef.current) {
        plotRef.current.destroy();
        plotRef.current = null;
      }
    };
  }, []); // Empty deps array - only run on mount

  // Set up the update cycle
  useEffect(() => {
    const updateData = () => {
      const progress = getProgress();
      const newPoint = progress[progress.length - 1];

      // Only update if values have changed
      if (
        newPoint.iteration !== lastValuesRef.current.iteration ||
        newPoint.metric !== lastValuesRef.current.metric
      ) {
        lastValuesRef.current = newPoint;

        setData((current) => {
          const updated = progress;
          return updated.slice(-maxPoints);
        });
      }

      // Schedule next update
      timeoutRef.current = setTimeout(updateData, updateRate);
    };

    // Start the update cycle
    updateData();

    // Cleanup
    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, [getProgress, updateRate, maxPoints]);

  // Update plot with new data
  useEffect(() => {
    if (!plotRef.current || data.length < 2) return;

    const uPlotData = [data.map((d) => d.iteration), data.map((d) => d.metric)];

    plotRef.current.setData(uPlotData);
  }, [data]);

  return (
    <div
      ref={containerRef}
      style={{
        width: "100%",
        height: "96px",
        backgroundColor: "#1e1e1e",
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        overflow: "hidden",
      }}
    />
  );
};
