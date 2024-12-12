import React, { useEffect, useRef } from "react";
import uPlot from "uplot";
import "uplot/dist/uPlot.min.css";

export const ProgressGraph = ({ data }) => {
  const containerRef = useRef(null);
  const plotRef = useRef(null);

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
      hooks: {
        draw: [
          (u) => {
            const ctx = u.ctx;
            if (u.data[1].length > 0) {
              const latestVal = u.data[1][u.data[1].length - 1].toFixed(3);
              ctx.fillStyle = "#fff";
              ctx.font = "12px sans-serif";
              ctx.fillText(`Latest: ${latestVal}`, 13, 25);
            }
          },
        ],
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

  // Update data
  useEffect(() => {
    if (!plotRef.current || !data || data.length < 2) return;

    // Transform and limit data
    const displayData = data.slice(-20000);
    const uPlotData = [
      displayData.map((d) => d.iteration),
      displayData.map((d) => d.metric),
    ];

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
