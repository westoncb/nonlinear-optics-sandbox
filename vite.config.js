import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  server: {
    open: true,
  },
  build: {
    outDir: "dist",
    sourcemap: true,
  },
  plugins: [react()],
  base: "/nonlinear-optics-sandbox/",
});
