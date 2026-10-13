import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// GitHub Pages serves the app under the repository name; Vercel and local runs serve from the root.
const GITHUB_PAGES_BASE = process.env.GITHUB_PAGES_BASE || "/Assured_QC/";

export default defineConfig(({ command }) => ({
  plugins: [react()],
  base: command === "build" && process.env.GITHUB_PAGES === "1" ? GITHUB_PAGES_BASE : "/",
  worker: { format: "es" },
  server: { host: true, port: 5174 },
}));
