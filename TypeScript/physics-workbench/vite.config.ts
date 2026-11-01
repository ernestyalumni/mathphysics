import { defineConfig } from "vite";
import { cpSync } from "node:fs";
import { resolve } from "node:path";
const outDir = process.env.WORKBENCH_DIST || "/tmp/mathphysics-workbench-dist";
const exploreOnly=process.env.WORKBENCH_APP === "explore";
export default defineConfig({
  build: { assetsInlineLimit: 0, outDir, emptyOutDir: true, rollupOptions: { input: exploreOnly ? { explore: resolve("explore.html") } : { notebook: resolve("index.html"), explore: resolve("explore.html") } } },
  plugins: [
    {
      name: "local-pdf-assets",
      closeBundle() {
        if(exploreOnly) return;
        for (const folder of ["cmaps", "standard_fonts", "wasm"])
          cpSync(
            resolve("node_modules/pdfjs-dist", folder),
            resolve(outDir, "pdfjs", folder),
            { recursive: true },
          );
        cpSync(
          resolve("node_modules/pdfjs-dist/LICENSE"),
          resolve(outDir, "pdfjs/LICENSE"),
        );
      },
    },
  ],
});
