import { defineConfig } from "@playwright/test";
export default defineConfig({
  testDir: "./tests",
  workers: 1,
  timeout: 45000,
  use: {
    baseURL: process.env.WORKBENCH_URL || "http://127.0.0.1:18876",
    viewport: { width: 1440, height: 1000 },
    launchOptions: {
      executablePath: process.env.CHROME_PATH || "/usr/bin/google-chrome",
      args: ["--no-sandbox", "--disable-dev-shm-usage"],
    },
  },
  outputDir:
    process.env.WORKBENCH_TEST_OUTPUT ||
    "/tmp/mathphysics-workbench-test-results",
  reporter: "list",
});
