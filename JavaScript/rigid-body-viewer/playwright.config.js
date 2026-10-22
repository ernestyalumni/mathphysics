import { defineConfig } from '@playwright/test';
export default defineConfig({
  testDir: './tests/browser', workers: 1, timeout: 60000,
  use: {
    baseURL: process.env.SIM_URL || 'http://127.0.0.1:8765',
    viewport: { width: 1600, height: 1000 },
    launchOptions: {
      executablePath: process.env.CHROME_PATH || '/usr/bin/google-chrome',
      args: ['--no-sandbox', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'],
    },
    trace: 'retain-on-failure',
  },
  outputDir: '/tmp/mathphysics-rigid-body-browser-tests',
  reporter: 'list',
});
