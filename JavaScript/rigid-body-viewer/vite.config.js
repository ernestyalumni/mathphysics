import { defineConfig } from 'vite';
export default defineConfig({
  server: { proxy: { '/api': { target: process.env.RIGID_BODY_API || 'http://127.0.0.1:8765' } } },
  build: { emptyOutDir: true, sourcemap: true },
});
