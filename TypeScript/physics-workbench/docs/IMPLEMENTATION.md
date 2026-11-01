> Superseded integration: Explore now launches independently. See [standalone separation](../../../docs/plans/rust-physics-platform/game-exploration/SEPARATION.md). Earlier shared-shell behavior/evidence below is historical.

# Browser slice status — 2026-10-31

Implemented in the TypeScript frontend: shared Notebook/Pixel state; autosaved cells/drafts/preferences; optional sound; file root/path/search/breadcrumbs, favorites, recent/frequent folders and history clearing; closeable tabs, paired TeX/PDF; safe local Markdown/TeX preview and export; actual local Rust Clifford trace with explicit exact-algebra/error status; responsive drawer and keyboard shortcuts. Text/math note source editors are collapsed by default; trace inputs stay visible.

Integrated PDF.js6.4.299 legacy compatibility build renders bounded pages into canvas with selectable text, local assets, page navigation/zoom, generation/cancellation guards and disposal. Style changes preserve the mounted page. Source inspection confirms these mechanisms; implementation presence alone is not independent acceptance evidence.

The durable browser log at `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/browser-tests.log` records **8 passed (6.4s)**, including nonwhite PDF canvas pixels, expected selectable text, PDF page preservation across styles and a local PDF rendered without uploading its bytes. These are maintained browser acceptance tests; independent PDF review remains pending. Earlier independent17/17 QA covered the preceding native-viewer slice; its iframe/byte assertions do not establish PDF.js rendering. The maintained tests supply the stronger canvas/text assertions. See ../README.md for synthetic fixture generation and direct fixture-root gateway setup.

Remaining: target-Linux/macOS desktop install/launch, real harness sessions/MCP tool adapters/job cancellation, complete Srednicki source-item-to-result provenance, general CAS/full TeX compilation, broader accessibility/zoom/font fidelity checks and audible sound review. Do not mark the umbrella UX story complete. PDD specifications describe manual implementation; full automated generation/sync/contracts are unverified.

Final independent PDF QA: the exact 665-page Srednicki textbook rendered its cover and page 10 (198 selectable text spans); next/previous, 125% zoom and Pixel switching passed. No page errors, console errors, failed requests or CSP violations were observed. Durable report/screenshots: `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/pdf-review/`. CSP is supplied in HTML meta, not an HTTP header. The real repository launcher also built and served the integrated frontend and gateway; all build output stayed outside the repository.

## Playable Explore extension

The first game-like SU(2) observatory and experiment journal are implemented. The combined browser suite now has 20 passing checks; the gateway has8 HTTP checks. See [Explore evidence](../../../docs/plans/rust-physics-platform/game-exploration/IMPLEMENTATION.md) for independent review, runtime limits and the QFT backlog. Canvas 2D is used for the first chart renderer; no GPU acceleration, RTS, desktop or live harness completion is claimed.

The scene-first redesign supersedes the lab worksheet: see [current evidence](../../../docs/plans/rust-physics-platform/game-exploration/SCENE-FIRST.md).
