# Unified mathematical physics notebook workspace

## Story

As a mathematical physics researcher using mathphysics, I want one attractive, readable notebook workspace that I can open in my main web browser or launch as a dedicated Linux or macOS application, with access to my available AI harnesses and Rust mathematical tools, so that I can study Srednicki, develop derivations, run calculations, coordinate agents and inspect their evidence in one consistent interface.

I want two interchangeable visual experiences: a professional notebook interface that looks and behaves as closely as practical to Mathematica, and a polished 16-bit pixel interface drawing on Stardew Valley, Super Mario World and SNES Legend of Zelda design language, with readable mathematical typography and matching interaction sounds.

## User requirements and acceptance

### One notebook workspace across browser and desktop

These shared frontend, document-model, theme and workflow requirements apply only to the existing Notebook across its browser and desktop versions. They preserve Notebook and MUST NOT constrain Explore’s design or interaction model. Explore Physics is separately launched and governed by its own story; it MUST NOT depend on the Notebook shell JavaScript, CSS or document model.

- The browser experience MUST be the primary version. Linux and macOS MUST have dedicated installable applications with icons, launchable by double-click from their application launcher or Applications folder.
- Browser and desktop MUST share the TypeScript frontend, notebook/document model, mathematical commands, agent connections and evidence presentation. A theme or platform change MUST preserve the notebook and its state.
- The existing Rust-first policy MUST remain: mathematical engines, orchestration and native services in Rust; frontend in TypeScript; justified C++/CUDA C++ and then Python remain fallback choices where needed.
- Users MUST be able to create/open/save notebooks, organize text and mathematical cells, evaluate supported calculations, inspect output, interrupt jobs and reopen their work. Notebook-format migration and tested persistence are release requirements.

### Mathematica style

- The professional theme MUST closely reproduce the familiar Mathematica notebook design language: menu/toolbar hierarchy, a spacious document canvas, input/output cell presentation, cell grouping/brackets, mathematical typesetting and restrained application chrome.
- Mathematical entry, text cells, results, source references and evidence MUST remain easy to read at ordinary desktop sizes. Keyboard evaluation and cell navigation MUST be first-class interactions.
- Supported mathematical syntax and operations MUST be discoverable. Visual resemblance does not itself imply a Wolfram Language runtime or `.nb` compatibility; those would require separately specified implementations if requested.

### Readable 16-bit pixel theme and sound

- The pixel theme MUST apply a coherent SNES-era language to windows, panels, buttons, borders, icons, selection/focus, palettes and interaction feedback, with the requested game references guiding its feel.
- Mathematical expressions and notebook body text MUST retain readability and comparable effective font sizes at the same zoom as the professional theme. Pixel styling MUST NOT force tiny text, distorted mathematical glyphs or loss of dotted indices, accents, subscripts and superscripts.
- Theme changes MUST preserve functionality, navigation, source/evidence visibility and unsaved notebook contents. Visual state must remain understandable in both themes.
- The pixel experience MUST have fitting interaction sounds for actions such as selection/activation, submitting a job, completion, cancellation and errors. Sound MUST have user-controlled enable/mute and volume, respect browser activation requirements, and never be the only indication of a result.
- Fonts, sound, animation and motion preferences MUST persist. Reduced-motion behavior, visible keyboard focus, contrast, zoom and keyboard-only use MUST be checked in both themes.
- The supplied snes.css and CSSwitch references MUST be downloaded/inspected for reuse where their source and permissions support it. Preserve source revisions, hashes and notices. Original or appropriately licensed artwork and sounds can supply the requested design language.

### All named agent harnesses and MCP-style tool access

- Both browser and desktop MUST offer connections to available Claude Code, Codex, Grok Build, Hermes Agent and Pi installations or authorized endpoints. Connection discovery, setup, reconnect, cancellation and errors MUST be visible.
- A shared capability registry and adapters MUST distinguish installed/available, configured, connected, unavailable and unsupported states. The UI MUST NOT display an integration as connected or verified without an actual working call.
- MCP/tool access MUST be supported through the appropriate client/server roles. Use supported CLI/API adapters where a harness does not expose the same MCP capabilities; a label saying MCP is not integration evidence.
- The frontend MUST let users select their harness/agent, inspect its supported tools, submit work, see streamed progress/results and cancel bounded jobs. Multiple-agent work MUST show ownership and state without conflating peer agreement with proof.
- Browser access to local harnesses MUST use a Rust companion/gateway or another explicitly configured service. The same shared services must be reachable through the desktop application; a browser page cannot silently start arbitrary local executables.
- Workspace connection credentials and execution permissions belong to the appropriate native/service boundary. Deliver a practical pairing/connection UX rather than exposing provider secrets in notebook or browser storage.

### Integrated with Srednicki tools and evidence

- The notebook MUST browse the existing chapter/item registry, source references, assumptions and pending work, and call actual Rust capabilities from the shared tool registry.
- Computed expressions, numerical support, exact algebra, replayed certificates, formal proofs, source-review states and inconclusive outcomes MUST remain visibly distinct. Metadata validation MUST NOT receive a mathematical-proof badge.
- A complete first vertical slice MUST open a real Srednicki-linked item, submit a supported calculation through a real connected route, display its actual result and assumptions, retain its provenance, and correctly present invalid/budget-exhausted outcomes in both themes.
- UX and Srednicki engine work MUST run as coordinated parallel streams, with shared contracts and integration gates. A frontend mockup or a source inventory MUST NOT be reported as a finished connected mathematical workstation.

## Delivery direction

Use one TypeScript Notebook workbench across browser and desktop with theme tokens and shared notebook state. This does not require Explore Physics to share that shell, state model or styling. Evaluate Tauri 2 as the Rust-native desktop shell for Linux/macOS, with a Rust companion/gateway supporting the primary browser experience. Keep local compute, job management, tool discovery and harness adapters in reusable services; desktop packaging should reuse those contracts rather than fork the product.

Proposed packages are `TypeScript/physics-workbench`, `Rust/physics-gateway` and a desktop wrapper containing Rust native code. The TypeScript browser workbench and bounded local Rust gateway are now implemented; the desktop wrapper remains a proposed boundary. Desktop bundles, browser/local-service connectivity and every named harness require their own executed acceptance gates, including macOS build/launch verification on a Mac.

## Independent acceptance scenarios

1. Open the same notebook in browser and the packaged desktop app; switch themes mid-edit. Inputs, outputs, source links, selected harness and unsaved changes remain intact.
2. Run the existing bounded Clifford trace tool through an actual adapter using new labels; the notebook shows a computed exact-algebra result with its metric convention and no formal-proof claim.
3. Submit an over-budget job and an unavailable harness request. Both themes display truthful error states; neither fabricates a completed result or a connected agent.
4. Inspect a Srednicki item containing dotted spinor indices, epsilon signs and nested scripts at 100%, 125% and 200% zoom. Both themes remain readable and operable by keyboard.
5. Enable, change and mute pixel-theme sounds after user activation. Actions retain visible feedback while muted, and preferences survive reopening.
6. Install and double-click the Linux and macOS bundles on their target systems, connect to a real available harness, evaluate a cell, cancel another job and reopen the saved notebook.

## Intent and ownership

Source: Ernest's two UX messages on 2026-10-31, preserved in `../docs/plans/rust-physics-platform/ux-request.txt`. This is the third program story, alongside the agent-accessible mathematics tools and complete Srednicki coverage. The root program plan and task graph coordinate all three.

This story is human-authoritative intent; it is not a generated contract. Versioned implementation prompts and generated contracts will link to actual frontend/gateway/desktop units as they are introduced. The implementation status below records the bounded browser delivery; desktop bundles and all-harness sessions remain pending.


## File finding and integrated reading acceptance

- The browser MUST offer approved roots, filename search, path/breadcrumb navigation, closeable document tabs and paired TeX/PDF actions. Favorites MUST survive clearing automatic recent/frequent history.
- PDF reading MUST render real page pixels and selectable text with page navigation and zoom using local PDF.js assets; theme changes MUST preserve the current mounted page. A successful PDF response or iframe URL alone is insufficient evidence.
- Local Markdown/TeX input MUST render safely; script/event/unsafe-URL payloads MUST stay inert. Source/math previews MUST disclose lack of full TeX compilation. Draft edits/export MUST NOT silently overwrite source files.

## Bounded delivery status — 2026-10-31

Browser slice implemented at ../TypeScript/physics-workbench: live Notebook/Pixel styles with shared persisted state, optional sound preferences, notebook cells, file search/favorites/history, paired source/PDF tabs, sanitized local preview/export and real Rust Clifford. Text/math note source editors are behind Edit source; trace inputs remain visible. PDF.js6.4.299 legacy compatibility build uses local assets for canvas/selectable-text pages with page/zoom controls; style changes preserve the page. Exact algebra remains distinct from formal proof; invalid labels/13-label oversize remain structured failures.

Earlier independent QA17/17 predates PDF.js and explicitly does not certify plugin page rendering. The durable browser log at `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/browser-tests.log` records **8 passed (6.4s)**, including nonwhite PDF canvas pixels, expected selectable text, PDF page preservation across styles and a local PDF rendered without uploading its bytes. These are maintained browser acceptance tests; independent PDF review remains pending. Desktop launch/install on Linux/macOS, actual harness sessions/MCP adapters, complete Srednicki item provenance and broader cross-platform acceptance remain pending. The umbrella story is not complete. Frontend PDD specs declare manual implementation; automated regeneration/sync/contracts are unverified.

## Additional exploration experience

Explore MUST “unshackle” itself from any notion of a notebook, beyond removing its shell. Interactivity MUST be the most prominent experience. A standalone page with a long SU(2) worksheet of controls, matrix and four plots does not satisfy this requirement. A static page stack of controls and plots MUST NOT be the primary interface.

The default lab MUST present a dominant interactive viewport, a directly manipulable mathematical object such as an SU(2) rotor, and a compact contextual HUD. Map/world locations and spatial tools MUST support immediately usable mathematical interactions. Detailed control/mathematics sheets and the field journal MUST be hidden until requested; dismissing them MUST restore the same scene, parameters, playback state and unsaved work. Contextual readouts MUST be legible and truthful. Precise equations, axes, conventions, assumptions, source references and evidence MUST remain inspectable on request. Numerical checks, exact computation and formal proof MUST remain distinct; animation or a discovery badge proves none of them.

Stardew Valley, SNES Super Mario World, The Legend of Zelda and Super Metroid guide navigation, discovery, tools, atmosphere and interactivity. Super Metroid does not require combat or difficult platforming. No commercial sprites, maps, fonts or sounds may be copied without permission. Available mathematical tools and evidence MUST NOT require traversal, unlocks or game progression. Direct access, keyboard alternatives, touch alternatives to hover/precision dragging, visible focus, reduced motion and user-controlled sound MUST remain available.

Notebook’s cell model, theme parity and Mathematica resemblance apply only to Notebook and MUST NOT govern Explore. See [scene-first acceptance](../docs/plans/rust-physics-platform/game-exploration/SCENE-FIRST.md).

The [game-like mathematical physics story](story__mathematical-physics-exploration.md) now requires Explore Physics to launch separately with its own game-first map, contextual tools and field journal. Notebook retains its established Notebook/Pixel styles, browser/desktop scope, agent connections and Srednicki requirements. Explore MUST work from direct start through refresh, experiment, evidence inspection and save/reopen without a user ever visiting Notebook or loading its shell JavaScript, CSS or document model. Explore MUST NOT inherit notebook tabs, cells, source-editor metaphors, Mathematica vocabulary, forced notebook workflows or “pin to notebook” controls. No mandatory Notebook button is required.

Shared Rust backend services and domain representations are allowed; future explicit data exchange is optional. Existing experiments, journal entries, notes, saved notebooks and unsaved work MUST be preserved without silent deletion. Earlier shared-shell Explore implementation and test evidence are historical, not acceptance evidence for the standalone target. The standalone entry and launcher are implemented; future RTS mode and broader game scope remain pending.

Standalone separation evidence: [current launch and checks](../docs/plans/rust-physics-platform/game-exploration/SEPARATION.md). The independent game launcher and 21 browser checks now pass; human usability iteration and native desktop packaging remain pending.

The approved scene-first correction remains pending parent runtime verification. Existing standalone launch/browser checks do not certify its viewport, direct manipulation, requested sheets or keyboard/touch accessibility. Existing saved and unsaved state MUST be preserved. The broader game remains iterative and incomplete.

Current scene-first evidence: [implementation and acceptance receipt](../docs/plans/rust-physics-platform/game-exploration/SCENE-FIRST.md). The first interactive scene, direct controls and requested detail sheets are implemented;24 browser checks pass. This does not mark the full game or human usability program complete.

Explore art-direction clarification: its earlier snes.css/CSSwitch/GitHub references are starting inspiration, not binding design libraries or a requirement to inherit Notebook’s styles. The separate Explore story permits a contemporary original pixel direction informed by River City Ransom, Stardew Valley and other useful game interfaces, with interaction clarity taking priority. Notebook’s own styles remain separate.
