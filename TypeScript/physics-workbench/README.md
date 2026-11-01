# Mathematical physics workbench

A local browser notebook and playable mathematical physics observatory, with Notebook/Pixel styles, file finding and real Rust mathematical tools.

From the mathphysics repository root:

```sh
./run_workbench.sh
```

Open http://127.0.0.1:8876. Node.js 24 is recommended; Rust/Cargo is required for the local gateway. The launcher snapshots the frontend into an external build cache and defaults to `npm ci --offline --ignore-scripts --no-audit --no-fund`; Rust builds also use `--offline --locked`. Cached npm/Cargo dependencies are required. For an explicitly authorized online npm bootstrap, use `WORKBENCH_NPM_OFFLINE=0 ./run_workbench.sh`. Alternatively set `WORKBENCH_NODE_MODULES` to matching installed dependencies. The npm opt-in does not enable online Cargo resolution. Use `--port N` to override the port or `--build-only` to build without serving.

Use roots, paths and search to find files; pin favorites or revisit recent/frequent folders. Clear visit history keeps favorites. Open paired TeX/PDF, switch among closeable tabs, edit drafts and export without overwriting source files. Local .md/.tex files use sanitized previews. Source+math preview supports common equations/headings, **not full TeX compilation**; unsupported document commands remain visible source.

PDF.js 6.4.299 uses its legacy compatibility build, local worker/font/CMap/WASM assets, canvas pages and selectable text. Page navigation and PDF zoom are available. Switching styles retains the current page and mounted PDF. This is integrated rendering rather than reliance on the browser's native PDF plugin. Assets remain same-origin; build.assetsInlineLimit=0 prevents inline font data URIs under strict font-src self.

Notebook text/math cells default to displayed content with Edit source controls; trace inputs remain visible. Shift+Enter runs actual Rust Clifford calculations. The two-label `mu, nu` result is `4*h(2:mu;2:nu)` in the formal symmetric h convention; exact algebra is not formal proof. Invalid labels and over-budget requests remain errors. Styles, drafts and sound mute/volume persist locally. Sound is optional and user activated. Ctrl/Cmd+K finds files, Escape closes the mobile drawer, Ctrl/Cmd+S exports.

Desktop Linux/macOS packaging, live harness sessions/MCP adapters, full CAS, complete Srednicki item workflows and formal proof are pending. Listed available harnesses are not connected sessions. Local PDFs must be reopened after browser restart.

See [implementation and evidence](docs/IMPLEMENTATION.md), [PDD scope](PDD.md) and [third-party notices](THIRD_PARTY_NOTICES.md).


## Reproduce browser acceptance tests

The normal launcher serves mathphysics and an optional books root. Browser tests additionally require an explicit synthetic `fixtures` root. From the repository root, create fixtures and build:

```sh
node TypeScript/physics-workbench/tests/create-fixtures.mjs /tmp/mathphysics-workbench-fixtures
./run_workbench.sh --build-only
```

The build prints `Frontend:` and `Gateway:` paths. In another terminal, substitute those exact paths below and launch the gateway directly with the fixture root:

```sh
/path/printed/Gateway --static-dir /path/printed/Frontend --root mathphysics="$PWD" --root fixtures=/tmp/mathphysics-workbench-fixtures --port 8876
```

With matching frontend dependencies available, run from `TypeScript/physics-workbench`:

```sh
WORKBENCH_URL=http://127.0.0.1:8876 npm run test:browser
```

The generator writes synthetic PDF/TeX/Markdown fixtures only; no source-corpus documents are required. The Playwright configuration uses `/usr/bin/google-chrome` by default (`CHROME_PATH` overrides it). A normal launcher without the `fixtures` root is insufficient for these tests.

The durable browser log at `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/browser-tests.log` records **8 passed (6.4s)**, including nonwhite PDF canvas pixels, expected selectable text, PDF page preservation across styles and a local PDF rendered without uploading its bytes. These are maintained browser acceptance tests; independent PDF review remains pending.

## Explore Physics is a separate application

Run `./run_explore.sh` from the repository root and open http://127.0.0.1:8877. It has its own full-window game shell, field guide and journal. There are no notebook tabs, cells, theme switch or notebook pinning. Mathematical notation remains precise. The shared Rust service supplies SU(2); it does not require Notebook to run. `--port N` is supported. The dedicated build contains only the game entry.

Notebook remains at8876 through `./run_workbench.sh`. Its external Explore link opens a separate page at /explore.html. Existing notebook notes are preserved. Existing exploration exports remain compatible; different ports have separate browser storage, so use full journey export/import to move existing discoveries. See [separation and migration](../../docs/plans/rust-physics-platform/game-exploration/SEPARATION.md).

World movement, direct lab entry, linked plots, playback, slices, journal and numerical evidence remain available. Discovery export replaces notebook pinning. Future QFT destinations remain labeled planned; no new amplitude engine is claimed. This is a browser application, not a delivered native desktop bundle.

## Scene-first observatory update

Explore now starts in a game viewport rather than a long lab worksheet. Drag the ring, use the Rotate/Slice tools, change axes, Play, Reset or Keep discovery. Scan opens precise settings, equations and linked plots only when requested; Field journal and World atlas are separate requested overlays. Escape returns to the scene. Restart the running Explore process after source updates, then refresh. See [scene-first acceptance](../../docs/plans/rust-physics-platform/game-exploration/SCENE-FIRST.md).

## Playful pixel polish

The original reference CSS/repos do not constrain Explore’s visual design. Current polish adds warmer original environments, whole-instrument grab cues, labeled pictogram tools and canvas-only R/S/Space/K shortcuts. Optional Tips offers real one-/two-turn experiments; keeping a moving point waits for its exact result. Original save tones stay opt-in. See [current design and evidence](../../docs/plans/rust-physics-platform/game-exploration/PLAYFUL-PIXEL.md).
