# Standalone Explore Physics — superseding correction

The user requested Explore be launched separately, leave notebook interactions behind, and use new game-first UX rather than Mathematica language. They explicitly asked to update all relevant stories. Both story__mathematical-physics-exploration.md and story__unified-notebook-workspace.md now carry that requirement. This supersedes the shared-shell integration in the first implementation, without deleting its saved work.

Run `./run_explore.sh`, then open http://127.0.0.1:8877. `--port N` overrides it. Explore has its own entry (explore.html, explore-main.ts), game base stylesheet, navigation, sound preference, field guide and journal. The dedicated build excludes notebook entry, KaTeX, PDF.js and notebook styles/model; the game runs at the root URL and refreshes there. It shares Rust APIs and mathematical data structures, not notebook UI/state. The small read-only field guide can attach and open references independently. Discovery export replaces notebook pinning.

Notebook remains `./run_workbench.sh` at8876 with its original Notebook/Pixel styles and saved data. An optional external link opens /explore.html independently; no Explore tab or embedded game exists in Notebook. An ordinary workbench build provides two separately bundled HTML entries for convenience. The dedicated Explore launch builds only the game. Browser app launching is delivered; native desktop bundles remain pending.

Existing exploration key mathphysics.explore.v1 and journal/convention/source records remain compatible. Existing notebook storage/notes are neither read nor changed by Explore. Browser storage is origin-specific: a new port has separate storage. To move a pre-separation journey from8876, open http://127.0.0.1:8876/explore.html directly with the workbench server running, export the exploration, then import it at8877. This does not load or require interacting with Notebook. Keep the export; no old storage is deleted. Exported discovery entries and full journey exports serve different purposes; full-journey import accepts the versioned exploration format.

Exact correction request and read-only intent plan are adjacent. The planner and manual implementation do not claim provider-backed apply, automatic regeneration or sync. Bounded prompt mappings register standalone entry/styles/launcher and remove the notebook dependency on exploration. Mathematical services are unchanged; new browser tests cover standalone isolation plus retained math and journal behavior. See the latest separation evidence receipt for actual results. Earlier shared-shell screenshots/test claims are historical and do not satisfy this correction.

This is the independent foundation, not the final game design: the walkable observatory, contextual lab controls and field journal continue to iterate. RTS, wider QFT simulations, GPU renderer and native packages are still future work. Feynman/scattering/Standard Model destination planning remains applicable.

## Verified receipt

21 browser acceptance tests pass. The actual run_explore.sh game-only root launch, refresh, Rust checkpoint and discovery reload pass. Evidence and screenshots: `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-explore-standalone`. Shared math was unchanged; the preceding eight Rust HTTP checks remain relevant to that unchanged kernel. The old 20-check shared-shell run is historical.

The subsequent scene-first redesign supersedes the default worksheet. See [current scene-first evidence](SCENE-FIRST.md).
