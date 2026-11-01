# Unified notebook UX delivery stream

The unified browser/Linux/macOS notebook is the program's third story. Its two themes share one notebook model and backend contracts; all named harnesses and MCP-style tools are integration targets. This stream runs alongside the Rust engines and Srednicki item work, using their actual capability and evidence records.

## Shared architecture

```text
TypeScript browser workbench        Linux/macOS desktop workbench
             |                         TypeScript + Rust shell
             +------------------------+
                          |
              Rust workspace gateway/companion
              notebook/jobs/capabilities/evidence
                /                     \
     Rust mathematics tools       Harness adapters
                              Claude Code / Codex / Grok Build
                                  Hermes Agent / Pi
                          |
            Appropriate MCP client/server or CLI/API routes
```

The gateway owns local process access and workspace connections. Browser and desktop use the same job and evidence protocol, with actual capability discovery. Existing `physics-tools` is a bounded JSON Clifford bridge, not yet a gateway, MCP server or general notebook language.

Tauri 2 is a candidate desktop shell because it combines a shared web frontend with Rust native application code and supports Linux/macOS distribution. Do not claim macOS readiness from a Linux-only build. A browser-only deployment cannot independently spawn local CLI harnesses; the companion connection is a real product requirement.

## Parallel scheduling

| Stream | Initial assignments | Dependency and integration gate |
|---|---|---|
| Mathematical engines | F06/F07 errors and exact coefficients; F08/F09 expressions; F10 replay | Stable typed capability/result/evidence contracts |
| Srednicki coverage | Refine CH003, independently review CH035, continue chapter inventories and bounded executable batches | Source-reviewed items, assumptions, negative controls |
| Notebook UX | UX01 interaction/format contracts; UX02 professional theme; UX03 pixel theme; UX04 sound/accessibility | Same notebook state, readable math, no invented result statuses |
| Services and harnesses | UX05 Rust gateway; UX06 adapters and live job states | Actual calls and cancellation; truthful availability states |
| Packaging | UX07 shared Tauri desktop shell; UX08 Linux/macOS packaging | Real target-OS launch/install tests |
| End-to-end verification | UX09 both-theme browser/native gates; UX10 Srednicki notebook release | Real sourced item → actual calculation → result/evidence persistence |

Task IDs and detailed acceptance live in the central tasks.json. UX specifications and CSS reference auditing can proceed alongside coefficient/expression work. The first connected notebook slice can use the already available Clifford operation without waiting for full CAS coverage, while clearly showing its current limits.

## References and assets

- User reference: https://github.com/devMiguelCarrero/snes.css/
- User reference: https://remiscan.fr/csswitch/
- Notebook reference: https://www.wolfram.com/mathematica/
- Candidate desktop architecture: https://v2.tauri.app/
- Desktop packaging: https://v2.tauri.app/distribute/macos-application-bundle/ and https://v2.tauri.app/distribute/appimage/

Download/audit records preserve exact source revisions and reusable files separately from product code. Use pixel chrome and expressive sound while retaining clear mathematical fonts. A source's existence does not establish that its assets, fonts or audio have all been reviewed for reuse or browser accessibility.

## Current implementation boundary — 2026-10-31

The first browser slice is implemented at TypeScript/physics-workbench with a local Rust gateway: live shared-state Notebook/Pixel styles, persisted notebook and sound preferences, file search/path navigation, favorites/history, paired TeX/PDF tabs, safe local previews/export and actual bounded Clifford calls. Root launch contract is ./run_workbench.sh → http://127.0.0.1:8876. Source+math preview is not TeX compilation.

Integrated PDF.js6.4.299 legacy canvas/selectable-text rendering replaces the headless blank native plugin, with local assets, page navigation/zoom and page preservation across styles. The durable browser log at `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/browser-tests.log` records **8 passed (6.4s)**, including nonwhite PDF canvas pixels, expected selectable text, PDF page preservation across styles and a local PDF rendered without uploading its bytes. These are maintained browser acceptance tests; independent PDF review remains pending. Earlier independent17/17 QA covered the preceding slice and cannot certify the new PDF renderer. The six bounded frontend prompts and nested architecture/config describe manual implementation, not successful automatic PDD regeneration.

UX01–04 and UX09 have bounded browser implementation/evidence, not complete cross-platform acceptance. UX05 provides a local file/Clifford gateway slice; UX06 live harness sessions/MCP adapters, UX07–08 desktop shell/target-OS bundles and UX10 complete Srednicki provenance workflow remain pending. Harness installation availability is explicitly not connectivity. General CAS, formal proof and full TeX compilation are not delivered. Do not mark the umbrella complete. Details: ../../../TypeScript/physics-workbench/docs/IMPLEMENTATION.md (repository-relative frontend evidence summary).

Final independent PDF QA: the exact 665-page Srednicki textbook rendered its cover and page 10 (198 selectable text spans); next/previous, 125% zoom and Pixel switching passed. No page errors, console errors, failed requests or CSP violations were observed. Durable report/screenshots: `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-workbench/pdf-review/`. CSP is supplied in HTML meta, not an HTTP header. The real repository launcher also built and served the integrated frontend and gateway; all build output stayed outside the repository.
