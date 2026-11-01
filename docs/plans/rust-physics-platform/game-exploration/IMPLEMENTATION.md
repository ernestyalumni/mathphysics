> Superseded integration: Explore now launches independently. See [standalone separation](SEPARATION.md). Earlier shared-shell behavior/evidence below is historical.

# Explore implementation evidence — 2026-10-31

The first playable mathematical physics experience is implemented in the existing browser workbench. Enter Explore from the sidebar/tab, walk or directly enter the SU(2) observatory, manipulate axis/angle/slice, play/pause/step/reset and compare 0/π/2π/4π. Four linked Canvas views show complex α/β, a labeled stereographic chart and a coordinate slice, plus numeric matrix/invariant/Bloch readouts. Camera interaction has pointer and keyboard alternatives. Mathematical values are returned by the Rust POST /api/su2 endpoint, with independently tested Pauli signs, unitary/determinant checks, chart and slice edge cases. U and -U remain distinct. Computed status is numerical, never formal proof.

Experience and style are independent: Explore has its own wider composition and keeps both Notebook and Pixel styles. Notebook drafts/file tabs survive transitions. The journal records inputs/convention/observations/source links, reopens and compares parameters, pins notebook notes and exports/imports versioned state. Ctrl/Cmd+S exports the active experience. Quota/blocked-data failure retains in-memory export with truthful feedback. Invalid/unknown stored state is preserved for explicit recovery. Async requests and animation are disposed on exit; background/hidden worlds skip drawing; BFCache restores the view. Playback paces Rust sampling instead of abort-starving slower responses, and only matching input/results receive checked status.

## Verification

- `npm run build`: strict TypeScript and production build pass, using cached pinned dependencies and external output.
- `WORKBENCH_URL=http://127.0.0.1:18878 npm run test:browser`: **20/20 pass**. Includes original 8 notebook/PDF checks; SU(2) matrix fixtures, sign, slices, custom invalid axis, journal/reload/pin, playback, delayed220ms response pacing, scrub during playback, export shortcut, storage quota and unknown-version recovery, mobile layout and style/draft compatibility.
- Rust `cargo test --offline --locked`: **8/8 HTTP tests pass**, including existing gateway cases and independently constructed Pauli/matrix/Bloch fixtures, unitary/determinant, noncommuting signs, mixed axes, chart poles/near-poles, invalid inputs, body/method/origin protections. `cargo fmt --check` and Clippy all-targets -D warnings pass.
- Independent browser QA: **23/23 functional + 13/13 targeted layout/lifecycle checks pass**, with source attachment/reopening and notebook provenance checked separately. Final screenshots cover desktop/mobile, Notebook/Pixel and world/lab. No runtime/CSP/HTTP errors observed; readback performance warnings come from test pixel inspection.
- Actual repository `./run_workbench.sh --port 18879` builds and serves from external caches. Launcher smoke verifies world, real 2π=-I calculation and Pixel switching.

Durable evidence: `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-explore/`. Source hashes and logs/screenshots live there. Local verification servers are stopped at handoff; rerun ./run_workbench.sh and open http://127.0.0.1:8876.

## Limits and next work

Canvas 2D renders the first 3D projections; no GPU compute or CUDA acceleration is claimed. This is a small playable lab, not a completed game world. Direct manipulation is currently through axes/time/slice controls and the plot camera; broader complex-analysis functions need their own specified domains. Journal comparison shows parameters, not a computed group distance. Source links are file references, not immutable item-level Srednicki provenance. Target-user usability sessions and broader accessibility/audio checks remain pending.

Chapter-linked QX01–QX07 in QFT-DESTINATIONS.md cover diagram building, scattering, Standard Model groups, field modes, spinors, RG and gauge/Higgs. Visible destinations are explicitly planned; no amplitude engine or Standard Model simulation is claimed. Future RTS mode, desktop bundles and live harness/MCP adapters remain pending. Existing Srednicki item inventories are not complete book implementation.

## PDD

Approved game story and exact request/read-only planning receipts are retained. Modules have explicit prompts/output mappings, including numerical su2.rs and exploration CSS; tests remain independent. Implementation is manual prompt-first work in the authorized existing checkout. No provider-backed intent apply, full automatic regeneration, synchronization or generated story-contract success is claimed. Original dirty files are preserved; no commit or publication.
