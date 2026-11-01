# Scene-first Explore Physics

## Accepted correction

The user sees a long SU(2) lab with controls, matrix and four plots and asks to “unshackle” Explore from any notion of a notebook, making interactivity most prominent. This goes beyond the independent launcher described in [SEPARATION.md](SEPARATION.md). A detached worksheet still fails the requested experience.

Explore MUST “unshackle” itself from any notion of a notebook, beyond removing its shell. Interactivity MUST be the most prominent experience. A standalone page with a long SU(2) worksheet of controls, matrix and four plots does not satisfy this requirement. A static page stack of controls and plots MUST NOT be the primary interface.

The default lab MUST present a dominant interactive viewport, a directly manipulable mathematical object such as an SU(2) rotor, and a compact contextual HUD. Map/world locations and spatial tools MUST support immediately usable mathematical interactions. Detailed control/mathematics sheets and the field journal MUST be hidden until requested; dismissing them MUST restore the same scene, parameters, playback state and unsaved work. Contextual readouts MUST be legible and truthful. Precise equations, axes, conventions, assumptions, source references and evidence MUST remain inspectable on request. Numerical checks, exact computation and formal proof MUST remain distinct; animation or a discovery badge proves none of them.

Stardew Valley, SNES Super Mario World, The Legend of Zelda and Super Metroid guide navigation, discovery, tools, atmosphere and interactivity. Super Metroid does not require combat or difficult platforming. No commercial sprites, maps, fonts or sounds may be copied without permission. Available mathematical tools and evidence MUST NOT require traversal, unlocks or game progression. Direct access, keyboard alternatives, touch alternatives to hover/precision dragging, visible focus, reduced motion and user-controlled sound MUST remain available.

## Bounded SU(2) loop

Open the observatory through the map or direct access. Manipulate the rotor in the dominant viewport immediately. See synchronized contextual feedback. Request a control/math sheet to inspect axis/path settings, matrix, complex coordinates, chart or slice; request the journal to compare, save and reopen. Close it and continue from the same scene/state. Detailed sheets are not a default worksheet. Essential chart identity, poles, invalid input and calculation status remain visible when relevant.

Retain the existing mathematical acceptance: SU(2) is S³; charts/slices are labeled, and Bloch S² is a distinct construction. Preserve U versus -U, declared generator/sign/units, reproducible parameter playback, singularities, source links and evidence provenance. Exact symbolic checkpoints must be distinguished from floating-point evaluation and formal proof. A visual rotor is an interaction tool, not a replacement for these guarantees.

## Acceptance gates

1. Capture direct-entry and lab-entry views at desktop and touch sizes: the interactive scene dominates, the rotor works immediately, and no static controls/matrix/four-plots stack or open journal is primary.
2. Demonstrate scene manipulation changing real linked mathematical state and contextual readouts. Open detailed controls, equations/evidence and journal only on request; dismiss each without losing parameters, playback state or unsaved work.
3. Complete manipulation, pause/step/reset, evidence inspection, save and reopen using keyboard only and touch. Check focus, readable text/zoom, reduced motion and muted feedback. No hover-only or precision movement dependency.
4. Reach available tools by map and direct access without traversal or unlocks. Confirm the standalone entry still refreshes without Notebook dependencies.
5. Reopen existing experiments and journal/source/evidence records; preserve Notebook documents and drafts. Unsupported data stays recoverable with an explanation. Recheck numerical/exact/proof labels and existing SU(2) geometry/path acceptance.

These are required observations, not executed results. Prior separation/browser checks do not verify the new layout or interactions. Parent implementation owns prompt/runtime changes and executed evidence. Human usability iteration, native packaging and a full game remain unverified; no broader game completion is claimed.

## Story ownership and planning receipt

Both [Explore](../../../../user_stories/story__mathematical-physics-exploration.md) and [Notebook](../../../../user_stories/story__unified-notebook-workspace.md) carry this correction. Notebook’s established cell model, themes and Mathematica resemblance preserve that existing experience only; they do not constrain Explore. Shared Rust services/domain data remain allowed.

This bounded documentation pass changes only those two stories and this document. No runtime, prompt mappings, generated contracts or parent-owned tasks.json changes are included. No commits are made. The stories are independent human-authoritative acceptance intent, not generated contracts or substitute implementation prompts.

Read-only PDD planning ran against TypeScript/physics-workbench before these edits, with exact request captured in `/tmp/mathphysics-scene-story-plan-request.txt`, stdout in `/tmp/mathphysics-scene-story-plan.json` and diagnostics in `/tmp/mathphysics-scene-story-plan.stderr`. Intent ID: `update-relevant-mathphysics-userstories-only-plu-a9a0d7c8`. Planning exited 0; stdout includes CLI text surrounding the JSON and a failed read-only debug-snapshot write. No intent apply or provider-backed regeneration is claimed. Planner extraction omitted the compressed constraints and proposed generic runtime targets; the explicit user scope and full preserved request govern this documentation correction.

## Verified first scene-first slice — 2026-10-31

The default laboratory now opens a dominant original observatory scene with a manipulable rotor and contextual tool dock. Drag the ring or use Left/Right/Home/End to rotate. Slice grips the plane interior, supports drag and Up/Down; change axes in the scene HUD. Play, Reset and Keep discovery act directly without a control form. Shape/color distinguish the a0 hemisphere; a compact legend identifies the needle as the Bloch direction, not the whole SU(2)=S3. Rust remains authoritative for coordinates/invariants. Pending/error state dims the device and is never marked checked.

Precision controls, equations, matrix and four plots are now inside native modal Scan, not the default page. Journal and future atlas open only on request. Field guide is also an overlay. Escape returns focus; live calculation/save feedback moves inside the active sheet and returns to the scene on close. This verifies structural placement and keyboard focus, not exhaustive screen-reader announcements.

The mobile layout fits the viewport with every tool visible. All 24 maintained browser tests pass, including 16 exploration/scene checks and the original8 notebook/PDF regressions. They cover scene-dominant defaults, direct real-Rust angle/slice/axis changes without Scan, Keep discovery, optional overlays/Escape, numerical signs, slow requests, persistence failures, source independence and preserved notebook behavior. Strict TypeScript/production build pass. Shared Rust math was unchanged from the previously tested8HTTP checks.

Independent Chrome QA verified actual touch slice-plane drag (c:0→0.56748…), expanded rotate rim touch (t:0→1.23806…), compact legend, visible mobile tool dock at390×844, and Scan/journal/atlas live-region placement/focus. Screenshots include desktop scene and mobile scene. No new page errors were observed. Native touch-tap activation of modal buttons is not fully qualified; the recorded modal checks used click/keyboard.

The actual run_explore.sh launcher at its root URL passed an additional smoke: the user’s existing5.021rad state was retained, direct4pi manipulation worked, discovery save/reload and journal/Escape worked without a visible worksheet or runtime errors. Same exploration schema/key remains; no notebook data or saved discoveries were removed.

Durable evidence and source hashes: `/home/propdev/.openclaw/workspace/outputs/mathphysics-platform/20261031-explore-scene/`. Verification servers are stopped at handoff. To load this update, stop/restart the running Explore process with `./run_explore.sh` and refresh http://127.0.0.1:8877.

The implementation is manual prompt-first work with explicit scene/layout output mappings. No provider-backed PDDapply, automatic regeneration or sync is claimed. Human playtesting, richer environments/challenges, general QFT simulations, GPU rendering and native desktop packages remain open. The updated stories require ongoing game-first iteration, not merely recoloring the former worksheet.

Latest interaction/art iteration: [playful contemporary pixel polish](PLAYFUL-PIXEL.md). Reference libraries are inspiration, not constraints.
