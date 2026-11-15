# YM-A02 — Comparator on Linux with real process isolation; second kernel if possible

Tier T2–T3 · Workstream A · Depends on: YM-A01 · Blocks: G1 · Review: T3 · Budget: one day

## Why
Assessment §9.2: the recorded Comparator runs used `fake-landrun.sh` on macOS, so they "do not test isolation of potentially hostile build code", and "an additional external kernel implementation was not used". We have Linux. Both gaps can be closed here.

## Do
1. Read `YMG/submissions/openai-20261028/verification/` and `docs/rerun.md` for the pinned Comparator and exporter revisions and configurations.
2. Build the pinned Comparator. Check that the kernel supports Landlock (`uname -r` is 6.8; confirm the LSM is enabled) and that real `landrun` is used, not the fake launcher.
3. Run the three recorded targets: `YangMills.jointExistence`, `YMReview.Acceptance.combined_acceptance`, `YMReview.OPE.constructed`. Run the four controls (valid proof accepted; changed statement, forbidden axiom, changed definition rejected).
4. Find out whether any independent Lean 4 kernel (for example `lean4lean` or `nanoda`) supports toolchain 4.32.0 exports. If one does, replay the export for the public theorem. If none does, record what was tried.

## Deliverable
`PLAN/notes/YM-A02-isolated-comparator.md`; logs under `PLAN/ledger/evidence/YM-A02/`.

## Acceptance
For each target: pass or fail with log path, and an explicit line stating that isolation was real, with the evidence. Item 4 has a definite outcome.

## Do not
Report anything upstream. That is Ernest's call after review.
