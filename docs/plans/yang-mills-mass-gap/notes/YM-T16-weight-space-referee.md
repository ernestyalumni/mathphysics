# YM-T16 — Weight-space edit of the local 9B referee: setup and status

Author: strategist (Claude Opus 5.5), 2026-11-18 UTC. Branch `ym/YM-T16-claude-opus55-wsp`. Code: claw-work `Rust/weight-space/` on branch `feat/weight-space-discovery` (a user of this repository is assumed to have claw-work as well).

## Result so far
No accuracy result yet: the GPU runs on the 9B have not started (the 3060 was held by another server). What exists:

- **Cases** (`bin/ym-micro export-cases`, rules in its docstring): 31 of 35 YM-T14 units have at most 25 steps; 29 of them produced cases. Discovery: 15 units, 45 mutants, 60 cases, SHA-256 `530c25ab54a20733826e418641f53c784e1e311f499ef148f1a9db291a7a97db`. Held-out: 14 units, 29 mutants, 43 cases, SHA-256 `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`. Held-out units include C23 and six proofs of the mono-notes topic (`ym-02`, `ym-03`, `ym-06`, `ym-07`, `ym-08`, `ym-09`), two topic sections and five N1 lemmas. Files are outside git (the N1/N2 text is not ours to publish).
- **Tools** (claw-work, verified there, see its `docs/PHASE2.md`): CUDA C++ dequantization of the 9B's Q5_K/Q6_K tensors equal to ggml's own routine in all 134 million values checked; an audited control/edited checkpoint pair on the real 9B; activation capture over reasoning spans with GPU sums.

## What will be measured
1. **Gate.** Within one mutant, are the samples that miss the planted error shorter than those that locate it? If not, "reasons too briefly" is the wrong target and the edit is not expected to help.
2. **Held-out.** Mutants located, mutants flagged, false alarms on originals, for control weights, edited weights (4% of the 896 attention and DeltaNet heads, projected off the short-reasoning direction), and control weights with an instruction to check every step at length.

## Context from the overnight run (not used to choose anything in this experiment)
The strategist session hand-checked six INCORRECT flags that the desktop 9B (Qwythos, a different model from the one edited here) raised on presumed-correct proofs: all six were false alarms (ym-01 step 3, ym-07 step 5, ym-09 step 5, ym-s02 step 14, ym-s09 step 14, N2 squaring step 1), against about 28% recall on planted errors. So the held-out comparison must be read for false alarms on originals as much as for mutants located: an edit that makes the model reason longer may also make it flag more. Of these six, ym-07 and ym-09 are held-out units here; ym-s02 is excluded (50 steps).

## Labels
All model output: `unreviewed (local model output; not evidence)`. Code and protocol: by the author, not reviewed by another family. Originals are presumed correct; a genuine error in one of our proofs would count against the model, and would be worth finding.

## Resume here
Run claw-work `Rust/weight-space/scripts/referee-experiment.sh` with `WSP_CASES` pointing at the exported cases; then fill in the two sections above from `gate.json`, `compare-control-vs-edited.json` and `compare-control-vs-prompt.json`.
