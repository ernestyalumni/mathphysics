# YM-T16 — Weight-space edit of the local 9B referee

Author: strategist (Claude Opus 5.5), 2026-11-18 UTC. Branch `ym/YM-T16-claude-opus55-wsp`. Code: claw-work `Rust/weight-space/` on branch `feat/weight-space-discovery` (a user of this repository is assumed to have claw-work as well).

## Result

**Held-out (2026-11-19, 43 cases × 2 samples per arm, 0 protocol errors): no difference established.** The edit made the 9B reason longer and stop passing planted errors as CORRECT. The gain in errors located is not significant, and the cost is twice as many answers cut off at the token limit.

| | control | edited | prompt |
|---|---|---|---|
| Located or passed (primary) | 39/86 | 46/86 | 48/86 |
| Mutants located | 29/58 | 37/58 | 30/58 |
| Mutants passed as CORRECT | 7 | 1 | 10 |
| Originals passed | 10/28 | 9/28 | 18/28 |
| Cut off at 12,288 tokens | 12 | 25 | 7 |
| Mean reasoning tokens, answers not cut off | 5,615 | 8,104 | 4,854 |

- **Edited against control:**
  - Primary: 17 attempts rescued, 10 regressed (McNemar p = 0.25).
  - Mutants located: 15 rescued, 7 regressed (p = 0.13).
  - By case: 9 cases better, 3 worse (sign test p = 0.15).
- **Prompt against control:**
  - Primary: 14 rescued, 5 regressed (p = 0.064).
  - The gain is in fewer false alarms on originals. Mutants located did not change.
- **Edited against prompt:** on mutants located, the edit alone succeeded on 13 attempts and the prompt alone on 6 (p = 0.17).

Samples of one case are not independent, so the attempt-level p-values overstate certainty.

Both effects are pointing in the expected direction but small. As a referee for this program, the 9B remains a flagger whose output needs reading, not a checker. Details and the follow-ups are in claw-work `Rust/weight-space/docs/PHASE2.md`. Each follow-up would need its own pre-registered run:
- A higher cap.
- Edit and prompt combined.
- A smaller alpha.

## Discovery
**Gate (discovery run, 2026-11-18, 180 answers, 0 protocol errors): not met as pre-registered.** Within the same mutant, a missed sample was shorter than a located one in 30 of 52 pairs (58%), which is not clearly more than half. Across all mutant samples, missed ones were on average *longer* (5,591 reasoning tokens) than located ones (5,177).

Discovery totals, base model, 3 samples per case:
- Mutants: 58 of 135 located, 92 flagged.
- Originals: 23 of 45 passed, 20 false alarms.
- 16 answers hit the 12,288-token cap and gave no verdict.
- Per mutant case, the planted error was located in 0, 1, 2 and 3 of 3 samples for 13, 12, 14 and 6 cases.

Exploratory, not pre-registered:
- A capped answer counts as missed and is long by construction. Without the 14 capped mutant answers, missed was shorter in 30 of 40 pairs (75%).
- Misses split three ways:
  - "CORRECT" verdict: 29 answers, mean 3,084 reasoning tokens.
  - Flagged at the wrong step: 34 answers, mean 4,972.
  - Capped: 14 answers.
  - Located answers average 5,177 tokens.
- So brevity goes with the model *passing* a mutant, not with it misplacing the error.
- This reading was chosen after seeing the data. Only the held-out run can test it.

**Mutants that are not errors.** The strategist found the first one; the second came from reading the rest of the exponent mutants by hand.
- `ym-04 mut4` (discovery) changes ‖E_H((0,Δ'])ψ‖² = 0 to ‖…‖³ = 0. These say the same thing, so the "mutant" is true.
  - Two of its three samples correctly answered CORRECT.
  - The third, counted as "located", flagged step 4 after 7,071 tokens. That was a false alarm.
- `ym-05 mut2` (discovery) changes r_j² ≤ r_j to r_j³ ≤ r_j under 0 ≤ r_j ≤ 1. This is still true; only the inference to q_n ≤ s_n gets weaker. All three samples answered CORRECT.
- Effect of dropping both:
  - Mutants located in discovery: 57 of 129.
  - Gate: 29 of 50 pairs (58%), so the gate verdict is unchanged.
- Held-out: none of the 11 exponent mutants is of this kind.
  - The syntactic check also catches `ym-09 mut1`. There, y¹ = y² = 0 becomes y¹ = y³ = 0.
  - The ² is a coordinate index, not an exponent, so that mutant is a real error and stays in.
- The generator now skips an exponent change where the power is one whole side of a comparison with 0 (`compared_with_zero` in `bin/ym-micro`). The exported case files were not changed: the protocol was fixed before the run.

The pipeline goes on to capture, discovery of the edit and the held-out arms unchanged. That was pre-registered: the edit's effect on held-out accuracy is the second, separate question.

Before the run, the following existed:

- **Cases** (`bin/ym-micro export-cases`, rules in its docstring): 31 of 35 YM-T14 units have at most 25 steps; 29 of them produced cases. Discovery: 15 units, 45 mutants, 60 cases, SHA-256 `530c25ab54a20733826e418641f53c784e1e311f499ef148f1a9db291a7a97db`. Held-out: 14 units, 29 mutants, 43 cases, SHA-256 `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`. Held-out units include C23 and six proofs of the mono-notes topic (`ym-02`, `ym-03`, `ym-06`, `ym-07`, `ym-08`, `ym-09`), two topic sections and five N1 lemmas. Files are outside git (the N1/N2 text is not ours to publish).
- **Tools** (claw-work, verified there, see its `docs/PHASE2.md`): CUDA C++ dequantization of the 9B's Q5_K/Q6_K tensors equal to ggml's own routine in all 134 million values checked; an audited control/edited checkpoint pair on the real 9B; activation capture over reasoning spans with GPU sums.

## What was measured
1. **Gate.** Within one mutant, are the samples that miss the planted error shorter than those that locate it? If not, "reasons too briefly" is the wrong target and the edit is not expected to help.
2. **Held-out.** Mutants located, mutants flagged, false alarms on originals, for control weights, edited weights (4% of the 896 attention and DeltaNet heads, projected off the short-reasoning direction), and control weights with an instruction to check every step at length.

## Context from the overnight run (not used to choose anything in this experiment)
The strategist session hand-checked six INCORRECT flags that the desktop 9B (Qwythos, a different model from the one edited here) raised on presumed-correct proofs: all six were false alarms (ym-01 step 3, ym-07 step 5, ym-09 step 5, ym-s02 step 14, ym-s09 step 14, N2 squaring step 1), against about 28% recall on planted errors. So the held-out comparison must be read for false alarms on originals as much as for mutants located: an edit that makes the model reason longer may also make it flag more. Of these six, ym-07 and ym-09 are held-out units here; ym-s02 is excluded (50 steps).

## Labels
All model output: `unreviewed (local model output; not evidence)`. Code and protocol: by the author, not reviewed by another family. Originals are presumed correct; a genuine error in one of our proofs would count against the model, and would be worth finding.

## Resume here
The run is complete. A follow-up is a new pre-registered run of claw-work `Rust/weight-space/scripts/referee-experiment.sh` with a new `WSP_DATA` directory.
