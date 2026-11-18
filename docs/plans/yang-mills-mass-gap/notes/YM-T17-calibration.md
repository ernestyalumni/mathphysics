# YM-T17 calibration: a 9B model as a numeric-test writer

Label: unreviewed (claude-opus55-strat, 2026-11-18). Model: Qwen3.8-9B-Distill Q5_K_M on the laptop (RTX 3070). Units: the 54 planted mutants from the 2026-11-16 overnight set (mono-notes units and C23), each paired with the same step of the unmodified unit. One test per step (seed 0), 200 random instances per test.

## Result

| Outcome over 54 mutants | No thinking | Thinking | Thinking, scipy in sandbox |
|---|---|---|---|
| Detected (mutant FAIL, original PASS) | 6 | 13 | 13 |
| Missed (mutant PASS) | 6 | 8 | 8 |
| FAIL on the unmodified step | 6 | 5 | 6 |
| Untestable (NONE, NOCODE, ERROR, CRASH on either side) | 36 | 28 | 27 |
| Model time per step, median | 30 s | about 60 s | (re-run of saved tests) |

For comparison, the 9B referee located 94 of 333 planted errors (28%) in the overnight run.

**Correction for non-errors.** Two of the 54 mutants are true statements: ym-04 mut4 (‖·‖² = 0 → ‖·‖³ = 0) and ym-05 mut2 (r_j² ≤ r_j → r_j³ ≤ r_j for 0 ≤ r_j ≤ 1). Both were confirmed by q7nd, who fixed the generator in YM-T14 fe6963f. Without them there are 52 real errors:

- No thinking: both were untestable, so the result stays 6 detected and 6 missed.
- Thinking: 12 detected and 7 missed. The ym-04 mut4 "detection" becomes a false alarm on a true statement, and the ym-05 mut2 "miss" is a correct PASS.

## Every FAIL on an unmodified step was read

All 12 (6 per run) are wrong tests, not errors in the notes:

- The test computes one side wrongly: ym-01 step 2, summing all entries where a trace is meant, in both runs.
- An inequality is tested as an equality: ym-05 step 3 (`≥ 0` tested as `= 0`), ym-07 step 2 (requires `c₁δ₀ = ½`).
- Quadrature or Monte Carlo is held to rtol 1e-5 to 1e-7: ym-03 step 1, s04 step 11, s06 step 3.
- The instance does not satisfy the hypotheses: s05 step 22 (takes `M_A = 100` without deriving it), s07 steps 9 and 23, s01 step 5 (also asserts the norm is 0).
- The lattice geometry is wrong: s10 step 19, ym-09 step 6.

## Detections are weaker than the count

Spot check of 6 of the 13 thinking-run detections:

- Genuine: s01 mut17 (3N = 2N), ym-02 mut2 (wrong indicator in Tonelli), ym-01 mut5.
- Plausibly genuine, on an arbitrary instance: ym-03 mut1.
- Right verdict from a wrong test: ym-04 mut1 (`≤` tested as equality), ym-07 mut4 (no-thinking run; same flaw).
- Not an error at all: ym-04 mut4 changes `‖E_H((0,Δ'])ψ‖² = 0` to `‖…‖³ = 0`, which is equivalent. The mutant is a true statement, so its FAIL is a false alarm. YM-T14 now skips such exponent changes (fe6963f). ym-05 mut2 is also true, and no text rule catches it.

## Harness defects found and fixed

1. A reply with an opening code fence but no closing one was run with the fence included (SyntaxError, reported as CRASH).
2. A test that returns `None` (truncated, or the model decided NONE inside its own code) counted as FAIL. It now counts as an error.
3. numpy 2 removed `np.trapz` and `np.math`; models also assume `np.linalg.expm`. The runner supplies these. `expm` is checked against `eigh` on skew-Hermitian matrices.
4. With thinking, the model imports `scipy.special` for Bessel functions (C23). `YM_NUMCHECK_IMAGE` selects an image with scipy; `resandbox` re-runs saved tests without asking the model again.

## Conclusion

The 9B is not a reliable test writer. Thinking roughly doubles the nominal detections (6 to 13 of 54), but the FAILs on correct steps stay at about one per detection, and they are 0 of 12 real. A FAIL here is not evidence about the mathematics until a person reads the test.

What is better than the referee: every verdict comes with 15 to 40 lines of code that a person can check in about a minute, and the false alarms are visibly wrong (an inequality tested as an equality, a Monte Carlo held to 1e-7). The referee's flags were prose that took a careful re-derivation to dismiss.

The mechanical check moves the problem from "is the step right" to "is the test right"; with a 9B writer the second question is as hard to trust as the first. Next steps, in order of cost:

1. Same 108 steps with a stronger writer (the Mac's 27B, or Codex) to see whether test quality scales with the model.
2. Give the writer verified samplers (su(N), SU(N), Haar, spectral measures from the Rust `lie-groups` crate) so that it only writes the two sides of the claim, not the instance.
3. Do not run the 295-step false-alarm sweep with the 9B: the 54 paired originals already show the rate (6 of 54, all wrong tests).

Raw output: `~/.local/state/mathphysics/ym-numcheck/{calib1,calib2-thinking,calib2-thinking-scipy}/` (outside the repo; `score.txt` in each).
