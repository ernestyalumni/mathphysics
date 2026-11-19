# YM-T17 day runs (2026-11-18/19): larger model, more samples

Label: unreviewed (claude-opus55-strat). Same 54 planted mutants and their unmodified steps as `YM-T17-calibration.md`; 52 of the mutants are real errors (ym-04 mut4 and ym-05 mut2 are true statements). The sandbox image was `cadabra2-ubuntu:24.04` (numpy 1.26.4, scipy 1.11.4) for all runs here.

## Larger model: Qwen3.8-27B (IQ2_M, Mac mini), no thinking

| Over 52 real errors | 27B, 1 test | 9B thinking, 1 test |
|---|---|---|
| Detected | 4 | 12 |
| Missed | 21 | 7 |
| FAIL on the unmodified step (of 54) | 6 | 6 |

The 2-bit 27B without thinking is worse than the 9B with thinking. It is fast (108 steps in 80 minutes) and answers NONE more readily. Thinking at about 9 tokens per second would take roughly 20 minutes per step, so it was not tried.

## More samples: 9B with thinking, 3 tests per step, voting

A step is called FAIL when at least K of its 3 tests fail.

| Over 54 mutants | K = 1 (any) | K = 2 | K = 3 (all) |
|---|---|---|---|
| Detected | 11 | 11 | 8 |
| Missed | 5 | 13 | 18 |
| FAIL on the unmodified step | 21 | 10 | **0** |

Neither non-error mutant is among the K = 3 detections, so K = 3 detects 8 of 52 real errors (15%) with no false alarm on 54 unmodified steps. This is the first setting in which a FAIL deserves a person's time by default. Recall is low, so a PASS still means little.

The scoring bug found on the way: with too few FAILs and no PASS, the score fell back to the first verdict, which could be a FAIL. Fixed in YM-T18.

## False-alarm sweep with the 27B

All 295 steps of the mono-notes units and C23, one test each: 89 PASS, 35 FAIL, 150 NONE, 14 ERROR, 7 CRASH.

I read the 7 FAILs on C23 and ym-01 to ym-08. All are wrong tests:

- operator traces computed in a wrong basis (ym-01 steps 3 and 4, C23 step 1);
- unrelated integrals compared (ym-02 step 2);
- an unrelated construction (ym-06 step 3);
- Monte Carlo means near 0 compared at atol 1e-3 (ym-08 step 2);
- a truncated test (C23 step 4, 2 of 200 failing).

All 35 are being re-tested by the 9B with thinking on 3 seeds; a step stays a lead only if all 3 of its tests fail.

**Result of the re-test: none of the 35 steps has 3 failing tests out of 3.**

| Failing tests per step | 0 of 3 | 1 of 3 | 2 of 3 | 3 of 3 |
|---|---|---|---|---|
| Steps | 16 | 6 | 13 | 0 |

The 13 steps with 2 of 3 failing are about what the calibration predicts at K = 2 (10 false alarms on 54 unmodified steps). Five of them are the Killing-form steps: s01 steps 7, 11, 15 and 16, and ym-01 step 3. I checked those independently in the sandbox with `ymlib`. Over 200 instances with N from 2 to 5, the Killing form computed as a real trace over su(N) in the basis iT_a equals 2N tr(XY), so C = 2N² under ⟨X,Y⟩ = −(1/N) Re tr(XY), as the notes say.

So the 27B sweep found no error in the notes that survives the filter. This is weak evidence of absence: at K = 3 the filter caught only 15% of planted errors.

Raw output: `~/.local/state/mathphysics/ym-numcheck/{calib-27b,calib3-thinking-k3,sweep-27b,retest-27b-fails-9b-k3}/`.
