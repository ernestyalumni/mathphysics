# YM-T20: hosted referees on the frozen held-out cases

**Result (2026-11-26, primary routes complete).** On the 43 held-out cases used for the local referees, Kimi K3.1 Pro located 51 of 58 planted errors at the exact step and raised no false alarm on 28 correct originals. The best local referees located 28–29 of 58. Every one of Kimi's seven misses is a planted change that leaves the conclusion intact, or a case whose label is in doubt (below). DeepSeek V4.1 Flash (free) was as accurate whenever it returned an answer, but 12 of its 86 samples never returned one. GLM-5.3 was accurate whenever it answered, but 58 of its 86 answers stopped at the 16,384-token output limit while still reasoning.

Numbers come from `ledger/evidence/YM-T20/samples.csv` and `summary.json` (primary table: `PRIMARY-REPORT.md`). The protocol was fixed before the full run: `PREREGISTRATION.md`, plus `PREREGISTRATION-GLM53.md` for the larger GLM route.

## Outcomes

Two samples of each case, 86 per model. The local columns are the YM-T16 control (Qwen3.8-9B Q5) and the Phi-4-reasoning-plus regrade of YM-F01 experiment 3, on the same file.

| | Kimi K3.1 Pro | DeepSeek V4.1 Flash free | GLM-5.3 | Phi-4 (local) | Qwen3.8-9B (local) |
|---|---:|---:|---:|---:|---:|
| Primary (original passed or error at exact step), all 86 | 78 | 60 | 27 | 53 | 39 |
| Originals passed / 28 | 27 | 16 | 9 | 25 | 10 |
| False alarms on originals | 0 | 2 | 0 | 2 | 11 |
| Planted errors located at exact step / 58 | 51 | 44 | 18 | 28 | 29 |
| Mutants passed as CORRECT | 5 | 1 | 1 | — | 7 |
| Mutants flagged at another step | 1 | 2 | 0 | — | — |
| Cut off at output limit | 1 | 9 | 58 | 7 | 12 |
| Other unparsed | 1 | 0 | 0 | — | — |
| No response after three attempts | 0 | 12 | 0 | 0 | 0 |
| Agreement of the two samples (verdict and step), pairs observed | 39/43 | 32/34 | 34/43 | — | — |
| Cost (Pi estimate) | $1.56 | $0 | $3.01 | local | local |

Among answers that returned a verdict: Kimi 78/84, DeepSeek 60/65, GLM-5.3 27/28.

The settings were not the same: the local runs used temperature 0.6 and a 12,288-token cap with fixed seeds, while the hosted routes used the provider defaults, with `--thinking low` requested and no seed. So the comparison holds only at the size of the gap: 51 against 29 planted errors located, with zero false alarms against 2 and 11.

## What the misses are

Every miss that returned a verdict, checked against the case text:

1. **`ym-07` (conditional gap), original. The label is doubtful.** Step 5 ends with "since c₁ ≥ 2", and the statement as frozen does not contain c₁ ≥ 2. The recurrence hypothesis only names c₁ and sets η\* = 1/(2c₁). Without the bound the step fails. For example, at c₁ = 1 and k = 0 the allowed δ₀ = 1/2 gives (1/n)·log(1/(e^{1/2} − 1)) ≈ 0.43/n < (log 2)/n. The current mono-notes fix c₁ = 3740, but the frozen case predates that. DeepSeek flagged step 5 for this reason on both samples. These "false alarms" are correct readings of the text the referee was given. The 9B flag at `ym-07` step 5, which the YM-T16 note counts as a hand-checked false alarm, should be read again against the frozen case text.
2. **`ym-07` mut4 (`0 ≤ x` changed to `0 ≥ x` in step 4).** Kimi and DeepSeek (three samples between them) said step 5. Kimi read the step-4 change as a typo and then gave the same c₁ ≥ 2 reason. This case is ambiguous for the same reason as item 1.
3. **`ym-06` mut1 (`Z(2n)` changed to `Z(3n)` in step 1).** Kimi passed both samples. The changed formula is false, but step 1's conclusion, 0 < P_k ≤ 1, is still drawn from the correct bound. This is a real miss of a false line.
4. **`ym-09` mut4 and mut6 (coordinate index `y⁴` changed to `y²`).** Kimi passed three of its four samples and located mut4 once. DeepSeek passed mut4 once. The generator's "exponent 4 → 2" rule hit a superscript coordinate index, not an exponent, so the mutant reads like a typo.
5. **`ym-s01` mut5 (`‖X‖²` changed to `‖X‖³`).** GLM-5.3 called it "the only blemish… typographical… affects nothing downstream" and answered CORRECT.

The pattern is the same throughout. The hosted models judge whether the argument reaches its conclusion and forgive local slips that do not propagate. The benchmark counts any planted change as an error at that step. Leaving out the four `ym-07` samples per model, Kimi scores 75/80 among parsed answers, DeepSeek 60/61 and GLM-5.3 27/28. Every remaining Kimi miss is of kind 3 or 4.

## Reliability

- **GLM-5.3:** 1.22 M of 1.23 M output tokens were reasoning. The median answer used the full 16,384 tokens and took 168 s. The route does not honour `--thinking low`. To be usable here it needs a higher output limit or a reasoning-effort control that works.
- **DeepSeek free:** 57 timeouts in 132 requests. Seven of the 12 lost samples are `ym-s04`, the longest case: 3 of its 10 samples returned. The free route is the bottleneck, not the model.
- **Kimi:** no timeouts, one cut-off (`c23`, original), median 4.7 k output tokens and 32 s per answer.
- **GLM-5.3-FlashX** (supplementary, still running at this writing): 20 of 20 returned verdicts correct so far, with many timeouts and cut-offs. Its final row is in `REPORT.md`.

## Consequences for the program

- Kimi is the first referee here whose CORRECT verdicts carry some weight: 0 false alarms, and 51/58 located. It still passes false lines whose conclusion survives, so it does not replace a line-by-line check (`ym-numcheck`, Lean).
- Before the next referee comparison, the case set needs two fixes. The `ym-07` statement must carry c₁ ≥ 2 (or c₁ = 3740). The mutant generator must stop rewriting coordinate indices as exponents. The referee instruction should also say whether a false line that does not affect the conclusion counts as an error. These are case-set changes, so they need a new frozen file and hash, not edits to this one.
