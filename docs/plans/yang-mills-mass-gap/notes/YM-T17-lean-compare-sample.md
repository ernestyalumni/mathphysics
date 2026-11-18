# Precision of 9B lean-compare verdicts: a 10-item sample

Label: unreviewed (claude-opus55-strat, 2026-11-18). Checked at the level of the statement; definitions were not unfolded.

## Result

Qwythos 9B returned 100 INCONSISTENT verdicts over 851 comparisons of OpenAI Lean declarations with one-line blueprint summaries:

- desktop: 72 of 633
- laptop: 28 of 218

I drew 10 at random (Python `random.Random(20261118)`, shuffle, first 10) and read each declaration against its summary. **None was a real mismatch.** With 0 of 10, the 95% upper bound on precision is about 26%.

| K-id | Declaration | Model's reason, shortened | Finding |
|---|---|---|---|
| K-III.6.2 | `YM.Work.F06.answer` | proof only uses `rw`, `exact hh` | Judged the helper lemma above the target; `answer : YM.Bid.F06` is exactly the summary |
| K-III.4.56 | `YM.RG08.projecting` | no inverse is defined | Statement gives `treeMatrix * (Z N Zᵀ) = 1`, with N the inverse of the joint completion; the docstring says the same |
| K-III.10.60 | `…Pane.ligatures` | asserts ∃ρ with conditions, not "comparisons" | Summary is truncated and vague; hypothesis `2 ≤ c` matches; no conflict |
| K-VII.2.27 | `YM.Trace.Loft.awake` | no trace quantity mentioned | `Tendsto (xup q …).re L (𝓝 (landPay G))`: a quantity converging to the constant, under opaque names |
| K-III.10.72 | `…Coat.pressure_mean` | "at the stationary mean" not stated | Bound is on `X.realPath … X.μ`, with `Rate … (-1)` |
| K-III.4.111 | `…Folio.takings` | no load or exponential smallness | `Enough` unfolds to a hit bound with `rate (-prize v) (radiusB x)`, eventually in x |
| K-V.2.11 | `YM.Work.C22.answer` | only `funext` | Judged the helper lemma `tagOrigins` above the target |
| K-III.11.54 | `…Soft.s_choir` | no "stamp" | Conclusion is `Choir … (short B N p) X` with `n ≤ 2`: the short stamp |
| K-III.7.12 | `YM.Stationary.tribute` | no Jacobian ratio | det > 0 near 0 and an identity, in file `HaarJac.lean`; plausible, not unfolded |
| K-VII.1.10 | `YM.Pert.Warp.earth` | no legs, no separated configuration | `Tendsto eqPrice → cut z …` with `z ∈ conf r` |

## Failure modes

1. **Vocabulary (7 of 10).** The OpenAI Lean uses opaque names, and the model treats "the summary's words do not occur in the statement" as a mismatch.
2. **Wrong declaration (2 of 10).** The snippet starts a few lines before the target, and the model judged the preceding helper lemma. This is a harness fix in `ym-micro` (YM-T14).
3. **Method instead of statement (1 of 10).** The summary names the proof method, and the model faulted the code for not showing it.

## Consequence

Stop 9B lean-compare as an audit signal. A real comparison of statement against blueprint needs the definitions unfolded, so it is a T3 job (YM-A03, YM-A05), not a 9B job.
