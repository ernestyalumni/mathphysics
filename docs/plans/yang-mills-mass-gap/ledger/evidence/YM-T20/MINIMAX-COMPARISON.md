# MiniMax-M3 versus Kimi — frozen 43-case replay

Generated UTC: 2026-11-26T08:18:02.874502+00:00. Each route has 86 scheduled samples; no expected label is supplied to the model. MiniMax is the only configured MiniMax route. MiniMax started later, so provider load and exact provider-side output caps are not controlled.

Primary success means an original passed or a mutant flagged at exactly its planted step. Mathematical rates count terminal responses; all-scheduled sensitivity also retains transport failures.

| Metric | MiniMax-M3 | Kimi |
|---|---:|---:|
| Terminal samples | 82 | 86 |
| Transport-failed samples | 4 | 0 |
| Requests / retries | 97 / 11 | 86 / 0 |
| Timeouts | 15 | 0 |
| Unparsed terminal answers | 16 | 2 |
| Token truncations | 30 | 1 |
| Primary success / terminal | 46/82 | 78/86 |
| Success / all 86 scheduled | 46 | 78 |
| Parseable verdicts | 66 | 84 |
| Success given parseable | 46 | 78 |
| Originals passed /28 | 15 | 27 |
| Original false alarms | 3 | 0 |
| Exact planted step /58 | 31 | 51 |
| Mutants passed | 11 | 5 |
| Wrong/missing step | 6 | 1 |
| Same verdict+step pairs | 23 | 39 |
| Both-sample-success cases | 17 | 37 |
| Mean seconds per terminal sample | 119.81 | 39.48 |
| Input tokens, all requests | 40605 | 72392 |
| Output tokens, all requests | 817725 | 505715 |
| Reasoning tokens, all requests | 813749 | 498328 |
| Cost estimate USD | $0.996576 | $1.561932 |
| Cost estimate per success USD | $0.021665 | $0.020025 |

## Matched terminal observations

Secondary, post-hoc subset: 82 observations; MiniMax 46/82, Kimi 76/82.
Discordances are recorded in `minimax-comparison.json`; these post-hoc common-response samples do not support an independent-case superiority claim.


## Matched parseable observations

Secondary, post-hoc subset: 66 observations; MiniMax 46/66, Kimi 61/66.
Discordances are recorded in `minimax-comparison.json`; these post-hoc common-response samples do not support an independent-case superiority claim.

## Error and research boundary

Receipt-level mathematical errors remain planted-label outcomes, not proof that each historical original has independent human certification. Exact per-case categories and planted/returned steps are reproduced below and in `samples.csv`. Different elapsed times are cumulative request time, not guaranteed wall-clock model latency. Usage/cost is Pi/TokenRouter accounting, not a provider invoice.

A double model pass is not a theorem certificate. Use every actual Yang–Mills flag only as input to independent strategist/human or Lean adjudication.

## Exact MiniMax per-case comparison

| Case | Expected | Planted step | MiniMax sample 1: category / verdict / step | MiniMax sample 2 | Kimi sample 1 | Kimi sample 2 |
|---|---|---:|---|---|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | unparsed / UNPARSED / — | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | mutant_exact_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | unparsed / UNPARSED / — | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | mutant_wrong_or_missing_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 8 | mutant_exact_step / INCORRECT / 8 | mutant_exact_step / INCORRECT / 8 |
| n1-clusters-04-lem-joint-association | CORRECT | — | original_false_alarm / INCORRECT / 7 | unparsed / UNPARSED / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | mutant_wrong_or_missing_step / INCORRECT / 7 | mutant_passed / CORRECT / — | mutant_exact_step / INCORRECT / 6 | mutant_exact_step / INCORRECT / 6 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| c23-transfer-positivity | CORRECT | — | transport_failure / — / — | transport_failure / — / — | original_passed / CORRECT / none | unparsed / UNPARSED / — |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | mutant_wrong_or_missing_step / INCORRECT / 4 | unparsed / UNPARSED / 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| n1-rotation-01-lem-gradient | CORRECT | — | original_passed / CORRECT / none | unparsed / UNPARSED / — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | unparsed / UNPARSED / — | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 4 | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-06-prop-crossings | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-07-ym-prop-conditional-gap | CORRECT | — | original_passed / CORRECT / none | original_false_alarm / INCORRECT / 5 | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_passed / CORRECT / 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | mutant_passed / CORRECT / 5 | mutant_wrong_or_missing_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 4 | mutant_wrong_or_missing_step / INCORRECT / 5 |
| ym-02-ym-lem-layer | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | original_false_alarm / INCORRECT / 13 | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_wrong_or_missing_step / INCORRECT / 12 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | mutant_passed / CORRECT / 5 | mutant_passed / CORRECT / 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_wrong_or_missing_step / INCORRECT / 12 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | unparsed / UNPARSED / none | unparsed / UNPARSED / — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | mutant_exact_step / INCORRECT / 6 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / none | mutant_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | transport_failure / — / — | transport_failure / — / — | mutant_exact_step / INCORRECT / 4 | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | mutant_passed / CORRECT / 4 | mutant_exact_step / INCORRECT / 4 | mutant_exact_step / INCORRECT / 4 | mutant_exact_step / INCORRECT / 4 |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_passed / CORRECT / none | mutant_exact_step / INCORRECT / 1 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_passed / CORRECT / 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-08-ym-lem-twist-move | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none | original_passed / CORRECT / none |
