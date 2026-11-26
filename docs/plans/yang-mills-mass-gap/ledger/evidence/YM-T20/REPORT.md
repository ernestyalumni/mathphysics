# Frozen hosted-model referee replay: deterministic report

Generated UTC: 2026-11-26T08:01:57.866707+00:00. **Status: IN PROGRESS.**

Primary requested comparison: Kimi, free DeepSeek, and `z-ai/glm-5.3`. FlashX receipts are supplementary and are not substituted for the larger GLM model.

Primary success means an original is passed or a planted error is flagged at exactly the planted step. Mathematical categories below count only terminal responses; transport failures and pending samples are separate.

| Model | Scheduled | Terminal | Transport-failed | Pending/retrying | Primary success / terminal | Pass+exact / all scheduled (sensitivity) | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| deepseek/deepseek-v4.1-flash-1017-free | 86 | 74 | 12 | 0 | 60/74 (81.1%) | 60/86 | 9 |
| moonshotai/kimi-k3.1-pro | 86 | 86 | 0 | 0 | 78/86 (90.7%) | 78/86 | 2 |
| z-ai/glm-5.3-flashx | 86 | 64 | 22 | 0 | 44/64 (68.8%) | 44/86 | 17 |
| z-ai/glm-5.3 | 86 | 86 | 0 | 0 | 27/86 (31.4%) | 27/86 | 58 |
| MiniMax-M3 | 86 | 70 | 2 | 14 | 37/70 (52.9%) | 37/86 | 16 |

## Mathematical outcomes (terminal responses)

| Model | Original passed | Original false alarm | Mutant exact step | Mutant passed | Mutant wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|
| deepseek/deepseek-v4.1-flash-1017-free | 16 | 2 | 44 | 1 | 2 | 9 |
| moonshotai/kimi-k3.1-pro | 27 | 0 | 51 | 5 | 1 | 2 |
| z-ai/glm-5.3-flashx | 19 | 0 | 25 | 3 | 0 | 17 |
| z-ai/glm-5.3 | 9 | 0 | 18 | 1 | 0 | 58 |
| MiniMax-M3 | 11 | 3 | 26 | 8 | 6 | 16 |

## Reliability, repeatability, usage, and cost

| Model | Requests | Retries | Timeouts | Token truncations | Fully observed pairs | Same verdict+step | Both samples successful | Input tokens | Output tokens | Reasoning tokens | Cost USD (Pi estimate) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| deepseek/deepseek-v4.1-flash-1017-free | 132 | 46 | 57 | 9 | 34 | 32 | 28 | 19008 | 615867 | 611374 | $0.000000 |
| moonshotai/kimi-k3.1-pro | 86 | 0 | 0 | 1 | 43 | 39 | 37 | 72392 | 505715 | 498328 | $1.561932 |
| z-ai/glm-5.3-flashx | 156 | 70 | 92 | 17 | 25 | 20 | 18 | 63313 | 660772 | 651778 | $0.849727 |
| z-ai/glm-5.3 | 87 | 1 | 1 | 58 | 43 | 34 | 9 | 35530 | 1227044 | 1222605 | $3.005188 |
| MiniMax-M3 | 80 | 7 | 10 | 28 | 35 | 20 | 14 | 33105 | 725806 | 722704 | $0.883511 |

## deepseek/deepseek-v4.1-flash-1017-free: proof groups

| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| c23-transfer-positivity | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| n1-clusters-01-lem-association | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| n1-clusters-02-lem-edge-correlation-inequalities | 8 | 0 | 0 | 8 | 0 | 0 | 0 |
| n1-clusters-04-lem-joint-association | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| n1-clusters-06-prop-crossings | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| n1-rotation-01-lem-gradient | 5 | 0 | 0 | 4 | 0 | 0 | 1 |
| ym-02-ym-lem-layer | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-03-ym-prop-inverse-square | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-06-ym-prop-defect | 10 | 2 | 0 | 7 | 0 | 0 | 1 |
| ym-07-ym-prop-conditional-gap | 10 | 0 | 2 | 6 | 0 | 2 | 0 |
| ym-08-ym-lem-twist-move | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| ym-09-ym-lem-thooft-class | 7 | 2 | 0 | 3 | 1 | 0 | 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | 10 | 2 | 0 | 8 | 0 | 0 | 0 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling | 3 | 0 | 0 | 2 | 0 | 0 | 1 |

### deepseek/deepseek-v4.1-flash-1017-free: exact per-case outcomes

| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |
|---|---|---:|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | transport_failure / — / — | transport_failure / — / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | transport_failure / — / — | transport_failure / — / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | transport_failure / — / — | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | transport_failure / — / — | mutant_exact_step / INCORRECT / 10 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | transport_failure / — / — | unparsed / UNPARSED / — |
| n1-clusters-04-lem-joint-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | transport_failure / — / — | transport_failure / — / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | mutant_exact_step / INCORRECT / 6 | mutant_exact_step / INCORRECT / 6 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| c23-transfer-positivity | CORRECT | — | transport_failure / — / — | unparsed / UNPARSED / — |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient | CORRECT | — | unparsed / UNPARSED / — | transport_failure / — / — |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | mutant_exact_step / INCORRECT / 4 |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-06-prop-crossings | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap | CORRECT | — | original_false_alarm / INCORRECT / 5 | original_false_alarm / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | mutant_wrong_or_missing_step / INCORRECT / 5 | mutant_wrong_or_missing_step / INCORRECT / 5 |
| ym-02-ym-lem-layer | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 6 |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | transport_failure / — / — | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | mutant_exact_step / INCORRECT / 4 |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | unparsed / UNPARSED / — |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-08-ym-lem-twist-move | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |

## moonshotai/kimi-k3.1-pro: proof groups

| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| c23-transfer-positivity | 4 | 1 | 0 | 2 | 0 | 0 | 1 |
| n1-clusters-01-lem-association | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| n1-clusters-02-lem-edge-correlation-inequalities | 10 | 2 | 0 | 8 | 0 | 0 | 0 |
| n1-clusters-04-lem-joint-association | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| n1-clusters-06-prop-crossings | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| n1-rotation-01-lem-gradient | 6 | 2 | 0 | 3 | 0 | 0 | 1 |
| ym-02-ym-lem-layer | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-03-ym-prop-inverse-square | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-06-ym-prop-defect | 10 | 2 | 0 | 6 | 2 | 0 | 0 |
| ym-07-ym-prop-conditional-gap | 10 | 2 | 0 | 7 | 0 | 1 | 0 |
| ym-08-ym-lem-twist-move | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| ym-09-ym-lem-thooft-class | 8 | 2 | 0 | 3 | 3 | 0 | 0 |
| ym-s01-gauge-group-inner-product-casimir-consta | 10 | 2 | 0 | 8 | 0 | 0 | 0 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling | 10 | 2 | 0 | 8 | 0 | 0 | 0 |

### moonshotai/kimi-k3.1-pro: exact per-case outcomes

| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |
|---|---|---:|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | mutant_exact_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | mutant_exact_step / INCORRECT / 8 | mutant_exact_step / INCORRECT / 8 |
| n1-clusters-04-lem-joint-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | mutant_exact_step / INCORRECT / 6 | mutant_exact_step / INCORRECT / 6 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| c23-transfer-positivity | CORRECT | — | original_passed / CORRECT / none | unparsed / UNPARSED / — |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| n1-rotation-01-lem-gradient | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-06-prop-crossings | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-07-ym-prop-conditional-gap | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | mutant_wrong_or_missing_step / INCORRECT / 5 |
| ym-02-ym-lem-layer | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_exact_step / INCORRECT / 5 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | mutant_exact_step / INCORRECT / 4 |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-08-ym-lem-twist-move | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |

## z-ai/glm-5.3-flashx: proof groups

| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| c23-transfer-positivity | 3 | 2 | 0 | 0 | 0 | 0 | 1 |
| n1-clusters-01-lem-association | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| n1-clusters-02-lem-edge-correlation-inequalities | 10 | 1 | 0 | 7 | 0 | 0 | 2 |
| n1-clusters-04-lem-joint-association | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| n1-clusters-06-prop-crossings | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| n1-rotation-01-lem-gradient | 4 | 0 | 0 | 0 | 0 | 0 | 4 |
| ym-02-ym-lem-layer | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-03-ym-prop-inverse-square | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-06-ym-prop-defect | 8 | 2 | 0 | 4 | 1 | 0 | 1 |
| ym-07-ym-prop-conditional-gap | 6 | 1 | 0 | 1 | 0 | 0 | 4 |
| ym-08-ym-lem-twist-move | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| ym-09-ym-lem-thooft-class | 3 | 2 | 0 | 0 | 1 | 0 | 0 |
| ym-s01-gauge-group-inner-product-casimir-consta | 8 | 1 | 0 | 6 | 1 | 0 | 0 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling | 5 | 0 | 0 | 1 | 0 | 0 | 4 |

### z-ai/glm-5.3-flashx: exact per-case outcomes

| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |
|---|---|---:|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | transport_failure / — / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | transport_failure / — / — | transport_failure / — / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | transport_failure / — / — | transport_failure / — / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-04-lem-joint-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | original_passed / CORRECT / none | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | mutant_exact_step / INCORRECT / 6 | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| c23-transfer-positivity | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | unparsed / UNPARSED / — | transport_failure / — / — |
| n1-rotation-01-lem-gradient | CORRECT | — | unparsed / UNPARSED / — | transport_failure / — / — |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | transport_failure / — / — | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-06-prop-crossings | CORRECT | — | unparsed / UNPARSED / — | transport_failure / — / — |
| ym-07-ym-prop-conditional-gap | CORRECT | — | original_passed / CORRECT / none | transport_failure / — / — |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 3 |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | unparsed / UNPARSED / — | transport_failure / — / — |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | transport_failure / — / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | unparsed / UNPARSED / — | transport_failure / — / — |
| ym-02-ym-lem-layer | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | transport_failure / — / — |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | transport_failure / — / — | mutant_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | transport_failure / — / — | transport_failure / — / — |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | transport_failure / — / — | transport_failure / — / — |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | transport_failure / — / — | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | transport_failure / — / — | mutant_exact_step / INCORRECT / 4 |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | transport_failure / — / — | mutant_passed / CORRECT / none |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | unparsed / UNPARSED / — |
| ym-08-ym-lem-twist-move | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |

## z-ai/glm-5.3: proof groups

| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| c23-transfer-positivity | 4 | 0 | 0 | 0 | 0 | 0 | 4 |
| n1-clusters-01-lem-association | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| n1-clusters-02-lem-edge-correlation-inequalities | 10 | 0 | 0 | 0 | 0 | 0 | 10 |
| n1-clusters-04-lem-joint-association | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| n1-clusters-06-prop-crossings | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| n1-rotation-01-lem-gradient | 6 | 0 | 0 | 1 | 0 | 0 | 5 |
| ym-02-ym-lem-layer | 4 | 1 | 0 | 2 | 0 | 0 | 1 |
| ym-03-ym-prop-inverse-square | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-06-ym-prop-defect | 10 | 0 | 0 | 2 | 0 | 0 | 8 |
| ym-07-ym-prop-conditional-gap | 10 | 0 | 0 | 2 | 0 | 0 | 8 |
| ym-08-ym-lem-twist-move | 2 | 2 | 0 | 0 | 0 | 0 | 0 |
| ym-09-ym-lem-thooft-class | 8 | 0 | 0 | 0 | 0 | 0 | 8 |
| ym-s01-gauge-group-inner-product-casimir-consta | 10 | 2 | 0 | 6 | 1 | 0 | 1 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling | 10 | 0 | 0 | 1 | 0 | 0 | 9 |

### z-ai/glm-5.3: exact per-case outcomes

| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |
|---|---|---:|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-04-lem-joint-association | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| c23-transfer-positivity | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | mutant_exact_step / INCORRECT / 4 | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-clusters-06-prop-crossings | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 2 |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-02-ym-lem-layer | CORRECT | — | unparsed / UNPARSED / — | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_exact_step / INCORRECT / 17 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_passed / CORRECT / none |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | unparsed / UNPARSED / — |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-06-ym-prop-defect | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 3 |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | unparsed / UNPARSED / — |
| ym-08-ym-lem-twist-move | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |

## MiniMax-M3: proof groups

| Group | Terminal | Original passed | False alarms | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|
| c23-transfer-positivity | 2 | 0 | 0 | 0 | 0 | 1 | 1 |
| n1-clusters-01-lem-association | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| n1-clusters-02-lem-edge-correlation-inequalities | 10 | 2 | 0 | 6 | 1 | 1 | 0 |
| n1-clusters-04-lem-joint-association | 2 | 0 | 1 | 0 | 0 | 0 | 1 |
| n1-clusters-06-prop-crossings | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| n1-rotation-01-lem-gradient | 6 | 1 | 0 | 1 | 0 | 0 | 4 |
| ym-02-ym-lem-layer | 4 | 2 | 0 | 0 | 2 | 0 | 0 |
| ym-03-ym-prop-inverse-square | 4 | 2 | 0 | 2 | 0 | 0 | 0 |
| ym-06-ym-prop-defect | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ym-07-ym-prop-conditional-gap | 10 | 1 | 1 | 5 | 2 | 1 | 0 |
| ym-08-ym-lem-twist-move | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ym-09-ym-lem-thooft-class | 6 | 0 | 0 | 3 | 1 | 0 | 2 |
| ym-s01-gauge-group-inner-product-casimir-consta | 10 | 1 | 1 | 4 | 2 | 2 | 0 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling | 10 | 0 | 0 | 3 | 0 | 1 | 6 |

### MiniMax-M3: exact per-case outcomes

| Case | Expected | Planted step | Sample 1: category / verdict / returned step | Sample 2: category / verdict / returned step |
|---|---|---:|---|---|
| ym-s04-asymptotic-freedom-and-the-bare-coupling | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut11-coefficient-2-to-3 | INCORRECT | 11 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut7-coefficient-2-to-3 | INCORRECT | 7 | mutant_exact_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 7 |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut10-exponent-2-to-3 | INCORRECT | 10 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-s04-asymptotic-freedom-and-the-bare-coupling--mut8-exponent-2-to-3 | INCORRECT | 8 | mutant_wrong_or_missing_step / INCORRECT / 7 | mutant_exact_step / INCORRECT / 8 |
| n1-clusters-04-lem-joint-association | CORRECT | — | original_false_alarm / INCORRECT / 7 | unparsed / UNPARSED / none |
| n1-clusters-02-lem-edge-correlation-inequalities | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-02-lem-edge-correlation-inequalities--mut9-coefficient-2-to-3 | INCORRECT | 9 | mutant_exact_step / INCORRECT / 9 | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut6-coefficient-2-to-3 | INCORRECT | 6 | mutant_wrong_or_missing_step / INCORRECT / 7 | mutant_passed / CORRECT / — |
| n1-clusters-02-lem-edge-correlation-inequalities--mut15-coefficient-2-to-3 | INCORRECT | 15 | mutant_exact_step / INCORRECT / 15 | mutant_exact_step / INCORRECT / 15 |
| n1-clusters-02-lem-edge-correlation-inequalities--mut10-exponent-2-to-3 | INCORRECT | 10 | mutant_exact_step / INCORRECT / 10 | mutant_exact_step / INCORRECT / 10 |
| c23-transfer-positivity | CORRECT | — | transport_failure / — / — | transport_failure / — / — |
| c23-transfer-positivity--mut3-coefficient-2-to-3 | INCORRECT | 3 | mutant_wrong_or_missing_step / INCORRECT / 4 | unparsed / UNPARSED / 3 |
| n1-rotation-01-lem-gradient | CORRECT | — | original_passed / CORRECT / none | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut4-exponent-2-to-3 | INCORRECT | 4 | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| n1-rotation-01-lem-gradient--mut9-exponent-2-to-3 | INCORRECT | 9 | unparsed / UNPARSED / — | mutant_exact_step / INCORRECT / 9 |
| n1-clusters-06-prop-crossings | CORRECT | — | unparsed / UNPARSED / — | unparsed / UNPARSED / — |
| ym-07-ym-prop-conditional-gap | CORRECT | — | original_passed / CORRECT / none | original_false_alarm / INCORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut3-flip-le-to-ge | INCORRECT | 3 | mutant_exact_step / INCORRECT / 3 | mutant_exact_step / INCORRECT / 3 |
| ym-07-ym-prop-conditional-gap--mut5-flip-le-to-ge | INCORRECT | 5 | mutant_exact_step / INCORRECT / 5 | mutant_passed / CORRECT / 5 |
| ym-07-ym-prop-conditional-gap--mut2-exponent-2-to-3 | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_exact_step / INCORRECT / 2 |
| ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge | INCORRECT | 4 | mutant_passed / CORRECT / 5 | mutant_wrong_or_missing_step / INCORRECT / 5 |
| ym-02-ym-lem-layer | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-02-ym-lem-layer--mut2-flip-le-to-ge | INCORRECT | 2 | mutant_passed / CORRECT / none | mutant_passed / CORRECT / 2 |
| ym-03-ym-prop-inverse-square | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| ym-03-ym-prop-inverse-square--mut1-flip-le-to-ge | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-s01-gauge-group-inner-product-casimir-consta | CORRECT | — | original_passed / CORRECT / none | original_false_alarm / INCORRECT / 13 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut17-coefficient-2-to-3 | INCORRECT | 17 | mutant_wrong_or_missing_step / INCORRECT / 12 | mutant_exact_step / INCORRECT / 17 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut5-exponent-2-to-3 | INCORRECT | 5 | mutant_passed / CORRECT / 5 | mutant_passed / CORRECT / 5 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut11-coefficient-2-to-3 | INCORRECT | 11 | mutant_exact_step / INCORRECT / 11 | mutant_exact_step / INCORRECT / 11 |
| ym-s01-gauge-group-inner-product-casimir-consta--mut2-flip-ge-to-le | INCORRECT | 2 | mutant_exact_step / INCORRECT / 2 | mutant_wrong_or_missing_step / INCORRECT / 12 |
| n1-clusters-01-lem-association | CORRECT | — | original_passed / CORRECT / none | original_passed / CORRECT / none |
| n1-clusters-01-lem-association--mut1-coefficient-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class | CORRECT | — | unparsed / UNPARSED / none | unparsed / UNPARSED / — |
| ym-09-ym-lem-thooft-class--mut6-exponent-4-to-2 | INCORRECT | 6 | mutant_exact_step / INCORRECT / 6 | mutant_passed / CORRECT / none |
| ym-09-ym-lem-thooft-class--mut1-exponent-2-to-3 | INCORRECT | 1 | mutant_exact_step / INCORRECT / 1 | mutant_exact_step / INCORRECT / 1 |
| ym-09-ym-lem-thooft-class--mut4-exponent-4-to-2 | INCORRECT | 4 | retry_pending / — / — | pending / — / — |
| ym-06-ym-prop-defect | CORRECT | — | pending / — / — | pending / — / — |
| ym-06-ym-prop-defect--mut3-flip-ge-to-le | INCORRECT | 3 | pending / — / — | pending / — / — |
| ym-06-ym-prop-defect--mut4-coefficient-2-to-3 | INCORRECT | 4 | pending / — / — | pending / — / — |
| ym-06-ym-prop-defect--mut1-coefficient-2-to-3 | INCORRECT | 1 | pending / — / — | pending / — / — |
| ym-06-ym-prop-defect--mut2-exponent-4-to-2 | INCORRECT | 2 | pending / — / — | pending / — / — |
| ym-08-ym-lem-twist-move | CORRECT | — | pending / — / — | pending / — / — |
