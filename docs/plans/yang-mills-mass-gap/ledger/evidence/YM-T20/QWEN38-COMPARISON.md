# Qwen 3.8 Max versus immutable primary models

Generated UTC: 2026-11-26T08:00:46.074902+00:00. All Qwen observations are fully resolved against the same 43-case set, two samples per case.
No expected label was supplied to any model. This is a separate comparison; the completed primary report was not overwritten.

| Metric | Qwen3.8-max | Kimi3.1-pro | DeepSeek4.1-free | GLM5.3 |
|---|---:|---:|---:|---:|
| Terminal responses /86 | 86 | 86 | 74 | 86 |
| Transport failures /86 | 0 | 0 | 12 | 0 |
| Primary success / terminal | 74/86 | 78/86 | 60/74 | 27/86 |
| Primary success / all 86 | 74 | 78 | 60 | 27 |
| Parseable verdicts | 83 | 84 | 65 | 28 |
| Success given parseable | 74/83 | 78/84 | 60/65 | 27/28 |
| Originals passed /28 | 26 | 27 | 16 | 9 |
| Original false alarms | 1 | 0 | 2 | 0 |
| Exact planted steps /58 | 48 | 51 | 44 | 18 |
| Mutants passed | 7 | 5 | 1 | 1 |
| Wrong/missing step | 1 | 1 | 2 | 0 |
| Unparsed | 3 | 2 | 9 | 58 |
| Token truncations | 3 | 1 | 9 | 58 |
| Requests | 86 | 86 | 132 | 87 |
| Retries | 0 | 0 | 46 | 1 |
| Timeouts | 0 | 0 | 57 | 1 |
| Input tokens | 50570 | 72392 | 19008 | 35530 |
| Output tokens | 689070 | 505715 | 615867 | 1227044 |
| Reasoning tokens | 678779 | 498328 | 611374 | 1222605 |
| Mean seconds per terminal sample | 26.65 | 39.48 | 192.11 | 159.07 |
| Estimated cost USD | $2.123748 | $1.561932 | 0 | $3.005188 |

## Secondary matched-observation views

These subsets condition on all four routes producing terminal or parseable responses; they are post-hoc, selected views, not independent-case superiority estimates.

- **all_terminal:** 74 observations. Successes: Qwen3.8-max 65, Kimi3.1-pro 66, DeepSeek4.1-free 60, GLM5.3 27.

- **all_parseable:** 28 observations. Successes: Qwen3.8-max 28, Kimi3.1-pro 28, DeepSeek4.1-free 28, GLM5.3 27.

## Evidence boundary

Exact per-case results are in `qwen38-case-results.csv` and `qwen38-comparison.json`. All raw Qwen requests and failures remain separately hashed in `QWEN38-RECEIPTS.sha256`. Usage and costs are Pi/TokenRouter estimates, not invoices; requests lacking usage cannot contribute unavailable token totals. Truncated or unparsed outputs, transport failures and mathematical failures remain distinct.
These mechanical historical labels are not independent human certification. No agreement among models, successful compilation or finite-volume numerical result establishes a continuum Yang–Mills mass-gap proof.
