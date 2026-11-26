# Requested three-model planted-error referee comparison

Generated UTC: 2026-11-26T06:04:44.542106+00:00. **Status: COMPLETE.**

Routes: Kimi, free DeepSeek, and larger `z-ai/glm-5.3` (not FlashX). Exact per-case planted/returned steps are in `samples.csv` and `REPORT.md`; all raw responses and failures are in `receipts/`. Mathematical outcomes count terminal responses only.

| Model | Terminal /86 | Transport-failed | Pending/retrying | Primary success | Original passed | Original false alarm | Exact mutant step | Mutant passed | Wrong/missing step | Unparsed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| moonshotai/kimi-k3.1-pro | 86 | 0 | 0 | 78/86 (90.7%) | 27 | 0 | 51 | 5 | 1 | 2 |
| deepseek/deepseek-v4.1-flash-1017-free | 74 | 12 | 0 | 60/74 (81.1%) | 16 | 2 | 44 | 1 | 2 | 9 |
| z-ai/glm-5.3 | 86 | 0 | 0 | 27/86 (31.4%) | 9 | 0 | 18 | 1 | 0 | 58 |

| Model | Requests | Retries | Timeouts | Token truncations | Input tokens | Output tokens | Reasoning tokens | Cost USD (Pi estimate) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| moonshotai/kimi-k3.1-pro | 86 | 0 | 0 | 1 | 72392 | 505715 | 498328 | $1.561932 |
| deepseek/deepseek-v4.1-flash-1017-free | 132 | 46 | 57 | 9 | 19008 | 615867 | 611374 | $0.000000 |
| z-ai/glm-5.3 | 87 | 1 | 1 | 58 | 35530 | 1227044 | 1222605 | $3.005188 |
