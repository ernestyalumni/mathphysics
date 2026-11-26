# Yang–Mills referee replay, 2026-11-25 Pacific

Four historical cases, two original arguments and two planted errors. I selected the layer identity and deficit-telescope arguments before generating responses. Each call received only the original numbered prompt in a fresh text-only Pi session. Verdict and first wrong step were scored against the archived labels. The original arguments are graded according to the historical benchmark. This is a selective pilot, with one completed response per case at most per hosted model. It is unsuitable for a general accuracy estimate.

| Case | Historical label | Kimi K3.1 Pro | DeepSeek V4.1 Flash 1017 free |
|---|---|---|---|
| Layer identity | Correct | Correct | Correct after retry |
| Layer, flipped inequality | Incorrect at 2 | Located step 2 | Located step 2 |
| Defect telescope | Correct | Correct | 90-second timeout |
| Defect, altered coefficient | Incorrect at 4 | Located step 4 | Located step 4 |

Kimi scored 4/4. DeepSeek scored 3/4 submitted cases, 3/3 returned verdicts. Its first layer-identity request also ended in a provider error. The pilot used Pi with low thinking requested, a short referee system prompt and a 90-second per-call timeout. Historical local arms used different model templates, sampling controls and caps. The two local samples per case are repeated observations of each proof. Historical Qwen control scored 2/4 and 3/4 on these same cases. Regraded Phi-4 scored 3/4 and 3/4. No GLM request was issued: the authenticated catalog contains z-ai/glm-5.3-flashx but does not contain z-ai/glm-5.3-flash.

## Evidence and replay

- Frozen test set SHA-256: b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe. Discovery set: 530c25ab54a20733826e418641f53c784e1e311f499ef148f1a9db291a7a97db. The frozen copies here reproduce the old export before exclusion of exponent changes compared to zero. The current exporter against the saved units produces 42 test cases instead of 43. Use frozen-test.json or the original archive for cross-model comparisons.
- Original cases and local output: the YM-T16 held-out case file and its control-arm results; Phi-4 scores from the YM-F01 experiment 3 regrade (all outside git).
- Hosted output: adjacent kimi-*.json and deepseek-*.json receipts. Each records case, route, verdict, step, response or transport outcome, elapsed time and returned usage. The timed-out call has no provider usage record.
- Next comparison: all 43 frozen cases twice per route, including the four pilot cases afresh. Count original passes, false alarms, planted errors located, timeouts, truncation and cost separately. Report paired results by proof group (14 groups). Include a sensitivity analysis for known underspecified originals. Keep any full-corpus external run within the source text's sharing permissions. For Lean proof generation, keep the statements frozen and grade with Lean.
