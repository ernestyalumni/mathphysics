# GLM-5.3 output/reasoning diagnostic — frozen before calls

Date: 2026-11-26. This does not overwrite or re-score the completed 86-sample GLM-5.3 baseline.

## Problem and diagnostic set

The baseline `z-ai/glm-5.3` route produced 58/86 token-truncated, unparsed samples at 16,384 tokens. Select four fixed case prompts on which baseline sample 1 exhausted with `stopReason=length` and no verdict, spanning original and mutant arguments and long/short proof groups:

1. `ym-s04-asymptotic-freedom-and-the-bare-coupling`
2. `ym-06-ym-prop-defect`
3. `ym-06-ym-prop-defect--mut1-coefficient-2-to-3`
4. `ym-07-ym-prop-conditional-gap--mut4-flip-le-to-ge`

Selection is explicitly based on output failure, not mathematical correctness. Four cases and one sample per arm are engineering diagnostics only.

## Fixed arms

Both use fresh text-only Pi calls through `tokenrouter`, model `z-ai/glm-5.3`, no saved session/tools/MCP/extensions/skills/context files, the unchanged case prompt and archived parser, concurrency two, 240-second timeout, and at most three requests only for transport-type failure. Never retry a returned wrong, malformed, or token-truncated response. Keep the provider-default output limit.

- `off_original`: `--thinking off`; original system prompt `Review the mathematical argument and follow the requested answer format.`
- `low_verdict_first`: `--thinking low`; system prompt adds `Budget your analysis so a final answer is always emitted. Emit the required verdict block promptly, with a one-sentence reason; do not return analysis without a verdict.`

Run every arm on every diagnostic case once, subject only to the transport retry rule. Preserve each receipt and usage record. No labels are sent to the model.

## Decision rule

Primary engineering statistic: number of parseable verdicts out of four; unparsed and transport failures do not count. Secondary: exact primary mathematical successes, token truncations, timeouts, output/reasoning tokens, elapsed seconds and Pi-estimated cost. Prefer the arm with more parseable answers; break ties first by mathematical successes then fewer output tokens. Do not infer general mathematical accuracy from this selected set.

If neither arm reaches at least 3/4 parseable answers, test a separately configured explicit larger output budget on the same four cases before any fresh validation. Any later configuration change gets an addendum and separate receipts.

A configuration that completes this diagnostic may be evaluated on a new adversarial set built with frozen mutation provenance. Targeted adversarial validation is not an estimate of performance on arbitrary proofs; original-proof labels require independent adjudication. All actual Yang–Mills flags remain subject to human/strategist or Lean review.
