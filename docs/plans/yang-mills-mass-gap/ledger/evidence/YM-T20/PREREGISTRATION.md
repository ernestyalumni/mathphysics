# Hosted Yang–Mills referee full replay — fixed before full-run requests

Date: 2026-11-25 Pacific. This extends, but does not replace, the four-case pilot in `../PILOT.md`.

## Question

How do three hosted TokenRouter models perform on the historical held-out mathematical-referee cases? Report mathematical performance separately from transport and response-format reliability. This is a benchmark of referee behavior, not evidence that any argument is correct and not a mass-gap proof.

## Frozen inputs and routes

- Cases: `../frozen-test.json`, 43 cases from 14 proof groups (14 originals, 29 planted-error mutants), SHA-256 `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`.
- The discovery set remains unused for this replay. The current exporter's 42-case output is not substituted for the historical 43-case set.
- Routes:
  - `deepseek`: `deepseek/deepseek-v4.1-flash-1017-free`;
  - `kimi`: `moonshotai/kimi-k3.1-pro`;
  - `flashx`: `z-ai/glm-5.3-flashx`.
- FlashX was chosen before the full run because it is the configured Flash route named by Ernest; on the small pilot its two returned responses were both correct (one original passed and one exact error location), while a third request timed out. The unavailable exact ID `z-ai/glm-5.3-flash` is not silently substituted. The larger `z-ai/glm-5.3` is not part of this run.

## Scheduled observations

- Two fresh, independently invoked samples for every case and every route: 86 scheduled samples per route, 258 total. Pilot receipts are not reused as full-run samples.
- Every call is a fresh text-only Pi invocation: no tools, MCP, extensions, skills, context files, or saved session; provider `tokenrouter`; requested thinking level `low`; system prompt `Review the mathematical argument and follow the requested answer format.`
- The case prompt, expected label, and planted step are exactly those in the frozen JSON. The model receives only the case prompt and system prompt, never the label.
- No sampling temperature or random seed is invented for routes whose Pi interface does not expose a pinned seed. Report this limitation.
- Maximum wall time per request: 240 seconds. A timeout kills that request's local process group.
- A scheduled sample receives at most three requests (initial plus two retries). Retries are allowed only when no terminal assistant response was obtained: local timeout, nonzero/no-message invocation, or assistant `stopReason=error`. Backoff is deterministic and recorded. A returned answer, including an unparseable or token-truncated one, is not retried. A wrong verdict or wrong step is never retried.
- If all allowed requests fail at transport, the scheduled sample is scored as a transport failure, not a mathematical mistake. Every request receives a separate immutable receipt.
- Concurrency is fixed by route: DeepSeek 8, Kimi 2, FlashX 2. Other routes do not consume DeepSeek's queue. Cases within each route run in frozen-file order subject to concurrency; sample 1 and sample 2 jobs are adjacent but independent.
- A `STOP` file prevents new requests and lets active requests finish. Existing receipts are resumed, not overwritten.

## Outcomes and analysis

The archived `ym-micro.parse_verdict` parser is used unchanged. For each scheduled sample with a terminal assistant response:

- Original passed: verdict `CORRECT`.
- Original false alarm: verdict `INCORRECT` (the asserted step is also reported).
- Mutant located: verdict `INCORRECT` and first wrong step exactly equals the planted step.
- Mutant passed as correct: verdict `CORRECT`.
- Mutant flagged at wrong step: `INCORRECT` with another step or no parsed step.
- Unparsed and token-truncated responses are reported separately and are not called mathematical errors.

For each route report all scheduled samples, terminal responses, transport failures, the above counts, two-sample agreement, per-proof-group outcomes, elapsed time, token usage, and Pi-estimated cost where returned. Also give results among terminal responses so transport reliability does not masquerade as mathematical performance. Count primary located-or-passed only when an original passes or a mutant is located exactly. In a sensitivity view, count transport/unparsed observations as not located-or-passed.

Historical Qwen/Phi figures remain context, not a perfectly controlled comparison: their inference settings, output caps, and local execution differed. No significance claim will treat two samples of one case as independent proofs. Any exclusions for underspecified or disputed original cases will be shown as a sensitivity analysis, with all-case results retained.

## Evidence and privacy

- All prompts, assistant text, per-request outcomes, timing, usage, run metadata, and summaries stay under this private directory.
- No credentials, private machine addresses, or private source text are copied into the public `mathphysics` repository.
- No upstream report, publication, or claim of a mass-gap proof follows from this benchmark.
