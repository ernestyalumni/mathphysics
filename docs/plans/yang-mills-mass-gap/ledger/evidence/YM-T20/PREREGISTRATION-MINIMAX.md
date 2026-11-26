# MiniMax-M3 versus Kimi — full frozen replay preregistration

Frozen before the first MiniMax benchmark call on 2026-11-26.

## Model choice

The configured TokenRouter catalog contains exactly one available MiniMax entry: **`MiniMax-M3`** (display name MiniMax M3, 128K context, 16.4K output limit, catalog input/output prices 0.30/1.20 USD per million tokens below its large-input tier). Thus this is the best and only available MiniMax route; no older, smaller, or high-speed substitute is made. Every MiniMax receipt records its provider-reported model for identity inspection.

## Set and samples

Use the historical 43-case held-out mathematical-referee set (14 original proofs and 29 planted-error mutants, 14 groups) exactly as in the completed Kimi run. Case file SHA-256: `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`. Run each case twice in a fresh text-only Pi process, for **86 MiniMax samples**. Models receive the case prompt but never the expected label. The same archived `ym-micro.parse_verdict`, prompts, labels and planted steps are used unchanged.

Kimi remains the existing immutable 86-receipt run; it is not rerun or selected anew. Compare the existing Kimi outcomes with `moonshotai/kimi-k3.1-pro` against MiniMax, with different launch time/provider load stated as a limitation.

## Protocol

- Provider `tokenrouter`, model `MiniMax-M3`, `--thinking low`, no saved session, tools, MCP, extensions, skills or context files.
- Fixed concurrency two; timeout 240 seconds per request.
- Initial request plus at most two retries, only for transport failure, timeout, or no terminal assistant response.
- Never retry a returned wrong, malformed, or token-truncated answer.
- Keep the provider-default output cap; do not change sampling, route or output limit mid-run.
- Preserve every raw response, usage record and failure as atomic per-attempt receipts. Resume existing receipts and support `STOP` without killing active calls.
- All cases stay in headline results. No case exclusion based on mathematical outcomes.

## Scoring and reporting

Primary success is an original passed or a mutant judged `INCORRECT` at exactly its planted first-wrong step. Separately report originals passed, false alarms, exact mutant locations, mutants incorrectly passed, wrong/missing-step flags, transport exhaustion, timeout requests, unparsed output and token truncation.

Report (1) primary success among terminal responses and (2) an all-scheduled sensitivity denominator where transport/unparsed observations are not successes. Include per-group and exact per-case planted/returned steps, two-sample agreement, request time, tokens and Pi-estimated cost. Compare MiniMax directly with the frozen Kimi totals, including matched/common-terminal and common-parseable sensitivity views if reported as secondary—not as replacements for all-scheduled scores.

Transport reliability and mathematical correctness are distinct. Two samples of a case are not independent proofs. No Kimi or MiniMax verdict certifies a Yang–Mills result; every actual research flag requires independent human/strategist or Lean adjudication.

Artifacts are private under this `full-replay` directory, using route `minimax`, `status-minimax.json`, `minimax-run.log`, `run-metadata-minimax.json` and a separate `MINIMAX-COMPARISON.md`.
