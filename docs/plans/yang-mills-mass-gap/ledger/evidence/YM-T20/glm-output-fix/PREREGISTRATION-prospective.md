# Prospective GLM validation on frozen adversarial candidates

Prepared before any call on the fresh set and before selection of the GLM diagnostic configuration.

## Admission gate

No inference on these cases is authorized unless `../selection.json` has status `selected`. If no original arm reaches ≥3/4 parseable answers and the separately preregistered explicit-budget arm also fails the gate, the launch service stops without calling the candidate set. The exact selected arm and its diagnostic hashes are recorded at startup and matched on resume; no post-selection switching is permitted.

A larger-client-budget selection does not prove a provider cap increase. It uses exactly the independently configured 32768-token requested output budget from `../larger-output/PREREGISTRATION.md`.

## Frozen inputs and calls

- `../../../adversarial-20261126/frozen-candidates.json`, SHA-256 `c156555e93303974c1a40f4dc141e412c68f680bba3e5df8c6a08ed747dc23b4`.
- All 21 cases: four originals and 17 planted-error mutants; four groups. Prompts, IDs, aliases and mechanical expected labels are not rewritten after a model answer.
- Two independent text-only calls per case, 42 scheduled terminal samples. The two calls of one case are not independent proofs.
- Fresh text-only Pi processes; no saved sessions/tools/MCP/extensions/skills/context files.
- Exactly the selected GLM thinking and system instruction and its output-budget configuration.
- Concurrency two; per-request timeout 240 seconds; at most three requests only for transport-type failure. Never retry a returned wrong, malformed, or token-truncated response.
- Archive the same `ym-micro.parse_verdict`; no manual reclassification of terminal model answers.
- Preserve every request receipt, stop reason, usage and timing. A scoped `STOP` prevents new calls without killing in-flight requests.

## Reported outcomes

Primary candidate-label success means original `CORRECT` or mutant `INCORRECT` at exactly its mechanically planted first-wrong step. Report terminal and all-scheduled denominators, original pass/false-alarm counts, mutant exact-step/wrong-or-missing-step/passed counts, unparsed/token-truncated outcomes, timeouts/retries/exhausted transport, request timing/tokens/cost, defect/group breakdowns and all per-case returned steps.

The freeze contains augmented historical sources and new self-contained material. Candidate labels are **not independently human-adjudicated**. This prospective result is targeted engineering evidence, not a randomly sampled proof-population estimate. Model agreement does not fill the blind human worksheet. Any later human adjudication is an addendum/sensitivity analysis and will preserve these outcomes rather than delete disputed cases.
