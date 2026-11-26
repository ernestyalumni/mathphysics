# Qwen 3.8 Max — full frozen replay preregistration

Frozen before the first Qwen mathematical benchmark request, following Ernest's request to add `qwen/qwen3.8-max`.

## Route and set

- Provider `tokenrouter`, route exactly **`qwen/qwen3.8-max`**. The TokenRouter catalog confirms this route with 128K context and a displayed 16.4K output limit. Do not substitute another Qwen model or change the provider output cap.
- Use the historical 43-case referee set: 14 originals and 29 planted-error mutants in 14 proof groups. Case file SHA-256: `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`.
- Run every case twice, **86 scheduled samples**, in fresh text-only Pi processes. Supply only the case prompt, never the expected label or planted step. Preserve the original per-prompt hashes and deterministic `ym-micro.parse_verdict` parser.

## Identical replay operating rules

- System instruction: `Review the mathematical argument and follow the requested answer format.` Request `--thinking low` as in the other full routes; no saved sessions, tools, MCP, extensions, skills or context files.
- Fixed provider-default output cap and concurrency two; 240 seconds per request.
- Initial request plus at most two retries, only for transport failure, timeout or no terminal assistant response. Never selectively retry a returned wrong, malformed or token-truncated response. Do not change prompts, route, sampling, parser, output cap or comparison rules during the run.
- Write an immutable atomic receipt for every request, retaining provider-reported model identity when available, parsed verdict/step, complete returned answer, stop reason, elapsed time, token usage and cost estimate when returned. Resume existing receipts; honor `STOP` without killing active calls.
- Use a supported runtime bootstrap (Node 24.18.1, Pi 1.6.0); record the runtime and runner/protocol hashes in the Qwen metadata. Receipts will use the distinct route directory `receipts/qwen38/` and model key `qwen38`.

## Scoring and reporting

Primary success is passing an original or flagging a mutant as `INCORRECT` at **exactly its planted step**. Retain all cases in the headline results. Report originals passed, false alarms, exact planted-step identifications, mutants incorrectly passed, wrong/missing-step flags, unparsed output, token truncation, timeouts, retry/transport accounting, request timing, token totals and cost estimates. Distinguish terminal-response counts from all-86-scheduled sensitivity counts. Parseable-response and matched-available-response views are secondary and post-hoc, not superiority tests over independent cases.

Produce separate `QWEN38-COMPARISON.md`, `qwen38-comparison.json`, per-case CSV data and `QWEN38-RECEIPTS.sha256` after all 86 samples reach terminal or exhausted-transport outcomes. Compare against immutable Kimi, free DeepSeek and larger GLM primary results. MiniMax remains a separately launched route and is not required to finish before the Qwen report; if included, its status must be stated explicitly. The completed three-model reports and `PRIMARY-RECEIPTS.sha256` must not be replaced by any Qwen/MiniMax result.

Different launch times/provider load are uncontrolled. No expected-label comparison is independently human-certified mathematical correctness. Two samples of a case are not two independent cases. Repeated transport failure or output truncation is not itself a demonstrated mathematical error. Model agreement and successful Lean compilation do not establish a continuum Yang–Mills mass-gap proof; actual research flags still require independent human/strategist or Lean adjudication.

All evidence remains under `<private-run>/full-replay/`. No credential or private proof source is copied into public `mathphysics`.
