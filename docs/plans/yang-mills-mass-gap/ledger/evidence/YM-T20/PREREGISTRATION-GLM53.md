# GLM-5.3 replay addendum — fixed before its first full-run call

Date: 2026-11-25 Pacific. Ernest has now requested the proof-checking comparison using **Kimi, free DeepSeek, and the larger GLM model `z-ai/glm-5.3` instead of `z-ai/glm-5.3-flashx`**.

## Primary comparison

- `moonshotai/kimi-k3.1-pro`
- `deepseek/deepseek-v4.1-flash-1017-free`
- `z-ai/glm-5.3`

The headline requested comparison will use these three routes. The already active FlashX run is retained as supplementary evidence rather than interrupted mid-request; its results are not substituted for GLM-5.3. This avoids losing in-flight evidence and preserves all existing receipts.

## GLM-5.3 protocol

Except for the model identity and the fact that this route runs in a separate process, every rule is inherited from `PREREGISTRATION.md`:

- Frozen `../frozen-test.json`, 43 historical held-out cases, SHA-256 `b90342fb05b0467c673bd7578018ad97c039d8fdfd668082b26cb6e3f06193fe`.
- Two completely fresh samples per case: 86 scheduled GLM-5.3 samples. No FlashX receipt or selective pilot response is reused.
- Same system prompt, numbered case prompts, provider `tokenrouter`, text-only fresh Pi calls, `--thinking low`, no tools/MCP/extensions/skills/context/session, and unchanged archived parser.
- Same 240-second request timeout; at most three requests per scheduled sample; retries only for timeout/no-terminal-response/assistant errors. Mathematical or format failures are terminal and are never retried.
- Concurrency two. Separate `receipts/glm53/` directory, immutable per-attempt receipts, and separate run-log/metadata/status files.
- Same success and error taxonomy: original passed/false alarm; exact planted step/mutant passed/wrong or missing step; unparsed and token truncation; transport failures, attempts, timeouts, token usage, Pi-estimated cost, and two-sample agreement.

The three-model report will not compare unreturned GLM responses as mathematical mistakes, and will show terminal-only and all-scheduled sensitivity denominators. No result-dependent exclusions or model substitution is allowed. The existing known limitation about the historical `ym-s04` case and any other label/case concerns remain secondary explicitly labelled sensitivity analyses.

The GLM run starts later than Kimi/DeepSeek/FlashX, at a different provider load. The same request parameters are used, but no claim is made that simultaneous execution or identical provider-side output caps were controlled.
