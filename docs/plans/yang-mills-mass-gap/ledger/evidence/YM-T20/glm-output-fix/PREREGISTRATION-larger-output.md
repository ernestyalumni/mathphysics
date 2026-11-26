# Conditional separate GLM output-budget diagnostic — frozen before calls

Created 2026-11-26 while the two provider-default diagnostic arms are still running. This file does not authorize calls until both original arms are resolved and neither reaches at least 3/4 parseable answers, as required by `../PREREGISTRATION.md`.

## Fixed configuration

- Model/provider: `z-ai/glm-5.3`, `tokenrouter`.
- Client model metadata in `pi-agent/models.json`: `maxTokens: 32768`, contextWindow 128K. Every other relevant request setting remains the same. This is an explicit requested client override, not a claim that TokenRouter permits 32,768 output tokens; provider enforcement is determined only from actual request receipts.
- Thinking: `low`.
- System instruction: exactly the frozen `low_verdict_first` text in the old diagnostic.
- Four fixed case IDs and unchanged prompts selected in the original diagnostic.
- `--thinking low`, 240-second request timeout, concurrency two, at most three attempts only for transport failure. Never selectively retry a returned wrong, malformed or token-truncated response.
- Fresh text-only Pi processes, `--no-session --no-tools --no-mcp --no-extensions --no-skills --no-context-files`.
- Separate private PI agent directory containing no copied credentials. It uses a new model-metadata JSON pointing to the same provider; its apiKey entry is the documented reference to the existing configured key command, not a copied credential.
- Same archived `ym-micro.parse_verdict`; immutable standalone receipts.

## Interpretation and selection

This is an engineering fallback, not a mathematical-performance retest. It is run only if neither original arm reaches ≥3/4 parseable answers. Its outcome does not overwrite the original diagnostic or primary report.

If it reaches ≥3/4 parseable answers, it is eligible for prospective evaluation on the frozen 21-case adversarial candidate set. Otherwise no GLM configuration has passed the engineering gate. Selection does not certify mathematical correctness, and new mechanical labels remain candidates until independent human adjudication. If the tests reveal that a further setting change is needed, that must be a separately preregistered addendum rather than silently extending this arm.
