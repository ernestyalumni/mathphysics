# One-time Node-24 recovery of the conditional 32768-client-budget diagnostic

Prepared after inspecting the failed fallback receipts, before any mathematical request in this recovery. This is a one-time bootstrap-repair addendum, not a fresh model/configuration search.

## Why the original fallback did not test the model

All 12 fallback attempts exited 1 in about 0.1 s, with no assistant message, no terminal model response, no usage and no estimated cost. Every receipt contains the same startup exception: `node:module` does not export `enableCompileCache`. System `/usr/bin/node` is v12.22.9; the known-working Node is `~/.nvm/versions/node/v24.18.1/bin/node` (v24.18.1), with Pi 1.6.0. No provider cap acceptance or mathematical performance can be inferred from these startup failures.

The four original fallback cases × three exhausted transport attempts and the original selector's blocked decision are preserved unchanged. This addendum is necessary because a successful transport retry budget was exhausted without ever starting Pi.

## Frozen recovery

- Re-run all four already-selected cases, exactly once each unless a transport failure triggers the existing at-most-three transport attempts. Do not select cases based on returned mathematics, and do not retry returned wrong/malformed/truncated responses.
- Exactly the original `low_verdict_first` prompt/instruction, `z-ai/glm-5.3`, concurrency two, 240-second timeout, provider route, parser, expected-label comparison logic and isolated `PI_CODING_AGENT_DIR` whose displayed client `maxTokens` is 32768.
- Only bootstrap difference: prefix Node 24 and the Pi executable directory to this recovery process's `PATH`, after verifying `node --version` == v24.18.1 and Pi `--version` == 1.6.0. No mathematical request may occur if this preflight fails.
- Import and execute the frozen original fallback functions, whose hash must match `configuration-hashes.json`. Do not edit that wrapper, the original arms, their receipts, the model listing or its credential reference. Redirect all new writes to this recovery directory.
- Preserve new receipts, raw stderr, stop reasons, usage/cost and an explicit `runtime_node_binary` command-option annotation. No claim that the provider accepts/enforces 32768 output tokens until returned request evidence supports it.
- The primary 43-case reports and primary receipt manifest remain frozen. This remains a four-case engineering diagnostic, not an independent mathematical benchmark.

## Gate and stop rule

Write a separate `effective-selection.json`. If this recovered arm has at least 3/4 parseable terminal answers, expose it under the existing selector-compatible schema (`low_verdict_first`, `explicit_requested_32768`). Include the immutable original-selection hash, recovery metadata and this preregistration hash. If it does not pass, select nothing and block prospective validation.

This one-time recovery does not authorize further caps, timeouts, system prompts, models, routes, Thinking settings, or configuration tuning. Those require another explicit preregistration. It authorizes only this repaired invocation environment and, on passage of the unchanged ≥3/4 parseable gate, the already frozen 21-case prospective validation under the resulting same model configuration and the same Node-24 bootstrap. If that gate fails, the 21-case set remains uncalled.
