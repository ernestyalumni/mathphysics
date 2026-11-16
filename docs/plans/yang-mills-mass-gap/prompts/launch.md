# Launch recipes

How to start a worker in each harness Ernest has. `PLAN` is this directory's parent (`docs/plans/yang-mills-mass-gap`). The prompt text is `prompts/worker.md` (or `reviewer.md`, `coordinator.md`) with the bracketed fields filled in; save the filled prompt to a file outside git, for example `/tmp/ym-prompt.txt`.

Every recipe starts in the canonical checkout (or a LAN clone of it, protocol §10):

    cd <MP>            # the mathphysics checkout
    P=docs/plans/yang-mills-mass-gap

| Harness | Tier | Command | Notes |
|---|---|---|---|
| Claude Code | T2–T4 | `claude -p "$(cat /tmp/ym-prompt.txt)"` or interactive `claude` | Use interactive for T3 mathematics. |
| Codex | T2–T3 | `codex exec "$(cat /tmp/ym-prompt.txt)"` | A different model family from Claude: preferred for reviews of Claude's work, and the reverse. |
| Grok Build | T2–T3 | `grok "$(cat /tmp/ym-prompt.txt)"` | Third family for red-team passes. |
| Pi | T1–T2 | `pi -p --model <provider/model> "$(cat /tmp/ym-prompt.txt)"` | With a local model: merge `claw-work/Deployments/Scripts/llama-cpp-server/examples/pi/` into Pi's configuration first, then `--model <local alias>`. |
| Hermes Agent | T1–T2 | `hermes -z "$(cat /tmp/ym-prompt.txt)"` | For a local server use `--provider custom:<name>`; see `TMP-HERMES-AGENT.md` at the repository root. |
| Local Qwen3.5, no tools | T0 | `$P/bin/ym-micro …` | Not written yet: packet YM-T14. Until then use Pi or Hermes with the local alias for T1 packets only. |

Check a command's flags with `<harness> --help` before scripting it; the table was written from the help output on 2026-11-15.

## Starting the local model server

From `claw-work/Deployments/Scripts/llama-cpp-server` (see its README):

    ./run-recommended.sh --dry-run
    ./run-recommended.sh                       # Qwythos 9B Q5: http://127.0.0.1:8080/v1, alias qwythos-9b-q5-k-m, context 65536, one slot

**As of 2026-11-15 this fails on the desktop**: the NVIDIA container toolkit was purged on 2026-11-11 (`notes/local-models.md` §3a). Ernest has to reinstall it.
    ./launch.sh --status
    ./launch.sh --stop

The server uses the RTX 3060. Do not start it while a GPU Monte Carlo run is using that card, and the reverse.

## What to give which tier (while Ernest is away)

| Tier | Packets that fit today |
|---|---|
| T0/T1 local 9B (use `qwythos-9b-q5`) | See `notes/local-models.md` §6, which is based on a measured probe: flagging referee, verbatim indexing, Lean signature cards, Rust tests with the expected **value** given, reruns. Not: digests, derivations, calculations, any review that closes a task. |
| T2 | YM-T04 (`lean-graph`), YM-T14 (micro-task runner), YM-N01 (SU(2) Monte Carlo baseline, on top of `Rust/lattice-gauge`; read its README first), YM-R06–R08, YM-W01–W02, reviews of YM-T01, T02, T05, O01. |
| T3 | YM-R01–R05, YM-A03 (twice, two families), YM-A04, YM-G02 (review the proofs in the mono-notes topic), YM-T11. |

A local 9B model given a T2 or T3 packet will produce fluent text that is wrong. The protocol's review step is what stops that from entering the ledger; do not skip it to save time.
