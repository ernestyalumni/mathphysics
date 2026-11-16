# YM-T14 — Micro-task runner for small local models (`bin/ym-micro`)

Tier T2 · Workstream B (tooling for the program itself) · Depends on: nothing · Blocks: use of tool-less 9B models · Review: T2, different family · Budget: half a day

## Why
Ernest's local models (Qwen3.5 9B, 32,768-token context, one server slot, about 40 tokens per second) cannot navigate a repository or hold a paper. They can transform one short snippet into one short structured answer. A runner that prepares the snippet, calls the local OpenAI-compatible endpoint, and checks the answer mechanically turns them into useful parallel labour without an agent harness.

## Inputs
`06-AGENT-PROTOCOL.md` §12; `ARM/data/key-statements.tsv` (990 rows: id, unit, Lean name, link with path and line, importance, one-line summary); `LEAN` sources; `claw-work/Deployments/Scripts/llama-cpp-server/README.md` for the endpoint.

## Do
Write `PLAN/bin/ym-micro` in Python 3 with the standard library only (`urllib`, `json`, `argparse`). No model-written code is ever executed.

1. `ym-micro prompt lean-reading <K-id>`: print a self-contained prompt of at most 1,500 tokens: the Lean text from three lines above the declaration's line to the end of its signature (stop at `:= by`, `:=`, `where`, or 30 lines), and a fixed answer template with three headed fields: `KIND` (definition / structure / theorem), `STATEMENT` (one or two sentences using only identifiers visible in the snippet), `IDENTIFIERS` (comma-separated identifiers copied from the snippet). Do **not** include Armstrong's summary in the prompt.
2. `ym-micro run lean-reading --count N [--endpoint URL] [--model NAME]`: for the next N rows without an output (main results first), call `POST {endpoint}/chat/completions`, temperature 0.2, `max_tokens` 400. Default endpoint `http://127.0.0.1:8080/v1`. Refuse any endpoint whose host is not loopback or a private-range address.
3. Mechanical acceptance of each answer: all three headers present; every item of `IDENTIFIERS` occurs verbatim in the snippet; `STATEMENT` under 80 words. Otherwise save under `rejected/`.
4. Save accepted answers to `PLAN/notes/lean-readings/<K-id>.md` with a header: model alias, UTC time, `LEAN` path and line, SHA-256 of the snippet, and the line `label: unreviewed (T0 model output; not evidence)`.
5. `ym-micro compare <K-id>`: print the model's `STATEMENT` next to Armstrong's summary for a human or T2 reviewer. No automatic verdict.
6. `ym-micro status`: counts of done, rejected, remaining.

## Added 2026-11-16, after the probe in `notes/local-models.md`
Two more task types, both justified by measurement there, and both more useful than `lean-reading`. Build them first.

7. `ym-micro run referee --file <text file of numbered steps> --seeds K`: send the steps with the fixed referee question of `ledger/evidence/local-model-probe/probe.py` (`REVIEW_QUESTION`) K times with different seeds. Output: a table of verdicts. Report **only** the INCORRECT flags, with step number and reason, to `PLAN/notes/referee-flags/<name>.md`. Never report "passed". Inputs to prepare by script: one file per proof of the mono-notes topic (split at `\begin{proof}`), with each displayed equation as one step.
8. `ym-micro run index --file <excerpt>`: ask for every numbered equation, definition and theorem statement, verbatim, with its number. Accept a line only if, after whitespace normalization, it occurs in the excerpt. Output: `PLAN/notes/indexes/<source>.tsv`.

## Deliverable
The script, a `--dry-run` mode that needs no server, unit tests under `PLAN/bin/tests/` (snippet extraction on five fixed rows; the acceptance check on three good and three bad answers), and a note `notes/YM-T14-micro-runner.md` reporting a real run of 20 rows: acceptance rate, tokens per second, three examples compared with Armstrong's summaries, and an honest judgment of whether the output is worth a reviewer's time.

## Acceptance
Tests pass. The 20-row run is reported with its numbers. If fewer than half of the accepted readings are judged usable by the reviewer, say so and recommend stopping; that is a valid result.

## Do not
Send any prompt to a non-local endpoint. Start or reconfigure the model server from inside the script. Treat the readings as evidence anywhere.
