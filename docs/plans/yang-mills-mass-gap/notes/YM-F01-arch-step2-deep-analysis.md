# YM-F01, architecture step 2 of 13: deep analysis

Intent: `lean-formalisation-of-load-bearing-steps-and-a-l-3be8135f` (`add`, `existing_pdd_change`, technology `docker`).
Author: Claude Opus 5, 2026-11-20. Input: `docs/PRODUCT_INTENT.md` + step-1 report.
Not posted to GitHub: the issue reference is `local-intent:...` and the repository argument is `/`; the
program's protocol forbids posting outside the local network without Ernest's approval. Report kept here.

**Concurrency note.** `docs/plans/yang-mills-mass-gap/formal/` is being written by another agent while this
analysis was produced (`Casimir.lean` and `check.sh` appeared mid-run). Nothing under `formal/` was modified here.

## Status

**Analysis Complete.** The stack is determined, not ambiguous. Two decisions must be taken before step 3
delivers anything runnable: the Lean image pin (D1) and whether prover weights are fetched (D2).

## Feature decomposition

| Feature (step 1) | Functional units | Candidate module |
|---|---|---|
| Pick 3–5 load-bearing steps | parse `yang-mills-lattice.tex` for `\label{ym:...}` inside theorem environments; extract statement body; resolve `\ref`/`\eqref` edges into a dependency graph; apply and record a stated load-bearing criterion; normalise and hash each statement | `step_corpus` |
| Lean project for our own statements | module list and dependency order; Mathlib supplied by the image (no `lake`, no `lakefile`, no `lean-toolchain` present or needed); skeleton with hypotheses explicit and `sorry` in the proof slot | `lean_project` (configuration + hand-authored sources; **not** prompt-generated) |
| Docker-only Lean execution | image selection; bind-mount `formal/` only; `--network none`; non-root `-u uid:gid`; per-module compile in dependency order; `LEAN_PATH` to the olean build dir; per-module timeout; structured result | `lean_compile_runner` (embryo: `formal/check.sh`) |
| Statement-fidelity audit | inventory declarations from the `.lean` sources; pair label ↔ declaration; record which parts enter as hypotheses rather than being proved; auditor + date; emit a LaTeX map section and flag orphans in both directions | `fidelity_map` |
| Local prover-model harness | endpoint guard (loopback/private range only); prompt assembly from statement + its stated hypotheses; N sampled attempts under a wall-clock budget; splice candidate proof into a scratch copy of the module; compile each in the offline container; append JSONL; stamp the unreviewed label | `prover_harness` |
| Lean-only acceptance judge | pure predicate over (expected statement hash, compile log, `#print axioms` output) | `lean_judge` |
| Evidence ledger and handoff | receipt with command + sha256 under `ledger/evidence/YM-F01/`; rows in `ledger/claims.md`; pass **and** fail logged; `status` / `report` markdown | `evidence_receipts` |
| Scope guard | no mass-gap claim; label vocabulary `reproduced` / `source-claims` / `unreviewed` enforced on every emitted artifact | review gate, enforced in `evidence_receipts` (not its own module) |
| (shared) paths, images, budgets | env-var resolution with no hardcoded paths | `formal_config` |

## Shared concerns

- **Auth:** none, and that is a real simplification — no API keys, no tokens, nothing to leak. The analogue
  that does bite is **network containment**, and it is a hard requirement, not a preference: every Lean
  *check* runs `--network none`; the model endpoint must be validated as loopback or RFC1918 before a request
  is made (the `ipaddress` check in `bin/ym-micro` is the precedent to reuse verbatim); network is permitted
  only for an image pull. Needed by `lean_compile_runner` and `prover_harness`.
- **Error handling:** a failed compile is *data*, not an exception. Each stage returns a structured verdict;
  only the top-level entry point maps "any FAIL" to a nonzero exit (`check.sh` already does this with
  `FAIL <module>` / `OK <module>` lines). Timeouts are mandatory on both the compile and the generation call,
  or the prover loop runs unbounded. Model output is never executed outside the bounded Lean compile.
- **Validation:** sha256 equality of the normalised statement text (the single most important check — it is
  what stops a model from "proving" a weakened goal); axiom set ⊆ `{propext, Classical.choice, Quot.sound}`;
  rejection of `sorry` / `sorryAx`, of `native_decide`, and of any new `axiom` or `@[implemented_by]`;
  JSONL schema validation on queue and result records; writes confined to the state directory;
  label-vocabulary validation on emitted claims.
- **Logging / observability:** append-only JSONL in a state directory outside git; regenerable markdown
  `report` / `status`; a receipt per run with the exact command and output checksums. Failures are logged
  under the same discipline as passes (06-AGENT-PROTOCOL.md).
- **Configuration:** environment variables only; no machine-local path in any tracked file (AGENTS.md).

## Tech stack (confirmed)

- **Lean 4 + Mathlib, in Docker, never on the host** — explicit in the intent. Two images exist locally:
  `ym-lean:4.32.0` (3.39 GB, built from `docs/plans/yang-mills-mass-gap/docker/lean`, and `bin/ym-lean` mounts
  only the audited external clone pinned at `3456e5e`), and `claw-work/lean-agent:lean4.35.0-mathlib` (9.87 GB,
  Mathlib prebuilt), which the in-flight `formal/check.sh` already uses via `YM_FORMAL_IMAGE`.
  **Recommendation (D1):** keep `ym-lean` for the external proof and use the prebuilt lean-agent image for
  `YMFormal`; it avoids a multi-hour cold Mathlib build and a second multi-GB working copy, and the Mathlib
  revision need not match the audited copy because we are stating our own propositions, not replaying theirs.
  Pin the image **digest** in the receipt, since the tag is mutable.
- **No `lake`.** There is no `lakefile`/`lean-toolchain` anywhere in the repo, and the chosen approach does not
  need one: `lean -R /work -o build/<Mod>.olean <Mod>.lean` with `LEAN_PATH` extended by the build dir. This
  corrects the step-1 expectation of a lake package; it also removes `lake exe cache get` and its network need.
- **Bash** for the container wrappers (precedents: `bin/ym-lean`, `bin/ym-local-llm`, `formal/check.sh`).
- **Python 3, standard library only** for corpus extraction, the judge, the harness and the ledger
  (precedents: `bin/ym-micro`, `bin/ym_micro_plan.py`, `watch/ym_watch.py`; tests as `tests/test_*.py`).
- **llama.cpp OpenAI-compatible server** on the RTX 3060 via `bin/ym-local-llm start <model> <alias> <ctx>`,
  serving loopback `http://127.0.0.1:${YM_LLM_PORT:-8081}/v1`. `nvidia-container-toolkit` is still purged on
  this host, so the hand-mounted-driver workaround remains load-bearing.
- **Storage:** flat files. No database. JSONL (outside git), markdown notes, `ledger/` receipts, LaTeX map.
- **Rust is not displaced, and must be cited rather than duplicated.** Two of the four statements already in
  flight have exact-arithmetic owners in the repo: `ym:prop:purity` → `Rust/lattice-gauge/src/spectral.rs`
  (+ `tests/purity.rs`, prompt `prompts/spectral_rust.prompt`), and `ym:prop:casimir-sun` →
  `Rust/lie-groups/tests/special_unitary_casimir.rs`. That pairing *is* the 2026-11-15 amendment's
  "Rust discovers and falsifies, Lean certifies", and the reusable-library policy forbids re-implementing the
  arithmetic inside the Lean tier or inside new Python.
- **Frontend:** none. `TypeScript/physics-workbench` is a different stream and stays out.

## Module candidates

1. **`formal_config`** — resolve images, the `formal/` directory, the state directory, attempt budget and
   per-compile timeout from environment variables with documented defaults. *Foundational.*
2. **`step_corpus`** — from `documents/mono-notes/topics/yang-mills-lattice.tex` produce one record per
   labelled statement: label, environment kind, statement text, referenced labels, normalised sha256, and
   the recorded reason it is load-bearing. *Foundational* (reads the notes only).
3. **`lean_judge`** — pure verdict function over (expected statement hash, compile log, axiom report) →
   `PASS | FAIL(reason)`. Deliberately has no I/O so it can be driven entirely by fixtures. *Foundational.*
4. **`lean_compile_runner`** — compile named modules in dependency order inside the offline container and
   return a structured `CompileResult` per module. Depends on `formal_config`. *Dependent.*
5. **`fidelity_map`** — pair labels with our declarations, carry the auditor and date, emit the LaTeX map
   section and flag unmapped labels and unmapped declarations. Depends on `step_corpus`. *Dependent.*
6. **`prover_harness`** — queue, endpoint guard, prompt assembly, sampling, splice, compile, record, label.
   Depends on `formal_config`, `step_corpus`, `lean_compile_runner`, `lean_judge`. *Dependent.*
7. **`evidence_receipts`** — receipts, claims-ledger rows, `report` / `status` markdown, pass-and-fail logs,
   label-vocabulary enforcement. Depends on `lean_judge` and `prover_harness` output. *Dependent.*

`lean_project` (the `.lean` sources themselves) is listed in the decomposition but is **not** a prompt-owned
module: each statement must be read against the LaTeX by hand, which is exactly what a generated artifact
cannot certify. The sources are hand-authored; only the tooling around them is prompt-owned.

## Inter-module interfaces

- `step_corpus` → `fidelity_map`, `prover_harness`: `Step{label, kind, statement_text, refs[], sha256, why}`.
- `formal_config` → `lean_compile_runner`, `prover_harness`: resolved paths, image reference, timeouts, budget.
- `lean_compile_runner` → `lean_judge`: `CompileResult{module, exit_code, stdout, warnings[], axioms[], seconds}`.
- `lean_judge` → `prover_harness`, `evidence_receipts`: `Verdict{module, decl, pass: bool, reasons[], hash_ok,
  axioms_ok, sorry_free, statement_unchanged}`.
- `prover_harness` → `evidence_receipts`: `ProofAttempt{step_label, decl, model, alias, prompt_sha256, sample_index,
  candidate, verdict, seconds, label: "unreviewed (local model output; not evidence)"}`.
- `fidelity_map` → `documents/mono-notes/topics/yang-mills-source-map.tex`: a **new, separate** section.
  The existing `ym:sec:source-map` maps our labels to *OpenAI's* declarations; mixing our own statements into it
  would blur "the group proved this in Lean" with "we stated this in Lean" — the one confusion this packet exists
  to prevent.
- `Rust/lattice-gauge::spectral`, `Rust/lie-groups` tests ↔ `YMFormal.Purity`, `YMFormal.Casimir`: both anchored
  to the same `\label`, exact-arithmetic falsification on one side and certification on the other; cross-reference
  by label, no shared code.
- `bin/ym-local-llm` → `prover_harness`: the loopback endpoint. The harness never starts or reconfigures a server
  (an `ym-micro` invariant worth keeping).

## Cross-cutting concerns (from the PRD)

| Concern | Modules that need it |
|---|---|
| Auth / secrets | none — no credentials anywhere in this tier |
| Network containment (the real "auth") | `lean_compile_runner` (`--network none`), `prover_harness` (endpoint allowlist) |
| Error handling as structured verdicts + mandatory timeouts | `lean_compile_runner`, `prover_harness`, `lean_judge` |
| Validation (hash, axioms, `sorry`, schema, path confinement) | `lean_judge` (core), `step_corpus`, `prover_harness`, `evidence_receipts` |
| Logging / receipts / reproducibility | `evidence_receipts` (owner), every other module as producer |
| Configuration via env only | `formal_config` (owner), all consumers |
| Untrusted-output handling and the honesty vocabulary | `prover_harness`, `evidence_receipts`, `fidelity_map` |

## Identified shared identifiers

- **Database collections/tables:** none. Flat files: `formal/YMFormal/*.lean`, `formal/check.sh`,
  `ledger/evidence/YM-F01/receipt.txt`, rows in `ledger/claims.md`, `notes/YM-F01-*.md`, and
  state-directory JSONL (`attempts.jsonl`, `verdicts.jsonl`) outside git.
- **API endpoints:** none built. One consumed, local only: `POST http://127.0.0.1:${YM_LLM_PORT:-8081}/v1/chat/completions`
  and `GET /v1/models`.
- **Shared types:** `Step`, `StatementHash`, `LeanStatement`, `CompileResult`, `Verdict`, `AxiomSet`,
  `ProofAttempt`, `Receipt`, `SourceMapEntry`.
- **Environment variables.** Existing: `YM_REPOS`, `YM_LEAN_WORK`, `YM_LEAN_IMAGE`, `YM_LEAN_TTY`,
  `YM_FORMAL_IMAGE` (already honoured by `formal/check.sh`), `YM_LLM_PORT`, `YM_CLAW_LLAMA`,
  `YM_LLM_REASONING_BUDGET`, `YM_NUMCHECK_IMAGE`, `CUDA_VISIBLE_DEVICES`. Implied new:
  `YM_FORMAL_STATE`, `YM_FORMAL_TIMEOUT`, `YM_FORMAL_ATTEMPTS`, `YM_FORMAL_MODEL`, `YM_FORMAL_ALIAS`.

## Decisions required before step 3 produces anything runnable

- **D1 — image and layout.** Adopt `claw-work/lean-agent:lean4.35.0-mathlib` (digest-pinned) and the
  no-`lake` `lean -R` layout, as `formal/check.sh` already does; leave `ym-lean:4.32.0` to the audited clone.
- **D2 — prover weights.** No prover-specific GGUF is on this machine; `notes/local-models.md` inventories
  general 9B chat distills only, and YM-T17 calibration measured them weak on multi-step bookkeeping. Either
  run a general 9B as the honest baseline (zero download, low expected yield) or fetch a prover model
  (DeepSeek-Prover / Goedel / Kimina). A weight download is an inbound fetch rather than a publication, but it
  is still off-network traffic — put it to Ernest rather than assume it.
- **D3 — interaction mode.** One-shot: sample a proof term or `by` block, splice into the `sorry` slot, compile.
  A tactic-feedback loop needs a Lean REPL service inside the container and should be its own packet.
- **D4 — map separation.** A new map section for our declarations, never appended to `ym:sec:source-map`.
- **D5 — tracker.** `BOARD.md` has no `YM-F*` packet; YM-F01 needs an entry, a tier and a noted dependency on
  YM-O01 (still in review).
- **D6 — success criterion, fixed before the run.** State now what counts as informative (e.g. at least one of
  the non-trivial statements closed with a clean axiom set) and record the zero-closure outcome as a result
  rather than as a packet failure.

## Judge fixtures the next step must ship

A statement with a correct proof (must PASS); the same with `sorry` (must FAIL); one with an added `axiom`
(must FAIL); one where the *statement* was edited (must FAIL on the hash, and this is the fixture that matters
most); one that exceeds the timeout (must FAIL); a `native_decide` proof (must FAIL); and a no-network
assertion inside the container. Plus a trivial Mathlib-free smoke module so the harness tests do not need
the 9.87 GB image.

---
*Proceeding to step 3: research.*
