# Statement checker, architecture step 1 of 13: PRD analysis

Intent: `automatic-lean-statement-checker-in-rust-3c8e1031` (`add`, `existing_pdd_change`, technology `Rust`).
Author: Claude Opus 5, 2026-11-20. Input: `docs/PRODUCT_INTENT.md` (the entry added by
`pdd intent apply`), `request-2026-11-21-statement-checker.txt`, `intent-plan-2026-11-21-statement-checker.json`.

**Not posted to GitHub.** The issue reference is `local-intent:automatic-lean-statement-checker-in-rust-3c8e1031`
and the repository argument is `/`; the workflow's own final override forbids `gh`, and the program protocol
forbids posting outside the local network without Ernest's approval. The report is kept here instead.

**No packet ID is allocated.** `BOARD.md` has no entry for this unit; the next free tooling number is YM-T18.
Allocating it is a step for the plan, not for this analysis.

## Status

**Analysis Complete.** The stack is fixed by the intent and by the existing crates, not ambiguous. Two premises
in the request do not survive contact with the repository and are recorded under gaps (G1, G2); neither blocks
architecture, but both change what the first unit should contain.

## Project summary

- **Title:** Automatic Lean statement checker in Rust (with a Rust/CUDA replacement for the SymPy/NumPy
  scripts the program currently leans on).
- **Goal:** decide mechanically whether a Lean statement says what the source claim says, so that a compiling
  statement is no longer taken as a faithful one, and so that the fidelity audit stops depending on Ernest
  reading every statement by hand.
- **Why now:** `notes/YM-F01-lean-formalisation.md`, experiment 2. Of 28 formalizer attempts, 9 compiled and
  **4 of those 9 were false statements**. Compiling is no evidence of faithfulness, and hand-grading is the
  bottleneck.

## Features identified

1. **Paper claim to an executable predicate.** A source claim (one `ym:*` label in
   `documents/mono-notes/topics/yang-mills-lattice.tex`) becomes a predicate over an instance, plus a generator
   of instances satisfying its hypotheses. Precedent, in Python/NumPy, is `bin/ym-numcheck`, where a local model
   writes `check(rng)` and the harness runs it in a sandboxed container.
2. **Lean statement to an executable predicate, mechanically.** Not from the Lean source text: from Lean's
   **elaborated term**, exported by a Lean metaprogram inside the Docker image and consumed in Rust. This is
   what makes the `(P₁ : ℝ := z₂ / z₁^2)` bug visible — an optional parameter is a universally quantified
   variable with a default, and three such `defect_telescope` statements compiled and were false.
3. **Differential comparison on random and edge-case instances.** Both predicates are evaluated on shared
   instances from a deterministic seeded stream; a disagreement is reported with the witnessing instance.
4. **Vacuity check.** Search for an instance satisfying the Lean statement's hypotheses. No witness found means
   the statement may be empty — the failure mode recorded as "a statement made vacuous by contradictory
   hypotheses counts as wrong".
5. **Counterexample search.** Edge-case and adversarial instances against the Lean conclusion alone, so a false
   statement is refuted rather than merely flagged as divergent.
6. **Hypothesis-necessity check.** Drop each hypothesis in turn and search for a counterexample to the weakened
   statement: a witness means the hypothesis is load-bearing, no witness means it is a candidate for removal.
   This is the check that would have classified YM-F01's own additions (`λ₂ > 0`, positive defects) and its
   omissions (`𝖯(n) > 1/2`, `Δ < log 2`).
7. **Three-valued verdict with an explicit supported fragment.** `AGREE` / `DISAGREE` + witness /
   `UNDECIDABLE-in-fragment`. Silence outside the fragment is the one outcome that must be impossible, because
   that is the current failure mode in a new costume.
8. **Pre-registered validation on two corpora.** (a) The graded formalizer attempts in
   `~/Data/ym-formal/formalizer-run1/grades.json`: 28 graded, of which 9 elaborate — 5 equivalent,
   4 wrong. (b) Mutants of the 18 reference theorems frozen at commit `1cdf83a`
   (`formal/YMFormal/{Purity,Defect,ConditionalGap,Casimir}.lean`), as false-alarm and detection controls.
9. **Rust numeric core to replace the SymPy/NumPy/SciPy use in the loop.** Deterministic counter-based RNG,
   certified interval arithmetic for `exp`/`log` (roadmap T12, currently deferred), closed-form `tsum` families,
   small exact dense matrices. 70 files in the repository import `sympy`, `numpy` or `scipy`; inside this
   program the live ones are `ledger/evidence/check_oneloop.py` (SymPy),
   `ledger/evidence/transfer-positivity/check_su2_character_coefficients.py` (SciPy `quad`, `iv`) and
   `bin/ym-numcheck` (NumPy, in the sandbox image).
10. **Evidence ledger and receipts.** Append-only JSONL in a state directory outside git, a receipt per run
    with the command, the image **digest** and output checksums, rows in `ledger/claims.md`, passes and
    failures logged alike, and the label vocabulary enforced on every artifact.

## Tech stack

- **Backend:** Rust, first (explicit in the intent and in the roadmap's durable requirements). Existing crates
  are `Rust/{qft-symbolics,lie-groups,lattice-gauge,physics-contracts,physics-gateway,physics-tools}`;
  `Rust/lean-graph` is planned for T04 and is the nearest existing home for the Lean-term side.
  Dependencies are pinned and minimal (`num-bigint`, `num-rational`, `num-traits`, `serde`, `serde_json`);
  builds are `--offline --locked`; no RNG, array or interval crate is present today.
- **Formal side:** Lean 4 + Mathlib, **in Docker, never on the host**. Images: `ym-lean:4.32.0` (built from
  `docker/lean`, mounts the audited external clone pinned at `3456e5e`) and
  `claw-work/lean-agent:lean4.35.0-mathlib`, which `formal/check.sh` uses via `YM_FORMAL_IMAGE`.
- **Frontend:** none. CLI subcommands plus regenerable markdown reports, following `prover.py`
  (`tasks`/`generate`/`judge`/`report`) and `ym-numcheck` (`run`/`score`/`resandbox`).
  `TypeScript/physics-workbench` is out of scope for this unit.
- **Database:** none. JSONL state outside git, committed fixtures and receipts.
- **Other:** Docker with `--network none` for every check; a llama.cpp server on loopback for any
  model-written predicate, behind the loopback/RFC1918 endpoint guard used by `bin/ym-micro`; GGUF local
  models (Goedel-Formalizer-V2-8B, Goedel-Prover-V2-8B, Phi-4-reasoning-plus). CUDA C++ only where a
  data-parallel workload is **measured** to be the bottleneck; no such workload is identified in this unit.

## Non-functional requirements

- **Offline and contained.** Every Lean invocation `--network none`, non-root `-u uid:gid`, `formal/` bind-mount
  only; network permitted for an image pull. The model endpoint is validated as loopback or RFC1918 before a
  request. Nothing is posted outside the local network without Ernest's approval.
- **Determinism and replay.** Seeded instance streams, recorded seeds, statement normalisation and sha256,
  image digest in the receipt (the tag is mutable). A run is re-derivable from its receipt.
- **Soundness direction.** A `DISAGREE` with a witness is a hard result, reproducible from the recorded
  instance. An `AGREE` is **not** a proof: random instances can miss, and quantification over all sequences,
  all matrices and all spectra is not exhaustible. The Lean kernel stays the certifier; this tool is the
  falsifier. Label: unreviewed tooling with measured rates.
- **Calibration as a shipping requirement.** The tool reports its detection rate on mutants and its
  false-alarm rate on known-faithful statements; the YM-T17 precedent (13 of 54 mutants detected by a 9B
  test-writer) is the measurement discipline to repeat.
- **Independence must be structural.** The two translations are only independent if the predicate writer
  never sees the Lean statement and the Lean reader never sees the paper prose, and if the model writing the
  predicate is from a different family than the one that produced the Lean statement. Blinding is an
  architectural constraint, not a convention.
- **Exactness.** Exact arithmetic wherever the mathematics is exact; interval arithmetic with controlled
  rounding for transcendentals; floating point only where marked. `overflow-checks = true` is already the
  profile in every crate.
- **Resource contention.** The RTX 3060 serves both the local model server and any GPU work; check
  `nvidia-smi` before running both.
- **No build products in the repository.** `CARGO_TARGET_DIR` outside the tree; no venvs, weights or run
  outputs committed.

## Development requirements

- **Local testing.** A mock Lean endpoint (the prover harness was tested against a mock server, commit
  `9f89e0d`) and a mock model endpoint, so the harness is testable with neither Docker nor a GPU. Fixtures:
  the 9 elaborating candidates from `formalizer-run1/candidates/*.lean` with their grades, the 18 frozen
  reference statements, and a mutant corpus with **per-mutant ground truth**. Every public function needs a
  positive test against an independent oracle, a negative control that must fail, and a doc comment stating
  conventions.
- **Dev-only UI.** None. The analogue is `report`/`status` markdown regenerated from the JSONL, plus a
  `resandbox` path to re-run saved CRASH/ERROR cases after a harness fix.
- **Utility modules implied.** normalise + hash statements; export (Lean elaborated term to JSON); mutate
  (planted-defect generator); sample (seeded instance streams); grade/score; import (lift `grades.json` into a
  committed fixture); replay (verify a receipt); guard (endpoint and write-path validation).
- **Migration/setup scripts.** No database. Needed: image-digest resolution and pinning; state-directory
  bootstrap outside git; the fixture import above; and a pre-registration file written **before** any
  generation, as YM-F01 did.
- **Environment variables.** `YM_FORMAL_IMAGE`, `YM_NUMCHECK_IMAGE`, `CARGO_TARGET_DIR`,
  `CUDA_VISIBLE_DEVICES`, the model endpoint URL and the state directory. **No secrets and no API keys** —
  there is nothing to leak, which is a real simplification; containment, not authentication, is the control.

## Data entities

`SourceClaim` (`ym:*` label, text, hash) · `Predicate` (expression + instance generator + fragment tag) ·
`LeanDeclaration` (module, name, elaborated term JSON, hypothesis list, binder kinds including defaults,
axioms) · `Instance` (seed, assignment, provenance: random or edge case) · `Verdict`
(AGREE / DISAGREE + witness / UNDECIDABLE + reason) · `VacuityReport` · `NecessityReport` (per hypothesis) ·
`Mutant` (base declaration, mutation kind, ground-truth label) · `GradedAttempt` (the 28) ·
`CalibrationRun` (detection and false-alarm rates) · `Receipt` (command, image digest, checksums) ·
`PreRegistration`.

## Roles

Ernest approves and is the only one who may authorise anything leaving the local network. A strategist agent
plans; worker tiers (Pi, Hermes Agent, Codex, Claude Code, Grok Build, local Qwen3.5) take one packet each; a
reviewer from a second model family accepts or rejects. The Lean kernel is the judge and is not a person.
Local model output is untrusted input, bounded by the compile and the sandbox.

## Gaps and ambiguities

- **G1. "We have a full Rust port already" is not accurate.** `Rust/qft-symbolics` self-describes as
  "exact, deliberately bounded QFT calculations; not a Cadabra-compatible CAS". There is no array layer, no
  RNG, no interval arithmetic (T12 is deferred, "do not start unasked") and no special functions — SciPy's
  `quad` and `iv` are used in the transfer-positivity check and have no Rust counterpart here. **Scope
  decision needed:** port what this checker needs (exact rationals, which exist; a counter-based RNG;
  interval `exp`/`log`; closed-form `tsum` families; small exact matrices) and register the general
  NumPy/SymPy replacement as its own program, rather than inferring it from this intent.
- **G2. "Banish SymPy/NumPy" collides with a durable requirement.** `05-TOOLING-ROADMAP.md` states both
  "Python (SymPy, the existing `.venv`) **stays as an independent oracle**" and "every tool MUST be validated
  against an independent oracle". Removing SymPy removes the oracle for the new Rust numerics.
  **Recommendation for Ernest:** keep SymPy and SciPy as test-only oracles that ship in no discovery path, so
  the loop becomes Rust while independent validation survives. Needs his call; it is his own requirement on
  both sides.
- **G3. The supported fragment is undefined, and it decides the tool's worth.** The statements use `∑'`
  (`tsum`), `Summable`, `Real.exp`, `Real.log`, `Matrix (Fin N) (Fin N) ℂ` and `LinearMap.trace`, with
  `noncomputable` definitions. Universal quantification over `ℕ → ℝ` is not enumerable; `exp` and `log` are
  not rational. A fragment must be pre-registered — plausibly: finitely-supported and geometric sequence
  families with closed-form `tsum`, interval arithmetic for transcendentals, exact rationals for matrices —
  with everything else returning `UNDECIDABLE`, never `AGREE`.
- **G4. The validation set is smaller than the intent assumes.** Only **9 of the 28** graded attempts
  elaborate, so only 9 can reach an elaborated-term comparison: 5 equivalent, 4 wrong. Nineteen are already
  caught by Lean and measure nothing. Nine cases cannot establish a false-alarm rate, so the mutants carry
  the statistical weight and the mutant design is the real experiment.
- **G5. The ground truth is one non-blind grader.** `grades.json` is Claude Opus 5.5's own grading of its own
  experiment, non-blind by its own admission, and it lives outside git at
  `~/Data/ym-formal/formalizer-run1/`. Validating a checker against it inherits that. Needs
  import into the repository as a fixture and a second independent grader, or re-grading by derivation in
  Lean as was done for three of the attempts.
- **G6. Mutants are not automatically false.** YM-T17 found 2 of 54 planted mutants were **true** statements
  (`‖·‖² = 0 → ‖·‖³ = 0`; `r² ≤ r → r³ ≤ r` on `[0,1]`), making their detections false alarms. Every mutant
  needs its own ground-truth label, not an assumed one.
- **G7. Who writes the paper-side predicate is unspecified.** A human writer contradicts the stated motive
  ("really can't rely upon me"); a model writer makes the tool's accuracy a measured property, per YM-T17.
  Blinding and model-family separation (see non-functionals) must be fixed before any run.
- **G8. The Lean export mechanism and the image pin are open.** The elaborated-term JSON format, the
  metaprogram's place in the build, and step 2's decision D1 (which image) all remain; the Mathlib revision
  changes `Expr` shapes and must be pinned by digest in the receipt.
- **G9. No acceptance thresholds.** Pre-registration needs numbers: the detection rate on mutants and the
  false-alarm rate on the 18 known-faithful statements at which the checker may gate the fidelity audit.
- **G10. Ownership and placement.** New crate versus extending `qft-symbolics` versus the planned
  `Rust/lean-graph`; and no `BOARD.md` packet exists. The reusable-library policy says the capability is
  implemented once behind a public API with Yang-Mills supplying fixtures only — so the checker belongs in a
  library, not in `docs/plans/`, where `formalizer.py` and `prover.py` currently sit.
- **G11. CUDA is not implied here.** "Rust / CUDA C++" in the request is the general porting direction. The
  rule is a Rust CPU reference path behind a public API, with CUDA added only for a measured data-parallel
  bottleneck. Predicate evaluation on a few thousand instances is not one. Do not start CUDA in this unit.

## For step 2

Decide G3 (fragment) and G7 (who writes the predicate, and blinding) before anything else: they determine
whether the tool can give a verdict at all and whether its verdict means anything. Put G1, G2, G9 and G10 to
Ernest as four questions with recommendations. G4, G5 and G6 are the experiment design and belong in the
pre-registration, which is written before any generation.
