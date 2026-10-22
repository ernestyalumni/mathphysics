# Current Execution Queue — TMP / MHV Discovery Engine

**Last updated:** 2026-10-21

## Implemented: Frankel SO(3)/SU(2) tooling (October 17 request)

- [Rust/lie-groups](Rust/lie-groups/README.md) and
  [CppMath/LieGroups](CppMath/LieGroups/README.md): generic matrix storage,
  distinct validated groups/algebras, f32/f64/software binary16 Rust backends,
  exact symbolic matrices/real charts and normalization-aware exponentials.
- Rust: 27 integration tests +4 doctests; C++: 13 named contract groups +epsilon
  regression, debug/release/ASan+UBSan. Actual SymPy/Sage/Cadabra differential
  checks and 25 independent numerical references; original QFT suite rerun.
- Original matrix implementations and mono-notes kept untouched.
  Scope, limitations and evidence: Rust/lie-groups/docs/VERIFICATION.md.
- Optional continuation: Rust/lie-groups/CONTINUE.md. No background job running.

## Implemented: Goldstein rotation dictionaries (October 21)

- Extended Rust/lie-groups and CppMath/LieGroups, added Python/rotations,
  and condensed mono-notes Goldstein §§4.6–4.7. Reusable numerical maps,
  explicit inverse branches, exact symbolic identities and independent
  cross-language/CAS regression all pass.
- Receipt: Rust/lie-groups/docs/GOLDSTEIN-VERIFICATION.md. No background test
  process is required to interpret this completed receipt.

## Implemented: interactive ideal six-DOF lab (October 21)

- `./run_rigid_body_lab.sh`: local Python backend plus WebGL viewer at
  <http://127.0.0.1:8765>, keyboard/touch force and torque, free coasting,
  tumbling/flip presets, live conservation telemetry, pause and single step.
- Reusable C++/Rust/Python dynamics share the body-to-inertial convention;
  implicit midpoint momentum with Cayley attitude, full COM inertia.
- Analytic, DOP853, cross-language, symbolic and actual-browser checks pass.
  Scope, commands and evidence:
  [six-DOF verification](JavaScript/rigid-body-viewer/VERIFICATION.md).
- Concise theory and engineering reading added to both mono-notes PDFs.
  The interactive server remains running by user request; test jobs finished.

## Active: Rust QFT corpus (Ernest's October 16 request)

- Two bounded slices implemented at [Rust/qft-symbolics](Rust/qft-symbolics/README.md):
  exact algebra, strict typed tensor/index kernel, 13 named textbook entries
  (partial chapters), now 55 integration tests plus two compile-fail doctests
  (including the shared polynomial-conjugation regression),
  and 200 exact SMGA recurrence comparisons. New: formal-d 47.2, angular 50.2a,
  evaluated 50.4, and live Cadabra positive/negative tensor differential checks.
- [Research-workflow findings](agent/research/2026-10-16-theophyseval.md) from
  the supplied TheoPhysEval paper are recorded with source/provenance limits.
- [Continuation units](Rust/qft-symbolics/CONTINUE.md): arbitrary precision,
  remaining Srednicki gaps, nonabelian vertices, ordinary Yang–Mills
  BG, then independently verified MHV/BCFW amplitudes.
- The Rust port is specifically requested; retain existing Python as an
  independent oracle until full parity. No chapter or whole-port completion
  is inferred from selected tests.

## Historical July application queue

**Deadline driving this historical queue:** YC application deadline 2026-07-27, 8: 00 PM PT
**Detailed task packets:** [`tasks/README.md`](tasks/README.md)

## Do now, in order

1. ~~**TASK-000 — baseline reproducibility**~~ DONE 2026-07-25. `uv run pytest` now works
   from a fresh checkout with no ad hoc flags: focused suite 43 passed, full suite 64 passed
   (the previously-noted WIP failures are gone — committed in `619b505`).

2. ~~**TASK-009 — two-minute verification demo**~~ DONE 2026-07-25. See
   [`amplitudes/DEMO.md`](amplitudes/DEMO.md) and
   `amplitudes/scripts/09_verification_demo.py`.

3. **TASK-002 / TASK-003 — durable result reports**
   Convert the already-working checks into human- and machine-readable artifacts with seeds, versions, and explicit scope.

4. **Application review**
   `yc/` contains public-safe product language and video rehearsals. Authenticated form state and private founder/company facts stay outside this public research repo.

## Parallel only if it cannot delay the critical path

- **TASK-005 — verifier wrapper:** useful if it can be demonstrated quickly with a known-true and known-false claim.
- **TASK-008 — public research/product README:** polish after the demo behavior is frozen.

## After submission

- **TASK-004:** amplitude consistency checks.
- **TASK-006:** adjacent-region exploration.
- **TASK-010:** spinor-helicity eval suite v0.
- **TASK-007:** NMHV pattern search.

## Current evidence

```text
Focused amplitude suite: 43 passed, no warnings
Full suite:              64 passed, 0 failures
Script 08:               OVERALL: PASS
Demo (script 09):        correct candidate PASS n=3..10; broken candidate correctly FAILs
```

Reproduction command, verified 2026-07-25 (no ad hoc `--with`/`--no-project` flags needed
anymore — see `tasks/TASK-000-baseline-reproducibility.md`):

```bash
cd amplitudes/MHVamplitudes
uv run pytest -q tests/unit_tests/amplitudes tests/integration_tests/test_smga_reproduction.py
```

Two-minute demo command, from repo root:

```bash
uv run --project amplitudes/MHVamplitudes python amplitudes/scripts/09_verification_demo.py
```

## Agent rules

- Claim one task at a time using the protocol in `tasks/README.md`.
- Never commit on `master`.
- Do not stage the pre-existing edits in `spinors/invariants.py` or `spinors/massless_spinors.py` unless Ernest explicitly assigns that work.
- A generated report must include command, git SHA, branch, seed, backend, and pass/fail criteria.
- Publishing, arXiv submission, application submission, and external messages require Ernest's approval.
