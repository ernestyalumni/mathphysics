# TASK-001 — Berends-Giele numeric recursion, hardened + tested

Status: TODO
Priority: P0 (blocks TASK-002, TASK-003)
Estimated size: one session (3–5 h)
Branch suggestion: `feat/task-001-berends-giele`

## Goal

A trustworthy, tested numerical implementation of Berends-Giele off-shell recursion for
color-ordered tree-level gluon amplitudes, arbitrary helicity configuration, n ≤ 10, in the
`amplitudes/MHVamplitudes/` package. This is the workhorse every downstream verification
task calls.

## Inputs / starting points

- Existing package: `amplitudes/MHVamplitudes/src/mhvamplitudes/` (spinors module exists;
  note `spinors/invariants.py` and `spinors/massless_spinors.py` have uncommitted edits on
  master — inspect and either adopt or reconcile them).
- Existing scripts: `amplitudes/scripts/02_spinor_construction.py` through
  `08_smga_stripped_amplitudes.py` — reuse, don't rewrite.
- Reference: Dixon arXiv:1310.5353 §2–3 (color decomposition, off-shell currents);
  Srednicki conventions per repo default (mostly-plus).

## Steps

1. Implement (or harden the existing) `berends_giele(momenta, helicities)` returning the
   color-ordered partial amplitude as a complex number, using complex spinor-helicity
   variables so complex/analytically-continued kinematics are supported.
2. Random-momentum generator: massless, on-shell, momentum-conserving n-point phase-space
   points (complex deformation optional flag), seedable for reproducibility.
3. Tests (pytest, in `amplitudes/MHVamplitudes/tests/`):
   - MHV configurations match Parke-Taylor ⟨ij⟩⁴/(⟨12⟩⟨23⟩…⟨n1⟩) to 1e-10 relative
     error for n = 4…8, ≥ 20 random points each.
   - All-plus and single-plus amplitudes vanish (to numerical zero) for n = 4…8.
   - Cyclic invariance of the color-ordered amplitude at random points.
   - 4- and 5-point known answers cross-checked against `amplitudes/scripts/04_*.py`, `06_*.py`.

## Definition of done

`cd amplitudes/MHVamplitudes && uv run pytest tests/ -k berends` passes, and a short note
in `TESTING.md` documents how to run it and observed tolerances.

## Out of scope

Gravitons, loops, symbolic (Cadabra2) versions — later tasks.
