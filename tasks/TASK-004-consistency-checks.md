# TASK-004 — Consistency checks: cyclicity, KK relations, soft & collinear limits

Status: TODO
Priority: P1
Depends on: TASK-002
Estimated size: one session (3–5 h)
Branch suggestion: `feat/task-004-consistency`

## Goal

A reusable `consistency` test module that any amplitude implementation in this repo must
pass. These are the field-theory identities that catch silently-wrong amplitudes — the
core of the "verifier" story.

## Checks to implement (each parametrized over n and helicity config)

1. **Cyclicity:** A(1,2,…,n) = A(2,…,n,1) for color-ordered partial amplitudes.
2. **Reflection:** A(1,2,…,n) = (−1)ⁿ A(n,…,2,1).
3. **U(1) decoupling / photon decoupling:** Σ_cyclic A(1,σ(2,…,n)) = 0.
4. **Kleiss-Kuijf relations** (at least the basic subcyclic ones for n = 5, 6).
5. **Soft limit:** as gluon s becomes soft, A_n → Soft(s) × A_{n−1} with the standard
   eikonal factor; verify the ratio numerically along a soft trajectory.
6. **Collinear limit:** leading splitting-function behavior as two adjacent momenta become
   collinear (compare against Dixon 1310.5353's split amplitudes).

## Steps

1. New module `amplitudes/MHVamplitudes/src/mhvamplitudes/consistency.py` with one function
   per check returning `{status, max_violation, details}`.
2. Pytest suite exercising all checks against the Berends-Giele engine (TASK-001) for
   n = 4…7, random seeds fixed.
3. Existing scripts `07_soft_gluon_limit.py` should be refactored to call the module
   (keep the script as a thin CLI).

## Definition of done

`uv run pytest tests/ -k consistency` passes; `amplitudes/results/consistency_report.md`
is generated summarizing max violations per check (all < 1e-8, soft/collinear ratios
converging at the expected rate as the limit parameter → 0).
