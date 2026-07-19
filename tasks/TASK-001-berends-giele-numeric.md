# TASK-001 — Harden the existing specialized SMGA recurrence

Status: PARTIAL — implementation and integration tests already exist
Priority: P1
Depends on: TASK-000
Estimated size: one session (3–5 h)
Branch suggestion: `feat/task-001-berends-giele`

## Goal

Audit and harden the existing **specialized SMGA half-collinear recurrence** for stripped
single-minus amplitudes through n ≤ 10. Do not turn it into a generic Lorentzian,
arbitrary-helicity Berends–Giele engine in this task; the current implementation explicitly
does something narrower.

## Inputs / starting points

- Existing implementation: `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/berends_giele.py`.
- Existing integration comparison: `tests/integration_tests/test_smga_reproduction.py`.
- Pre-existing uncommitted spinor edits are unrelated. Preserve and do not stage them.
- Existing scripts: `amplitudes/scripts/02_spinor_construction.py` through
  `08_smga_stripped_amplitudes.py` — reuse, don't rewrite.
- Primary reference: arXiv:2602.12176v2, especially the recurrence and Appendix B.
- Background reference: Dixon arXiv:1310.5353 for standard amplitude conventions;
  Srednicki mostly-plus remains the repo default where applicable.

## Steps

1. Map every implemented recurrence component to the exact SMGA equation and notation.
2. Add boundary/degeneracy behavior tests for zero brackets and chamber walls.
3. Confirm deterministic results for fixed seeds through n = 10.
4. Profile n = 8..10 and document practical runtime and recurrence growth.
5. Make error messages distinguish invalid kinematics, chamber-boundary points, and code failure.

## Definition of done

The focused recurrence/integration suite passes; equation mapping and runtime notes are in
`TESTING.md`; and no unrelated spinor work is included in the commit.

## Out of scope

Generic Lorentzian BG currents, arbitrary helicities, gravitons, loops, and symbolic
Cadabra2 recursion.
