# YM-A01 — Independent clean rebuild on Linux x86-64

Tier T2 · Workstream A · Depends on: YM-O01 · Blocks: G1 · Review: T2 · Budget: one day wall time (the group's build took 184 minutes on Apple Silicon; expect longer on 12 cores, 31 GB RAM)

## Do
1. Fresh copy of `LEAN` with **no** unpacked project artifacts. Only the official Mathlib cache.
2. `lake build WilsonReference WilsonConstruction` from source. Record peak memory; if 31 GB is not enough, lower parallelism (`-j`) and say so.
3. Run `YMG/submissions/openai-20261028/scripts/verify.py integrity`, then the full `verify.py all` as far as it runs on Linux (read `docs/rerun.md` first).
4. Compare: number of modules built, warning counts by group, declaration counts, against assessment §9.1 (9,470 modules; 2,453 warnings; 250,319 declarations; 100,168 reached).

## Deliverable
`PLAN/notes/YM-A01-linux-rebuild.md`; logs under `PLAN/ledger/evidence/YM-A01/`.

## Acceptance
A table of our numbers against the assessment's, every difference explained or flagged `FINDING`. Zero errors, or the first error quoted in full.
