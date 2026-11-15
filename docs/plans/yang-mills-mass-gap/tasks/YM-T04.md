# YM-T04 — `lean-graph`: make the million-line proof queryable

Tier T2 · Workstream B · Depends on: PDD intent approved; phase 2 on YM-O01 · New crate `Rust/lean-graph` · Review: T2

## Behaviour wanted (phase 1: published tables only)
Load `ARM/data/units.tsv`, `declaration-units.tsv`, `key-statements.tsv` and `LEAN/verification/declaration-usage.tsv`, `module-usage.tsv`. CLI with JSON output:
- `unit-of <declaration>`; `declarations-in <unit>`; `key-statements <unit>`;
- `path-from-headline <declaration>` using the dependency path column of `declaration-usage.tsv`;
- `where <declaration>` → file and line span; `stats` → counts.

## Fixtures
Counts stated in the sources: 32 units; 989 key statements; 250,319 imported project declarations of which 100,168 reached by the 17 headline results. Negative control: an unknown declaration returns a structured "not found", exit code non-zero.

## Phase 2 (after YM-O01)
Direct dependency edges from a Lean export, and `depends-on <declaration> --depth k`.

## Note
Read the header rows of the TSV files before designing types; column meanings are documented in `ARM/data/README.md`.
