# YM-T05 — Lattice geometry with exact arithmetic (Rust)

Tier T2 · Workstream B · Depends on: PDD intent approved · New crate `Rust/lattice-gauge` · Review: T2

## Behaviour wanted
A d-dimensional periodic lattice with side lengths (L₀,…,L_{d−1}): typed indices for sites, oriented links, plaquettes, cubes; boundary and coboundary maps over the integers and over Z₂; iteration over the links of a plaquette in a fixed documented orientation; a spanning-tree (maximal tree) gauge; centre-twisted boundary conditions specified by an antisymmetric integer matrix n_{μν} mod N.

## Fixtures
- Counts: V sites, dV links, d(d−1)/2·V plaquettes.
- coboundary∘coboundary = 0 and boundary∘boundary = 0 on random integer cochains.
- A maximal tree on a connected periodic lattice has V − 1 links, leaving (d−1)V + 1 free: in d = 4, **3V + 1**, the count used in the proof's determinant budget (`ARM/summary.md`, "No loss of log β per site").
- Twisted boundary conditions: the product of plaquette twist signs over any closed 2-surface of the torus equals the prescribed flux.
- Negative control: a non-closed surface gives a flux that depends on its boundary.

Conventions go in the crate README and in `PLAN/notes/conventions.md`.
