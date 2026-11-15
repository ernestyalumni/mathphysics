# lattice-gauge

Exact lattice gauge theory on periodic tori, in integer arithmetic. No dependencies.

Part of the Yang–Mills program (`docs/plans/yang-mills-mass-gap/`). The mathematics, with labels cited below, is in `documents/mono-notes/topics/yang-mills-lattice.tex`.

| Module | What it does | Checks which statement |
|---|---|---|
| `torus` | Sites, links, plaquettes, cubes of `(Z/n_0) x … x (Z/n_{d-1})`; stable index conventions | Cell counts; δδ = 0; first Betti number of the torus |
| `cochain` | Sets of plaquettes mod 2: the 't Hooft twist surface, monopole counts, membership in the span of link coboundaries by Gaussian elimination over GF(2) | Lemma `ym:lem:thooft-class` on 4D tori, by two independent methods |
| `z2_gauge` | Exact density of states of Z₂ gauge theory with twists on tori with at most 28 links | 2D closed form; Lemma `ym:lem:twist-move` |
| `spectral` | Purity bounds for a nonnegative transfer operator on integer spectra | Proposition `ym:prop:purity` |

## Run

    export CARGO_TARGET_DIR=/tmp/mathphysics-rust-target
    cargo test --offline --locked          # 19 tests; the 2x2x2 enumeration takes a few seconds in debug mode

## What this is not

A test on finitely many lattices or spectra is evidence, not a proof. The proofs are in the notes. Nothing here is a Monte Carlo simulation, and nothing here treats SU(2) yet: see roadmap units T06–T08.

## Conventions

Directions are numbered from 0 (the notes number them from 1). A link `(site, direction)` runs from `site` to `site + e_direction`. A plaquette `(site, first < second)` has the four links of `U(y,μ) U(y+e_μ,ν) U(y+e_ν,μ)⁻¹ U(y,ν)⁻¹`. Index formulas are in the module documentation of `torus` and must not change without updating every consumer.
