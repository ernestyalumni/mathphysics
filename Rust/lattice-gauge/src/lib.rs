//! Exact lattice gauge theory on periodic tori.
//!
//! Everything in this crate is integer arithmetic. Conventions follow the Wilson lattice
//! gauge theory of `documents/mono-notes/topics/yang-mills-lattice.tex`: sites of
//! `(Z/n_0) x ... x (Z/n_{d-1})`, links `(site, direction)` from `site` to
//! `site + e_direction`, plaquettes `(site, first < second)`. Directions are numbered from 0.
//! Computed evidence is not a formal proof; see README.md.
/// Sets of plaquettes modulo 2: twist surfaces, monopole counts, coboundary classes.
pub mod cochain;
/// Purity of the trace of a nonnegative transfer operator, in exact integer arithmetic.
pub mod spectral;
/// Sites, links, plaquettes and cubes of a periodic torus.
pub mod torus;
/// Exact enumeration of Z2 lattice gauge theory on small tori.
pub mod z2_gauge;
