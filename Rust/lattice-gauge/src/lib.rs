//! Exact lattice gauge theory on periodic tori.
//!
//! The Z2 and geometry modules use integer arithmetic; the SU(2) modules use f64. Conventions follow the Wilson lattice
//! gauge theory of `documents/mono-notes/topics/yang-mills-lattice.tex`: sites of
//! `(Z/n_0) x ... x (Z/n_{d-1})`, links `(site, direction)` from `site` to
//! `site + e_direction`, plaquettes `(site, first < second)`. Directions are numbered from 0.
//! Computed evidence is not a formal proof; see README.md.
/// Modified Bessel functions and the SU(2) character expansion; exact 2D results.
pub mod character;
/// GPU path for the SU(2) update (feature `cuda`).
#[cfg(feature = "cuda")]
pub mod cuda;
/// Sets of plaquettes modulo 2: twist surfaces, monopole counts, coboundary classes.
pub mod cochain;
/// Zero-momentum 0++ glueball operators (APE smearing) and correlators: the lattice mass gap.
pub mod glueball;
/// Maximal abelian gauge, abelian projection and DeGrand-Toussaint monopoles.
pub mod monopole;
/// Philox4x32-10 counter-based random numbers, reproducible on CPU and GPU.
pub mod philox;
/// Purity of the trace of a nonnegative transfer operator, in exact integer arithmetic.
pub mod spectral;
/// Jackknife errors and integrated autocorrelation times.
pub mod stats;
/// SU(2) as unit quaternions: Haar measure, Kennedy-Pendleton heat bath, over-relaxation.
pub mod su2;
/// Topological charge (clover) and cooling for SU(2) in d = 4: instantons.
pub mod topology;
/// Sites, links, plaquettes and cubes of a periodic torus.
pub mod torus;
/// Exact enumeration of Z2 lattice gauge theory on small tori.
pub mod z2_gauge;
/// Compact U(1) lattice gauge theory: the counter-model (Coulomb phase in 4D).
pub mod u1;
/// Centre vortices: direct maximal centre gauge, centre projection, P-vortices.
pub mod vortex;
/// SU(2) Wilson lattice gauge theory: checkerboard Monte Carlo and Wilson loops.
pub mod wilson;
