//! Algebraic topology for the Yang–Mills mass-gap program (YM-T19), std only, exact arithmetic.
//!
//! - `smith`: Smith normal form invariant factors of integer matrices, and rank modulo a prime.
//! - `complex`: finite chain complexes with sparse integer boundary maps; homology over Z (Betti numbers and
//!   torsion) and over Z_p.
//! - `cubical`: the cubical complex of a periodic d-torus with the site conventions of `Rust/lattice-gauge`.
//! - `simplicial`: simplicial complexes from facets, Z_p cohomology bases, the Alexander–Whitney cup product,
//!   and the Freudenthal triangulation of a d-torus with integer seam cocycles.
//! - `surface`: components, Euler characteristic and orientability of square-complex surfaces (P-vortex sheets).
//!
//! Computed evidence is not a proof. Label: unreviewed.
pub mod complex;
pub mod cubical;
pub mod simplicial;
pub mod smith;
pub mod surface;
