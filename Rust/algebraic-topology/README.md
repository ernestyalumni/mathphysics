# algebraic-topology

Exact algebraic topology for the Yang–Mills mass-gap program (task YM-T19). Std only, no dependencies.

- `smith`: invariant factors (Smith normal form) of integer matrices; rank mod p.
- `complex`: chain complexes with sparse integer boundaries; homology over Z (Betti numbers and torsion) and Z_p.
- `cubical`: cubical complex of the periodic d-torus, with the site conventions of `Rust/lattice-gauge`.
- `simplicial`: simplicial complexes, coboundary, Alexander–Whitney cup product, Z_p cohomology bases, and the
  Freudenthal triangulation of T^d with integer seam cocycles.
- `surface`: components, Euler characteristic, closedness, branching and orientability of square surfaces
  (P-vortex sheets).

Tests are independent oracles from known topology: H_*(T^d), S^2, RP^2 (torsion Z_2), Klein bottle (Z + Z_2),
the 7-vertex torus; w cup w != 0 on RP^2 and a cup b != 0 on T^2 mod 2; and 't Hooft's twisted-torus formula
int n cup n = 2 Pf(n) on T^4, computed from first principles.

Computed evidence is not a formal proof. Label: unreviewed.
