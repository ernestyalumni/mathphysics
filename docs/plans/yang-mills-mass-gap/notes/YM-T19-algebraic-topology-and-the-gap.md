# YM-T19: algebraic topology and the mass gap — theory map and first computations

Label: reviewed with corrections by Codex (YM-O03, `notes/YM-O03-codex-review-2026-11-21.md` §5); corrections applied 2026-11-21. Author claude-opus55-strat. Code: `Rust/algebraic-topology` and `Rust/lattice-gauge`. Citations are from memory except where the review checked them.

## Where algebraic topology enters the 4D Yang–Mills gap problem

**1. Twist sectors.** On T⁴ a class w ∈ H²(T⁴; Z_N), with six components n_{μν}, is the obstruction to lifting a PSU(N) bundle to SU(N); 't Hooft's twisted boundary conditions fix it. Fixing w still leaves an integer instanton label: w = 0 allows different second Chern numbers.
- Take an integral lift u of w on T⁴ with an oriented basis. Then ∫u² = 2(n₀₁n₂₃ − n₀₂n₁₃ + n₀₃n₁₂) = 2 Pf(n), and in the charge-sign convention used here the fractional charge is −Pf(n)/N mod 1.
- Changing the lift by N v changes (∫u²)/(2N) by an integer, because the intersection form of T⁴ is even. This also fixes the lift convention for odd N.
- The relation of this square to the gauge characteristic number is bundle theory ('t Hooft, Commun. Math. Phys. 81, 1981). A test of the cup product alone cannot derive that normalization.
- Flux free energies: the sector partition functions (Boltzmann weights, not their logarithms) undergo a discrete Fourier projection over temporal twists at fixed spatial magnetic twist, and the logarithms taken afterwards define electric and magnetic flux free energies. Spatial twists need separate interpretation.
- In a confining phase electric flux can cost energy proportional to length. Relating the approach Z_tw/Z → 1 to the lightest glueball correlation length needs an additional argument relating distinct sectors and observables; the estimator (`examples/twist_z.rs`) supplies none (see von Smekal and de Forcrand, hep-lat/0209149).

**2. Generalized symmetry and anomalies.** The centre symmetry is a 1-form Z_N symmetry with a 2-form background B. In pure SU(N), an area law for lines charged under it signals its unbroken realization (Gaiotto, Kapustin, Seiberg, Willett, JHEP 2015).

At θ = π, SU(2) Yang–Mills has a mixed anomaly between time reversal and the centre symmetry, with background phase exp(±iπ∫P(B)/2) and P(B) taken mod 4 (Gaiotto, Kapustin, Komargodski, Seiberg, JHEP 2017, §2). Anomaly matching excludes a *unique, symmetry-preserving, trivially gapped* vacuum at θ = π. Two time-reversed vacua, each with a positive excitation gap, are allowed. So a proof that bounds the gap at both angles can be consistent.

The resulting filter is narrower than an earlier version of this note claimed. It applies to arguments that also fix the symmetry realization and rule out vacuum degeneracy and topological order. The conditional-gap chain rests on reflection positivity, which fails at θ ≠ 0 (`notes/YM-G02-filters.md`).

**3. Topological charge on the lattice.** M. Lüscher, *Topology of lattice gauge fields* (Commun. Math. Phys. 85, 1982), shows that a small-plaquette (admissibility) condition is *sufficient* for a geometric lattice charge that is an integer characteristic class. It is not a necessity theorem for every continuum construction.

**4. Vortices and monopoles as cycles.** After centre projection with periodic links z, the plaquette cochain is δz. The full dual P-vortex sheet is therefore a closed mod-2 2-cycle whose ambient homology class vanishes, even when individual components wind. Intrinsic homology (orientability, H₁ torsion, Euler characteristic) is a different measurement.

After abelian projection, monopole currents are integer 1-cycles closed by lattice current conservation. Göpfert–Mack (Commun. Math. Phys. 82, 1982) concerns 3D compact U(1) with the Villain action. It does not establish a gauge-projected SU(2) percolation criterion.

## Computations (reproducible with the code named; replayed by the review unless noted)

**'t Hooft's formula as a cochain identity** (`tests/thooft.rs`). On a Freudenthal triangulation of T⁴ (sides 3), the seam 1-cocycles x_μ are integer cocycles, x₀x₁x₂x₃ evaluates to ±1, and for n = Σ n_{μν} x_μ ∪ x_ν the integral of n ∪ n equals 2 Pf(n) for every tested n. This checks the cup-product algebra; the gauge-theoretic normalization of the charge is the bundle-theory statement above.

**Pontryagin square** (`tests/steenrod.rs`). Cup-1 satisfies its coboundary identity on random cochains. For an integer lift of a **mod-2** cocycle, ∫P(w) mod 4 is invariant under changes of representative on T⁴ and equals 2 Pf(n) for constant twists. The routine does not implement general Z_N backgrounds.

**Vortex-sheet homology** (`examples/vortex_surfaces.rs`, data from lattice-gauge `vortex_scan` after direct maximal centre gauge).
- Every component is closed, as the Bianchi identity requires: an end-to-end check.
- What is measured: for a finite 2-complex H₂(Z) is free. By universal coefficients, b₂(F_p) = rank H₂(Z) + t_p(H₁), where t_p counts the cyclic factors of H₁ of p-power order. The code compares F₂ with F_p for p = 1 000 003, so the reported difference is t₂ − t_p.
  - A difference of 1 proves at least one even-order torsion factor in H₁ of the whole sheet.
  - It does not prove exactly one Z₂: torsion at p is not excluded, and Z/2, Z/4 and Z/8 are not distinguished.
  - With branching edges it does not imply a surface crosscap. An earlier version said "exactly one Z₂ obstruction, a crosscap"; that was an overstatement.
  - Settling the algebra needs an integral Smith form, component by component.
- At β = 2.4 on 8⁴ (10 configurations): 9–17 components, the largest carrying 0.894–0.962 of the area, and difference 1 in all 10. A largest-area fraction describes concentration of the sheet; it does not by itself test spanning or winding.
- Control: random Z₂ links at matched density on 8⁴ gave 115–120 components and difference 0 in 4/4 samples. The review did not replay these.
- Configurations with difference 1:

  | β | 8⁴ | 12⁴ | 16⁴ |
  |---|---|---|---|
  | 2.2 | 20/20 | | |
  | 2.3 | 20/20 | | |
  | 2.4 | 10/10 (default start); 18/20 (random start, different configurations) | | |
  | 2.5 | 7/20 (largest share 0.60) | 20/20 (0.90) | |
  | 2.6 | 4/20 (0.42) | 14/20 (0.68) | 19/20 (0.82) |
  | 2.7 | | | 15/20 (0.80) |
  | 2.8 | 1/20 | | 1/20 (0.33) |
  | 2.9 | | | 3/20 (0.39) |

  At β = 2.8 the scan processed 20 configurations and configuration 19 had no pierced plaquettes. The dump writes only pierced plaquettes, so it does not appear; an empty sheet has difference 0.
- The default-start and random-start runs at β = 2.4 are not paired configuration by configuration, so no per-configuration comparison is claimed.
- **Reading.** At fixed β the count rises with the volume (β = 2.6: 4/20, 14/20, 19/20 on 8⁴, 12⁴, 16⁴). On 16⁴ the drop between β = 2.7 and 2.8 brackets the N_t = 16 thermal coupling β_c(16) = 2.7395(100); on 8⁴ the drop between 2.4 and 2.5 is near, but does not bracket, β_c(8) = 2.5115(40) (Teper, Table 14). Thermal criticality is defined at large spatial volume with fixed temporal extent; symmetric boxes have no such singularity, and the small rebound at 2.9 argues against fitting a critical point (YM-O04 §4). The fall on 8⁴ sits near the N_t = 8 thermal coupling (≈ 2.51), a finite-size crossover on a symmetric torus. This is an *association* with confinement, supported by three box sizes, four random controls and one gauge-fixing prescription. It does not isolate confinement, define an infinite-volume order parameter, or establish a per-component or per-configuration lemma.

**Monopole worldlines** (lattice-gauge `monopole::worldline_topology`). Current conservation holds on every slice. The largest-cluster fraction of the current measures concentration; calling it percolation needs a winding or finite-size test, and no windings occurred in these runs.

| | largest-cluster fraction of the current |
|---|---|
| U(1) on 6⁴ | 0.99 at β = 0.8 (confined), 0.22 at β = 1.2 (Coulomb) |
| SU(2) after MAG, 8⁴ | 0.87, 0.73, 0.45, 0.26, 0.23 at β = 2.2–2.6 |
| SU(2) after MAG, 12⁴ | 0.55, 0.35, 0.18, 0.15 at β = 2.4–2.7 |

On 12⁴ the fraction drops below one half at β = 2.5, where the vortex count is 20/20. As measured here it is a weaker association than the vortex count.

**Instantons by persistent homology** (`examples/instanton_persistence.rs`, cooled charge densities). The code compares counts with round(1.15 Q_clover), using a persistence threshold (default 0.002) or a basin-charge cut. The agreement is 12/12 peaks on 8⁴ at β = 2.4, but on 12⁴ only 8/30 by peaks and 12/30 by basins. This tests a counting procedure. The failure may come from segmentation, cooling or overlapping lumps, and it cannot decide what generates the mass gap.

**Centre dominance.** The centre-projected Creutz ratio at β = 2.3 is χ_cp(3,3) = 0.1357(28) from the 200-configuration scan; the 20-configuration homology scan at the same parameters gives 0.1348(58). Teper's Table 10 gives a√σ = 0.3687(22), so a²σ ≈ 0.1359. At β = 2.4, 0.0584(32) is about 81% of a²σ = 0.0721 (a√σ = 0.2686(21)). A finite-loop ratio does not determine the asymptotic tension.

**Topological susceptibility.** After 20 cooling sweeps at β = 2.4 on 12⁴: aχ_t^{1/4} = 0.1075(7) from the raw clover charge, and 0.1154(8) after rescaling by 1.149 and rounding to integers. Teper's 0.1163(13) is one charge prescription in his Table 15; it is closer to the rescaled value, and the comparison is not normalization-independent.

## What this suggests for the program

- The conditional-gap chain is analytic, and its recurrence is now a theorem (YM-O03 §1, `notes/YM-N03-defect-measurement.md`). Its open content is a seed at a bounded physical size. Twist sectors offer a topological handle on that seed: the defect compares trivial-sector partition functions only.
- The vortex count is a candidate *association* with confinement. Before any lemma it needs integral invariant factors, an infinite-volume formulation, and control of the gauge-fixing prescription.

## Next

1. Integral Smith normal form of the vortex-sheet complexes, component by component, to replace t₂ − t_p.
2. Winding and finite-size tests for vortex and monopole clusters.
3. Twist-sector free energies (YM-N02, waiting on G04), with the Fourier projection done on Boltzmann weights.
4. Lüscher admissibility: the fraction of admissible configurations against β, and the geometric charge compared with the cooled charge.
