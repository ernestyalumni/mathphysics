# YM-T19: algebraic topology and the mass gap — theory map and first computations

Label: unreviewed (claude-opus55-strat, 2026-11-21). Code: `Rust/algebraic-topology` (this branch) and `Rust/lattice-gauge` (branch `ym/YM-N01-claude-opus55-su2`). Citations are from memory; check before quoting.

## Where algebraic topology enters the 4D Yang–Mills gap problem

**1. Twist sectors are cohomology classes.** On the torus T⁴, SU(N)/Z_N bundles, i.e. 't Hooft's twisted boundary conditions, are classified by w ∈ H²(T⁴; Z_N), with six components n_{μν}.
- 't Hooft (Nucl. Phys. B153, 1979; Commun. Math. Phys. 81, 1981): the instanton number of a twisted bundle is fractional, Q ≡ −P(w)/(2N) mod 1, where P is the Pontryagin square. For translation-invariant twists P(w)[T⁴] = 2 Pf(n) with Pf(n) = n₀₁n₂₃ − n₀₂n₁₃ + n₀₃n₁₂.
- Electric and magnetic flux free energies are Fourier duals over these sectors. Confinement is the statement that magnetic flux costs no free energy at large volume, Z_tw/Z → 1, while electric flux costs energy proportional to the length.
- How fast Z_tw/Z → 1 is controlled by the same correlation length as the gap. This is the YM-N02 measurement, and the estimator exists (`examples/twist_z.rs`).

**2. Generalized symmetry and anomalies.** In modern language the centre symmetry is a 1-form Z_N symmetry with a 2-form background B ∈ H²(X; Z_N).
- Gaiotto, Kapustin, Seiberg, Willett, *Generalized global symmetries* (JHEP 2015): confinement means the 1-form symmetry is unbroken.
- Gaiotto, Kapustin, Komargodski, Seiberg, *Theta, time reversal and temperature* (JHEP 2017): at θ = π, SU(2) Yang–Mills has a mixed anomaly between time reversal and the 1-form symmetry, proportional to the Pontryagin square ∫ P(B)/2. The theory therefore cannot be trivially gapped at θ = π. Either the gap closes, or the vacuum breaks T, or there is topological order.

These are rigorous constraints on which kinds of gapped phase are possible. Any proof of a gap at θ = 0 must use something that fails at θ = π: a concrete filter for arguments, complementing the U(1) counter-model filter.

**3. Topological charge on the lattice.** M. Lüscher, *Topology of lattice gauge fields* (Commun. Math. Phys. 85, 1982), reconstructs a principal bundle from lattice fields whose plaquettes are close enough to 1 ("admissible"), so the integer Q is a genuine characteristic class. A rigorous continuum limit must pass through such admissible configurations. Instantons are not the gap mechanism, but the bundle topology they measure is part of the continuum object whose spectrum is in question.

**4. Vortices and monopoles are cycles.** After centre projection, P-vortex sheets are closed Z₂ 2-cycles on the dual lattice. After abelian projection, monopole worldlines are integer 1-cycles. Two kinds of invariant matter:
- the homology classes in the ambient torus, tied to twist sectors and winding;
- the intrinsic topology of the sheets (orientability, H₁ torsion, Euler characteristic).

The vortex picture of confinement (Engelhardt, Faber, Greensite, Langfeld, Olejník, Reinhardt) says confinement is the percolation of vortex sheets, which are non-orientable. The rigorous analogue in 3D compact U(1) is Göpfert–Mack, where the gap comes from a monopole plasma.

## First computations (all reproducible with the code named)

**'t Hooft's formula from first principles** (`tests/thooft.rs`). On a Freudenthal triangulation of T⁴ (sides 3), the seam 1-cocycles x_μ are integer cocycles, and x₀x₁x₂x₃ evaluates to ±1 on the fundamental cycle. For the twist class n = Σ n_{μν} x_μ ∪ x_ν, the integral of n ∪ n equals 2 Pf(n) for every tested n, including non-unit entries. The fractional instanton number −Pf(n)/N mod 1 thus follows from our own cup product, not from a textbook.

**Vortex-sheet homology** (`examples/vortex_surfaces.rs`, data from lattice-gauge `vortex_scan` after direct maximal centre gauge).
- At β = 2.4 on 8⁴ (10 configurations), the P-vortex sheet has 9–17 components, and one percolating component carries 89–96% of the area.
- Every component is closed, as the Bianchi identity requires; this is an end-to-end check.
- Branching edges are frequent, so naive orientation fails. Homology decides anyway: b₂(Z₂) − b₂(Q) = 1 in all 10 configurations. In a 2-complex this difference counts the Z₂-torsion summands of H₁. The sheet contains many orientable closed sub-surfaces plus exactly one Z₂ obstruction, a homological "crosscap".
- **Control.** Sheets from independent random Z₂ links at the same density on 8⁴ have 115–120 small components and b₂(Z₂) − b₂(Q) = 0 (4 of 4 samples). The obstruction is therefore a property of the dynamics, not of the construction.
- **Phase and volume dependence.** Configurations with the obstruction b₂(Z₂) − b₂(Q) = 1:

  | β | 8⁴ | 12⁴ |
  |---|---|---|
  | 2.2 | 20/20 | |
  | 2.3 | 20/20 | |
  | 2.4 | 10/10 | |
  | 2.5 | 7/20 (largest share 0.60) | **20/20** (largest share 0.90) |
  | 2.6 | 4/20 (largest share 0.42) | **14/20** (largest share 0.68) |
  | 2.8 | 1/19 | |

  On 8⁴ the obstruction switches off across the finite-volume deconfinement crossover (β_c ≈ 2.51 for N_t = 8), together with percolation. At β = 2.6 it returns on 12⁴, which is still confined (β_c ≈ 2.64 for N_t = 12). The invariant therefore tracks confinement, not the coupling. That is the behaviour a topological order parameter for confinement should have.
- Candidate lemma to formalize: in the confined phase the centre-projected vortex sheet has a percolating component whose H₁ carries Z₂ torsion; in the deconfined phase it does not.

**Monopole worldlines** (lattice-gauge `monopole::worldline_topology`): integer 1-cycles on the dual lattice. Current conservation holds on every slice, which is a test. Compact U(1) on 6⁴ percolates in its confined phase (largest cluster 0.99 of the current at β = 0.8) and not in its Coulomb phase (0.22 at β = 1.2). MAG-projected SU(2) at β = 2.3 percolates (0.72). Full scans are chained on the Mac (SU(2)) and the laptop (U(1)).

**Instantons by persistent homology** (`persistence`, `examples/instanton_persistence.rs`, cooled charge densities):
- On 8⁴ at β = 2.4, the persistent peaks of q and −q give N_I − N_A = Q in 12/12 configurations, including one instanton–anti-instanton pair at Q = 0.
- On 12⁴ this holds in only 8–13/30, both for peak-height persistence and for steepest-ascent basins with integrated charge. Lumps overlap and partly cancel, so the dilute-instanton-gas picture fails at this density. This agrees with instantons not being the gap mechanism.

**Steenrod cup-1 and the Pontryagin square** (`tests/steenrod.rs`): cup-1 satisfies its coboundary identity on random cochains, and ∫P(w) mod 4 on T⁴ is a class invariant equal to 2 Pf(n) for constant twists. The plain cup square is not invariant for non-cocycle lifts, so the test is non-vacuous. Twist sectors of any representative can now be evaluated.

**Centre dominance** (same configurations). The centre-projected Creutz ratio χ_cp(3,3) = 0.1357(28) at β = 2.3 reproduces the full string tension σa² = 0.1359 (Teper, hep-th/9812187, Table 10). At β = 2.4 it gives 82% of σ.

**Instantons** (lattice-gauge `topology`). After 20 cooling sweeps the clover charge clusters near 0.86 × integers. The topological susceptibility at β = 2.4 on 12⁴ is aχ_t^{1/4} = 0.1154(8), against Teper's 0.1163(13) (Table 15).

## What this suggests for the program

- The conditional-gap chain (YM-G02, YM-F01) is purely analytic. The twist sectors give it a topological handle: the defect Δ_β(n) compares Z on n⁴ and (2n)⁴ tori in the trivial sector only. A version that sums or compares twist sectors would test the same physics, the 1-form symmetry being unbroken, that the anomaly argument constrains.
- A gap argument must distinguish θ = 0 from θ = π. A step that does not see θ cannot be the whole argument.
- The vortex obstruction (b₂(Z₂) − b₂(Q) = 1, plus percolation) is a candidate topological order parameter. If it switches off where the finite torus deconfines, it is worth trying to state as a lemma: in a confining phase, the centre-projected sheet carries a non-orientable percolating component.

## Next

1. Vortex obstruction on 16⁴ at β = 2.6 and 2.7 (running; sparse rank now makes this cheap).
2. Monopole-worldline scans across β for SU(2) and U(1) (chained on the Mac and laptop).
3. Twist-sector free energies (YM-N02) with P(w) to label sectors.
4. Lüscher admissibility: the fraction of admissible configurations against β, and Q from bundle reconstruction compared with Q from cooling.
