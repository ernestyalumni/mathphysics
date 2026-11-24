# YM-T10: exact Z₂ lattice gauge theory on tori with sides 1 and 2, and a GPU-versus-CPU test of exact enumeration

Label: unreviewed (claude-opus55-strat, 2026-11-23). Code: `Rust/lattice-gauge/src/z2_code.rs`, `cuda/z2_enum.cu`, and the examples `z2_exact`, `z2_defect`, `z2_duality`.

## The computation is a weight enumerator

With Z₂ links written additively in F₂:
- The plaquette field is δU, where δ: F₂^links → F₂^plaquettes is the coboundary. On a side of length 1, links that repeat in a plaquette cancel mod 2.
- Every plaquette configuration that occurs lies in the binary linear code C = im δ, of rank r, and each one is hit by 2^{links − r} link configurations.
- With the handbook energy Σ_p (1 − σ_p) and normalized counting measure:

  Z(β) = 2^{−r} Σ_w A_w e^{−2βw},   A_w = #{c ∈ C : |c| = w}.

- A twist t (reversed plaquette signs) replaces C by the coset t + C.

The integers A_w give Z exactly at every β. They are computed by a Gray-code walk over C: each step XORs one basis vector into a 128-bit word and counts its popcount.

**Checks.**
- Against the independent link-by-link enumerator `z2_gauge` on 2×2, 2×4, 3×3, 2×6 and 2×2×2.
- Against a direct walk over all link configurations on boxes with sides 1, including the twisted cosets.
- The MacWilliams identity between C and its dual code, verified exactly on 2²×1², 2³ and 2³×1.
- CPU and GPU histograms identical on every box tested and on the 2⁴ subcodes used for benchmarking (`tests/z2_code.rs`).

## GPU against CPU for exact enumeration

The question was whether GPUs help with symbolic or exact computation, given that they are poor at branching. This workload has no data-dependent branching: XOR, popcount, and a private histogram per thread.

| Hardware | Codewords per second | 2⁴ torus (2⁴⁵ codewords) |
|---|---|---|
| RTX 3060 (desktop), CUDA C++ | 3.1–3.5 × 10¹⁰ | 17 min per twist class |
| RTX 3070 Laptop, CUDA C++ | 2.7–3.0 × 10¹⁰ | 19–21 min |
| Mac mini, 14 CPU threads (Rust) | 8.8 × 10⁹ (2³² subcode) | about 70 min |
| Desktop, 12 CPU threads (Rust), while a Lean build used all cores | 1.1 × 10⁹ | not a fair measurement |

- Results are bit-identical across all backends tested.
- The GPU advantage here is a factor of about 3–4 over a 14-thread Apple CPU. It is not orders of magnitude, because XOR+popcount is cheap on modern CPUs too.
- What the GPU cannot do is the branchy part of symbolic work: proof search, kernel type checking, term rewriting with pattern matching. Those stay on CPUs.

## Exact results at n = 1

**Ranks.** 1⁴: 0; 2·1³: 3; 2²·1²: 9; 2³·1: 21; 2⁴: 45.

**The defect.**
- On 1⁴ every plaquette is a commutator of commuting elements, so Z(1⁴) = 1 and Δ_β(1) = −log Z(2⁴).
- The mean frustration on 2⁴ at β = 0 is 48 = P/2, so Δ_β(1) = 96β + O(β²).
- The exact seed threshold, Δ(1) = 1/7480, is at β = 1.3926034652 × 10⁻⁶.
- The elementary bound Δ ≤ 192βn⁴ gives β = 6.9630124777 × 10⁻⁷: a ratio of **exactly 2**.
- At n = 1 the bound is tight up to the factor 2, because the abelian side-1 torus is degenerate. This says nothing about n ≥ 2, where 4⁴ is out of reach of enumeration.
- As β → ∞, Δ(1) → 45 log 2 = 31.19162: only the flat configuration survives.

**Purities.** The telescope 8 log P₁ + 4 log P₂ + 2 log P₃ + log P₄ = −Δ holds to machine precision at every β tested, with P_k ≤ 1 and P₄ ≥ e^{−Δ}.

**Twist free energies on 2⁴**, F = −log(Z_twist/Z), for a single plane:

| β | 0.01 | 0.1 | 0.3 | 0.4 | 0.5 | 1.0 | 2.0 |
|---|---|---|---|---|---|---|---|
| F | 8.0 × 10⁻⁸ | 7.9 × 10⁻⁴ | 0.078 | 0.565 | 2.08 | 6.61 | 14.61 |

The rise is steepest between 0.4 and 0.5, near the self-dual coupling β_c = ½ log(1 + √2) ≈ 0.4407 of four-dimensional Z₂ gauge theory. On a box this small it is a crossover.
