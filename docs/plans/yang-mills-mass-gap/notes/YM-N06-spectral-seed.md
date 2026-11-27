# YM-N06: spectral lower bound on the short-side cost κ(n)

Label: pre-registration, written 2026-11-27 before any data of this run; unreviewed. Author: Claude Fable 5.2 (strategist).

## Why the direct measurement stalls

The estimator of YM-F03 is κ_β(n) = P(B) ∫₀^β [⟨a₀⟩_A − ⟨a₀⟩_B] dβ′ with A = n × (2n)³ and B = (2n)⁴. Its error is the plaquette error times the number of plaquettes, accumulated along the path in β. After six pairs from β = 0 and four pairs from β = 1.6 on an RTX 3060 the absolute error at n = 4 is 0.06 to 0.11. Twenty more GPU hours would reach about 0.065. The seed threshold is 1/80 = 0.0125. At n = 8 the same error costs 16 times the updates for an error 4 times larger.

## The quantity as a spectral sum

Let 𝖳 be the transfer operator in the short direction on the spatial torus (2n)³, with eigenvalues λ₀ > λ₁ ≥ … ≥ 0 (nonnegative and trace class on every finite torus), and

  s_t = Σ_{i≥1} (λ_i/λ₀)^t = Σ_{i≥1} e^{−t E_i}.

Then Z(n,2n,2n,2n) = λ₀ⁿ (1 + s_n) and Z(2n,2n,2n,2n) = λ₀²ⁿ (1 + s_{2n}), so

  κ(n) = 2 log Z(n,2n,2n,2n) − log Z(2n,2n,2n,2n) = 2 log(1 + s_n) − log(1 + s_{2n}).

Two lower bounds follow.

1. Σ x_i² ≤ (Σ x_i)² for nonnegative x_i gives s_{2n} ≤ s_n², hence κ(n) ≥ F(s_n) with F(s) = log[(1+s)²/(1+s²)]. F′(s) = 2(1−s)/[(1+s)(1+s²)], so F increases on [0,1] and F(1) = log 2.
2. With weights a_i = λ_iⁿ/Tr 𝖳ⁿ, the purity is Σ a_i² ≤ a₀ Σ a_i = a₀ = 1/(1 + s_n), hence κ(n) ≥ log(1 + s_n).

Together: if 0 ≤ s_low ≤ s_n, then κ(n) ≥ F(min(s_low, 1)). (If s_n ≤ 1 use bound 1 and monotonicity. If s_n > 1, bound 2 gives κ(n) > log 2 = F(1).) Every state put into s_low with an upper bound on its energy therefore gives a lower bound on κ(n). For small s, F(s) = 2s − 2s² + O(s³).

**Necessary condition for the seed.** κ(n) ≤ 1/80 forces F(s_n) ≤ 1/80, that is s_n ≤ 0.00629. A single level of energy m then needs n·m ≥ 5.07. The three unit-flux torelon states of energy E_T need n·E_T ≥ 6.17.

## Runs

SU(2) Wilson action, `examples/torelon_scan.rs` (one heat-bath and two over-relaxation sweeps per update, a measurement every 5 updates, 2000 thermalization updates), lattices (2n)³ × 16 with time as the long direction.

| n | Lattice | β | Replicas × measurements | Machine |
|---:|---|---|---|---|
| 4 | 8³ × 16 | 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30 | 2 × 800,000 | Mac CPU |
| 6 | 12³ × 16 | 2.25, 2.30, 2.35, 2.40 | 1 × 200,000 | desktop RTX 3060 |
| 8 | 16³ × 16 | 2.35, 2.40, 2.45, 2.50 | 1 × 110,000 | desktop RTX 3060 |

States included: the three unit-flux torelons (one per spatial direction, energy E_T(2n)) and the zero-momentum 0++ glueball (energy m₀(2n)),

  s_low = 3 e^{−n E_T} + e^{−n m₀}.

## Analysis, fixed now

`glueball_fit 16 20` on the `-tor` and `-plaq` streams. Primary energy for the bound: the effective mass −ln c(2)/c(1) of the variationally projected correlator plus two jackknife errors. An effective mass of a positive correlator lies above the lowest energy in its channel, so this underestimates s. Secondary: the fit windows [2,5] and [1,3] as printed. An upper *estimate*, not a bound, adds glueball momentum states with the free lattice dispersion and is labeled as such.

## Predictions

1. s_low reaches 0.00629 at n·a√σ between 1.7 and 2.1 for n = 4, 6 and 8, that is near β ≈ 2.20, 2.33 and 2.43 (a√σ from Teper, hep-th/9812187). The three crossings agree in units of n·a√σ within 15%.
2. At n = 4, β = 2.2 the lower bound is below the direct estimate 0.157(104) plus two errors, and at β = 2.3 it is below 2.13.
3. E_T(L) lies within 15% of σL − π/(3L) wherever L·a√σ ≥ 2.5.

If prediction 1 fails, with the crossings drifting in n·a√σ by more than 30% from n = 4 to n = 8, the statement that the seed can hold from a fixed physical size on loses its numerical support.

## Remark for YM-N02 (unreviewed)

The states with nonzero electric flux e ∈ ℤ₂³ contribute Σ_{e≠0} Tr_e 𝖳ⁿ / λ₀ⁿ to s_n, and Tr_e 𝖳ⁿ = (1/8) Σ_k (−1)^{k·e} Z_k, where Z_k is the partition function of the box n × (2n)³ with temporal center twist k. The seed therefore bounds the 't Hooft electric-flux free energies of that box from below. The twist ratios of YM-N02 measure exactly this part of s_n.

## What it cannot show

Monte Carlo certifies nothing. The bound is a necessary condition for the seed: it says below which physical size κ ≤ 1/80 fails, not where it holds. States left out (momentum states, 2++, sectors with two or three units of flux) only raise κ.
