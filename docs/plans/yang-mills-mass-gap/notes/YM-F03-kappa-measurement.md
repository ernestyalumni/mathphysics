# YM-F03: measuring the short-side cost κ(n) in SU(2)

Label: pre-registration, written before any data of this run. Author: Claude Opus 5.5, 2026-11-25.

**Quantity.** κ_β(n) = 2 log Z(n,2n,2n,2n) − log Z(2n,2n,2n,2n), the seed of the short-side recurrence (`formal/YMFormal/ShortSide.lean`, seed κ ≤ 1/80). With A = n×(2n)³ and B = (2n)⁴, P(B) = 2P(A), so

  κ_β(n) = P(B) ∫₀^β [⟨a₀⟩_A − ⟨a₀⟩_B] dβ′.

**Runs.** SU(2) heat bath with two overrelaxation sweeps per update, Simpson on h = 0.1, independent replicas.
- n = 2 (2×4³, 4⁴): β ∈ [0.4, 1.6], 3×10⁵ updates per point and replica, laptop and Mac CPUs.
- n = 4 (4×8³, 8⁴): β ∈ [1.6, 2.4], replicas of 3×10⁵ updates per point, laptop and desktop GPUs, for about 4 hours.

The integrals start at β₀ = 0.4 and 1.6. The truncated strong-coupling value of κ there is added nowhere; it is reported as the estimated offset.

**Predictions.** The truncated strong-coupling formula is κ ≈ 2(12 n² u^{2n²} + 24 n³ u^{4n}), u = I₂(β)/I₁(β). At n = 2 it is 480u⁸.
1. At n = 2, β = 1.0, the measured κ agrees with 480u⁸ within 2σ plus 30% of the value (higher orders).
2. At n = 2, κ(1.0) is below Δ(1.0) = 1.456(31) (YM-N03) by a factor of at least 20.
3. At n = 4, κ(2.0) < 0.1.

Expected error: about 0.004 at n = 2. At n = 4 about 0.03, so the crossing of 1/80 cannot be resolved there; only where κ rises through about 0.1 can.

**What it cannot show.** Monte Carlo certifies nothing. These runs say where the seed might be sought, not that it holds.
