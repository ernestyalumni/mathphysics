# YM-T01 — Exact sphere averages and the one-loop coefficient (Rust)

Tier T2 · Workstream B · Depends on: **PDD intent approved by Ernest** · Crate: `Rust/qft-symbolics` · Review: T2 different family

## Behaviour wanted
- Exact average of any monomial n₀^{a₀}…n_{d−1}^{a_{d−1}} over the unit sphere S^{d−1}: zero if any exponent is odd, else the rational Π(a_i − 1)!! / (d (d+2) … (d + |a| − 2)).
- Area of S^{d−1} as rational × power of π. Radial shell integral ∫_{r₁<|k|≤r₂} |k|^{−s} dᵈk in closed form including the logarithmic case s = d.
- A function that, given a rational function of k written as Σ c_j (monomial)/|k|^{2m_j} homogeneous of degree −d, returns its shell integral as coefficient × log(r₂/r₁).

## Fixtures (independent oracle: `PLAN/ledger/evidence/check_oneloop.py`, SymPy)
- d = 4: ⟨n₁²⟩ = 1/4, ⟨n₁² n₀²⟩ = 1/24, area 2π².
- h(k) = 4C (k₁²/|k|⁶ − 4k₁²k₀²/|k|⁸ − 1/|k|⁴): angular average −(11/3) C |k|⁻⁴; (2π)⁻⁴ ∫_{1/L<|k|≤1} h = −(11 C / 24π²) log L.
- Negative controls: a deliberately wrong double factorial must fail; an odd exponent must return exactly zero; s ≠ d must not produce a logarithm.

## Process
Read `MP/PDD.md` and `Rust/qft-symbolics/PDD.md`. Prompt first, then code, then tests. `CARGO_TARGET_DIR=/tmp/mathphysics-rust-target`, `--offline --locked`. No new dependencies.

## Acceptance
`cargo test` passes for the new tests; the κ_G fixture is asserted exactly (rational coefficient of log L / π²), not in floating point.
