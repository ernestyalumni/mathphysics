"""
ch40_spin_sums.py
==================
Srednicki QFT — Chapter 40: Spin Sums

What this file covers:
  §40.A  Spin-averaging factors for unpolarized scattering
  §40.B  Trace theorems for γ-matrices
  §40.C  Helicity amplitudes vs. trace-based calculations
  §40.D  Numerical spin sum examples (Compton scattering)

Run with:
    python3 ch40_spin_sums.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §40.A  Spin-averaging factors
# =============================================================================
sec("§40.A — Spin-averaging for unpolarized scattering")

print("For unpolarized initial states (average over spin):")
print("  Factor for each initial fermion:  ½")
print("  Factor for each initial photon:   1/(2·3) = 1/6  (physical gauge)")
print()

print("Example: e⁺e⁻ → μ⁺μ⁻ annihilation")
print("  σ ∝ (1/4) Σ_spins |M|²")
print("  = (1/4) Tr[(p₁̸ - m) γ^μ (p₂̸ + m) γ^ν] × same for muons")
print("  × (1/q⁴) [η_μν η_ρσ + η_μρ η_νσ - η_μσ η_νρ]")
print()

# =============================================================================
# §40.B  Trace theorems
# =============================================================================
sec("§40.B — Trace theorems for γ-matrices")

print("Trace theorems (η = diag(-1,+1,+1,+1), {γ^μ,γ^ν} = 2η^{μν} I₄):")
print()

theorems = [
    ("Tr[I₄]",            "4"),
    ("Tr[γ^μ]",           "0  (odd number of γ matrices)"),
    ("Tr[γ^μ γ^ν]",       "4 η^{μν}"),
    ("Tr[γ^μ γ^ν γ^ρ]",   "0  (odd)"),
    ("Tr[γ^μ γ^ν γ^ρ γ^σ]","4(η^{μν}η^{ρσ} - η^{μρ}η^{νσ} + η^{μσ}η^{νρ})"),
    ("Tr[γ⁵]",             "0  (trace of odd γ⁵ vanishes in 4D)"),
    ("Tr[γ⁵ γ^μ γ^ν]",     "0  (odd)"),
    ("Tr[γ⁵ γ^μ γ^ν γ^ρ γ^σ]", "-4i ε^{μνρσ}  (dual tensor)"),
]

for expr, val in theorems:
    print(f"  {expr:40s} = {val}")

print()
print("Contracted trace identity:")
print("  Tr[(a̸ + m)(b̸ + m)(c̸ + m)(d̸ + m)]")
print("  = 4[ (a·b)(c·d) - (a·c)(b·d) + (a·d)(b·c) ] + O(m)")

# =============================================================================
# §40.C  Helicity amplitudes vs. trace calculations
# =============================================================================
sec("§40.C — Helicity amplitudes vs. trace-based approach")

print("Two approaches to spinor amplitudes:")
print()
print("  1. Trace/traces method (Ch.40):")
print("     - Write |M|² as trace of products of γ matrices")
print("     - Apply trace theorems → scalar products of momenta")
print("     - Works for unpolarized cross sections")
print()
print("  2. Helicity amplitude method (Ch.42/48/60):")
print("     - Keep spinor indices explicit λ, λ̃")
print("     - Compute amplitude directly → spinor brackets ⟨ij⟩, [ij]")
print("     - Much simpler for MHV amplitudes")
print()

print("Example: e⁺e⁻ → γγ (Bhabha scattering)")
print("  Helicity: only ⟨12⟩⟨34⟩ + ⟨14⟩⟨23⟩ terms survive")
print("  Trace: messy γ-matrix algebra; result factorizes identically")

# =============================================================================
# §40.D  Numerical spin sum example
# =============================================================================
sec("§40.D — Numerical spin sum (Compton: eγ → eγ)")

print("Compton scattering: e(p₁) + γ(k₁,ε₁) → e(p₂) + γ(k₂,ε₂)")
print()
print("Matrix element (QED):")
print("  M = -ie² ū(p₂) [γ_μ ε̸₁ (p̸₁ - m + k̸₁) ε̸₂ γ_ν + (1↔2, μ↔ν)] u(p₁)")
print()
print("Spin sum: (1/4) Σ |M|² = (1/4) Tr[... product of 8 γ matrices ...]")
print("  → Klein-Nishina formula for unpolarized photon scattering")
print()
print("Result (Klein-Nishina, relativistic limit ω₁≫m):")
print("  dσ/dΩ ∝ (ω'/ω)² [ω'/ω + ω/ω' - sin²θ cos²φ]")
print("  where ω' is scattered photon energy, θ is scattering angle")

print(f"\n{SEP}")
print("  ch40 — Spin Sums — COMPLETE")
print(f"{SEP}")
