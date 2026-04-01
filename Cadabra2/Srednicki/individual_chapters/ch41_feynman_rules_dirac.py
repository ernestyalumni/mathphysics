"""
ch41_feynman_rules_dirac.py
===========================
Srednicki QFT — Chapter 41: The Feynman Rules for Dirac Fields

What this file covers:
  §39.A  Dirac propagator S_F(x-y) — position and momentum space
  §39.B  External fermion wave functions u(p,s), v(p,s)
  §39.C  Dirac trace and spin sum theorems
  §39.D  QED vertex factor and photon propagator
  §39.E  Feynman rules summary for spinor QED

Run with:
    python3 ch41_feynman_rules_dirac.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §39.A  Dirac propagator
# =============================================================================
sec("§39.A — Dirac propagator S_F(x-y)")

cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

print("Position-space Feynman propagator for a Dirac field:")
print("  S_F(x-y) = ⟨0| T{ψ(x)ψ̄(y)} |0⟩")
print("  = ∫ d⁴k e^{-ik·(x-y)} / [(2π)⁴ (k² - m² + iε)]")
print()

print("Momentum-space propagator:")
print("  S_F(k) = i(k̸ + m) / (k² - m² + iε)")
print("          = i / (k̸ - m + iε)")
print()

print("Index structure (Weyl components):")
print("  S_F^{α α̇}(k) = i(-σ^μ_{α α̇} k_μ + m δ^α_{α̇}) / (k² - m² + iε)")

# =============================================================================
# §39.B  External fermion wave functions
# =============================================================================
sec("§39.B — External fermion wave functions u(p,s), v(p,s)")

print("Incoming fermion u(p,s):")
print("  u_s(p) = √(E+m) [ χ_s , (σ·p/(E+m)) χ_s ]^T")
print()
print("Outgoing antifermion v(p,s):")
print("  v_s(p) = √(E+m) [ (σ·p/(E+m)) χ_s , χ_s ]^T")
print("  where χ_+ = (1,0)^T, χ_- = (0,1)^T (spin along z)")
print()

print("Massless limit (E = |p|, m→0):")
print("  u_+(p) ∝ √(2E) [ 1, 0 ]^T   (spin along +z)")
print("  u_-(p) ∝ √(2E) [ 0, 1 ]^T   (spin along -z)")
print("  v_+(p) ∝ √(2E) [ 0, 1 ]^T   (antifermion, opposite)")
print("  v_-(p) ∝ √(2E) [ 1, 0 ]^T   (antifermion, opposite)")
print()

print("Completeness relations:")
print("  Σ_s u_s(p) ū_s(p) = k̸ + m")
print("  Σ_s v_s(p) v̄_s(p) = k̸ - m")

# =============================================================================
# §39.C  Dirac trace and spin sum theorems
# =============================================================================
sec("§39.C — Trace theorems and spin sums")

print("Spin sum (averaging over initial spins, summing over final):")
print("  Σ_s u_s(p) ū_s(p) = (p̸ + m)")
print("  Σ_s v_s(p) v̄_s(p) = (p̸ - m)")
print()

print("Trace theorems (handles closed fermion loops):")
print("  Tr[1] = 4")
print("  Tr[p̸] = 4(p·q)")
print("  Tr[p̸ q̸] = 4[(p·q)² - p² q²]")
print("  Tr[p̸ q̸ r̸] = 4[(p·q)(q·r)(r·p) - (p·q)² r² - ...]")
print()

print("Useful: Tr[γ⁵] = 0, Tr[γ⁵ γ^μ γ^ν] = 0, Tr[γ⁵ γ^μ γ^ν γ^ρ γ^σ] = -4i ε^{μνρσ}")

# =============================================================================
# §39.D  QED vertex factor and photon propagator
# =============================================================================
sec("§39.D — QED vertex factor and photon propagator")

print("QED Lagrangian (Dirac + photon):")
print("  L = ψ̄(iγ^μ D_μ - m)ψ - ¼ F_μν F^{μν}")
print("  D_μ = ∂_μ + ieA_μ")
print()

print("Feynman rules for spinor QED:")
print("  Photon propagator: -i η_μν / q²   (Feynman gauge)")
print("  Fermion propagator: i(p̸ + m)/(p² - m² + iε)")
print("  QED vertex: -ie γ_μ")
print()

print("Compton scattering e⁺e⁻ → γγ: two diagrams (s-channel, u-channel)")
print("  Amplitude ∝ ū(p₂)[-ieγ_μ](i(p₁̸ + m)/(p₁-k)²-m²)(-ieγ_ν)u(p₁) ε^μ(k) ε^ν(k')")

# =============================================================================
# §39.E  Feynman rules summary for spinor QED
# =============================================================================
sec("§39.E — Complete Feynman rules for spinor QED")

rules = [
    ("Photon propagator (Feynman)", "-i η_μν / q²"),
    ("Fermion propagator", "i(p̸ + m) / (p² - m² + iε)"),
    ("QED vertex", "-ie γ_μ"),
    ("External fermion", "u(p,s) or ū(p,s)"),
    ("External antifermion", "v(p,s) or v̄(p,s)"),
    ("External photon", "ε_μ(k,λ)"),
    ("Closed fermion loop", "Tr[ ... ] × spin sum"),
    ("Prop. vertex correction", "ie³ γ_μ Γ^μ(p,q)"),
]
for name, rule in rules:
    print(f"  {name:35s} : {rule}")

print(f"\n{SEP}")
print("  ch41 — Feynman Rules for Dirac Fields — COMPLETE")
print(f"{SEP}")
