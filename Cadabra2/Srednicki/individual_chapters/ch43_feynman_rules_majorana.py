"""
ch43_feynman_rules_majorana.py
===============================
Srednicki QFT — Chapter 43: The Feynman Rules for Majorana Fields

What this file covers:
  §43.A  Majorana condition ψ = ψ^c = C ψ̄^T
  §43.B  Majorana propagator (same as Dirac with ½ factor)
  §43.C  Majorana Feynman rules
  §43.D  Four-fermion contact interactions
  §43.E  Connection to Weyl theory

Run with:
    python3 ch43_feynman_rules_majorana.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §43.A  Majorana condition
# =============================================================================
sec("§43.A — Majorana condition ψ = ψ^c = C ψ̄^T")

print("Majorana field: ψ = ψ^c = C ψ̄^T")
print("  where C = iγ²γ⁰ (charge-conjugation matrix)")
print()
print("Consequences:")
print("  - Particle = antiparticle (neutrino, neutralino, gravitino)")
print("  - Only 2 degrees of freedom (not 4 like Dirac)")
print("  - Charge-neutral, so no electromagnetic interactions (unless Weyl + heavy)")
print()
print("Majorana vs. Dirac:")
print("  Dirac:  ψ = ψ_L + ψ_R  (4 dof, charged)")
print("  Majorana: ψ = ψ_L + ψ_L^c  (2 dof, neutral)")

# =============================================================================
# §43.B  Majorana propagator
# =============================================================================
sec("§43.B — Majorana propagator")

print("Propagator is same as Dirac propagator:")
print("  S_F(k) = i(k̸ + m)/(k² - m² + iε)")
print()
print("BUT: Wick contraction gives ½ factor when two ends are same field:")
print("  ⟨0| T{ψ(x) ψ̄(y)} |0⟩ = ½ S_F(x-y)   (Majorana)")
print("  ⟨0| T{ψ(x) ψ̄(y)} |0⟩ =     S_F(x-y)   (Dirac)")
print()
print("This ½ factor is crucial for loop calculations with Majorana fermions.")

# =============================================================================
# §43.C  Majorana Feynman rules
# =============================================================================
sec("§43.C — Majorana Feynman rules")

print("Feynman rules for Majorana theory (ψ⁴ theory or SUSY):")
print()
rules = [
    ("Majorana propagator",   "i(p̸ + m)/(p² - m² + iε)   [½ factor in contractions]"),
    ("Majorana vertex",       "-iy γ_μ (for Yukawa coupling)"),
    ("4-Majorana contact",    "-iy (for ψ⁴ theory)"),
    ("Closed Majorana loop",  "½ Tr[S_F]   (extra ½ from Wick pairing)"),
]
for name, rule in rules:
    print(f"  {name:35s} : {rule}")
print()
print("Ghost fields: Majorana ghosts appear in nonabelian gauge theories (Ch.80).")

# =============================================================================
# §43.D  Four-fermion contact interactions
# =============================================================================
sec("§43.D — Four-fermion contact interactions")

print("Effective 4-fermion operator (Fermi theory, or effective SUSY):")
print("  L_4f = -(G_F/√2) [ψ̄ γ^μ (1-γ⁵) ψ] [ψ̄ γ_μ (1-γ⁵) ψ]")
print()
print("Feynman diagram: 4-external-leg contact term (no propagator)")
print("  Amplitude: M ∝ -i G_F [ū γ^μ (1-γ⁵) u] [ū γ_μ (1-γ⁵) u]")
print()
print("In spinor-helicity form:")
print("  Contact terms have no propagator singularities → no BCFW poles")

# =============================================================================
# §43.E  Connection to Weyl theory
# =============================================================================
sec("§43.E — Majorana ↔ Weyl connection")

print("Majorana in 2-component Weyl notation:")
print("  ψ_α = ψ̄^α̇  (Majorana: left-handed only)")
print()
print("Mass term (Majorana mass):")
print("  L_mass = -½ m (ψᵀ C ψ + ψ̄ C γ⁵ ψ)    [SUSY-like]")
print()
print("In MSSM: Majoranaino mass terms = gaugino Majorana masses")
print("  L ⊃ ½ M_1 λ̃ λ̃ + ½ M_2 Ẽ Ẽ + h.c.")
print("  These are spinor-helicity vertices with |λ⟩|λ⟩ form")

print(f"\n{SEP}")
print("  ch43 — Feynman Rules for Majorana Fields — COMPLETE")
print(f"{SEP}")
