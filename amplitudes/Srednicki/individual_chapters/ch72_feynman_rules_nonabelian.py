"""
ch72_feynman_rules_nonabelian.py
================================
Srednicki QFT — Chapter 72: The Feynman Rules for Nonabelian Gauge Theory

What this file covers:
  §72.A  Gluon propagator in covariant gauges
  §72.B  Ghost propagator and ghost-gluon vertex
  §72.C  Three-gluon vertex
  §72.D  Four-gluon vertex
  §72.E  Quark-gluon vertex and color factors
  §72.F  Color-ordered Feynman rules
  §72.G  Sample: ggg and qqb→gg color-ordered subamplitudes

Run with:
    python3 ch72_feynman_rules_nonabelian.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §72.A  Gluon propagator in covariant gauges
# =============================================================================
sec("§72.A — Gluon propagator in covariant gauges")

print("Gluon propagator (Feynman gauge ξ = 1):")
print("  D^{ab}_μν(k) = -i δ^{ab} η_μν / k²")
print()
print("General covariant gauge (gauge parameter ξ):")
print("  D^{ab}_μν(k) = -i δ^{ab} [η_μν - (1-ξ) k_μ k_ν/k²] / k²")
print()
print("Physical observables are ξ-independent (Ward identities).")
print("Popular gauges:")
print("  ξ = 1 : Feynman gauge (simplest propagator)")
print("  ξ = 0 : Landau gauge (transverse propagator)")
print("  ξ → ∞ : axial gauge (no ghosts)")

# =============================================================================
# §72.B  Ghost propagator and ghost-gluon vertex
# =============================================================================
sec("§72.B — Ghost propagator and ghost-gluon vertex")

print("Fadeev-Popov ghosts c̄^a, c^b (complex Grassmann fields):")
print()
print("Ghost propagator:")
print("  ⟨0| T{c^a(x) c̄^b(y)} |0⟩ = i δ^{ab} / ∂²   (same as scalar)")
print()
print("Ghost-gluon vertex:")
print("  -g f^{abc} k̄_μ   (momentum k̄ flowing INTO ghost vertex)")
print()
print("Important: ghosts are unphysical — they cancel unphysical gluon polarizations.")
print("In axial gauge (n·A = 0): no ghost diagrams at all.")

# =============================================================================
# §72.C  Three-gluon vertex
# =============================================================================
sec("§72.C — Three-gluon vertex")

print("Three-gluon vertex (color-stripped tensor structure):")
print()
print("  V^{abc}_μνρ(k₁, k₂, k₃) =")
print("    g f^{abc} [ (k₁ - k₂)_ρ η_μν + (k₂ - k₃)_μ η_νρ + (k₃ - k₁)_ν η_ρμ ]")
print()
print("Color factor: g f^{abc} (structure constant coupling)")
print()
print("The vertex tensor satisfies:")
print("  k₁^μ V^{abc}_μνρ = 0   (Ward identity, ghost cancels)")
print("  Momentum conservation: k₁ + k₂ + k₃ = 0")

# =============================================================================
# §72.D  Four-gluon vertex
# =============================================================================
sec("§72.D — Four-gluon vertex")

print("Four-gluon vertex (two equivalent forms):")
print()
print("Form 1 (direct):")
print("  -i g² [ f^{abe} f^{cde} (η_μρ η_νσ - η_μσ η_νρ)")
print("           + f^{ace} f^{bde} (η_μν η_ρσ - η_μρ η_νσ)")
print("           + f^{ade} f^{bce} (η_μσ η_νρ - η_μν η_ρσ) ]")
print()
print("Form 2 (dual):")
print("  -ig² [ C_1^{abcd} η_μν η_ρσ + C_2^{abcd} η_μρ η_νσ + C_3^{abcd} η_μσ η_νρ ]")
print()
print("Both forms related by renaming (μ,ν,a) ↔ (ρ,σ,d) and Jacobi identity.")
print("No momentum factors in the 4-gluon vertex (pure contact interaction).")

# =============================================================================
# §72.E  Quark-gluon vertex and color factors
# =============================================================================
sec("§72.E — Quark-gluon vertex and color factors")

print("Quark-gluon vertex (color-stripped):")
print("  -ig γ_μ   (same as QED, but color matrix T^a acts on quark color)")
print()
print("Color factor for q → q g:")
print("  T^a_ij  (quark i → quark j, gluon color a)")
print()
print("Color factors:")
print("  C_F = (N_c²-1)/(2N_c) = 4/3 for N_c = 3")
print("  C_A = N_c = 3")
print("  T_R = ½   (fundamental trace normalization)")

# =============================================================================
# §72.F  Color-ordered Feynman rules
# =============================================================================
sec("§72.F — Color-ordered Feynman rules (subamplitudes)")

print("Color-ordered amplitudes A_n (gluons only):")
print("  M_n = Σ_{σ∈S_n/Z_n} Tr[T^{a_σ(1)}...T^{a_σ(n)}] A_n(σ(1),...,σ(n))")
print()
print("Primitive (single-trace) subamplitudes A_n obey:")
print("  - Gauge invariance: sum over permutations")
print("  - Color decomposition: multi-trace at subleading order in 1/N_c")
print()
print("Leading-color (Nc → ∞) basis: n!/(n-2)! = n(n-1) traces")
print("Subleading corrections: 1/N_c² suppressed")

# =============================================================================
# §72.G  Sample: ggg and qqb→gg subamplitudes
# =============================================================================
sec("§72.G — Sample: ggg and q̄q→gg color-ordered subamplitudes")

print("Three-gluon amplitude (MHV, all helicities):")
print("  A_3(1^+, 2^+, 3^+) = 0   (all-plus vanishes at tree)")
print("  A_3(1^-, 2^-, 3^+) = i ⟨12⟩³ / ⟨23⟩⟨31⟩   (Parke-Taylor)")
print()
print("Quark antiquark → two gluons:")
print("  A_4(q̄^-, q^+, g^+, g^+) = 0   (vanishing helicity configuration)")
print("  A_4(q̄^-, q^+, g^-, g^+) = i ⟨q 1⟩² ⟨q 2⟩² / (⟨q 1⟩⟨12⟩⟨2q⟩⟨q q⟩)")
print()
print("BCFW recursion builds all n-gluon from 3-particle seeds (MHV amplitudes).")

print(f"\n{SEP}")
print("  ch72 — Feynman Rules for Nonabelian Gauge Theory — COMPLETE")
print(f"{SEP}")
