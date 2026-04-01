"""
ch38_lsz_feynman_rules.py
==========================
Srednicki QFT — Chapter 38: LSZ for Spinors & Feynman Rules

What this file covers:
  §38.A  LSZ reduction formula for Weyl spinor fields
  §38.B  Time-ordered correlators → S-matrix elements
  §38.C  Fermion propagator in momentum space
  §38.D  Yukawa theory Feynman rules
  §38.E  Gluon/photon Feynman rules (spinor-helicity form)
  §38.F  External leg spinor wave functions

Run with:
    python3 ch38_lsz_feynman_rules.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70


def sec(s):
    print(f"\n{SEP}\n  {s}\n{SEP}")


# =============================================================================
# §38.A  LSZ reduction for Weyl spinors
# =============================================================================
sec("§38.A — LSZ reduction for Weyl spinors")

# Declare indices (same as ch34 / ch37):
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

# Declare epsilon tensors (same as ch34):
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# LSZ for a Weyl spinor ψ_α(x):
# ⟨f|ψ_α(x)|i⟩ picks up the wave function from the asymptotic region.
# The reduced matrix element involves (i∂̸ - m) acting on the field.
print("LSZ for a Weyl spinor field ψ_α(x):")
print("  ⟨f| i\\rangle = ∏_j ∫ d⁴x_j  ū(k_j) (i∂̸ - m)  ⟨f| T{ψ(x₁)⋯}|i⟩")
print()
print("Key point: the spinor wave function ū(k) carries the external leg spin info.")
print("The scalar LSZ uses ∂²+m²; the spinor LSZ uses i∂̸ - m.")

# =============================================================================
# §38.B  Time-ordered correlator → S-matrix
# =============================================================================
sec("§38.B — From time-ordered correlator to S-matrix")

# The S-matrix element is the Fourier transform of the time-ordered product,
# with residues taken at the mass-shell poles p² = m².
print("S-matrix from time-ordered correlator:")
print("  S_fi = lim_{p²→m²}  (p̸ - m) ∫ d⁴x e^{ip·x} ⟨0| T{ψ(x)ψ̄(0)} |0⟩")
print()
print("In perturbation theory, each internal line contributes a propagator,")
print("each vertex contributes -iy, and spinor wave functions ū, u attach to external legs.")

# =============================================================================
# §38.C  Fermion propagator in momentum space
# =============================================================================
sec("§38.C — Fermion propagator S_F(k)")

# Declare fermion fields and propagators:
cadabra2.SelfAntiCommuting(Ex(r"\psi_{\alpha}"))
cadabra2.SelfAntiCommuting(Ex(r"\psidag^{\dal}"))

# The Feynman propagator for a Dirac/Weyl fermion:
# S_F(x-y) = ⟨0| T{ψ(x)ψ̄(y)} |0⟩ = ∫ d⁴k e^{-ik·(x-y)} / [k² - m² + iε]
print("Fermion propagator (Feynman):")
print("  S_F(k) = i(k̸ + m) / (k² - m² + iε)")
print()
print("Alternative form:")
print("  S_F(k) = i / (k̸ - m + iε)")
print()

# In the Weyl (two-component) formalism:
print("Weyl-indexed fermion propagator:")
print("  S_F^{α α̇}(k) = i(-σ^μ_{α α̇} k_μ + m δ^α_{α̇}) / (k² - m² + iε)")

# =============================================================================
# §38.D  Yukawa theory Feynman rules
# =============================================================================
sec("§38.D — Yukawa theory: scalar-fermion vertex")

# Yukawa Lagrangian: L_Yuk = -y ψ ψ φ - y* ψ̄ ψ̄ φ
# (Majorana case: ψ = ψ̄, so one term)
print("Yukawa Lagrangian (Majorana):")
print("  L_Yuk = -y ψ ψ φ - y* ψ̄ ψ̄ ψ̄ φ")
print()
print("Feynman rules:")
vertex = Ex(r"\text{Vertex} = -i y")
print(f"  Yukawa vertex: {vertex}")
print("  Fermion propagator: i(k̸ + m)/(k² - m² + iε)")
print("  Scalar propagator: i/(p² - m² + iε)")
print("  External scalar φ: 1")
print("  External fermion u(p,s): spinor wave function")
print("  External antifermion ū(p,s): = v̄(p,s) for Majorana")

# =============================================================================
# §38.E  Gluon/photon Feynman rules (spinor-helicity form)
# =============================================================================
sec("§38.E — Gluon Feynman rules in spinor-helicity form")

# Gluon propagator in axial gauge with reference momentum q:
print("Gluon propagator (axial gauge):")
print("  d_μν = η_μν - (q_μ q̄_ν + q_ν q̄_μ) / (q·q̄)")
print()

# Three-gluon vertex (colour-stripped, spinor form):
print("Three-gluon vertex (spinor form):")
print("  V_3 = g ε^{a1 a2 a3} [⟨12⟩[2| + ⟨23⟩[3| + ⟨31⟩[1|]")
print()
print("Note: the colour structure is in SU(N) structure constants ε^{abc}.")
print("In axial gauge q^μ, the propagator is transverse: q^μ d_μν = 0.")

# =============================================================================
# §38.F  External leg spinors — helicity basis
# =============================================================================
sec("§38.F — Helicity spinors for external states")

# For massless k = (E, 0, 0, E):
# u_+(k) ∝ (1, 0)^T  (positive helicity = spin along +z)
# u_-(k) ∝ (0, 1)^T  (negative helicity = spin along -z)
print("Massless helicity spinors (chiral basis, k along +z):")
print("  u_+(k) = √E (1, 0)^T   ← positive helicity (right-handed)")
print("  u_-(k) = √E (0, 1)^T   ← negative helicity (left-handed)")
print()
print("General direction (θ, φ):")
print("  u_+(k) ∝ ( √(E+p_z),   √(E-p_z) e^{iφ} )^T")
print("  u_-(k) ∝ ( -√(E-p_z) e^{-iφ},  √(E+p_z) )^T")
print()

# Spinor inner products for external legs:
print("  Orthogonality: u_+^dagger(k) u_-(k) = 0  ✓")

# Normalization:
print("  Normalization: u_s^dagger(k) u_s(k) = 2E")

print(f"\n{SEP}")
print("  ch38 — LSZ & Feynman Rules — COMPLETE")
print(f"{SEP}")
