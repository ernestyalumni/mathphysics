"""
01_weyl_spinors.py
==================
Introduction to 2-component (Weyl) spinors in cadabra2.

Reference: Srednicki QFT, Chapters 34-35.
           van der Waerden notation: undotted indices α,β for (1/2,0) rep,
           dotted indices α̇,β̇ for (0,1/2) rep of SL(2,C).

Run with:
    python3 01_weyl_spinors.py
or inside Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/01_weyl_spinors.py
"""

import cadabra2
from cadabra2 import Ex, ExNode, __cdbkernel__

print("=" * 60)
print("Weyl spinors in cadabra2 (Srednicki Ch. 34-35)")
print("=" * 60)

# -------------------------------------------------------------------------
# 1. Declare spinor indices (undotted: SL(2,C) fundamental rep)
#    α, β, γ, δ  — undotted lower indices  (spinor indices)
#    In cadabra2, we declare indices with their range/properties.
# -------------------------------------------------------------------------
__cdbkernel__ = cadabra2.create_scope()

# Undotted spinor indices α,β,γ,δ
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))

# Dotted spinor indices \dot{\alpha} etc. — represented as \dal, \dbe, etc.
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))

print("\n[1] Declared spinor indices:")
print("    Undotted: alpha, beta, gamma, delta  (SL(2,C) fundamental)")
print("    Dotted:   dal, dbe, dga, dde          (SL(2,C) conjugate)")

# -------------------------------------------------------------------------
# 2. Epsilon tensors for raising/lowering spinor indices.
#    ε_{αβ} and ε^{αβ} with ε_{12} = 1, ε^{12} = -1 (Srednicki convention)
#    ε_{αβ} ε^{βγ} = δ_α^γ
# -------------------------------------------------------------------------
# Declare epsilon as antisymmetric
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))

print("\n[2] Epsilon tensors (SL(2,C) invariant):")
print("    eps_{alpha beta}: lowers undotted indices")
print("    eps^{alpha beta}: raises undotted indices")
print("    Convention: eps_{12} = 1, eps_{21} = -1  (Srednicki)")

# -------------------------------------------------------------------------
# 3. Define Weyl spinors
#    ψ_α : left-handed Weyl spinor (undotted lower index)
#    χ^α : left-handed with raised index
#    Raising: χ^α = ε^{αβ} χ_β
# -------------------------------------------------------------------------
# These are just symbols/expressions in cadabra2
psi_lower = Ex(r"\psi_{\alpha}")
chi_lower = Ex(r"\chi_{\beta}")
chi_upper = Ex(r"\chi^{\alpha}")

print("\n[3] Weyl spinor fields:")
print("    psi_{alpha} =", psi_lower)
print("    chi_{beta}  =", chi_lower)
print("    chi^{alpha} = eps^{alpha beta} chi_{beta}  (raised index)")

# -------------------------------------------------------------------------
# 4. Spinor inner products (Lorentz-invariant contractions)
#    Angle bracket:  <psi chi> = ε^{αβ} ψ_α χ_β  =  ψ^α χ_α
#    Square bracket: [psi chi] = ε_{α̇β̇} ψ̄^{α̇} χ̄^{β̇}
#
#    The angle bracket is antisymmetric: <psi chi> = -<chi psi>
# -------------------------------------------------------------------------
# Represent as expressions
angle_product = Ex(r"\epsilon^{\alpha\beta} \psi_{\alpha} \chi_{\beta}")
print("\n[4] Lorentz-invariant spinor products:")
print("    <psi chi> = eps^{alpha beta} psi_{alpha} chi_{beta}")
print("             =", angle_product)

# Show antisymmetry symbolically: <psi chi> = -<chi psi>
# This follows from antisymmetry of epsilon and Grassmann nature of spinors.
# eps^{ab} psi_a chi_b = eps^{ab} psi_a chi_b
# swap a <-> b: eps^{ba} chi_a psi_b = -eps^{ab} chi_a psi_b
# For Grassmann spinors: -eps^{ab} chi_a psi_b = +eps^{ab} psi_b chi_a ... wait
# More carefully: spinors anticommute (Grassmann) so psi chi = -chi psi
# Hence <psi chi> = -<chi psi>

print("\n[5] Antisymmetry of angle bracket:")
print("    <psi chi> = eps^{ab} psi_a chi_b")
print("    <chi psi> = eps^{ab} chi_a psi_b")
print("    Because eps is antisymmetric and spinors are Grassmann:")
print("    <psi chi> = -<chi psi>  [antisymmetric under swap]")

# -------------------------------------------------------------------------
# 5. Dotted spinors (right-handed, complex conjugate representation)
#    psibar^{alpha-dot} = complex conjugate of psi_alpha
#    psibar_{alpha-dot} = eps_{alpha-dot beta-dot} psibar^{beta-dot}
# -------------------------------------------------------------------------
psibar_upper = Ex(r"\bar{\psi}^{\dal}")
psibar_lower = Ex(r"\bar{\psi}_{\dal}")

print("\n[6] Dotted (right-handed) spinors:")
print("    psibar^{dotalpha} — conjugate rep (0,1/2)")
print("    psibar_{dotalpha} = eps_{dotalpha dotbeta} psibar^{dotbeta}")

square_product = Ex(r"\epsilon_{\dal\dbe} \bar{\psi}^{\dal} \bar{\chi}^{\dbe}")
print("    [psi chi] = eps_{dala dalb} psibar^{dala} chibar^{dalb}")
print("             =", square_product)

# -------------------------------------------------------------------------
# 6. Summary of conventions (Srednicki Ch. 34)
# -------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Summary (Srednicki Ch. 34 conventions)")
print("=" * 60)
print("""
  Undotted index α:   left-handed (1/2, 0) Weyl spinor
  Dotted index α̇:    right-handed (0, 1/2) Weyl spinor

  Raising:    ψ^α = ε^{αβ} ψ_β
  Lowering:   ψ_α = ε_{αβ} ψ^β

  Epsilon normalization (Srednicki):
    ε_{12} = ε^{21} = +1
    ε_{21} = ε^{12} = -1
    ε^{αγ} ε_{γβ} = δ^α_β

  Angle bracket (left-handed product):
    <ψ χ> = ε^{αβ} ψ_α χ_β = ψ^α χ_α = -ψ_α χ^α

  Square bracket (right-handed product):
    [ψ̄ χ̄] = ε_{α̇β̇} ψ̄^{α̇} χ̄^{β̇}

  Antisymmetry:
    <ψ χ> = -<χ ψ>    [Grassmann + antisymmetric ε]
    [ψ̄ χ̄] = -[χ̄ ψ̄]
""")

print("Done: 01_weyl_spinors.py")
