"""
05_fierz_identities.py
======================
Fierz rearrangement identities for 2-component spinors.

Reference: Srednicki QFT, Appendix B; Chapter 36.
           Dreiner, Haber, Martin (2010) — "Two-component spinor techniques"
           arXiv:0812.1594

Fierz identities arise from the completeness of the sigma matrices
as a basis for 2x2 matrices. They are used extensively in:
  - Loop calculations (Srednicki Ch. 51, 62)
  - SUSY algebra
  - Spinor index gymnastics in gauge theories

Run with:
    python3 05_fierz_identities.py
or inside Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/05_fierz_identities.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__

print("=" * 60)
print("Fierz identities for 2-component spinors (Srednicki App. B)")
print("=" * 60)

__cdbkernel__ = cadabra2.create_scope()

# -------------------------------------------------------------------------
# 1. Completeness relation for sigma matrices
#    The 2x2 matrices {I, σ^i} form a complete basis.
#    In spinor index notation (Srednicki notation):
#
#    (σ^μ)_{αα̇} (σ̄_μ)^{β̇β} = -2 δ_α^β δ^{β̇}_{α̇}
#
#    This is the fundamental Fierz identity from which all others follow.
# -------------------------------------------------------------------------
print("\n[1] Fundamental Fierz identity (sigma completeness):")
print("""
    sigma^mu_{alpha dotalpha} sigmabar_{mu}^{dotbeta beta} = -2 delta_alpha^beta delta^{dotbeta}_{dotalpha}

    Proof: The sigma matrices {I, sigma^i} form a complete orthonormal basis
    for 2x2 matrices under the trace inner product Tr[A B] = 2 delta_{AB}.
    The completeness relation in matrix form:

    sum_mu sigma^mu_{aa'} sigmabar_mu^{b'b} = -2 delta_a^b delta^{b'}_{a'}

    (factor -2 from normalization convention Tr[sigma^mu sigmabar^nu] = -2 eta^{mu nu})
""")

# -------------------------------------------------------------------------
# 2. Clifford algebra for sigma matrices
#    (σ^μ σ̄^ν + σ^ν σ̄^μ)_α^β = -2 η^{μν} δ_α^β
#    (σ̄^μ σ^ν + σ̄^ν σ^μ)^{α̇}_{β̇} = -2 η^{μν} δ^{α̇}_{β̇}
# -------------------------------------------------------------------------
print("[2] Clifford algebra (2-component form):")
print("""
    (sigma^mu sigmabar^nu + sigma^nu sigmabar^mu)_alpha^beta = -2 eta^{mu nu} delta_alpha^beta
    (sigmabar^mu sigma^nu + sigmabar^nu sigma^mu)^{dala}_{dalb} = -2 eta^{mu nu} delta^{dala}_{dalb}

    Antisymmetric combinations define Lorentz generators:
    sigma^{mu nu}_alpha^beta = (i/4)(sigma^mu sigmabar^nu - sigma^nu sigmabar^mu)_alpha^beta
    sigmabar^{mu nu}^{dala}_{dalb} = (i/4)(sigmabar^mu sigma^nu - sigmabar^nu sigma^mu)^{dala}_{dalb}

    These satisfy [sigma^{mu nu}, sigma^{rho sigma}] = Lorentz algebra
""")

# -------------------------------------------------------------------------
# 3. Standard Fierz rearrangement for spinor bilinears
#    For any four 2-component Weyl spinors χ, ψ, φ, ξ:
#
#    (χ^α ψ_α)(φ^β ξ_β) = -(χ^α φ_α)(ψ^β ξ_β) - (1/2)(χ σ^μ ξ̄)(φ σ_μ ψ̄)
#
#    Or in angle-bracket notation: <χψ><φξ> = ...
#
#    Simplest 2-component Fierz:
#    (θ_α θ_β) = -(1/2) ε_{αβ} (θ^γ θ_γ)
#    (θ_α θ̄_{α̇}) = -(1/2) σ^μ_{αα̇} (θ σ_μ θ̄)
# -------------------------------------------------------------------------
print("[3] Fierz rearrangement identities:")
print("""
    Four-spinor Fierz identity:
    (chi_alpha psi_beta) = -(1/2) eps_{alpha beta} (chi^gamma psi_gamma)

    [This is the antisymmetrization in 2D: any rank-2 tensor is proportional
     to epsilon in the antisymmetric part]

    More general Fierz (Srednicki notation):
    (chi psi)(phi xi) = -(chi phi)(psi xi)  [using eps to raise/lower]

    Sigma Fierz:
    (chi sigma^mu xi_bar)(phi sigma_mu psi_bar)
    = -2 (chi phi)(xi_bar psi_bar)
      - 2 (chi psi)(xi_bar phi_bar)  [incomplete without proper index structure]

    In practice the most used form (Srednicki App. B):
    (xi_bar sigmabar^mu chi)(psi sigma_mu xi_bar)
    = -2 (chi psi)(xi_bar xi_bar)   [only schematic]
""")

# -------------------------------------------------------------------------
# 4. Specific Fierz used in loop integrals (Srednicki Ch. 51, 62)
#    In scalar QED / Yukawa theories, propagator numerators produce
#    terms like (σ^μ σ̄^ν)_α^β which need to be simplified via:
#
#    σ^μ σ̄^ν = η^{μν} I - i σ^{μν}   (Cayley-Hamilton for 2-component)
#    σ̄^μ σ^ν = η^{μν} I - i σ̄^{μν}
#
#    And sigma^{mu nu} = (i/4)(sigma^mu sigmabar^nu - sigma^nu sigmabar^mu)
# -------------------------------------------------------------------------
print("[4] Product formula for sigma matrices:")
print("""
    sigma^mu sigmabar^nu = eta^{mu nu} I_{2x2} - i sigma^{mu nu}
    sigmabar^mu sigma^nu = eta^{mu nu} I_{2x2} - i sigmabar^{mu nu}

    Trace:
    Tr[sigma^mu sigmabar^nu] = -2 eta^{mu nu}
    Tr[sigma^mu sigmabar^nu sigma^rho sigmabar^sigma]
    = 2(eta^{mu nu} eta^{rho sigma} - eta^{mu rho} eta^{nu sigma} + eta^{mu sigma} eta^{nu rho})
    - 2i epsilon^{mu nu rho sigma}  [epsilon is Levi-Civita]

    These traces appear in loop calculations.
""")

# -------------------------------------------------------------------------
# 5. SUSY context: Fierz for Grassmann variables
#    In SUSY, θ^α are Grassmann-valued spinors.
#    Key identities:
#    θ_α θ_β = -(1/2) ε_{αβ} θθ      where θθ = θ^α θ_α
#    θ_α θ̄_{α̇} = -(1/2) σ^μ_{αα̇} θ σ_μ θ̄
#    θ^2 θ̄^2 = (1/4)|θθ|^2
# -------------------------------------------------------------------------
print("[5] Fierz for Grassmann spinors (SUSY context):")
print("""
    theta_alpha theta_beta = -(1/2) eps_{alpha beta} (theta theta)
    where (theta theta) = theta^alpha theta_alpha = eps^{alpha beta} theta_alpha theta_beta

    Important: theta_alpha theta_beta is symmetric in alpha<->beta
    (two Grassmann variables anticommute: theta_a theta_b = -theta_b theta_a)
    so theta_a theta_b = (1/2)(theta_a theta_b - theta_b theta_a) = (1/2) antisymmetrized
    => theta_a theta_b = -(1/2) eps_{ab} (theta^2)

    Similarly for dotted:
    theta_bar_{dala} theta_bar_{dalb} = -(1/2) eps_{dala dalb} (theta_bar)^2

    Mixed:
    theta_alpha theta_bar_{dala} = -(1/4) sigma^mu_{alpha dala} (theta sigma_mu theta_bar)

    These are the building blocks for superspace Feynman rules (Srednicki Ch. 85+)
""")

# -------------------------------------------------------------------------
# 6. Vanishing identities from antisymmetry
#    In 2 dimensions (2-component spinor space):
#    - Any antisymmetric rank-3 tensor T_{αβγ} = 0 (only 2 components)
#    - ε_{αβ} ε_{γδ} + ε_{αγ} ε_{δβ} + ε_{αδ} ε_{βγ} = 0  [Schouten]
# -------------------------------------------------------------------------
print("[6] Schouten identity from Fierz perspective:")
print("""
    In 2-component spinor space, there are only 2 basis vectors.
    Any antisymmetric combination of 3 or more indices vanishes.

    This gives the Schouten identity:
    eps_{alpha beta} eps_{gamma delta} + eps_{alpha gamma} eps_{delta beta} + eps_{alpha delta} eps_{beta gamma} = 0

    In spinor-product notation:
    <ab><cd> + <ac><db> + <ad><bc> = 0  [Schouten]

    This identity is heavily used in:
    - Simplifying n-point amplitudes
    - Proving BCJ relations
    - Loop integral reduction
""")

print("Done: 05_fierz_identities.py")
