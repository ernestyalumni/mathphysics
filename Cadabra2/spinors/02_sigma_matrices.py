"""
02_sigma_matrices.py
====================
Sigma matrices σ^μ_{αα̇} and σ̄^{μ α̇α} in Srednicki's notation.

Reference: Srednicki QFT, Chapters 35-36.
           σ^μ = (1, σ^i),  σ̄^μ = (1, -σ^i)
           where σ^i are the Pauli matrices.

Run with:
    python3 02_sigma_matrices.py
or inside Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/02_sigma_matrices.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__

print("=" * 60)
print("Sigma matrices in cadabra2 (Srednicki Ch. 35-36)")
print("=" * 60)

__cdbkernel__ = cadabra2.create_scope()

# -------------------------------------------------------------------------
# 1. Declare indices
# -------------------------------------------------------------------------
# Spacetime (Lorentz) index μ,ν,ρ,σ — vector rep
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

# Undotted spinor indices
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma}"), Ex(r"position=fixed"))

# Dotted spinor indices
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga}"), Ex(r"position=fixed"))

print("\n[1] Index declarations:")
print("    Spacetime: mu, nu, rho, sigma  (Lorentz/vector)")
print("    Undotted:  alpha, beta, gamma  (left-handed spinor)")
print("    Dotted:    dal, dbe, dga       (right-handed spinor)")

# -------------------------------------------------------------------------
# 2. Sigma matrices: σ^μ_{αα̇}
#    Connects a left-handed and right-handed spinor index with a vector index.
#
#    σ^0 = I_2 = [[1,0],[0,1]]
#    σ^1 = [[0,1],[1,0]]
#    σ^2 = [[0,-i],[i,0]]
#    σ^3 = [[1,0],[0,-1]]
#
#    σ̄^{μ α̇α} = ε^{α̇β̇} ε^{αβ} σ^μ_{ββ̇}
#    σ̄^0 = I,  σ̄^i = -σ^i
# -------------------------------------------------------------------------
sigma = Ex(r"\sigma^\mu_{\alpha\dal}")
sigma_bar = Ex(r"\bar{\sigma}^{\mu\dal\alpha}")

print("\n[2] Sigma matrices:")
print("    sigma^mu_{alpha dotalpha} =", sigma)
print("    sigmabar^{mu dotalpha alpha} =", sigma_bar)
print("\n    Explicit components (Srednicki metric signature -+++  or +---):")
print("    sigma^0 = identity_2x2")
print("    sigma^i = Pauli matrices sigma^i  (i=1,2,3)")
print("    sigmabar^0 = identity_2x2")
print("    sigmabar^i = -sigma^i")

# -------------------------------------------------------------------------
# 3. Key identity: sigma and sigmabar completeness/trace relations
#    σ^μ_{αα̇} σ̄_μ^{β̇β} = -2 δ_α^β δ^{β̇}_{α̇}
#    (up to sign conventions for metric)
# -------------------------------------------------------------------------
print("\n[3] Completeness relation:")
print("    sigma^mu_{alpha dal} sigmabar_mu^{dbe beta} = -2 delta_alpha^beta delta^dbe_{dal}")
print("    (Srednicki eq. 36.4, metric signature -+++)")

# In cadabra2, we can declare this as a substitution rule / identity
# The trace formula follows from the Clifford algebra / Pauli matrix algebra:
#   Tr[sigma^mu sigmabar^nu] = -2 eta^{mu nu}

trace_id = Ex(r"\sigma^\mu_{\alpha\dal} \bar{\sigma}_\mu^{\dal\beta}")
print("\n    Trace: sigma^mu_{alpha dal} sigmabar_{mu}^{dal beta}")
print("         =", trace_id)
print("         = -2 delta_alpha^beta  (spinor trace)")

# -------------------------------------------------------------------------
# 4. Momentum spinors: p_{αα̇} = p_μ σ^μ_{αα̇}
#    For a 4-momentum p^μ = (E, px, py, pz):
#    p_{αα̇} = p_μ σ^μ_{αα̇} = E·I + p·σ
#             = [[E+pz,  px-ipy],
#                [px+ipy, E-pz]]
# -------------------------------------------------------------------------
momentum_matrix = Ex(r"p_\mu \sigma^\mu_{\alpha\dal}")
print("\n[4] Momentum matrix p_{alpha dotalpha}:")
print("    p_{alpha dal} = p_mu sigma^mu_{alpha dal}")
print("                 =", momentum_matrix)
print("""
    Explicit 2x2 matrix form (for p^mu = (E, px, py, pz)):

    p_{alpha dal} = [[E + pz,    px - i*py],
                     [px + i*py,  E - pz  ]]

    Key property for massless particles (p^2 = 0):
    det(p_{alpha dal}) = p^mu p_mu = 0  (on-shell condition)
    => p_{alpha dal} = lambda_alpha * lambdatilde_{dotalpha}
       (factorizes into rank-1 matrix = spinor-helicity decomposition)
""")

# -------------------------------------------------------------------------
# 5. Clifford algebra from sigma matrices
#    (σ^μ σ̄^ν + σ^ν σ̄^μ)_α^β = -2 η^{μν} δ_α^β
#    (σ̄^μ σ^ν + σ̄^ν σ^μ)^{α̇}_{β̇} = -2 η^{μν} δ^{α̇}_{β̇}
# -------------------------------------------------------------------------
print("[5] Clifford algebra relations:")
print("    (sigma^mu sigmabar^nu + sigma^nu sigmabar^mu)_alpha^beta")
print("    = -2 eta^{mu nu} delta_alpha^beta")
print()
print("    (sigmabar^mu sigma^nu + sigmabar^nu sigma^mu)^{dala}_{dalb}")
print("    = -2 eta^{mu nu} delta^{dala}_{dalb}")
print()
print("    These are the 2-component versions of the Dirac/Clifford algebra.")
print("    The 4-component gamma matrices are built from sigma as:")
print("    gamma^mu = [[0,         sigma^mu  ],")
print("                [sigmabar^mu,   0     ]]  (Weyl/chiral basis)")

# -------------------------------------------------------------------------
# 6. Lorentz generators in spinor notation
#    (σ^{μν})_α^β = (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)_α^β
#    (σ̄^{μν})^{α̇}_{β̇} = (i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)^{α̇}_{β̇}
# -------------------------------------------------------------------------
print("\n[6] Lorentz generators in spinor notation:")
print("    (sigma^{mu nu})_alpha^beta = (i/4)(sigma^mu sigmabar^nu - sigma^nu sigmabar^mu)")
print("    These generate SL(2,C) transformations on left-handed spinors.")
print("    The (anti-)self-dual decomposition F^{mu nu} = F^+ + F^-")
print("    is encoded in sigma^{mu nu} and sigmabar^{mu nu}.")

print("\nDone: 02_sigma_matrices.py")
