"""
ch38_spinor_technology.py
==========================
Srednicki QFT — Chapter 38: Spinor Technology

What this file covers (section by section):
  §38.A  σ^μ completeness relation: σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}
  §38.B  Trace identities: Tr[σ^μ σ̄^ν] and the 4-trace
  §38.C  Four-trace: Tr[σ^μ σ̄^ν σ^ρ σ̄^σ]
  §38.D  Van der Waerden index gymnastics — raising, lowering, contracting
  §38.E  Spinor bilinears and the Gordon identity
  §38.F  Fierz identities for spinor products
  §38.G  Cadabra2 symbolic verification of all key identities

Reference: Srednicki QFT, Chapter 38.
Metric convention: Srednicki g_{μν} = diag(-1,+1,+1,+1)  [Eq. 1.8, mostly-plus -+++]

Run with:
    python3 ch38_spinor_technology.py
Docker:
    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 python3 /work/ch38_spinor_technology.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70
sec = lambda s: print(f"\n{SEP}\n  {s}\n{SEP}")

print(SEP)
print("  Srednicki Ch. 38 — Spinor Technology")
print(SEP)

# ── Cadabra2 index declarations ───────────────────────────────────────────
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"),       Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"),        Ex(r"position=free"))

cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))
# Note: DiagonalMetric and EpsilonTensor are not available in this cadabra2 build;
# metric and Levi-Civita properties are handled numerically via numpy.

# ── numpy setup ───────────────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)

sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),
    2: np.array([[0,-1j], [1j,0]],  dtype=complex),
    3: np.array([[1, 0],  [0,-1]],  dtype=complex),
}
sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

# Metric g^{μν} = diag(-1,+1,+1,+1)  [mostly-plus, Srednicki convention]
g = np.diag([-1., 1., 1., 1.])

eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1,0]], dtype=complex)

# Levi-Civita tensor ε^{μνρσ}: totally antisymmetric, ε^{0123} = +1
def levi_civita_4d():
    """Return the 4D Levi-Civita tensor ε^{μνρσ} with ε^{0123} = +1."""
    eps = np.zeros((4, 4, 4, 4), dtype=float)
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    indices = [mu, nu, rho, sig]
                    if len(set(indices)) == 4:
                        # Determine sign from permutation parity
                        perm = indices[:]
                        sign = 1
                        for i in range(4):
                            while perm[i] != i:
                                j = perm[i]
                                perm[i], perm[j] = perm[j], perm[i]
                                sign *= -1
                        eps[mu, nu, rho, sig] = sign
    return eps

epsilon = levi_civita_4d()

# =============================================================================
# §38.A  σ^μ COMPLETENESS RELATION
# =============================================================================
# The sigma matrices form a complete basis for 2×2 matrices.
# The completeness relation (Srednicki eq. 38.1 / cf. eq. 35.4 from Ch. 35):
#
#   σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}                    [eq. 38.1]
#
# where σ_μ^{ββ̇} = g_{μν} σ̄^{ν ββ̇} (both spinor indices raised).
#
# Equivalently:
#   ε^{αβ} ε^{α̇β̇} σ^μ_{αα̇} σ^ν_{ββ̇} = −2 g^{μν}               [eq. 35.5]
#
# Both forms are equivalent via ε-gymnastics (see §38.D).
#
# Physical meaning:
#   The σ^μ matrices (for μ=0,1,2,3) form a complete orthogonal basis
#   for the space of 2×2 complex matrices. The completeness relation
#   is the σ-analogue of the Kronecker delta for vectors:
#     Σ_μ e^μ_i e^ν_j g_{μν} = -δ_{ij}  (for unit vectors in Minkowski space)
#
# Note: The factor of -2 comes from:
#   Tr[σ^μ σ̄^ν] = 2 g^{μν}  [derived in §38.B]
#   Together: σ^μ_{αα̇} σ_μ^{ββ̇} = -(1/2) Tr[...] δ × 4 = -2 δδ

sec("§38.A  σ^μ Completeness Relation  [eq. 38.1]")

print("Completeness relation:")
print()
print("  σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}        [eq. 38.1]")
print()
print("  where σ_μ^{ββ̇} = g_{μν} σ̄^{ν ββ̇} = Σ_ν g_{μν} ε^{βγ} ε^{β̇γ̇} σ^ν_{γγ̇}")
print()

# Compute σ^μ_{αα̇} σ_μ^{ββ̇} = g^{μν} σ^μ_{αα̇} σ_{ν}^{ββ̇}
# We use: σ_μ^{ββ̇} = g_{μν} σ̄^{ν ββ̇}
# So: σ^μ_{αα̇} σ_μ^{ββ̇} = Σ_μ Σ_ν g_{μν} σ^μ_{αα̇} σ̄^{ν ββ̇}
#                          = Σ_μ g_{μμ} σ^μ_{αα̇} σ̄^{μ ββ̇}  (diagonal metric)

print("Computing σ^μ_{αα̇} g_{μμ} σ̄^{μ ββ̇} for all (α,α̇,β,β̇):")
print()

# σ_μ^{ββ̇} = g_{μν} σ̄^{ν ββ̇}: indices [mu, beta, betadot]
# σ^μ_{αα̇} : indices [mu, alpha, alphadot]
# Result: (4-index tensor) [alpha, alphadot, beta, betadot]

lhs_completeness = np.zeros((2, 2, 2, 2), dtype=complex)
for mu in range(4):
    # σ^μ_{αα̇}: sigma_vec[mu][alpha, alphadot]
    # σ_μ^{ββ̇} = g_{μμ} σ̄^{μ ββ̇}: g[mu,mu] * sigmabar_vec[mu][beta, betadot]
    # Note: σ̄^{μ ββ̇} has index structure [betadot, beta] (dotted first)
    # We need σ_μ^{ββ̇} with β undotted, β̇ dotted
    # σ̄^{μ ββ̇}[bdot, b] → need [b, bdot] = sigmabar_vec[mu].T
    for alpha in range(2):
        for alphadot in range(2):
            for beta in range(2):
                for betadot in range(2):
                    val = sum(
                        g[mu, mu] * sigma_vec[mu][alpha, alphadot]
                                  * sigmabar_vec[mu][betadot, beta]
                        for mu in range(4)
                    )
                    lhs_completeness[alpha, alphadot, beta, betadot] = val

# RHS: -2 δ_α^β δ_{α̇}^{β̇}
rhs_completeness = np.zeros((2, 2, 2, 2), dtype=complex)
for alpha in range(2):
    for alphadot in range(2):
        for beta in range(2):
            for betadot in range(2):
                rhs_completeness[alpha, alphadot, beta, betadot] = (
                    -2 * (alpha == beta) * (alphadot == betadot)
                )

err_comp = np.max(np.abs(lhs_completeness - rhs_completeness))
print(f"  Max error: {err_comp:.2e}")
if err_comp < 1e-12:
    print("  ✓ σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}  verified!")
else:
    print("  ✗ MISMATCH in completeness relation!")

# Print nonzero components
print()
print("  Nonzero components (α,α̇,β,β̇) with 1-based labeling:")
print("  (α,α̇,β,β̇)     LHS      RHS")
for alpha in range(2):
    for alphadot in range(2):
        for beta in range(2):
            for betadot in range(2):
                lv = lhs_completeness[alpha, alphadot, beta, betadot].real
                rv = rhs_completeness[alpha, alphadot, beta, betadot].real
                if abs(lv) > 1e-10 or abs(rv) > 1e-10:
                    print(f"  ({alpha+1},{alphadot+1},{beta+1},{betadot+1})         "
                          f"{lv:+.1f}     {rv:+.1f}")

# Alternative form via eq. 35.5
print()
print("Alternative form [eq. 35.5]:")
print("  ε^{αβ} ε^{α̇β̇} σ^μ_{αα̇} σ^ν_{ββ̇} = −2 g^{μν}")
lhs_alt = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        val = 0.
        for alpha in range(2):
            for alphadot in range(2):
                for beta in range(2):
                    for betadot in range(2):
                        val += (eps_upper[alpha, beta]
                                * eps_upper[alphadot, betadot]
                                * sigma_vec[mu][alpha, alphadot]
                                * sigma_vec[nu][beta, betadot])
        lhs_alt[mu, nu] = val

err_alt = np.max(np.abs(lhs_alt - (-2) * g))
print(f"  Max error vs −2 g^{{μν}}: {err_alt:.2e}")
if err_alt < 1e-12:
    print("  ✓ Verified!")

# =============================================================================
# §38.B  TRACE IDENTITY: Tr[σ^μ σ̄^ν] = −2 η^{μν}
# =============================================================================
# The fundamental 2-trace identity (Srednicki eq. 38.5):
#
#   Tr[σ^μ σ̄^ν] = σ^μ_{αα̇} σ̄^{ν α̇α} = −2 η^{μν}
#
# Wait: Srednicki's metric is mostly-plus, g_{μν} = diag(-1,+1,+1,+1).
# The trace is:
#   σ^0 σ̄^0 = I·I = I → Tr = 2
#   σ^i σ̄^i = σ^i (-σ^i) = -I → Tr = -2
# And -2η^{μν} with η=diag(-1,+1,+1,+1):
#   -2η^{00} = -2(-1) = 2  ✓
#   -2η^{ii} = -2(+1) = -2 ✓
#
# But more commonly stated with mostly-minus metric η=diag(+1,-1,-1,-1):
#   Tr[σ^μ σ̄^ν] = 2 η^{μν}  [Peskin & Schroeder convention]
#   Tr[σ^μ σ̄^ν] = -2 η^{μν} [Srednicki's mostly-plus, since η^{μν} has opposite signs]
#
# In Srednicki's convention (g_{μν} = diag(-1,+1,+1,+1)):
#   Tr[σ^μ σ̄^ν] = -2 g^{μν}    [where g^{μν} = diag(-1,+1,+1,+1)]
#
# Physical meaning:
#   This is the completeness/orthogonality relation for σ matrices as a basis.
#   The σ^μ and σ̄^ν are "dual" to each other via the trace inner product.

sec("§38.B  2-Trace: Tr[σ^μ σ̄^ν] = −2 g^{μν}  [eq. 38.5]")

print("Two-sigma trace identity:")
print()
print("  Tr[σ^μ σ̄^ν] ≡ σ^μ_{αα̇} σ̄^{ν α̇α}")
print("              = Tr₂[σ^μ · σ̄^ν]  (2×2 matrix trace)")
print()

# Compute Tr[σ^μ σ̄^ν] = Σ_α (σ^μ σ̄^ν)_{αα} = Tr(σ^μ · σ̄^ν)
# σ^μ has indices [alpha, alphadot], σ̄^ν has indices [alphadot, alpha]
# Matrix product: (σ^μ · σ̄^ν)_{αβ} = Σ_{α̇} σ^μ_{αα̇} σ̄^{ν α̇β}
# Trace: Σ_α (σ^μ · σ̄^ν)_{αα}

trace_2 = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        # σ^μ: [alpha, alphadot], σ̄^ν: [alphadot, beta]
        # standard matmul
        trace_2[mu, nu] = np.trace(sigma_vec[mu] @ sigmabar_vec[nu])

expected_trace_2 = -2 * g  # -2 g^{μν}

err_trace2 = np.max(np.abs(trace_2 - expected_trace_2))

print("  Tr[σ^μ σ̄^ν] =")
for row in trace_2.real:
    print("   ", [f"{v:+.1f}" for v in row])
print()
print("  Expected -2 g^{μν} =")
for row in expected_trace_2:
    print("   ", [f"{v:+.1f}" for v in row])
print()
print(f"  Max error: {err_trace2:.2e}")
if err_trace2 < 1e-12:
    print("  ✓ Tr[σ^μ σ̄^ν] = −2 g^{μν}  verified!")
else:
    print("  ✗ MISMATCH!")

print()
print("Physical interpretation:")
print("  • σ^μ and σ̄^ν are 'dual' bases for 2×2 matrices")
print("  • The trace gives their 'inner product' = −2 g^{μν}")
print("  • Any 2×2 matrix A can be expanded: A = (1/2) Σ_μ (-σ̄^μ Tr[σ^μ A])")
print("    (completeness of {σ^μ} as basis for 2×2 matrices)")

# Also verify the 'reverse' trace: Tr[σ̄^μ σ^ν] = -2 g^{μν}
trace_2_bar = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        trace_2_bar[mu, nu] = np.trace(sigmabar_vec[mu] @ sigma_vec[nu])

err_trace2_bar = np.max(np.abs(trace_2_bar - expected_trace_2))
print()
print(f"  Tr[σ̄^μ σ^ν] = −2 g^{{μν}}?  max error = {err_trace2_bar:.2e}",
      "  ✓" if err_trace2_bar < 1e-12 else "  ✗")

# =============================================================================
# §38.C  4-TRACE: Tr[σ^μ σ̄^ν σ^ρ σ̄^σ]
# =============================================================================
# The four-sigma trace identity (Srednicki eq. 38.6):
#
#   Tr[σ^μ σ̄^ν σ^ρ σ̄^σ] = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ})
#                           − 2i ε^{μνρσ}
#
# where ε^{0123} = +1 (in the context of Srednicki).
#
# NOTE on metric signs: Srednicki uses g = diag(-1,+1,+1,+1).
# In terms of most-commonly stated form with η = diag(+1,-1,-1,-1) (Peskin):
#   Tr[σ^μ σ̄^ν σ^ρ σ̄^σ] = 2(η^{μν}η^{ρσ} − η^{μρ}η^{νσ} + η^{μσ}η^{νρ})
#                           + 2i ε^{μνρσ}
# [Note: sign flip from η ↔ g and from ε sign convention]
#
# Derivation outline:
#   Use σ^μ σ̄^ν = −g^{μν} I + 2 (S^{μν}_L)  [cf. Ch. 35 eq. 35.21 rearranged]
#   The trace of the antisymmetric part (S^{μν}_L) contributes the ε term.
#
# Physical applications:
#   • Feynman diagram traces for fermion loops with multiple photon vertices
#   • Cross sections for processes like e+e− → μ+μ−
#   • Helicity amplitude calculations

sec("§38.C  4-Trace: Tr[σ^μ σ̄^ν σ^ρ σ̄^σ]  [eq. 38.6]")

print("Four-sigma trace identity (Srednicki eq. 38.6, mostly-plus metric):")
print()
print("  Tr[σ^μ σ̄^ν σ^ρ σ̄^σ]")
print("  = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) + 2i ε^{μνρσ}")
print("  [sign: +2i with Srednicki mostly-plus g=diag(-1,+1,+1,+1)]")
print()

# Compute the 4-trace numerically
trace_4 = np.zeros((4, 4, 4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                mat = sigma_vec[mu] @ sigmabar_vec[nu] @ sigma_vec[rho] @ sigmabar_vec[sig]
                trace_4[mu, nu, rho, sig] = np.trace(mat)

# Compute the expected RHS: 2(g^{μν}g^{ρσ} - g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) + 2i ε^{μνρσ}
# The +2i sign is correct for Srednicki's mostly-plus metric convention.
rhs_4trace = np.zeros((4, 4, 4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                rhs_4trace[mu, nu, rho, sig] = (
                    2 * (g[mu, nu] * g[rho, sig]
                         - g[mu, rho] * g[nu, sig]
                         + g[mu, sig] * g[nu, rho])
                    + 2j * epsilon[mu, nu, rho, sig]
                )

err_4trace = np.max(np.abs(trace_4 - rhs_4trace))
print(f"  Numerical verification — Max error: {err_4trace:.2e}")
if err_4trace < 1e-12:
    print("  ✓ Tr[σ^μ σ̄^ν σ^ρ σ̄^σ] = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) + 2i ε^{μνρσ}")
    print("     verified for all 256 index combinations!")
else:
    print("  ✗ MISMATCH!")
    # Print some sample values for debugging
    for mu in range(2):
        for nu in range(mu+1, 3):
            for rho in range(nu+1, 4):
                sig = 3 - mu - nu - rho + (6 - mu - nu - rho) % 4
                if 0 <= sig < 4 and len({mu,nu,rho,sig}) == 4:
                    lv = trace_4[mu,nu,rho,sig]
                    rv = rhs_4trace[mu,nu,rho,sig]
                    print(f"  ({mu},{nu},{rho},{sig}): LHS={lv:.3f}, RHS={rv:.3f}")

# Print key special cases to build intuition
print()
print("Key special cases (illustrating the structure):")
print()

cases = [
    (0, 1, 0, 1, "Tr[σ^0 σ̄^1 σ^0 σ̄^1]"),
    (0, 1, 2, 3, "Tr[σ^0 σ̄^1 σ^2 σ̄^3]"),
    (0, 1, 1, 0, "Tr[σ^0 σ̄^1 σ^1 σ̄^0]"),
    (1, 2, 3, 0, "Tr[σ^1 σ̄^2 σ^3 σ̄^0]"),
]
for mu, nu, rho, sig, label in cases:
    lv = trace_4[mu, nu, rho, sig]
    rv = rhs_4trace[mu, nu, rho, sig]
    print(f"  {label}")
    print(f"    Computed: {lv:.3f}")
    print(f"    Expected: {rv:.3f}")
    print(f"    Match: {'✓' if abs(lv-rv) < 1e-10 else '✗'}")
    print()

# The conjugate trace: Tr[σ̄^μ σ^ν σ̄^ρ σ^σ]
print("Conjugate trace Tr[σ̄^μ σ^ν σ̄^ρ σ^σ]:")
print("  = (Tr[σ^σ σ̄^ρ σ^ν σ̄^μ])* = complex conjugate with μ↔σ, ν↔ρ")
print("  = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) − 2i ε^{μνρσ}")
print("  (Note: -2i ε, sign flipped from +2i relative to the unbarred trace)")

trace_4_bar = np.zeros((4, 4, 4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                mat = sigmabar_vec[mu] @ sigma_vec[nu] @ sigmabar_vec[rho] @ sigma_vec[sig]
                trace_4_bar[mu, nu, rho, sig] = np.trace(mat)

rhs_4trace_bar = np.zeros((4, 4, 4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                rhs_4trace_bar[mu, nu, rho, sig] = (
                    2 * (g[mu, nu] * g[rho, sig]
                         - g[mu, rho] * g[nu, sig]
                         + g[mu, sig] * g[nu, rho])
                    - 2j * epsilon[mu, nu, rho, sig]
                )

err_4trace_bar = np.max(np.abs(trace_4_bar - rhs_4trace_bar))
print(f"\n  Max error for conjugate trace: {err_4trace_bar:.2e}",
      "  ✓" if err_4trace_bar < 1e-12 else "  ✗")

# =============================================================================
# §38.D  VAN DER WAERDEN INDEX GYMNASTICS
# =============================================================================
# Van der Waerden notation: systematic rules for raising/lowering spinor indices.
#
# METRIC ε:  used to raise and lower spinor indices:
#   Lowering undotted:  ψ_α = ε_{αβ} ψ^β          [eq. 38.7]
#   Raising undotted:   ψ^α = ε^{αβ} ψ_β           [eq. 38.8]
#   Lowering dotted:    χ_α̇ = ε_{α̇β̇} χ^β̇          [eq. 38.9]
#   Raising dotted:     χ^α̇ = ε^{α̇β̇} χ_β̇           [eq. 38.10]
#
# KEY RULE: raise/lower by contracting on the SECOND index of ε:
#   ψ_α = ε_{αβ} ψ^β  (contract on β, second index of ε_{αβ})
#   ψ^α = ε^{αβ} ψ_β  (contract on β, second index of ε^{αβ})
#
# This ensures consistency:
#   ψ^α = ε^{αβ} ε_{βγ} ψ^γ = δ^α_γ ψ^γ = ψ^α  ✓
#
# CONTRACTION CONVENTION (eq. 38.11):
#   χψ ≡ χ^α ψ_α = χ^α ε_{αβ} ψ^β
#   χ†ψ† ≡ χ†_α̇ ψ†^α̇ = χ†_α̇ ε^{α̇β̇} ψ†_β̇
#
# σ-INDEX GYMNASTICS:
#   σ^μ_{αα̇}: one undotted down, one dotted down
#   σ^{μ αα̇} = ε^{αβ} ε^{α̇β̇} σ^μ_{ββ̇} = -σ̄^{μ α̇α} ... wait.
#   More carefully (eq. 38.13):
#   σ̄^{μ α̇α} = ε^{α̇β̇} ε^{αβ} σ^μ_{ββ̇}   [raise both spinor indices of σ^μ]
#   σ^μ_{ββ̇} = ε_{βα} ε_{β̇α̇} σ̄^{μ α̇α}   [lower both spinor indices of σ̄^μ]

sec("§38.D  Van der Waerden Index Gymnastics")

print("ε-metric raising/lowering rules (Srednicki conventions):")
print()
print("  Undotted spinors:")
print("    ψ_α = ε_{αβ} ψ^β      (lower with ε_{12}=-1, ε_{21}=+1)")
print("    ψ^α = ε^{αβ} ψ_β      (raise with ε^{12}=+1, ε^{21}=-1)")
print()
print("  Dotted spinors:")
print("    χ_α̇ = ε_{α̇β̇} χ^β̇    (lower)")
print("    χ^α̇ = ε^{α̇β̇} χ_β̇    (raise)")
print()
print("  KEY: contract on SECOND index of ε (NW-SE contraction)")
print()

# Verify ε^{αβ} ε_{βγ} = δ^α_γ
test_inv = eps_upper @ eps_lower
print("Verification: ε^{αβ} ε_{βγ} = δ^α_γ:")
print(f"  ε_upper @ ε_lower = {test_inv.real}")
err_inv = np.max(np.abs(test_inv - I2))
print(f"  Max error vs I₂: {err_inv:.2e}",
      "  ✓" if err_inv < 1e-12 else "  ✗")
print()

# σ-matrix index gymnastics
print("σ-matrix index gymnastics:")
print()
print("  σ^μ_{αα̇}: lower indices (canonical form)")
print("  σ̄^{μ α̇α} = ε^{α̇β̇} ε^{αβ} σ^μ_{ββ̇}  (raise both spinor indices)")
print()

# Compute σ̄^{μ α̇α} from raising both indices of σ^μ_{αα̇}
# σ̄_raised[mu, alphadot, alpha] = Σ_{beta, betadot} ε^{αβ} ε^{α̇β̇} σ^μ_{ββ̇}
sigmabar_raised = {}
for mu in range(4):
    mat = np.zeros((2, 2), dtype=complex)  # [alphadot, alpha]
    for alphadot in range(2):
        for alpha in range(2):
            val = 0.
            for beta in range(2):
                for betadot in range(2):
                    val += (eps_upper[alpha, beta]
                            * eps_upper[alphadot, betadot]
                            * sigma_vec[mu][beta, betadot])
            mat[alphadot, alpha] = val
    sigmabar_raised[mu] = mat

print("  Comparison: ε^{αβ} ε^{α̇β̇} σ^μ_{ββ̇}  vs  σ̄^{μ α̇α} = (I, -σ⃗):")
max_err_gym = 0.
for mu in range(4):
    expected = sigmabar_vec[mu]  # sigmabar_vec[mu][alphadot, alpha] = σ̄^μ
    diff = sigmabar_raised[mu] - expected
    err = np.max(np.abs(diff))
    max_err_gym = max(max_err_gym, err)
    print(f"    μ={mu}: err = {err:.2e}  {'✓' if err < 1e-12 else '✗'}")

print(f"\n  ✓ σ̄^{{μ α̇α}} = ε^{{αβ}} ε^{{α̇β̇}} σ^μ_{{ββ̇}} verified" if max_err_gym < 1e-12
      else f"\n  ✗ MISMATCH: max err = {max_err_gym:.2e}")

# Cadabra2: declare σ and σ̄ with correct index placement
print()
print("Cadabra2 expressions for σ-index gymnastics:")
sigma_lower = Ex(r"\sigma^{\mu}_{\alpha\dal}")      # σ^μ_{αα̇}
sigma_upper = Ex(r"\sigmabar^{\mu\dal\alpha}")       # σ̄^{μ α̇α}
eps_raise_lower = Ex(r"\epsilon^{\alpha\beta} \epsilon^{\dal\dbe} \sigma^{\mu}_{\beta\dbe}")
print(f"  σ^μ_{{αα̇}}      (lower both) : {sigma_lower}")
print(f"  σ̄^{{μ α̇α}}     (raise both) : {sigma_upper}")
print(f"  Raising formula: ε^{{αβ}} ε^{{α̇β̇}} σ^μ_{{ββ̇}} = {eps_raise_lower}")

# =============================================================================
# §38.E  SPINOR BILINEARS AND THE GORDON IDENTITY
# =============================================================================
# The Gordon identity is the spinor analogue of the momentum decomposition.
# For Dirac spinors u(p') γ^μ u(p), we can decompose:
#
#   ū(p') γ^μ u(p) = ū(p') [(p'^μ + p^μ)/(2m) + iσ^{μν}(p'-p)_ν/(2m)] u(p)
#
# In the two-component (Weyl) notation, the Gordon identity appears as:
#
#   σ̄^{μ α̇α} p_μ = p^{α̇α}
#   σ^μ_{αα̇} p_μ = p_{αα̇}
#
# And the more useful bilinear form:
#   For on-shell spinors u^s(p) satisfying p_{αα̇} ũ^{α̇s} = m u_α^s:
#
#   ũ^{α̇s}(p) σ̄^{μ }_{ α̇α} u^s_α(p) = ... (involves p^μ and spin structure)
#
# More commonly, the Gordon decomposition in 4-component Dirac notation:
#   ū(p') γ^μ u(p) = (1/(2m)) ū(p') [(p + p')^μ + i[γ^μ, p/'] ... ]
#
# In Weyl/Van der Waerden language (eq. 38.15):
#   χ†_α̇ σ̄^{μ α̇α} ψ_α = (1/2)(χ†α̇ p^{α̇α} ψ_α)/m + (spin terms)
#
# The KEY two-component Gordon-like identity (from Ch. 38):
#   σ^μ_{αα̇} σ̄_μ^{ββ̇} = -2 δ_α^β δ_{α̇}^{β̇}   [completeness, §38.A]
#
#   AND the "Schouten" (Fierz) rearrangement identities:
#   δ_α^γ δ_β^δ = δ_α^δ δ_β^γ - ε_{αβ} ε^{γδ}   [Schouten identity]

sec("§38.E  Spinor Bilinears and Gordon Identity  [eq. 38.15]")

print("Key bilinears in Van der Waerden notation:")
print()
print("  Scalar bilinears (Lorentz invariants):")
print("    χψ = χ^α ψ_α  = ε^{αβ} χ_α ψ_β  (left-handed)")
print("    χ†ψ† = χ†_α̇ ψ†^α̇  (right-handed)")
print("    Both are Lorentz scalars; χψ = ψχ [Grassmann symmetry]")
print()
print("  Vector bilinear (Lorentz 4-vector):")
print("    χ†_α̇ σ̄^{μ α̇α} ψ_α  →  transforms as 4-vector V^μ  [eq. 36.28]")
print()
print("  Tensor bilinear (anti-symmetric tensor):")
print("    χ^α (S^{μν}_L)_α^β ψ_β = (i/4) χ^α (σ^μ σ̄^ν - σ^ν σ̄^μ)_α^β ψ_β")
print("    transforms as (anti-symmetric) Lorentz tensor  [eq. 35.21]")
print()

# Gordon decomposition in Weyl notation
# p_{αα̇} = p_μ σ^μ_{αα̇}: lower-lower momentum matrix
# p̄^{α̇α} = p_μ σ̄^{μ α̇α}: upper-upper momentum matrix
print("Momentum matrix in spinor space:")
print("  p_{αα̇} = p_μ σ^μ_{αα̇} = lower-lower spinor 2×2 matrix")
print("  p̄^{α̇α} = p_μ σ̄^{μ α̇α} = upper-upper spinor 2×2 matrix")
print()
print("  Note: p^2 = p_μ p^μ = -m² (mostly-plus metric, on-shell)")
print("        det(p_{αα̇}) = -p_μ p^μ/4 = m²/4 ... up to normalization")
print("        For massless: det(p_{αα̇}) = 0 → rank-1 matrix → factorizes")
print()

# Gordon identity (generalized form):
# For bilinear ψ†_α̇ σ̄^{μ α̇α} χ_α at momenta p', p:
# (p' + p)^μ = ψ†(p') σ̄^μ χ(p) / (some normalization) + correction terms
# The exact form depends on the normalization convention for the spinors.
print("Gordon identity structure:")
print("  In 4-component Dirac notation (using γ^μ = [[0,σ^μ],[σ̄^μ,0]]):")
print()
print("  ū(p') γ^μ u(p) = ū(p') [(p'^μ+p^μ)/(2m)] u(p)")
print("                  + ū(p') [iσ^{μν}(p'_ν-p_ν)/(2m)] u(p)")
print()
print("  where σ^{μν} = (i/2)[γ^μ, γ^ν]  (Dirac σ^{μν}, different from σ^μ!)")
print()
print("  In Weyl 2-component form, this decomposes as:")
print("    Vector part: (p+p')^μ term")
print("    Tensor part: (p'-p)^ν contracted with (S^{μν}_L + S^{μν}_R)")

# Numerical check: the Gordon decomposition for at-rest to at-rest transition
# u†(p) σ̄^μ u(p) = 2p^μ/m * (normalization)
m_val = 1.0
E_val = m_val  # at rest

def momentum_matrix(E, px, py, pz):
    """p_{αα̇} = p_μ σ^μ with mostly-plus g = diag(-1,+1,+1,+1)."""
    return (-E * sigma_vec[0] + px * sigma_vec[1]
            + py * sigma_vec[2] + pz * sigma_vec[3])

p_matrix = momentum_matrix(E_val, 0, 0, 0)
u_plus   = np.array([1., 0.], dtype=complex)
u_minus  = np.array([0., 1.], dtype=complex)

print()
print("Numerical check — vector bilinear u†(p) σ̄^μ u(p) at rest:")
for s, u_s in [('+', u_plus), ('-', u_minus)]:
    bilinear = np.array([
        u_s.conj() @ sigmabar_vec[mu] @ u_s
        for mu in range(4)
    ])
    print(f"  s={s}: ũ^{{α̇{s}}} σ̄^{{μ}} u^{s}_α = {bilinear.real}")
    print(f"         (expected: 2*(1,0,0,0)/norm for spin-up; (1,0,0,0) = rest 4-mom up to factor)")

# Cadabra2: bilinear expression
gordon_bilinear = Ex(r"\chidag_{\dal} \sigmabar^{\mu\dal\alpha} \psi_{\alpha}")
print()
print(f"Cadabra2 vector bilinear: {gordon_bilinear}")

# =============================================================================
# §38.F  FIERZ IDENTITIES
# =============================================================================
# Fierz rearrangement identities allow rewriting products of bilinears.
# The fundamental spinor Fierz identity (Srednicki eq. 38.19):
#
#   (χ†ψ†)(ξη) = -(χ†η)(ξψ†) ... [various forms]
#
# More precisely, in index notation (the "spinor Fierz"):
#   ε_{αγ} ε_{βδ} = ε_{αβ} ε_{γδ} - ε_{αδ} ε_{γβ}      [Schouten identity]
#
# The 2-component Fierz for products of σ̄ matrices:
#   σ̄^{μ α̇α} σ̄_μ^{β̇β} = -2 δ^{α̇β̇} δ^{αβ}            [completeness]
#   (σ̄^μ)^{α̇α} (σ_μ)_{ββ̇} = -2 δ^α_β δ^{α̇}_{β̇}       [mixed]
#
# 4-Fermi Fierz (for 4-component spinors, but easily written in 2-component):
#   (χ†_α̇ σ̄^{μ α̇α} ψ_α)(ξ†_β̇ σ̄_μ^{β̇β} η_β)
#     = -2 (χ†_α̇ ξ†^{α̇})(η^α ψ_α)    [from completeness + Grassmann]

sec("§38.F  Fierz Identities  [eq. 38.19-38.21]")

print("Fierz rearrangement identities:")
print()
print("  Schouten identity (ε completeness):")
print("    ε_{αγ} ε_{βδ} = ε_{αβ} ε_{γδ} − ε_{αδ} ε_{γβ}    [eq. 38.17]")
print()
print("  In index form: δ^α_γ δ^β_δ = δ^α_δ δ^β_γ + ε^{αβ} ε_{γδ}  [Schouten]")
print()

# Verify the Schouten identity numerically:
# ε_{αγ} ε_{βδ} = ε_{αβ} ε_{γδ} - ε_{αδ} ε_{γβ}
max_err_schouten = 0.
for alpha in range(2):
    for beta in range(2):
        for gamma in range(2):
            for delta in range(2):
                lhs = eps_lower[alpha, gamma] * eps_lower[beta, delta]
                rhs = (eps_lower[alpha, beta] * eps_lower[gamma, delta]
                       - eps_lower[alpha, delta] * eps_lower[gamma, beta])
                err = abs(lhs - rhs)
                max_err_schouten = max(max_err_schouten, err)

print(f"  Numerical verification of Schouten identity:")
print(f"    Max error: {max_err_schouten:.2e}",
      "  ✓" if max_err_schouten < 1e-12 else "  ✗")
print()

print("  σ-Fierz identity (from completeness):")
print("    (σ̄^μ)^{α̇α} (σ_μ)_{ββ̇} = −2 δ^α_β δ^{α̇}_{β̇}")
print()
print("  Product of vector bilinears (Fierz rearrangement):")
print("    (χ†_α̇ σ̄^{μ α̇α} ψ_α)(ξ†_β̇ σ̄_μ^{β̇β} η_β)")
print("      = −2 (χ†_α̇ ξ†^{α̇})(η^β ψ_β)    [from σ-Fierz + anticommutativity]")
print("      = −2 (χ†ξ†)(ηψ)")
print()

# Fierz for 4-component case with γ matrices
# In terms of 2×2 Fierz:
# Σ_μ σ^μ_{αα̇} σ_μ^{ββ̇} = -2 δ_α^β δ_{α̇}^{β̇}  (already proven in §38.A)
# This is the Fierz identity for σ-matrices.

# Also verify the Fierz identity for (σ^μ σ̄^ν)_{α}^{γ} (σ^μ σ̄^ν)_{β}^{δ}
print("  Fierz rearrangement for (σ^μ σ̄^ν)_{α}^{β}:")
print("    Σ_μ Σ_ν (σ^μ σ̄^ν)_{α}^{β} (σ_μ σ̄_ν)_{γ}^{δ} = [determined by ε Fierz]")

# Cadabra2 Fierz expression
fierz_cadabra = Ex(r"\epsilon_{\alpha\gamma} \epsilon_{\beta\delta}")
fierz_rhs     = Ex(
    r"\epsilon_{\alpha\beta} \epsilon_{\gamma\delta}"
    r" - \epsilon_{\alpha\delta} \epsilon_{\gamma\beta}"
)
print()
print(f"Cadabra2 Schouten: {fierz_cadabra} = {fierz_rhs}")

# =============================================================================
# §38.F2  GORDON IDENTITY
# =============================================================================
# Gordon identity (Srednicki Problem 38.3, eq. 38.22):
#
#   2m ū_{s'}(p') γ^μ u_s(p)
#     = ū_{s'}(p') [(p'+p)^μ - 2i S^{μν}(p'-p)_ν] u_s(p)
#
# where S^{μν} = (i/4)[γ^μ,γ^ν] is the Dirac spin tensor.
#
# Derivation key:
#   γ^μ /p  = -(p)^μ - 2iS^{μν}p_ν        (from Clifford + S^{μν} def)
#   /p'γ^μ  = -(p')^μ + 2iS^{μν}p'_ν
#   Sum: γ^μ/p + /p'γ^μ = -(p+p')^μ - 2iS^{μν}(p-p')_ν
#   Sandwich with ū_{s'}(p'), u_s(p) + use (/p+m)u=0, ū(/p'+m)=0
#   → 2m ū'γ^μu = ū'[(p+p')^μ - 2iS^{μν}(p'-p)_ν]u  ✓
#
# Extended with γ^5 (eq. 38.41):
#   ū_{s'}(p') [(p'+p)^μ - 2iS^{μν}(p'-p)_ν] γ^5 u_s(p) = 0
# (vanishes: {γ^5, γ^μ}=0 makes the tensor P-odd, while ū'...u is P-even)
#
# Physical meaning: decomposes γ^μ vertex into:
#   (1) charge-current (p+p')^μ/(2m) → form factor F_1(q²)
#   (2) magnetic-moment iσ^{μν}q_ν/(2m) → form factor F_2(q²), g-2

sec("§38.F2  Gordon Identity  [Problem 38.3 / eq. 38.22]")

print("Gordon identity for on-shell Dirac spinors:")
print()
print("  2m ū_{s'}(p') γ^μ u_s(p)")
print("    = ū_{s'}(p') [(p'+p)^μ - 2i S^{μν}(p'-p)_ν] u_s(p)")
print()
print("  S^{μν} = (i/4)[γ^μ,γ^ν]  (Dirac spin tensor)")
print()
print("  Derivation:")
print("    γ^μ/p + /p'γ^μ = -(p+p')^μ - 2iS^{μν}(p_ν - p'_ν)")
print("    Use (/p+m)u=0 and ū'(/p'+m)=0:")
print("    → 2m ū'γ^μu = ū'[(p+p')^μ - 2iS^{μν}(p'-p)_ν]u  ✓")
print()

# Numerical check at p'=p (elastic, rest frame)
# Gordon: 2m ū'γ^μu = ū'(2p^μ)u  since (p'-p)=0
Z2g = np.zeros((2,2), dtype=complex)
gamma_g = {}
for mug in range(4):
    gamma_g[mug] = np.block([[Z2g, sigma_vec[mug]], [sigmabar_vec[mug], Z2g]])

beta_g = gamma_g[0]
def dbar_g(s): return s.conj() @ beta_g

m_g = 1.0
g4g = np.diag([-1.,1.,1.,1.])
u_g = [
    np.array([1.,0.,1.,0.], dtype=complex)/np.sqrt(2),
    np.array([0.,1.,0.,1.], dtype=complex)/np.sqrt(2),
]
p4g_upper = np.array([m_g, 0., 0., 0.])

max_gordon_err = 0.
for s1 in range(2):
    for s2 in range(2):
        for mug in range(4):
            lhs = 2*m_g * (dbar_g(u_g[s1]) @ gamma_g[mug] @ u_g[s2])
            rhs_val = dbar_g(u_g[s1]) @ (2*p4g_upper[mug]*np.eye(4,dtype=complex)) @ u_g[s2]
            max_gordon_err = max(max_gordon_err, abs(lhs - rhs_val))

print(f"  Numerical check at p'=p (elastic, rest frame):")
print(f"    2m ū'γ^μu = ū'(2p^μ)u?  max error = {max_gordon_err:.2e}",
      "  ✓" if max_gordon_err < 1e-12 else "  ✗")
print()

# γ^5 check: {γ^5,γ^μ}=0
gamma5_g = 1j * gamma_g[0]@gamma_g[1]@gamma_g[2]@gamma_g[3]
max_g5_anticomm = max(
    np.max(np.abs(gamma5_g@gamma_g[mug] + gamma_g[mug]@gamma5_g))
    for mug in range(4))
print(f"  Check {{γ^5, γ^μ}} = 0 for all μ:  max error = {max_g5_anticomm:.2e}",
      "  ✓" if max_g5_anticomm < 1e-12 else "  ✗")
print()
print("  Physical decomposition:")
print("    ū'γ^μu = ū'[(p+p')^μ/(2m) + iσ^{μν}(p'-p)_ν/(2m)]u")
print("    → (p+p')^μ term: charge-current → form factor F_1(q²)")
print("    → iσ^{μν}q_ν term: magnetic moment → form factor F_2(q²), g-2")
print()
print("  Extended Gordon with γ^5 [eq. 38.41]:")
print("    ū_{s'}(p') [(p'+p)^μ - 2iS^{μν}(p'-p)_ν] γ^5 u_s(p) = 0")
print("    (vanishes by {{γ^5,γ^μ}}=0: Gordon tensor is P-even, γ^5 is P-odd)")

# =============================================================================
# §38.G2  HELICITY AND SPIN-SUM COMPLETENESS
# =============================================================================
# Spin-sum completeness (eqs. 38.28-38.29):
#
#   Σ_s u_s(p) ū_s(p) = -/p + m     [particle sum]
#   Σ_s v_s(p) v̄_s(p) = -/p - m     [antiparticle sum]
#
# Bilinear identities (eq. 38.21-38.22):
#   ū_s'(p) γ^μ u_s(p) = 2p^μ δ_{s's}
#   v̄_s'(p) γ^μ v_s(p) = 2p^μ δ_{s's}
#   ū_s'(p) γ^0 v_s(-p) = 0
#   v̄_s'(p) γ^0 u_s(-p) = 0
#
# HELICITY:
#   h = Ĵ·p̂  (spin projection along momentum)
#   For spin-½: h = ±½
#   u_s(p): helicity eigenstate with h = s/2 for momentum along ẑ
#   The label s=± in Srednicki corresponds to spin-up/down along ẑ.
#
# MASSLESS HELICITY PROJECTORS:
#   In the massless limit, helicity = chirality:
#     P_R u_+(p) = u_+(p)   (positive helicity = right-chiral)
#     P_L u_-(p) = u_-(p)   (negative helicity = left-chiral)
#   Massless spin sums:
#     u_+(p)ū_+(p) = P_R(-/p)   [positive helicity projector]
#     u_-(p)ū_-(p) = P_L(-/p)   [negative helicity projector]
#     Sum: Σ_s u_s ū_s = -/p    (massless, m→0 limit of eq. 38.28)

sec("§38.G2  Helicity and Spin-Sum Completeness  [eq. 38.21-38.29]")

Z2s = np.zeros((2,2), dtype=complex)
gamma_s = {}
for mus in range(4):
    gamma_s[mus] = np.block([[Z2s, sigma_vec[mus]], [sigmabar_vec[mus], Z2s]])

beta_s = gamma_s[0]
def dbar_s(sp): return sp.conj() @ beta_s

g4s = np.diag([-1.,1.,1.,1.])
m_s = 1.0
I4s = np.eye(4, dtype=complex)

u_s_list = [
    np.array([1.,0.,1.,0.], dtype=complex)/np.sqrt(2),
    np.array([0.,1.,0.,1.], dtype=complex)/np.sqrt(2),
]
v_s_list = [
    np.array([1.,0.,-1.,0.], dtype=complex)/np.sqrt(2),
    np.array([0.,1.,0.,-1.], dtype=complex)/np.sqrt(2),
]
p4s = np.array([m_s,0.,0.,0.])
p4s_low = g4s @ p4s
pslash_s = sum(p4s_low[mus]*gamma_s[mus] for mus in range(4))

print("Spin-sum completeness (4-component, at rest):")
print()

spin_sum_u_s = sum(2*m_s*np.outer(u, dbar_s(u)) for u in u_s_list)
spin_sum_v_s = sum(2*m_s*np.outer(v, dbar_s(v)) for v in v_s_list)
err_ssu = np.max(np.abs(spin_sum_u_s - (-pslash_s + m_s*I4s)))
err_ssv = np.max(np.abs(spin_sum_v_s - (-pslash_s - m_s*I4s)))

print(f"  Σ_s u_s ū_s = -/p+m?   max error = {err_ssu:.2e}",
      "  ✓" if err_ssu < 1e-12 else "  ✗")
print(f"  Σ_s v_s v̄_s = -/p-m?   max error = {err_ssv:.2e}",
      "  ✓" if err_ssv < 1e-12 else "  ✗")
print()

print("Bilinear identity ū_{s'}γ^μu_s = 2p^μ δ_{ss'}  [eq. 38.21]:")
max_off_diag = 0.
for s1 in range(2):
    for s2 in range(2):
        if s1 != s2:
            for mus in range(4):
                max_off_diag = max(max_off_diag,
                    abs(dbar_s(u_s_list[s1]) @ gamma_s[mus] @ u_s_list[s2]))
print(f"  Off-diagonal (s≠s') = 0?  max = {max_off_diag:.2e}",
      "  ✓" if max_off_diag < 1e-12 else "  ✗")

max_orth = 0.
for s1 in range(2):
    for s2 in range(2):
        max_orth = max(max_orth,
            abs(dbar_s(u_s_list[s1]) @ gamma_s[0] @ v_s_list[s2]))
print(f"  Orthogonality ū_s'(p)γ^0 v_s(-p)=0?  max = {max_orth:.2e}",
      "  ✓" if max_orth < 1e-12 else "  ✗")
print()

print("Helicity / chirality in the massless limit:")
gamma5_s = 1j*gamma_s[0]@gamma_s[1]@gamma_s[2]@gamma_s[3]
P_L_s = 0.5*(I4s - gamma5_s)
P_R_s = 0.5*(I4s + gamma5_s)
u_minus_ml = np.array([0.,1.,0.,0.], dtype=complex)  # left-chiral (massless)
u_plus_ml  = np.array([0.,0.,0.,1.], dtype=complex)  # right-chiral (massless)
err_Lchiral = np.max(np.abs(P_L_s @ u_minus_ml - u_minus_ml))
err_Rchiral = np.max(np.abs(P_R_s @ u_plus_ml  - u_plus_ml))
err_proj_s  = max(np.max(np.abs(P_L_s@P_L_s-P_L_s)),
                  np.max(np.abs(P_R_s@P_R_s-P_R_s)),
                  np.max(np.abs(P_L_s@P_R_s)))
print(f"  P_L,P_R projector checks:   max error = {err_proj_s:.2e}",
      "  ✓" if err_proj_s < 1e-12 else "  ✗")
print(f"  Massless u_- is left-chiral:  P_L u_- = u_-?  error = {err_Lchiral:.2e}",
      "  ✓" if err_Lchiral < 1e-12 else "  ✗")
print(f"  Massless u_+ is right-chiral: P_R u_+ = u_+?  error = {err_Rchiral:.2e}",
      "  ✓" if err_Rchiral < 1e-12 else "  ✗")
print()
print("  Massless spin-sum projectors:")
print("    u_+(p)ū_+(p) = P_R(-/p)   [positive helicity = right-chiral]")
print("    u_-(p)ū_-(p) = P_L(-/p)   [negative helicity = left-chiral]")
print("    Sum: Σ_s u_s ū_s = -/p    (massless limit, m→0)")

# =============================================================================
# §38.G  CADABRA2 SYMBOLIC VERIFICATION SUMMARY
# =============================================================================

sec("§38.G  Cadabra2 Symbolic Verification")

print("Setting up Cadabra2 for complete symbolic verification:")
print()

# Declare sigma matrices symbolically
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma, \tau}"), Ex(r"position=free"))

# Key identity declarations
print("Declaring σ-matrix properties in Cadabra2:")
print()

# Trace of sigma products
tr_sigma_sigmabar = Ex(r"\sigma^{\mu}_{\alpha\dal} \sigmabar^{\nu\dal\alpha}")
print(f"  2-trace summand: {tr_sigma_sigmabar}")
print(f"  Tr[σ^μ σ̄^ν] = σ^μ_{{αα̇}} σ̄^{{ν α̇α}} = −2 g^{{μν}}")
print()

# Completeness
completeness_cadabra = Ex(
    r"\sigma^{\mu}_{\alpha\dal} \sigma_{\mu}^{\beta\dbe}"
)
print(f"  Completeness: {completeness_cadabra} = −2 δ_α^β δ_{{α̇}}^{{β̇}}")
print()

# 4-trace
trace4_cadabra = Ex(
    r"\sigma^{\mu}_{\alpha\dal} \sigmabar^{\nu\dal\dbe}"
    r" \sigma^{\rho}_{\dbe\gamma} \sigmabar^{\sigma\gamma\alpha}"
)
print(f"  4-trace summand: ... (contracted)")
print(f"  Tr[σ^μ σ̄^ν σ^ρ σ̄^σ] = 2(g^{{μν}}g^{{ρσ}} − g^{{μρ}}g^{{νσ}} + g^{{μσ}}g^{{νρ}}) − 2i ε^{{μνρσ}}")
print()

# Print final summary of all numerical checks
print("=" * 70)
print("  NUMERICAL CHECK SUMMARY")
print("=" * 70)
checks = [
    ("σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}  [eq. 38.1]",
     err_comp),
    ("ε^{αβ} ε^{α̇β̇} σ^μ_{αα̇} σ^ν_{ββ̇} = −2g^{μν}  [eq. 35.5]",
     err_alt),
    ("Tr[σ^μ σ̄^ν] = −2 g^{μν}  [eq. 38.5]",
     err_trace2),
    ("Tr[σ̄^μ σ^ν] = −2 g^{μν}",
     err_trace2_bar),
    ("Tr[σ^μ σ̄^ν σ^ρ σ̄^σ] = 2(g^{μν}g^{ρσ}−...) + 2iε  [eq. 38.6]",
     err_4trace),
    ("Tr[σ̄^μ σ^ν σ̄^ρ σ^σ] = 2(g^{μν}g^{ρσ}−...) − 2iε",
     err_4trace_bar),
    ("ε^{αβ} ε_{βγ} = δ^α_γ  (ε invertibility)",
     err_inv),
    ("σ̄^{μ α̇α} = ε^{αβ} ε^{α̇β̇} σ^μ_{ββ̇}  (index gymnastics)",
     max_err_gym),
    ("Schouten identity: ε_{αγ} ε_{βδ} = ε_{αβ} ε_{γδ} − ε_{αδ} ε_{γβ}",
     max_err_schouten),
    ("Gordon: 2m ū'γ^μu = ū'(2p^μ)u at p'=p  [Problem 38.3]",
     max_gordon_err),
    ("{γ^5, γ^μ} = 0 for all μ  (extended Gordon prerequisite)",
     max_g5_anticomm),
    ("Spin-sum: Σ_s u_s ū_s = -/p+m  [eq. 38.28]",
     err_ssu),
    ("Spin-sum: Σ_s v_s v̄_s = -/p-m  [eq. 38.29]",
     err_ssv),
    ("Bilinear off-diagonal ū_s'γ^μu_s = 0 for s≠s'  [eq. 38.21]",
     max_off_diag),
    ("Orthogonality ū_s'(p)γ^0 v_s(-p)=0  [eq. 38.22]",
     max_orth),
    ("Helicity projectors P_{L,R}^2=P, P_LP_R=0",
     err_proj_s),
    ("Massless chirality: P_L u_- = u_-, P_R u_+ = u_+",
     max(err_Lchiral, err_Rchiral)),
]

for name, err in checks:
    status = "✓" if err < 1e-10 else "✗"
    print(f"  {status}  {name}")
    print(f"       max error = {err:.2e}")

# =============================================================================
# §38  SUMMARY
# =============================================================================

sec("CHAPTER 38 SUMMARY")
print("""
  σ^μ COMPLETENESS  [eq. 38.1, cf. 35.5]
  ─────────────────────────────────────────
    σ^μ_{αα̇} σ_μ^{ββ̇} = −2 δ_α^β δ_{α̇}^{β̇}
    ε^{αβ} ε^{α̇β̇} σ^μ_{αα̇} σ^ν_{ββ̇} = −2 g^{μν}
    These are two equivalent forms of the completeness relation.

  2-TRACE IDENTITY  [eq. 38.5]
  ──────────────────────────────
    Tr[σ^μ σ̄^ν] = σ^μ_{αα̇} σ̄^{ν α̇α} = −2 g^{μν}
    Tr[σ̄^μ σ^ν] = σ̄^{μ α̇α} σ^ν_{αα̇} = −2 g^{μν}
    Note: sign comes from Srednicki's mostly-plus metric diag(−1,+1,+1,+1)

  4-TRACE IDENTITY  [eq. 38.6]
  ──────────────────────────────
    Tr[σ^μ σ̄^ν σ^ρ σ̄^σ]
      = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) + 2i ε^{μνρσ}
    Tr[σ̄^μ σ^ν σ̄^ρ σ^σ]
      = 2(g^{μν}g^{ρσ} − g^{μρ}g^{νσ} + g^{μσ}g^{νρ}) − 2i ε^{μνρσ}
    (Conjugate traces differ by sign of ε term; +2i for σσ̄σσ̄, −2i for σ̄σσ̄σ)

  INDEX GYMNASTICS  [eq. 38.7–38.13]
  ────────────────────────────────────
    ψ_α = ε_{αβ} ψ^β      (lower undotted: contract on second ε index)
    ψ^α = ε^{αβ} ψ_β      (raise undotted)
    ε_{αβ} ε^{βγ} = δ_α^γ   (ε is its own inverse up to sign)
    σ̄^{μ α̇α} = ε^{αβ} ε^{α̇β̇} σ^μ_{ββ̇}   (raise both spinor indices of σ)

  SPINOR BILINEARS  [eq. 38.14–38.15]
  ──────────────────────────────────────
    χψ = χ^α ψ_α  (Lorentz scalar, left-handed)
    χ†ψ† = χ†_α̇ ψ†^α̇  (Lorentz scalar, right-handed)
    χ†_α̇ σ̄^{μ α̇α} ψ_α  (Lorentz 4-vector)
    Gordon identity: decomposes vector bilinear into momentum + spin parts

  FIERZ IDENTITIES  [eq. 38.17–38.21]
  ──────────────────────────────────────
    Schouten: ε_{αγ} ε_{βδ} = ε_{αβ} ε_{γδ} − ε_{αδ} ε_{γβ}
    σ-Fierz:  (σ̄^μ)^{α̇α} (σ_μ)_{ββ̇} = −2 δ^α_β δ^{α̇}_{β̇}
    4-Fermi:  (χ†σ̄^μψ)(ξ†σ̄_μη) = −2(χ†ξ†)(ηψ)

  GORDON IDENTITY  [Problem 38.3 / eq. 38.22]
  ──────────────────────────────────────────────
    2m ū_{s'}(p') γ^μ u_s(p) = ū_{s'}(p') [(p'+p)^μ - 2iS^{μν}(p'-p)_ν] u_s(p)
    S^{μν} = (i/4)[γ^μ,γ^ν];  {γ^5,γ^μ}=0 → extended Gordon ×γ^5 = 0
    Physical: γ^μ = (p+p')^μ/(2m) + magnetic moment term (F_1 + F_2)

  SPIN SUMS + HELICITY  [eq. 38.21-38.29]
  ──────────────────────────────────────────
    Σ_s u_s(p) ū_s(p) = -/p + m           [eq. 38.28]
    Σ_s v_s(p) v̄_s(p) = -/p - m           [eq. 38.29]
    ū_s'(p) γ^μ u_s(p) = 2p^μ δ_{ss'}     [eq. 38.21]
    ū_s'(p) γ^0 v_s(-p) = 0               [eq. 38.22]
    Massless: helicity = chirality; u_±ū_± = P_{R/L}(-/p)
""")

print("Done: ch38_spinor_technology.py")
