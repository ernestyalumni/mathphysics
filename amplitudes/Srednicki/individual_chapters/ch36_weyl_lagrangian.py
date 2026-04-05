"""
ch36_weyl_lagrangian.py
========================
Srednicki QFT — Chapter 36: Lagrangians for Spinor Fields
(Weyl Field Lagrangian and Majorana/Dirac Fermions)

Prerequisites: Ch. 22, Ch. 35

What this file covers (section by section):
  §36.A  The Weyl Lagrangian — all terms and hermiticity
  §36.B  Equations of motion — from δS/δψ† = 0
  §36.C  γ matrices in Weyl rep — 4×4 block matrices from σ^μ, σ̄^μ
  §36.D  Verify Clifford algebra {γ^μ, γ^ν} = -2g^{μν} for ALL 10 pairs
  §36.E  Majorana spinor structure — 4-component EOM in matrix form
  §36.F  Dirac fermion from two Weyl fields — the χ, ξ construction

Run with:
    python3 ch36_weyl_lagrangian.py
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

SEP = "=" * 70
sec = lambda s: print(f"\n{SEP}\n  {s}\n{SEP}")

print(SEP)
print("  Srednicki Ch. 36 — Lagrangians for Spinor Fields")
print(SEP)

# ── Cadabra2 index declarations (mirror ch34) ─────────────────────────────
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"),       Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"),        Ex(r"position=free"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# ── numpy setup ───────────────────────────────────────────────────────────
sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),
    2: np.array([[0,-1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0],  [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)

# σ^μ_{aȧ} = (I, σ⃗),  σ̄^{μ ȧa} = (I, -σ⃗)
sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

# Metric g^{μν} = diag(-1,+1,+1,+1)   [Srednicki convention, "mostly plus"]
# Srednicki defines g_{μν} explicitly in eq. (2.4): diag(-1,+1,+1,+1).
# k·p = k_μ p^μ = -k^0 p^0 + k^i p^i  (spatial dot product dominates)
g = np.diag([-1., 1., 1., 1.])

# =============================================================================
# §36.A  THE WEYL LAGRANGIAN
# =============================================================================
# We seek a hermitian, Lorentz-invariant Lagrangian that is quadratic in
# ψ_a and ψ†^ȧ, leading to linear equations of motion.
#
# KINETIC TERM: iψ†σ̄^μ ∂_μ ψ
#   Even though this expression is not individually hermitian, its hermitian
#   conjugate equals itself up to a total derivative (eq. 36.1):
#       (iψ†σ̄^μ ∂_μ ψ)† = iψ†σ̄^μ ∂_μ ψ − i∂_μ(ψ†σ̄^μ ψ)
#   So the kinetic term is real in the action S = ∫d⁴x ℒ.
#
# MASS TERM: -½m ψψ − ½m* ψ†ψ†
#   ψψ = ψ^a ψ_a = ε^{ab} ψ_b ψ_a  (Lorentz scalar, nonzero for Grassmann ψ)
#   The phase of m is unphysical — can always be rotated away by ψ → e^{-iα/2}ψ
#   So take m real and positive: m* = m.
#
# FULL LAGRANGIAN (eq. 36.2):
#   ℒ = iψ†σ̄^μ ∂_μ ψ − ½m ψψ − ½m ψ†ψ†

sec("§36.A  The Weyl Lagrangian")

# Symbolic cadabra2 expressions for all terms
L_kinetic = Ex(r"i \psidag^{\dal} \sigmabar^{\mu}_{\dal\alpha} \partial_{\mu}(\psi^{\alpha})")
L_mass1   = Ex(r"-\frac{1}{2} m \epsilon^{\alpha\beta} \psi_{\alpha} \psi_{\beta}")
L_mass2   = Ex(r"-\frac{1}{2} m \epsilon_{\dal\dbe} \psidag^{\dal} \psidag^{\dbe}")

print("Weyl Lagrangian (eq. 36.2):")
print("  ℒ = ℒ_kinetic + ℒ_mass + ℒ_mass†")
print(f"\n  ℒ_kinetic = {L_kinetic}")
print(f"  ℒ_mass    = {L_mass1}")
print(f"  ℒ_mass†   = {L_mass2}")

print("""
Physical interpretation:
  • Kinetic term: iψ†_{ȧ} σ̄^{μ ȧα} ∂_μ ψ_α
      – Mixes left-handed ψ and right-handed ψ† via σ̄^μ (the "bridge")
      – Real in the action up to boundary terms [eq. 36.1]
  • Mass term: -½m ε^{αβ} ψ_α ψ_β = -½m ψ^α ψ_α
      – Lorentz scalar from antisymmetric ε^{αβ}
      – Nonzero because ψ_α are Grassmann-valued fields (anticommute)
      – Called a "Majorana mass" term
  • For m real: the Lagrangian is fully hermitian (m* = m)
""")

# =============================================================================
# §36.B  EQUATIONS OF MOTION
# =============================================================================
# Varying ℒ with respect to ψ†^ȧ (treating ψ and ψ† as independent):
#
#   0 = -δS/δψ†^ȧ = -iσ̄^μ_{ȧα} ∂_μ ψ^α + m ψ†_ȧ   [eq. 36.3, index-free]
#
# With full spinor indices (eq. 36.4):
#   0 = -i σ̄^{μ ȧc} ∂_μ ψ_c + m ψ†^ȧ
#
# Taking the hermitian conjugate (or varying w.r.t ψ):
#   0 = -i σ^μ_{aċ} ∂_μ ψ†^ċ + m ψ_a             [eq. 36.5 combined]
#
# These two can be combined into a 4-component matrix equation [eq. 36.6]:
#   [ m δ_a^c        -iσ^μ_{aċ}∂_μ  ] [ ψ_c   ] = 0
#   [ -iσ̄^{μȧc}∂_μ  m δ^ȧ_ċ        ] [ ψ†^ċ  ]

sec("§36.B  Equations of Motion")

eom_up_t1   = Ex(r"-i \sigmabar^{\mu\dal\alpha} \partial_{\mu}(\psi_{\alpha})")
eom_up_t2   = Ex(r"m \psidag^{\dal}")
eom_down_t1 = Ex(r"-i \sigma^{\mu}_{\alpha\dal} \partial_{\mu}(\psidag^{\dal})")
eom_down_t2 = Ex(r"m \psi_{\alpha}")

print("EOM from δS/δψ†: (eq. 36.3-36.4)")
print(f"  0 = {eom_up_t1} + {eom_up_t2}")
print("\nEOM from δS/δψ (hermitian conjugate): (eq. 36.5)")
print(f"  0 = {eom_down_t1} + {eom_down_t2}")

print("""
Combined 4-component form (eq. 36.6):
  ⎡  m δ_a^c          -iσ^μ_{aċ}∂_μ  ⎤ ⎡ ψ_c   ⎤   ⎡ 0 ⎤
  ⎢                                    ⎥ ⎢       ⎥ = ⎢   ⎥
  ⎣  -iσ̄^{μȧc}∂_μ    m δ^ȧ_ċ        ⎦ ⎣ ψ†^ċ  ⎦   ⎣ 0 ⎦

Key structural observation:
  • The 2×2 diagonal blocks are proportional to the identity (mass term)
  • The 2×2 off-diagonal blocks mix left/right sectors via σ^μ and σ̄^μ
  • This is the hallmark structure of the Dirac/Majorana equation
""")

# =============================================================================
# §36.C  γ MATRICES IN THE WEYL REPRESENTATION
# =============================================================================
# Define 4×4 gamma matrices via the block structure (eq. 36.7):
#
#   γ^μ = [ 0             σ^μ_{aċ}    ]
#          [ σ̄^{μ ȧc}    0            ]
#
# In terms of numpy arrays with block structure:
#   gamma[mu][0:2, 0:2] = 0₂    (upper-left)
#   gamma[mu][0:2, 2:4] = σ^μ   (upper-right)
#   gamma[mu][2:4, 0:2] = σ̄^μ  (lower-left)
#   gamma[mu][2:4, 2:4] = 0₂    (lower-right)
#
# This is the WEYL (chiral) representation of the γ matrices.
# (Other common representations: Dirac, Majorana — related by unitary transforms)

sec("§36.C  γ Matrices in Weyl Representation")

gamma = {}
for mu in range(4):
    gamma[mu] = np.block([
        [np.zeros((2,2), dtype=complex), sigma_vec[mu]   ],
        [sigmabar_vec[mu],               np.zeros((2,2), dtype=complex)]
    ])

mu_label = {0: '0', 1: '1', 2: '2', 3: '3'}

print("γ^μ = [[0, σ^μ], [σ̄^μ, 0]]  (Weyl/chiral representation, eq. 36.7)")
for mu in range(4):
    print(f"\n  γ^{mu_label[mu]} =")
    for row in gamma[mu]:
        # Format each entry cleanly
        def fmt(z):
            r, i = z.real, z.imag
            if abs(r) < 1e-10 and abs(i) < 1e-10: return " 0 "
            if abs(i) < 1e-10: return f"{int(round(r)):+d} " if abs(r-round(r))<1e-10 else f"{r:+.2f}"
            if abs(r) < 1e-10: return f"{'+i' if i>0 else '-i'}"
            return f"{r:+.1f}{i:+.1f}i"
        print("   [", "  ".join(fmt(z) for z in row), "]")

print("""
The Weyl representation is "chiral" — γ^5 is diagonal (see §36.D).
This makes it natural for massless particles where L and R decouple.
""")

# =============================================================================
# §36.D  VERIFY CLIFFORD ALGEBRA {γ^μ, γ^ν} = -2g^{μν}
# =============================================================================
# The gamma matrices must satisfy the Clifford algebra (eq. 36.9):
#   {γ^μ, γ^ν} = γ^μ γ^ν + γ^ν γ^μ = -2g^{μν} I₄
#
# This follows from the sigma-matrix relations (eq. 36.8):
#   (σ^μ σ̄^ν + σ^ν σ̄^μ)_a^c = -2g^{μν} δ_a^c
#   (σ̄^μ σ^ν + σ̄^ν σ^μ)^ȧ_ċ = -2g^{μν} δ^ȧ_ċ
#
# With γ^μ γ^ν = [[0,σ^μ],[σ̄^μ,0]] · [[0,σ^ν],[σ̄^ν,0]]
#              = [[σ^μσ̄^ν, 0], [0, σ̄^μσ^ν]]
# So {γ^μ,γ^ν} = [[σ^μσ̄^ν+σ^νσ̄^μ, 0], [0, σ̄^μσ^ν+σ̄^νσ^μ]]
#             = [[-2g^{μν}I₂, 0], [0, -2g^{μν}I₂]]
#             = -2g^{μν} I₄  ✓
#
# We verify ALL 10 pairs (μ≤ν) numerically.

sec("§36.D  Verify Clifford Algebra {γ^μ, γ^ν} = -2g^{μν}  [eq. 36.9]")

# With Srednicki's mostly-plus metric g = diag(-1,+1,+1,+1):
#   {γ^0, γ^0} = 2(γ^0)^2 = 2I   and  -2g^{00} = -2(-1) = +2  ✓
#   {γ^i, γ^i} = 2(γ^i)^2 = -2I  and  -2g^{ii} = -2(+1) = -2  ✓
# NOTE: if you (incorrectly) use g = diag(+1,-1,-1,-1), the check fails!
print("Checking {γ^μ, γ^ν} = -2g^{μν} I₄  [Srednicki mostly-plus g = diag(-1,+1,+1,+1)]\n")

max_clifford_err = 0.0
all_ok = True
for mu in range(4):
    for nu in range(mu, 4):
        anticomm = gamma[mu] @ gamma[nu] + gamma[nu] @ gamma[mu]
        expected = -2 * g[mu, nu] * I4
        err = np.max(np.abs(anticomm - expected))
        status = "✓" if err < 1e-12 else "✗"
        print(f"  {status}  {{γ^{mu_label[mu]}, γ^{mu_label[nu]}}} = -2×{g[mu,nu]:+.0f}×I₄  "
              f"  [err={err:.1e}]")
        if err > max_clifford_err:
            max_clifford_err = err
        if err > 1e-12:
            all_ok = False

print(f"\n  Max error: {max_clifford_err:.2e}")
print("  ✓ Clifford algebra verified for all 10 pairs" if all_ok
      else "  ✗ MISMATCH in Clifford algebra!")

# Also compute and verify γ^5
# γ^5 = i γ^0 γ^1 γ^2 γ^3  [eq. 36.46]
# In Weyl rep: γ^5 = diag(-I₂, +I₂)  [eq. 36.43]
print("\n--- γ^5 = i γ^0 γ^1 γ^2 γ^3 ---")
gamma5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]
expected_g5 = np.block([
    [-I2,                    np.zeros((2,2), dtype=complex)],
    [np.zeros((2,2), dtype=complex), I2]
])
err5 = np.max(np.abs(gamma5 - expected_g5))
print(f"\n  γ^5 =")
for row in gamma5.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")
print(f"\n  γ^5 = diag(-I₂, +I₂)?  err = {err5:.2e}")
print("  ✓ γ^5 = diag(-I₂, +I₂) confirmed" if err5 < 1e-12 else "  ✗ MISMATCH!")

# Projection operators P_L = ½(1 - γ^5), P_R = ½(1 + γ^5)  [eq. 36.44]
P_L = 0.5 * (I4 - gamma5)
P_R = 0.5 * (I4 + gamma5)

print("\n--- Projection operators ---")
print("  P_L = ½(I - γ^5) = diag(I₂, 0₂):")
for row in P_L.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")
print("  P_R = ½(I + γ^5) = diag(0₂, I₂):")
for row in P_R.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")

# Verify P_L² = P_L, P_R² = P_R, P_L·P_R = 0, P_L + P_R = I
err_PL  = np.max(np.abs(P_L @ P_L - P_L))
err_PR  = np.max(np.abs(P_R @ P_R - P_R))
err_PLP = np.max(np.abs(P_L @ P_R))
err_sum = np.max(np.abs(P_L + P_R - I4))
print(f"\n  P_L² = P_L:          err = {err_PL:.1e}")
print(f"  P_R² = P_R:          err = {err_PR:.1e}")
print(f"  P_L·P_R = 0:         err = {err_PLP:.1e}")
print(f"  P_L + P_R = I₄:      err = {err_sum:.1e}")
print("  ✓ Projectors verified" if max(err_PL, err_PR, err_PLP, err_sum) < 1e-12
      else "  ✗ Projector error!")

# =============================================================================
# §36.E  MAJORANA SPINOR STRUCTURE
# =============================================================================
# The 4-component Majorana spinor is defined (eq. 36.10):
#   Ψ = ( ψ_c   )
#       ( ψ†^ċ  )
# where both components come from the SAME left-handed Weyl field ψ.
#
# The EOM (eq. 36.6) in terms of γ matrices:
#   (-i γ^μ ∂_μ + m) Ψ = 0    [eq. 36.11, Dirac equation]
#
# To verify structure: the matrix M = m I₄ - i γ^μ p_μ (Fourier space)
# must annihilate any on-shell Majorana spinor.
#
# For a spinor at rest (p^μ = (m,0,0,0)), the EOM becomes:
#   (-i γ^0 ∂_0 + m) Ψ = 0
# In Fourier space: (-i γ^0 (-im) + m) u = (- m γ^0 + m) u = m(I - γ^0) u = 0
# So (I - γ^0) u = 0, meaning u is an eigenstate of γ^0 with eigenvalue 1.

sec("§36.E  Majorana Spinor and 4-Component EOM")

print("Majorana spinor: Ψ = (ψ_c, ψ†^ċ)^T  [eq. 36.10]")
print("Dirac equation:  (-iγ^μ ∂_μ + m)Ψ = 0  [eq. 36.11]")

print("\n--- At-rest Dirac operator (Fourier space, p = (m,0,0,0)) ---")
# In Fourier space: ∂_μ → -i p_μ, p^μ = (E, p⃗) with p^2 = -m²
# EOM: (γ^μ p_μ + m) u = 0  [Srednicki metric: p_μ = g_{μν}p^ν = (-E, p⃗)]
# At rest p^μ = (m, 0, 0, 0): p_μ = (-m, 0, 0, 0)
# Dirac operator: γ^μ p_μ + m = -m γ^0 + m = m(I₄ - γ^0)
m_val = 1.0
p_mu = np.array([-m_val, 0., 0., 0.])   # p_μ (lower index, Srednicki)
dirac_op = sum(gamma[mu] * p_mu[mu] for mu in range(4)) + m_val * I4
print(f"  Dirac op at rest = γ^μ p_μ + m =")
for row in dirac_op.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")

# At rest, solutions satisfy (I - γ^0)u = 0, i.e. γ^0 u = u
# γ^0 has eigenvalues ±1; its +1 eigenspace gives physical solutions
eigenvalues, eigenvectors = np.linalg.eig(gamma[0])
print(f"\n  γ^0 eigenvalues: {eigenvalues.real}")
# Find the +1 eigenstates
pos_eig_idx = [i for i, ev in enumerate(eigenvalues) if abs(ev - 1.0) < 1e-10]
print(f"  +1 eigenstates (physical at rest): indices {pos_eig_idx}")
for idx in pos_eig_idx:
    u = eigenvectors[:, idx]
    residual = np.max(np.abs(dirac_op @ u))
    print(f"    u_{idx} = [{', '.join(f'{v:.3f}' for v in u.real)}]  "
          f"  ||Dirac·u|| = {residual:.1e}")

print("""
Physical interpretation of Majorana spinor:
  • Ψ = (ψ_c, ψ†^ċ)^T has BOTH components from the SAME field ψ
  • This means: the particle IS its own antiparticle (Majorana condition: Ψ^C = Ψ)
  • Analogous to a real scalar field (φ = φ†)
  • Majorana mass term breaks any U(1) charge symmetry
  • Candidate for neutrino masses in physics beyond the Standard Model
""")

# Verify the 4-component EOM matrix structure [eq. 36.6]
print("--- 4-component EOM matrix M = m I₄ - i γ^μ ∂_μ  structure ---")
print("  At rest (Fourier space), M = m(I - γ^0) = m·diag(0,0,2,2) in Weyl rep")
M_rest = m_val * (I4 - gamma[0])
print("  M_rest =")
for row in M_rest.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")

# =============================================================================
# §36.F  DIRAC FERMION FROM TWO WEYL FIELDS
# =============================================================================
# Start with TWO left-handed Weyl fields ψ_1, ψ_2 with SO(2) symmetry [eq.36.12]
# Rewrite in complex basis (eq. 36.14-36.15):
#   χ = (ψ_1 + iψ_2)/√2
#   ξ = (ψ_1 - iψ_2)/√2
#
# The resulting Lagrangian (eq. 36.16):
#   ℒ = iχ†σ̄^μ ∂_μ χ + iξ†σ̄^μ ∂_μ ξ − m χξ − m ξ†χ†
#
# This has U(1) symmetry: χ → e^{-iα}χ, ξ → e^{+iα}ξ  [eq. 36.17]
#
# The 4-component DIRAC spinor (eq. 36.19):
#   Ψ = ( χ_c   )   ← left-handed component (charge -1 under U(1))
#       ( ξ†^ċ  )   ← right-handed component (charge +1 under U(1))
#
# The Dirac Lagrangian (eq. 36.28):
#   ℒ = i Ψ̄ γ^μ ∂_μ Ψ − m Ψ̄Ψ
# where Ψ̄ = Ψ†β with β = γ^0 (numerically)
#
# KEY DISTINCTIONS:
#   Majorana: built from ONE Weyl field, Ψ^C = Ψ, no U(1) symmetry
#   Dirac:    built from TWO Weyl fields, U(1) charge symmetry, particle ≠ antiparticle

sec("§36.F  Dirac Fermion from Two Weyl Fields")

# Cadabra2 expressions for the Dirac Lagrangian terms
L_Dirac_kin_chi = Ex(r"i \chidag^{\dal} \sigmabar^{\mu}_{\dal\alpha} \partial_{\mu}(\chi^{\alpha})")
L_Dirac_kin_xi  = Ex(r"i \xidag^{\dal} \sigmabar^{\mu}_{\dal\alpha} \partial_{\mu}(\xi^{\alpha})")
L_Dirac_mass    = Ex(r"-m \epsilon^{\alpha\beta} \chi_{\alpha} \xi_{\beta}")
L_Dirac_massdg  = Ex(r"-m \epsilon_{\dal\dbe} \xidag^{\dal} \chidag^{\dbe}")

print("Two-Weyl-field Lagrangian (eq. 36.16):")
print(f"  ℒ_kin_χ = {L_Dirac_kin_chi}")
print(f"  ℒ_kin_ξ = {L_Dirac_kin_xi}")
print(f"  ℒ_mass  = {L_Dirac_mass}")
print(f"  ℒ_mass† = {L_Dirac_massdg}")

print("\nDirac spinor Ψ = (χ_c, ξ†^ċ)^T  [eq. 36.19]")
print("Majorana spinor Ψ = (ψ_c, ψ†^ċ)^T  [eq. 36.10]  (same field ψ)")

# Define the beta matrix numerically [eq. 36.21]
# β = [[0, δ^ȧ_ċ], [δ_a^c, 0]] = γ^0 numerically
beta = gamma[0].copy()
print(f"\n  β (= γ^0 numerically) =")
for row in beta.real:
    print("   [", "  ".join(f"{int(round(v)):+d}" for v in row), "]")

# Verify Ψ̄Ψ structure: Ψ̄ = Ψ†β
# For Dirac Ψ = (χ_c, ξ†^ċ)^T:
#   Ψ† = (χ†_ȧ, ξ^a)
#   Ψ̄ = Ψ†β = (ξ^a, χ†_ȧ)   [eq. 36.22]
print("""
Ψ̄Ψ structure (eq. 36.23):
  Ψ  = (χ_c, ξ†^ċ)^T
  Ψ† = (χ†_ȧ, ξ^a)
  Ψ̄  = Ψ†β = (ξ^a, χ†_ȧ)
  Ψ̄Ψ = ξ^a χ_a + χ†_ȧ ξ†^ȧ    (mass term in Dirac Lagrangian)
""")

# Dirac vs Majorana comparison table
print("=" * 60)
print("  DIRAC vs MAJORANA comparison")
print("=" * 60)
print(f"  {'Property':<30} {'Dirac':<15} {'Majorana'}")
print(f"  {'-'*58}")
print(f"  {'Weyl fields needed':<30} {'2 (χ, ξ)':<15} {'1 (ψ)'}")
print(f"  {'U(1) charge symmetry':<30} {'Yes':<15} {'No'}")
print(f"  {'Particle = antiparticle':<30} {'No':<15} {'Yes'}")
print(f"  {'Ψ^C = Ψ (Majorana cond.)':<30} {'No':<15} {'Yes'}")
print(f"  {'Lagrangian':<30} {'iΨ̄γ∂Ψ - mΨ̄Ψ':<15} {'(i/2)Ψ̄γ∂Ψ - (m/2)Ψ̄Ψ'}")
print(f"  {'Scalar field analogue':<30} {'complex φ':<15} {'real φ'}")
print(f"  {'Neutrino candidate':<30} {'(Dirac ν)':<15} {'Yes (seesaw)'}")
print("=" * 60)

# Charge conjugation matrix C [eq. 36.32]
# C = [[ε_{ac}, 0], [0, ε^{ȧċ}]]
eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)   # ε_{12}=-1, ε_{21}=+1
eps_upper = np.array([[0,  1], [-1, 0]], dtype=complex)   # ε^{12}=+1, ε^{21}=-1
C_matrix = np.block([
    [eps_lower, np.zeros((2,2), dtype=complex)],
    [np.zeros((2,2), dtype=complex), eps_upper]
])

print("\n--- Charge conjugation matrix C [eq. 36.32] ---")
print("  C =")
for row in C_matrix:
    def fmt(z):
        r, i = z.real, z.imag
        if abs(r) < 1e-10 and abs(i) < 1e-10: return " 0"
        if abs(i) < 1e-10: return f"{int(round(r)):+d}"
        return f"{r:+.0f}{i:+.0f}i"
    print("   [", "  ".join(fmt(z) for z in row), "]")

# Verify C properties (eq. 36.36):
# C^T = C† = C^{-1} = -C
CT  = C_matrix.T
Cdg = C_matrix.conj().T
Cinv = np.linalg.inv(C_matrix)
err_CT  = np.max(np.abs(CT  - C_matrix))          # C^T = C  or -C?
err_Cdg = np.max(np.abs(Cdg - Cinv))               # C† = C^{-1}
err_neg = np.max(np.abs(C_matrix + CT))             # C^T = -C

print(f"\n  C^T = -C:   err = {err_neg:.1e}")
print(f"  C† = C^{{-1}}: err = {err_Cdg:.1e}")
C_inv_gamma_C = Cinv @ gamma[0] @ C_matrix
print(f"\n  C^{{-1}}γ^0 C = -γ^0^T?  "
      f"err = {np.max(np.abs(C_inv_gamma_C + gamma[0].T)):.1e}")

# Verify C^{-1} γ^μ C = -(γ^μ)^T for all μ  [eq. 36.40]
print("\n--- Verify C^{-1}γ^μ C = -(γ^μ)^T  [eq. 36.40] ---")
max_cc_err = 0.0
for mu in range(4):
    lhs = Cinv @ gamma[mu] @ C_matrix
    rhs = -gamma[mu].T
    err = np.max(np.abs(lhs - rhs))
    status = "✓" if err < 1e-12 else "✗"
    print(f"  {status}  C^{{-1}}γ^{mu_label[mu]}C = -(γ^{mu_label[mu]})^T  [err={err:.1e}]")
    max_cc_err = max(max_cc_err, err)

print(f"\n  Max error: {max_cc_err:.2e}")
print("  ✓ C^{-1}γ^μ C = -(γ^μ)^T verified for all μ" if max_cc_err < 1e-12
      else "  ✗ MISMATCH!")

# =============================================================================
# §36  SUMMARY
# =============================================================================

sec("CHAPTER 36 SUMMARY")
print("""
  WEYL LAGRANGIAN (eq. 36.2)
  ───────────────────────────
    ℒ = iψ†σ̄^μ ∂_μ ψ − ½m ψψ − ½m ψ†ψ†
    • Kinetic term real in action (total derivative cancels)
    • Majorana mass m can always be taken real and positive

  EQUATIONS OF MOTION (eq. 36.3-36.4)
  ────────────────────────────────────
    −iσ̄^{μȧc}∂_μ ψ_c + m ψ†^ȧ = 0
    −iσ^μ_{aċ}∂_μ ψ†^ċ + m ψ_a = 0

  γ MATRICES — WEYL REPRESENTATION (eq. 36.7)
  ─────────────────────────────────────────────
    γ^μ = [[0, σ^μ], [σ̄^μ, 0]]   (4×4 block matrices)
    {γ^μ, γ^ν} = −2g^{μν} I₄     (Clifford algebra, eq. 36.9;  g = diag(−1,+1,+1,+1))
    γ^5 = iγ^0γ^1γ^2γ^3 = diag(−I₂, +I₂)

  MAJORANA SPINOR (eq. 36.10-36.11)
  ────────────────────────────────────
    Ψ = (ψ_c, ψ†^ċ)^T    [one Weyl field]
    (−iγ^μ∂_μ + m)Ψ = 0  [Dirac equation]
    Ψ^C = Ψ              [Majorana condition: particle = antiparticle]

  DIRAC SPINOR (eq. 36.19)
  ─────────────────────────
    Ψ = (χ_c, ξ†^ċ)^T    [two Weyl fields]
    ℒ = iΨ̄γ^μ∂_μΨ − mΨ̄Ψ
    U(1): χ → e^{-iα}χ, ξ → e^{+iα}ξ  (conserved charge = particle number)

  CHARGE CONJUGATION (eq. 36.40)
  ────────────────────────────────
    C^{-1}γ^μ C = −(γ^μ)^T
    Majorana: Ψ^C = Ψ  (self-conjugate, no U(1) charge)
    Dirac: Ψ^C ≠ Ψ    (distinct particle and antiparticle)
""")

print("Done: ch36_weyl_lagrangian.py")
