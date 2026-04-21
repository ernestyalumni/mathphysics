from mhvamplitudes.spinors.GammaMatrices import GammaMatrices
from mhvamplitudes.spinors import massless_spinors
# =============================================================================
# §60.B  EXPLICIT SPINOR CONSTRUCTION (WEYL REPRESENTATION)
# =============================================================================
# For massless momentum k^μ = (ω, 0, 0, ω) along the +z axis, the
# explicit Dirac spinors are given by Srednicki eq. 60.10:
#
#   |k]  = u_-(k) = √(2ω) (0, 1, 0, 0)^T
#   |k⟩  = u_+(k) = √(2ω) (0, 0, 1, 0)^T
#
# These satisfy the massless Dirac equation:  γ^μ k_μ |k] = k/ |k] = 0
# and the helicity eigenvalue equation:       h |k] = -½ |k],  h |k⟩ = +½ |k⟩
#
# For the bra spinors (row vectors):
#   [k|  = ū_+(k) = √(2ω) (0, 1, 0, 0)^† = √(2ω) (0, 1, 0, 0) · γ^0
#   ⟨k|  = ū_-(k) = √(2ω) (0, 0, 0, 1) · ... (see below)
#
# The bra-ket relation in Srednicki's convention:
#   [k| = (|k])† γ^0   and   ⟨k| = (|k⟩)† γ^0

def test_verify_massless_Dirac_equation():

    # arbitrary massless energy
    omega = 2.5
    ket_sq_k, ket_an_k, bra_sq_k, bra_an_k = massless_spinors(omega)

    print(f"Massless momentum k^μ = ({omega}, 0, 0, {omega})  [along +z axis]")
    print(f"  |k] = u_-(k) = {ket_sq_k}")
    print(f"  |k⟩ = u_+(k) = {ket_an_k}")
    print(f"  [k| = ū_+(k) = {bra_sq_k}")
    print(f"  ⟨k| = ū_-(k) = {bra_an_k}")

    # Verify: massless Dirac equation k/ |k] = 0
    k_mu = np.array([omega, 0., 0., omega])

    gamma = GammaMatrices()
    gam = gamma.gam()
    g = gamma.metric

    slash_k = sum(g[mu, mu] * k_mu[mu] * gam[mu] for mu in range(4))
    # NOTE: with Srednicki's (-,+,+,+) metric, p/ = k^μ γ_μ = Σ_μ g_{μμ} k^μ γ^μ
    #       = -k^0 γ^0 + k^i γ^i  (time component picks up minus sign)
    # Lowering: p_μ = g_{μν} p^ν → p_0 = +p^0, p_i = -p^i
    # So p/ = γ^0 p_0 + γ^i p_i = γ^0 p^0 - γ^i p^i (Einstein convention)
    # Direct: p/ = γ^μ p_μ with p_μ = (ω, 0, 0, -ω) for k^μ = (ω,0,0,ω)
    k_lower = np.array([omega, 0., 0., -omega])   # k_μ = g_{μν} k^ν
    slash_k = sum(k_lower[mu] * gam[mu] for mu in range(4))

    dirac_sq = slash_k @ ket_sq_k
    dirac_an = slash_k @ ket_an_k
    print(f"\nDirac equation check: k/ |k] = {np.max(np.abs(dirac_sq)):.2e}  (should be 0)")
    print(f"Dirac equation check: k/ |k⟩ = {np.max(np.abs(dirac_an)):.2e}  (should be 0)")
    assert np.allclose(dirac_sq, 0), "k/ |k] ≠ 0!"
    assert np.allclose(dirac_an, 0), "k/ |k⟩ ≠ 0!"
    print("  ✓ Both spinors satisfy massless Dirac equation k/ |ψ⟩ = 0")

    # Verify: k^2 = 0
    k2 = k_mu @ g @ k_mu
    print(f"\nMasslessness check: k^2 = k^μ k_μ = {k2:.4f}  (should be 0)")
    assert abs(k2) < 1e-14, "k^2 ≠ 0!"
    print("  ✓ k^2 = 0 verified")

def test_helicity_check():
    omega = 2.5
    ket_sq_k, ket_an_k, _, _ = massless_spinors(omega)

    # Helicity operator h = ½ Σ·p̂  (Σ = spin matrix, p̂ = unit 3-momentum)
    # For k along z:  h = ½ γ^1 γ^2 γ^3 γ^0 ... or more simply use γ^5 γ^3
    # Helicity = ½ (momentum · σ) / |momentum| acting on upper/lower components
    # In Weyl rep, γ^5 = diag(-I_2, +I_2), so helicity eigenstates are:
    #   positive h = +½: lower components nonzero → |k⟩ = u_+(k)
    #   negative h = -½: upper components nonzero → |k]  = u_-(k)
    print("\nHelicity check (for +z momentum): lower 2 components = + helicity:")
    print(f"  |k] upper (a=0,1): {ket_sq_k[:2]}  (should be 0 for negative helicity)")
    print(f"  |k] lower (a=2,3): {ket_sq_k[2:]}  (upper 2 in Weyl block)")
    print(f"  |k⟩ lower (a=2,3): {ket_an_k[2:]}  (should be 0 for positive helicity)")
    print(f"  |k⟩ upper (a=0,1): {ket_an_k[:2]}  (lower 2 in Weyl block)")

    assert np.allclose(ket_sq_k[:2], 0), "k/ |k] ≠ 0!"
    assert np.allclose(ket_sq_k[2:], 0), "k/ |k] ≠ 0!"
    assert np.allclose(ket_an_k[2:], 0), "k/ |k⟩ ≠ 0!"
    assert np.allclose(ket_an_k[:2], 0), "k/ |k⟩ ≠ 0!"
