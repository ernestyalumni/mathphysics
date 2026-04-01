# Patch to replace pol_vector_plus/minus in ch60_spinor_helicity.py
# This uses 2-component spinor formalism per Srednicki eq. 60.13-60.16

replacement_text = '''sec("§60.D Polarization vectors in spinor-helicity formalism")

def pol_vector_plus(lam_q, lamtil_q, lam_k, lamtil_k):
    """
    Compute ε_+^μ(k; q) = ⟨q|σ^μ|k] / (√2 ⟨qk⟩) using 2-component spinors.
    
    Per Srednicki eq. 60.13: ε_+^μ = ⟨q|σ^μ|k] / (√2 ⟨qk⟩)
    where ⟨q|σ^μ|k] = λ̃_q^* · σ^μ · λ_k (with σ^μ = (I, σ_x, σ_y, σ_z))
    """
    ang_qk = angle_bracket(lam_q, lam_k)
    denom = np.sqrt(2) * ang_qk
    eps_plus = np.zeros(4, dtype=complex)
    # σ^μ = (I, σ_x, σ_y, σ_z)
    sigma_list = [I2, sigma[1], sigma[2], sigma[3]]
    for mu in range(4):
        # ⟨q|σ^μ|k] = lamtil_q.conj() · σ^μ · lam_k
        val = lamtil_q.conj() @ sigma_list[mu] @ lam_k
        eps_plus[mu] = val / denom
    return eps_plus

def pol_vector_minus(lam_q, lamtil_q, lam_k, lamtil_k):
    """
    Compute ε_-^μ(k; q) = [q|σ̄^μ|k⟩ / (√2 [qk]) using 2-component spinors.
    
    Per Srednicki eq. 60.14: ε_-^μ = [q|σ̄^μ|k⟩ / (√2 [qk])
    where σ̄^μ = (I, -σ_x, -σ_y, -σ_z)
    """
    sq_qk = square_bracket(lamtil_q, lamtil_k)
    denom = np.sqrt(2) * sq_qk
    eps_minus = np.zeros(4, dtype=complex)
    sigmabar_list = [I2, -sigma[1], -sigma[2], -sigma[3]]
    for mu in range(4):
        # [q|σ̄^μ|k⟩ = lam_q.conj() · σ̄^μ · lamtil_k
        val = lam_q.conj() @ sigmabar_list[mu] @ lamtil_k
        eps_minus[mu] = val / denom
    return eps_minus
'''

print(replacement_text)
