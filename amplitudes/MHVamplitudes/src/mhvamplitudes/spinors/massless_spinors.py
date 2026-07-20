import numpy as np

from .GammaMatrices import GammaMatrices


def massless_spinors(omega: float):
    """
    Returns |k], |k⟩ (4-component Dirac spinors in Weyl rep)
    for momentum k^μ = (ω, 0, 0, ω) along the +z axis.
    Also returns [k| = (|k])†γ^0 and ⟨k| = (|k⟩)†γ^0.
    """
    gamma_zero = GammaMatrices().gamma_mu(0)
    sq2om = np.sqrt(2 * omega)
    # From Srednicki eq. 60.10:
    ket_sq = sq2om * np.array([0, 1, 0, 0], dtype=complex)   # |k]
    ket_an = sq2om * np.array([0, 0, 1, 0], dtype=complex)   # |k⟩
    # Bra spinors: ū = ψ† γ^0
    bra_sq = ket_sq.conj() @ gamma_zero    # [k| = (|k])† γ^0
    bra_an = ket_an.conj() @ gamma_zero    # ⟨k| = (|k⟩)† γ^0
    return ket_sq, ket_an, bra_sq, bra_an
