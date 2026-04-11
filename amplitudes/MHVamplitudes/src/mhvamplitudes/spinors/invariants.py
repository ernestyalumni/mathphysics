"""
spinors.invariants — ε tensor and index raising/lowering.
Follows Srednicki’s sign convention exactly (eq. 34.22).
"""

import numpy as np


class EpsilonTensor:
    """
    SL(2,C) invariant antisymmetric tensor ε_{ab} and ε^{ab}.
    Srednicki convention: ε_{12} = -1, ε_{21} = +1, ε^{12} = +1, ε^{21} = -1.
    """

    def __init__(self):
        self.eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
        self.eps_upper = np.array([[0, 1], [-1, 0]], dtype=complex)

    def raise_index(self, psi_lower: np.ndarray) -> np.ndarray:
        """ψ^a = ε^{ab} ψ_b"""
        return self.eps_upper @ psi_lower

    def lower_index(self, psi_upper: np.ndarray) -> np.ndarray:
        """ψ_a = ε_{ab} ψ^b"""
        return self.eps_lower @ psi_upper

    def contract_left(self, psi, chi) -> complex:
        """ψ^a χ_a = ε^{ab} ψ_b χ_a (includes the minus sign relative to ψ_a χ^a)"""
        return np.dot(self.raise_index(psi), chi)

    def verify(self):
        """Check ε_{ab} ε^{bc} = δ_a^c and invariance."""
        product = self.eps_lower @ self.eps_upper
        assert np.allclose(product, np.eye(2)), "ε normalization failed"
        print("✓ ε_{ab} ε^{bc} = δ_a^c verified")

        # Simple invariance test under SL(2,C) rotation (matches original
        # script)
        theta = 0.7
        L = (np.cos(theta/2) * np.eye(2, dtype=complex) +
             1j * np.sin(theta/2) * np.array([[1, 0], [0, -1]], dtype=complex))
        # Correct transformation for left-handed: L^T ε L = ε (with Srednicki
        # sign convention)
        check = L.T @ self.eps_lower @ L
        assert np.allclose(check, self.eps_lower, atol=1e-14), "ε not invariant"
        print("✓ ε_{ab} is SL(2,C) invariant")
        print(f"  det(L) = {np.linalg.det(L):.6f}")
        return True
