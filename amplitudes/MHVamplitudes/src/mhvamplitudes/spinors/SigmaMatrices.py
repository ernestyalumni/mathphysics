"""
Srednicki.spinors.sigma — σ^μ_{aȧ} and σ̄^μ_{ȧa} matrices.

Implements the spinor-vector dictionary and key identities.
"""

import numpy as np
from .pauli_matrices import I2, pauli_matrices

class SigmaMatrices:
    """σ^μ and σ̄^μ in the Weyl basis."""

    def __init__(self):
        self.sigma = {
            0: I2,
            1: pauli_matrices[1],
            2: pauli_matrices[2],
            3: pauli_matrices[3],
        }
        self.sigmabar = {
            0: I2,
            1: -pauli_matrices[1],
            2: -pauli_matrices[2],
            3: -pauli_matrices[3],
        }

    def sigma_mu(self, mu: int):
        return self.sigma[mu]

    def sigmabar_mu(self, mu: int):
        return self.sigmabar[mu]

    def trace_sigma_sigmabar(self, mu: int, nu: int) -> complex:
        """tr(σ^μ σ̄^ν) = -2 g^{μν} in Srednicki mostly-plus metric."""
        return np.trace(self.sigma[mu] @ self.sigmabar[nu])

    def check_trace_identity(self, metric) -> bool:
        """Verify tr(σ^μ σ̄^ν) = -2 g^{μν}."""
        expected = -2.0 * metric.g
        max_err = 0.0
        for mu in range(4):
            for nu in range(4):
                val = self.trace_sigma_sigmabar(mu, nu)
                err = abs(val - expected[mu, nu])
                if err > max_err:
                    max_err = err
        print(f"  Max error in tr(σ^μ σ̄^ν) = -2g^{{μν}}: {max_err:.2e}")
        return max_err < 1e-12

    def print_sigma(self):
        print("σ^μ_{aȧ} = (I, σ⃗):")
        for mu in range(4):
            print(f"  σ^{mu} =")
            print("   ", self.sigma[mu][0])
            print("   ", self.sigma[mu][1])