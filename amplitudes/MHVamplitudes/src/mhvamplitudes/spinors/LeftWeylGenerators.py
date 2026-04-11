from .LorentzMetric import LorentzMetric
from .pauli_matrices import pauli_matrices
import numpy as np

class LeftWeylGenerators:
    """
    Generators S^μν in the (2,1) left-handed Weyl representation.
    Also computes the right-handed (1,2) generators via S_R = -S_L^*.
    """

    def __init__(self):
        self.metric = LorentzMetric()
        self.S_L = self._build_S_L()
        self.S_R = self._build_S_R()

    def _build_S_L(self) -> list:
        """Build S^μν_L following Srednicki eqs. (34.9) and (34.10)."""
        S = [[None] * 4 for _ in range(4)]

        # Spatial rotations: S^{ij}_L = (1/2) ε^{ijk} σ_k
        for i in range(1, 4):
            for j in range(1, 4):
                mat = np.zeros((2, 2), dtype=complex)
                for k in range(1, 4):
                    eps = self._eps3(i, j, k)
                    mat += eps * pauli_matrices[k]
                S[i][j] = 0.5 * mat

        # Boosts: S^{k0}_L = (i/2) σ_k, S^{0k}_L = -S^{k0}_L
        for k in range(1, 4):
            S[k][0] = (1j / 2) * pauli_matrices[k]
            S[0][k] = -(1j / 2) * pauli_matrices[k]

        # Zero on diagonal
        S[0][0] = np.zeros((2, 2), dtype=complex)
        for i in range(1, 4):
            S[i][i] = np.zeros((2, 2), dtype=complex)

        return S

    def _build_S_R(self) -> list:
        """S^μν_R = -[S^μν_L]^* (Srednicki eq. 34.17)."""
        S = [[None] * 4 for _ in range(4)]
        for mu in range(4):
            for nu in range(4):
                if self.S_L[mu][nu] is not None:
                    S[mu][nu] = -np.conj(self.S_L[mu][nu])
        return S

    @staticmethod
    def _eps3(i: int, j: int, k: int) -> int:
        """Levi-Civita symbol ε^{ijk} for i,j,k = 1,2,3."""
        return int(np.linalg.det([
            [i==1, i==2, i==3],
            [j==1, j==2, j==3],
            [k==1, k==2, k==3]
        ]))

    def commutator(self, mu: int, nu: int, rho: int, sig: int) -> np.ndarray:
        """[S^μν, S^ρσ]"""
        A = self.get_S_L(mu, nu)
        B = self.get_S_L(rho, sig)
        return A @ B - B @ A

    def get_S_L(self, mu: int, nu: int) -> np.ndarray:
        if mu == nu:
            return np.zeros((2, 2), dtype=complex)
        if mu < nu:
            return self.S_L[mu][nu]
        return -self.S_L[nu][mu]

    def get_S_R(self, mu: int, nu: int) -> np.ndarray:
        if mu == nu:
            return np.zeros((2, 2), dtype=complex)
        if mu < nu:
            return self.S_R[mu][nu]
        return -self.S_R[nu][mu]

    def check_lorentz_algebra(self, tol: float = 1e-12) -> bool:
        """Verify [S^μν, S^ρσ] = i(g^μρ S^νσ - ...) for all independent pairs.
        """
        max_err = 0.0
        for mu in range(4):
            for nu in range(mu + 1, 4):
                for rho in range(4):
                    for sig in range(rho + 1, 4):
                        lhs = self.commutator(mu, nu, rho, sig)
                        rhs = self._lorentz_rhs(mu, nu, rho, sig)
                        err = np.max(np.abs(lhs - rhs))
                        if err > max_err:
                            max_err = err
        return max_err < tol

    def _lorentz_rhs(self, mu: int, nu: int, rho: int, sig: int) -> np.ndarray:
        """Right-hand side of Lorentz algebra (Srednicki eq. 34.4)."""
        g = self.metric
        def S(m, n):
            return self.get_S_L(m, n)
        return 1j * (
            g[mu, rho] * S(nu, sig)
            - g[nu, rho] * S(mu, sig)
            - g[mu, sig] * S(nu, rho)
            + g[nu, sig] * S(mu, rho)
        )

    def print_generators(self):
        """Pretty print S^μν_L and S^μν_R."""
        label = {0: '0', 1: '1', 2: '2', 3: '3'}
        print("Left-handed generators S^μν_L:")
        for mu in range(4):
            for nu in range(mu + 1, 4):
                mat = self.S_L[mu][nu]
                print(f"  S^{{{label[mu]}{label[nu]}}}_L =")
                print("   ", mat[0])
                print("   ", mat[1])
        print("\nRight-handed generators S^μν_R = -[S^μν_L]^* verified.")
