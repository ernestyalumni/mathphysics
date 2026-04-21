from .LorentzMetric import LorentzMetric
from .SigmaMatrices import SigmaMatrices

import numpy as np

class GammaMatrices:
    Z2 = np.zeros((2, 2), dtype=complex)
    def __init__(self, sigma: SigmaMatrices):
        self.sigma = SigmaMatrices()
        self.metric = LorentzMetric()
    # 4×4 γ matrices in Weyl representation: γ^μ =
    # [[  0    ,  σ^μ  ], [ σ̄^μ  ,   0   ]]
    def gamma_mu(self, mu: int):
        return np.block(
            [
                [self.Z2, self.sigma.sigma_mu(mu)],
                [self.sigma.sigmabar_mu(mu), self.Z2]])

    def gamma_bar_mu(self, mu: int):
        return self.sigma.sigmabar_mu(mu)

    def gam(self):
        return [self.gamma_mu(mu) for mu in range(4)]

    def gamma5(self):
        return 1j * self.gam[0] @ self.gam[1] @ self.gam[2] @ self.gam[3]

    def gam_mu_lower(self, mu: int):
        return self.gamma_mu(mu) @ self.metric[(mu, mu)]