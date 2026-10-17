"""theta-logistics (TheoPhysEval, Stochastic processes and biological physics).

    dX = X [1 - (X/K)^theta] dt + sqrt(sigma) X dW     (Ito)

Stationary Fokker-Planck density, with drift a(x) and diffusion b(x)^2 = sigma x^2:

    p(x) ∝ b^-2 exp( ∫ 2a/b^2 dx ) ∝ x^(2/sigma - 2) exp( -(2/(sigma theta)) (x/K)^theta ),

normalizable at 0 iff 2/sigma - 2 > -1, i.e. sigma < 2. Put a = 2/sigma - 1 > 0. Then

    E[X^m] = (K c)^m Gamma((a+m)/theta) / Gamma(a/theta),   c = (sigma theta / 2)^(1/theta)

(substitute y = (2/(sigma theta)) (x/K)^theta in the Gamma integral).

theta = 1: Gamma(shape a, scale sigma K / 2). Mean 1 and variance 1 force shape 1, so
           sigma = 1, K = 2, and X ~ Exp(1).
theta = 2: mean 1 and E[X^2] = 2 fix a (one equation after eliminating K) and then K.

Detection with Poisson sampling: P_detect(V) = 1 - E[exp(-V X)]; V_theta solves
E[exp(-V X)] = 0.05. For Exp(1), E[e^{-VX}] = 1/(1+V), so V_1 = 19 exactly.

Run: python theta_logistic.py   (analytic answer + Euler-Maruyama Monte Carlo check)
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammaln


def moment(m: float, a: float, K: float, sigma: float, theta: float) -> float:
    c = (sigma * theta / 2) ** (1 / theta)
    return (K * c) ** m * np.exp(gammaln((a + m) / theta) - gammaln(a / theta))


def calibrate(theta: float):
    """Find (K, sigma) with stationary mean 1 and variance 1."""
    # E[X^2]/E[X]^2 = Gamma((a+2)/t) Gamma(a/t) / Gamma((a+1)/t)^2 must equal 2; K drops out.
    def ratio(a):
        return np.exp(gammaln((a + 2) / theta) + gammaln(a / theta) - 2 * gammaln((a + 1) / theta)) - 2.0
    a = brentq(ratio, 1e-6, 1e3, xtol=1e-15)
    sigma = 2.0 / (a + 1.0)
    K = 1.0 / moment(1, a, 1.0, sigma, theta)
    return K, sigma, a


def density(x, K, sigma, theta, a):
    c = (sigma * theta / 2) ** (1 / theta)
    y = x / (K * c)
    log_norm = np.log(theta) - np.log(K * c) - gammaln(a / theta)
    return np.exp(log_norm + (a - 1) * np.log(y) - y ** theta)


def laplace(V, K, sigma, theta, a):
    """E[exp(-V X)]. In z = (x/(K c))^theta the stationary law is Gamma(s), s = a/theta:
        (1/Gamma(s)) ∫ z^(s-1) e^-z exp(-V K c z^(1/theta)) dz.
    The z^(s-1) endpoint singularity (s ~ 0.29 for theta = 2) defeats generic quadrature
    (mpmath at 15-20 digits and scipy quad both drift at 1e-4). Substituting z = w^(1/s)
    gives z^(s-1) dz = dw/s and a smooth integrand:
        (1/Gamma(s+1)) ∫ exp(-w^(1/s) - V K c w^(1/(s theta))) dw."""
    c = (sigma * theta / 2) ** (1 / theta)
    s = a / theta
    f = lambda w: np.exp(-w ** (1 / s) - V * K * c * w ** (1 / (s * theta)))
    val, _ = quad(f, 0, np.inf, epsabs=1e-15, epsrel=1e-13, limit=500)
    return val / np.exp(gammaln(s + 1))


def detection_volume(theta: float, target: float = 0.95) -> float:
    K, sigma, a = calibrate(theta)
    return brentq(lambda V: 1 - laplace(V, K, sigma, theta, a) - target, 1e-6, 1e6, xtol=1e-13)


def monte_carlo(theta: float, K: float, sigma: float, n_paths: int = 20000, t_end: float = 60.0,
                dt: float = 2e-3, seed: int = 0):
    """Log-space Euler-Maruyama (keeps X > 0); returns stationary samples at t_end."""
    rng = np.random.default_rng(seed)
    y = np.zeros(n_paths)                     # y = ln X, start at X = 1
    for _ in range(int(t_end / dt)):
        x = np.exp(y)
        y += (1 - (x / K) ** theta - sigma / 2) * dt + np.sqrt(sigma * dt) * rng.standard_normal(n_paths)
    return np.exp(y)


if __name__ == "__main__":
    results = {}
    for theta in (1.0, 2.0):
        K, sigma, a = calibrate(theta)
        m1, m2 = moment(1, a, K, sigma, theta), moment(2, a, K, sigma, theta)
        V = detection_volume(theta)
        results[theta] = V
        print(f"theta={theta:g}: K={K:.12f} sigma={sigma:.12f} (a={a:.12f}); mean={m1:.12f} var={m2 - m1**2:.12f}; "
              f"V={V:.10f}")
    print(f"V2/V1 = {results[2.0] / results[1.0]:.10f}  ->  {results[2.0] / results[1.0]:.1f}")

    print("\nMonte Carlo check (Euler-Maruyama in ln X, 20k paths to t=60):")
    for theta in (1.0, 2.0):
        K, sigma, _ = calibrate(theta)
        xs = monte_carlo(theta, K, sigma)
        p_detect = 1 - np.mean(np.exp(-results[theta] * xs))
        print(f"theta={theta:g}: sample mean={xs.mean():.3f} var={xs.var():.3f} "
              f"P_detect(V_theta)={p_detect:.4f} (target 0.95)")
