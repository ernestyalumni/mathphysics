#!/usr/bin/env python3
"""Character coefficients of the SU(2) Wilson plaquette weight are positive (two independent methods).

Class functions of SU(2): g ~ diag(e^{it}, e^{-it}), tr g = 2 cos t; normalized Haar measure on class
functions is (2/pi) sin^2 t dt on [0, pi]; chi_j(t) = sin((2j+1)t)/sin t.
Weight per plaquette with m = 2: exp(-beta (1 - Re tr g / 2)) = e^{-beta} exp(beta cos t).
a_j = <chi_j, exp(beta cos t)> by quadrature, compared with the closed form 2 (2j+1) I_{2j+1}(beta) / beta.
The convolution operator by exp(beta cos t) acts on the isotypic block V_j (x) V_j^* by a_j / (2j+1).
Run with the mathphysics .venv (needs scipy).
"""
import math
from scipy.integrate import quad
from scipy.special import iv

for beta in (0.5, 2.0, 8.0):
    for two_j in range(0, 7):
        n = two_j + 1
        numeric, _ = quad(lambda t: (2 / math.pi) * math.exp(beta * math.cos(t)) * math.sin(n * t) * math.sin(t),
                          0, math.pi, limit=200)
        closed = 2 * n * iv(n, beta) / beta
        assert numeric > 0, (beta, two_j, numeric)
        assert abs(numeric - closed) < 1e-9 * max(1.0, abs(closed)), (beta, two_j, numeric, closed)
print("OK: a_j > 0 and a_j = 2(2j+1) I_{2j+1}(beta)/beta for beta in {0.5, 2, 8}, j = 0..3")
