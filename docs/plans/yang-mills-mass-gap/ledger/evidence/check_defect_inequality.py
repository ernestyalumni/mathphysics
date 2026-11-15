"""Numerical test of claim C9: for eigenvalues l0 >= l1 >= ... >= 0 of a positive operator,
D(n) = (sum l_i^n)^2 / sum l_i^(2n) - 1 >= sum_{i>=1} (l_i/l0)^n.
Random spectra, including near-degenerate and heavy-tailed ones. Exact rational arithmetic."""
import random
from fractions import Fraction

random.seed(20261114)
worst = None
violations = 0
trials = 0
for kind in ("uniform", "near_degenerate", "geometric", "many_small"):
    for _ in range(2000):
        size = random.randint(2, 40)
        if kind == "uniform":
            spectrum = [Fraction(random.randint(0, 10**6), 10**6) for _ in range(size)]
        elif kind == "near_degenerate":
            spectrum = [Fraction(10**6 - random.randint(0, 50), 10**6) for _ in range(size)]
        elif kind == "geometric":
            ratio = Fraction(random.randint(1, 99), 100)
            spectrum = [ratio**k for k in range(size)]
        else:
            spectrum = [Fraction(1)] + [Fraction(random.randint(0, 1000), 10**6) for _ in range(size - 1)]
        spectrum.sort(reverse=True)
        if spectrum[0] == 0:
            continue
        for power in (1, 2, 3, 7):
            trace_n = sum(value**power for value in spectrum)
            trace_2n = sum(value ** (2 * power) for value in spectrum)
            defect = trace_n**2 / trace_2n - 1
            excited = sum((value / spectrum[0]) ** power for value in spectrum[1:])
            trials += 1
            slack = defect - excited
            if slack < 0:
                violations += 1
            if worst is None or slack < worst:
                worst = slack
print("trials", trials, "violations", violations, "smallest slack", float(worst))
