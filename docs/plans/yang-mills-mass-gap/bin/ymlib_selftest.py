"""Checks for ymlib (YM-T18). Run in the ym-numcheck sandbox: `ym-numcheck selftest`.

Exact identities are checked to 1e-10. Haar moments are checked by Monte Carlo against known values,
with a tolerance of 6 standard errors estimated from the same sample:
  U(N):  E|tr U|^2 = 1,  E|tr U^2|^2 = min(2, N),  E|U_11|^2 = 1/N,  E (tr U)^N = 0.
  SU(N): E|tr U|^2 = 1 (N >= 2),  E (tr U)^N = 1  (one singlet in the N-fold tensor power).
The last pair separates Haar measure on SU(N) from Haar measure on U(N).
"""
import sys
import numpy as np

sys.path.insert(0, "/w")
import ymlib as ym  # noqa: E402

FAILS = []


def check(name, ok):
    print(("ok   " if ok else "FAIL ") + name)
    if not ok:
        FAILS.append(name)


def mc(name, samples, expected):
    x = np.asarray(samples)
    se = max(np.std(x) / np.sqrt(len(x)), 1e-12)
    m = np.mean(x)
    check("%s = %.4f (expected %s, %.1f se)" % (name, abs(m) if np.iscomplexobj(m) else m, expected,
                                                abs(m - expected) / se), abs(m - expected) <= 6 * se)


rng = np.random.default_rng(20261118)
I = lambda n: np.eye(n)

for n in (2, 3, 4, 5):
    x = ym.su_algebra(rng, n)
    check("su_algebra(%d) anti-Hermitian, traceless" % n, np.allclose(x, -x.conj().T, atol=1e-12) and abs(np.trace(x)) < 1e-12)

    g = ym.gell_mann(n)
    gram = np.array([[np.trace(a @ b) for b in g] for a in g])
    check("gell_mann(%d): %d generators, Hermitian, traceless, tr(TaTb) = delta/2" % (n, n * n - 1),
          len(g) == n * n - 1 and all(np.allclose(t, t.conj().T) and abs(np.trace(t)) < 1e-12 for t in g)
          and np.allclose(gram, np.eye(n * n - 1) / 2, atol=1e-12))
    comp = sum(np.einsum("ij,kl->ijkl", t, t) for t in g)
    d = np.eye(n)
    want = 0.5 * (np.einsum("il,jk->ijkl", d, d) - np.einsum("ij,kl->ijkl", d, d) / n)
    check("gell_mann(%d) completeness sum_a T_ij T_kl = (d_il d_jk - d_ij d_kl / N) / 2" % n, np.allclose(comp, want, atol=1e-12))
    cas = sum(t @ t for t in g)
    check("gell_mann(%d) Casimir sum_a T_a^2 = (N^2-1)/(2N) I" % n, np.allclose(cas, (n * n - 1) / (2 * n) * I(n), atol=1e-12))

    u = ym.haar_unitary(rng, n)
    s = ym.haar_su(rng, n)
    check("haar_unitary(%d) unitary" % n, np.allclose(u @ u.conj().T, I(n), atol=1e-10))
    check("haar_su(%d) unitary with det 1" % n, np.allclose(s @ s.conj().T, I(n), atol=1e-10) and abs(np.linalg.det(s) - 1) < 1e-10)

    e = ym.expm(x)
    w, v = np.linalg.eigh(-1j * x)
    check("expm(su(%d)) matches the spectral exponential and lies in SU(%d)" % (n, n),
          np.allclose(e, v @ np.diag(np.exp(1j * w)) @ v.conj().T, atol=1e-10) and abs(np.linalg.det(e) - 1) < 1e-10)

    m = 20000
    us = [ym.haar_unitary(rng, n) for _ in range(m)]
    tr = np.array([np.trace(a) for a in us])
    mc("U(%d)  E|tr U|^2" % n, np.abs(tr) ** 2, 1.0)
    mc("U(%d)  E|tr U^2|^2" % n, [abs(np.trace(a @ a)) ** 2 for a in us], float(min(2, n)))
    mc("U(%d)  E|U_11|^2" % n, [abs(a[0, 0]) ** 2 for a in us], 1.0 / n)
    mc("U(%d)  Re E (tr U)^N" % n, (tr ** n).real, 0.0)
    ss = np.array([np.trace(ym.haar_su(rng, n)) for _ in range(m)])
    mc("SU(%d) E|tr U|^2" % n, np.abs(ss) ** 2, 1.0)
    mc("SU(%d) Re E (tr U)^N" % n, (ss ** n).real, 1.0)
    mc("SU(%d) Im E (tr U)^N" % n, (ss ** n).imag, 0.0)

d = np.array([0.3, -1.2, 2.0])
check("expm(diag) = diag(exp)", np.allclose(ym.expm(np.diag(d)), np.diag(np.exp(d)), rtol=1e-12))

p = ym.prob_vector(rng, 6, sort_desc=True)
check("prob_vector: nonnegative, sums to 1, decreasing", (p >= 0).all() and abs(p.sum() - 1) < 1e-12 and (np.diff(p) <= 0).all())

for n in (2, 4, 6):
    h, p0, psi, en, v = ym.spectral_instance(rng, n, gap=0.7)
    ev = np.linalg.eigvalsh(h)
    check("spectral_instance(%d): Hermitian, ground energy 0, gap >= 0.7" % n,
          np.allclose(h, h.conj().T) and abs(ev[0]) < 1e-10 and ev[1] >= 0.7 - 1e-10)
    check("spectral_instance(%d): P0 is the projection onto ker H" % n,
          np.allclose(p0 @ p0, p0) and np.allclose(h @ p0, 0, atol=1e-10) and abs(np.trace(p0) - 1) < 1e-10)
    lam, wt = ym.spectral_measure(h, psi)
    s_ = 0.37
    lhs = (psi.conj() @ (ym.expm(-s_ * h) - p0) @ psi).real
    rhs = np.sum(wt * np.exp(-s_ * lam)) - np.sum(wt[np.abs(lam) < 1e-9])
    check("spectral_measure(%d): the vacuum eigenvalue is exactly 0 and nothing else lies in (0, gap)" % n,
          np.sum(lam == 0.0) == 1 and not np.any((lam > 0) & (lam < 0.7 - 1e-10)))
    check("spectral_measure(%d): weights sum to |psi|^2, Laplace transform matches <psi,(e^-sH - P0)psi>" % n,
          abs(wt.sum() - 1) < 1e-10 and abs(lhs - rhs) < 1e-10)

a = ym.psd(rng, 5, rank=2)
ev = np.linalg.eigvalsh(a)
check("psd(5, rank 2): PSD with rank 2", ev.min() > -1e-10 and int((ev > 1e-8).sum()) == 2)

print("\n%d failed" % len(FAILS))
sys.exit(1 if FAILS else 0)
