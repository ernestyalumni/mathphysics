"""ymlib: verified random instances for numeric tests of proof steps (YM-T18).

ym-numcheck mounts this file in its sandbox and exposes it to model-written tests as `ym`, so a test
only has to state the claim; the instance comes from code that `ym-numcheck selftest` checks.
Every sampler takes a numpy Generator first. Conventions follow the mono-notes: su(N) is the
traceless anti-Hermitian matrices, generators T_a are Hermitian with tr(T_a T_b) = delta_ab / 2.

Label: unreviewed tooling.
"""
import numpy as np


def complex_gaussian(rng, n, m=None):
    """Matrix (or vector, if m is None and n is an int) of i.i.d. standard complex Gaussians."""
    shape = (n,) if m is None else (n, m)
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def hermitian(rng, n, traceless=False):
    """Random Hermitian n x n matrix (GUE-like); traceless if asked."""
    a = complex_gaussian(rng, n, n)
    h = (a + a.conj().T) / 2
    if traceless:
        h = h - np.trace(h) / n * np.eye(n)
    return h


def su_algebra(rng, n):
    """Random element of su(n): traceless and anti-Hermitian."""
    return 1j * hermitian(rng, n, traceless=True)


def u_algebra(rng, n):
    """Random element of u(n): anti-Hermitian."""
    return 1j * hermitian(rng, n)


def haar_unitary(rng, n):
    """Haar-random U(n) (QR of a complex Gaussian with the phases of R's diagonal removed; Mezzadri 2007)."""
    q, r = np.linalg.qr(complex_gaussian(rng, n, n))
    d = np.diag(r)
    return q * (d / np.abs(d))


def haar_su(rng, n):
    """Haar-random SU(n): U det(U)^(-1/n) for Haar U in U(n). The map commutes with left multiplication
    by SU(n), so it carries Haar measure on U(n) to Haar measure on SU(n)."""
    u = haar_unitary(rng, n)
    return u * np.linalg.det(u) ** (-1.0 / n)


def gell_mann(n):
    """Hermitian traceless generators T_a of su(n), a = 1..n^2-1, with tr(T_a T_b) = delta_ab / 2."""
    gens = []
    for j in range(n):
        for k in range(j + 1, n):
            s = np.zeros((n, n), complex); s[j, k] = s[k, j] = 0.5; gens.append(s)
            a = np.zeros((n, n), complex); a[j, k] = -0.5j; a[k, j] = 0.5j; gens.append(a)
    for l in range(1, n):
        d = np.zeros((n, n), complex)
        d[np.arange(l), np.arange(l)] = 1.0
        d[l, l] = -l
        gens.append(d / np.sqrt(2 * l * (l + 1)))
    return gens


def expm(a):
    """Matrix exponential by scaling and squaring with a Taylor series (exact enough for small matrices)."""
    a = np.asarray(a, dtype=complex if np.iscomplexobj(a) else float)
    s = max(0, int(np.ceil(np.log2(max(np.linalg.norm(a, 1), 1e-300)))) + 1)
    b = a / 2.0 ** s
    term = np.eye(a.shape[-1], dtype=a.dtype)
    out = term.copy()
    for k in range(1, 30):
        term = term @ b / k
        out = out + term
    for _ in range(s):
        out = out @ out
    return out


def prob_vector(rng, n, sort_desc=False):
    """Random probability vector (Dirichlet(1,...,1)); sorted in decreasing order if asked."""
    p = rng.exponential(size=n)
    p = p / p.sum()
    return np.sort(p)[::-1] if sort_desc else p


def spectral_instance(rng, n, gap=None):
    """A finite model of a Hamiltonian with vacuum: H = V diag(0, E_2..E_n) V^dagger with E_i >= gap > 0,
    the vacuum projection P0, and a random unit vector psi. Returns (H, P0, psi, energies, V)."""
    gap = rng.uniform(0.2, 2.0) if gap is None else gap
    e = np.concatenate([[0.0], gap + rng.exponential(size=n - 1)])
    v = haar_unitary(rng, n)
    h = v @ np.diag(e) @ v.conj().T
    p0 = np.outer(v[:, 0], v[:, 0].conj())
    psi = complex_gaussian(rng, n)
    return h, p0, psi / np.linalg.norm(psi), e, v


def spectral_measure(h, psi):
    """Atoms of the spectral measure nu_psi of a Hermitian h: (eigenvalues, weights |<v_i, psi>|^2)."""
    e, v = np.linalg.eigh(h)
    return e, np.abs(v.conj().T @ psi) ** 2


def psd(rng, n, rank=None):
    """Random positive semidefinite n x n matrix of the given rank (default full)."""
    a = complex_gaussian(rng, n, rank or n)
    return a @ a.conj().T


def close(a, b, rtol=1e-7, atol=1e-9):
    """Elementwise closeness for scalars or arrays (complex allowed)."""
    return bool(np.allclose(a, b, rtol=rtol, atol=atol))
