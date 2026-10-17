"""Maple-leaf problems from TheoPhysEval (arXiv:2610.16733, Appendix A, Quantum phases of matter).

Three problems share one six-sublattice lattice:

  dirac   -- maple-leaf-dirac-nodes: nodes of the middle bands, their energy, count,
             symmetry status, and the Berry phase (holonomy of the band-3 Berry connection)
  chern   -- maple-leaf-chern-number: Chern number of the occupied bands of the chiral
             deformation h_phi, direct and indirect gaps
  ed      -- maple-leaf-exact-diagonalization: 18-site Heisenberg cluster, fixed-S^z
             Lanczos, total spin, translation eigenvalue, spin gap

Run:  python maple_leaf.py [dirac|chern|ed|all]

Every answer is computed from the Hamiltonian in the prompt; nothing is looked up. The
official answer keys live in the gated Hugging Face bundle, so these are self-checked
values (convergence, independent methods, symmetry checks), not graded ones.
"""

from __future__ import annotations

import itertools
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

# ---------------------------------------------------------------------------
# Bloch Hamiltonians
# ---------------------------------------------------------------------------

T_D, T_H, T_T = 1.0, -0.73, 0.19


def bloch_hamiltonian(u, v, phi: float = 0.0, lam: float = 0.0,
                      t_d: float = T_D, t_h: float = T_H, t_t: float = T_T) -> np.ndarray:
    """h_phi(k) at fractional momenta (u, v); u, v may be arrays of equal shape.

    phi = 0 is the Dirac-node problem's h(k) (with lambda_z = lam); phi != 0 puts
    t_h e^{+i phi} on every bond of the oriented cycle 1->2->...->6->1.
    Returns an array of shape u.shape + (6, 6).
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    P = np.exp(2j * np.pi * u)
    Q = np.exp(2j * np.pi * v)
    Pi, Qi = 1 / P, 1 / Q
    one = np.ones_like(P)
    f, b = t_h * np.exp(1j * phi), t_h * np.exp(-1j * phi)  # forward a->a+1, backward
    rows = [
        [lam * one, f * one, t_t * P * Q, t_d * Q, t_t * Q, b * one],
        [b * one, lam * one, f * one, t_t * Q, t_d * Pi, t_t * Pi],
        [t_t * Pi * Qi, b * one, lam * one, f * one, t_t * Pi, t_d * Pi * Qi],
        [t_d * Qi, t_t * Qi, b * one, lam * one, f * one, t_t * Pi * Qi],
        [t_t * Qi, t_d * P, t_t * P, b * one, lam * one, f * one],
        [f * one, t_t * P, t_d * P * Q, t_t * P * Q, b * one, lam * one],
    ]
    h = np.stack([np.stack(r, axis=-1) for r in rows], axis=-2)
    return h


def check_hermitian(phi: float) -> float:
    rng = np.random.default_rng(0)
    u, v = rng.random(50), rng.random(50)
    h = bloch_hamiltonian(u, v, phi)
    return float(np.max(np.abs(h - np.conj(np.swapaxes(h, -1, -2)))))


def bands(u, v, phi: float = 0.0) -> np.ndarray:
    return np.linalg.eigvalsh(bloch_hamiltonian(u, v, phi))


def c6(u, v):
    """Sixfold rotation on fractional momenta: R b1 = b2, R b2 = -b1 + b2."""
    return -v, u + v


def c6_spectrum_defect(phi: float, n: int = 400) -> float:
    """max |eps(Rk) - eps(k)| over random k; zero if the spectrum is C6 symmetric."""
    rng = np.random.default_rng(1)
    u, v = rng.random(n), rng.random(n)
    ru, rv = c6(u, v)
    return float(np.max(np.abs(bands(u, v, phi) - bands(ru, rv, phi))))


# ---------------------------------------------------------------------------
# Problem: Dirac nodes
# ---------------------------------------------------------------------------

def middle_gap(x, phi: float = 0.0) -> float:
    e = np.linalg.eigvalsh(bloch_hamiltonian(x[0], x[1], phi))
    return float(e[3] - e[2])


def find_nodes(phi: float = 0.0, n_grid: int = 400, tol_gap: float = 1e-9):
    g = np.linspace(0.0, 1.0, n_grid, endpoint=False)
    U, V = np.meshgrid(g, g, indexing="ij")
    e = bands(U, V, phi)
    gap = e[..., 3] - e[..., 2]
    # Local minima on the periodic grid.
    is_min = np.ones_like(gap, dtype=bool)
    for du, dv in itertools.product((-1, 0, 1), repeat=2):
        if du == 0 and dv == 0:
            continue
        is_min &= gap <= np.roll(np.roll(gap, du, 0), dv, 1)
    seeds = np.argwhere(is_min & (gap < 0.05))
    nodes = []
    for i, j in seeds:
        x0 = np.array([g[i], g[j]])
        res = minimize(middle_gap, x0, args=(phi,), method="Nelder-Mead",
                       options=dict(xatol=1e-13, fatol=1e-15, maxiter=20000,
                                    initial_simplex=[x0, x0 + [2e-3, 0], x0 + [0, 2e-3]]))
        if res.fun < tol_gap:
            x = np.mod(res.x, 1.0)
            x[np.isclose(x, 1.0, atol=1e-9)] = 0.0
            if not any(np.max(np.abs(((x - y + 0.5) % 1.0) - 0.5)) < 1e-6 for y, _ in nodes):
                e_node = np.linalg.eigvalsh(bloch_hamiltonian(x[0], x[1], phi))
                nodes.append((x, 0.5 * (e_node[2] + e_node[3])))
    nodes.sort(key=lambda t: (round(t[0][0], 6), round(t[0][1], 6)))
    return nodes, float(gap.min())


def band_berry_phase(center, band: int, radius: float = 1e-4, n: int = 400, phi: float = 0.0) -> float:
    """Berry phase of one band around a counterclockwise circle in Cartesian k.

    gamma = -Arg prod_j <psi(k_j)|psi(k_{j+1})>, gauge invariant for a closed loop.
    The circle is taken in Cartesian momentum so 'counterclockwise' has its plane meaning.
    """
    b1 = 2 * np.pi / 21 * np.array([5 * np.sqrt(3), 3.0])
    b2 = 2 * np.pi / 21 * np.array([np.sqrt(3), 9.0])
    B = np.column_stack([b1, b2])          # k = B @ (u, v)
    k0 = B @ np.asarray(center)
    s = np.linspace(0, 2 * np.pi, n, endpoint=False)
    k = k0[:, None] + radius * np.vstack([np.cos(s), np.sin(s)])
    uv = np.linalg.solve(B, k)
    _, vecs = np.linalg.eigh(bloch_hamiltonian(uv[0], uv[1], phi))
    psi = vecs[:, :, band]                 # (n, 6)
    overlaps = np.einsum("ji,ji->j", psi.conj(), np.roll(psi, -1, axis=0))
    return float(-np.angle(np.prod(overlaps / np.abs(overlaps))))


def on_gamma_m_line(x, tol: float = 1e-6) -> bool:
    """Is (u,v) on a line through Gamma and an M point (b1/2, b2/2, (b1-b2)/2 or a C6 image),
    up to a reciprocal-lattice vector?"""
    orbit = []
    for m in ([0.5, 0.0], [0.0, 0.5], [0.5, -0.5]):
        p = np.array(m)
        for _ in range(6):
            orbit.append(p.copy())
            p = np.array(c6(p[0], p[1]))
    for d in orbit:
        for shift in itertools.product((-1, 0, 1), repeat=2):
            y = np.asarray(x) + shift
            if abs(y[0] * d[1] - y[1] * d[0]) < tol:
                return True
    return False


def solve_dirac():
    print("== maple-leaf-dirac-nodes ==")
    print(f"Hermiticity defect: {check_hermitian(0.0):.1e}")
    print(f"C6 spectrum defect: phi=0: {c6_spectrum_defect(0.0):.1e}, phi=0.15: {c6_spectrum_defect(0.15):.1e}")
    nodes, grid_min = find_nodes(0.0)
    print(f"grid minimum of eps4-eps3 on 400^2 mesh: {grid_min:.3e}")
    for x, e in nodes:
        gamma3 = band_berry_phase(x, band=2)
        gamma_small = band_berry_phase(x, band=2, radius=3e-5, n=800)
        print(f"node (u,v)=({x[0]:.8f}, {x[1]:.8f})  E={e:+.10f}  on Gamma-M line: {on_gamma_m_line(x)}  "
              f"Berry phase band 3: {gamma3:+.6f} (r=1e-4), {gamma_small:+.6f} (r=3e-5)")
    print(f"number of inequivalent nodes: {len(nodes)}")
    energies = [e for _, e in nodes]
    print(f"spread of nodal energies: {max(energies) - min(energies):.2e}")

    # Symmetry test: p6 contains C6 but no mirrors. Turn on the chiral deformation
    # (C6-symmetric per the defect printed above, mirror-breaking) and see whether the
    # nodes survive on the Gamma-M lines.
    x0 = nodes[0][0]
    for phi in (1e-4, 1e-3, 1e-2):
        res = min((minimize(middle_gap, x0 + d, args=(phi,), method="Nelder-Mead",
                            options=dict(xatol=1e-12, fatol=1e-15, maxiter=20000))
                   for d in (np.zeros(2), np.array([1e-2, 0]), np.array([0, 1e-2]))), key=lambda r: r.fun)
        print(f"phi={phi:g}: minimum middle gap near the node {res.fun:.3e} (gap/phi = {res.fun / phi:.4f})")
    return nodes


# ---------------------------------------------------------------------------
# Problem: Chern number and gaps
# ---------------------------------------------------------------------------

def chern_occupied(phi: float, n: int, n_occ: int = 3) -> float:
    g = np.arange(n) / n
    U, V = np.meshgrid(g, g, indexing="ij")
    _, vecs = np.linalg.eigh(bloch_hamiltonian(U, V, phi))
    occ = vecs[..., :n_occ]                               # (n, n, 6, 3)

    def link(shift_axis):
        nb = np.roll(occ, -1, axis=shift_axis)            # V(k + delta_mu), periodic
        d = np.linalg.det(np.einsum("...ia,...ib->...ab", occ.conj(), nb))
        return d / np.abs(d)

    Uu, Uv = link(0), link(1)
    F = np.angle(Uu * np.roll(Uv, -1, axis=0) / (np.roll(Uu, -1, axis=1) * Uv))
    return float(F.sum() / (2 * np.pi))


def band_extrema(phi: float, n_grid: int = 300):
    g = np.linspace(0, 1, n_grid, endpoint=False)
    U, V = np.meshgrid(g, g, indexing="ij")
    e = bands(U, V, phi)

    def refine(fun, idx):
        best = None
        for i, j in idx:
            res = minimize(fun, np.array([g[i], g[j]]), method="Nelder-Mead",
                           options=dict(xatol=1e-12, fatol=1e-14, maxiter=20000))
            if best is None or res.fun < best.fun:
                best = res
        return best

    def k_smallest(a, k=8):
        flat = np.argsort(a, axis=None)[:k]
        return [np.unravel_index(f, a.shape) for f in flat]

    gap = e[..., 3] - e[..., 2]
    direct = refine(lambda x: np.diff(np.linalg.eigvalsh(bloch_hamiltonian(x[0], x[1], phi))[2:4])[0],
                    k_smallest(gap))
    e4min = refine(lambda x: np.linalg.eigvalsh(bloch_hamiltonian(x[0], x[1], phi))[3],
                   k_smallest(e[..., 3]))
    e3max = refine(lambda x: -np.linalg.eigvalsh(bloch_hamiltonian(x[0], x[1], phi))[2],
                   k_smallest(-e[..., 2]))
    return dict(direct=direct.fun, direct_at=np.mod(direct.x, 1),
                e4_min=e4min.fun, e4_min_at=np.mod(e4min.x, 1),
                e3_max=-e3max.fun, e3_max_at=np.mod(e3max.x, 1),
                grid_direct=float(gap.min()),
                grid_indirect=float(e[..., 3].min() - e[..., 2].max()))


def solve_chern(phi: float = 0.15):
    print("== maple-leaf-chern-number ==")
    print(f"Hermiticity defect: {check_hermitian(phi):.1e}")
    for n in (24, 48, 96, 192):
        print(f"mesh {n:4d}^2: C_occ = {chern_occupied(phi, n):+.10f}")
    ex = band_extrema(phi)
    print(f"direct gap   = {ex['direct']:.10f} at (u,v)={ex['direct_at']}  (grid: {ex['grid_direct']:.6f})")
    ind = ex["e4_min"] - ex["e3_max"]
    print(f"min eps4     = {ex['e4_min']:.10f} at {ex['e4_min_at']}")
    print(f"max eps3     = {ex['e3_max']:.10f} at {ex['e3_max_at']}")
    print(f"indirect gap = {ind:.10f}  (grid: {ex['grid_indirect']:.6f})")
    return ex


# ---------------------------------------------------------------------------
# Problem: 18-site exact diagonalization
# ---------------------------------------------------------------------------

J_H, J_T, J_D = 1.0, -0.58, -0.31
RULES = {
    J_H: [(0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 0, 0)],
    J_T: [(0, 2, 0), (0, 4, 2), (1, 3, 2), (1, 5, 2), (2, 4, 2), (3, 5, 0)],
    J_D: [(0, 3, 2), (1, 4, 2), (2, 5, 0)],
}
N_CELLS, N_SUB = 3, 6
N_SITES = N_CELLS * N_SUB


def site(c: int, a: int) -> int:
    return 6 * (c % N_CELLS) + a


def bonds():
    out, seen = [], set()
    for J, rules in RULES.items():
        for c in range(N_CELLS):
            for a, b, d in rules:
                i, j = site(c, a), site(c + d, b)
                key = (min(i, j), max(i, j))
                assert i != j, "self bond"
                assert key not in seen, f"duplicate bond {key}"
                seen.add(key)
                out.append((i, j, J))
    return out


def sector_states(n_up: int) -> np.ndarray:
    states = [sum(1 << k for k in combo) for combo in itertools.combinations(range(N_SITES), n_up)]
    return np.array(sorted(states), dtype=np.int64)


def heisenberg_sector(n_up: int):
    st = sector_states(n_up)
    index = {int(s): k for k, s in enumerate(st)}
    bl = bonds()
    rows, cols, vals = [], [], []
    diag = np.zeros(len(st))
    for k, s in enumerate(st):
        s = int(s)
        for i, j, J in bl:
            si, sj = (s >> i) & 1, (s >> j) & 1
            diag[k] += J * (0.25 if si == sj else -0.25)
            if si != sj:                       # S+S- + S-S+ flips the pair, amplitude J/2
                t = s ^ ((1 << i) | (1 << j))
                rows.append(index[t]); cols.append(k); vals.append(0.5 * J)
    rows += list(range(len(st))); cols += list(range(len(st))); vals += list(diag)
    H = csr_matrix((vals, (rows, cols)), shape=(len(st), len(st)))
    return H, st, index


def translate_state(s: int) -> int:
    out = 0
    for c in range(N_CELLS):
        for a in range(N_SUB):
            if (s >> site(c, a)) & 1:
                out |= 1 << site(c + 1, a)
    return out


def solve_ed(k_eigs: int = 6):
    print("== maple-leaf-exact-diagonalization ==")
    bl = bonds()
    counts = {J: sum(1 for *_, x in bl if x == J) for J in RULES}
    print(f"bond counts: J_h={counts[J_H]}, J_t={counts[J_T]}, J_d={counts[J_D]}")
    H0, st0, idx0 = heisenberg_sector(9)
    H1, _, _ = heisenberg_sector(10)
    print(f"dim(S^z=0) = {H0.shape[0]}, dim(S^z=1) = {H1.shape[0]}")
    w0, v0 = eigsh(H0, k=k_eigs, which="SA", tol=1e-13)
    w1, _ = eigsh(H1, k=k_eigs, which="SA", tol=1e-13)
    o0, o1 = np.argsort(w0), np.argsort(w1)
    w0, v0, w1 = w0[o0], v0[:, o0], w1[o1]
    print("lowest S^z=0:", np.array2string(w0, precision=10))
    print("lowest S^z=1:", np.array2string(w1, precision=10))
    E0 = w0[0]
    # SU(2) consistency: each S^z=1 level must also appear in the S^z=0 spectrum.
    for e in w1[:3]:
        print(f"  S^z=1 level {e:.10f}: closest S^z=0 level differs by {np.min(np.abs(w0 - e)):.2e}"
              f"{' (beyond computed S^z=0 window)' if e > w0[-1] else ''}")
    S = 0 if E0 < w1[0] - 1e-9 else None
    print(f"E0 = {E0:.10f}, E0/18 = {E0 / 18:.10f}, ground-state degeneracy in S^z=0: "
          f"{int(np.sum(np.abs(w0 - E0) < 1e-8))}")
    print(f"total spin S = {S}  (E0 below every S^z=1 level -> no S>=1 multiplet contains it)")
    perm = np.array([idx0[translate_state(int(s))] for s in st0])
    psi = v0[:, 0]
    Tpsi = np.zeros_like(psi)
    Tpsi[perm] = psi                       # (T psi)(T s) = psi(s)
    lam = np.vdot(psi, Tpsi)
    print(f"<psi0|T|psi0> = {lam.real:+.10f}{lam.imag:+.10f}i  -> k = {np.angle(lam) % (2 * np.pi):.6f}")
    print(f"spin gap = {w1[0] - E0:.10f}")
    return dict(E0=E0, S=S, T=lam, gap=w1[0] - E0)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("dirac", "all"):
        solve_dirac()
    if which in ("chern", "all"):
        solve_chern()
    if which in ("ed", "all"):
        solve_ed()
