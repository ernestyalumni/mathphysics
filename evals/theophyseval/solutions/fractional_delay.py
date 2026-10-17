"""autonomous-de (TheoPhysEval, Classical physics): the inverse problem no agent solved.

The prompt shows a curve f(t), 0 <= t <= 5, on a log axis and asks for "an autonomous
differential equation where no derivative is of order greater than 1". The authors'
answer (paper, Appendix "Solutions") is a fractional *delay* equation,

    f'(t) = D^{1/5} f(t) - 5 D^{1/3} f(t) + 500 f(t-1) + 550,   t > 0,
    f(t) = -11/10 on [-1, 0),   f(0) = 1.0073e-5,

with Riemann-Liouville derivatives D^a f(t) = (1/Gamma(1-a)) d/dt ∫_0^t (t-s)^-a f(s) ds.
The fractional orders are below 1, so "no derivative of order greater than 1" holds.

This script
  1. solves the equation (implicit Euler + Grunwald-Letnikov weights, first order in h),
  2. digitizes the published figure fig1.png (frame-calibrated, log-y),
  3. reports the misfit in decades and a convergence table,
  4. reads the structure off the curve the way a human would, and fits the four
     coefficients by least squares to show the inverse problem is well posed once the
     ansatz is right.

Run: python fractional_delay.py [--figure PATH] [--plot OUT.png]
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

DEFAULT_FIGURE = Path("/media/propdev/Expansion/openclaw/.openclaw/workspace/Data/Public/papers/AI/"
                      "arXiv-2610.16733v1/fig1.png")


@dataclass
class FDE:
    """f' = c1 D^a1 f + c2 D^a2 f + c3 f(t - tau) + c4, history h0 on [-tau, 0), f(0) = f0."""
    c1: float = 1.0
    a1: float = 1 / 5
    c2: float = -5.0
    a2: float = 1 / 3
    c3: float = 500.0
    c4: float = 550.0
    tau: float = 1.0
    h0: float = -1.1
    f0: float = 1.0073e-5


def gl_weights(alpha: float, n: int) -> np.ndarray:
    """Grunwald-Letnikov weights w_j = (-1)^j binom(alpha, j)."""
    w = np.empty(n + 1)
    w[0] = 1.0
    for j in range(1, n + 1):
        w[j] = w[j - 1] * (1 - (alpha + 1) / j)
    return w


def solve(p: FDE, t_end: float = 5.0, steps_per_tau: int = 2000):
    """Implicit Euler for f' with the fractional terms at t_{n+1}:
        D^a f(t_{n+1}) ~ h^-a sum_{j=0}^{n+1} w_j f_{n+1-j},
    which is linear in f_{n+1}. The delay lands on the grid because h = tau / steps_per_tau."""
    h = p.tau / steps_per_tau
    n_steps = int(round(t_end / h))
    w1, w2 = gl_weights(p.a1, n_steps), gl_weights(p.a2, n_steps)
    k1, k2 = h ** -p.a1, h ** -p.a2
    f = np.empty(n_steps + 1)
    f[0] = p.f0
    lhs = 1 / h - p.c1 * k1 - p.c2 * k2
    for n in range(n_steps):
        m = n + 1                                    # computing f[m]
        past = f[m - 1::-1]                          # f[m-1], f[m-2], ..., f[0]
        s1 = w1[1:m + 1] @ past
        s2 = w2[1:m + 1] @ past
        j = m - steps_per_tau
        f_delay = p.h0 if j < 0 else f[j]
        rhs = f[n] / h + p.c1 * k1 * s1 + p.c2 * k2 * s2 + p.c3 * f_delay + p.c4
        f[m] = rhs / lhs
    return np.linspace(0, n_steps * h, n_steps + 1), f


# ---------------------------------------------------------------------------
# Digitizing the published figure
# ---------------------------------------------------------------------------

def digitize(path: Path, t_range=(0.0, 5.0), log10_range=(-7.0, 8.0), tick_margin: int = 32):
    """Return (t, log10 f) sampled from the black curve in fig1.png.

    The plot frame is found as the rows/columns with the longest dark runs; its edges are
    t = 0, 5 and f = 1e-7, 1e8 (the axis labels sit on the frame). Columns within
    tick_margin pixels of the left/right frame and rows near the top/bottom frame are
    skipped because the tick marks live there.
    """
    from PIL import Image

    img = np.asarray(Image.open(path).convert("L"), dtype=float) / 255.0
    dark = img < 0.5
    row_runs, col_runs = dark.sum(axis=1), dark.sum(axis=0)
    rows = np.argsort(row_runs)[::-1]
    cols = np.argsort(col_runs)[::-1]
    top, bottom = sorted([rows[0], next(r for r in rows if abs(r - rows[0]) > 100)])
    left, right = sorted([cols[0], next(c for c in cols if abs(c - cols[0]) > 100)])
    ts, ys = [], []
    y_prev = None
    for x in range(left + tick_margin, right - tick_margin):
        hits = np.nonzero(dark[top + 20:bottom - 20, x])[0] + top + 20
        if hits.size == 0:
            continue
        # Split into contiguous runs; a major tick mark is a second run in the same column.
        runs = np.split(hits, np.nonzero(np.diff(hits) > 1)[0] + 1)
        run = runs[0] if y_prev is None else min(runs, key=lambda r: abs(r.mean() - y_prev))
        if run.size > 12:                            # the near-vertical jump at t = 1
            continue
        y = run.mean()
        y_prev = y
        ts.append(t_range[0] + (x - left) / (right - left) * (t_range[1] - t_range[0]))
        ys.append(log10_range[1] - (y - top) / (bottom - top) * (log10_range[1] - log10_range[0]))
    return np.array(ts), np.array(ys), dict(left=left, right=right, top=top, bottom=bottom)


def jump_mask(t_dig, t_jump: float = 1.0, width: float = 0.01):
    """Drop points on the near-vertical jump at t = 1: there a one-pixel horizontal offset
    is a ~6 decade vertical error, so a vertical misfit is meaningless."""
    return np.abs(t_dig - t_jump) > width


def misfit_decades(t_sol, f_sol, t_dig, y_dig):
    """RMS and max |log10 f_model - log10 f_digitized| at the digitized abscissae
    (jump window excluded)."""
    keep = jump_mask(t_dig)
    t_dig, y_dig = t_dig[keep], y_dig[keep]
    y_model = np.log10(np.interp(t_dig, t_sol, f_sol))
    d = y_model - y_dig
    return float(np.sqrt(np.mean(d ** 2))), float(np.max(np.abs(d)))


# ---------------------------------------------------------------------------
# The inverse step: structure from the curve, coefficients by least squares
# ---------------------------------------------------------------------------

def fit_coefficients(t_dig, y_dig, start: FDE, steps_per_tau: int = 1000):
    """Fit (c1, c2, c3) with c4 = -h0 c3 and the orders, delay and history fixed.

    Reading the figure fixes the rest: kinks at t = 1, 2, 3, 4 give tau = 1; on [0, 1) the
    curve only drifts from 1e-5 to 1e-6, so the delay and constant terms must cancel against
    the constant history, c3 h0 + c4 = 0. Without that constraint a generic start grows to
    O(10) on [0, 1), seven decades off, and least squares stalls in a far basin.
    """
    keep = jump_mask(t_dig)
    t_dig, y_dig = t_dig[keep], y_dig[keep]

    def model(x):
        return FDE(**{**start.__dict__, "c1": x[0], "c2": x[1], "c3": x[2], "c4": -start.h0 * x[2]})

    def residual(x):
        t, f = solve(model(x), steps_per_tau=steps_per_tau)
        f = np.where(f > 0, f, 1e-30)
        return np.log10(np.interp(t_dig, t, f)) - y_dig

    x0 = np.array([start.c1, start.c2, start.c3])
    res = least_squares(residual, x0, x_scale=np.abs(x0), diff_step=1e-4, max_nfev=300)
    return res, model(res.x)


def identifiability(t_dig, c=(1.0, -5.0, 500.0), h0: float = -1.1, steps_per_tau: int = 1000):
    """Singular values of d log10 f(t_i) / d log c_k at the published point (c4 = -h0 c3).
    A small singular value is a coefficient combination the figure cannot see."""
    t_keep = t_dig[jump_mask(t_dig)]
    c0 = np.asarray(c, float)

    def y(cc):
        t, f = solve(FDE(c1=cc[0], c2=cc[1], c3=cc[2], c4=-h0 * cc[2]), steps_per_tau=steps_per_tau)
        return np.log10(np.interp(t_keep, t, f))

    cols = []
    for k in range(3):
        d = np.zeros(3)
        d[k] = 1e-3 * abs(c0[k])
        cols.append((y(c0 + d) - y(c0 - d)) / (2 * d[k]) * abs(c0[k]))
    _, S, Vt = np.linalg.svd(np.array(cols).T, full_matrices=False)
    return S, Vt, t_keep.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    ap.add_argument("--plot", type=Path, default=None)
    ap.add_argument("--fit", action="store_true", help="also run the coefficient fit")
    args = ap.parse_args()

    p = FDE()
    t_dig, y_dig, frame = digitize(args.figure)
    print(f"digitized {t_dig.size} columns; frame px {frame}")

    print("\nforward solve of the published equation (implicit Euler + Grunwald-Letnikov):")
    sols = {}
    for m in (250, 500, 1000, 2000, 4000):
        t, f = solve(p, steps_per_tau=m)
        sols[m] = (t, f)
        rms, mx = misfit_decades(t, f, t_dig, y_dig)
        at = {tt: f[int(round(tt * m))] for tt in (2.0, 3.0, 4.0, 5.0)}
        at[1.0] = f[m - 1]                           # f(1-): the delay switches on at t = 1
        print(f"  h=1/{m:<5d} f(1-)={at[1.0]:.4e} f(2)={at[2.0]:.4e} f(3)={at[3.0]:.4e} "
              f"f(4)={at[4.0]:.4e} f(5)={at[5.0]:.4e}   misfit vs figure: rms {rms:.3f}, max {mx:.3f} decades")

    if args.fit:
        print("\ninverse step: fit (c1, c2, c3, c4) from a perturbed start, orders/delay/history fixed")
        start = FDE(c1=1.4, c2=-3.5, c3=420.0, c4=1.1 * 420.0)
        res, q = fit_coefficients(t_dig, y_dig, start)
        print(f"  constraint c4 = -h0 c3 = 1.1 c3 (flat segment on [0, 1))")
        print(f"  start  (c1, c2, c3, c4) = ({start.c1}, {start.c2}, {start.c3}, {start.c4:.1f})")
        print(f"  fitted (c1, c2, c3, c4) = ({q.c1:.4f}, {q.c2:.4f}, {q.c3:.3f}, {q.c4:.3f})  "
              f"rms {np.sqrt(np.mean(res.fun ** 2)):.4f} decades, {res.nfev} solves")
        print("  published               = (1, -5, 500, 550)")
        S, Vt, n = identifiability(t_dig)
        print(f"  Jacobian singular values (decades per unit relative change): {np.round(S, 3)}")
        print(f"  weakest direction in (dc1/c1, dc2/c2, dc3/c3): {np.round(Vt[-1], 3)}; a 10% move along it "
              f"changes the curve by rms {0.1 * S[-1] / np.sqrt(n):.4f} decades (one pixel ~ 0.017)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        t, f = sols[2000]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.semilogy(t, np.where(f > 0, f, np.nan), lw=1.5, label="solution, h = 1/2000")
        ax.semilogy(t_dig, 10 ** y_dig, ".", ms=2, alpha=0.5, label="digitized fig1.png")
        ax.set_xlabel("t"); ax.set_ylabel("f(t)"); ax.set_ylim(1e-7, 1e8); ax.legend()
        ax.set_title("f' = D^{1/5}f - 5 D^{1/3}f + 500 f(t-1) + 550")
        fig.tight_layout(); fig.savefig(args.plot, dpi=150)
        print(f"\nwrote {args.plot}")


if __name__ == "__main__":
    main()
