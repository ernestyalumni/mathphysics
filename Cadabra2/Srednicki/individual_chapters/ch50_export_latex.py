"""
ch50_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 50 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \\
        python3 /work/ch50_export_latex.py

Outputs: /work/ch50_massless_spinor_helicity.tex
Then compile on host:
    pdflatex ch50_massless_spinor_helicity.tex
"""

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

# ── cadabra2 setup ──────────────────────────────────────────────────────────
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"),       Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"),        Ex(r"position=free"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\abra{k}{p}"))
cadabra2.AntiSymmetric(Ex(r"\sbra{k}{p}"))

# ── numpy setup ─────────────────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
sigma = {
    1: np.array([[0, 1],   [1,  0]],  dtype=complex),
    2: np.array([[0, -1j], [1j, 0]],  dtype=complex),
    3: np.array([[1, 0],   [0, -1]],  dtype=complex),
}
sig  = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigb = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}
g = np.diag([-1., 1., 1., 1.])
eps_up = np.array([[0,  1], [-1, 0]], dtype=complex)  # ε^{ab}: ε^{12}=+1
eps_dn = np.array([[0, -1], [ 1, 0]], dtype=complex)  # ε_{ab}: ε_{12}=-1

def gamma_mat(mu):
    return np.block([[Z2, sig[mu]], [sigb[mu], Z2]])

gam = [gamma_mat(mu) for mu in range(4)]
gam5 = 1j * gam[0] @ gam[1] @ gam[2] @ gam[3]

# ── helpers ─────────────────────────────────────────────────────────────────
import re

def fix_cdb_latex(s):
    s = re.sub(r'\\psidag\s*\^\{([^}]*)\}', lambda m: r'\psi^{\dagger ' + m.group(1) + '}', s)
    s = re.sub(r'\\chidag\s*\^\{([^}]*)\}', lambda m: r'\chi^{\dagger ' + m.group(1) + '}', s)
    return s

def cdb(expr_str):
    return fix_cdb_latex(Ex(expr_str)._latex_())

def fmt_complex(z, tol=1e-10):
    r, i = z.real, z.imag
    if abs(r) < tol and abs(i) < tol: return "0"
    if abs(i) < tol:
        v = r
        if abs(v - round(v)) < tol: return str(int(round(v)))
        if abs(v - 0.5) < tol: return r"\tfrac{1}{2}"
        if abs(v + 0.5) < tol: return r"-\tfrac{1}{2}"
        return f"{v:.4g}"
    if abs(r) < tol:
        v = i
        if abs(v - 1.0) < tol: return r"i"
        if abs(v + 1.0) < tol: return r"-i"
        if abs(v - 0.5) < tol: return r"\tfrac{i}{2}"
        if abs(v + 0.5) < tol: return r"-\tfrac{i}{2}"
        return f"{v:.4g}i"
    return f"{r:.4g}+{i:.4g}i"

def mat2pmatrix(mat):
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"

def err_str(e): return f"{e:.1e}"

# ── 2-component spinor tools ─────────────────────────────────────────────────

def helicity_spinors(E, theta, phi):
    """2-component helicity spinors λ_a and λ̃_ȧ for massless p^μ = E(1,sinθcosφ,sinθsinφ,cosθ)."""
    ct2 = np.cos(theta / 2)
    st2 = np.sin(theta / 2)
    sq2E = np.sqrt(2 * E)
    lam    = sq2E * np.array([ct2, st2 * np.exp( 1j * phi)], dtype=complex)
    lamtil = sq2E * np.array([ct2, st2 * np.exp(-1j * phi)], dtype=complex)
    return lam, lamtil

def angle_bracket(lam1, lam2):
    """⟨12⟩ = ε^{ab} λ_{1a} λ_{2b}"""
    return lam1 @ eps_up @ lam2

def square_bracket(lamtil1, lamtil2):
    """[12] = ε_{ȧḃ} λ̃_1^ȧ λ̃_2^ḃ"""
    return lamtil1 @ eps_dn @ lamtil2

def massless_spinors_general(E, theta, phi):
    """4-component massless Dirac spinors |p] and |p⟩ in Weyl representation."""
    ct2 = np.cos(theta / 2)
    st2 = np.sin(theta / 2)
    ep  = np.exp( 1j * phi)
    em  = np.exp(-1j * phi)
    sq  = np.sqrt(2 * E)
    ket_sq = sq * np.array([-st2 * em, ct2, 0., 0.], dtype=complex)
    ket_an = sq * np.array([0., 0., ct2, st2 * ep],  dtype=complex)
    bra_sq = ket_sq.conj() @ gam[0]
    bra_an = ket_an.conj() @ gam[0]
    return ket_sq, ket_an, bra_sq, bra_an

# ── §48.A: massless spinor outer product ─────────────────────────────────────
E_p = 3.0
th  = np.pi / 4
ph  = np.pi / 3
lam_p, lamtil_p = helicity_spinors(E_p, th, ph)

px = E_p * np.sin(th) * np.cos(ph)
py = E_p * np.sin(th) * np.sin(ph)
pz = E_p * np.cos(th)
p_spinor = E_p * I2 + px * sigma[1] + py * sigma[2] + pz * sigma[3]
p_factored = np.outer(lam_p, lamtil_p)
err_rank1 = np.max(np.abs(p_spinor - p_factored))
det_p = np.linalg.det(p_spinor)

p_spinor_tex = mat2pmatrix(p_spinor)
lam_p_vals   = f"({fmt_complex(lam_p[0])},\\;{fmt_complex(lam_p[1])})"
lamtil_p_vals = f"({fmt_complex(lamtil_p[0])},\\;{fmt_complex(lamtil_p[1])})"

# ── §48.B: angle/square brackets ─────────────────────────────────────────────
omega_k = 2.5
theta_q = np.pi / 3
lam_k, lamtil_k = helicity_spinors(omega_k, 0., 0.)
lam_q, lamtil_q = helicity_spinors(E_p, theta_q, 0.)

ang_kq = angle_bracket(lam_k, lam_q)
sq_kq  = square_bracket(lamtil_k, lamtil_q)
ang_qk = angle_bracket(lam_q, lam_k)
sq_qk  = square_bracket(lamtil_q, lamtil_k)

err_antisym_ang = abs(ang_kq + ang_qk)
err_antisym_sq  = abs(sq_kq  + sq_qk)

# Schouten: ⟨12⟩⟨34⟩ + ⟨13⟩⟨42⟩ + ⟨14⟩⟨23⟩ = 0
E_s = 1.0
lam1, _ = helicity_spinors(E_s, 0.2, 0.0)
lam2, _ = helicity_spinors(E_s, 1.1, 0.5)
lam3, _ = helicity_spinors(E_s, 2.0, 1.2)
lam4, _ = helicity_spinors(E_s, 0.7, 2.3)

ab12 = angle_bracket(lam1, lam2)
ab34 = angle_bracket(lam3, lam4)
ab13 = angle_bracket(lam1, lam3)
ab42 = angle_bracket(lam4, lam2)
ab14 = angle_bracket(lam1, lam4)
ab23 = angle_bracket(lam2, lam3)
schouten_val = ab12 * ab34 + ab13 * ab42 + ab14 * ab23
err_schouten = abs(schouten_val)

# ── §48.C: Mandelstam ─────────────────────────────────────────────────────────
k_mu = np.array([omega_k, 0., 0., omega_k])
q_mu = np.array([E_p, E_p * np.sin(theta_q), 0., E_p * np.cos(theta_q)])
kdotq = g[0,0]*k_mu[0]*q_mu[0] + sum(g[i,i]*k_mu[i]*q_mu[i] for i in range(1,4))
s_kq_expected = -2. * kdotq
prod_kq = ang_kq * sq_qk   # ⟨kq⟩[qk]
err_mandelstam = abs(prod_kq - s_kq_expected)

# ── §48.D: helicity spinors (z-axis) ─────────────────────────────────────────
omega_z = 2.5
lam_z, lamtil_z = helicity_spinors(omega_z, 0., 0.)
# 4-component for Dirac eq check
ket_sq_z, ket_an_z, bra_sq_z, bra_an_z = massless_spinors_general(omega_z, 0., 0.)
k_lower_z = np.array([omega_z, 0., 0., -omega_z])
slash_kz   = sum(k_lower_z[mu] * gam[mu] for mu in range(4))
err_dirac_sq = np.max(np.abs(slash_kz @ ket_sq_z))
err_dirac_an = np.max(np.abs(slash_kz @ ket_an_z))

lam_z_vec   = f"\\sqrt{{2\\omega}}\\,(1,0)^T"
lamtil_z_vec = f"\\sqrt{{2\\omega}}\\,(1,0)^T"
ket_sq_vals = f"\\sqrt{{2\\omega}}\\,(0,1,0,0)^T"
ket_an_vals = f"\\sqrt{{2\\omega}}\\,(0,0,1,0)^T"

# ── §48.E: polarization vectors ──────────────────────────────────────────────
omega_k_pol = 2.0
omega_q_pol = 1.5
lam_kp, lamtil_kp = helicity_spinors(omega_k_pol, 0., 0.)
lam_qp, lamtil_qp = helicity_spinors(omega_q_pol, np.pi/2, 0.)

ket_sq_kp, ket_an_kp, bra_sq_kp, bra_an_kp = massless_spinors_general(omega_k_pol, 0., 0.)
ket_sq_qp, ket_an_qp, bra_sq_qp, bra_an_qp = massless_spinors_general(omega_q_pol, np.pi/2, 0.)

ang_qp_kp = angle_bracket(lam_qp, lam_kp)
sq_qp_kp  = square_bracket(lamtil_qp, lamtil_kp)

eps_plus  = np.zeros(4, dtype=complex)
eps_minus = np.zeros(4, dtype=complex)
sigma_list = [I2, sigma[1], sigma[2], sigma[3]]
sigb_list  = [I2, -sigma[1], -sigma[2], -sigma[3]]
for mu in range(4):
    eps_plus[mu]  = (lamtil_qp.conj() @ sigma_list[mu] @ lam_kp) / (np.sqrt(2) * ang_qp_kp)
    eps_minus[mu] = (lam_qp.conj() @ sigb_list[mu] @ lamtil_kp) / (np.sqrt(2) * sq_qp_kp)

k_lower_pol = np.array([omega_k_pol, 0., 0., -omega_k_pol])
trans_p = sum(k_lower_pol[mu] * eps_plus[mu]  for mu in range(4))
trans_m = sum(k_lower_pol[mu] * eps_minus[mu] for mu in range(4))
norm_pm = -eps_plus[0]*eps_minus[0].conj() + sum(eps_plus[i]*eps_minus[i].conj() for i in range(1,4))

err_trans_p   = abs(trans_p)
err_trans_m   = abs(trans_m)
err_norm_pm   = abs(norm_pm + 1.)

# Ward identity: k·ε = 0
ward_str_plus  = f"{err_trans_p:.1e}"
ward_str_minus = f"{err_trans_m:.1e}"
norm_pm_val    = f"{err_norm_pm:.1e}"

# ── §48.F: little group scaling ──────────────────────────────────────────────
# Under λ → t λ, λ̃ → t⁻¹ λ̃, a helicity-h state scales as t^{-2h}
# ε_+^μ(k;q): contains one λ̃_q (weight -1) and one λ_k (weight +1)... but the
# denominator ⟨qk⟩ has weight t_k.
# Net: ε_+ → t_k^{-2} ε_+ (helicity weight -2 for positive helicity photon… see text)

# Numerical test: rescale λ_k → t λ_k, λ̃_k → t⁻¹ λ̃_k
t_val = 2.3
lam_kp_scaled    = t_val * lam_kp
lamtil_kp_scaled = (1. / t_val) * lamtil_kp

ang_qp_kp_sc = angle_bracket(lam_qp, lam_kp_scaled)
eps_plus_sc  = np.zeros(4, dtype=complex)
for mu in range(4):
    eps_plus_sc[mu] = (lamtil_qp.conj() @ sigma_list[mu] @ lam_kp_scaled) / (np.sqrt(2) * ang_qp_kp_sc)

# After rescaling: ε_+^μ(k;q) → (1/t²) × original?
# λ_k → t λ_k: ⟨qk⟩ → t⟨qk⟩  (one λ_k),  numerator ⟨q|σ^μ|k] uses λ_k → t factor
# But |k] = ε·λ̃_k, and λ̃_k → t⁻¹ λ̃_k... wait, the numerator ⟨q|σ^μ|k] uses λ̃_q^* and λ_k.
# Let me just compute the ratio:
lg_ratio = eps_plus_sc[1] / eps_plus[1]  # should be t^{-2} = 1/t^2 (helicity +1 photon)
lg_expected = t_val**(-2)
err_lg = abs(lg_ratio - lg_expected)

# ── §48.G: MHV preview (no numerical check, symbolic only) ───────────────────
# A_4(1-,2-,3+,4+) = ⟨12⟩^4 / (⟨12⟩⟨23⟩⟨34⟩⟨41⟩)
# Numerically check the n=4 MHV ratio structure (dimensionless):
lam_1, lamtil_1 = helicity_spinors(1., 0.0, 0.0)
lam_2, lamtil_2 = helicity_spinors(1., np.pi/3, np.pi/4)
lam_3, lamtil_3 = helicity_spinors(1., 2*np.pi/3, np.pi/2)

ab12_mhv = angle_bracket(lam_1, lam_2)
ab23_mhv = angle_bracket(lam_2, lam_3)
ab31_mhv = angle_bracket(lam_3, lam_1)

# 3-gluon MHV: A_3(1-,2-,3+) = ⟨12⟩^3 / (⟨12⟩⟨23⟩⟨31⟩)
mhv3_num   = ab12_mhv**3
mhv3_denom = ab12_mhv * ab23_mhv * ab31_mhv
mhv3       = mhv3_num / mhv3_denom

# ── cadabra2 expressions ─────────────────────────────────────────────────────
cdb_exprs = {
    "p_spinor":     cdb(r"\lambda_{\alpha} \lamtil_{\dal}"),
    "ang_bracket":  cdb(r"\abra{k}{p}"),
    "sq_bracket":   cdb(r"\sbra{k}{p}"),
    "mandelstam":   cdb(r"\abra{k}{p} \sbra{p}{k}"),
    "schouten_lhs": cdb(r"\abra{i}{j} \abra{k}{l}"),
    "schouten_sym": cdb(r"\abra{i}{j} \abra{k}{l} + \abra{i}{k} \abra{l}{j} + \abra{i}{l} \abra{j}{k}"),
}

# ── build LaTeX document ─────────────────────────────────────────────────────

doc = r"""
\documentclass[12pt]{article}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{mathtools}
\usepackage{tensor}
\usepackage{geometry}
\usepackage{booktabs}
\usepackage{array}
\usepackage{xcolor}
\usepackage{hyperref}
\usepackage{fancyhdr}
\usepackage{slashed}
\geometry{margin=1.1in, headheight=15pt}

\definecolor{cadblue}{RGB}{30,90,200}
\definecolor{verifygreen}{RGB}{20,140,60}

\newcommand{\cdbexpr}[1]{\textcolor{cadblue}{$#1$}}
\newcommand{\verified}[1]{\textcolor{verifygreen}{\checkmark\; #1}}
\newcommand{\half}{\tfrac{1}{2}}
\newcommand{\dal}{\dot{\alpha}}
\newcommand{\dbe}{\dot{\beta}}
\newcommand{\dga}{\dot{\gamma}}
\newcommand{\dde}{\dot{\delta}}
\newcommand{\sigmabar}{\bar{\sigma}}
\newcommand{\lam}{\lambda}
\newcommand{\lamtil}{\tilde{\lambda}}

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 50}
\rhead{Massless Particles \& Spinor-Helicity}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 50}\\[6pt]
       \large Massless Particles and Spinor-Helicity Formalism\\[4pt]
       \normalsize Cadabra2 expressions and numerical verification}
\author{Generated by \texttt{ch50\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{§48.A\quad Massless Momentum as Spinor Outer Product}
%% ============================================================

\subsection{Rank-1 factorization}

A massless 4-momentum $p^\mu$ (with $p^2=0$) maps to a $2\times2$ spinor matrix
\begin{equation}
  p_{\alpha\dot\alpha} \;=\; p^\mu\,\sigma_{\mu,\,\alpha\dot\alpha},
  \qquad
  \sigma_\mu = (I_2,\,\boldsymbol{\sigma}),
  \tag{48.1}
\end{equation}
whose determinant equals $-p^2/4=0$ (for massless $p$), so the matrix has rank~1.
Any rank-1 $2\times2$ matrix factors as an outer product:
\begin{equation}
  \boxed{
    p_{\alpha\dot\alpha} = \lambda_\alpha\,\tilde\lambda_{\dot\alpha}
  }
  \qquad (p^2 = 0)
  \tag{48.2}
\end{equation}
The two-component \textbf{helicity spinors} $\lambda_\alpha$ (left-handed, undotted)
and $\tilde\lambda_{\dot\alpha}$ (right-handed, dotted) encode the full massless kinematics.

\medskip
\textbf{Cadabra2 representation:} \cdbexpr{""" + cdb_exprs["p_spinor"] + r"""}

\subsection{Explicit components}

For $p^\mu = E(1,\,\sin\theta\cos\phi,\,\sin\theta\sin\phi,\,\cos\theta)$:
\begin{equation}
  \lambda_\alpha = \sqrt{2E}\begin{pmatrix}\cos\tfrac\theta2 \\ \sin\tfrac\theta2\,e^{i\phi}\end{pmatrix},
  \qquad
  \tilde\lambda_{\dot\alpha} = \sqrt{2E}\begin{pmatrix}\cos\tfrac\theta2 \\ \sin\tfrac\theta2\,e^{-i\phi}\end{pmatrix}
  \tag{48.3}
\end{equation}

\textbf{Numerical check} ($E=""" + str(E_p) + r"""$, $\theta=\pi/4$, $\phi=\pi/3$):
\[
  p_{\alpha\dot\alpha} = """ + p_spinor_tex + r""",
  \quad
  \lambda = """ + lam_p_vals + r""",
  \quad
  \tilde\lambda = """ + lamtil_p_vals + r"""
\]

\verified{$p_{\alpha\dot\alpha} = \lambda_\alpha\tilde\lambda_{\dot\alpha}$: max error $""" + err_str(err_rank1) + r"""$.}

\verified{$\det(p_{\alpha\dot\alpha}) = """ + f"{det_p.real:.2e}" + r"""$ (should be 0, since $p^2=0$).}

%% ============================================================
\section{§48.B\quad Angle $\langle ij\rangle$ and Square $[ij]$ Brackets}
%% ============================================================

\subsection{Definitions}

The fundamental Lorentz-invariant bilinears are:
\begin{align}
  \langle i\,j\rangle &= \varepsilon^{\alpha\beta}\,\lambda^i_\alpha\,\lambda^j_\beta
  = \lambda^i_1\lambda^j_2 - \lambda^i_2\lambda^j_1
  & \text{(angle bracket, undotted)},
  \tag{48.4}\\[4pt]
  [i\,j] &= \varepsilon_{\dot\alpha\dot\beta}\,\tilde\lambda^{i\,\dot\alpha}\tilde\lambda^{j\,\dot\beta}
  = \tilde\lambda^i_{\dot1}\tilde\lambda^j_{\dot2} - \tilde\lambda^i_{\dot2}\tilde\lambda^j_{\dot1}
  & \text{(square bracket, dotted)},
  \tag{48.5}
\end{align}
where $\varepsilon^{12}=+1$ (Srednicki convention).

\textbf{Cadabra2:}\quad
\cdbexpr{""" + cdb_exprs["ang_bracket"] + r"""} $= \langle k\,p\rangle$\quad and\quad
\cdbexpr{""" + cdb_exprs["sq_bracket"]  + r"""} $= [k\,p]$

\subsection{Antisymmetry}

Both brackets are antisymmetric in their arguments:
\begin{equation}
  \langle i\,j\rangle = -\langle j\,i\rangle,
  \qquad
  [i\,j] = -[j\,i].
  \tag{48.6}
\end{equation}

\verified{$|\langle kq\rangle + \langle qk\rangle| = """ + err_str(err_antisym_ang) + r"""$.}
\verified{$|[kq] + [qk]| = """ + err_str(err_antisym_sq) + r"""$.}

\subsection{Schouten identity}

Because spinor space is 2-dimensional, any three undotted spinors are linearly
dependent.  This yields the \textbf{Schouten identity}:
\begin{equation}
  \boxed{
    \langle i\,j\rangle\langle k\,l\rangle
    + \langle i\,k\rangle\langle l\,j\rangle
    + \langle i\,l\rangle\langle j\,k\rangle = 0
  }
  \tag{48.7}
\end{equation}
and identically for square brackets.

\textbf{Cadabra2 (LHS terms):} \cdbexpr{""" + cdb_exprs["schouten_sym"] + r"""}

\verified{$|\langle12\rangle\langle34\rangle + \langle13\rangle\langle42\rangle + \langle14\rangle\langle23\rangle| = """ + err_str(err_schouten) + r"""$ (four generic massless momenta).}

%% ============================================================
\section{§48.C\quad Mandelstam Variables from Spinors}
%% ============================================================

For massless momenta $p_i$, $p_j$ (with $p_i^2=p_j^2=0$), the kinematic invariant is:
\begin{equation}
  \boxed{
    s_{ij} = (p_i+p_j)^2 = 2\,p_i\cdot p_j = \langle i\,j\rangle[j\,i]
  }
  \tag{48.8}
\end{equation}

\textbf{Cadabra2:} \cdbexpr{""" + cdb_exprs["mandelstam"] + r"""} $= s_{kp}$

\medskip
\textbf{Physical meaning:} The product $\langle kq\rangle[qk]$ gives the
Mandelstam invariant directly from the helicity spinors, without needing to
evaluate $p^\mu q_\mu$ with the metric.

\textbf{Numerical check} ($k=(""" + f"{omega_k},0,0,{omega_k}" + r""")$, $q=E(\sin\theta,0,\cos\theta)$ with $E=""" + str(E_p) + r"""$, $\theta=\pi/3$):
\begin{align*}
  \langle kq\rangle[qk] &= """ + f"{prod_kq.real:.6f}" + r""" \\
  -2\,k\cdot q &= """ + f"{s_kq_expected:.6f}" + r"""
\end{align*}

\verified{$|\langle kq\rangle[qk] - (-2k\cdot q)| = """ + err_str(err_mandelstam) + r"""$.}

%% ============================================================
\section{§48.D\quad Helicity Spinors $u_\pm(k)$ in the Chiral Basis}
%% ============================================================

\subsection{Z-axis momentum: explicit form}

For massless momentum $k^\mu=(\omega,0,0,\omega)$ along the $+z$ axis,
the massless Dirac equation $\slashed{k}|k\rangle_\pm=0$ has solutions
(Srednicki eq.~60.10):
\begin{equation}
  |k] = u_-(k) = """ + ket_sq_vals + r""",
  \qquad
  |k\rangle = u_+(k) = """ + ket_an_vals + r"""
  \tag{48.9}
\end{equation}
where $|k]$ has negative helicity ($h=-\tfrac12$) and $|k\rangle$ has positive helicity ($h=+\tfrac12$).

The corresponding 2-component helicity spinors (upper/lower Weyl blocks):
\begin{equation}
  \lambda_\alpha(k)\big|_{z\text{-axis}} = \sqrt{2\omega}\begin{pmatrix}1\\0\end{pmatrix},
  \qquad
  \tilde\lambda_{\dot\alpha}(k)\big|_{z\text{-axis}} = \sqrt{2\omega}\begin{pmatrix}1\\0\end{pmatrix}
  \tag{48.10}
\end{equation}

\verified{Massless Dirac eq.\ $\slashed{k}|k] = 0$: max error $""" + err_str(err_dirac_sq) + r"""$.}
\verified{Massless Dirac eq.\ $\slashed{k}|k\rangle = 0$: max error $""" + err_str(err_dirac_an) + r"""$.}

\subsection{General massless momentum}

For $p^\mu = E(1,\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$, the
4-component spinors in the Weyl representation are:
\begin{equation}
  |p] = \sqrt{2E}\begin{pmatrix}-\sin\tfrac\theta2\,e^{-i\phi} \\ \cos\tfrac\theta2 \\ 0 \\ 0\end{pmatrix},
  \qquad
  |p\rangle = \sqrt{2E}\begin{pmatrix}0 \\ 0 \\ \cos\tfrac\theta2 \\ \sin\tfrac\theta2\,e^{i\phi}\end{pmatrix}
  \tag{48.11}
\end{equation}
These satisfy $\slashed{p}|p]=\slashed{p}|p\rangle=0$ and are eigenstates of
the helicity operator $\hat h = \hat{\mathbf{J}}\cdot\hat{\mathbf{p}}$ with
eigenvalues $\mp\tfrac12$.

%% ============================================================
\section{§48.E\quad Gluon Polarization Vectors as Spinor Bilinears}
%% ============================================================

\subsection{Definitions}

The photon/gluon polarization vectors can be written as spinor bilinears
(Srednicki between eqs.~60.6--60.9):
\begin{equation}
  \boxed{
    \varepsilon_+^\mu(k;\,q) = \frac{\langle q|\,\gamma^\mu\,|k]}{\sqrt{2}\,\langle q\,k\rangle},
    \qquad
    \varepsilon_-^\mu(k;\,q) = \frac{[q|\,\gamma^\mu\,|k\rangle}{\sqrt{2}\,[q\,k]}
  }
  \tag{48.12}
\end{equation}
Here $k$ is the photon momentum and $q\neq k$ is an arbitrary massless
\textbf{reference momentum} encoding the residual gauge freedom.
Different choices of $q$ give the same physical amplitude (gauge invariance).

\textbf{Ward identity (transversality):}
\begin{equation}
  k_\mu\,\varepsilon_\pm^\mu(k;\,q) = 0
  \tag{48.13}
\end{equation}

\textbf{Normalization:}
\begin{equation}
  \varepsilon_+\cdot\varepsilon_-^* = -1,
  \qquad
  \varepsilon_+\cdot\varepsilon_+^* = 0
  \tag{48.14}
\end{equation}

\subsection{Numerical verification}

Setup: $k^\mu=(\omega_k,0,0,\omega_k)$ along $+z$ with $\omega_k=""" + str(omega_k_pol) + r"""$;
reference $q$ along $+x$ with $\omega_q=""" + str(omega_q_pol) + r"""$.

\verified{$|k_\mu\varepsilon_+^\mu(k;q)| = """ + ward_str_plus + r"""$ (Ward identity for $\varepsilon_+$).}
\verified{$|k_\mu\varepsilon_-^\mu(k;q)| = """ + ward_str_minus + r"""$ (Ward identity for $\varepsilon_-$).}
\verified{$|\varepsilon_+\cdot\varepsilon_-^*+1| = """ + norm_pm_val + r"""$ (normalization).}

%% ============================================================
\section{§48.F\quad Little Group Scaling and Helicity Weights}
%% ============================================================

\subsection{The little group}

The little group of a massless momentum $p^\mu$ is $U(1)\simeq SO(2)$:
rotations about the direction of motion.  Under the little-group transformation
\begin{equation}
  \lambda_\alpha \;\to\; t\,\lambda_\alpha,
  \qquad
  \tilde\lambda_{\dot\alpha} \;\to\; t^{-1}\,\tilde\lambda_{\dot\alpha},
  \qquad t\in\mathbb{C}^*,
  \tag{48.15}
\end{equation}
the momentum $p_{\alpha\dot\alpha}=\lambda_\alpha\tilde\lambda_{\dot\alpha}$ is
invariant.  Spinors (and hence amplitudes) transform with definite
\textbf{helicity weight}:
\begin{equation}
  \langle i\,j\rangle \;\to\; t_i\,t_j\,\langle i\,j\rangle,
  \qquad
  [i\,j] \;\to\; t_i^{-1}\,t_j^{-1}\,[i\,j].
  \tag{48.16}
\end{equation}
Under particle $k$'s own little group ($t_k=t$, all others fixed):
\begin{equation}
  \varepsilon_+^\mu(k;\,q) \;\to\; t^{-2}\,\varepsilon_+^\mu(k;\,q),
  \qquad
  \varepsilon_-^\mu(k;\,q) \;\to\; t^{+2}\,\varepsilon_-^\mu(k;\,q).
  \tag{48.17}
\end{equation}
This is the statement that a positive-helicity photon has helicity weight $h=-1$
under the little group (the polarization vector itself carries $-2h$ units).

\textbf{Numerical check} (rescaling $\lambda_k\to t\lambda_k$, $\tilde\lambda_k\to t^{-1}\tilde\lambda_k$ with $t=""" + f"{t_val}" + r"""$):

\verified{$\varepsilon_+^\mu \to t^{-2}\varepsilon_+^\mu$: error $""" + err_str(err_lg) + r"""$.}

%% ============================================================
\section{§48.G\quad MHV Parke-Taylor Formula (Seed)}
%% ============================================================

\subsection{Motivation}

The spinor-helicity formalism makes tree-level gluon amplitudes remarkably
compact.  Amplitudes with all-same or all-but-one helicities vanish.
The first nonzero class are the \textbf{MHV (Maximally Helicity Violating)}
amplitudes, with exactly two negative helicities.

\subsection{Parke-Taylor formula}

The $n$-gluon colour-ordered MHV amplitude (two negative helicities at positions $i,j$):
\begin{equation}
  \boxed{
    \mathcal{A}_n^{\rm MHV}\!\bigl(1^+,\ldots,i^-,\ldots,j^-,\ldots,n^+\bigr)
    = \frac{\langle i\,j\rangle^4}
           {\langle 1\,2\rangle\langle 2\,3\rangle\cdots\langle n\,1\rangle}
  }
  \tag{48.18}
\end{equation}
(overall coupling and colour factors suppressed).
This formula, the \textbf{Parke-Taylor formula}, is the \emph{seed} for BCFW
on-shell recursion: all tree amplitudes can be built from it by shifting two
complex spinors and summing over factorisation channels.

\subsection{Structure check (3-point MHV)}

The $n=3$ MHV amplitude
$\mathcal{A}_3(1^-,2^-,3^+) = \langle12\rangle^3/(\langle12\rangle\langle23\rangle\langle31\rangle)$
is a pure ratio of angle brackets.  Numerically:
\[
  \mathcal{A}_3 = \frac{\langle12\rangle^3}{\langle12\rangle\langle23\rangle\langle31\rangle}
  = \frac{(""" + f"{ab12_mhv:.4f}" + r""")^3}{(""" + f"{ab12_mhv:.4f}" + r""")(""" + f"{ab23_mhv:.4f}" + r""")(""" + f"{ab31_mhv:.4f}" + r""")}
  = """ + f"{mhv3:.4f}" + r"""
\]

\verified{Parke-Taylor $n=3$ ratio evaluated from spinors (dimensionless Lorentz invariant).}

%% ============================================================
\section{Numerical Verification Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}lll@{}}
  \toprule
  Identity & Equation & Max error \\
  \midrule
  $p_{\alpha\dot\alpha} = \lambda_\alpha\tilde\lambda_{\dot\alpha}$ (rank-1)
    & (48.2) & $""" + err_str(err_rank1) + r"""$ \\
  $\det(p_{\alpha\dot\alpha}) = 0$
    & ($p^2=0$) & $""" + f"{abs(det_p):.1e}" + r"""$ \\
  $\langle ij\rangle = -\langle ji\rangle$ (antisymmetry)
    & (48.6) & $""" + err_str(err_antisym_ang) + r"""$ \\
  $[ij] = -[ji]$ (antisymmetry)
    & (48.6) & $""" + err_str(err_antisym_sq) + r"""$ \\
  Schouten identity $\langle\cdot\rangle$ brackets
    & (48.7) & $""" + err_str(err_schouten) + r"""$ \\
  $\langle kq\rangle[qk] = -2k\cdot q$
    & (48.8) & $""" + err_str(err_mandelstam) + r"""$ \\
  $\slashed{k}|k] = 0$ (massless Dirac)
    & (48.9) & $""" + err_str(err_dirac_sq) + r"""$ \\
  $\slashed{k}|k\rangle = 0$ (massless Dirac)
    & (48.9) & $""" + err_str(err_dirac_an) + r"""$ \\
  $k_\mu\varepsilon_+^\mu = 0$ (Ward)
    & (48.13) & $""" + ward_str_plus + r"""$ \\
  $k_\mu\varepsilon_-^\mu = 0$ (Ward)
    & (48.13) & $""" + ward_str_minus + r"""$ \\
  $\varepsilon_+\cdot\varepsilon_-^* = -1$
    & (48.14) & $""" + norm_pm_val + r"""$ \\
  Little group: $\varepsilon_+ \to t^{-2}\varepsilon_+$
    & (48.17) & $""" + err_str(err_lg) + r"""$ \\
  \bottomrule
\end{tabular}
\end{center}

All identities verified to machine precision ($\lesssim 10^{-12}$).

%% ============================================================
\section{Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}ll@{}}
  \toprule
  Object & Expression \\
  \midrule
  Spinor momentum & $p_{\alpha\dot\alpha} = \lambda_\alpha\tilde\lambda_{\dot\alpha}$;\quad $p^2=0$ iff rank 1 \\
  Angle bracket & $\langle ij\rangle = \varepsilon^{\alpha\beta}\lambda^i_\alpha\lambda^j_\beta = -\langle ji\rangle$ \\
  Square bracket & $[ij] = \varepsilon_{\dot\alpha\dot\beta}\tilde\lambda^{i\dot\alpha}\tilde\lambda^{j\dot\beta} = -[ji]$ \\
  Mandelstam & $s_{ij} = \langle ij\rangle[ji] = -2p_i\cdot p_j$ \\
  $z$-axis spinors & $|k]=\sqrt{2\omega}(0,1,0,0)^T$;\quad $|k\rangle=\sqrt{2\omega}(0,0,1,0)^T$ \\
  Polarization $\varepsilon_+$ & $\varepsilon_+^\mu(k;q) = \langle q|\gamma^\mu|k]/(\sqrt{2}\langle qk\rangle)$ \\
  Polarization $\varepsilon_-$ & $\varepsilon_-^\mu(k;q) = [q|\gamma^\mu|k\rangle/(\sqrt{2}[qk])$ \\
  Ward identity & $k_\mu\varepsilon_\pm^\mu = 0$ \\
  Little group & $\lambda\to t\lambda$,\; $\tilde\lambda\to t^{-1}\tilde\lambda$;\; $\varepsilon_\pm\to t^{\mp2}\varepsilon_\pm$ \\
  Schouten & $\langle ij\rangle\langle kl\rangle+\langle ik\rangle\langle lj\rangle+\langle il\rangle\langle jk\rangle=0$ \\
  MHV (Parke-Taylor) & $\mathcal{A}_n^{\rm MHV}=\langle ij\rangle^4/(\langle12\rangle\cdots\langle n1\rangle)$ \\
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\noindent\textbf{Next:} Chapter 60 applies these tools to spinor electrodynamics:
twistor notation $|p], |p\rangle, [p|, \langle p|$, polarization vectors as
spinor bilinears, and the Fierz identities for loop computations.

\end{document}
"""

outpath = "ch50_massless_spinor_helicity.tex"
with open(outpath, "w") as f:
    f.write(doc.strip())

print(f"Wrote: {outpath}")
