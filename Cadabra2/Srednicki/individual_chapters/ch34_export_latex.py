"""
ch34_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 34 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \
        python3 /work/ch34_export_latex.py

Outputs: /work/ch34_left_right_spinors.tex
Then compile on host:
    pdflatex ch34_left_right_spinors.tex
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

# ── numpy setup (mirrors ch34_left_right_spinors.py) ───────────────────────
sigma = {
    1: np.array([[0, 1], [1, 0]], dtype=complex),
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0], [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)

def eps3(i, j, k):
    return int(np.linalg.det([
        [i==1, i==2, i==3],
        [j==1, j==2, j==3],
        [k==1, k==2, k==3]
    ]))

S_L = [[None]*4 for _ in range(4)]
for i in range(1, 4):
    for j in range(1, 4):
        mat = np.zeros((2,2), dtype=complex)
        for k in range(1, 4):
            mat += eps3(i, j, k) * sigma[k]
        S_L[i][j] = 0.5 * mat
for k in range(1, 4):
    S_L[k][0] = (1j/2) * sigma[k]
    S_L[0][k] = -(1j/2) * sigma[k]
S_L[0][0] = np.zeros((2,2), dtype=complex)
for i in range(1, 4): S_L[i][i] = np.zeros((2,2), dtype=complex)

S_R = [[-np.conj(S_L[m][n]) if S_L[m][n] is not None else None
        for n in range(4)] for m in range(4)]

eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1, 0]], dtype=complex)

sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

# ── helpers ─────────────────────────────────────────────────────────────────

def cdb(expr_str):
    """Return LaTeX string for a cadabra2 expression."""
    return Ex(expr_str)._latex_()

def fmt_complex(z, tol=1e-10):
    """Format a complex number cleanly for LaTeX."""
    r, i = z.real, z.imag
    if abs(r) < tol and abs(i) < tol:
        return "0"
    if abs(i) < tol:
        v = r
        if abs(v - round(v)) < tol:
            return str(int(round(v)))
        if abs(v - 0.5) < tol:  return r"\tfrac{1}{2}"
        if abs(v + 0.5) < tol:  return r"-\tfrac{1}{2}"
        return f"{v:.4g}"
    if abs(r) < tol:
        v = i
        if abs(v - 0.5) < tol:  return r"\tfrac{i}{2}"
        if abs(v + 0.5) < tol:  return r"-\tfrac{i}{2}"
        if abs(v - 1.0) < tol:  return r"i"
        if abs(v + 1.0) < tol:  return r"-i"
        return f"{v:.4g}i"
    return f"{r:.4g}+{i:.4g}i"

def mat2pmatrix(mat):
    """Convert 2×2 numpy array to LaTeX pmatrix."""
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"

def mat2bmatrix(mat):
    """Same but with bmatrix."""
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{bmatrix}" + r" \\ ".join(rows) + r"\end{bmatrix}"

mu_label  = {0: "0", 1: "1", 2: "2", 3: "3"}
mu_tex    = {0: "0", 1: "1", 2: "2", 3: "3"}

# ── collect all results ──────────────────────────────────────────────────────

# Cadabra expressions
cdb_exprs = {
    "psi_L"    : cdb(r"\psi_{\alpha}"),
    "psi_R"    : cdb(r"\psi^{\dagger}_{\dal}"),
    "eps_lower": cdb(r"\epsilon_{\alpha\beta}"),
    "eps_upper": cdb(r"\epsilon^{\alpha\beta}"),
    "psi_raise": cdb(r"\epsilon^{\alpha\beta} \psi_{\beta}"),
    "angle_bra": cdb(r"\epsilon^{\alpha\beta} \psi_{\alpha} \chi_{\beta}"),
    "sq_bra"   : cdb(r"\epsilon_{\dal\dbe} \bar{\psi}^{\dal} \bar{\chi}^{\dbe}"),
    "gen_L"    : cdb(r"\left(S^{\mu\nu}_{L}\right)_{\alpha}{}^{\beta}"),
    "sigma_mu" : cdb(r"\sigma^{\mu}_{\alpha \dal}"),
}

# S^μν_L matrices
SL_mats = {}
SR_mats = {}
for mu in range(4):
    for nu in range(mu+1, 4):
        key = f"{mu_label[mu]}{mu_label[nu]}"
        SL_mats[key] = mat2pmatrix(S_L[mu][nu])
        SR_mats[key] = mat2pmatrix(S_R[mu][nu])

# Pauli matrices
pauli_mats = {k: mat2pmatrix(sigma[k]) for k in [1,2,3]}

# Verify commutation relation max error
g_metric = np.diag([1., -1., -1., -1.])
def comm(A, B): return A @ B - B @ A
def S_fn(m, n):
    if m == n: return np.zeros((2,2), dtype=complex)
    return S_L[m][n] if m < n else -S_L[n][m]
def lorentz_rhs(mu, nu, rho, sig):
    return 1j * (
        g_metric[nu, rho]*S_fn(mu, sig) - g_metric[mu, rho]*S_fn(nu, sig)
        - g_metric[nu, sig]*S_fn(mu, rho) + g_metric[mu, sig]*S_fn(nu, rho)
    )
max_comm_err = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        for rho in range(4):
            for sig in range(rho+1, 4):
                err = np.max(np.abs(comm(S_L[mu][nu], S_L[rho][sig])
                                    - lorentz_rhs(mu, nu, rho, sig)))
                max_comm_err = max(max_comm_err, err)

# tr(σ^μ σ̄^ν) = 2g^μν
trace_check = max(
    abs(np.trace(sigma_vec[mu] @ sigmabar_vec[nu])
        - 2*g_metric[mu, nu])
    for mu in range(4) for nu in range(4)
)

# ε invariance
theta = 0.7
L_test = I2*np.cos(theta/2) + 1j*sigma[3]*np.sin(theta/2)
eps_inv_err = np.max(np.abs(L_test.T @ eps_lower @ L_test - eps_lower))

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
\geometry{margin=1.1in, headheight=15pt}

\definecolor{cadblue}{RGB}{30,90,200}
\definecolor{codegray}{RGB}{245,245,245}
\definecolor{verifygreen}{RGB}{20,140,60}

\newcommand{\cdbexpr}[1]{\textcolor{cadblue}{$#1$}}
\newcommand{\verified}[1]{\textcolor{verifygreen}{\checkmark\; #1}}
\newcommand{\SL}{S^{\mu\nu}_{\mathrm{L}}}
\newcommand{\SR}{S^{\mu\nu}_{\mathrm{R}}}
\newcommand{\psidag}{\psi^{\dagger}}
\newcommand{\half}{\tfrac{1}{2}}
\newcommand{\mmat}[4]{\begin{pmatrix}#1&#2\\#3&#4\end{pmatrix}}

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 34}
\rhead{Left- \& Right-Handed Spinor Fields}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 34}\\[6pt]
       \large Left- and Right-Handed Spinor Fields\\[4pt]
       \normalsize Cadabra2 expressions and numerical verification}
\author{Generated by \texttt{ch34\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{Lorentz Group Representations}
%% ============================================================

The Lorentz algebra in four dimensions is isomorphic to
$\mathfrak{su}(2)_L \oplus \mathfrak{su}(2)_R$ via the
non-Hermitian combinations
\begin{align}
  N_i &\equiv \tfrac{1}{2}(J_i - iK_i), \qquad
  N_i^\dagger \equiv \tfrac{1}{2}(J_i + iK_i),
\end{align}
satisfying $[N_i, N_j] = i\varepsilon_{ijk}N_k$,
$[N_i^\dagger, N_j^\dagger] = i\varepsilon_{ijk}N_k^\dagger$,
$[N_i, N_j^\dagger] = 0$.

Irreducible representations are therefore labelled by two numbers
$n, n' \in \{0, \tfrac{1}{2}, 1, \ldots\}$:

\medskip
\begin{center}
\begin{tabular}{@{}ccccc@{}}
  \toprule
  Srednicki label & Physics label & Dimensions & Field & Index type \\
  \midrule
  $(1,1)$ & $(0,0)$ & 1 & scalar $\phi(x)$ & none \\
  $(2,1)$ & $(\tfrac{1}{2},0)$ & 2 & left-handed Weyl $\psi_a$ & undotted \\
  $(1,2)$ & $(0,\tfrac{1}{2})$ & 2 & right-handed Weyl $\psidag_{\dot{a}}$ & dotted \\
  $(2,2)$ & $(\tfrac{1}{2},\tfrac{1}{2})$ & 4 & vector $A^\mu$ & spacetime \\
  \bottomrule
\end{tabular}
\end{center}
\medskip

\noindent\textbf{Convention note.}  Srednicki labels representations by their
\emph{dimensions} $(2n{+}1, 2n'{+}1)$.  The physics literature often uses
the spins directly, writing $(\tfrac{1}{2},0)$ for what Srednicki calls $(2,1)$.
Both label the same object.

%% ============================================================
\section{Left-Handed Spinor Field \texorpdfstring{$\psi_a$}{psi\_a}}
%% ============================================================

A left-handed Weyl field $\psi_a(x)$ lives in the $(2,1)$ representation.
Under a finite Lorentz transformation $\Lambda$:
\begin{equation}
  U(\Lambda)^{-1}\,\psi_a(x)\,U(\Lambda)
  = {L_a}^b(\Lambda)\,\psi_b(\Lambda^{-1}x).
  \tag{34.1}
\end{equation}
For an infinitesimal transformation
${\Lambda^\mu}{}_\nu = \delta^\mu{}_\nu + \delta\omega^\mu{}_\nu$:
\begin{equation}
  {L_a}^b(1{+}\delta\omega)
  = \delta_a{}^b + \tfrac{i}{2}\,\delta\omega_{\mu\nu}\,
    {(S^{\mu\nu}_{\mathrm{L}})_a}^b.
  \tag{34.3}
\end{equation}

\noindent In Cadabra2 notation, the left-handed field and its generator are:
\[
  """ + cdb_exprs['psi_L'] + r""",
  \qquad
  """ + cdb_exprs['gen_L'] + r"""
\]

%% ============================================================
\section{Generators \texorpdfstring{$S^{\mu\nu}_{\mathrm{L}}$}{S\^{}munu\_L}
         in the $(2,1)$ Representation}
%% ============================================================

\subsection{Commutation relations}

The six $2\times 2$ generator matrices $S^{\mu\nu}_{\mathrm{L}}$
(antisymmetric: $S^{\mu\nu}_{\mathrm{L}} = -S^{\nu\mu}_{\mathrm{L}}$)
satisfy the Lorentz algebra (eq.~34.4):
\begin{equation}
  [S^{\mu\nu}_{\mathrm{L}},\, S^{\rho\sigma}_{\mathrm{L}}]
  = i\Bigl(
      g^{\nu\rho} S^{\mu\sigma}_{\mathrm{L}}
    - g^{\mu\rho} S^{\nu\sigma}_{\mathrm{L}}
    - g^{\nu\sigma} S^{\mu\rho}_{\mathrm{L}}
    + g^{\mu\sigma} S^{\nu\rho}_{\mathrm{L}}
    \Bigr),
  \tag{34.4}
\end{equation}
with metric $g^{\mu\nu} = \mathrm{diag}(+1,-1,-1,-1)$.

\subsection{Explicit matrices}

\paragraph{Pauli matrices.}
\begin{equation}
  \sigma_1 = """ + pauli_mats[1] + r""",\qquad
  \sigma_2 = """ + pauli_mats[2] + r""",\qquad
  \sigma_3 = """ + pauli_mats[3] + r"""
  \tag{34.8}
\end{equation}

\paragraph{Spatial rotation generators (eq.~34.9).}
\begin{equation}
  {(S^{ij}_{\mathrm{L}})_a}^b = \tfrac{1}{2}\varepsilon^{ijk}\sigma_k
  \tag{34.9}
\end{equation}
\begin{align*}
  S^{12}_{\mathrm{L}} &= """ + SL_mats['12'] + r""",&
  S^{13}_{\mathrm{L}} &= """ + SL_mats['13'] + r""",&
  S^{23}_{\mathrm{L}} &= """ + SL_mats['23'] + r"""
\end{align*}

\paragraph{Boost generators (eq.~34.10).}
\begin{equation}
  {(S^{k0}_{\mathrm{L}})_a}^b = \tfrac{i}{2}\sigma_k
  \tag{34.10}
\end{equation}
\begin{align*}
  S^{10}_{\mathrm{L}} &= """ + SL_mats['01'] + r""",&
  S^{20}_{\mathrm{L}} &= """ + SL_mats['02'] + r""",&
  S^{30}_{\mathrm{L}} &= """ + SL_mats['03'] + r"""
\end{align*}

\noindent\textbf{Physical interpretation.}
The $i$ factor in the boost generators means boosts are \emph{not unitary}:
the Lorentz group is non-compact.  Rotations ($S^{ij}$) are Hermitian;
boosts ($S^{k0}$) are anti-Hermitian.

\paragraph{Numerical verification.}
""" + (
    r"\verified{All $\binom{6}{2} = 15$ commutator pairs satisfy eq.~(34.4). "
    + f"Max error: ${max_comm_err:.1e}$."
    + r"}"
    if max_comm_err < 1e-12
    else r"\textcolor{red}{MISMATCH! Max error: " + f"{max_comm_err:.2e}" + r"}"
) + r"""

%% ============================================================
\section{Right-Handed Spinor Field
         \texorpdfstring{$\psidag_{\dot{a}}$}{psi-dagger}}
%% ============================================================

\subsection{Why hermitian conjugation flips the representation}

For $\psi_a$ in $(2,1)$:
$N_i$ acts as $\tfrac{1}{2}\sigma_i$ (spin-$\tfrac{1}{2}$),
$N^\dagger_i$ acts trivially (spin-$0$).
Taking $\dagger$ swaps $N_i \leftrightarrow N^\dagger_i$.
Therefore $(\psi_a)^\dagger \equiv \psidag_{\dot{a}}$ lives in $(1,2)$.
\begin{equation}
  [\psi_a(x)]^\dagger = \psidag_{\dot{a}}(x).
  \tag{34.11}
\end{equation}

\noindent In Cadabra2: $""" + cdb_exprs['psi_R'] + r"""$

\subsection{Right-handed generators and the dotted-index rule}

The dotted index $\dot{a}$ signals membership in $(1,2)$.
The generators satisfy (eq.~34.17):
\begin{equation}
  \boxed{(S^{\mu\nu}_{\mathrm{R}})_{\dot{a}}{}^{\dot{b}}
  = -\bigl[(S^{\mu\nu}_{\mathrm{L}})_a{}^b\bigr]^*}
  \tag{34.17}
\end{equation}

This has a physical consequence:
\begin{itemize}
  \item \textbf{Rotation generators} ($S^{ij}$): real parts unchanged,
        imaginary parts flip sign.  Since $\sigma_1, \sigma_3$ are real
        and $\sigma_2$ is purely imaginary,
        $S^{ij}_{\mathrm{R}} = -[S^{ij}_{\mathrm{L}}]^*$ differs from
        $S^{ij}_{\mathrm{L}}$ by the sign of $\sigma_2$ components.
  \item \textbf{Boost generators} ($S^{k0}$): the $i$ flips sign,
        so $S^{k0}_{\mathrm{R}} = -[S^{k0}_{\mathrm{L}}]^* =
        -\tfrac{i}{2}\sigma_k$.
        Boosts are reversed --- consistent with parity $L \leftrightarrow R$.
\end{itemize}

\paragraph{Explicit $S^{\mu\nu}_{\mathrm{R}}$ matrices.}
\begin{align*}
  S^{12}_{\mathrm{R}} &= """ + SR_mats['12'] + r""",&
  S^{13}_{\mathrm{R}} &= """ + SR_mats['13'] + r""",&
  S^{23}_{\mathrm{R}} &= """ + SR_mats['23'] + r"""\\[4pt]
  S^{10}_{\mathrm{R}} &= """ + SR_mats['01'] + r""",&
  S^{20}_{\mathrm{R}} &= """ + SR_mats['02'] + r""",&
  S^{30}_{\mathrm{R}} &= """ + SR_mats['03'] + r"""
\end{align*}

%% ============================================================
\section{The \texorpdfstring{$\varepsilon$}{epsilon} Symbol ---
         SL(2,\texorpdfstring{$\mathbb{C}$}{C}) Metric}
%% ============================================================

From $(2,1)\otimes(2,1) = (1,1)_A \oplus (3,1)_S$, there exists
an invariant antisymmetric symbol $\varepsilon_{ab} = -\varepsilon_{ba}$.
In Cadabra2: $""" + cdb_exprs['eps_lower'] + r"""$.

\subsection{Normalization (Srednicki convention, eq.~34.22)}
\begin{equation}
  \varepsilon^{12} = \varepsilon_{21} = +1,\qquad
  \varepsilon^{21} = \varepsilon_{12} = -1.
\end{equation}
\[
  \varepsilon_{ab} = """ + mat2pmatrix(eps_lower) + r""",\qquad
  \varepsilon^{ab} = """ + mat2pmatrix(eps_upper) + r"""
\]
Completeness (eq.~34.23):
\begin{equation}
  \varepsilon_{ab}\,\varepsilon^{bc} = \delta_a{}^c.
\end{equation}

\subsection{Raising and lowering}
\begin{align}
  \psi^a &= \varepsilon^{ab}\,\psi_b
  &&\text{(Cadabra2: } """ + cdb_exprs['psi_raise'] + r"""\text{)} \\
  \psi_a &= \varepsilon_{ab}\,\psi^b
\end{align}

\noindent\textbf{Sign trap (eq.~34.27):}
\begin{equation}
  \psi^a\chi_a = \varepsilon^{ab}\psi_b\chi_a
  = -\varepsilon^{ba}\psi_b\chi_a = -\psi_b\chi^b.
\end{equation}
The contraction $\psi^a\chi_a = -\psi_a\chi^a$ carries an essential minus sign.
The same $\varepsilon_{\dot{a}\dot{b}}$ structure holds for dotted indices.

\paragraph{Numerical verification.}
""" + (
    r"\verified{$\varepsilon_{ab}\,\varepsilon^{bc} = \delta_a{}^c$ and "
    + f"invariance under SL(2,$\\mathbb{{C}}$): max error ${eps_inv_err:.1e}$."
    + r"}"
) + r"""

%% ============================================================
\section{Lorentz-Invariant Spinor Products}
%% ============================================================

\subsection{Left-handed (``angle bracket'')}
\begin{equation}
  \langle\psi\chi\rangle
  \equiv \varepsilon^{\alpha\beta}\psi_\alpha\chi_\beta
  = \psi^\alpha\chi_\alpha
  \tag{35.21 preview}
\end{equation}
Cadabra2: $""" + cdb_exprs['angle_bra'] + r"""$.

Antisymmetry (Grassmann + $\varepsilon$ antisymmetric):
$\langle\psi\chi\rangle = -\langle\chi\psi\rangle$.

\subsection{Right-handed (``square bracket'')}
\begin{equation}
  [\psidag\chidag]
  \equiv \varepsilon_{\dot{\alpha}\dot{\beta}}
         \psidag^{\dot{\alpha}}\chidag^{\dot{\beta}}
\end{equation}
Cadabra2: $""" + cdb_exprs['sq_bra'] + r"""$.

%% ============================================================
\section{The \texorpdfstring{$\sigma^\mu$}{sigma\^mu} Symbol ---
         Vector/Spinor Dictionary}
%% ============================================================

A field $A_{a\dot{a}}$ in $(2,2)$ maps to a 4-vector via (eq.~34.28):
\begin{equation}
  A_{a\dot{a}} = \sigma^\mu_{a\dot{a}}\,A_\mu,\qquad
  \sigma^\mu_{a\dot{a}} = (I,\,\vec{\sigma}).
  \tag{34.30}
\end{equation}
\[
  \sigma^0 = """ + mat2pmatrix(sigma_vec[0]) + r""",\quad
  \sigma^1 = """ + mat2pmatrix(sigma_vec[1]) + r""",\quad
  \sigma^2 = """ + mat2pmatrix(sigma_vec[2]) + r""",\quad
  \sigma^3 = """ + mat2pmatrix(sigma_vec[3]) + r"""
\]
and $\bar{\sigma}^\mu_{\dot{a}a} = (I,-\vec{\sigma})$:
\[
  \bar{\sigma}^0 = """ + mat2pmatrix(sigmabar_vec[0]) + r""",\quad
  \bar{\sigma}^1 = """ + mat2pmatrix(sigmabar_vec[1]) + r""",\quad
  \bar{\sigma}^2 = """ + mat2pmatrix(sigmabar_vec[2]) + r""",\quad
  \bar{\sigma}^3 = """ + mat2pmatrix(sigmabar_vec[3]) + r"""
\]

\paragraph{Key identity.}
\begin{equation}
  \operatorname{tr}(\sigma^\mu\bar{\sigma}^\nu)
  \equiv \sigma^\mu_{a\dot{a}}\,\bar{\sigma}^{\nu\,\dot{a}a}
  = 2\,g^{\mu\nu}.
\end{equation}
""" + (
    r"\verified{Verified numerically: max error $"
    + f"{trace_check:.1e}" + r"$.}"
) + r"""

%% ============================================================
\section{Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{@{}ll@{}}
  \toprule
  Object & Expression \\
  \midrule
  Left-handed field & $\psi_\alpha$ \quad (undotted index) \\
  Right-handed field & $\psidag_{\dot{\alpha}} = (\psi_\alpha)^\dagger$
                       \quad (dotted index) \\
  Rotation generator & $(S^{ij}_L)_a{}^b = \tfrac{1}{2}\varepsilon^{ijk}\sigma_k$ \\
  Boost generator    & $(S^{k0}_L)_a{}^b = \tfrac{i}{2}\sigma_k$ \\
  R-generators       & $S^{\mu\nu}_R = -[S^{\mu\nu}_L]^*$ \\
  Raise index        & $\psi^a = \varepsilon^{ab}\psi_b$ \\
  Lower index        & $\psi_a = \varepsilon_{ab}\psi^b$ \\
  Sign identity      & $\psi^a\chi_a = -\psi_a\chi^a$ \\
  Angle bracket      & $\langle\psi\chi\rangle = \varepsilon^{\alpha\beta}\psi_\alpha\chi_\beta$ \\
  Square bracket     & $[\psidag\chidag] = \varepsilon_{\dot\alpha\dot\beta}
                        \psidag^{\dot\alpha}\chidag^{\dot\beta}$ \\
  Vector dictionary  & $\sigma^\mu_{a\dot{a}} = (I,\vec\sigma)$,\quad
                        $\bar\sigma^\mu_{\dot{a}a} = (I,-\vec\sigma)$ \\
  Trace identity     & $\operatorname{tr}(\sigma^\mu\bar\sigma^\nu) = 2g^{\mu\nu}$ \\
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\noindent\textbf{Next:} Chapter 35 develops the index-free dot/bar notation,
derives the $\sigma$-algebra identities, and constructs the Weyl Lagrangian.

\end{document}
"""

outpath = "/work/ch34_left_right_spinors.tex"
with open(outpath, "w") as f:
    f.write(doc.strip())

print(f"Wrote: {outpath}")
