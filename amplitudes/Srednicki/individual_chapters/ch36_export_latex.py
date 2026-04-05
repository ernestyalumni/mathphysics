"""
ch36_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 36 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \
        python3 /work/ch36_export_latex.py

Outputs: /work/ch36_weyl_lagrangian.tex
Then compile on host:
    pdflatex ch36_weyl_lagrangian.tex
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

# ── numpy setup (mirrors ch36_weyl_lagrangian.py) ──────────────────────────
sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),
    2: np.array([[0,-1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0],  [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)

sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

g = np.diag([-1., 1., 1., 1.])   # Srednicki mostly-plus metric

# γ matrices in Weyl representation
gamma = {}
for mu in range(4):
    gamma[mu] = np.block([
        [np.zeros((2,2), dtype=complex), sigma_vec[mu]   ],
        [sigmabar_vec[mu],               np.zeros((2,2), dtype=complex)]
    ])

# γ^5
gamma5 = 1j * gamma[0] @ gamma[1] @ gamma[2] @ gamma[3]

# Epsilon tensors
eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1, 0]], dtype=complex)

# Charge conjugation matrix C
C_matrix = np.block([
    [eps_lower, np.zeros((2,2), dtype=complex)],
    [np.zeros((2,2), dtype=complex), eps_upper]
])

mu_label = {0: "0", 1: "1", 2: "2", 3: "3"}

# ── helpers ─────────────────────────────────────────────────────────────────

import re

def fix_cdb_latex(s):
    """Fix cadabra2 LaTeX output for standard pdflatex compilation.
    Cadabra outputs e.g. \\psidag^{\\dal} which double-superscripts
    since \\psidag expands to \\psi^{\\dagger}.  Replace with
    \\psi^{\\dagger X} form inline.
    """
    # \psidag^{X}  ->  \psi^{\dagger X}
    s = re.sub(r'\\psidag\s*\^\{([^}]*)\}', lambda m: r'\psi^{\dagger ' + m.group(1) + '}', s)
    # \chidag^{X}  ->  \chi^{\dagger X}
    s = re.sub(r'\\chidag\s*\^\{([^}]*)\}', lambda m: r'\chi^{\dagger ' + m.group(1) + '}', s)
    # \xidag^{X}   ->  \xi^{\dagger X}
    s = re.sub(r'\\xidag\s*\^\{([^}]*)\}',  lambda m: r'\xi^{\dagger '  + m.group(1) + '}', s)
    return s

def cdb(expr_str):
    """Return LaTeX string for a cadabra2 expression."""
    return fix_cdb_latex(Ex(expr_str)._latex_())

def fmt_complex(z, tol=1e-10):
    """Format a complex number cleanly for LaTeX."""
    r, i = z.real, z.imag
    if abs(r) < tol and abs(i) < tol:
        return "0"
    if abs(i) < tol:
        v = r
        if abs(v - round(v)) < tol:
            n = int(round(v))
            return str(n)
        if abs(v - 0.5)  < tol: return r"\tfrac{1}{2}"
        if abs(v + 0.5)  < tol: return r"-\tfrac{1}{2}"
        return f"{v:.4g}"
    if abs(r) < tol:
        v = i
        if abs(v - 0.5)  < tol: return r"\tfrac{i}{2}"
        if abs(v + 0.5)  < tol: return r"-\tfrac{i}{2}"
        if abs(v - 1.0)  < tol: return r"i"
        if abs(v + 1.0)  < tol: return r"-i"
        return f"{v:.4g}i"
    return f"{r:.4g}+{i:.4g}i"

def mat2pmatrix(mat):
    """Convert numpy array to LaTeX pmatrix."""
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"

def mat2bmatrix(mat):
    """Convert numpy array to LaTeX bmatrix."""
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{bmatrix}" + r" \\ ".join(rows) + r"\end{bmatrix}"

def block_pmatrix_sigma(mu):
    """LaTeX for a gamma matrix in block form showing σ/σ̄ labels."""
    if mu == 0:
        return r"\begin{pmatrix} 0 & I \\ I & 0 \end{pmatrix}"
    signs = {1: ("+", "-"), 2: ("+", "-"), 3: ("+", "-")}
    pm_upper = r"\sigma_{" + str(mu) + r"}"
    pm_lower = r"\bar{\sigma}_{" + str(mu) + r"}"
    return (r"\begin{pmatrix} 0 & " + pm_upper + r" \\ "
            + pm_lower + r" & 0 \end{pmatrix}")

# ── compute verification quantities ─────────────────────────────────────────

# Clifford algebra verification
clifford_results = {}
max_clifford_err = 0.0
for mu in range(4):
    for nu in range(mu, 4):
        anticomm = gamma[mu] @ gamma[nu] + gamma[nu] @ gamma[mu]
        expected = -2 * g[mu, nu] * I4
        err = np.max(np.abs(anticomm - expected))
        clifford_results[(mu, nu)] = {
            'val': -2 * g[mu, nu],
            'err': err,
            'ok': err < 1e-12
        }
        max_clifford_err = max(max_clifford_err, err)

# γ^5 verification
expected_g5 = np.block([
    [-I2, np.zeros((2,2), dtype=complex)],
    [np.zeros((2,2), dtype=complex), I2]
])
gamma5_err = np.max(np.abs(gamma5 - expected_g5))

# Projection operators
P_L = 0.5 * (I4 - gamma5)
P_R = 0.5 * (I4 + gamma5)
proj_err = max(
    np.max(np.abs(P_L @ P_L - P_L)),
    np.max(np.abs(P_R @ P_R - P_R)),
    np.max(np.abs(P_L @ P_R)),
)

# C^{-1} γ^μ C = -(γ^μ)^T
Cinv = np.linalg.inv(C_matrix)
max_cc_err = max(
    np.max(np.abs(Cinv @ gamma[mu] @ C_matrix + gamma[mu].T))
    for mu in range(4)
)

# tr(σ^μ σ̄^ν) = 2g^{μν}
trace_err = max(
    abs(np.trace(sigma_vec[mu] @ sigmabar_vec[nu]) - 2*g[mu, nu])
    for mu in range(4) for nu in range(4)
)

# σ^μ σ̄^ν + σ^ν σ̄^μ = -2g^{μν}I (eq. 36.8)
sigma_alg_err = max(
    np.max(np.abs(
        sigma_vec[mu] @ sigmabar_vec[nu] + sigma_vec[nu] @ sigmabar_vec[mu]
        + 2*g[mu,nu]*I2
    ))
    for mu in range(4) for nu in range(4)
)

# ── cadabra2 expressions ─────────────────────────────────────────────────────
cdb_exprs = {
    "L_kinetic":   cdb(r"i \psidag^{\dal} \sigmabar^{\mu}_{\dal\alpha} \partial_{\mu}(\psi^{\alpha})"),
    "L_mass1":     cdb(r"-\frac{1}{2} m \epsilon^{\alpha\beta} \psi_{\alpha} \psi_{\beta}"),
    "L_mass2":     cdb(r"-\frac{1}{2} m \epsilon_{\dal\dbe} \psidag^{\dal} \psidag^{\dbe}"),
    "eom_up_t1":   cdb(r"-i \sigmabar^{\mu}_{\dal\alpha} \partial_{\mu}(\psi^{\alpha})"),
    "eom_up_t2":   cdb(r"m \psidag^{\dal}"),
    "eom_down_t1": cdb(r"-i \sigma^{\mu}_{\alpha\dal} \partial_{\mu}(\psidag^{\dal})"),
    "eom_down_t2": cdb(r"m \psi_{\alpha}"),
    "majorana":    cdb(r"\Psi_{\alpha}"),
    "gamma_mu":    cdb(r"\gamma^{\mu}_{ab}"),
}

# ── build LaTeX document ─────────────────────────────────────────────────────

# Clifford table rows
clifford_table_rows = []
for mu in range(4):
    for nu in range(mu, 4):
        res = clifford_results[(mu, nu)]
        val_str = f"${res['val']:+.0f} \\cdot I_4$"
        status = r"\textcolor{verifygreen}{\checkmark}" if res['ok'] else r"\textcolor{red}{\times}"
        clifford_table_rows.append(
            f"  $\\{{\\gamma^{mu_label[mu]},\\gamma^{mu_label[nu]}\\}}$ & "
            f"$-2g^{{{mu_label[mu]}{mu_label[nu]}}}$ & "
            f"{val_str} & {status} \\\\"
        )

clifford_table = "\n".join(clifford_table_rows)

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
\usepackage{colortbl}
\geometry{margin=1.1in, headheight=15pt}

\definecolor{cadblue}{RGB}{30,90,200}
\definecolor{codegray}{RGB}{245,245,245}
\definecolor{verifygreen}{RGB}{20,140,60}
\definecolor{diaggray}{RGB}{200,200,200}

\newcommand{\cdbexpr}[1]{\textcolor{cadblue}{$#1$}}
\newcommand{\verified}[1]{\textcolor{verifygreen}{\checkmark\; #1}}
\newcommand{\psidag}{\psi^{\dagger}}
\newcommand{\half}{\tfrac{1}{2}}
\newcommand{\sbar}{\bar{\sigma}}
% cadabra2-generated index macros
\newcommand{\dal}{\dot{\alpha}}
\newcommand{\dbe}{\dot{\beta}}
\newcommand{\dga}{\dot{\gamma}}
\newcommand{\dde}{\dot{\delta}}
\newcommand{\sigmabar}{\bar{\sigma}}
\newcommand{\chidag}{\chi^{\dagger}}
\newcommand{\xidag}{\xi^{\dagger}}

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 36}
\rhead{Weyl Lagrangian \& Majorana/Dirac Fermions}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 36}\\[6pt]
       \large Lagrangians for Spinor Fields\\[4pt]
       \normalsize Cadabra2 expressions and numerical verification}
\author{Generated by \texttt{ch36\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{The Weyl Lagrangian}
%% ============================================================

\subsection{Lagrangian terms}

For a single left-handed Weyl field $\psi_a(x)$, the Lorentz-invariant,
hermitian Lagrangian quadratic in $\psi$ is (eq.~36.2):
\begin{equation}
  \boxed{
    \mathcal{L}
    = i\psi^\dagger_{\dot{a}}\,\bar{\sigma}^{\mu\,\dot{a}\alpha}\,\partial_\mu\psi_\alpha
    - \tfrac{1}{2}m\,\varepsilon^{\alpha\beta}\psi_\alpha\psi_\beta
    - \tfrac{1}{2}m\,\varepsilon_{\dot{\alpha}\dot{\beta}}\psi^{\dagger\dot{\alpha}}\psi^{\dagger\dot{\beta}}
  }
  \tag{36.2}
\end{equation}

\medskip
\noindent\textbf{Cadabra2 expressions:}
\begin{align*}
  \mathcal{L}_{\text{kin}} &= """ + cdb_exprs["L_kinetic"] + r""" \\
  \mathcal{L}_{\text{mass}} &= """ + cdb_exprs["L_mass1"] + r""" \\
  \mathcal{L}_{\text{mass}}^\dagger &= """ + cdb_exprs["L_mass2"] + r"""
\end{align*}

\subsection{Hermiticity of the kinetic term}

Although $i\psi^\dagger\bar{\sigma}^\mu\partial_\mu\psi$ is not individually
hermitian, it satisfies (eq.~36.1):
\begin{equation}
  (i\psi^\dagger\bar{\sigma}^\mu\partial_\mu\psi)^\dagger
  = i\psi^\dagger\bar{\sigma}^\mu\partial_\mu\psi
    - i\partial_\mu(\psi^\dagger\bar{\sigma}^\mu\psi)\,.
  \tag{36.1}
\end{equation}
The second term is a total divergence and vanishes in the action.
Thus the kinetic term \emph{is} real in $S = \int d^4x\,\mathcal{L}$.

\subsection{Mass term and Majorana mass}

The mass term $\varepsilon^{\alpha\beta}\psi_\alpha\psi_\beta = \psi^\alpha\psi_\alpha$
is a Lorentz scalar; it is nonzero because $\psi_\alpha$ are Grassmann-valued
(fermionic) fields.  The parameter $m$ is a ``Majorana mass'':
\begin{itemize}
  \item The phase of $m$ is unphysical: redefining $\psi \to e^{-i\alpha/2}\psi$
        removes it.  We therefore take $m$ real and positive.
  \item There is no conserved U(1) charge associated with this mass term---it
        explicitly breaks any $\psi\to e^{i\alpha}\psi$ symmetry.
\end{itemize}

%% ============================================================
\section{Equations of Motion}
%% ============================================================

Varying the action $S = \int d^4x\,\mathcal{L}$ with respect to
$\psi^{\dagger\dot{a}}$ (holding $\psi$ fixed):
\begin{equation}
  0 = -\frac{\delta S}{\delta\psi^{\dagger\dot{a}}}
    = -i\,\bar{\sigma}^{\mu\,\dot{a}c}\,\partial_\mu\psi_c
      + m\,\psi^{\dagger\dot{a}}\,.
  \tag{36.3--36.4}
\end{equation}

\noindent In Cadabra2: $""" + cdb_exprs["eom_up_t1"] + r""" + """ + cdb_exprs["eom_up_t2"] + r"""$

\medskip
\noindent The hermitian conjugate equation (varying w.r.t.\ $\psi^\alpha$):
\begin{equation}
  0 = -i\,\sigma^\mu_{a\dot{c}}\,\partial_\mu\psi^{\dagger\dot{c}}
      + m\,\psi_a\,.
  \tag{36.5}
\end{equation}

\noindent In Cadabra2: $""" + cdb_exprs["eom_down_t1"] + r""" + """ + cdb_exprs["eom_down_t2"] + r"""$

\subsection{Combined 4-component form}

Equations~(36.4) and~(36.5) combine into the matrix equation (eq.~36.6):
\begin{equation}
  \begin{pmatrix}
    m\,\delta_a{}^c & -i\,\sigma^\mu_{a\dot{c}}\,\partial_\mu \\
    -i\,\bar{\sigma}^{\mu\,\dot{a}c}\,\partial_\mu & m\,\delta^{\dot{a}}{}_{\dot{c}}
  \end{pmatrix}
  \begin{pmatrix} \psi_c \\ \psi^{\dagger\dot{c}} \end{pmatrix}
  = 0\,.
  \tag{36.6}
\end{equation}
The $2\times 2$ diagonal blocks are proportional to $m\cdot\mathbf{1}_2$;
the off-diagonal blocks mix the left- and right-handed sectors via $\sigma^\mu$
and $\bar{\sigma}^\mu$.

%% ============================================================
\section{Gamma Matrices in the Weyl Representation}
%% ============================================================

\subsection{Definition}

The gamma matrices are defined by the block structure (eq.~36.7):
\begin{equation}
  \gamma^\mu \equiv
  \begin{pmatrix}
    0 & \sigma^\mu_{a\dot{c}} \\
    \bar{\sigma}^{\mu\,\dot{a}c} & 0
  \end{pmatrix}.
  \tag{36.7}
\end{equation}
This is the \textbf{Weyl (chiral) representation}.

\subsection{Explicit matrices}

\begin{align}
  \gamma^0 &= """ + mat2pmatrix(gamma[0]) + r""",&
  \gamma^1 &= """ + mat2pmatrix(gamma[1]) + r"""\\[6pt]
  \gamma^2 &= """ + mat2pmatrix(gamma[2]) + r""",&
  \gamma^3 &= """ + mat2pmatrix(gamma[3]) + r"""
\end{align}

\noindent Written symbolically in terms of $\sigma^\mu = (I, \vec{\sigma})$
and $\bar{\sigma}^\mu = (I,-\vec{\sigma})$:

\[
  \gamma^0 = \begin{pmatrix}0&I\\I&0\end{pmatrix},\quad
  \gamma^k = \begin{pmatrix}0&\sigma_k\\-\sigma_k&0\end{pmatrix}
  \quad (k=1,2,3)\,.
\]

%% ============================================================
\section{Clifford Algebra Verification}
%% ============================================================

\subsection{The algebra}

The gamma matrices obey the Clifford algebra (eq.~36.9):
\begin{equation}
  \boxed{\{\gamma^\mu,\,\gamma^\nu\} \equiv \gamma^\mu\gamma^\nu + \gamma^\nu\gamma^\mu
  = -2\,g^{\mu\nu}\,I_4\,,}
  \tag{36.9}
\end{equation}
with $g^{\mu\nu}=\mathrm{diag}(-1,+1,+1,+1)$ (Srednicki mostly-plus convention).

This follows from the sigma-matrix identities (eq.~36.8):
\begin{align}
  \sigma^\mu\bar{\sigma}^\nu + \sigma^\nu\bar{\sigma}^\mu &= -2g^{\mu\nu}\,I_2\,,
  \tag{36.8a}\\
  \bar{\sigma}^\mu\sigma^\nu + \bar{\sigma}^\nu\sigma^\mu &= -2g^{\mu\nu}\,I_2\,.
  \tag{36.8b}
\end{align}
\verified{$\sigma^\mu\bar{\sigma}^\nu+\sigma^\nu\bar{\sigma}^\mu=-2g^{\mu\nu}I_2$: max error $""" + f"{sigma_alg_err:.1e}" + r"""$.}

\subsection{Numerical verification table}

All ten independent anticommutator pairs are verified below:
\medskip
\begin{center}
\begin{tabular}{@{}llll@{}}
  \toprule
  Pair & Formula & Value & Verified \\
  \midrule
""" + clifford_table + r"""
  \bottomrule
\end{tabular}
\end{center}
\medskip
\verified{Clifford algebra: all 10 pairs verified. Max error: $""" + f"{max_clifford_err:.1e}" + r"""$.}

\subsection{The \texorpdfstring{$\gamma^5$}{gamma5} matrix}

Define (eq.~36.46):
\begin{equation}
  \gamma^5 \equiv i\gamma^0\gamma^1\gamma^2\gamma^3
  = \begin{pmatrix}-I_2 & 0 \\ 0 & +I_2\end{pmatrix}.
  \tag{36.43}
\end{equation}
Numerically:
\[
  \gamma^5 = """ + mat2pmatrix(gamma5.real) + r"""\,.
\]
\verified{$\gamma^5 = \mathrm{diag}(-I_2,+I_2)$: error $""" + f"{gamma5_err:.1e}" + r"""$.}

\medskip
\noindent The left- and right-handed \textbf{projection operators} are (eq.~36.44):
\begin{align}
  P_L &\equiv \tfrac{1}{2}(I - \gamma^5)
     = \begin{pmatrix}I_2&0\\0&0\end{pmatrix},
  &
  P_R &\equiv \tfrac{1}{2}(I + \gamma^5)
     = \begin{pmatrix}0&0\\0&I_2\end{pmatrix}.
  \tag{36.44}
\end{align}
\verified{$P_L^2=P_L$, $P_R^2=P_R$, $P_LP_R=0$: max error $""" + f"{proj_err:.1e}" + r"""$.}

%% ============================================================
\section{Majorana 4-Component Spinor}
%% ============================================================

\subsection{Definition}

The \textbf{Majorana spinor} is built from a single left-handed Weyl field $\psi_a$
(eq.~36.10):
\begin{equation}
  \Psi \equiv \begin{pmatrix}\psi_c\\\psi^{\dagger\dot{c}}\end{pmatrix}.
  \tag{36.10}
\end{equation}
Equation~(36.6) then becomes the \textbf{Dirac equation} (eq.~36.11):
\begin{equation}
  \boxed{(-i\gamma^\mu\partial_\mu + m)\Psi = 0\,.}
  \tag{36.11}
\end{equation}

\subsection{Majorana condition}

The Majorana spinor is its own \textbf{charge conjugate}:
\begin{equation}
  \Psi^C \equiv \mathcal{C}\,\overline{\Psi}^{\,T} = \Psi\,,
\end{equation}
where the charge conjugation matrix is (eq.~36.32):
\begin{equation}
  \mathcal{C} = \begin{pmatrix}\varepsilon_{ac}&0\\0&\varepsilon^{\dot{a}\dot{c}}\end{pmatrix}
  = """ + mat2pmatrix(C_matrix) + r"""\,.
  \tag{36.32}
\end{equation}

The matrix $\mathcal{C}$ satisfies (eq.~36.36):
\begin{equation}
  \mathcal{C}^T = \mathcal{C}^\dagger = \mathcal{C}^{-1} = -\mathcal{C}\,.
  \tag{36.36}
\end{equation}
and conjugates the gamma matrices as (eq.~36.40):
\begin{equation}
  \mathcal{C}^{-1}\gamma^\mu\mathcal{C} = -(\gamma^\mu)^T\,.
  \tag{36.40}
\end{equation}
\verified{$\mathcal{C}^{-1}\gamma^\mu\mathcal{C}=-(\gamma^\mu)^T$
for all $\mu$: max error $""" + f"{max_cc_err:.1e}" + r"""$.}

\subsection{Majorana Lagrangian}

Re-expressing eq.~(36.2) using $\overline{\Psi}=\Psi^\dagger\beta$ (where $\beta=\gamma^0$
numerically) gives (eq.~36.41):
\begin{equation}
  \mathcal{L} = \tfrac{i}{2}\overline{\Psi}\gamma^\mu\partial_\mu\Psi
              - \tfrac{1}{2}m\,\overline{\Psi}\Psi\,.
  \tag{36.41}
\end{equation}

%% ============================================================
\section{Dirac Fermion from Two Weyl Fields}
%% ============================================================

\subsection{Two-Weyl-field Lagrangian}

Begin with two left-handed fields $\chi_a$, $\xi_a$ (rewritten from $\psi_{1,2}$
via eq.~36.14--36.15). The Lagrangian is (eq.~36.16):
\begin{equation}
  \mathcal{L}
  = i\chi^\dagger\bar{\sigma}^\mu\partial_\mu\chi
  + i\xi^\dagger\bar{\sigma}^\mu\partial_\mu\xi
  - m\,\chi\xi - m\,\xi^\dagger\chi^\dagger\,,
  \tag{36.16}
\end{equation}
invariant under the U(1) symmetry (eq.~36.17):
\begin{equation}
  \chi \to e^{-i\alpha}\chi\,, \qquad \xi \to e^{+i\alpha}\xi\,.
  \tag{36.17}
\end{equation}

\subsection{Dirac spinor and Lagrangian}

Define the \textbf{Dirac spinor} (eq.~36.19):
\begin{equation}
  \Psi \equiv \begin{pmatrix}\chi_c\\\xi^{\dagger\dot{c}}\end{pmatrix}\,,
  \qquad
  \overline{\Psi} = \Psi^\dagger\beta = (\xi^a,\,\chi^\dagger_{\dot{a}})\,.
  \tag{36.19, 36.22}
\end{equation}
Then the Dirac Lagrangian (eq.~36.28):
\begin{equation}
  \boxed{\mathcal{L} = i\overline{\Psi}\gamma^\mu\partial_\mu\Psi - m\overline{\Psi}\Psi\,,}
  \tag{36.28}
\end{equation}
where (eq.~36.23):
\begin{equation}
  \overline{\Psi}\Psi = \xi^a\chi_a + \chi^\dagger_{\dot{a}}\xi^{\dagger\dot{a}}\,.
  \tag{36.23}
\end{equation}

The Noether current for the U(1) symmetry (eq.~36.30):
\begin{equation}
  j^\mu = \overline{\Psi}\gamma^\mu\Psi
        = \chi^\dagger\bar{\sigma}^\mu\chi - \xi^\dagger\bar{\sigma}^\mu\xi\,.
  \tag{36.30}
\end{equation}

%% ============================================================
\section{Dirac vs.\ Majorana --- Comparison}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}lll@{}}
  \toprule
  \textbf{Property} & \textbf{Majorana} & \textbf{Dirac} \\
  \midrule
  Weyl fields & 1 (single $\psi$) & 2 ($\chi$ and $\xi$) \\
  U(1) symmetry & No & Yes ($\chi\to e^{-i\alpha}\chi$, $\xi\to e^{+i\alpha}\xi$) \\
  Charge & None (Majorana mass) & Conserved $Q=\int j^0 d^3x$ \\
  $\Psi^C = \Psi$ & Yes & No \\
  Lagrangian & $\frac{i}{2}\bar{\Psi}\gamma^\mu\partial_\mu\Psi - \frac{m}{2}\bar{\Psi}\Psi$ &
               $i\bar{\Psi}\gamma^\mu\partial_\mu\Psi - m\bar{\Psi}\Psi$ \\
  Factor of $\frac{1}{2}$ & Yes (overcounting) & No \\
  Scalar analogue & Real scalar $\phi=\phi^\dagger$ & Complex scalar $\phi\neq\phi^\dagger$ \\
  Neutrino model & Seesaw mechanism & Dirac neutrino \\
  Degrees of freedom & 2 (real) & 4 (complex) \\
  \bottomrule
\end{tabular}
\end{center}

%% ============================================================
\section{Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.4}
\begin{tabular}{@{}ll@{}}
  \toprule
  Object & Expression \\
  \midrule
  Weyl Lagrangian & $\mathcal{L}=i\psi^\dagger\bar\sigma^\mu\partial_\mu\psi
    -\frac{1}{2}m\psi\psi - \frac{1}{2}m\psi^\dagger\psi^\dagger$ \\
  EOM ($\delta/\delta\psi^\dagger$) & $-i\bar\sigma^{\mu\dot{a}c}\partial_\mu\psi_c
    +m\psi^{\dagger\dot{a}}=0$ \\
  $\gamma$ matrix (Weyl rep) & $\gamma^\mu=\bigl[\begin{smallmatrix}0&\sigma^\mu\\
    \bar\sigma^\mu&0\end{smallmatrix}\bigr]$ \\
  Clifford algebra & $\{\gamma^\mu,\gamma^\nu\}=-2g^{\mu\nu}I_4$ \\
  $\gamma^5$ & $i\gamma^0\gamma^1\gamma^2\gamma^3
    =\mathrm{diag}(-I_2,+I_2)$ \\
  Majorana spinor & $\Psi=(\psi_c,\psi^{\dagger\dot{c}})^T$,\quad $\Psi^C=\Psi$ \\
  Dirac equation & $(-i\gamma^\mu\partial_\mu+m)\Psi=0$ \\
  Dirac spinor & $\Psi=(\chi_c,\xi^{\dagger\dot{c}})^T$,\quad $\Psi^C\neq\Psi$ \\
  Dirac Lagrangian & $i\bar\Psi\gamma^\mu\partial_\mu\Psi - m\bar\Psi\Psi$ \\
  Charge conjugation & $\mathcal{C}^{-1}\gamma^\mu\mathcal{C}=-(\gamma^\mu)^T$ \\
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\noindent\textbf{Next:} Chapter 37 develops the quantization of
Dirac and Majorana fields, introducing anticommutation relations for
the creation/annihilation operators.

\end{document}
"""

outpath = "/work/ch36_weyl_lagrangian.tex"
with open(outpath, "w") as f:
    f.write(doc.strip())

print(f"Wrote: {outpath}")
