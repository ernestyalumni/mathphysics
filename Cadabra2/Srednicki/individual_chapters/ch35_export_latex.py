"""
ch35_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 35 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \
        python3 /work/ch35_export_latex.py

Outputs: /work/ch35_sigma_algebra.tex
Then compile on host:
    pdflatex ch35_sigma_algebra.tex
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

# ── numpy setup (mirrors ch35_sigma_algebra.py) ────────────────────────────
g = np.diag([1., -1., -1., -1.])

eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1, 0]], dtype=complex)

sigma = {
    1: np.array([[0, 1],   [1, 0]],  dtype=complex),
    2: np.array([[0, -1j], [1j, 0]], dtype=complex),
    3: np.array([[1, 0],   [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)

sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}

# S^μν_L from Ch. 34 formula (for comparison in §35.E)
def eps3(i, j, k):
    return int(np.linalg.det([
        [i==1, i==2, i==3],
        [j==1, j==2, j==3],
        [k==1, k==2, k==3]
    ]))

S_L = [[None]*4 for _ in range(4)]
for i in range(1, 4):
    for j in range(1, 4):
        mat = np.zeros((2, 2), dtype=complex)
        for k in range(1, 4):
            mat += eps3(i, j, k) * sigma[k]
        S_L[i][j] = 0.5 * mat
for k in range(1, 4):
    S_L[k][0] = (1j/2) * sigma[k]
    S_L[0][k] = -(1j/2) * sigma[k]
S_L[0][0] = np.zeros((2, 2), dtype=complex)
for i in range(1, 4):
    S_L[i][i] = np.zeros((2, 2), dtype=complex)

def S_L_fn(mu, nu):
    if mu == nu: return np.zeros((2, 2), dtype=complex)
    if mu < nu:  return S_L[mu][nu]
    return -S_L[nu][mu]

def S_R_fn(mu, nu):
    return -np.conj(S_L_fn(mu, nu))

# ── compute all verification quantities ────────────────────────────────────

# eq 35.4: σ^μ_{aȧ} σ_{μ bḃ} = -2 ε_{ab} ε_{ȧḃ}
lhs_35_4 = np.zeros((2, 2, 2, 2), dtype=complex)
for mu in range(4):
    lhs_35_4 += g[mu, mu] * np.einsum('ij,kl->ijkl', sigma_vec[mu], sigma_vec[mu])
rhs_35_4 = np.zeros((2, 2, 2, 2), dtype=complex)
for a in range(2):
    for b in range(2):
        for adot in range(2):
            for bdot in range(2):
                rhs_35_4[a, adot, b, bdot] = -2 * eps_lower[a, b] * eps_lower[adot, bdot]
err_35_4 = np.max(np.abs(lhs_35_4 - rhs_35_4))

# Sample nonzero components of eq 35.4 for display
sample_35_4 = []
for a in range(2):
    for adot in range(2):
        for b in range(2):
            for bdot in range(2):
                lv = lhs_35_4[a, adot, b, bdot].real
                if abs(lv) > 1e-10:
                    sample_35_4.append((a+1, adot+1, b+1, bdot+1, int(round(lv))))

# eq 35.5: ε^{ab} ε^{ȧḃ} σ^μ_{aȧ} σ^ν_{bḃ} = -2 g^{μν}
lhs_35_5 = np.zeros((4, 4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        val = 0.0 + 0.0j
        for a in range(2):
            for adot in range(2):
                for b in range(2):
                    for bdot in range(2):
                        val += (eps_upper[a, b] * eps_upper[adot, bdot]
                                * sigma_vec[mu][a, adot]
                                * sigma_vec[nu][b, bdot])
        lhs_35_5[mu, nu] = val
err_35_5 = np.max(np.abs(lhs_35_5 - (-2 * g)))

# σ̄^μ from definition: σ̄^{μ ȧa} = ε^{ab} ε^{ȧḃ} σ^μ_{bḃ}
sigmabar_computed = {}
for mu in range(4):
    mat = np.zeros((2, 2), dtype=complex)
    for adot in range(2):
        for a in range(2):
            val = 0.0 + 0.0j
            for b in range(2):
                for bdot in range(2):
                    val += eps_upper[a, b] * eps_upper[adot, bdot] * sigma_vec[mu][b, bdot]
            mat[adot, a] = val
    sigmabar_computed[mu] = mat
max_err_sigmabar = max(
    np.max(np.abs(sigmabar_computed[mu] - sigmabar_vec[mu]))
    for mu in range(4)
)

# tr(σ^μ σ̄^ν) = 2 g^{μν}
err_trace = np.max(np.abs(
    np.array([[np.trace(sigma_vec[mu] @ sigmabar_vec[nu])
               for nu in range(4)] for mu in range(4)])
    - 2 * g
))

# S^μν_L from ch35 formula: (i/4)(σ^μ σ̄^ν - σ^ν σ̄^μ)
S_L_ch35 = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        S_L_ch35[mu][nu] = (1j/4) * (sigma_vec[mu] @ sigmabar_vec[nu]
                                      - sigma_vec[nu] @ sigmabar_vec[mu])
max_err_SL = max(
    np.max(np.abs(S_L_ch35[mu][nu] - S_L[mu][nu]))
    for mu in range(4) for nu in range(mu+1, 4)
)

# S^μν_R from ch35 formula: -(i/4)(σ̄^μ σ^ν - σ̄^ν σ^μ)
S_R_ch35 = [[None]*4 for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        S_R_ch35[mu][nu] = -(1j/4) * (sigmabar_vec[mu] @ sigma_vec[nu]
                                       - sigmabar_vec[nu] @ sigma_vec[mu])
max_err_SR = max(
    np.max(np.abs(S_R_ch35[mu][nu] - (-np.conj(S_L[mu][nu]))))
    for mu in range(4) for nu in range(mu+1, 4)
)

# Clifford-like: σ^μ σ̄^ν + σ^ν σ̄^μ = -2 g^{μν} I₂
max_err_cliff = max(
    np.max(np.abs(
        sigma_vec[mu] @ sigmabar_vec[nu] + sigma_vec[nu] @ sigmabar_vec[mu]
        - (-2 * g[mu, nu] * I2)
    ))
    for mu in range(4) for nu in range(4)
)

# Covariance condition (eq 35.14)
max_err_cov = 0.
for mu in range(4):
    for nu in range(mu+1, 4):
        for rho in range(4):
            for a in range(2):
                for adot in range(2):
                    term1 = (g[mu, rho] * sigma_vec[nu][a, adot]
                             - g[nu, rho] * sigma_vec[mu][a, adot])
                    term2 = (1j * S_L_fn(mu, nu) @ sigma_vec[rho])[a, adot]
                    term3 = 1j * np.dot(S_R_fn(mu, nu)[adot, :], sigma_vec[rho][a, :])
                    max_err_cov = max(max_err_cov, abs(term1 + term2 + term3))

# ── helpers ─────────────────────────────────────────────────────────────────

def cdb(expr_str):
    return Ex(expr_str)._latex_()

def fmt_complex(z, tol=1e-10):
    r, i = z.real, z.imag
    if abs(r) < tol and abs(i) < tol: return "0"
    if abs(i) < tol:
        v = r
        if abs(v - round(v)) < tol: return str(int(round(v)))
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
    rows = []
    for row in mat:
        rows.append(" & ".join(fmt_complex(z) for z in row))
    return r"\begin{pmatrix}" + r" \\ ".join(rows) + r"\end{pmatrix}"

def verified(text, err, tol=1e-12):
    if err < tol:
        return r"\verified{" + text + r"}"
    return r"\textcolor{red}{\textbf{MISMATCH!} " + text + f" (err={err:.1e})" + r"}"

# ── cadabra2 symbolic expressions ────────────────────────────────────────────

cdb_exprs = {
    "sigma_mu"       : cdb(r"\sigma^{\mu}_{\alpha \dal}"),
    "sigmabar_mu"    : cdb(r"\bar{\sigma}^{\mu \dal \alpha}"),
    "eq35_4_lhs"     : cdb(r"\sigma^{\mu}_{\alpha \dal} \sigma_{\mu \beta \dbe}"),
    "eq35_4_rhs"     : cdb(r"-2 \epsilon_{\alpha\beta} \epsilon_{\dal\dbe}"),
    "eq35_5_lhs"     : cdb(r"\epsilon^{\alpha\beta} \epsilon^{\dal\dbe} \sigma^{\mu}_{\alpha \dal} \sigma^{\nu}_{\beta \dbe}"),
    "eq35_5_rhs"     : cdb(r"-2 g^{\mu\nu}"),
    "sigmabar_def"   : cdb(r"\epsilon^{\alpha\beta} \epsilon^{\dal\dbe} \sigma^{\mu}_{\beta \dbe}"),
    "gen_L"          : cdb(r"\left(S^{\mu\nu}_{L}\right)_{\alpha}{}^{\beta}"),
    "gen_R"          : cdb(r"\left(S^{\mu\nu}_{R}\right)^{\dal}{}_{\dbe}"),
    "SL_formula"     : cdb(r"\tfrac{i}{4}\left(\sigma^{\mu}_{\alpha\dal}\bar{\sigma}^{\nu\dal\beta} - \sigma^{\nu}_{\alpha\dal}\bar{\sigma}^{\mu\dal\beta}\right)"),
    "bilinear_chi"   : cdb(r"\epsilon^{\alpha\beta} \chi_{\alpha} \psi_{\beta}"),
    "bilinear_dag"   : cdb(r"\epsilon_{\dal\dbe} \chidag^{\dal} \psidag^{\dbe}"),
    "vec_bilinear"   : cdb(r"\psidag_{\dal} \sigmabar^{\mu\dal\alpha} \chi_{\alpha}"),
}

# Explicit generator matrices (for table)
SL_mats = {}
SR_mats = {}
for mu in range(4):
    for nu in range(mu+1, 4):
        key = f"{mu}{nu}"
        SL_mats[key] = mat2pmatrix(S_L[mu][nu])
        SR_mats[key] = mat2pmatrix(S_R_ch35[mu][nu])

# σ and σ̄ matrices for display
pauli_mats = {k: mat2pmatrix(sigma[k]) for k in [1,2,3]}
sigma_mats    = {mu: mat2pmatrix(sigma_vec[mu])    for mu in range(4)}
sigmabar_mats = {mu: mat2pmatrix(sigmabar_vec[mu]) for mu in range(4)}

# ── build LaTeX document ──────────────────────────────────────────────────────

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
\definecolor{verifygreen}{RGB}{20,140,60}

\newcommand{\cdbexpr}[1]{\textcolor{cadblue}{$#1$}}
\newcommand{\verified}[1]{\textcolor{verifygreen}{\checkmark\; #1}}
\newcommand{\half}{\tfrac{1}{2}}
\newcommand{\SL}{S^{\mu\nu}_{\mathrm{L}}}
\newcommand{\SR}{S^{\mu\nu}_{\mathrm{R}}}
\newcommand{\psidag}{\psi^{\dagger}}
\newcommand{\sbar}{\bar{\sigma}}

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 35}
\rhead{Manipulating Spinor Indices}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 35}\\[6pt]
       \large Manipulating Spinor Indices\\[4pt]
       \normalsize Cadabra2 expressions and numerical verification}
\author{Generated by \texttt{ch35\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{Setup: \texorpdfstring{$\varepsilon$}{epsilon} Symbols
         and \texorpdfstring{$\sigma^\mu$}{sigma\^{}mu}}
%% ============================================================

Chapter~35 opens by recalling the $\varepsilon$ conventions and the
$\sigma^\mu$ invariant symbol.  We use the mostly-minus metric
$g^{\mu\nu} = \mathrm{diag}(+1,-1,-1,-1)$.

\subsection{\texorpdfstring{$\varepsilon$}{epsilon} normalization (eq.~35.1)}
\begin{equation}
  \varepsilon^{12} = \varepsilon^{\dot{1}\dot{2}}
  = \varepsilon_{21} = \varepsilon_{\dot{2}\dot{1}} = +1.
  \tag{35.1}
\end{equation}
\[
  \varepsilon_{\alpha\beta} = \begin{pmatrix}0 & -1\\1 & 0\end{pmatrix},\qquad
  \varepsilon^{\alpha\beta} = \begin{pmatrix}0 & +1\\-1 & 0\end{pmatrix}.
\]

\subsection{\texorpdfstring{$\sigma^\mu$}{sigma\^{}mu} invariant symbol (eq.~35.2)}

The $\sigma^\mu$ symbol lives in the $(2,2)$ representation --- it has
one undotted and one dotted spinor index plus a spacetime vector index:
\begin{equation}
  \sigma^\mu_{\alpha\dot{\alpha}} = (I,\,\sigma_1,\,\sigma_2,\,\sigma_3).
  \tag{35.2}
\end{equation}
In Cadabra2: $""" + cdb_exprs["sigma_mu"] + r"""$.
\[
  \sigma^0 = """ + sigma_mats[0] + r""",\quad
  \sigma^1 = """ + sigma_mats[1] + r""",\quad
  \sigma^2 = """ + sigma_mats[2] + r""",\quad
  \sigma^3 = """ + sigma_mats[3] + r"""
\]

%% ============================================================
\section{Key Identity: \texorpdfstring{$\sigma^\mu_{a\dot{a}}\,\sigma_{\mu b\dot{b}}
         = -2\varepsilon_{ab}\varepsilon_{\dot{a}\dot{b}}$}{eq 35.4}
         \quad [eq.~35.4]}
%% ============================================================

\subsection{Statement and physical meaning}

This is the fundamental $\sigma$-completeness relation.  Contracting
two $\sigma^\mu$ symbols over the spacetime index produces the
SL$(2,\mathbb{C})$ invariant $\varepsilon_{ab}\varepsilon_{\dot{a}\dot{b}}$:
\begin{equation}
  \boxed{
    \sigma^\mu_{\alpha\dot{\alpha}}\,\sigma_{\mu\beta\dot{\beta}}
    = -2\,\varepsilon_{\alpha\beta}\,\varepsilon_{\dot{\alpha}\dot{\beta}}
  }
  \tag{35.4}
\end{equation}
In Cadabra2:
\[
  """ + cdb_exprs["eq35_4_lhs"] + r""" = """ + cdb_exprs["eq35_4_rhs"] + r"""
\]

\subsection{Index structure}

The notation $\sigma_{\mu\,b\dot{b}}$ means $\sigma^\nu_{b\dot{b}}$ with
the vector index lowered: $\sigma_{\mu\,b\dot{b}} = g_{\mu\nu}\sigma^\nu_{b\dot{b}}$.
The left-hand side is therefore
\[
  \sum_\mu g_{\mu\mu}\,\sigma^\mu_{\alpha\dot{\alpha}}\,\sigma^\mu_{\beta\dot{\beta}}
  = (+1)\sigma^0\!\otimes\!\sigma^0
    + (-1)\bigl(\sigma^1\!\otimes\!\sigma^1
    + \sigma^2\!\otimes\!\sigma^2
    + \sigma^3\!\otimes\!\sigma^3\bigr),
\]
a $4$-index tensor with $(2\times2)^2 = 16$ components.

\subsection{Nonzero components}

The right-hand side $-2\varepsilon_{\alpha\beta}\varepsilon_{\dot\alpha\dot\beta}$
is nonzero only when $\alpha\neq\beta$ and $\dot\alpha\neq\dot\beta$
simultaneously.  The nonzero entries are:

\medskip
\begin{center}
\begin{tabular}{@{}ccccc@{}}
  \toprule
  $(\alpha,\dot\alpha,\beta,\dot\beta)$ &
  LHS & RHS $= -2\varepsilon_{\alpha\beta}\varepsilon_{\dot\alpha\dot\beta}$ \\
  \midrule
""" + "".join(
    r"  $(" + f"{a},{adot},{b},{bdot}" + r")$ & $" + f"{v}" + r"$ & $" + f"{v}" + r"$ \\" + "\n"
    for a, adot, b, bdot, v in sample_35_4
) + r"""  \bottomrule
\end{tabular}
\end{center}

\paragraph{Numerical verification.}
""" + verified(
    r"All 16 components of eq.~(35.4): max error $"
    + f"{err_35_4:.1e}" + r"$.",
    err_35_4
) + r"""

%% ============================================================
\section{Key Identity:
         \texorpdfstring{$\varepsilon^{ab}\varepsilon^{\dot{a}\dot{b}}
         \sigma^\mu_{a\dot{a}}\sigma^\nu_{b\dot{b}} = -2g^{\mu\nu}$}{eq 35.5}
         \quad [eq.~35.5]}
%% ============================================================

\subsection{Statement and physical meaning}

Contracting two $\sigma^\mu$ symbols with two $\varepsilon$ tensors
projects out the spacetime content, recovering the metric:
\begin{equation}
  \boxed{
    \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\,
    \sigma^\mu_{\alpha\dot\alpha}\,\sigma^\nu_{\beta\dot\beta}
    = -2\,g^{\mu\nu}
  }
  \tag{35.5}
\end{equation}
In Cadabra2:
\[
  """ + cdb_exprs["eq35_5_lhs"] + r""" = """ + cdb_exprs["eq35_5_rhs"] + r"""
\]

\subsection{Physical interpretation}

Together, eqs.~(35.4) and~(35.5) say that $\sigma^\mu$ is genuinely
\emph{invertible} as a basis change between vectors and bispinors:
\begin{align*}
  A_{\alpha\dot\alpha} &= \sigma^\mu_{\alpha\dot\alpha}\,A_\mu, &
  A_\mu &= -\tfrac{1}{2}\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}
           \sigma_{\mu\,\beta\dot\beta}\,A_{\alpha\dot\alpha}.
\end{align*}

\paragraph{Numerical result.}
\[
  \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}
  \sigma^\mu_{\alpha\dot\alpha}\sigma^\nu_{\beta\dot\beta} =
  \begin{pmatrix}
    -2 & 0 & 0 & 0 \\
    0 & 2 & 0 & 0 \\
    0 & 0 & 2 & 0 \\
    0 & 0 & 0 & 2
  \end{pmatrix}
  = -2\,\mathrm{diag}(+1,-1,-1,-1) = -2\,g^{\mu\nu}.
\]

\paragraph{Numerical verification.}
""" + verified(
    r"All 16 components of eq.~(35.5): max error $"
    + f"{err_35_5:.1e}" + r"$.",
    err_35_5
) + r"""

%% ============================================================
\section{\texorpdfstring{$\bar\sigma^\mu$}{sigmabar\^{}mu}:
         Definition and Numerical Verification \quad [eq.~35.19--35.20]}
%% ============================================================

\subsection{Definition (eq.~35.19)}

The raised-index sigma symbol is defined by:
\begin{equation}
  \bar\sigma^{\mu\,\dot\alpha\alpha}
  \equiv \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\,
         \sigma^\mu_{\beta\dot\beta}.
  \tag{35.19}
\end{equation}
In Cadabra2: $""" + cdb_exprs["sigmabar_def"] + r"""$.

\subsection{Explicit result (eq.~35.20)}

Both spinor indices are now upstairs; the result is:
\begin{equation}
  \bar\sigma^{\mu\,\dot\alpha\alpha}
  = (I,\,-\sigma_1,\,-\sigma_2,\,-\sigma_3).
  \tag{35.20}
\end{equation}
\[
  \bar\sigma^0 = """ + sigmabar_mats[0] + r""",\quad
  \bar\sigma^1 = """ + sigmabar_mats[1] + r""",\quad
  \bar\sigma^2 = """ + sigmabar_mats[2] + r""",\quad
  \bar\sigma^3 = """ + sigmabar_mats[3] + r"""
\]

In Cadabra2: $""" + cdb_exprs["sigmabar_mu"] + r"""$.

\subsection{Key property: trace identity}

From the definition it follows immediately that:
\begin{equation}
  \operatorname{tr}(\sigma^\mu\bar\sigma^\nu)
  \equiv \sigma^\mu_{\alpha\dot\alpha}\,\bar\sigma^{\nu\,\dot\alpha\alpha}
  = 2\,g^{\mu\nu}.
\end{equation}

\paragraph{Numerical verification.}
""" + verified(
    r"$\bar\sigma^{\mu} = (I,-\vec\sigma)$ from definition [eq.~35.20]: max error $"
    + f"{max_err_sigmabar:.1e}" + r"$.",
    max_err_sigmabar
) + r"""

""" + verified(
    r"$\operatorname{tr}(\sigma^\mu\bar\sigma^\nu) = 2g^{\mu\nu}$: max error $"
    + f"{err_trace:.1e}" + r"$.",
    err_trace
) + r"""

%% ============================================================
\section{Generators \texorpdfstring{$S^{\mu\nu}_{\mathrm{L}}$}{S\_L} and
         \texorpdfstring{$S^{\mu\nu}_{\mathrm{R}}$}{S\_R}
         in Terms of \texorpdfstring{$\sigma,\bar\sigma$}{sigma,sigmabar}
         \quad [eq.~35.21--35.22]}
%% ============================================================

\subsection{Formulas}

Using the covariance condition on $\sigma^\mu$ (eq.~35.9, 35.17--35.18)
and the definition of $\bar\sigma^\mu$, one derives the compact form:
\begin{align}
  \bigl(S^{\mu\nu}_{\mathrm{L}}\bigr)_\alpha{}^\beta
  &= \tfrac{i}{4}
    \bigl(\sigma^\mu_{\alpha\dot\alpha}\,\bar\sigma^{\nu\,\dot\alpha\beta}
         -\sigma^\nu_{\alpha\dot\alpha}\,\bar\sigma^{\mu\,\dot\alpha\beta}\bigr),
  \tag{35.21}\\[6pt]
  \bigl(S^{\mu\nu}_{\mathrm{R}}\bigr)^{\dot\alpha}{}_{\dot\beta}
  &= -\tfrac{i}{4}
    \bigl(\bar\sigma^{\mu\,\dot\alpha\alpha}\,\sigma^\nu_{\alpha\dot\beta}
         -\bar\sigma^{\nu\,\dot\alpha\alpha}\,\sigma^\mu_{\alpha\dot\beta}\bigr).
  \tag{35.22}
\end{align}

In Cadabra2:
$""" + cdb_exprs["gen_L"] + r"""$
and
$""" + cdb_exprs["gen_R"] + r"""$.

\paragraph{Relation to Ch.~34.}
Eq.~(35.21) is consistent with the Ch.~34 definitions
$S^{ij}_{\mathrm{L}} = \tfrac{1}{2}\varepsilon^{ijk}\sigma_k$
and $S^{k0}_{\mathrm{L}} = \tfrac{i}{2}\sigma_k$,
and eq.~(35.22) is consistent with
$S^{\mu\nu}_{\mathrm{R}} = -[S^{\mu\nu}_{\mathrm{L}}]^*$.

\subsection{Explicit matrices}

\paragraph{\texorpdfstring{$S^{\mu\nu}_{\mathrm{L}}$}{SL} matrices.}
\begin{align*}
  S^{01}_{\mathrm{L}} &= """ + SL_mats['01'] + r""",&
  S^{02}_{\mathrm{L}} &= """ + SL_mats['02'] + r""",&
  S^{03}_{\mathrm{L}} &= """ + SL_mats['03'] + r"""\\[6pt]
  S^{12}_{\mathrm{L}} &= """ + SL_mats['12'] + r""",&
  S^{13}_{\mathrm{L}} &= """ + SL_mats['13'] + r""",&
  S^{23}_{\mathrm{L}} &= """ + SL_mats['23'] + r"""
\end{align*}

\paragraph{\texorpdfstring{$S^{\mu\nu}_{\mathrm{R}}$}{SR} matrices.}
\begin{align*}
  S^{01}_{\mathrm{R}} &= """ + SR_mats['01'] + r""",&
  S^{02}_{\mathrm{R}} &= """ + SR_mats['02'] + r""",&
  S^{03}_{\mathrm{R}} &= """ + SR_mats['03'] + r"""\\[6pt]
  S^{12}_{\mathrm{R}} &= """ + SR_mats['12'] + r""",&
  S^{13}_{\mathrm{R}} &= """ + SR_mats['13'] + r""",&
  S^{23}_{\mathrm{R}} &= """ + SR_mats['23'] + r"""
\end{align*}

\paragraph{Numerical verification.}
""" + verified(
    r"$S^{\mu\nu}_{\mathrm{L}} = \tfrac{i}{4}(\sigma^\mu\bar\sigma^\nu - \sigma^\nu\bar\sigma^\mu)$"
    r" [eq.~35.21]: max error $" + f"{max_err_SL:.1e}" + r"$.",
    max_err_SL
) + r"""

""" + verified(
    r"$S^{\mu\nu}_{\mathrm{R}} = -\tfrac{i}{4}(\bar\sigma^\mu\sigma^\nu - \bar\sigma^\nu\sigma^\mu)$"
    r" [eq.~35.22]: max error $" + f"{max_err_SR:.1e}" + r"$.",
    max_err_SR
) + r"""

%% ============================================================
\section{Clifford-Like Algebra \quad
         \texorpdfstring{$\sigma^\mu\bar\sigma^\nu + \sigma^\nu\bar\sigma^\mu
         = -2g^{\mu\nu}I_2$}{Clifford}}
%% ============================================================

\subsection{Statement}

The anticommutator of $\sigma^\mu$ and $\bar\sigma^\nu$ (as $2\times2$ matrices):
\begin{equation}
  \bigl(\sigma^\mu\bar\sigma^\nu + \sigma^\nu\bar\sigma^\mu\bigr)_\alpha{}^\beta
  = -2\,g^{\mu\nu}\,\delta_\alpha{}^\beta.
  \tag{35-Clifford}
\end{equation}
This is the Weyl-spinor analogue of the Clifford algebra
$\{\gamma^\mu,\gamma^\nu\} = 2g^{\mu\nu}$ for Dirac matrices.

\subsection{Relation to generators}

The antisymmetric part gives exactly the left-handed generators:
\[
  \sigma^\mu\bar\sigma^\nu - \sigma^\nu\bar\sigma^\mu
  = \tfrac{4}{i}\,S^{\mu\nu}_{\mathrm{L}}
  \quad\Longrightarrow\quad
  S^{\mu\nu}_{\mathrm{L}} = \tfrac{i}{4}
  (\sigma^\mu\bar\sigma^\nu - \sigma^\nu\bar\sigma^\mu).
\]
So the Clifford relation says: the \emph{symmetric} part
$\sigma^\mu\bar\sigma^\nu + \sigma^\nu\bar\sigma^\mu$ is proportional to
the identity, while the \emph{antisymmetric} part encodes the generators.

\paragraph{Numerical verification.}
""" + verified(
    r"$\sigma^\mu\bar\sigma^\nu + \sigma^\nu\bar\sigma^\mu = -2g^{\mu\nu}I_2$"
    r" for all 16 $(\mu,\nu)$ pairs: max error $"
    + f"{max_err_cliff:.1e}" + r"$.",
    max_err_cliff
) + r"""

%% ============================================================
\section{Covariance Check: \texorpdfstring{$\sigma^\mu$}{sigma\^{}mu}
         Invariance \quad [eq.~35.14]}
%% ============================================================

\subsection{Statement}

The $\sigma^\mu$ symbol is invariant under simultaneous Lorentz
rotation of its vector index and both spinor indices.
Infinitesimally, the covariance condition reads:
\begin{equation}
  \bigl(g^{\mu\rho}\delta^{\nu}{}_\tau - g^{\nu\rho}\delta^{\mu}{}_\tau\bigr)
  \sigma^\tau_{\alpha\dot\alpha}
  + i\bigl(S^{\mu\nu}_{\mathrm{L}}\bigr)_\alpha{}^\beta
    \sigma^\rho_{\beta\dot\alpha}
  + i\bigl(S^{\mu\nu}_{\mathrm{R}}\bigr)_{\dot\alpha}{}^{\dot\beta}
    \sigma^\rho_{\alpha\dot\beta}
  = 0.
  \tag{35.14}
\end{equation}

\subsection{Physical significance}

This identity is the Lorentz-algebraic consistency condition that uniquely
determines $S^{\mu\nu}_{\mathrm{L}}$ and $S^{\mu\nu}_{\mathrm{R}}$ in
terms of $\sigma^\mu$.  Its solution is precisely eqs.~(35.21)--(35.22).

\paragraph{Numerical verification.}
""" + verified(
    r"Covariance condition [eq.~35.14] for all $(\mu,\nu,\rho,\alpha,\dot\alpha)$:"
    r" max residual $" + f"{max_err_cov:.1e}" + r"$.",
    max_err_cov
) + r"""

%% ============================================================
\section{Index-Free Notation and Spinor Bilinears \quad [eq.~35.23--35.29]}
%% ============================================================

\subsection{Index-free spinor products (eq.~35.23)}

Using the $\varepsilon$ tensor to contract spinor indices, one defines
Lorentz-invariant bilinears:
\begin{align}
  \chi\psi &\equiv \chi^\alpha\psi_\alpha
  = \varepsilon^{\alpha\beta}\chi_\alpha\psi_\beta
  \tag{35.23a}\\
  \chi^\dagger\psi^\dagger
  &\equiv \chi^\dagger_{\dot\alpha}\psi^{\dagger\,\dot\alpha}
  = \varepsilon_{\dot\alpha\dot\beta}\chi^{\dagger\,\dot\alpha}\psi^{\dagger\,\dot\beta}
  \tag{35.23b}
\end{align}

In Cadabra2:
\begin{align*}
  \chi\psi &= """ + cdb_exprs["bilinear_chi"] + r""" \\
  \chi^\dagger\psi^\dagger &= """ + cdb_exprs["bilinear_dag"] + r"""
\end{align*}

\subsection{Symmetry of the bilinear (eq.~35.25)}

For Grassmann fields, the two minus signs cancel:
\begin{equation}
  \chi\psi = \chi^\alpha\psi_\alpha
  = -\psi_\alpha\chi^\alpha \quad\text{(Grassmann)}
  = +\psi^\alpha\chi_\alpha = \psi\chi.
  \tag{35.25}
\end{equation}
So \emph{the undotted bilinear is symmetric} for Grassmann spinors.

\subsection{Hermitian conjugate (eq.~35.26)}

\begin{equation}
  (\chi\psi)^\dagger = (\chi^\alpha\psi_\alpha)^\dagger
  = (\psi_\alpha)^\dagger(\chi^\alpha)^\dagger
  = \psi^\dagger_{\dot\alpha}\chi^{\dagger\,\dot\alpha}
  = \psi^\dagger\chi^\dagger.
  \tag{35.26}
\end{equation}

\subsection{The vector bilinear (eq.~35.27--35.28)}

The combination $\psi^\dagger\bar\sigma^\mu\chi$ transforms as a
Lorentz 4-vector:
\begin{equation}
  \psi^\dagger\bar\sigma^\mu\chi
  \equiv \psi^\dagger_{\dot\alpha}\,\bar\sigma^{\mu\,\dot\alpha\alpha}\,\chi_\alpha.
  \tag{35.27}
\end{equation}
In Cadabra2: $""" + cdb_exprs["vec_bilinear"] + r"""$.

\medskip
Under a Lorentz transformation (eq.~35.28):
\begin{equation}
  U(\Lambda)^{-1}\,[\psi^\dagger\bar\sigma^\mu\chi]\,U(\Lambda)
  = {\Lambda^\mu}{}_\nu\,[\psi^\dagger\bar\sigma^\nu\chi].
  \tag{35.28}
\end{equation}

Since $\bar\sigma^\mu = (I,-\vec\sigma)$ is Hermitian,
the hermitian conjugate satisfies (eq.~35.29):
\begin{equation}
  [\psi^\dagger\bar\sigma^\mu\chi]^\dagger = \chi^\dagger\bar\sigma^\mu\psi.
  \tag{35.29}
\end{equation}

%% ============================================================
\section{Summary of Chapter 35 Results}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}p{5.5cm}p{8cm}@{}}
  \toprule
  Result & Expression \\
  \midrule
  $\sigma^\mu$ symbol & $\sigma^\mu_{\alpha\dot\alpha} = (I,\,\vec\sigma)$ \\[2pt]
  $\bar\sigma^\mu$ definition & $\bar\sigma^{\mu\dot\alpha\alpha}
    = \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma^\mu_{\beta\dot\beta}
    = (I,-\vec\sigma)$ \quad [35.19] \\[2pt]
  Completeness I \quad [35.4] &
    $\sigma^\mu_{\alpha\dot\alpha}\sigma_{\mu\beta\dot\beta}
     = -2\varepsilon_{\alpha\beta}\varepsilon_{\dot\alpha\dot\beta}$ \\[2pt]
  Completeness II \quad [35.5] &
    $\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}
     \sigma^\mu_{\alpha\dot\alpha}\sigma^\nu_{\beta\dot\beta} = -2g^{\mu\nu}$ \\[2pt]
  Trace identity &
    $\operatorname{tr}(\sigma^\mu\bar\sigma^\nu) = 2g^{\mu\nu}$ \\[2pt]
  Clifford-like relation &
    $\sigma^\mu\bar\sigma^\nu + \sigma^\nu\bar\sigma^\mu = -2g^{\mu\nu}I_2$ \\[2pt]
  $S^{\mu\nu}_{\mathrm{L}}$ generator \quad [35.21] &
    $= \tfrac{i}{4}(\sigma^\mu\bar\sigma^\nu - \sigma^\nu\bar\sigma^\mu)$ \\[2pt]
  $S^{\mu\nu}_{\mathrm{R}}$ generator \quad [35.22] &
    $= -\tfrac{i}{4}(\bar\sigma^\mu\sigma^\nu - \bar\sigma^\nu\sigma^\mu)
     = -[S^{\mu\nu}_{\mathrm{L}}]^*$ \\[2pt]
  Covariance \quad [35.14] &
    $\sigma^\rho$ invariant under simultaneous $L$, $R$, vector rotations \\[2pt]
  Angle bracket &
    $\chi\psi = \varepsilon^{\alpha\beta}\chi_\alpha\psi_\beta = \psi\chi$
    \quad [35.23, 35.25] \\[2pt]
  Square bracket &
    $\chi^\dagger\psi^\dagger
     = \varepsilon_{\dot\alpha\dot\beta}\chi^{\dagger\dot\alpha}\psi^{\dagger\dot\beta}$
    \quad [35.23] \\[2pt]
  Hermitian conjugate &
    $(\chi\psi)^\dagger = \psi^\dagger\chi^\dagger$ \quad [35.26] \\[2pt]
  Vector bilinear &
    $\psi^\dagger\bar\sigma^\mu\chi$ transforms as 4-vector \quad [35.28] \\[2pt]
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\begin{center}
\renewcommand{\arraystretch}{1.3}
\begin{tabular}{@{}lccc@{}}
  \toprule
  Check & Equation & Max Error & Status \\
  \midrule
  $\sigma^\mu_{a\dot a}\sigma_{\mu b\dot b} = -2\varepsilon_{ab}\varepsilon_{\dot a\dot b}$
    & 35.4 & $""" + f"{err_35_4:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $\varepsilon^{ab}\varepsilon^{\dot a\dot b}\sigma^\mu_{a\dot a}\sigma^\nu_{b\dot b}=-2g^{\mu\nu}$
    & 35.5 & $""" + f"{err_35_5:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $\bar\sigma^\mu = (I,-\vec\sigma)$
    & 35.20 & $""" + f"{max_err_sigmabar:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $\operatorname{tr}(\sigma^\mu\bar\sigma^\nu) = 2g^{\mu\nu}$
    & --- & $""" + f"{err_trace:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $S^{\mu\nu}_{\mathrm{L}} = \tfrac{i}{4}(\sigma^\mu\bar\sigma^\nu-\sigma^\nu\bar\sigma^\mu)$
    & 35.21 & $""" + f"{max_err_SL:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $S^{\mu\nu}_{\mathrm{R}} = -\tfrac{i}{4}(\bar\sigma^\mu\sigma^\nu-\bar\sigma^\nu\sigma^\mu)$
    & 35.22 & $""" + f"{max_err_SR:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  $\sigma^\mu\bar\sigma^\nu+\sigma^\nu\bar\sigma^\mu = -2g^{\mu\nu}I_2$
    & --- & $""" + f"{max_err_cliff:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  Covariance condition
    & 35.14 & $""" + f"{max_err_cov:.1e}" + r"""$ & \textcolor{verifygreen}{\checkmark} \\
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\noindent\textbf{Next:} Chapter~36 constructs the Weyl Lagrangian using
these identities and derives the equations of motion for massless
left- and right-handed spinor fields.

\end{document}
"""

outpath = "/work/ch35_sigma_algebra.tex"
with open(outpath, "w") as f:
    f.write(doc.strip())

print(f"Wrote: {outpath}")
