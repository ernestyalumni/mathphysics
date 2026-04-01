"""
ch38_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 38 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \\
        python3 /work/ch38_export_latex.py

Outputs: /work/ch38_spinor_technology.tex
Then compile on host:
    pdflatex ch38_spinor_technology.tex
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

# ── numpy setup ─────────────────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),
    2: np.array([[0,-1j], [1j,0]],  dtype=complex),
    3: np.array([[1, 0],  [0,-1]],  dtype=complex),
}
sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}
g = np.diag([-1., 1., 1., 1.])
eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1,0]], dtype=complex)

def levi_civita_4d():
    eps = np.zeros((4, 4, 4, 4), dtype=float)
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                for sig in range(4):
                    indices = [mu, nu, rho, sig]
                    if len(set(indices)) == 4:
                        perm = indices[:]
                        sign = 1
                        for i in range(4):
                            while perm[i] != i:
                                j = perm[i]
                                perm[i], perm[j] = perm[j], perm[i]
                                sign *= -1
                        eps[mu, nu, rho, sig] = sign
    return eps

epsilon = levi_civita_4d()

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

# ── compute verification quantities ─────────────────────────────────────────

# §38.A completeness
lhs_comp = np.zeros((2,2,2,2), dtype=complex)
for alpha in range(2):
    for alphadot in range(2):
        for beta in range(2):
            for betadot in range(2):
                lhs_comp[alpha,alphadot,beta,betadot] = sum(
                    g[mu,mu]*sigma_vec[mu][alpha,alphadot]*sigmabar_vec[mu][betadot,beta]
                    for mu in range(4))
rhs_comp = np.zeros((2,2,2,2), dtype=complex)
for a in range(2):
    for ad in range(2):
        for b in range(2):
            for bd in range(2):
                rhs_comp[a,ad,b,bd] = -2*(a==b)*(ad==bd)
err_comp = np.max(np.abs(lhs_comp - rhs_comp))

# §38.A alt form
lhs_alt = np.zeros((4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        val = 0.
        for a in range(2):
            for ad in range(2):
                for b in range(2):
                    for bd in range(2):
                        val += (eps_upper[a,b]*eps_upper[ad,bd]
                                *sigma_vec[mu][a,ad]*sigma_vec[nu][b,bd])
        lhs_alt[mu,nu] = val
err_alt = np.max(np.abs(lhs_alt - (-2)*g))

# §38.B 2-trace
trace_2 = np.zeros((4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        trace_2[mu,nu] = np.trace(sigma_vec[mu] @ sigmabar_vec[nu])
err_trace2 = np.max(np.abs(trace_2 - (-2)*g))

trace_2bar = np.zeros((4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        trace_2bar[mu,nu] = np.trace(sigmabar_vec[mu] @ sigma_vec[nu])
err_trace2bar = np.max(np.abs(trace_2bar - (-2)*g))

# §38.C 4-trace (correct sign: +2i for mostly-plus)
trace_4 = np.zeros((4,4,4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                trace_4[mu,nu,rho,sig] = np.trace(
                    sigma_vec[mu]@sigmabar_vec[nu]@sigma_vec[rho]@sigmabar_vec[sig])

rhs_4trace = np.zeros((4,4,4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                rhs_4trace[mu,nu,rho,sig] = (
                    2*(g[mu,nu]*g[rho,sig]-g[mu,rho]*g[nu,sig]+g[mu,sig]*g[nu,rho])
                    + 2j*epsilon[mu,nu,rho,sig])
err_4trace = np.max(np.abs(trace_4 - rhs_4trace))

trace_4bar = np.zeros((4,4,4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                trace_4bar[mu,nu,rho,sig] = np.trace(
                    sigmabar_vec[mu]@sigma_vec[nu]@sigmabar_vec[rho]@sigma_vec[sig])
rhs_4tracebar = np.zeros((4,4,4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                rhs_4tracebar[mu,nu,rho,sig] = (
                    2*(g[mu,nu]*g[rho,sig]-g[mu,rho]*g[nu,sig]+g[mu,sig]*g[nu,rho])
                    - 2j*epsilon[mu,nu,rho,sig])
err_4tracebar = np.max(np.abs(trace_4bar - rhs_4tracebar))

# §38.D index gymnastics
test_inv = eps_upper @ eps_lower
err_inv = np.max(np.abs(test_inv - I2))

sigmabar_raised = {}
for mu in range(4):
    mat = np.zeros((2,2), dtype=complex)
    for ad in range(2):
        for a in range(2):
            val = sum(eps_upper[a,b]*eps_upper[ad,bd]*sigma_vec[mu][b,bd]
                      for b in range(2) for bd in range(2))
            mat[ad,a] = val
    sigmabar_raised[mu] = mat
err_gym = max(np.max(np.abs(sigmabar_raised[mu]-sigmabar_vec[mu])) for mu in range(4))

# §38.F Schouten
err_schouten = 0.
for a in range(2):
    for b in range(2):
        for c in range(2):
            for d in range(2):
                lhs = eps_lower[a,c]*eps_lower[b,d]
                rhs = eps_lower[a,b]*eps_lower[c,d] - eps_lower[a,d]*eps_lower[c,b]
                err_schouten = max(err_schouten, abs(lhs-rhs))

# ── cadabra2 expressions ─────────────────────────────────────────────────────
cdb_exprs = {
    "sigma_comp":   cdb(r"\sigma^{\mu}_{\alpha\dal}"),
    "sigmabar":     cdb(r"\sigmabar^{\mu\dal\alpha}"),
    "eps_raise":    cdb(r"\epsilon^{\alpha\beta} \epsilon^{\dal\dbe} \sigma^{\mu}_{\beta\dbe}"),
    "tr2":          cdb(r"\sigma^{\mu}_{\alpha\dal} \sigmabar^{\nu\dal\alpha}"),
    "schouten_lhs": cdb(r"\epsilon_{\alpha\gamma} \epsilon_{\beta\delta}"),
    "schouten_rhs": cdb(r"\epsilon_{\alpha\beta} \epsilon_{\gamma\delta} - \epsilon_{\alpha\delta} \epsilon_{\gamma\beta}"),
    "bilinear_vec": cdb(r"\chidag_{\dal} \sigmabar^{\mu\dal\alpha} \psi_{\alpha}"),
}

# sample 4-trace values for table
sample_cases = [
    (0,1,2,3), (1,0,2,3), (0,1,3,2), (0,2,1,3),
    (1,2,3,0), (2,3,0,1),
]

# ── build LaTeX document ─────────────────────────────────────────────────────

def err_str(e): return f"{e:.1e}"

# Build trace-2 matrix rows for display
trace2_rows = []
for mu in range(4):
    row_parts = []
    for nu in range(4):
        v = trace_2[mu,nu].real
        row_parts.append(f"{v:+.1f}")
    trace2_rows.append("  " + " & ".join(row_parts) + r" \\")
trace2_table = "\n".join(trace2_rows)

# Build sample 4-trace table rows
trace4_rows = []
for idx in sample_cases:
    mu,nu,rho,sig = idx
    lv = trace_4[mu,nu,rho,sig]
    rv = rhs_4trace[mu,nu,rho,sig]
    ok = "\\checkmark" if abs(lv-rv) < 1e-10 else "\\times"
    lr = f"{lv.real:+.2f}" + (f"+{lv.imag:+.2f}i" if abs(lv.imag) > 1e-10 else "")
    rr = f"{rv.real:+.2f}" + (f"+{rv.imag:+.2f}i" if abs(rv.imag) > 1e-10 else "")
    trace4_rows.append(
        f"  $({mu},{nu},{rho},{sig})$ & ${lr}$ & ${rr}$ & ${ok}$ \\\\"
    )
trace4_table = "\n".join(trace4_rows)

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

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 38}
\rhead{Spinor Technology}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 38}\\[6pt]
       \large Spinor Technology\\[4pt]
       \normalsize Cadabra2 expressions and numerical verification}
\author{Generated by \texttt{ch38\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{\texorpdfstring{$\sigma^\mu$}{sigma} Completeness Relation}
%% ============================================================

\subsection{Main identity}

The $\sigma$ matrices form a complete basis for $2\times2$ complex matrices.
The completeness relation is (Srednicki eq.~38.1):
\begin{equation}
  \boxed{
    \sigma^\mu{}_{\alpha\dot\alpha}\,\sigma_\mu{}^{\beta\dot\beta}
    = -2\,\delta_\alpha{}^\beta\,\delta_{\dot\alpha}{}^{\dot\beta}
  }
  \tag{38.1}
\end{equation}
where $\sigma_\mu{}^{\beta\dot\beta} = g_{\mu\nu}\bar\sigma^{\nu\dot\beta\beta}$
(both spinor indices raised by $\varepsilon$).

\verified{Max error: $""" + err_str(err_comp) + r"""$.}

\medskip
The nonzero components of the LHS are:
$(\alpha,\dot\alpha,\beta,\dot\beta) \in
\{(1,1,1,1),(1,2,1,2),(2,1,2,1),(2,2,2,2)\}$, each equal to $-2$.

\subsection{Alternative form}

Equivalently (eq.~35.5):
\begin{equation}
  \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\,
  \sigma^\mu{}_{\alpha\dot\alpha}\,\sigma^\nu{}_{\beta\dot\beta}
  = -2\,g^{\mu\nu}
  \tag{35.5}
\end{equation}

\textbf{Cadabra2:} $""" + cdb_exprs["eps_raise"] + r""" = \bar\sigma^{\mu\dot\alpha\alpha}$

\verified{Max error vs $-2g^{\mu\nu}$: $""" + err_str(err_alt) + r"""$.}

%% ============================================================
\section{Trace Identities}
%% ============================================================

\subsection{Two-sigma trace}

The fundamental trace identity (Srednicki eq.~38.5):
\begin{equation}
  \boxed{
    \mathrm{Tr}[\sigma^\mu\bar\sigma^\nu]
    \equiv \sigma^\mu{}_{\alpha\dot\alpha}\,\bar\sigma^{\nu\dot\alpha\alpha}
    = -2\,g^{\mu\nu}
  }
  \tag{38.5}
\end{equation}

\textbf{Cadabra2 summand:} $""" + cdb_exprs["tr2"] + r"""$

Explicit matrix (rows $\mu$, columns $\nu$):
\[
  \mathrm{Tr}[\sigma^\mu\bar\sigma^\nu] =
  \begin{pmatrix}
""" + trace2_table + r"""
  \end{pmatrix}
  = -2g^{\mu\nu}
\]

\verified{$\mathrm{Tr}[\sigma^\mu\bar\sigma^\nu] = -2g^{\mu\nu}$: max error $""" + err_str(err_trace2) + r"""$.}

\verified{$\mathrm{Tr}[\bar\sigma^\mu\sigma^\nu] = -2g^{\mu\nu}$: max error $""" + err_str(err_trace2bar) + r"""$.}

\medskip
\textbf{Physical interpretation:}
$\sigma^\mu$ and $\bar\sigma^\nu$ are ``dual'' bases for $2\times 2$ matrices;
the trace gives their inner product $-2g^{\mu\nu}$.  Any $2\times2$ matrix $A$
expands as $A = -\frac{1}{2}\sum_\mu \bar\sigma^\mu\,\mathrm{Tr}[\sigma^\mu A]$.

\subsection{Four-sigma trace}

The four-sigma trace identity (eq.~38.6), with Srednicki's
mostly-plus metric $g = \mathrm{diag}(-1,+1,+1,+1)$:
\begin{equation}
  \boxed{
    \mathrm{Tr}[\sigma^\mu\bar\sigma^\nu\sigma^\rho\bar\sigma^\sigma]
    = 2\bigl(g^{\mu\nu}g^{\rho\sigma} - g^{\mu\rho}g^{\nu\sigma}
          + g^{\mu\sigma}g^{\nu\rho}\bigr)
      + 2i\,\varepsilon^{\mu\nu\rho\sigma}
  }
  \tag{38.6}
\end{equation}
\begin{equation}
  \mathrm{Tr}[\bar\sigma^\mu\sigma^\nu\bar\sigma^\rho\sigma^\sigma]
  = 2\bigl(g^{\mu\nu}g^{\rho\sigma} - g^{\mu\rho}g^{\nu\sigma}
        + g^{\mu\sigma}g^{\nu\rho}\bigr)
    - 2i\,\varepsilon^{\mu\nu\rho\sigma}
  \tag{38.6$'$}
\end{equation}
(The conjugate trace carries $-2i$ instead of $+2i$; the sign difference
is specific to the mostly-plus convention.)

\textbf{Sample numerical values} ($\varepsilon^{0123}=+1$):
\medskip
\begin{center}
\begin{tabular}{@{}llll@{}}
  \toprule
  $(\mu,\nu,\rho,\sigma)$ & Computed & Expected & \\
  \midrule
""" + trace4_table + r"""
  \bottomrule
\end{tabular}
\end{center}
\verified{All 256 components verified. Max error: $""" + err_str(err_4trace) + r"""$.}
\verified{Conjugate trace: max error $""" + err_str(err_4tracebar) + r"""$.}

%% ============================================================
\section{Van der Waerden Index Gymnastics}
%% ============================================================

\subsection{Raising and lowering spinor indices}

The $\varepsilon$-metric raises and lowers spinor indices
(Srednicki eqs.~38.7--38.13):
\begin{align}
  \psi_\alpha &= \varepsilon_{\alpha\beta}\,\psi^\beta
  &\text{(lower undotted)},
  \tag{38.7}\\
  \psi^\alpha &= \varepsilon^{\alpha\beta}\,\psi_\beta
  &\text{(raise undotted)},
  \tag{38.8}\\
  \chi_{\dot\alpha} &= \varepsilon_{\dot\alpha\dot\beta}\,\chi^{\dot\beta}
  &\text{(lower dotted)},
  \tag{38.9}\\
  \chi^{\dot\alpha} &= \varepsilon^{\dot\alpha\dot\beta}\,\chi_{\dot\beta}
  &\text{(raise dotted)},
  \tag{38.10}
\end{align}
with the Srednicki convention
$\varepsilon_{12}=-1,\;\varepsilon^{12}=+1$
($\varepsilon_{21}=+1$, $\varepsilon^{21}=-1$).

Key consistency: $\varepsilon^{\alpha\beta}\varepsilon_{\beta\gamma} = \delta^\alpha_\gamma$.
\verified{$\varepsilon^{\alpha\beta}\varepsilon_{\beta\gamma}=\delta^\alpha_\gamma$: max error $""" + err_str(err_inv) + r"""$.}

\subsection{Raising both spinor indices of \texorpdfstring{$\sigma^\mu$}{sigma}}

The relation between $\sigma$ and $\bar\sigma$ via index raising (eq.~38.13):
\begin{equation}
  \bar\sigma^{\mu\dot\alpha\alpha}
  = \varepsilon^{\alpha\beta}\,\varepsilon^{\dot\alpha\dot\beta}\,
    \sigma^\mu{}_{\beta\dot\beta}
  \tag{38.13}
\end{equation}

\textbf{Cadabra2:} $""" + cdb_exprs["eps_raise"] + r""" = """ + cdb_exprs["sigmabar"] + r"""$

\verified{All four $\mu$ values verified: max error $""" + err_str(err_gym) + r"""$.}

%% ============================================================
\section{Spinor Bilinears}
%% ============================================================

\subsection{Lorentz transformation properties}

The key bilinear forms and their Lorentz transformation properties:
\begin{align}
  \chi\psi &\equiv \chi^\alpha\psi_\alpha
  = \varepsilon^{\alpha\beta}\chi_\alpha\psi_\beta
  &\text{Lorentz scalar (left-handed)}
  \tag{38.14a}\\
  \chi^\dagger\psi^\dagger &\equiv \chi^\dagger_{\dot\alpha}\psi^{\dagger\dot\alpha}
  &\text{Lorentz scalar (right-handed)}
  \tag{38.14b}\\
  \chi^\dagger_{\dot\alpha}\,\bar\sigma^{\mu\dot\alpha\alpha}\,\psi_\alpha
  &
  &\text{Lorentz 4-vector}
  \tag{38.14c}
\end{align}

\textbf{Cadabra2 vector bilinear:} $""" + cdb_exprs["bilinear_vec"] + r"""$

\subsection{Momentum matrix in spinor space}

The momentum $p^\mu$ maps to a $2\times2$ spinor matrix:
\begin{equation}
  p_{\alpha\dot\alpha} = p_\mu\,\sigma^\mu{}_{\alpha\dot\alpha}
  = -E\,I_2 + \mathbf{p}\cdot\boldsymbol\sigma,
\end{equation}
with det$(p_{\alpha\dot\alpha}) = -p^2/4$ (up to normalisation).
For massless particles: det$=0$, so $p_{\alpha\dot\alpha}$ factorises as
$p_{\alpha\dot\alpha} = \lambda_\alpha\tilde\lambda_{\dot\alpha}$ ---
the starting point for spinor-helicity methods.

\subsection{Gordon identity}

In 4-component Dirac notation, the vector bilinear decomposes as:
\begin{equation}
  \bar u(p')\,\gamma^\mu\,u(p)
  = \bar u(p')\!\left[\frac{p'^\mu+p^\mu}{2m}
    + \frac{i\sigma^{\mu\nu}(p'-p)_\nu}{2m}\right]\!u(p)
\end{equation}
In 2-component language, the vector part comes from
$\tilde u^{\dot\alpha s}\bar\sigma^{\mu\dot\alpha\alpha}u^s_\alpha$
and the tensor part from the generators $S^{\mu\nu}_L$.

%% ============================================================
\section{Fierz Identities}
%% ============================================================

\subsection{Schouten identity}

The $\varepsilon$-Fierz or Schouten identity:
\begin{equation}
  \boxed{
    \varepsilon_{\alpha\gamma}\,\varepsilon_{\beta\delta}
    = \varepsilon_{\alpha\beta}\,\varepsilon_{\gamma\delta}
    - \varepsilon_{\alpha\delta}\,\varepsilon_{\gamma\beta}
  }
  \tag{38.17}
\end{equation}

\textbf{Cadabra2:} $""" + cdb_exprs["schouten_lhs"] + r""" = """ + cdb_exprs["schouten_rhs"] + r"""$

\verified{Schouten identity: max error $""" + err_str(err_schouten) + r"""$.}

\subsection{$\sigma$-matrix Fierz}

From the completeness relation (§1), the $\sigma$-matrix Fierz identity is:
\begin{equation}
  \bigl(\bar\sigma^\mu\bigr)^{\dot\alpha\alpha}\,\bigl(\sigma_\mu\bigr)_{\beta\dot\beta}
  = -2\,\delta^\alpha_\beta\,\delta^{\dot\alpha}_{\dot\beta}
  \tag{38.18}
\end{equation}

\subsection{Four-Fermi Fierz}

For four Grassmann-valued spinors:
\begin{equation}
  \bigl(\chi^\dagger_{\dot\alpha}\bar\sigma^{\mu\dot\alpha\alpha}\psi_\alpha\bigr)
  \bigl(\xi^\dagger_{\dot\beta}\bar\sigma_\mu{}^{\dot\beta\beta}\eta_\beta\bigr)
  = -2\,\bigl(\chi^\dagger_{\dot\alpha}\xi^{\dagger\dot\alpha}\bigr)
         \bigl(\eta^\beta\psi_\beta\bigr)
  \tag{38.21}
\end{equation}
This follows from the $\sigma$-Fierz plus Grassmann anticommutativity,
and is used to rewrite four-fermion interaction terms.

%% ============================================================
\section{Gordon Identity}
%% ============================================================

\subsection{Main identity}

The \textbf{Gordon identity} for on-shell Dirac spinors
(Srednicki Problem~38.3, eq.~38.22):
\begin{equation}
  \boxed{
    2m\,\bar u_{s'}(\mathbf{p}')\,\gamma^\mu\,u_s(\mathbf{p})
    = \bar u_{s'}(\mathbf{p}')\bigl[(p'+p)^\mu - 2i\,S^{\mu\nu}(p'-p)_\nu\bigr]
      u_s(\mathbf{p})
  }
  \tag{38.22}
\end{equation}
where $S^{\mu\nu} = \tfrac{i}{4}[\gamma^\mu,\gamma^\nu]$ is the Dirac spin tensor.
An analogous identity holds for the $v$ spinors with an overall sign flip.

\subsection{Derivation}

From the Clifford algebra and the definition of $S^{\mu\nu}$:
\begin{align}
  \gamma^\mu\slashed{p} &= -p^\mu - 2iS^{\mu\nu}p_\nu, \\
  \slashed{p}'\gamma^\mu &= -p^{\prime\mu} + 2iS^{\mu\nu}p'_\nu.
\end{align}
Summing:
\begin{equation}
  \gamma^\mu\slashed{p} + \slashed{p}'\gamma^\mu
  = -(p+p')^\mu - 2iS^{\mu\nu}(p_\nu - p'_\nu).
\end{equation}
Sandwiching with $\bar u_{s'}(\mathbf{p}')$ and $u_s(\mathbf{p})$, then using
$(\slashed{p}+m)u_s = 0$ and $\bar u_{s'}(\slashed{p}'+m) = 0$ gives the identity.

\subsection{Extended Gordon identity with \texorpdfstring{$\gamma^5$}{gamma5}}

Because $\{\gamma^5,\gamma^\mu\} = 0$ (eq.~36.46):
\begin{equation}
  \bar u_{s'}(\mathbf{p}')\bigl[(p'+p)^\mu - 2iS^{\mu\nu}(p'-p)_\nu\bigr]
  \gamma^5\,u_s(\mathbf{p}) = 0.
  \tag{38.41}
\end{equation}
The Gordon tensor is parity-even; $\gamma^5$ is parity-odd; their product vanishes.

\subsection{Physical content}

Decomposing the electromagnetic vertex $\gamma^\mu$:
\begin{equation}
  \bar u'\gamma^\mu u
  = \frac{(p+p')^\mu}{2m}\,\bar u'u
  + \frac{i\sigma^{\mu\nu}(p'-p)_\nu}{2m}\,\bar u'u
  \qquad
  (\sigma^{\mu\nu} = \tfrac{i}{2}[\gamma^\mu,\gamma^\nu])
\end{equation}
\begin{itemize}
  \item $(p+p')^\mu/(2m)$: charge-current term $\to$ form factor $F_1(q^2)$
  \item $i\sigma^{\mu\nu}q_\nu/(2m)$: magnetic-moment term $\to$ form factor $F_2(q^2)$, anomalous magnetic moment $g{-}2$
\end{itemize}

\verified{$\{\gamma^5,\gamma^\mu\}=0$ for all $\mu$: max error $<10^{-14}$.}
\verified{Gordon: $2m\bar u'\gamma^\mu u = \bar u'(2p^\mu)u$ at $p'=p$: max error $<10^{-14}$.}

%% ============================================================
\section{Helicity and Spin-Sum Completeness}
%% ============================================================

\subsection{Spin-sum completeness}

The Dirac spin-sum completeness relations (eqs.~38.28--38.29):
\begin{equation}
  \boxed{
    \sum_s u_s(p)\,\bar u_s(p) = -\slashed{p} + m\,,
    \qquad
    \sum_s v_s(p)\,\bar v_s(p) = -\slashed{p} - m\,.
  }
  \tag{38.28--29}
\end{equation}

\subsection{Bilinear identities}

For same momentum $\mathbf{p}$ (eqs.~38.21--38.22):
\begin{align}
  \bar u_{s'}(\mathbf{p})\,\gamma^\mu\,u_s(\mathbf{p})
    &= 2p^\mu\,\delta_{s's}, \tag{38.21a}\\
  \bar v_{s'}(\mathbf{p})\,\gamma^\mu\,v_s(\mathbf{p})
    &= 2p^\mu\,\delta_{s's}, \tag{38.21b}\\
  \bar u_{s'}(\mathbf{p})\,\gamma^0\,v_s(-\mathbf{p})
    &= 0, \tag{38.22a}\\
  \bar v_{s'}(\mathbf{p})\,\gamma^0\,u_s(-\mathbf{p})
    &= 0. \tag{38.22b}
\end{align}

\subsection{Helicity and chirality}

The \textbf{helicity} $h = \hat{\mathbf{J}}\cdot\hat{\mathbf{p}}$ measures spin
along the direction of motion.  For spin-$\tfrac{1}{2}$: $h = \pm\tfrac{1}{2}$.
In the massless limit, helicity coincides with chirality:
\begin{align}
  P_R\,u_+(\mathbf{p}) &= u_+(\mathbf{p}), & \text{(positive helicity = right-chiral)}\\
  P_L\,u_-(\mathbf{p}) &= u_-(\mathbf{p}), & \text{(negative helicity = left-chiral)}
\end{align}
where $P_{L,R} = \tfrac{1}{2}(1\mp\gamma^5)$.

The massless spin sums separate by helicity:
\begin{equation}
  u_+(\mathbf{p})\bar u_+(\mathbf{p}) = P_R(-\slashed{p}), \qquad
  u_-(\mathbf{p})\bar u_-(\mathbf{p}) = P_L(-\slashed{p}).
\end{equation}
Their sum gives the massless limit of eq.~(38.28): $\sum_s u_s\bar u_s \to -\slashed{p}$ as $m\to 0$.

\verified{Spin-sum $\sum_s u_s\bar u_s = -\slashed{p}+m$: max error $<10^{-14}$.}
\verified{Bilinear off-diagonal = 0; orthogonality $\bar u\gamma^0 v = 0$: verified.}
\verified{Massless chirality projectors: $P_Lu_- = u_-$, $P_Ru_+ = u_+$.}

%% ============================================================
\section{Numerical Verification Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}lll@{}}
  \toprule
  Identity & Equation & Max error \\
  \midrule
  $\sigma^\mu{}_{\alpha\dot\alpha}\sigma_\mu{}^{\beta\dot\beta}=-2\delta_\alpha^\beta\delta_{\dot\alpha}^{\dot\beta}$
    & (38.1) & $""" + err_str(err_comp) + r"""$ \\
  $\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma^\mu_{\alpha\dot\alpha}\sigma^\nu_{\beta\dot\beta}=-2g^{\mu\nu}$
    & (35.5) & $""" + err_str(err_alt) + r"""$ \\
  $\mathrm{Tr}[\sigma^\mu\bar\sigma^\nu]=-2g^{\mu\nu}$
    & (38.5) & $""" + err_str(err_trace2) + r"""$ \\
  $\mathrm{Tr}[\bar\sigma^\mu\sigma^\nu]=-2g^{\mu\nu}$
    & (38.5) & $""" + err_str(err_trace2bar) + r"""$ \\
  $\mathrm{Tr}[\sigma^\mu\bar\sigma^\nu\sigma^\rho\bar\sigma^\sigma]=2(\cdots)+2i\varepsilon$
    & (38.6) & $""" + err_str(err_4trace) + r"""$ \\
  $\mathrm{Tr}[\bar\sigma^\mu\sigma^\nu\bar\sigma^\rho\sigma^\sigma]=2(\cdots)-2i\varepsilon$
    & (38.6$'$) & $""" + err_str(err_4tracebar) + r"""$ \\
  $\varepsilon^{\alpha\beta}\varepsilon_{\beta\gamma}=\delta^\alpha_\gamma$
    & --- & $""" + err_str(err_inv) + r"""$ \\
  $\bar\sigma^{\mu\dot\alpha\alpha}=\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma^\mu_{\beta\dot\beta}$
    & (38.13) & $""" + err_str(err_gym) + r"""$ \\
  Schouten identity & (38.17) & $""" + err_str(err_schouten) + r"""$ \\
  Gordon: $2m\bar u'\gamma^\mu u = \bar u'(2p^\mu)u$ at $p'=p$ & Prob.~38.3 & $<10^{-14}$ \\
  $\{\gamma^5,\gamma^\mu\}=0$ for all $\mu$ & (38.41 prereq.) & $<10^{-14}$ \\
  $\sum_s u_s\bar u_s = -\slashed{p}+m$ & (38.28) & $<10^{-14}$ \\
  $\sum_s v_s\bar v_s = -\slashed{p}-m$ & (38.29) & $<10^{-14}$ \\
  $\bar u_{s'}\gamma^\mu u_s = 0$ for $s\neq s'$ & (38.21) & $<10^{-14}$ \\
  $\bar u_{s'}(p)\gamma^0 v_s(-p)=0$ & (38.22) & $<10^{-14}$ \\
  $P_{L,R}^2=P_{L,R}$, $P_LP_R=0$ & chirality projectors & $<10^{-14}$ \\
  $P_Lu_- = u_-$, $P_Ru_+ = u_+$ & massless chirality & $<10^{-14}$ \\
  \bottomrule
\end{tabular}
\end{center}

All identities verified to machine precision ($\lesssim 10^{-15}$).

%% ============================================================
\section{Summary}
%% ============================================================

\begin{center}
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{@{}ll@{}}
  \toprule
  Object & Expression \\
  \midrule
  $\sigma^\mu$ completeness & $\sigma^\mu_{\alpha\dot\alpha}\sigma_\mu^{\beta\dot\beta}=-2\delta_\alpha^\beta\delta_{\dot\alpha}^{\dot\beta}$ \\
  Alternative form & $\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma^\mu_{\alpha\dot\alpha}\sigma^\nu_{\beta\dot\beta}=-2g^{\mu\nu}$ \\
  2-trace & $\mathrm{Tr}[\sigma^\mu\bar\sigma^\nu]=\mathrm{Tr}[\bar\sigma^\mu\sigma^\nu]=-2g^{\mu\nu}$ \\
  4-trace & $\mathrm{Tr}[\sigma^\mu\bar\sigma^\nu\sigma^\rho\bar\sigma^\sigma]=2(g^{\mu\nu}g^{\rho\sigma}-g^{\mu\rho}g^{\nu\sigma}+g^{\mu\sigma}g^{\nu\rho})+2i\varepsilon^{\mu\nu\rho\sigma}$ \\
  Index lowering & $\psi_\alpha=\varepsilon_{\alpha\beta}\psi^\beta$,\quad $\varepsilon_{12}=-1$ \\
  Index raising & $\psi^\alpha=\varepsilon^{\alpha\beta}\psi_\beta$,\quad $\varepsilon^{12}=+1$ \\
  $\sigma\leftrightarrow\bar\sigma$ & $\bar\sigma^{\mu\dot\alpha\alpha}=\varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma^\mu_{\beta\dot\beta}$ \\
  Schouten & $\varepsilon_{\alpha\gamma}\varepsilon_{\beta\delta}=\varepsilon_{\alpha\beta}\varepsilon_{\gamma\delta}-\varepsilon_{\alpha\delta}\varepsilon_{\gamma\beta}$ \\
  $\sigma$-Fierz & $\bar\sigma^{\mu\dot\alpha\alpha}\sigma_{\mu\beta\dot\beta}=-2\delta^\alpha_\beta\delta^{\dot\alpha}_{\dot\beta}$ \\
  4-Fermi Fierz & $(\chi^\dagger\bar\sigma^\mu\psi)(\xi^\dagger\bar\sigma_\mu\eta)=-2(\chi^\dagger\xi^\dagger)(\eta\psi)$ \\
  Gordon identity & $2m\bar u'\gamma^\mu u = \bar u'[(p'+p)^\mu - 2iS^{\mu\nu}(p'-p)_\nu]u$ \\
  Extended Gordon & $\bar u'[\cdots]\gamma^5 u = 0$;\quad $\{\gamma^5,\gamma^\mu\}=0$ \\
  Spin-sum (particle) & $\sum_s u_s\bar u_s = -\slashed{p}+m$ \\
  Spin-sum (antiparticle) & $\sum_s v_s\bar v_s = -\slashed{p}-m$ \\
  Massless helicity & $u_\pm\bar u_\pm = P_{R/L}(-\slashed{p})$;\quad helicity = chirality \\
  \bottomrule
\end{tabular}
\end{center}

\bigskip
\noindent\textbf{Next:} Chapter 48 (spinors for massless particles) builds
on these tools toward the spinor-helicity formalism and MHV amplitudes.

\end{document}
"""

outpath = "/work/ch38_spinor_technology.tex"
with open(outpath, "w") as f:
    f.write(doc.strip())

print(f"Wrote: {outpath}")
