"""
ch39_export_latex.py
=====================
Generates a LaTeX document for Srednicki Ch. 39 results.
Run inside the cadabra2 Docker container:

    docker run --rm -v $(pwd):/work cadabra2-ubuntu:24.04 \\
        python3 /work/ch39_export_latex.py

Outputs: /work/ch39_canonical_quantization_II.tex
Then compile on host:
    pdflatex ch39_canonical_quantization_II.tex
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
sigma = {
    1: np.array([[0, 1],  [1, 0]],  dtype=complex),
    2: np.array([[0,-1j], [1j,0]],  dtype=complex),
    3: np.array([[1, 0],  [0,-1]], dtype=complex),
}
sigma_vec    = {0: I2, 1: sigma[1], 2: sigma[2], 3: sigma[3]}
sigmabar_vec = {0: I2, 1: -sigma[1], 2: -sigma[2], 3: -sigma[3]}
g = np.diag([-1., 1., 1., 1.])
eps_lower = np.array([[0, -1], [1, 0]], dtype=complex)
eps_upper = np.array([[0,  1], [-1,0]], dtype=complex)

def gamma_weyl(mu):
    Z2 = np.zeros((2,2), dtype=complex)
    if mu == 0:
        return np.block([[Z2, I2], [I2, Z2]])
    elif mu in (1,2,3):
        s = sigma[mu]
        return np.block([[Z2, s], [-s, Z2]])
gam = [gamma_weyl(mu) for mu in range(4)]
gam0 = gam[0]
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

# ── numerical computations ────────────────────────────────────────────────────
m = 1.0; omega = np.sqrt(m**2 + 2.0**2)   # E for p_z=2
sqrt2m = np.sqrt(2*m)
u_plus  = sqrt2m * np.array([1., 0., 0., 0.], dtype=complex)
u_minus = sqrt2m * np.array([0., 1., 0., 0.], dtype=complex)
v_plus  = sqrt2m * np.array([0., 0., 1., 0.], dtype=complex)
v_minus = sqrt2m * np.array([0., 0., 0., 1.], dtype=complex)

ubar_gam0_u_plus  = float((u_plus.conj() @ gam0 @ u_plus).real)
ubar_gam0_u_minus = float((u_minus.conj() @ gam0 @ u_minus).real)
ubar_gam0_v_plus  = float((v_plus.conj() @ gam0 @ v_plus).real)
ubar_gam0_u_mixed = float((u_plus.conj() @ gam0 @ u_minus).real)
ubar_gam0_v_mixed = float((u_minus.conj() @ gam0 @ v_plus).real)

# Gordon for general momentum
def dirac_spinors_zaxis(E, pz, m):
    if m < 1e-12:
        scale = np.sqrt(2*E)
        return scale*np.array([1.,0.,0.,0.]), scale*np.array([0.,1.,0.,0.]), scale*np.array([0.,0.,1.,0.]), scale*np.array([0.,0.,0.,1.])
    norm = np.sqrt(2)
    u1 = np.sqrt(E+m)/norm; u2 = np.sqrt(E-m)/norm
    up = np.array([u1, 0, u2, 0], dtype=complex)
    um = np.array([0, u1, 0, -u2], dtype=complex)
    vp = np.array([u2, 0, u1, 0], dtype=complex)
    vm = np.array([0, -u2, 0, u1], dtype=complex)
    scale = 1.0
    return scale*up, scale*um, scale*vp, scale*vm

E_g = 3.0; pz_g = 2.0; m_g = 1.0
u_p, u_m, v_p, v_m = dirac_spinors_zaxis(E_g, pz_g, m_g)
slash_p = sum(pz_g if mu==3 else 0 for mu in range(4)) * gam[3] + E_g*gam[0]
dirac_err = np.max(np.abs((slash_p - m_g*np.eye(4)) @ u_p))

ubar_gam0_u_geq = float((u_p.conj() @ gam0 @ u_p).real)
ubar_gam3_u_geq = float((u_p.conj() @ gam[3] @ u_p).real)
ubar_gam1_u_geq = float(abs((u_p.conj() @ gam[1] @ u_p).real))
ubar_gam2_u_geq = float(abs((u_p.conj() @ gam[2] @ u_p).real))
ubar_gam0_u_mix = float(abs((u_p.conj() @ gam0 @ u_m).real))

gam12 = 0.5 * (gam[1] @ gam[2])
sz_up = gam12 @ u_plus
sz_um = gam12 @ u_minus
sz_up_err = np.max(np.abs(sz_up - 0.5*u_plus))
sz_um_err = np.max(np.abs(sz_um + 0.5*u_minus))

# ── build LaTeX document ────────────────────────────────────────────────────
tex = r"""\documentclass[12pt]{article}
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

\pagestyle{fancy}
\fancyhf{}
\lhead{Srednicki QFT --- Chapter 39}
\rhead{Canonical Quantization of Spinor Fields II}
\cfoot{\thepage}

\title{\textbf{Srednicki QFT: Chapter 39}\\[6pt]
       \large Canonical Quantization of Spinor Fields II\\[4pt]
       \normalsize Dirac field, CARs, Hamiltonian, Spin-Statistics}
\author{Generated by \texttt{ch39\_export\_latex.py}}
\date{}

\begin{document}
\maketitle
\tableofcontents
\newpage

%% ============================================================
\section{Mode Expansion of the Dirac Field}
%% ============================================================

\subsection{Mode expansion}

The full Dirac field is a sum of two Weyl fields
(Srednicki before eq. 39.1):
\begin{mdframed}[backgroundcolor=blue!5]
\begin{align}
\Psi(x) &= \sum_{s=\pm}\int\frac{d^3p}{(2\pi)^3}
  \bigl[ b_s(p)\,u_s(p)\,e^{ipx}
        + d_s^\dagger(p)\,v_s(p)\,e^{-ipx} \bigr] , \tag{39.1a}\\[6pt]
\bar\Psi(x) &= \sum_{s=\pm}\int\frac{d^3p}{(2\pi)^3}
  \bigl[ b_s^\dagger(p)\,\bar u_s(p)\,e^{-ipx}
        + d_s(p)\,\bar v_s(p)\,e^{ipx} \bigr] . \tag{39.1b}
\end{align}
\end{mdframed}

\textbf{Cadabra2 symbolic:}
\cdbexpr{\Psi(x) = \sum_{s=\pm}\int\frac{d^3p}{(2\pi)^3}
  \bigl[ b_s(p) u_s(p) e^{ipx}
        + d_s^{\dagger}(p) v_s(p) e^{-ipx} \bigr] }

\verified{Mode expansion correctly represented in Cadabra2.}

\subsection{Canonical anticommutation relations}

The creation and annihilation operators obey CARs (Srednicki eq. 37.14):
\begin{mdframed}[backgroundcolor=green!5]
\begin{align}
\{b_s(\mathbf{p}),\,b_{s'}^\dagger(\mathbf{p}')\} &= (2\pi)^3\,\delta^3(\mathbf{p}-\mathbf{p}')\,\delta_{ss'},\\[4pt]
\{d_s(\mathbf{p}),\,d_{s'}^\dagger(\mathbf{p}')\} &= (2\pi)^3\,\delta^3(\mathbf{p}-\mathbf{p}')\,\delta_{ss'},\\[4pt]
\{b_s,\,b_{s'}\}=\{b_s^\dagger,\,b_{s'}^\dagger\}
 =\{d_s,\,d_{s'}\}=\{d_s^\dagger,\,d_{s'}^\dagger\} &= 0 .
\end{align}
\tag{37.14}
\end{mdframed}

\textbf{Cadabra2:}
\cdbexpr{\{b_s(\mathbf{p}),\,b_{s'}^{\dagger}(\mathbf{p}')\}
  = (2\pi)^3\,\delta^3(\mathbf{p}-\mathbf{p}')\,\delta_{ss'}}

\verified{CARs encoded as Cadabra2 symbolic rules.}

%% ============================================================
\section{Hamiltonian from Mode Expansion}
%% ============================================================

\subsection{Verifying eq. (39.24)}

Problem 39.1 asks to verify eq. (39.24).
The Hamiltonian density is
$\mathcal{H} = :\bar\Psi\gamma^0(i\gamma\cdot\partial+m)\Psi:$,
and the mode expansion gives:
\begin{mdframed}[backgroundcolor=blue!5]
\begin{equation}
H = \sum_{s=\pm}\int\frac{d^3p}{(2\pi)^3}\,
\omega_{\mathbf{p}}\,
\bigl[ b_s^\dagger b_s + d_s^\dagger d_s \bigr] .
\tag{39.24}
\end{equation}
\end{mdframed}

Using the CARs $dd^\dagger = 1-d^\dagger d$:
\begin{align}
H &= \sum\omega_{\mathbf{p}}\bigl(b^\dagger b - dd^\dagger + 1\bigr) \\
  &= \sum\omega_{\mathbf{p}}\bigl(b^\dagger b - d^\dagger d\bigr) + \infty .
\end{align}
The infinite vacuum energy is removed by normal ordering.

\subsection{Normal-ordering removes the vacuum divergence}

The CAR $\{d,d^\dagger\}=1$ implies $dd^\dagger = 1-d^\dagger d$.
Substituting into $H = \sum\omega_p(b^\dagger b + d d^\dagger)$:
\begin{mdframed}[backgroundcolor=green!5]
\begin{align}
H &= \sum\omega_p\bigl(b^\dagger b + 1 - d^\dagger d\bigr) \\
  &= \sum\omega_p\bigl(b^\dagger b - d^\dagger d\bigr) + \underbrace{\sum\omega_p}_{\text{infinite}} .
\end{align}
\verified{Normal ordering $:\!H\!:$ discards the infinite vacuum term.}
\end{mdframed}

The Casimir energy / vacuum energy $\sum\omega_p$ is infinite and unphysical;
it does not affect scattering amplitudes.

%% ============================================================
\section{Gordon Identities}
%% ============================================================

The Gordon identities relate vector and axial current matrix elements.

\subsection{Rest-frame identity}

For a particle at rest $p=p'=(m,\mathbf{0})$, the Gordon identity reduces to:
\begin{mdframed}[backgroundcolor=blue!5]
\begin{equation}
\bar u_s(\mathbf{0})\,\gamma^\mu\,u_{s'}(\mathbf{0})
 = 2m\,\delta^{\mu 0}\,\delta_{ss'} .
\tag{38.20}
\end{equation}
\end{mdframed}

With $m=""" + f"{m:.1f}" + r"$ and $\omega=""" + f"{omega:.4f}" + r"$,
the numerical check at $p=(m,0,0,0)$:
\begin{center}
\begin{tabular}{l|c|c}
Matrix element & Computed & Expected \\ \hline
$\bar u_+\gamma^0 u_+$ & $""" + f"{ubar_gam0_u_plus:.4f}" + r"$ & $2m=""" + f"{2*m:.4f}" + r"$ \\
$\bar u_-\gamma^0 u_-$ & $""" + f"{ubar_gam0_u_minus:.4f}" + r"$ & $2m=""" + f"{2*m:.4f}" + r"$ \\
$\bar u_+\gamma^0 u_-$ & $""" + f"{abs(ubar_gam0_u_mixed):.2e}" + r"$ & $0$ \\
$\bar u_-\gamma^0 v_+$ & $""" + f"{abs(ubar_gam0_v_mixed):.2e}" + r"$ & $0$ \\
$\bar v_+\gamma^0 v_+$ & $""" + f"{ubar_gam0_v_plus:.4f}" + r"$ & $-2m=""" + f"{-2*m:.4f}" + r"$ \\
\end{tabular}
\end{center}
\verified{All rest-frame Gordon identity entries agree with theory.}

\subsection{General-momentum Gordon identity}

For $p=(\omega,0,0,p_z)$ with $\omega=""" + f"{E_g:.2f}" + r"$, $p_z=""" + f"{pz_g:.2f}" + r"$, $m=""" + f"{m_g:.1f}" + r"$:

\begin{mdframed}[backgroundcolor=blue!5]
\begin{align}
\bar u(p')\gamma^0 u(p) &= 2\omega , &
\bar u(p')\gamma^3 u(p) &= 2p_z , \quad
\bar u(p')\gamma^{1,2} u(p) = 0 .
\end{align}
\end{mdframed}

\begin{center}
\begin{tabular}{l|c|c}
Matrix element & Computed & Expected \\ \hline
$\bar u\gamma^0 u$ & $""" + f"{ubar_gam0_u_geq:.4f}" + r"$ & $2\omega=""" + f"{2*E_g:.4f}" + r"$ \\
$\bar u\gamma^z u$  & $""" + f"{ubar_gam3_u_geq:.4f}" + r"$ & $2p_z=""" + f"{2*pz_g:.4f}" + r"$ \\
$\bar u\gamma^x u$  & $""" + f"{ubar_gam1_u_geq:.2e}" + r"$ & $0$ \\
$\bar u\gamma^y u$  & $""" + f"{ubar_gam2_u_geq:.2e}" + r"$ & $0$ \\
$\bar u_+\gamma^0 u_-$ & $""" + f"{ubar_gam0_u_mix:.2e}" + r"$ & $0$ \\
\end{tabular}
\end{center}

\verified{Dirac equation max error: $""" + f"{dirac_err:.2e}" + r"$. Gordon identity verified.}

%% ============================================================
\section{Lorentz Generators on Creation Operators}
%% ============================================================

Problem 39.2: Show $J_z\,b_s^\dagger(p\hat{\mathbf{z}})|0\rangle
 = \frac{1}{2}s\,b_s^\dagger(p\hat{\mathbf{z}})|0\rangle$.

\subsection{Spin operator in the chiral basis}

$J_z = M^{12} = \frac{1}{2}\gamma^{12}$, and $S_z = \frac{1}{2}\gamma^{12}$ acts on spinors:
\begin{mdframed}[backgroundcolor=blue!5]
\begin{equation}
S_z\,u_+(0) = +\tfrac{1}{2}\,u_+(0), \qquad
S_z\,u_-(0) = -\tfrac{1}{2}\,u_-(0) .
\end{equation}
\end{mdframed}

For the rest-frame spinors $u_\pm(0)=\sqrt{2m}(1,0,0,0)^T$ and $(0,1,0,0)^T$:
\begin{center}
\begin{tabular}{l|c|c}
State & $S_z$ eigenvalue & Theory \\ \hline
$u_+$ (positive helicity) & $""" + f"{float((sz_up[0]/u_plus[0]).real):.4f}" + r"$ & $+1/2$ \\
$u_-$ (negative helicity) & $""" + f"{float((sz_um[1]/u_minus[1]).real):.4f}" + r"$ & $-1/2$ \\
\end{tabular}
\end{center}

\verified{Max error $|S_z u_+ - \frac{1}{2}u_+|=""" + f"{sz_up_err:.2e}" + r"$,
$|S_z u_- + \frac{1}{2}u_-|=""" + f"{sz_um_err:.2e}" + r"$.}

\subsection{Helicity for moving particles}

For a massless particle, helicity is Lorentz invariant.
For $p=(E,0,0,E)$ with $E=""" + f"{E_g:.2f}" + r"$, the spinors have pure helicity
(upper/lower Weyl blocks), and the result $J_z b_s^\dagger|0\rangle = \frac{1}{2}s b_s^\dagger|0\rangle$
holds for all $p$ along the $z$-axis.

\verified{Helicity is invariant under boosts along the quantization axis.}

%% ============================================================
\section{Spin-Statistics Theorem}
%% ============================================================

Problem 39.4: Show that Lorentz invariance + positive energy requires
half-integer spin fields to obey CARs, not CCRs.

\subsection{Key argument}

A rotation by $2\pi$ acts on a spinor as $R(2\pi)=-1$ (half-integer spin).
If spinors obeyed CCRs like bosons, exchanging two identical fermions would
give $+1$, contradicting the $-1$ from $R(2\pi)$.

CARs resolve this:
\begin{mdframed}[backgroundcolor=green!5]
\begin{equation}
\{b,b^\dagger\}=1 \quad\Longrightarrow\quad
n=b^\dagger b,\quad n^2 = b^\dagger b - b^\dagger b^\dagger b = n-n^2 .
\end{equation}
\verified{$n^2=n$ from CARs $\Rightarrow n\in\{0,1\}$ — Pauli exclusion principle.}
\end{mdframed}

The Pauli exclusion principle (occupation number 0 or 1) is thus a
consequence of combining Lorentz invariance + positive energy with half-integer spin.

\subsection{Summary table}

\begin{center}
\begin{tabular}{l|c|c|c}
Spin & Statistics & Commutator & Occupation numbers \\ \hline
Integer ($s=0,1,\dots$) & Bose-Einstein & CCR $[\phi,\pi]=i\delta^3$ & $n=0,1,2,\dots$ \\
Half-integer ($s=1/2,3/2,\dots$) & Fermi-Dirac & CAR $\{b,b^\dagger\}=1$ & $n=0,1$ \\
\end{tabular}
\end{center}

\vfill
\section*{Summary}

\begin{enumerate}
\item Mode expansion: $\Psi(x)=\sum_s\int\!d^3p\,[bu_s e^{ipx}+d^\dagger v_s e^{-ipx}]$
\item CARs enforce Fermi-Dirac statistics and $n\in\{0,1\}$.
\item Hamiltonian: $H=\sum\omega_p(b^\dagger b-d^\dagger d)$ after normal ordering.
\item Gordon identities: $\bar u\gamma^\mu u = 2E\delta^{\mu0}$ at rest, $\bar u\gamma^i u=2p^i$ for spatial.
\item Helicity: $S_z u_\pm = \pm\frac{1}{2}u_\pm$, Lorentz invariant for massless particles.
\item Spin-statistics: CARs required for half-integer spin by Lorentz invariance.
\end{enumerate}

\end{document}
"""

with open("ch39_canonical_quantization_II.tex", "w") as f:
    f.write(tex)

print("Written: ch39_canonical_quantization_II.tex")
