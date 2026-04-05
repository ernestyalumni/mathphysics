r"""\documentclass[12pt]{article}
\usepackage{amsmath,amssymb,amsthm,fullpage,hyperref,xcolor,mdframed}
\usepackage[a4paper,margin=1in]{geometry}

% ── Custom commands ────────────────────────────────────────────────────────────
\DeclareDocumentCommand\bra{ m }{\langle{#1}\rvert}
\DeclareDocumentCommand\ket{ m }{\lvert{#1}\rangle}
\DeclareDocumentCommand\avg{ m }{\langle{#1}\rangle}

% =============================================================================
\title{Chapter 39: The Feynman Rules for Dirac Fields}
\subtitle{Srednicki QFT --- Cadabra2 Computation and Notes}
\author{Srednicki QFT Notes Series}
\date{\today}

\begin{document}
\maketitle

\begin{center}
\textbf{Context:} Ch.39 applies the canonical quantization of Ch.37 to Dirac
fields, deriving the \emph{Feynman propagator} $S_F(x-y)$ and the complete set
of Feynman rules for spinor QED.  These are the rules you will use in every
tree-level and loop calculation involving fermions.
\end{center}

\tableofcontents
\newpage

% ─────────────────────────────────────────────────────────────────────────────
\section{Dirac Propagator}\label{sec:propagator}
% ─────────────────────────────────────────────────────────────────────────────

The Dirac field $\psi(x)$ was promoted to a quantum operator in Ch.37 via
canonical anticommutation relations.  The \emph{Feynman propagator} is the
time-ordered vacuum expectation value:
\begin{equation}
S_F(x-y) \equiv \langle 0| T\{\,\psi(x)\bar\psi(y)\,]|0\rangle .
\tag{39.1}\label{eq:sf_def}
\end{equation}

For a free Dirac field the mode expansion is
\begin{equation}
\psi(x) = \int \frac{d^3k}{(2\pi)^3}
\frac{1}{\sqrt{2E_\mathbf{k}}}
\sum_s \bigl[ u_s(\mathbf{k})\, a_\mathbf{k}^s\,
e^{-ik\cdot x} + v_s(\mathbf{k})\, b_{\mathbf{k}}^{s\dagger}\,
e^{+ik\cdot x}\bigr] .
\tag{39.2}\label{eq:dirac_mode}
\end{equation}
where $a_\mathbf{k}^s$, $b_{\mathbf{k}}^{s\dagger}$ are the fermion and
antifermion annihilation/creation operators.

Evaluating \eqref{eq:sf_def} using \eqref{eq:dirac_mode} and the
anticommutation algebra gives the famous result:
\begin{mdframed}[backgroundcolor=blue!5]
\begin{equation}
S_F(x-y) = \int\frac{d^4k}{(2\pi)^4}\,
\frac{i}{k\!\!\!/-m+i\varepsilon}\,
e^{-ik\cdot(x-y)}
= \frac{i}{k\!\!\!/-m+i\varepsilon}\;(2\pi)^4\delta^{(4)}(k_\mu-\text{变量})\dots
\tag{39.3}\label{eq:sf_momentum}
\end{equation}
In position space:
\begin{equation}
S_F(x-y) = \left(i\gamma^\mu\partial_\mu + m\right)
\frac{1}{4\pi^2(x-y)^2}
\quad\text{(in 4D, massless case)}.
\end{equation}
\end{mdframed}

The minus-$i\varepsilon$ prescription is the Feynman $i\varepsilon$ that tells
you how to handle the poles at $k^0 = \pm E_\mathbf{k}$: you close the
 contour in the lower half-plane for $t > t'$ (positive-energy propagation)
and in the upper half-plane for $t < t'$ (negative-energy = antiparticle).

% ─────────────────────────────────────────────────────────────────────────────
\section{External Fermion Wave Functions}\label{sec:external}
% ─────────────────────────────────────────────────────────────────────────────

When a fermion line appears as an \emph{external} leg in a Feynman diagram,
it is replaced by a \emph{spinor wave function} $u_s(p)$ or $v_s(p)$.

\subsection*{Incoming fermion $u_s(p)$}
\begin{equation}
u_s(p) = \begin{pmatrix}
\sqrt{E_\mathbf{p}+\mathbf{p}\cdot\boldsymbol{\sigma}} \,\chi_s \\[2mm]
\sqrt{E_\mathbf{p}-\mathbf{p}\cdot\boldsymbol{\sigma}} \,\chi_s
\end{pmatrix},
\qquad
\chi_+ = \begin{pmatrix}1\\0\end{pmatrix},\;
\chi_- = \begin{pmatrix}0\\1\end{pmatrix} .
\tag{39.5}\label{eq:u_spinor}
\end{equation}
In the massless limit $E=|\mathbf{p}|$ and for the special case
$\mathbf{p} = (0,0,E)$ (momentum along the $+z$ axis):
\begin{equation}
u_+(k) = \sqrt{2E}\,\begin{pmatrix}1\\0\\1\\0\end{pmatrix},
\qquad
u_-(k) = \sqrt{2E}\,\begin{pmatrix}0\\1\\0\\-1\end{pmatrix}.
\tag{39.6}\label{eq:u_massless}
\end{equation}
The subscripts $\pm$ label \emph{helicity}, not spin along $z$: $h = +\tfrac12$
for $u_+$ and $h=-\tfrac12$ for $u_-$.

\subsection*{Outgoing antifermion $\bar v_s(p)$}
The antifermion wave function is
\begin{equation}
\bar v_s(p) = \bar u_s(p)\,\gamma^0
\quad\text{(Dirac adjoint)} .
\tag{39.7}\label{eq:v_bar}
\end{equation}
In the massless helicity basis: $\bar v_+(k) = \sqrt{2E}(0,-1,0,1)$.

% ─────────────────────────────────────────────────────────────────────────────
\section{Spin Sums and the Optical Theorem}\label{sec:spin_sums}
% ─────────────────────────────────────────────────────────────────────────────

When computing $|M|^2$ for processes with external fermions, we must sum over
the spins of the intermediate fermion lines.  The two essential identities are:

\begin{mdframed}[backgroundcolor=green!5]
\begin{align}
\sum_{s} u_s(p)\,\bar u_s(p) &= p\!\!\!/ + m , \tag{39.8}\label{eq:spin_sum_u}\\
\sum_{s} v_s(p)\,\bar v_s(p) &= p\!\!\!/ - m . \tag{39.9}\label{eq:spin_sum_v}
\end{align}
\end{mdframed}

These are completeness relations: the sum over all spin states of the outer
product of a spinor with its Dirac adjoint gives the identity operator in the
spinor space (the $4\times4$ matrix $p\!\!\!/ + m$ projects onto the
positive-energy subspace).

For an \emph{unpolarrized} process we average over the spins of each incoming
fermion and sum over the spins of each outgoing fermion.  Each external
fermion therefore contributes a factor of $\tfrac12$ when averaging over its
initial spin.

% ─────────────────────────────────────────────────────────────────────────────
\section{Trace Theorems}\label{sec:traces}
% ─────────────────────────────────────────────────────────────────────────────

Closed fermion loops in perturbation theory generate traces of products of
$\gamma$ matrices.  The complete set of trace identities (in 4D Minkowski
space with $\eta_{\mu\nu} = \mathrm{diag}(-1,+1,+1,+1)$) is:

\begin{mdframed}[backgroundcolor=yellow!5]
\begin{align}
\mathrm{Tr}[\mathbf{1}] &= 4, \tag{39.10}\label{eq:trace1} \\
\mathrm{Tr}[\gamma^\mu] &= 0, \tag{39.11}\label{eq:trace_gamma} \\
\mathrm{Tr}[\gamma^\mu\gamma^\nu] &= 4\,\eta^{\mu\nu}, \tag{39.12}\label{eq:trace_2} \\
\mathrm{Tr}[\gamma^\mu\gamma^\nu\gamma^\rho] &= 0, \tag{39.13}\label{eq:trace_3} \\
\mathrm{Tr}[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4\bigl(\eta^{\mu\nu}\eta^{\rho\sigma}
- \eta^{\mu\rho}\eta^{\nu\sigma}
+ \eta^{\mu\sigma}\eta^{\nu\rho}\bigr), \tag{39.14}\label{eq:trace_4}\\
\mathrm{Tr}[\gamma^5] &= 0, \tag{39.15}\label{eq:trace_5} \\
\mathrm{Tr}[\gamma^5\gamma^\mu\gamma^\nu] &= 0, \tag{39.16}\label{eq:trace_5_2} \\
\mathrm{Tr}[\gamma^5\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= -4i\,\varepsilon^{\mu\nu\rho\sigma}. \tag{39.17}\label{eq:trace_5_4}
\end{align}
\end{mdframed}

The last identity, \eqref{eq:trace_5_4}, is the source of the Adler-Bell-Jackiw
(chiral) anomaly in QED.

% ─────────────────────────────────────────────────────────────────────────────
\section{QED Feynman Rules}\label{sec:qed_rules}
% ─────────────────────────────────────────────────────────────────────────────

Combining the propagator and vertex from the QED Lagrangian
$\mathcal{L} = \bar\psi\,(i\gamma^\mu D_\mu - m)\,\psi - \tfrac14 F_{\mu\nu}F^{\mu\nu}$
gives the complete set of Feynman rules for spinor QED:

\begin{mdframed}[backgroundcolor=red!5]
\begin{center}
\begin{tabular}{|l|c|}
\hline
\textbf{Object} & \textbf{Feynman rule} \\
\hline
Photon propagator & $\displaystyle\frac{-i\eta_{\mu\nu}}{q^2+i\varepsilon}$ \\[2mm]
Fermion propagator & $\displaystyle\frac{i(p\!\!\!/+m)}{p^2-m^2+i\varepsilon}$ \\[2mm]
QED vertex & $-ie\,\gamma_\mu$ \\[2mm]
External fermion & $u_s(p)$ or $\bar u_s(p)$ \\[2mm]
External antifermion & $v_s(p)$ or $\bar v_s(p)$ \\[2mm]
External photon & $\varepsilon_\mu(k,\lambda)$ \\
\hline
\end{tabular}
\end{center}
\end{mdframed}

Each vertex contributes a factor $-ie\gamma_\mu$ and conserves 4-momentum at the
vertex.  The photon propagator in Feynman gauge has the simple form $-i\eta_{\mu\nu}/q^2$.
Note: there is \emph{no} photon momentum factor at the vertex (unlike nonabelian
gauge theory, where there are three-gluon vertices with momentum-dependent tensors).

% ─────────────────────────────────────────────────────────────────────────────
\section{Compton Scattering as a Worked Example}\label{sec:compton}
% ─────────────────────────────────────────────────────────────────────────────

Consider Compton scattering $e^-(p) + \gamma(k,\varepsilon) \to e^-(p') + \gamma(k',\varepsilon')$.
The amplitude from QED is (two Feynman diagrams: $s$-channel and $u$-channel):

\begin{align}
i\mathcal{M}
&= (-ie\bar u(p')\gamma^\nu u(p))\,
\frac{-i\eta_{\nu\rho}}{s-m_e^2}\,
(-ie\bar\varepsilon'^\mu\varepsilon_\mu)\,
+\;(p,k,\varepsilon,\varepsilon' \leftrightarrow p',k',-\varepsilon',-\varepsilon) .
\end{align}

After spin-averaging and applying the trace theorems, the unpolarized squared
amplitude simplifies to the Klein-Nishina formula:

\begin{equation}
\frac{d\sigma}{d\Omega}
= \frac{\alpha^2}{2m_e^2}
\left(\frac{\omega'}{\omega}\right)^2
\left[\frac{\omega'}{\omega} + \frac{\omega}{\omega'} - \sin^2\theta\cos^2\phi\right] .
\tag{39.$\star$}\label{eq:klein_nishina}
\end{equation}

This formula is a direct prediction of QED and was historically important as
a test of the Dirac equation's prediction of the electron's magnetic moment.

% ─────────────────────────────────────────────────────────────────────────────
\section{Significance for the MHV Program}\label{sec:mhv_significance}
% ─────────────────────────────────────────────────────────────────────────────

Ch.39 is foundational for the MHV program in two distinct ways:

\subsection*{Fermion lines in gauge theories}
In QCD (and in the electroweak sector), quarks are Dirac fermions.  Their
propagators and spin sums appear in every process with internal or external
quark lines.  The spin-sum formula $\sum u\bar u = p\!\!\!/ + m$ is used
constantly when squaring amplitudes.

\subsection*{Helicity spinors from massless limit}
In the massless limit (relevant for high-energy scattering, the regime of
interest for MHV amplitudes), the spinor wave functions \eqref{eq:u_massless}
become the helicity spinors $\lambda$, $\tilde\lambda$ of Ch.42.  The
$\gamma$-matrix trace calculations of Ch.39 reduce, in the massless chiral
basis, to the angle/square-bracket algebra of the spinor-helicity formalism.

The trace theorems \eqref{eq:trace_1}--\eqref{eq:trace_4} have direct analogues
in the spinor-helicity formalism: they are replaced by Schouten identities
($\langle ij\rangle\langle kl\rangle + \langle jk\rangle\langle il\rangle
+ \langle ki\rangle\langle jl\rangle = 0$) which encode the same algebraic
structure in a much more compact notation.

\vfill
\bibliographystyle{plain}
\bibliography{refs}

\end{document}
"""
# %s placeholder for dynamic content - we'll write a simpler version without
# the heavy latex document string since we compute in python

import cadabra2
from cadabra2 import Ex, __cdbkernel__
import numpy as np

__cdbkernel__ = cadabra2.create_scope()

# ── Setup ─────────────────────────────────────────────────────────────────────
cadabra2.Indices(Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
cadabra2.AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

# ── Dirac gamma matrices (Weyl/chiral basis) ───────────────────────────────────
gamma = {
    0: np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=complex),
    1: np.array([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]], dtype=complex),
    2: np.array([[0, 0, 0, -1j], [0, 0, 1j, 0], [0, 1j, 0, 0], [-1j, 0, 0, 0]], dtype=complex),
    3: np.array([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]], dtype=complex),
}
gamma5 = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=complex)
eta = np.diag([-1, 1, 1, 1])

# ── Trace calculations ────────────────────────────────────────────────────────
tr1 = np.trace(np.eye(4, dtype=complex))
tr_gamma = {mu: np.trace(gamma[mu]) for mu in range(4)}
tr_gamma_gamma = {mu: np.trace(gamma[mu] @ gamma[nu]) for mu in range(4) for nu in range(4)}

# Full trace of 4 gamma matrices -- numerical
tr4 = np.zeros((4,4,4,4), dtype=complex)
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                tr4[mu,nu,rho,sig] = np.trace(gamma[mu] @ gamma[nu] @ gamma[rho] @ gamma[sig])

# Verify trace identities
trace_theorem4 = {
    "tr(gamma^mu gamma^nu) = 4 eta^mu_nu": [],
}
# Check: tr(gm gn) = 4 eta_mn
tr_gmgn_check = {}
for mu in range(4):
    for nu in range(4):
        val = np.trace(gamma[mu] @ gamma[nu]).real
        expected = 4 * eta[mu, nu]
        tr_gmgn_check[(mu, nu)] = round(val, 6) == round(expected, 6)

# Check trace of 4 gammas: tr(gm gn gr gs) = 4(d_mn drs - dm_r dn_s + dm_s dn_r)
tr4_check = {}
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            for sig in range(4):
                val = np.trace(gamma[mu] @ gamma[nu] @ gamma[rho] @ gamma[sig])
                expected = 4 * (
                    eta[mu,nu]*eta[rho,sig]
                    - eta[mu,rho]*eta[nu,sig]
                    + eta[mu,sig]*eta[nu,rho]
                )
                tr4_check[(mu,nu,rho,sig)] = abs(val - expected) < 1e-10

all_tr4_ok = all(tr4_check.values())
all_tr2_ok = all(tr_gmgn_check.values())

# ── Dirac spinors (massless, along +z) ─────────────────────────────────────────
# For massless particle along +z: p = (E, 0, 0, E)
# u_+(p) = sqrt(2E) * (1, 0, 1, 0)^T  (positive helicity)
# u_-(p) = sqrt(2E) * (0, 1, 0, -1)^T  (negative helicity)
E0 = 1.0
sqrt2E = np.sqrt(2 * E0)
u_plus = np.array([sqrt2E, 0, sqrt2E, 0], dtype=complex)
u_minus = np.array([0, sqrt2E, 0, -sqrt2E], dtype=complex)

# Dirac adjoint
ubar_plus = u_plus.conj() @ gamma[0]
ubar_minus = u_minus.conj() @ gamma[0]

# Spin sum: sum_s u_s(p) * ubar_s(p) = p-slash + m (for m=0: p-slash)
p4 = np.array([E0, 0, 0, E0], dtype=complex)
p_slash = sum(p4[mu] * gamma[mu] for mu in range(4))

spin_sum_plus = np.outer(u_plus, ubar_plus)
spin_sum_minus = np.outer(u_minus, ubar_minus)
spin_sum_total = spin_sum_plus + spin_sum_minus

# For massless: should equal p-slash
spin_sum_diff = spin_sum_total - p_slash

# ── Compton scattering kinematics ───────────────────────────────────────────────
# e(p) + gamma(k) -> e(p') + gamma(k')
# s = (p+k)^2 = 2p·k (in massless limit for both)
alpha_em = 1/137.0
m_e = 0.511e-3  # MeV (for scale)
# For high-energy limit, cross-section ratio omega'/omega

# ── Write LaTeX ────────────────────────────────────────────────────────────────
OUT = r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 39: Feynman Rules for Dirac Fields}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Note:} This chapter applies Ch.37's canonical quantization to Dirac
fields. The Feynman propagator and spin sums here are used in every fermion-loop
and external-leg calculation throughout the rest of the course.
\end{mdframed}

\section{Dirac Propagator}\label{sec:propagator}

The Feynman propagator is the time-ordered vev:
\begin{equation}
S_F(x-y) = \langle0| T\{\psi(x)\bar\psi(y)\}|0\rangle .
\end{equation}
Evaluating via the mode expansion gives, in momentum space:
\begin{mdframed}[backgroundcolor=green!5]
\begin{equation}
S_F(k) = \frac{i(k\!\!\!/+m)}{k^2-m^2+i\varepsilon}
      = \frac{i}{k\!\!\!/-m+i\varepsilon} .
\tag{39.3}
\end{equation}
\end{mdframed}

The Feynman $i\varepsilon$ prescription places poles at $k^0 = \pm E_\mathbf{k} \mp i\varepsilon$,
specifying which contour deformation gives the correct time ordering.

\section{External Fermion Wave Functions}\label{sec:external}

For a fermion with momentum $p^\mu=(E,\mathbf{p})$ and spin label $s=\pm$:
\begin{equation}
u_s(p) = \begin{pmatrix}
\sqrt{E+\mathbf{p}\cdot\boldsymbol{\sigma}}\,\chi_s \\
\sqrt{E-\mathbf{p}\cdot\boldsymbol{\sigma}}\,\chi_s
\end{pmatrix},
\qquad
\chi_+ = \begin{pmatrix}1\\0\end{pmatrix},\;
\chi_- = \begin{pmatrix}0\\1\end{pmatrix} .
\tag{39.5}
\end{equation}
In the massless limit and for $\mathbf{p}=(0,0,E)$ (along $+z$):
\begin{equation}
u_+(p) = \sqrt{2E}\begin{pmatrix}1\\0\\1\\0\end{pmatrix},
\qquad
u_-(p) = \sqrt{2E}\begin{pmatrix}0\\1\\0\\-1\end{pmatrix} .
\tag{39.6}
\end{equation}

\section{Spin Sums}\label{sec:spin_sums}

The completeness (spin-sum) identities are:
\begin{mdframed}[backgroundcolor=yellow!5]
\begin{align}
\sum_s u_s(p)\,\bar u_s(p) &= p\!\!\!/+m , \tag{39.8}\\
\sum_s v_s(p)\,\bar v_s(p) &= p\!\!\!/-m . \tag{39.9}
\end{align}
\end{mdframed}

For unpolarized cross sections, each initial fermion spin contributes a factor $\frac12$.

\section{Trace Theorems}\label{sec:traces}

The fundamental trace identities (4D Minkowski, $\eta=\mathrm{diag}(-1,+1,+1,+1)$):
\begin{mdframed}[backgroundcolor=red!5]
\begin{align}
\Tr[\mathbf{1}] &= 4 , \tag{39.10}\\
\Tr[\gamma^\mu\gamma^\nu] &= 4\eta^{\mu\nu} , \tag{39.12}\\
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4\bigl(\eta^{\mu\nu}\eta^{\rho\sigma}
          -\eta^{\mu\rho}\eta^{\nu\sigma}
          +\eta^{\mu\sigma}\eta^{\nu\rho}\bigr) . \tag{39.14}
\end{align}
\end{mdframed}

\vfill
\ centraunted
\textbf{Numpy verification of trace theorems:}
\end{center}
\begin{center}
\begin{tabular}{|l|c|}
\hline
Check & Result \\
\hline
$\Tr[\gamma^\mu\gamma^\nu] = 4\eta^{\mu\nu}$ & """ + ("PASS" if all_tr2_ok else "FAIL") + r""" \\
$\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]$ identity & """ + ("PASS" if all_tr4_ok else "FAIL") + r""" \\
$\Tr[\gamma^5\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma] = -4i\varepsilon^{\mu\nu\rho\sigma}$ & See Ch.41 \\
\hline
\end{tabular}
\end{center}

\section{QED Feynman Rules}\label{sec:qed_rules}

\begin{center}
\begin{tabular}{|l|c|}
\hline
Object & Rule \\
\hline
Photon propagator & $\displaystyle\frac{-i\eta_{\mu\nu}}{q^2+i\varepsilon}$ \\
Fermion propagator & $\displaystyle\frac{i(p\!\!\!/+m)}{p^2-m^2+i\varepsilon}$ \\
QED vertex & $-ie\,\gamma_\mu$ \\
External fermion & $u_s(p)$, $\bar u_s(p)$ \\
External antifermion & $v_s(p)$, $\bar v_s(p)$ \\
External photon & $\varepsilon_\mu(k,\lambda)$ \\
\hline
\end{tabular}
\end{center}

\section{Spinor Completeness (Massless)}\label{sec:spin_completeness}

For massless $p^\mu=(E,0,0,E)$ along $+z$, the spin-sum gives:
\begin{mdframed}
\begin{align*}
\sum_s u_s\bar u_s
&= p\!\!\!/ = \begin{pmatrix}
0&0&E&0\\0&0&0&-E\\E&0&0&0\\0&-E&0&0
\end{pmatrix} .
\end{align*}
Numpy difference from $p\!\!\!/$: $\|$diff$\|_F = $""" + f"{np.linalg.norm(spin_sum_diff):.2e}" + r""" \\
(Zero = PASS)
\end{mdframed}

\section{Significance for the MHV Program}

Ch.39 is foundational because:
\begin{enumerate}
\item Every QCD amplitude with external quarks uses the spin-sum formula
  $\sum u\bar u = p\!\!\!/+m$.
\item The trace theorems reduce loop calculations to scalar products of momenta.
\item In the massless chiral limit, the Dirac spinor basis becomes the
  helicity-spinor basis of Ch.42: $u_+(k)\to\lambda$, $u_-(k)\to\tilde\lambda$.
\end{enumerate}

\vfill\bibliographystyle{plain}\bibliography{refs}
\end{document}
"""

outpath = "/work/ch41_feynman_rules_dirac.tex"
with open(outpath, "w") as f:
    f.write(OUT)
print(f"Wrote {outpath}  ({len(OUT)} chars)")
