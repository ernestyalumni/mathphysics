r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 40: Spin Sums}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Chapter 40 connects Ch.39's Feynman rules to observable cross sections.}
Every experimental prediction requires summing over final spins and averaging over
initial spins.  This chapter shows how to do that efficiently using trace theorems,
and introduces the crucial distinction between the \emph{trace method} (what we do
here) and the \emph{helicity amplitude method} (Ch.42, the heart of MHV).
\end{mdframed}

\section{Spin-Averaging Factors}\label{sec:averaging}

For an unpolarized initial state, we must average over the spins of each
incoming particle.  For a:
\begin{itemize}
\item Fermion (spin-$\tfrac12$): average over 2 spin states $\to \tfrac12$
\item Photon/gluon (spin-1): average over 2 physical polarizations (transverse)
      $\to \tfrac12$ in QED/QCD gauge;
      in general covariant gauge $\to \tfrac13$ of the 4 gauge-field polarizations
\end{itemize}

For $e^+e^- \to \mu^+\mu^-$ at high energy:
\begin{equation}
\frac{d\sigma}{d\Omega}\bigg|_{
\substack{\text{unpolarized}\\\text{partons}}}
= \frac{\alpha^2}{4s}\left[
1 + \cos^2\theta
+ \frac{s-4m_\mu^2}{s}\sin^2\theta\right] .
\tag{40.$\star$}\label{eq:unpolarized_cross}
\end{equation}
This is a genuine prediction: no free parameters after fixing $\alpha$ and $m_\mu$.

\section{Trace Method}\label{sec:trace_method}

The trace method converts a squared amplitude with closed fermion loops into a
trace of $\gamma$ matrices:
\begin{mdframed}[backgroundcolor=green!5]
\begin{equation}
\frac{1}{4}\sum_{\mathrm{spins}}|M|^2
= \frac{1}{4}\,\Tr\bigl[\ldots\gamma^\mu\ldots\gamma^\nu\bigr]\,
\Tr\bigl[\ldots\gamma_\mu\ldots\gamma_\nu\bigr] .
\tag{40.1}\label{eq:trace_method}
\end{equation}
\end{mdframed}
Each closed fermion loop contributes one trace; external fermion lines contribute
to the spin-sum as in Ch.39.

\section{Trace Theorems}\label{sec:trace_theorems}

From Ch.39 (verified by numpy), the complete set:
\begin{align}
\Tr[\mathbf{1}] &= 4, \\
\Tr[\gamma^\mu] &= 0, \\
\Tr[\gamma^\mu\gamma^\nu] &= 4\eta^{\mu\nu}, \\
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho] &= 0, \\
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4\bigl(\eta^{\mu\nu}\eta^{\rho\sigma}
          -\eta^{\mu\rho}\eta^{\nu\sigma}
          +\eta^{\mu\sigma}\eta^{\nu\rho}\bigr) .
\end{align}

A key identity follows from contracting with $p_{1\mu}p_{2\nu}$:
\begin{equation}
\Tr[(\p1\!\!\!/+m)\gamma^\mu(\p2\!\!\!/+m)\gamma^\nu]
= 2\Bigl[p_1^\mu p_2^\nu + p_1^\nu p_2^\mu
- \eta^{\mu\nu}(p_1\cdot p_2 - m^2)\Bigr] .
\tag{40.2}\label{eq:contracted_trace}
\end{equation}
This is the workhorse of QED cross-section calculations.

\section{Helicity Amplitudes vs. Trace Method}\label{sec:vs_helicity}

The trace method and the helicity method (Ch.42) give the same answers but
are computationally complementary:

\begin{center}
\begin{tabular}{|l|p{0.4\textwidth}|p{0.4\textwidth}|}
\hline
Feature & Trace Method & Helicity Method \\
\hline
Strength & Systematic; works for any process & Extremely compact for MHV \\
Best for & Loops, many fermions & Gluon amplitudes, gauge theory \\
Output & Scalar products $p_i\cdot p_j$ & Spinor brackets $\langle ij\rangle$, $[ij]$ \\
Symmetry & Manifest Lorentz invariance & Manifest helicity structure \\
\hline
\end{tabular}
\end{center}

For MHV amplitudes (Ch.60), the helicity method is overwhelmingly superior:
the Parke-Taylor formula $A_n \propto \langle 12\rangle^4/\prod\langle i\,i+1\rangle$
is a \emph{single algebraic expression} whereas the trace method would give
$O(n!)$ terms.

\section{Compton Scattering: dsigma/dOmega}\label{sec:compton}

For $e^-(p) + \gamma(k,\varepsilon) \to e^-(p') + \gamma(k',\varepsilon')$,
the unpolarized Klein-Nishina cross-section is:

\begin{mdframed}[backgroundcolor=red!5]
\begin{equation}
\frac{d\sigma}{d\Omega}
= \frac{\alpha^2}{2m_e^2}
\left(\frac{\omega'}{\omega}\right)^2
\left[\frac{\omega'}{\omega} + \frac{\omega}{\omega'} - \sin^2\theta\cos^2\phi\right] .
\tag{40.$\star\star$}\label{eq:klein_nishina}
\end{equation}
\end{mdframed}

In the high-energy limit $\omega\gg m_e$, $\omega'/\omega \approx 1$ and
this reduces to the Thomson cross-section $\sigma_T = \frac{8\pi r_e^2}{3}$
integrated over angles.

\section{Significance for the MHV Program}\label{sec:mhv}

Ch.40's trace method is how you compute $|M|^2$ for processes with closed
fermion loops (which appear in NLO QCD calculations).  However, for the
tree-level MHV program (Ch.42, 48, 60, BCFW), the helicity method replaces
the trace theorems with Schouten identities.  Understanding both gives you
the full picture of how amplitudes are built.

\end{document}
"""
import numpy as np

OUT = __import__('textwrap').dedent(OUT) if 'OUT' in dir() else ""

# Simple write
with open("/work/ch40_spin_sums.tex", "w") as f:
    f.write(OUT if OUT else r"""\documentclass[12pt]{article}
\usepackage[margin=1in]{geometry}\usepackage{amsmath,amssymb,xcolor,mhchem,mdframed}
\begin{document}\title{Chapter 40: Spin Sums}\date{\today}\maketitle
\section{Spin-Averaging}
For unpolarized scattering, average over initial spins: $\tfrac12$ per fermion.
\section{Trace Method}
Squared amplitudes with closed fermion loops:
$\frac{1}{4}\sum|M|^2 = \frac{1}{4}\Tr[\cdots\gamma^\mu\cdots\gamma_\mu]$.
\section{Trace Theorems}
$\Tr[\gamma^\mu\gamma^\nu]=4\eta^{\mu\nu}$;
$\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
=4(\eta^{\mu\nu}\eta^{\rho\sigma}-\eta^{\mu\rho}\eta^{\nu\sigma}+\eta^{\mu\sigma}\eta^{\nu\rho})$.
\section{Compton: Klein-Nishina}
$d\sigma/d\Omega = (\alpha^2/2m_e^2)(\omega'/\omega)^2
[(\omega'/\omega)+(\omega/\omega')-\sin^2\theta\cos^2\phi]$.
\section{Helicity vs. Trace}
Trace method: Lorentz-invariant but $O(n!)$ terms.
Helicity method (Ch.42): compact $\langle ij\rangle^4/\prod\langle i\,i+1\rangle$ expressions.
\vfill\bibliographystyle{plain}\end{document}""")

# We'll write the full version directly
full_out = r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 40: Spin Sums}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Chapter 40 bridges the gap between Feynman rules (Ch.39) and
observable cross sections.} Every experimental prediction requires
summing over final spins and averaging over initial spins. This chapter
shows how to do that efficiently using trace theorems.
\end{mdframed}

\section{Spin-Averaging Factors}\label{sec:averaging}
For an unpolarized initial state:
\begin{itemize}
\item Fermion (spin-$\frac12$): average over 2 spin states $\to \frac12$
\item Photon/gluon: average over 2 physical polarizations $\to \frac12$ (QED/QCD gauge)
\end{itemize}
For $e^+e^-\to\mu^+\mu^-$ at high energy ($s\gg m_\mu^2$):
\begin{equation}
\frac{d\sigma}{d\Omega}
= \frac{\alpha^2}{4s}\left[1+\cos^2\theta\right] .
\tag{40.$\star$}
\end{equation}

\section{Trace Method}\label{sec:trace_method}
Squared amplitudes with closed fermion loops become traces:
\begin{mdframed}[backgroundcolor=green!5]
\begin{equation}
\frac{1}{4}\sum_{\rm spins}|M|^2
= \frac{1}{4}\,\Tr\bigl[(\p1\!\!\!/+m)\gamma^\mu(\p2\!\!\!/+m)\gamma^\nu\bigr]\;
\Tr\bigl[(\p3\!\!\!/+m')\gamma_\mu(\p4\!\!\!/+m')\gamma_\nu\bigr] .
\tag{40.1}
\end{equation}
\end{mdframed}
Each closed fermion loop contributes one trace.

\section{Trace Theorems}\label{sec:trace_theorems}
Fundamental identities ($\eta=\mathrm{diag}(-1,+1,+1,+1)$):
\begin{align}
\Tr[\mathbf{1}] &= 4, &
\Tr[\gamma^\mu\gamma^\nu] &= 4\eta^{\mu\nu}, \\[2mm]
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4\bigl(\eta^{\mu\nu}\eta^{\rho\sigma}
          -\eta^{\mu\rho}\eta^{\nu\sigma}
          +\eta^{\mu\sigma}\eta^{\nu\rho}\bigr) .
\end{align}

Contracted trace identity (used in every QED cross-section):
\begin{equation}
\Tr[(\p1\!\!\!/+m)\gamma^\mu(\p2\!\!\!/+m)\gamma^\nu]
= 2\Bigl[p_1^\mu p_2^\nu + p_1^\nu p_2^\mu
         -\eta^{\mu\nu}(p_1\cdot p_2-m^2)\Bigr] .
\tag{40.2}
\end{equation}

\section{Helicity vs. Trace Method}\label{sec:vs_helicity}
\begin{center}
\begin{tabular}{|l|p{0.35\linewidth}|p{0.35\linewidth}|}
\hline
Feature & Trace Method & Helicity Method (Ch.42) \\
\hline
Output & Scalar products $p_i\cdot p_j$ & $\langle ij\rangle$, $[ij]$ \\
Best for & Loops, fermions & Gluon tree amplitudes \\
Result size & $O(n!)$ terms & Single expression for MHV \\
\hline
\end{tabular}
\end{center}

\section{Klein-Nishina (Compton Scattering)}\label{sec:compton}
\begin{mdframed}[backgroundcolor=red!5]
\begin{equation}
\frac{d\sigma}{d\Omega}
= \frac{\alpha^2}{2m_e^2}
\left(\frac{\omega'}{\omega}\right)^2
\left[\frac{\omega'}{\omega}+\frac{\omega}{\omega'}
-\sin^2\theta\cos^2\phi\right] .
\tag{40.$\star\star$}
\end{equation}
\end{mdframed}
In the Thomson limit $\omega\gg m_e$: $\sigma\to\sigma_T=\frac{8\pi r_e^2}{3}$.

\section{Significance for MHV Research}
Ch.40 is the last chapter before Ch.42 that uses only Lorentz-covariant
methods.  The trace theorems here are replaced by Schouten identities in the
helicity formalism: the same algebraic structure, but packaged in spinor
brackets that make helicity selection rules transparent.

\end{document}
"""

with open("/work/ch40_spin_sums.tex", "w") as f:
    f.write(full_out)
print(f"Wrote /work/ch40_spin_sums.tex ({len(full_out)} chars)")
