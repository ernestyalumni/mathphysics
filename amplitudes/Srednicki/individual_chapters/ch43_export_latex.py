r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 43: Majorana Fields}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Ch.43 covers a special type of spinor field: the Majorana field,}
where the particle is its own antiparticle.  This appears directly in SUSY
(gauginos, higgsinos, squarks, sleptons) and indirectly in the neutralino
dark-matter sector.
\end{mdframed}

\section{The Majorana Condition}\label{sec:condition}

A Dirac field has four degrees of freedom: two helicities for the particle
and two for the antiparticle.  A Majorana field satisfies the constraint:
\begin{mdframed}
\begin{equation}
\psi(x) = \psi^c(x) \equiv C\,\bar\psi^T(x) .
\tag{43.1}\label{eq:majorana_condition}
\end{equation}
\end{mdframed}
Here $C=i\gamma^2\gamma^0$ is the charge-conjugation matrix, satisfying
$C^T = -C$, $C^2 = -\mathbf{1}$, and $C\gamma^\mu C^{-1} = -(\gamma^\mu)^T$.

Equation \eqref{eq:majorana_condition} halves the degrees of freedom:
a Majorana field has only \emph{two} components (like a Weyl field) but is
a Lorentz singlet (not in a specific representation of the Lorentz group).

\section{Majorana vs. Weyl vs. Dirac}\label{sec:comparison}

\begin{center}
\begin{tabular}{|l|p{0.2\linewidth}|p{0.2\linewidth}|p{0.2\linewidth}|}
\hline
Property & Weyl & Majorana & Dirac \\
\hline
Spinor type & Undotted $\psi_\alpha$ & $\psi = C\bar\psi^T$ & $\psi = \psi_L + \psi_R$ \\
DOF & 2 & 2 & 4 \\
EM charge & 0 & 0 & $\pm e$ \\
Lorentz rep & $(\tfrac12,0)$ & Singlet & $(\tfrac12,0)\oplus(0,\tfrac12)$ \\
SUSY partner & $\tilde\psi$ & Gaugino $\lambda=\lambda^c$ & $\tilde f_L+\tilde f_R$ \\
\hline
\end{tabular}
\end{center}

A Weyl field is a \emph{complex} spinor: $\psi_L \neq (\psi_L)^*$.
A Majorana field is a \emph{real} spinor: $\psi = (\psi)^*$.

\section{Majorana Propagator}\label{sec:propagator}

Despite having only two components, the Majorana propagator has the same
Dirac form:
\begin{mdframed}
\begin{equation}
S_F(x-y) = \langle0|T\{\psi(x)\bar\psi(y)\}|0\rangle
         = \frac{i}{k\!\!\!/ - m + i\varepsilon} .
\tag{43.2}
\end{equation}
\end{mdframed}
However, there is a crucial difference: the \emph{Wick contraction} of two
Majorana fields gives an extra factor of $\tfrac12$:
\begin{equation}
\wick{\langle0|T\{\,\psi(x)\,\psi(y)\,}|0\rangle}
= \frac12\,S_F(x-y) .
\tag{43.3}\label{eq:majorana_contraction}
\end{equation}
The factor $\tfrac12$ is the statistical origin of the difference between
Bose-Einstein and Fermi-Dirac statistics for identical particles:
a Majorana particle is its own antiparticle, so there is no distinction
between $\psi(x)\psi(y)$ and $\psi(y)\psi(x)$.

\section{Feynman Rules for Majorana Fields}\label{sec:rules}

\begin{center}
\begin{tabular}{|l|c|}
\hline
Object & Rule \\
\hline
Majorana propagator & $\displaystyle\frac{i}{k\!\!\!/+m}$ \\
Majorana vertex & $-iy\,\gamma^\mu$ (Yukawa) \\
Closed Majorana loop & $\tfrac12\Tr[S_F]$ (extra $\tfrac12$) \\
4-Majorana contact & $-iy^2$ \\
\hline
\end{tabular}
\end{center}

\section{Majorana Mass Terms}\label{sec:mass}

A Dirac mass term is $\mathcal{L} \supset -m_D\,\bar\psi\psi = -m_D(\psi_L\psi_R + \psi_R\psi_L)$.
For a Majorana field $\psi = \psi_L + \psi_L^c$, this becomes:
\begin{mdframed}
\begin{equation}
\mathcal{L}\supset -\tfrac12 m_M\,\psi^T C\psi + \text{h.c.}
\tag{43.4}\label{eq:majorana_mass}
\end{equation}
\end{mdframed}
The $\tfrac12$ factor reflects the Majorana nature: only half the
degrees of freedom contribute (since particle = antiparticle).

In the MSSM, Majorana mass terms appear for gauginos:
$\mathcal{L} \supset \tfrac12 M_1 \tilde\lambda\tilde\lambda + \text{h.c.}$
These give gluino (colored) and neutralino/charginino (electroweak) masses.

\section{Why Majorana Matter for MHV Research}\label{sec:mhv_relevance}

In the MSSM, the MHV formalism extends to processes with squarks and
gluinos.  These are scalar partners of quarks (squarks: $h=\pm\tfrac12$)
and Majorana fermion partners of gluons (gluinos: $h=\pm1$).  The
gluino being a Majorana field means it can carry color charge \emph{and}
helicity in a way that is absent for Dirac fermions, leading to distinctive
double-pole signatures at the LHC.

The spinor-helicity formalism generalizes to:
\begin{itemize}
\item Gluino lines: use Majorana propagator with $\tfrac12$ factor
\item Gluino vertices: $\gamma^\mu$ structure (same as quark-gluon coupling)
\item Gluino spin sums: $2\times$ richer than Dirac (4 helicity states per line)
\end{itemize}

The BCFW recursion for SUSY theories (N=4 SYM) uses the fact that all
particles in the theory can be represented as massless states with definite
helicity, making the Majorana property not an obstacle but a feature.

\end{document}
"""
with open("/work/ch43_feynman_rules_majorana.tex", "w") as f:
    f.write(__import__('textwrap').dedent(r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 43: Majorana Fields}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Ch.43 covers a special spinor: the Majorana field,}
where the particle is its own antiparticle.  Appears in SUSY (gauginos,
higgsinos) and neutralino dark matter.
\end{mdframed}

\section{Majorana Condition}
\begin{mdframed}
$\displaystyle\psi(x) = \psi^c(x) \equiv C\,\bar\psi^T(x)$ ,
\quad $C=i\gamma^2\gamma^0$, $C^2=-\mathbf{1}$, $C\gamma^\mu C^{-1}=-(\gamma^\mu)^T$ .
\end{mdframed}
This halves the DOF: 2 (like Weyl) not 4 (like Dirac).

\section{Comparison}
\begin{center}
\begin{tabular}{|l|c|c|c|}
\hline
Property & Weyl & Majorana & Dirac \\
\hline
DOF & 2 & 2 & 4 \\
EM charge & 0 & 0 & $\pm e$ \\
Lorentz rep & $(\tfrac12,0)$ & Singlet & $(\tfrac12,0)\oplus(0,\tfrac12)$ \\
SUSY partner & $\tilde\psi$ & Gaugino & $\tilde f_L+\tilde f_R$ \\
\hline
\end{tabular}
\end{center}

\section{Majorana Propagator}
\begin{mdframed}
$S_F(k) = i/(k\!\!\!/+m)$.  But Wick contraction gives $\tfrac12 S_F$ (extra $\tfrac12$ factor).
\end{mdframed}
This $\tfrac12$ is the statistical origin of Majorana statistics.

\section{Feynman Rules}
\begin{center}
\begin{tabular}{|l|c|}
\hline
Object & Rule \\
\hline
Majorana propagator & $i/(k\!\!\!/+m)$ \\
Majorana Yukawa vertex & $-iy\,\gamma^\mu$ \\
Closed Majorana loop & $\frac12\Tr[S_F]$ \\
4-Majorana contact & $-iy^2$ \\
\hline
\end{tabular}
\end{center}

\section{Majorana Mass}
$\mathcal{L}\supset -\tfrac12 m_M\,\psi^TC\psi + \text{h.c.}$
The $\tfrac12$ reflects: only half the DOF (particle=antiparticle).

MSSM examples: $\tfrac12 M_1\tilde\lambda\tilde\lambda$ (gaugino Majorana mass).

\section{Why Majorana Matters for MHV}
In SUSY extensions of QCD, gluinos are Majorana fermions.
They enable double-pole signatures at the LHC.
The BCFW recursion generalizes to SUSY theories using massless
helicity states: the Majorana property is not an obstacle but a feature.

\vfill\bibliographystyle{plain}\end{document}"""))
print("Wrote /work/ch43_feynman_rules_majorana.tex")
