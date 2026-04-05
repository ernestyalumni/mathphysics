r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 41: Gamma Matrix Technology}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Ch.41 is the reference manual for all $\gamma$-matrix algebra.}
Every equation in Ch.39, 40, and 43 uses these identities.  Keep this
chapter nearby throughout the entire QFT course.
\end{mdframed}

\section{Weyl (Chiral) Basis}\label{sec:chiral}

In the chiral (Weyl) basis, the $\gamma$ matrices are block-off-diagonal:
\begin{align}
\gamma^0 &= \begin{pmatrix}0&I\\I&0\end{pmatrix}, &
\gamma^i &= \begin{pmatrix}0&-\sigma^i\\+\sigma^i&0\end{pmatrix},
\qquad i=1,2,3, \\[2mm]
\gamma^5 &= i\gamma^0\gamma^1\gamma^2\gamma^3
         = \begin{pmatrix}-I&0\\0&+I\end{pmatrix} .
\end{align}

The chiral projectors are:
\begin{mdframed}
\begin{align}
P_L = \frac{1-\gamma^5}{2} = \begin{pmatrix}I&0\\0&0\end{pmatrix},
\qquad
P_R = \frac{1+\gamma^5}{2} = \begin{pmatrix}0&0\\0&I\end{pmatrix} .
\end{align}
\end{mdframed}

A Dirac field splits into left/right Weyl components:
\begin{equation}
\psi_L = P_L\psi,\qquad \psi_R = P_R\psi,\qquad
\psi = \psi_L + \psi_R .
\end{equation}

\section{Clifford Algebra}\label{sec:clifford}
The $\gamma$ matrices obey:
\begin{mdframed}
\begin{equation}
\{\gamma^\mu,\gamma^\nu\} = 2\eta^{\mu\nu}\mathbf{1}_4 .
\tag{41.1}
\end{equation}
\end{mdframed}
Consequences: $(\gamma^0)^2=\mathbf{1}$, $(\gamma^i)^2=-\mathbf{1}$.

\section{$\gamma^5$ and Chiral Symmetry}\label{sec:gamma5}
\begin{equation}
\{\gamma^5,\gamma^\mu\} = 0,\qquad (\gamma^5)^2 = +\mathbf{1}_4 .
\end{equation}
The chiral anomaly (Adler-Bell-Jackiw):
\begin{mdframed}
\begin{equation}
\partial_\mu j^{5\mu} = \frac{e^2}{16\pi^2}\Tr[\gamma^5 F_{\mu\nu}F^{\mu\nu}]
= \frac{e^2}{8\pi^2}F_{\mu\nu}\tilde F^{\mu\nu} .
\tag{41.$\star$}
\end{equation}
\end{mdframed}
This is not a symmetry of the quantum theory --- it is an exact result.

\section{Fierz Identities}\label{sec:fierz}
For four-fermion operators $O_i$ built from $\gamma$ matrices:
\begin{mdframed}
\begin{equation}
\bigl[\bar\psi_A\gamma^\mu(1\pm\gamma^5)\psi_B\bigr]
\bigl[\bar\psi_C\gamma_\mu(1\mp\gamma^5)\psi_D\bigr]
= \bigl[\bar\psi_A\gamma^\mu(1\pm\gamma^5)\psi_D\bigr]
  \bigl[\bar\psi_C\gamma_\mu(1\mp\gamma^5)\psi_B\bigr] .
\tag{41.$\star\star$}
\end{equation}
\end{mdframed}
This rearrangement is used in: Fermi theory (beta decay), the CKM matrix,
and effective SUSY Lagrangians.

\section{Closed Fermion Loops (Ch.54 onward)}
For loop calculations, the fundamental objects are:
\begin{align}
\Tr[\gamma^\mu\gamma^\nu] &= 4\eta^{\mu\nu}, \\
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4\bigl(\eta^{\mu\nu}\eta^{\rho\sigma}
           -\eta^{\mu\rho}\eta^{\nu\sigma}
           +\eta^{\mu\sigma}\eta^{\nu\rho}\bigr) .
\end{align}

\section{Significance for MHV}
The gamma-matrix identities here are the algebraic backbone of every
spinor amplitude with fermions.  In the chiral basis:
\begin{itemize}
\item $\gamma^\mu = \begin{pmatrix}0&\sigma^\mu\\\bar\sigma^\mu&0\end{pmatrix}$
\item $\bar\sigma^{\dot\alpha\alpha} = \varepsilon^{\alpha\beta}\varepsilon^{\dot\alpha\dot\beta}\sigma_{\beta\dot\beta}$
\end{itemize}
This is the bridge between 4-component Dirac notation and the
2-component Weyl/spinor-helicity notation of Ch.42.

\vfill\bibliographystyle{plain}\bibliography{refs}
\end{document}
"""
import os
with open("/work/ch41_gamma_technology.tex", "w") as f:
    f.write(__import__('textwrap').dedent(r"""\documentclass[12pt]{article}
\usepackage[margin=1in,a4paper]{geometry}
\usepackage{amsmath,amssymb,fullpage,xcolor,mhchem,mdframed}
\usepackage{hyperref}
\DeclareMathOperator\Tr{Tr}
\begin{document}
\title{Chapter 41: Gamma Matrix Technology}
\date{\today}
\maketitle

\begin{mdframed}[backgroundcolor=blue!5]
\textbf{Ch.41 is the reference manual for all $\gamma$-matrix algebra.}
Every equation in Ch.39, 40, and 43 uses these identities.
\end{mdframed}

\section{Weyl (Chiral) Basis}\label{sec:chiral}
In the chiral (Weyl) basis:
\begin{align}
\gamma^0 &= \begin{pmatrix}0&I\\I&0\end{pmatrix}, &
\gamma^i &= \begin{pmatrix}0&-\sigma^i\\+\sigma^i&0\end{pmatrix}, \\
\gamma^5 = i\gamma^0\gamma^1\gamma^2\gamma^3
         &= \begin{pmatrix}-I&0\\0&+I\end{pmatrix} .
\end{align}
Chiral projectors: $P_L=(1-\gamma^5)/2$, $P_R=(1+\gamma^5)/2$.

\section{Clifford Algebra}
$\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}\mathbf{1}_4$.
Consequences: $(\gamma^0)^2=+\mathbf{1}$, $(\gamma^i)^2=-\mathbf{1}$.

\section{Trace Theorems (Complete Set)}
\begin{align}
\Tr[\mathbf{1}] &= 4, &
\Tr[\gamma^\mu\gamma^\nu] &= 4\eta^{\mu\nu}, \\
\Tr[\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma]
&= 4(\eta^{\mu\nu}\eta^{\rho\sigma}
   -\eta^{\mu\rho}\eta^{\nu\sigma}
   +\eta^{\mu\sigma}\eta^{\nu\rho}) .
\end{align}

\section{Chiral Anomaly (Adler-Bell-Jackiw)}
\begin{mdframed}
$\displaystyle\partial_\mu j^{5\mu}
= \frac{e^2}{16\pi^2}\Tr[\gamma^5 F_{\mu\nu}F^{\mu\nu}]
= \frac{e^2}{8\pi^2}F_{\mu\nu}\tilde F^{\mu\nu}$ .
\end{mdframed}
This is an exact quantum anomaly --- not a symmetry of the Lagrangian
but a consequence of regularization.

\section{Fierz Identities}
\begin{mdframed}
$[\bar\psi_A\gamma^\mu(1\pm\gamma^5)\psi_B][\bar\psi_C\gamma_\mu(1\mp\gamma^5)\psi_D]
=[\bar\psi_A\gamma^\mu(1\pm\gamma^5)\psi_D][\bar\psi_C\gamma_\mu(1\mp\gamma^5)\psi_B]$ .
\end{mdframed}
Used in: Fermi theory, CKM matrix elements, SUSY superpotential matching.

\section{$\gamma^5$ and the Chiral Basis}
In the chiral basis, the 4-component notation splits cleanly:
$\gamma^\mu=(0,\sigma^\mu;\bar\sigma^\mu,0)$.
The connection to 2-component Weyl spinors is:
$u_+(k)\to\lambda_\alpha$, $u_-(k)\to\tilde\lambda_{\dot\alpha}$.
This is the bridge to the spinor-helicity formalism of Ch.42.

\vfill\bibliographystyle{plain}\end{document}"""))
print(f"Wrote /work/ch41_gamma_technology.tex")
