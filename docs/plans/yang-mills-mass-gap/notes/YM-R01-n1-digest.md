# YM-R01 — Digest of N1: "Exponential decay in two-dimensional classical O(n) models"

- Author: Claude Opus 5.5 subagent (claude-opus55-r01)
- Date: 2026-11-17
- Source pin: openai/math 4eaee87, directory `preprints/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026/` (files `build/main.tex`, `build/sections/rotation.tex`, `build/sections/annulus.tex`, `build/sections/clusters.tex`; `paper.pdf`, 20 pages, used for equation and page numbers)
- Status: unreviewed; needs review by a different model family

## 0. Conventions

- Every fenced `latex` block below is copied verbatim from the named source file and line range (extracted mechanically with `sed -n`, not retyped). Prose outside the fences is mine; inline formulas in prose are paraphrase and carry no claim of verbatim match.
- Equation numbers are the printed numbers in `paper.pdf`. The source sets `\numberwithin{equation}{section}`, and all theorem-like environments share one counter per section (main.tex lines 14–21, 30). Section 1 is the introduction, so the rotation section is §2, the annulus §3, sign clusters §4.
- Page numbers are `paper.pdf` page numbers, read from `pdftotext -layout`.
- Comparison with the packet's assumptions: the packet's (2.1), (2.8), (2.10), (2.11), (2.12) and "Proposition 2.2" coincide with the PDF numbering. (2.1) is the definition of H, W; (2.8) is the inequality of Proposition 2.2; (2.10), (2.11), (2.12) are the non-abelian gauge identity, the quadratic expansion of B, and the fourth-order cancellation. The "second-variation identity" has no number of its own: it is the equality in (2.6) together with the unnumbered display in the proof of Lemma 2.1 (p. 5).

### Label concordance

| Printed | Source label | Page | Content |
|---|---|---|---|
| (1.1) | `eq:model` | 1 | free O(n) Gibbs measure |
| (1.2) | `eq:main-bound` | 1 | main two-point bound |
| (1.3) | `eq:limsup` | 1 | infinite-volume limsup |
| Thm 1.1 / Cor 1.2 | `thm:main` / `cor:limits` | 1 | main results |
| (2.1) | `eq:rotation-partition` | 4 | H(A), W(A) |
| (2.2) | `eq:rotation-local-coefficients` | 4 | J_e |
| (2.3) | `eq:rotation-local-expansion` | 4 | H_k expansion |
| (2.4) | `eq:rotation-R-definition` | 4 | R(h,k) |
| Lemma 2.1 | `lem:gradient` | 4 | gradient estimates |
| (2.5) | `eq:rotation-gradient-bound` | 4 | |
| (2.6) | `eq:rotation-R-bound` | 4 | second variation of R |
| (2.7) | `eq:rotation-S-definition` | 5 | S(f,g) |
| Prop 2.2 | `prop:fourth` | 5 | fourth-order rotation inequality |
| (2.8) | `eq:rotation-family-bound` | 5 | family inequality |
| (2.9) | `eq:rotation-explicit-constant` | 5 | K_rot |
| (2.10) | `eq:rotation-nonabelian-gauge` | 5 | W(dF)=W(B) |
| (2.11) | `eq:rotation-B-quadratic` | 5 | B quadratic |
| (2.12) | `eq:rotation-B-fourth` | 6 | fourth-order cancellation |
| (2.13) | `eq:rotation-four-directions` | 6 | |
| (2.14) | `eq:rotation-fourth-identity` | 6 | |
| (2.15) | `eq:rotation-response-balance` | 6 | |
| (2.16) | `eq:rotation-gradient-adjustment` | 6 | |
| (2.17) | `eq:rotation-cross-bound` | 6 | |
| (2.18) | `eq:rotation-fourth-lower` | 6 | |
| (2.19) | `eq:rotation-fifteen-terms` | 7 | |
| (2.20) | `eq:rotation-positive-square` | 7 | |
| (2.21) | `eq:rotation-row-gradient` | 7 | |
| (2.22) | `eq:rotation-row-bound` | 7 | |
| Prop 3.1 | `prop:annulus` | 8 | |
| (3.1), (3.2) | `eq:annulus-minmax`, `eq:annulus-ratio` | 8 | |
| (3.3)–(3.7) | `eq:annulus-h`, `-boundary-gauge`, `-profiles`, `-modes`, `-htilde` | 9–10 | |
| Lemma 3.2 | `lem:modes` | 10 | |
| (3.8)–(3.10) | `eq:annulus-commutator`, `-midpoint`, `-test-cost` | 10 | |
| (3.11)–(3.15) | `eq:annulus-one-mode`, `-mode-norms`, `-minimum-sum`, `-fourier-norm`, `-operator-factorization` | 10–11 | |
| (3.16)–(3.20) | `eq:annulus-tested-response`, `-current-control`, `-current-error`, `-local-error`, `-ode` | 12 | |
| (4.1)–(4.6) | `eq:spin-disintegration`, `amplitude-couplings`, `amplitude-density`, `bond-parameter`, `ising-bond-expansion`, `spin-connection` | 13 | |
| Lemmas 4.1–4.4 | `lem:association`, `lem:edge-correlation-inequalities`, `lem:amplitude-association`, `lem:joint-association` | 14–15 | |
| (4.7), (4.8) | `eq:doubled-covariance`, `eq:conditional-fk-weight` | 14, 15 | |
| Props 4.5–4.7 | `prop:pinning`, `prop:crossings`, `prop:connection` | 16–17 | |
| (4.9) | `eq:crossing-partition-ratio` | 17 | |
| (4.10), (4.11) | `eq:coarse-path-alpha`, `eq:coarse-path-length` | 17, 18 | |
| Lemma 4.8 | `lem:conditioning` | 18 | |

## 1. Main result and architecture (§1, pp. 1–3)

The model (1.1):

`build/main.tex` lines 64–70:

```latex
\begin{equation}
 \dd\mu_{G,b}^{(n)}(\sigma)
 =\frac{1}{Z_{G,b}^{(n)}}
   \exp\!\left(\sum_{\{x,y\}\in E}b_{xy}\sigma_x\cdot\sigma_y\right)
   \prod_{x\in V}\dd\omega_{n-1}(\sigma_x),
 \label{eq:model}
\end{equation}
```

Theorem 1.1 (p. 1), quoted:

`build/main.tex` lines 76–87:

```latex
\begin{theorem}
\label{thm:main}
For every integer $n\ge3$ and every finite $\beta>0$, there are constants
$A=A(n,\beta)<\infty$ and $m=m(n,\beta)>0$ such that, for every finite
nearest-neighbor square-lattice subgraph $G$, every strength array
$b\in[0,\beta]^E$, and all $x,y\in V$,
\begin{equation}
 0\le \langle\sigma_x\cdot\sigma_y\rangle_{G,b}^{(n)}
 \le A\exp(-m\|x-y\|_2).
 \label{eq:main-bound}
\end{equation}
\end{theorem}
```

Corollary 1.2 (p. 1), quoted:

`build/main.tex` lines 89–101:

```latex
\begin{corollary}[Free-boundary limits]\label{cor:limits}
Let $\Lambda_L=[-L,L]^2\cap\Z^2$, with all its nearest-neighbor edges,
free boundary conditions, and constant strength $\beta>0$. For every
$n\ge3$ and $x\in\Z^2$,
\begin{equation}
 \limsup_{L\to\infty}
 \left|\langle\sigma_0\cdot\sigma_x\rangle_{\Lambda_L,\beta}^{(n)}\right|
 \le A e^{-m\|x\|_2},
 \label{eq:limsup}
\end{equation}
with the constants of Theorem~\ref{thm:main}. Every subsequential local
weak limit of these laws satisfies the same two-point bound.
\end{corollary}
```

Architecture as the paper describes it (§1.3, p. 3): (i) for n = 3, on any finite graph with arbitrary pinned spins, a fourth derivative of the twisted partition function in two noncommuting rotation directions yields an upper bound on the second variation along the commutator direction (Prop. 2.2); (ii) on a square annulus of inner radius N with boundary layers pinned at relative angle t, a multiscale test family turns that into a curvature bound giving Z(π)/Z(0) ≥ 1 − O_β(1/log N) (Prop. 3.1); (iii) Z(π)/Z(0) is identified with 1 − (crossing probability) in a sign-cluster representation, pinning comparison makes crossings of disjoint annuli independent and individually small, and a coarse-path sum gives exponential decay of connectivity, hence of ⟨s^1_x s^1_y⟩ (Props. 4.5–4.7); (iv) conditioning on coordinates 4..n reduces O(n) to O(3) with weakened couplings (Lemma 4.8). The paper stresses that all estimates are uniform over finite subgraphs of Z² and over strength arrays in [0,β]^E; this uniformity is what makes step (iv) work (p. 2, last paragraph of §1.1).

The paper's stated limits (p. 2): the result concerns only the spin two-point function; a full transfer gap "requires a separate argument"; no correlation-length asymptotics or uniqueness of infinite-volume Gibbs states; constants may deteriorate as β → ∞. The paper says Theorem 1.1 cannot extend to n = 2 by Fröhlich–Spencer (p. 2).

## 2. §2 — A fourth-order rotation estimate (pp. 3–8)

### 2.1 Setting and definitions (p. 4)

`build/sections/rotation.tex` lines 17–23:

```latex
Let \(G=(V,E)\) be a finite graph, with one orientation \(e=(x,y)\)
chosen for each edge, and let \(0\le b_e\le\beta\). Fix a set
\(P\subseteq V\) of vertices and arbitrary spins \(s_x\in S^2\) for
\(x\in P\). At the remaining vertices, integrate against uniform
probability measure on \(S^2\), with normalized density proportional to
\(\exp(\sum_{e=(x,y)}b_e s_x^Ts_y)\). Write \(\langle\cdot\rangle\)
for this expectation.
```

Definitions of the generators and of H(A), W(A), eq. (2.1):

`build/sections/rotation.tex` lines 25–41:

```latex
Let \(a,b,c\in\mathfrak{so}(3)\) denote cross product with the first,
second, and third coordinate vectors, respectively. Thus
\([a,b]=c\), and their operator norms are one. For an edge field
\(A=(A_e)_{e\in E}\) of skew matrices, put
\begin{equation}
\label{eq:rotation-partition}
\begin{split}
H(A)&=\sum_{e=(x,y)}b_e s_x^T\bigl(\exp(A_e)-I\bigr)s_y,\\
W(A)&=\bigl\langle\exp H(A)\bigr\rangle.
\end{split}
\end{equation}
Thus \(W\) is the ratio of the twisted partition function to the
untwisted one. All derivatives below are evaluated at \(A=0\). Set
\(X=DH\) and \(H_k=D^kH\) for \(k=2,3,4\). Compactness of the
finite spin space permits differentiation under the integral. We extend
all multilinear derivatives complex linearly in their arguments; changes
of spin variables themselves will always use real rotations.
```

So W(A) is the ratio of the partition function with every edge interaction b_e s_x^T s_y replaced by b_e s_x^T exp(A_e) s_y to the untwisted partition function; X = DH and H_k = D^kH are evaluated at A = 0. All derivatives with complex arguments are the complex-multilinear extensions of the real derivatives.

Local coefficients (2.2) and the expansion (2.3):

`build/sections/rotation.tex` lines 48–63:

```latex
\begin{equation}
\label{eq:rotation-local-coefficients}
J_e(g_1,\ldots,g_k)
=b_e s_x^T\left(\frac{1}{k!}\sum_{\pi\in\mathfrak S_k}
g_{\pi(1)}\cdots g_{\pi(k)}\right)s_y,
\qquad |J_e(g_1,\ldots,g_k)|\le\beta,
\end{equation}
where each \(g_r\) is one of \(a,b,c\). Indeed, every matrix product
in this average has operator norm at most one. Consequently
\begin{equation}
\label{eq:rotation-local-expansion}
H_k(v_1g_1,\ldots,v_kg_k)
=\sum_e J_e(g_1,\ldots,g_k)\prod_{r=1}^k v_r(e).
\end{equation}
The coefficients depend on the spins and the generators, but not on the
scalar test functions. The same formula with \(k=1\) describes \(X\).
```

(Reason for (2.3), in my words: the coefficient of t_1⋯t_k in exp(Σ_r t_r A_r) is the symmetrized product (1/k!) Σ_π A_{π(1)}⋯A_{π(k)}, so the mixed k-th derivative of H at 0 is Σ_e b_e s_x^T (that average) s_y; with A_r = v_r g_r edgewise this is (2.3). The bound |J_e| ≤ β uses |s_x| = |s_y| = 1, operator norm one of every product of a, b, c, and 0 ≤ b_e ≤ β.)

The response form (2.4):

`build/sections/rotation.tex` lines 65–71:

```latex
Define the complex bilinear form
\begin{equation}
\label{eq:rotation-R-definition}
R(h,k)=-D^2W(hc,kc)
\end{equation}
on scalar edge forms. It is symmetric, but is not assumed to be positive
semidefinite.
```

Notation used later: for a scalar edge function v and g ∈ {a,b,c}, X_g(v) = X(vg) (introduced p. 7). Because J_e(g) = b_e s_x^T g s_y is real, X(v̄g) is the complex conjugate of X(vg), and likewise H_2(v̄g, w̄g') is the conjugate of H_2(vg, wg'). This conjugation rule is used, without being stated as a separate step, in (2.6), (2.20) and the table on p. 7.

### 2.2 Statements

Lemma 2.1 (p. 4), quoted:

`build/sections/rotation.tex` lines 73–90:

```latex
\begin{lemma}[Gradient estimates]\label{lem:gradient}
If \(f:V\to\mathbb C\) vanishes on \(P\), then, for
\(g\in\{a,b,c\}\),
\begin{equation}
\label{eq:rotation-gradient-bound}
\bigl\langle |X((df)g)|^2\bigr\rangle
\le\beta\|df\|_2^2.
\end{equation}
For every complex edge form \(h\),
\begin{equation}
\label{eq:rotation-R-bound}
R(h,\bar h)
=-\bigl\langle H_2(hc,\bar h c)+|X(hc)|^2\bigr\rangle
\le\beta\|h\|_2^2.
\end{equation}
Moreover, adding \(df\) in either argument of \(R\) does not change
its value whenever \(f\) vanishes on \(P\).
\end{lemma}
```

The form S, eq. (2.7) (p. 5):

`build/sections/rotation.tex` lines 119–124:

```latex
For vertex functions \(f,g\), write
\begin{equation}
\label{eq:rotation-S-definition}
S(f,g)_e=\tfrac12(f_xg_y-g_xf_y),\qquad e=(x,y).
\end{equation}
This form is bilinear and antisymmetric in \(f,g\).
```

Proposition 2.2 (p. 5), quoted in full, containing (2.8) and (2.9):

`build/sections/rotation.tex` lines 126–152:

```latex
\begin{proposition}[Fourth-order rotation inequality]\label{prop:fourth}
Let \((f_i)_{i\in I}\) be a finite family of complex vertex functions
vanishing on \(P\), and put
\[
u_i=df_i,\qquad
T_{el}=\sum_{i\in I}u_i(e)\overline{u_i(l)},\qquad
Q=\sum_{i\in I}S(f_i,\bar f_i).
\]
Then
\begin{equation}
\label{eq:rotation-family-bound}
\begin{split}
R(Q,\bar Q)
&\le K_{\mathrm{rot}}(\beta)\|T\|_{\mathrm{HS}}^2\\
&\quad+\beta\sum_{i,j\in I}
\min\bigl\{\|f_i\|_\infty^2\|df_j\|_2^2,
             \|f_j\|_\infty^2\|df_i\|_2^2\bigr\},
\end{split}
\end{equation}
where
\begin{equation}
\label{eq:rotation-explicit-constant}
K_{\mathrm{rot}}(\beta)=6\beta^2+4\beta^{3/2}+\beta.
\end{equation}
In particular, the constants are independent of the graph, the family
size, and the values of the fixed spins.
\end{proposition}
```

Quantifiers, spelled out: G any finite graph with an orientation, any 0 ≤ b_e ≤ β, any P ⊆ V, any pinned values s_x ∈ S² on P, any finite family of complex vertex functions vanishing on P. K_rot depends only on β.

### 2.3 Mechanism

**(a) Second-variation identity and abelian gauge invariance (Lemma 2.1, pp. 4–5).** Differentiating W(A) = ⟨exp H(A)⟩ at A = 0 gives DW(A) = ⟨X(A)⟩ and D²W(A,A') = ⟨H_2(A,A') + X(A)X(A')⟩ (this is the two-index case of the set-partition expansion (2.19)). Hence R(h,k) = −⟨H_2(hc,kc) + X(hc)X(kc)⟩ and, with k = h̄ and the conjugation rule above, the equality in (2.6):

(that equality is the second display of the Lemma 2.1 quote above, `sections/rotation.tex` lines 83–87).

The abelian gauge step: for real f vanishing on P and a fixed generator g, the substitution s_x^old = exp(−t f_x g) s_x^new at integrated vertices leaves the uniform measure and the pinned spins unchanged, and turns each edge matrix exp(t(f_y − f_x)g) into exp(t f_x g) exp(t(f_y−f_x)g) exp(−t f_y g) = I. So W(t(df)g) = 1 for all t. Its second derivative is the unnumbered second-variation identity on p. 5:

`build/sections/rotation.tex` lines 98–100:

```latex
\[
0=\bigl\langle H_2((df)g,(df)g)+X((df)g)^2\bigr\rangle.
\]
```

Since |H_2((df)g,(df)g)| ≤ β‖df‖²₂ by (2.2)–(2.3), this gives (2.5) for real f; real and imaginary parts add for complex f. The inequality in (2.6) follows because −⟨|X(hc)|²⟩ ≤ 0 and |H_2(hc,h̄c)| ≤ β‖h‖²₂. Gradient invariance: W((h+df)c) = W(hc) for real h, f by the same substitution with generator c, so R(h, df) = R(df, h) = 0 after differentiating and extending complex-linearly (p. 5). The paper stresses (p. 4, after (2.4)) that R is symmetric but not assumed positive semidefinite, and (p. 8, last sentence of the proof) that only the upper bound (2.6) and gradient invariance of R are used.

**(b) Non-abelian gauge identity (2.10)–(2.12) (pp. 5–6).** Now take a real matrix potential F_x in span{a,b}, zero on P, and substitute s_x^old = exp(−F_x) s_x^new. The edge matrix exp(A_e) with A = dF becomes exp(F_x) exp(F_y − F_x) exp(−F_y), which is no longer I because a and b do not commute:

`build/sections/rotation.tex` lines 161–181:

```latex
\emph{The noncommuting gauge identity.}
Let \(F_x\) be a real matrix potential in the span of \(a,b\), zero
on \(P\). The substitution
\(s_x^{\mathrm{old}}=\exp(-F_x)s_x^{\mathrm{new}}\) gives
\begin{equation}
\label{eq:rotation-nonabelian-gauge}
W(dF)=W(B),\qquad
B_e=\log\bigl(\exp(F_x)\exp(F_y-F_x)\exp(-F_y)\bigr)
\end{equation}
for \(F\) sufficiently small. The product inside the logarithm is
orthogonal and near the identity; hence its local logarithm is skew
symmetric. The logarithm also commutes with conjugation. Multiplying
the exponential series to degree two, and then taking the logarithm,
gives
\begin{equation}
\label{eq:rotation-B-quadratic}
B_e=\tfrac12[F_x,F_y]+O(F^3),\qquad DB_e(0)=0.
\end{equation}
For directions \(fa\) and \(gb\), the mixed second derivative of
\(B_e\) is therefore \(S(f,g)_e c\); it is zero for two directions
using the same generator.
```

Check of the mixed second derivative (my computation, consistent with the text above): with F = εf a + δg b, the ε δ coefficient of ½[F_x,F_y] is ½(f_x g_y [a,b] + g_x f_y [b,a]) = ½(f_x g_y − g_x f_y) c = S(f,g)_e c. For two directions with the same generator the commutator vanishes.

Fourth-order cancellation, eq. (2.12):

`build/sections/rotation.tex` lines 183–201:

```latex
We also need a fourth-order cancellation. Take four independent
parameters whose matrix directions have generator pattern \((a,b,a,b)\),
with arbitrary spatial potentials. The corresponding mixed fourth
derivative of \(B_e\) is zero. To see this, conjugate by the half-turn
about the first axis. This fixes \(a\) and negates \(b\), so the
coefficient involving all four parameters is fixed by conjugation:
its two \(b\)-directions contribute two minus signs. Conjugation by
the half-turn about the second axis likewise fixes this coefficient.
On \(\mathfrak{so}(3)\), the fixed spaces of these two conjugations
are respectively \(\mathbb R a\) and \(\mathbb R b\), with zero
intersection. Thus
\begin{equation}
\label{eq:rotation-B-fourth}
D_{1234}B_e=0
\quad\text{for the generator pattern }(a,b,a,b).
\end{equation}
Both \eqref{eq:rotation-B-quadratic} and
\eqref{eq:rotation-B-fourth} extend complex multilinearly. No complex
change of spin variables is involved in this extension.
```

Mechanism in my words: B is equivariant under conjugation, B(θFθ⁻¹) = θB(F)θ⁻¹, because exp and log commute with conjugation. The half-turn θ₁ about the first axis acts on so(3) by a ↦ a, b ↦ −b, c ↦ −c. Replacing each of the four directions (a,b,a,b) by its conjugate multiplies the all-four mixed coefficient by (+1)(−1)(+1)(−1) = 1, so that coefficient (an element of so(3), since B is skew) is fixed by Ad θ₁, i.e. lies in ℝa. The same with θ₂ puts it in ℝb. Hence it is zero.

**(c) The fourth-derivative identity (2.13)–(2.15) (p. 6).** For fixed i, j take the four real-parameter potentials f_i a, f̄_i b, f̄_j a, f_j b (complex-multilinear extension), whose edge gradients are:

`build/sections/rotation.tex` lines 203–208:

```latex
For fixed \(i,j\), take the four edge directions
\begin{equation}
\label{eq:rotation-four-directions}
A_1=u_i a,\qquad A_2=\bar u_i b,\qquad
A_3=\bar u_j a,\qquad A_4=u_j b.
\end{equation}
```

The left side of (2.10) is W of a linear function of the parameters, so its mixed fourth derivative is D⁴W(A_1,A_2,A_3,A_4). The right side is W∘B; its mixed fourth derivative is the sum over set partitions of {1,2,3,4} of D^{#blocks}W applied to the block derivatives of B (Faà di Bruno). Any block of size one carries DB(0) = 0 by (2.11), so only {1234} and the three pairings survive; {1234} gives DW[D_{1234}B] = 0 by (2.12); the pairing 13|24 pairs a with a and b with b, so both D_{13}B and D_{24}B vanish. The remaining pairings:

`build/sections/rotation.tex` lines 209–230:

```latex
The fourth derivative of the composition in
\eqref{eq:rotation-nonabelian-gauge} has no term containing \(DB\),
by \eqref{eq:rotation-B-quadratic}. Its term \(DW[D_{1234}B]\)
also vanishes, by \eqref{eq:rotation-B-fourth}. Only pairings of two
second derivatives remain. The pairing \(13|24\) is zero because
both pairs have the same generator. For the other two pairings, put
\(S_i=S(f_i,\bar f_i)\) and \(S_{ij}=S(f_i,f_j)\). Their arguments
are
\[
12|34:\quad(S_i c,\overline{S_j}c),\qquad
14|23:\quad(S_{ij}c,-\overline{S_{ij}}c).
\]
Here \(S(\bar f_j,f_j)=\overline{S_j}\) and
\(S(\bar f_j,\bar f_i)=-\overline{S_{ij}}\). Each labeled pairing
occurs once in the composition rule, so
\begin{equation}
\label{eq:rotation-fourth-identity}
D^4W(A_1,A_2,A_3,A_4)
=-R(S_i,\overline{S_j})+R(S_{ij},\overline{S_{ij}}).
\end{equation}
This also holds when some spatial functions coincide: the derivative
identity is multilinear in the four labeled directions.
```

(Check, mine: D_{23}B comes from directions f̄_i b and f̄_j a, giving ½(f̄_{i,x} f̄_{j,y}[b,a] + f̄_{j,x} f̄_{i,y}[a,b]) = S(f̄_j, f̄_i)c = −conj(S_{ij}) c; and S(f̄_j, f_j) = conj(S(f_j, f̄_j)). Each term D²W(·c,·c) is −R(·,·) by (2.4).)

Summing (2.14) over i, j and using bilinearity Σ_{i,j} R(S_i, S̄_j) = R(Q, Q̄):

`build/sections/rotation.tex` lines 232–241:

```latex
Summing over \(i,j\) and rearranging gives the identity to be estimated:
\begin{equation}
\label{eq:rotation-response-balance}
R(Q,\bar Q)
=\sum_{i,j}R(S_{ij},\overline{S_{ij}})
 -\sum_{i,j}D^4W(u_i a,\bar u_i b,\bar u_j a,u_j b).
\end{equation}
We bound the first sum from above using gradient invariance, and the
second from below by retaining a nonnegative square.

```

**(d) Cross terms (2.16)–(2.17) (p. 6).**

`build/sections/rotation.tex` lines 243–260:

```latex
\emph{The cross terms.}
For functions \(f,g\) vanishing on \(P\), the product \(fg\)
also vanishes there, and the edge identities
\begin{equation}
\label{eq:rotation-gradient-adjustment}
\begin{split}
S(f,g)+\tfrac12d(fg)&=\tfrac12(f_x+f_y)\,dg,\\
S(f,g)-\tfrac12d(fg)&=-\tfrac12(g_x+g_y)\,df
\end{split}
\end{equation}
hold on \(e=(x,y)\). Applying gradient invariance in both arguments
of \(R\), followed by \eqref{eq:rotation-R-bound}, gives
\begin{equation}
\label{eq:rotation-cross-bound}
R(S(f,g),\overline{S(f,g)})
\le\beta\min\{\|f\|_\infty^2\|dg\|_2^2,
                 \|g\|_\infty^2\|df\|_2^2\}.
\end{equation}
```

In my words: S(f,g) differs from ½(f_x+f_y) dg by the gradient −½d(fg) of a function vanishing on P (and S(f,g)‾ from its conjugate counterpart by −½d(f̄ḡ)), so gradient invariance in both arguments of R lets one replace S(f,g) by ½(f_x+f_y) dg, after which (2.6) gives β‖½(f_x+f_y)dg‖²₂ ≤ β‖f‖²_∞‖dg‖²₂; symmetrically with the second line of (2.16). Applied with f = f_i, g = f_j this bounds Σ_{i,j} R(S_{ij}, S̄_{ij}) by the second line of (2.8).

**(e) Lower bound on the fourth-derivative sum (2.18)–(2.22) (pp. 6–8).** Claim (2.18):

`build/sections/rotation.tex` lines 263–269:

```latex
\emph{The fourth-derivative lower bound.}
We claim that
\begin{equation}
\label{eq:rotation-fourth-lower}
\sum_{i,j}D^4W(u_i a,\bar u_i b,\bar u_j a,u_j b)
\ge-K_{\mathrm{rot}}(\beta)\|T\|_{\mathrm{HS}}^2.
\end{equation}
```

Expansion of the fourth derivative into the fifteen set partitions, eq. (2.19):

`build/sections/rotation.tex` lines 270–284:

```latex
To enumerate its terms, set \(X_r=X(A_r)\), and write
\(H_{r_1\cdots r_k}=H_k(A_{r_1},\ldots,A_{r_k})\). Direct
differentiation of \(\exp H\) gives all fifteen set partitions:
\begin{equation}
\label{eq:rotation-fifteen-terms}
\begin{aligned}
D^4W(A_1,A_2,A_3,A_4)=\bigl\langle
&X_1X_2X_3X_4\\
&+H_{12}X_3X_4+H_{13}X_2X_4+H_{14}X_2X_3\\
&+H_{23}X_1X_4+H_{24}X_1X_3+H_{34}X_1X_2\\
&+H_{12}H_{34}+H_{13}H_{24}+H_{14}H_{23}\\
&+H_{123}X_4+H_{124}X_3+H_{134}X_2+H_{234}X_1\\
&+H_{1234}\bigr\rangle.
\end{aligned}
\end{equation}
```

The nonnegative square, eq. (2.20):

`build/sections/rotation.tex` lines 285–294:

```latex
The four terms comprising
\((X_1X_2+H_{12})(X_3X_4+H_{34})\), after summation over \(i,j\),
have expectation
\begin{equation}
\label{eq:rotation-positive-square}
\left\langle\left|\sum_i
\bigl[X(u_i a)X(\bar u_i b)+H_2(u_i a,\bar u_i b)\bigr]
\right|^2\right\rangle\ge0.
\end{equation}
We bound the other eleven terms in absolute value.
```

(Why it is a square, my words: X_1X_2 + H_{12} = X(u_i a)X(ū_i b) + H_2(u_i a, ū_i b), while X_3X_4 + H_{34} = X(ū_j a)X(u_j b) + H_2(ū_j a, u_j b) is the complex conjugate of the same expression with index j, by the conjugation rule in §2.1. Summing over i, j gives |Σ_i[…]|².)

The key structural fact for the other eleven terms, eqs. (2.21)–(2.22):

`build/sections/rotation.tex` lines 296–310:

```latex
For an edge \(e\), let \(t_e(l)=T_{el}\) and
\(\rho_e=(\sum_l|T_{el}|^2)^{1/2}\), and abbreviate
\(X_g(v)=X(vg)\). Crucially,
\begin{equation}
\label{eq:rotation-row-gradient}
t_e=d\left(\sum_i u_i(e)\overline{f_i}\right).
\end{equation}
The potential on the right, and its conjugate, vanish on \(P\).
Lemma~\ref{lem:gradient} therefore gives
\begin{equation}
\label{eq:rotation-row-bound}
\langle|X_g(t_e)|^2\rangle\le\beta\rho_e^2,
\qquad
\langle|X_g(\bar t_e)|^2\rangle\le\beta\rho_e^2.
\end{equation}
```

In my words: T_{el} = Σ_i u_i(e) ū_i(l) = Σ_i u_i(e) (d f̄_i)(l), so for fixed e the row l ↦ T_{el} is the gradient of the vertex function Σ_i u_i(e) f̄_i, which vanishes on P. Lemma 2.1 then bounds the variance of the current X_g along that row by β times its ℓ² norm ρ_e. This is how all currents X that appear outside the square are controlled by gauge invariance rather than estimated directly.

The table of the eleven remaining contractions (p. 7) and the worked first row:

`build/sections/rotation.tex` lines 312–361:

```latex
The following table records every remaining contraction. Expand each
\(H_k\) using \eqref{eq:rotation-local-expansion}, sum over \(i,j\),
and multiply the last two columns. One then sums over \(e\), except
in the two \(H_2H_2\) rows, where one sums over \(e,l\).
\begin{center}
\renewcommand{\arraystretch}{1.18}
\begin{tabular}{@{}lll@{}}
\toprule
Term & Local coefficient & Contracted expression\\
\midrule
\(H_{13}X_2X_4\) & \(J_e(a,a)\)
  & \(X_b(t_e)X_b(\bar t_e)\)\\
\(H_{14}X_2X_3\) & \(J_e(a,b)\)
  & \(X_b(t_e)X_a(t_e)\)\\
\(H_{23}X_1X_4\) & \(J_e(b,a)\)
  & \(X_a(\bar t_e)X_b(\bar t_e)\)\\
\(H_{24}X_1X_3\) & \(J_e(b,b)\)
  & \(X_a(\bar t_e)X_a(t_e)\)\\
\midrule
\(H_{13}H_{24}\) & \(J_e(a,a)J_l(b,b)\)
  & \(T_{el}\overline{T_{el}}\)\\
\(H_{14}H_{23}\) & \(J_e(a,b)J_l(b,a)\)
  & \(T_{el}^{\,2}\)\\
\midrule
\(H_{123}X_4\) & \(J_e(a,b,a)\)
  & \(T_{ee}X_b(\bar t_e)\)\\
\(H_{124}X_3\) & \(J_e(a,b,b)\)
  & \(T_{ee}X_a(t_e)\)\\
\(H_{134}X_2\) & \(J_e(a,a,b)\)
  & \(T_{ee}X_b(t_e)\)\\
\(H_{234}X_1\) & \(J_e(b,a,b)\)
  & \(T_{ee}X_a(\bar t_e)\)\\
\midrule
\(H_{1234}\) & \(J_e(a,b,a,b)\) & \(T_{ee}^{\,2}\)\\
\bottomrule
\end{tabular}
\end{center}
For example, the first row follows from
\[
\begin{split}
&\sum_{i,j}u_i(e)\bar u_j(e)
             X_b(\bar u_i)X_b(u_j)\\
&\qquad=
X_b\left(\sum_i u_i(e)\bar u_i\right)
X_b\left(\sum_j\bar u_j(e)u_j\right)
=X_b(t_e)X_b(\bar t_e).
\end{split}
\]
The other rows follow by the same finite summation; their explicit
expressions distinguish the conjugated and unconjugated contractions.
```

The bounds (p. 8):

`build/sections/rotation.tex` lines 363–393:

```latex
For each of the first four rows, the pointwise coefficient bound and
Cauchy--Schwarz, followed by \eqref{eq:rotation-row-bound}, give an
absolute expectation at edge \(e\) of at most
\(\beta^2\rho_e^2\). The coefficient may depend on the spins;
its pointwise absolute value is bounded before applying
Cauchy--Schwarz. These four rows contribute at most
\(4\beta^2\|T\|_{\mathrm{HS}}^2\).
The next two rows contribute at most
\(2\beta^2\|T\|_{\mathrm{HS}}^2\), since both contractions have
absolute value \(|T_{el}|^2\).

For each of the four \(H_3X\) rows,
\eqref{eq:rotation-row-bound} bounds the total absolute expectation by
\[
\beta^{3/2}\sum_e T_{ee}\rho_e
\le\beta^{3/2}
\left(\sum_e T_{ee}^2\right)^{1/2}
\left(\sum_e\rho_e^2\right)^{1/2}
\le\beta^{3/2}\|T\|_{\mathrm{HS}}^2.
\]
Here \(T_{ee}=\sum_i|u_i(e)|^2\ge0\), and the squared diagonal
entries form part of the Hilbert--Schmidt sum. The last row is bounded
by \(\beta\sum_e T_{ee}^2\le\beta\|T\|_{\mathrm{HS}}^2\).

Finally, the sum on the left of \eqref{eq:rotation-fourth-lower} is
real: complex conjugation interchanges \(i,j\) and permutes the
arguments of the symmetric fourth derivative. Keeping the nonnegative
square \eqref{eq:rotation-positive-square}, and subtracting the bounds
for the remaining terms, proves \eqref{eq:rotation-fourth-lower} with
exactly \(K_{\mathrm{rot}}(\beta)\) from
\eqref{eq:rotation-explicit-constant}.
```

Assembly (p. 8): (2.15) = [sum bounded above by (2.17)] − [sum bounded below by (2.18)], which is (2.8).

### 2.4 Packet item 3: square, controlled local terms, and their bounds

- **Nonnegative square.** The four set-partition terms X_1X_2X_3X_4 + H_{12}X_3X_4 + H_{34}X_1X_2 + H_{12}H_{34} of (2.19), i.e. the product (X_1X_2 + H_{12})(X_3X_4 + H_{34}), which after Σ_{i,j} equals the expectation (2.20) ⟨|Σ_i [X(u_i a)X(ū_i b) + H_2(u_i a, ū_i b)]|²⟩ ≥ 0. It is kept, not estimated, and enters (2.18) only through its sign.
- **Controlled local terms.** The other eleven terms of (2.19), listed in the table on p. 7, in four groups:
  1. H_{13}X_2X_4, H_{14}X_2X_3, H_{23}X_1X_4, H_{24}X_1X_3 (coefficients J_e(a,a), J_e(a,b), J_e(b,a), J_e(b,b)): each ≤ β²‖T‖²_HS, using |J_e| ≤ β pointwise (2.2), Cauchy–Schwarz, and (2.22). Total 4β²‖T‖²_HS.
  2. H_{13}H_{24}, H_{14}H_{23} (coefficients J_e(a,a)J_l(b,b), J_e(a,b)J_l(b,a); contractions T_{el}T̄_{el}, T_{el}²): each ≤ β²‖T‖²_HS. Total 2β²‖T‖²_HS.
  3. H_{123}X_4, H_{124}X_3, H_{134}X_2, H_{234}X_1 (coefficients J_e(a,b,a), J_e(a,b,b), J_e(a,a,b), J_e(b,a,b); contractions T_{ee}X_·(t_e or t̄_e)): each ≤ β^{3/2} Σ_e T_{ee}ρ_e ≤ β^{3/2}‖T‖²_HS via (2.22) and Cauchy–Schwarz over e. Total 4β^{3/2}‖T‖²_HS.
  4. H_{1234} (coefficient J_e(a,b,a,b); contraction T_{ee}²): ≤ β Σ_e T_{ee}² ≤ β‖T‖²_HS.
  
  Sum: K_rot(β) = 6β² + 4β^{3/2} + β, eq. (2.9).
- **Cross terms** (not part of D⁴W): Σ_{i,j} R(S_{ij}, S̄_{ij}) from (2.15), bounded by (2.17) using gradient invariance (2.16) and (2.6); this yields the β Σ min{…} term of (2.8).
- **What is not used.** No positivity of R; no lower bound on R; the term DW[D_{1234}B] and the pairing 13|24 vanish identically, by (2.12) and same-generator commutation.

### 2.5 Re-derivation path for O(3), condensed

1. Define H, W by (2.1); X, H_k, J_e by (2.2)–(2.3); R by (2.4).
2. Abelian gauge (Lemma 2.1): W(t(df)g) = 1 ⇒ p. 5 display ⇒ (2.5); differentiate (2.1) twice ⇒ equality in (2.6); |H_2| ≤ β‖h‖² ⇒ inequality in (2.6); W((h+df)c) = W(hc) ⇒ R(h,df) = R(df,h) = 0.
3. Non-abelian gauge: (2.10) with B from (2.11); mixed second derivative S(f,g)c; (2.12) by the two half-turn conjugations.
4. Choose (2.13); Faà di Bruno on W∘B leaves pairings 12|34 and 14|23 ⇒ (2.14); sum ⇒ (2.15).
5. Upper bound for Σ R(S_{ij}, S̄_{ij}): (2.16) + gradient invariance + (2.6) ⇒ (2.17).
6. Lower bound for Σ D⁴W: expand (2.19); keep the square (2.20); rewrite rows of T as gradients (2.21); Lemma 2.1 ⇒ (2.22); bound the eleven table terms ⇒ (2.18) with (2.9). Realness of the sum: conjugation swaps i and j and permutes the arguments of the symmetric D⁴W.
7. Insert (2.17) and (2.18) into (2.15) ⇒ (2.8).

### 2.6 What §2 uses from earlier sections

Only the model (1.1) specialized to n = 3, with the generalization to an arbitrary finite graph and arbitrary pinned spins on P. Nothing from §1 beyond notation. The graph is arbitrary: no lattice structure, planarity or dimension enters §2.

### 2.7 What is specific to the model or dimension in §2

- Target S² and group SO(3): the change of variables s ↦ exp(−F_x)s preserves the uniform measure on S² (used in Lemma 2.1 and (2.10)); [a,b] = c (used in (2.11), (2.14)); the half-turn structure of so(3) (used in (2.12)); unit spins and unit-norm generators (used in |J_e| ≤ β, (2.2)).
- Nearest-neighbor pair interaction b_e s_x^T s_y: the twisted energy (2.1) and the local expansion (2.3) use that the energy is a sum of edge terms bilinear in the spins.
- No dimension dependence: the statement is for every finite graph.

### 2.8 GAPs in §2

- GAP (minor, not reached in detail): the paper asserts (p. 5) that the identity (2.10) holds "for F sufficiently small" and then differentiates at 0 and extends complex-multilinearly. I accepted that the real-analytic identity in real potentials passes to the complex-multilinear derivative identity (2.14) without writing out the argument; the paper states "No complex change of spin variables is involved in this extension" (p. 6).
- I checked the table rows H_{13}X_2X_4 (worked in the paper), H_{23}-type pairing signs, H_{123}X_4 and H_{1234} by hand. I did not recompute the remaining seven table rows individually; they follow "by the same finite summation" (p. 7).

## 3. §3 — A uniform estimate for an annular twist (pp. 8–12)

### 3.1 Geometry and pinning (p. 8)

`build/sections/annulus.tex` lines 12–29:

```latex
For an integer $N\geq 2$, let
\[
 D_N=\{x\in\Z^2:N\leq\|x\|_\infty\leq 2N\}.
\]
Let $G=(V,E)$ be any subgraph of the nearest-neighbor graph on $D_N$,
and orient each edge in a positive coordinate direction.  Give its edges
strengths $0\leq b_e\leq\beta$.  Put $s_*=(1,0,0)$ and prescribe
\[
 s_x=\exp(tc)s_*\quad(\|x\|_\infty=N),
 \qquad
 s_x=s_*\quad(\|x\|_\infty=2N)
 \quad (x\in V).
\]
All other spins are integrated against uniform probability measure on
$S^2$.  Write $Z(t)$ for the resulting partition integral, including the
interactions between any prescribed spins.  Missing vertices on either
boundary require no additional convention: only vertices of $V$ are
pinned.
```

So: the square annulus D_N = {N ≤ ‖x‖_∞ ≤ 2N} ⊂ Z², N ≥ 2, any subgraph G of its nearest-neighbor graph, edges oriented in a positive coordinate direction, strengths in [0,β]. Pinned: every vertex of V on the inner layer ‖x‖_∞ = N is set to exp(tc)s_* (s_* = (1,0,0) rotated by angle t about the third axis), every vertex of V on the outer layer ‖x‖_∞ = 2N to s_*. Z(t) includes interactions among pinned spins. Z is 2π-periodic in t.

### 3.2 Statements

Proposition 3.1 (p. 8), quoted, containing (3.1) and (3.2):

`build/sections/annulus.tex` lines 31–44:

```latex
\begin{proposition}[Vanishing annular twist cost]\label{prop:annulus}
For every finite $\beta>0$ there is a finite constant $C_\beta$ such that,
for all $N\geq2$ and all graphs and strengths just described,
\begin{equation}\label{eq:annulus-minmax}
 \frac{\min_{t\in\R}Z(t)}{\max_{t\in\R}Z(t)}
 \geq 1-\frac{C_\beta}{\log N}.
\end{equation}
In particular,
\begin{equation}\label{eq:annulus-ratio}
 \frac{Z(\pi)}{Z(0)}\geq1-\frac{C_\beta}{\log N}.
\end{equation}
The constant is independent of the subgraph, the strength array, and
the twist parameter.
\end{proposition}
```

Lemma 3.2 (p. 10), quoted, containing (3.8), (3.9), (3.10):

`build/sections/annulus.tex` lines 158–186:

```latex
\begin{lemma}[The multiscale family]\label{lem:modes}
There is a constant $C$, depending only on the fixed profiles, with the
following properties for every $N\geq2$ and every such subgraph.
First,
\begin{equation}\label{eq:annulus-commutator}
 \sum_i S(f_i,\overline{f_i})=-i\lambda_N\widetilde h,
 \qquad
 \lambda_N
 =N\!\!\sum_{\substack{1\leq r_p\leq N\\p_1>0}}
      r_p^{-3}\sin\left(\frac{2\pi p_1}{16N}\right)
 \geq\frac12\sum_{r=1}^N\frac1r.
\end{equation}
Second,
\begin{equation}\label{eq:annulus-midpoint}
 \|h\|_2\leq C,\qquad
 \|h+\widetilde h\|_\infty\leq\frac CN,\qquad
 \|h-\widetilde h\|_1\leq\frac CN.
\end{equation}
Finally,
\begin{equation}\label{eq:annulus-test-cost}
 \|T\|_{\mathrm{HS}}^2+
 \sum_{i,j}\min\left\{
   \|f_i\|_\infty^2\|df_j\|_2^2,
   \|f_j\|_\infty^2\|df_i\|_2^2
 \right\}
 \leq C(1+\log N).
\end{equation}
All edge norms use counting measure on $E$.
\end{lemma}
```

### 3.3 Mechanism

**Boundary gauge (3.3)–(3.4) (p. 9).** The twist is realized as a gauge field: Θ is a fixed smooth bump equal to 1 near [−1,1]² and supported in (−2,2)², and h = d[Θ(·/N)].

`build/sections/annulus.tex` lines 48–53:

```latex
Fix once and for all a real smooth function $\Theta$ equal to one in a
neighborhood of $[-1,1]^2$ and compactly supported in $(-2,2)^2$.  Define
the real edge form
\begin{equation}\label{eq:annulus-h}
 h_e=\Theta(y/N)-\Theta(x/N),\qquad e=(x,y).
\end{equation}
```

`build/sections/annulus.tex` lines 54–70:

```latex
Use a subscript $t$ on $W,X,H_k,R$ for the quantities from the preceding
section evaluated with these boundary values, and write
$\langle\,\cdot\,\rangle_t$ for the corresponding expectation.  Although
$\Theta$ does not vanish on the inner pinned layer, its gauge action is
still explicit:
\begin{equation}\label{eq:annulus-boundary-gauge}
 W_t(shc)=\frac{Z(t+s)}{Z(t)},
 \qquad
 R_t(h,h)=-\frac{Z''(t)}{Z(t)}.
\end{equation}
Indeed, substituting
$s_x^{\mathrm{old}}=\exp(-s\Theta(x/N)c)s_x^{\mathrm{new}}$
cancels every edge matrix.  It changes the inner boundary value to
$\exp((t+s)c)s_*$ and leaves the outer value unchanged.  Thus the gauge
changes the prescribed boundary values, rather than treating them as
integrated variables.

```

In my words: h is a gradient, but of a function that equals 1 on the inner pinned layer and 0 on the outer one, so it is not gauge-trivial in the sense of Lemma 2.1; instead the gauge transformation moves the inner boundary angle from t to t + s. So R_t(h,h) is minus the logarithmic second derivative Z''/Z, and an upper bound δ on R_t(h,h) uniform in t is a lower bound Z'' ≥ −δZ.

The plan (p. 9, quoted):

`build/sections/annulus.tex` lines 71–80:

```latex
The analytic target is now an upper bound of order $1/\log N$ for
$R_t(h,h)$.  We will construct functions whose collective form
$\sum_i S(f_i,\overline{f_i})$ equals $-i\lambda_N\widetilde h$,
where $\lambda_N\geq\tfrac12\log N$ and $\widetilde h$ approximates $h$.
The two family costs in Proposition~\ref{prop:fourth} will total
$O(\log N)$, whereas bilinearity multiplies
$R_t(\widetilde h,\widetilde h)$ by $\lambda_N^2$.
This yields the required logarithmic gain; after constructing the family,
we will transfer the bound from $\widetilde h$ to $h$ without assuming
that $R_t$ is positive.
```

**Profiles (3.5) (p. 9).** g_m = ∂_mΘ; χ_m a smooth cutoff supported in the open annulus 𝒜 = {1 < ‖z‖_∞ < 2} and equal to 1 on supp g_m.

`build/sections/annulus.tex` lines 82–96:

```latex
For $m=1,2$, put $g_m=\partial_m\Theta$.  Choose a real smooth cutoff
$\chi_m$ compactly supported in
\[
 \mathcal A=\{z\in\R^2:1<\|z\|_\infty<2\}
\]
and equal to one on $\operatorname{supp}g_m$.  Such a cutoff exists
because $g_m$ vanishes near both boundary squares.  The real profiles
\begin{equation}\label{eq:annulus-profiles}
 \phi_{m,+}=\frac{\chi_m+g_m}{2},\qquad
 \phi_{m,-}=\frac{\chi_m-g_m}{2}
 \quad\hbox{satisfy}\quad
 \phi_{m,+}^2-\phi_{m,-}^2=g_m.
\end{equation}
This construction uses no square root of a signed smooth function.
All four profiles have compact support in $\mathcal A$.
```

**Test-function family (3.6) (p. 10).**

`build/sections/annulus.tex` lines 127–156:

```latex
For $p\in\Z^2$ let $r_p=\|p\|_\infty$.  Introduce four groups of
vertex functions, restricted to $V$:
\begin{equation}\label{eq:annulus-modes}
 f_{m,\pm,p}(x)
 =r_p^{-3/2}\phi_{m,\pm}(x/N)
       \exp\left(\frac{2\pi i\,p\cdot x}{16N}\right),
 \qquad
 1\leq r_p\leq N,\quad \pm p_m>0.
\end{equation}
They vanish at every pinned vertex.  The power $r_p^{-3/2}$ is chosen
to produce a harmonic sum: the shell face $p_m=\pm r$ contains $2r+1$ modes, each with
squared normalization $r^{-3}$ and coordinate-$m$ phase increment of order
$r/N$.  Its summed frequency factor in the commutator is therefore of
order $1/(Nr)$.  The opposite frequency signs give the adjacent-point
products of $\phi_{m,+}$ and $\phi_{m,-}$ opposite coefficients.
Index this finite family by $i$, write $u_i=df_i$, and, as in
Proposition~\ref{prop:fourth}, put
\[
 T_{e\ell}=\sum_i u_i(e)\overline{u_i(\ell)},
 \qquad
 S(f,g)_e=\tfrac12(f_xg_y-g_xf_y).
\]
For an edge $e=(x,x+e_m)$ define
\begin{equation}\label{eq:annulus-htilde}
 \widetilde h_e=\frac1N
 \left[
 \phi_{m,+}(x/N)\phi_{m,+}((x+e_m)/N)
 -\phi_{m,-}(x/N)\phi_{m,-}((x+e_m)/N)
 \right].
\end{equation}
```

So the family is four groups (m ∈ {1,2}, sign ±) of plane waves of frequency p/(16N) with 1 ≤ ‖p‖_∞ ≤ N and ±p_m > 0, amplitude r_p^{−3/2}, modulated by the profile φ_{m,±}(x/N). All vanish on both pinned layers because the profiles are compactly supported in 𝒜. h̃ (3.7) is the discrete version of (φ_{m,+}² − φ_{m,−}²)/N = g_m/N ≈ h on direction-m edges.

**Where log N comes from (gain).** One mode's commutator, eq. (3.11) (p. 10):

`build/sections/annulus.tex` lines 188–209:

```latex
\begin{proof}
For a real profile $\phi$ and one mode on shell $r_p$, direct substitution
on a positive coordinate-$j$ edge gives
\begin{equation}\label{eq:annulus-one-mode}
 S(f_p,\overline{f_p})_e
 =-i r_p^{-3}\phi(x/N)\phi((x+e_j)/N)
       \sin\left(\frac{2\pi p_j}{16N}\right).
\end{equation}
If the group is indexed by $m\ne j$, its frequency set is invariant
under $p_j\mapsto-p_j$, so its sum vanishes.  For $m=j$, the negative
half-plane gives the negative of the positive-half-plane sine sum.
Interchanging coordinates identifies the positive sums for $j=1,2$.
This proves the identity in \eqref{eq:annulus-commutator}.

On shell $r$, retain only the $2r+1$ frequencies with $p_1=r$ and
$|p_2|\leq r$.  Since
\[
 \sin\left(\frac{\pi r}{8N}\right)\geq\frac r{4N}
 \qquad(1\leq r\leq N),
\]
their contribution to $\lambda_N$ is at least
$(2r+1)/(4r^2)\geq1/(2r)$.  This gives its asserted lower bound.
```

So Σ_i S(f_i, f̄_i) = −iλ_N h̃ with λ_N given in (3.8): on a direction-j edge, groups with m ≠ j cancel by the symmetry p_j ↦ −p_j; for m = j, the + and − groups give φ_{j,+}φ_{j,+} and −φ_{j,−}φ_{j,−} times the same positive sine sum. Retaining on shell r only the 2r+1 frequencies with p_1 = r and using sin(πr/(8N)) ≥ r/(4N), shell r contributes at least N·(2r+1)·r^{−3}·r/(4N) = (2r+1)/(4r²) ≥ 1/(2r) to λ_N. So λ_N ≥ ½ Σ_{r=1}^N 1/r ≥ ½ log N: the log N in the gain is a harmonic sum over shells, produced by the choice r_p^{−3/2} (p. 10: "The power r_p^{-3/2} is chosen to produce a harmonic sum").

**Where (log N)² comes from.** Prop. 2.2 is applied with Q = Σ_i S(f_i, f̄_i) = −iλ_N h̃; by bilinearity and reality of h̃, R_t(Q, Q̄) = λ_N² R_t(h̃, h̃) (p. 12). The left side of (2.8) therefore carries λ_N² ≥ ¼(log N)².

**Approximation h̃ ≈ h, eq. (3.9) (pp. 10–11).** Centered Taylor expansion at z = (x+y)/(2N): h_e = N^{−1}∂_mΘ(z) + O(N^{−3}) and φ(x/N)φ(y/N) = φ(z)² + O(N^{−2}), so with (3.5) |h_e − h̃_e| ≤ CN^{−3}. With O(N²) edges and |h_e| + |h̃_e| ≤ C/N this gives ‖h‖_2 ≤ C, ‖h + h̃‖_∞ ≤ C/N, ‖h − h̃‖_1 ≤ C/N (p. 11).

**Where log N comes from (cost), eq. (3.10).** Mode norms (3.12) and the min-sum (3.13) (p. 11):

`build/sections/annulus.tex` lines 222–240:

```latex
For a mode on shell $r$, the product rule for a discrete difference
and $|e^{iv}-1|\leq|v|$ give
\begin{equation}\label{eq:annulus-mode-norms}
 \|f_i\|_\infty\leq Cr^{-3/2},\qquad
 \|df_i\|_2\leq Cr^{-3/2}(1+r)\leq C' r^{-1/2}.
\end{equation}
There are $O(r)$ modes on shell $r$, including all four groups.
For shells $s\leq r$, choose the bound in the minimum that uses the
supremum norm of the shell-$r$ mode.  The cost of a pair is at most
$Cr^{-3}s^{-1}$; all $O(rs)$ pairs therefore cost at most $Cr^{-2}$.
Summing both shell orderings gives
\begin{equation}\label{eq:annulus-minimum-sum}
 \sum_{i,j}\min\left\{
   \|f_i\|_\infty^2\|df_j\|_2^2,
   \|f_j\|_\infty^2\|df_i\|_2^2
 \right\}
 \leq C\sum_{r=1}^N\sum_{s=1}^r r^{-2}
 \leq C'(1+\log N).
\end{equation}
```

The Hilbert–Schmidt bound uses Fourier orthogonality on the torus (Z/16NZ)² rather than the triangle inequality, which would give only O(N) (p. 11):

`build/sections/annulus.tex` lines 242–300:

```latex
It remains to bound $T$.  Applying the triangle inequality mode by mode
would give only $\|T\|_{\mathrm{HS}}\leq\sum_i\|df_i\|_2^2=O(N)$;
we use Fourier orthogonality for the sharper bound.  Work first in a
single profile group and let $U$ be the matrix whose $p$-th column is
$r_p^{1/2}df_p$.
Embed the annulus in the discrete torus
$Q_N=(\Z/(16N)\Z)^2$.  Its representatives may be chosen in
$[-8N,8N-1]^2$; all profiles vanish near the boundary of this square,
so they also define torus functions.  Let $F$ have columns
\[
 F_{x,p}=\exp\left(\frac{2\pi i\,p\cdot x}{16N}\right),
 \qquad x\in Q_N,
\]
using the frequencies of this group.  They are distinct modulo $16N$,
and orthogonality of these unnormalized waves gives
\begin{equation}\label{eq:annulus-fourier-norm}
 F^*F=(16N)^2I,\qquad \|F\|_{\mathrm{op}}=16N.
\end{equation}
For direction $j$, let $P_j$ restrict functions on $Q_N$ to the
origins of the coordinate-$j$ edges in $E$.  Let $M_v$ denote
multiplication by $v$ on $Q_N$, and set
\[
 \Delta_j\phi(x)=\phi((x+e_j)/N)-\phi(x/N),
 \qquad \phi^{(j)}(x)=\phi((x+e_j)/N).
\]
The directional block of $U$ has the exact factorization
\begin{align}\label{eq:annulus-operator-factorization}
 U_j={}&P_jM_{\Delta_j\phi}F
              \operatorname{diag}(r_p^{-1})\notag\\
 &+P_jM_{\phi^{(j)}}F
   \operatorname{diag}\left(
      \frac{\exp(2\pi i p_j/(16N))-1}{r_p}\right).
\end{align}
Here $\|P_j\|_{\mathrm{op}}\leq1$,
$\|\Delta_j\phi\|_\infty\leq C/N$, and the last diagonal factor
has operator norm at most $C/N$, since $|p_j|\leq r_p$.
Equation~\eqref{eq:annulus-fourier-norm} now gives
$\|U_j\|_{\mathrm{op}}\leq C$ and hence
$\|U\|_{\mathrm{op}}\leq C$ after combining the two directions.

The contribution of this group is
\[
 T_{\mathrm{group}}=U\operatorname{diag}(r_p^{-1})U^*.
\]
The finite-matrix inequality
$\|ADB\|_{\mathrm{HS}}\leq
\|A\|_{\mathrm{op}}\|D\|_{\mathrm{HS}}\|B\|_{\mathrm{op}}$
therefore implies
\[
 \|T_{\mathrm{group}}\|_{\mathrm{HS}}^2
 \leq C\sum_p r_p^{-2}
 \leq C'\sum_{r=1}^N r\,r^{-2}
 \leq C''(1+\log N).
\]
The triangle inequality over the four groups gives the same order
for $\|T\|_{\mathrm{HS}}^2$.  Together with
\eqref{eq:annulus-minimum-sum}, this proves
\eqref{eq:annulus-test-cost}.
\end{proof}
```

So both costs in (2.8) are sums Σ_r (number of modes on shell r ≈ r) × r^{−2}, i.e. harmonic sums ≲ 1 + log N.

**Proof of Prop. 3.1 (p. 12).** Apply Prop. 2.2 at twist t with the family (3.6):

`build/sections/annulus.tex` lines 305–317:

```latex
Apply Proposition~\ref{prop:fourth} in the law with twist $t$, using
the family \eqref{eq:annulus-modes}.  The left side is
\[
 R_t(-i\lambda_N\widetilde h,i\lambda_N\widetilde h)
 =\lambda_N^2R_t(\widetilde h,\widetilde h),
\]
because $R_t$ is bilinear and $\widetilde h$ is real.
Lemma~\ref{lem:modes} gives, uniformly in $t$,
\begin{equation}\label{eq:annulus-tested-response}
 R_t(\widetilde h,\widetilde h)
 \leq\frac{C_\beta(1+\log N)}{(\log N)^2}
 \leq\frac{C'_\beta}{\log N}.
\end{equation}
```

(Here, in my words: λ_N² R_t(h̃,h̃) ≤ (K_rot + β)·C(1+log N) by (2.8) and (3.10), and λ_N ≥ ½ log N.)

Transfer from h̃ to h without positivity of R_t, eqs. (3.17)–(3.20):

`build/sections/annulus.tex` lines 319–357:

```latex
We must transfer this upper bound to $h$ without treating $R_t$ as a
positive form.  If $R_t(h,h)\leq0$, the desired positive upper bound
already holds.  Otherwise, the defining second-derivative identity
and the pointwise bound on $H_{2,t}$ give
\begin{equation}\label{eq:annulus-current-control}
 \langle X_t(hc)^2\rangle_t
 <-\langle H_{2,t}(hc,hc)\rangle_t
 \leq\beta\|h\|_2^2\leq C\beta.
\end{equation}
Put $d=h-\widetilde h$.  The current is a real random variable on
real edge forms, and
\[
 |X_t(dc)|\leq\beta\|d\|_1\leq\frac{C\beta}{N}.
\]
Consequently ordinary Cauchy--Schwarz for random variables, followed
by \eqref{eq:annulus-current-control}, yields
\begin{align}\label{eq:annulus-current-error}
 \left|\langle X_t(hc)^2-X_t(\widetilde h c)^2\rangle_t\right|
 &\leq
 2\langle X_t(hc)^2\rangle_t^{1/2}
       \|X_t(dc)\|_\infty+\|X_t(dc)\|_\infty^2\notag\\
 &\leq\frac{C_\beta}{N}.
\end{align}
The local quadratic term satisfies the stronger pointwise bound
\begin{align}\label{eq:annulus-local-error}
 |H_{2,t}(hc,hc)-H_{2,t}(\widetilde h c,\widetilde h c)|
 &\leq\beta\sum_{e\in E}|d_e|\,|h_e+\widetilde h_e|\notag\\
 &\leq\frac{C\beta}{N^2}.
\end{align}
Equations~\eqref{eq:annulus-tested-response}--\eqref{eq:annulus-local-error}
and the formula for $R_t$ imply in both cases that
\begin{equation}\label{eq:annulus-ode}
 -\frac{Z''(t)}{Z(t)}=R_t(h,h)
 \leq\frac{C_\beta}{\log N}+\frac{C_\beta}{N}
 \leq\delta_N,\qquad
 \delta_N=\frac{C'_\beta}{\log N}.
\end{equation}
All constants are independent of the boundary twist, the retained
vertices and edges, and the admissible strengths.
```

Integration to the min/max ratio (p. 12):

`build/sections/annulus.tex` lines 359–375:

```latex
Finally, $Z$ is positive, smooth, and $2\pi$-periodic.  Set
$M=\max_t Z(t)$ and choose $t_0$ with $Z(t_0)=M$.  Then
$Z'(t_0)=0$, and \eqref{eq:annulus-ode} implies
$Z''(t)\geq-\delta_N Z(t)\geq-\delta_N M$.
For any $t$, choose an orientation $\epsilon\in\{-1,1\}$ and
$0\leq\ell\leq\pi$ with $t=t_0+\epsilon\ell$ modulo $2\pi$.
Twice integrating the second-derivative bound for
$v\mapsto Z(t_0+\epsilon v)$ gives
\[
 Z(t)\geq M-\frac{\delta_N M\ell^2}{2}
 \geq M\left(1-\frac{\pi^2\delta_N}{2}\right).
\]
Taking the minimum and enlarging $C_\beta$ proves
\eqref{eq:annulus-minmax}.  Since
$Z(\pi)/Z(0)\geq\min Z/\max Z$, it also proves
\eqref{eq:annulus-ratio}.
\end{proof}
```

### 3.4 Packet item 4, summarized

- Geometry: square annulus D_N = {x ∈ Z²: N ≤ ‖x‖_∞ ≤ 2N}, N ≥ 2, any subgraph, any strengths in [0,β]; outer side 4N.
- Pinned: only vertices of V; inner layer ‖x‖_∞ = N at exp(tc)s_*, outer layer ‖x‖_∞ = 2N at s_*; s_* = (1,0,0).
- Test family: (3.6), f_{m,±,p}(x) = r_p^{−3/2} φ_{m,±}(x/N) exp(2πi p·x/(16N)), 1 ≤ r_p ≤ N, ±p_m > 0, with profiles (3.5).
- log N: in the gain, λ_N ≥ ½ Σ_{r≤N} 1/r (3.8); in the cost, ‖T‖²_HS ≤ C Σ_p r_p^{−2} and the min-sum ≤ C Σ_r Σ_{s≤r} r^{−2}, both ≤ C(1+log N) (3.10).
- (log N)²: λ_N² multiplying R_t(h̃,h̃) by bilinearity, left side of (2.8).
- Final inequality: Z(π)/Z(0) ≥ 1 − C_β/log N, (3.2), from min Z/max Z ≥ 1 − C_β/log N, (3.1), via −Z''/Z = R_t(h,h) ≤ δ_N = C'_β/log N (3.20).

### 3.5 What §3 uses from earlier sections

Prop. 2.2 (2.8)–(2.9), applied with P = the two pinned layers; the definitions (2.1)–(2.4) at fixed twist t; the equality in (2.6) and the pointwise bound |H_2(hc,hc)| ≤ β‖h‖²₂ from (2.2)–(2.3) (used in (3.17), (3.19)); the bound |X_t(dc)| ≤ β‖d‖_1 from (2.2)–(2.3) with k = 1.

### 3.6 What is specific to the model or dimension in §3

- Dimension 2: the counts of modes per shell (2r+1 per face, O(r) per shell) in (3.8), (3.13) and in Σ_p r_p^{−2}; the O(N²) edge count in (3.9); the torus (Z/16NZ)² in (3.14). The exponent −3/2 in (3.6) is chosen so that these 2D counts give harmonic sums (p. 10).
- Spin model: rotation about the third axis c moves the inner boundary (3.4); 2π-periodicity of Z (p. 12); exp(πc)s_* = −s_* makes t = π the antipodal twist used in §4.

### 3.7 GAPs in §3

- GAP (minor): the Taylor estimates behind (3.9) are stated in two lines on p. 11 ("Centered Taylor expansion, with uniform bounds on the fixed smooth functions"). I verified the orders (the first-order terms cancel in the centered expansion) but not explicit constants.
- No other step in §3 that I could not follow.

## 4. §4 — Sign clusters and exponential decay (pp. 13–18)

Until §4.4, n = 3 and spins are written s_x; the law is the free law (1.1) on an arbitrary finite graph with b_e ∈ [0,β] (p. 13).

### 4.1 The amplitude/sign/angle decomposition (p. 13)

`build/sections/clusters.tex` lines 15–36:

```latex
Write
\begin{equation}\label{eq:spin-disintegration}
 s_x=(r_x\tau_x,q_x\cos\theta_x,q_x\sin\theta_x),
 \qquad q_x=(1-r_x^2)^{1/2}.
\end{equation}
Under the reference measure, the variables $r_x$, $\tau_x$, and
$\theta_x$ are independent and uniform on $[0,1]$, $\{-1,1\}$, and
$\R/(2\pi\Z)$, respectively.  Indeed, the first coordinate of a uniform
point on $S^2$ is uniform on $[-1,1]$, independently of its azimuthal
angle.  Conditional on all the amplitudes $r=(r_x)_{x\in V}$, the signs
and angles are independent Gibbs systems with respective couplings
\begin{equation}\label{eq:amplitude-couplings}
 K_{xy}=b_{xy}r_xr_y,
 \qquad L_{xy}=b_{xy}q_xq_y.
\end{equation}
Their edge energies are $\tau_x\tau_y$ and
$\cos(\theta_x-\theta_y)$.  Let $Z_I(K)$ and $Z_P(L)$ denote their
partition functions, using uniform probability reference measures.
The amplitude density is therefore proportional to
\begin{equation}\label{eq:amplitude-density}
 \rho(r)=Z_I(K(r))Z_P(L(r)).
\end{equation}
```

`build/sections/clusters.tex` lines 38–61:

```latex
We now use the embedded Ising representation underlying Wolff's
reflection clusters~\cite{Wolff}, coupled to bonds as in
Edwards--Sokal~\cite{EdwardsSokal}. Given the spins, sample bonds
$\eta_e\in\{0,1\}$ independently. For $e=\{x,y\}$, set
\begin{equation}\label{eq:bond-parameter}
 p_e=1-e^{-2K_e},
 \qquad
 \mathbb P(\eta_e=1\mid s)
   =p_e\mathbf1_{\{\tau_x=\tau_y\}}.
\end{equation}
The elementary identity
\begin{equation}\label{eq:ising-bond-expansion}
 e^{K_e\tau_x\tau_y}
  =e^{K_e}\big[(1-p_e)+p_e\mathbf1_{\{\tau_x=\tau_y\}}\big]
\end{equation}
shows that, conditional on $r,\eta$, the signs are independent fair
signs on the open connected components. In particular, we obtain the
weighted connection identity of Campbell and Chayes~\cite[Corollary~4,
Equation~(7)]{CC}:
\begin{equation}\label{eq:spin-connection}
 \langle s_x^1s_y^1\rangle_{G,b}^{(3)}
 =\E\big[r_xr_y\mathbf1_{\{x\leftrightarrow y\}}\big],
\end{equation}
where the connection uses only open edges of $G$.
```

In my words: the first coordinate of a uniform point of S² is uniform on [−1,1], independent of its azimuth; write it as r_x τ_x with amplitude r_x ∈ [0,1] and sign τ_x. Given all amplitudes, s_x·s_y = r_xr_y τ_xτ_y + q_xq_y cos(θ_x − θ_y), so the signs form a free Ising model with couplings K and the angles a free plane rotator with couplings L, conditionally independent, and the amplitude marginal has density ∝ Z_I(K(r)) Z_P(L(r)) (4.3). Edwards–Sokal bonds on the Ising part (4.4)–(4.5) make signs i.i.d. fair per open cluster given (r, η), hence (4.6).

### 4.2 Statements (pp. 14–17)

Lemma 4.1 (p. 14), quoted:

`build/sections/clusters.tex` lines 81–88:

```latex
\begin{lemma}\label{lem:association}
Let $w$ be a positive log-supermodular density on a finite product of
compact intervals or finite chains, with respect to the product of
Lebesgue or counting measures.  In the interval case, suppose $w$ is
continuous.  The probability law proportional to $w$ is associated.
Moreover, its conditional law on all but one coordinate increases
stochastically as the remaining coordinate increases.
\end{lemma}
```

Lemma 4.2 (p. 14), quoted:

`build/sections/clusters.tex` lines 124–128:

```latex
\begin{lemma}\label{lem:edge-correlation-inequalities}
In either the free Ising or the free planar system on a finite graph
with nonnegative couplings, every edge energy has nonnegative mean,
and any two edge energies have nonnegative covariance.
\end{lemma}
```

Its doubled-covariance identity (4.7):

`build/sections/clusters.tex` lines 137–143:

```latex
For covariances, take two independent copies of the Gibbs system.
If $T_e$ and $T_f$ are two edge energies and primes denote the second
copy, then
\begin{equation}\label{eq:doubled-covariance}
 2\operatorname{Cov}(T_e,T_f)
 =\E\big[(T_e-T'_e)(T_f-T'_f)\big].
\end{equation}
```

Lemma 4.3 (p. 15), quoted:

`build/sections/clusters.tex` lines 189–191:

```latex
\begin{lemma}\label{lem:amplitude-association}
The amplitude law with density \eqref{eq:amplitude-density} is associated.
\end{lemma}
```

Lemma 4.4 (p. 15), quoted:

`build/sections/clusters.tex` lines 229–232:

```latex
\begin{lemma}\label{lem:joint-association}
The joint free law of $r$ and $\eta$ is associated, with both amplitudes
and bond indicators ordered coordinatewise.
\end{lemma}
```

Proposition 4.5 (p. 16), quoted, with the pinned law it refers to:

`build/sections/clusters.tex` lines 263–273:

```latex
Let $s_*=(1,0,0)$.  For any set $B\subseteq V$, the pinned spin law
means that $s_x=s_*$ for $x\in B$, while all other spins retain uniform
reference measures and the same edge interactions.  Bonds are sampled
by the same rule \eqref{eq:bond-parameter}.

\begin{proposition}[Pinning comparison]\label{prop:pinning}
For the three-component free spin law on any finite graph with
nonnegative strengths, and bonds sampled by \eqref{eq:bond-parameter},
pinning any set of spins to $s_*$ cannot decrease the probability of
any increasing bond event $A$.
\end{proposition}
```

Annulus and crossing definitions, and Proposition 4.6 (p. 17), quoted:

`build/sections/clusters.tex` lines 306–324:

```latex
We now specialize to finite subgraphs of the nearest-neighbor graph
on $\Z^2$, with arbitrary strengths in $[0,\beta]$.  For $z\in\Z^2$
and an integer $N\geq2$, put
\[
 D_N(z)=\{v\in\Z^2:N\leq\|v-z\|_\infty\leq2N\}.
\]
An \emph{annulus crossing} is an open path in the restriction of $G$
to $D_N(z)$, joining its layers of radii $N$ and $2N$.
Edges or vertices outside this annulus are not allowed in the path.

\begin{proposition}[Uniform crossing suppression]\label{prop:crossings}
For each $\beta<\infty$ there are numbers $p_N=p_N(\beta)>0$ tending
to zero as $N\to\infty$ with the following property.  For every finite
lattice subgraph, every strength array in $[0,\beta]$, and every list
of centers whose closed boxes of sup radius $2N$ are pairwise disjoint,
the free probability that all the corresponding annuli are crossed
is at most $p_N$ to the number of annuli.  The boxes need not lie
inside the finite graph.
\end{proposition}
```

Proposition 4.7 (p. 17), quoted:

`build/sections/clusters.tex` lines 376–386:

```latex
\begin{proposition}\label{prop:connection}
For every finite $\beta>0$ there are constants $A_0<\infty$ and
$m_0>0$, depending only on $\beta$, such that the three-component free
bond law on every finite square-lattice subgraph with strengths in
$[0,\beta]$ satisfies
\[
 \mathbb P(x\leftrightarrow y)
 \leq A_0\exp(-m_0\|x-y\|_2)
 \qquad(x,y\in V).
\]
\end{proposition}
```

### 4.3 Mechanism

**Lemma 4.1 (FKG-type association).** Induction on the number of coordinates: log-supermodularity gives a likelihood ratio w(u,t')/w(u,t) increasing in u for t' ≥ t, so conditional laws increase stochastically in the last coordinate; the covariance decomposition Cov(F,H) = E[Cov(F,H|t)] + Cov(E[F|t], E[H|t]) finishes (p. 14).

**Lemma 4.2 (Ginibre, Ising and plane-rotator only).** Mean: expand the Gibbs exponential; Ising sign moments are 0 or 1; products of cosines of angle differences have nonnegative constant Fourier coefficient. Covariance: duplicate the system, (4.7); for the plane rotator substitute θ = φ+ψ, θ' = φ−ψ (a surjective homomorphism of tori, preserving Haar probability), making the doubled weight exp(2Σ L_g cos φ_g cos ψ_g) and T_e − T'_e = −2 sin φ_e sin ψ_e, so the numerator is a sum of nonnegative coefficients times squares; for Ising use u = (d+d')/2, v = (d−d')/2 (pp. 14–15). The paper states its scope: "no such inequality for the full vector-spin law is assumed" (p. 14).

**Lemma 4.3 (amplitude association, = [CC, Lemma 2]).** Mixed partials of log Z(J(r)) in distinct sites:

`build/sections/clusters.tex` lines 193–222:

```latex
\begin{proof}
For either partition function, with couplings $J_e(r)$ and edge
energies $T_e$, differentiation in distinct site variables $r_x,r_y$
gives
\begin{align*}
 \partial_x\partial_y\log Z(J(r))
 &=\sum_e\langle T_e\rangle\,\partial_x\partial_y J_e\\
 &\quad+\sum_{e,f}\operatorname{Cov}(T_e,T_f)
              (\partial_xJ_e)(\partial_yJ_f).
\end{align*}
For $J=K$, the first derivatives are nonnegative.  For $J=L$, they are
nonpositive, so the products in the second sum are again nonnegative.
The distinct-site mixed derivatives of individual couplings vanish
unless $e=\{x,y\}$, in which case
\[
 \partial_x\partial_y K_{xy}=b_{xy},
 \qquad
 \partial_x\partial_y L_{xy}
   =b_{xy}\frac{r_xr_y}{q_xq_y}\geq0
\]
in the open amplitude cube.  Lemma~\ref{lem:edge-correlation-inequalities}
therefore shows that every distinct-site mixed derivative of
$\log\rho$ is nonnegative there.  An increment of $\log\rho$ in one
coordinate is thus nondecreasing in every other coordinate.  For two
points $u,v$, compare the successive increments from $u\wedge v$ to
$u$ with the same coordinate increments from $v$ to $u\vee v$.
Their sum gives the log-supermodular inequality in the open cube,
and continuity extends it to the closed cube.  The possible singularity
of derivatives at $r_x=1$ is immaterial.  The density itself is positive
and continuous, so Lemma~\ref{lem:association} applies.
```

So the sign structure needed is: ∂K/∂r_x ≥ 0 and ∂L/∂r_x ≤ 0 (since q = (1−r²)^{1/2} decreases), making the products of first derivatives in the covariance sum nonnegative in both factors Z_I and Z_P; plus nonnegative edge-energy means and covariances from Lemma 4.2; plus ∂_x∂_y K_{xy} = b_{xy} ≥ 0, ∂_x∂_y L_{xy} ≥ 0. The paper flags that "both its sign structure and its boundary-free hypothesis matter" (p. 15).

**Lemma 4.4 (joint association of (r, η)).** Given r, the bond weight is the q = 2 random-cluster weight

`build/sections/clusters.tex` lines 234–258:

```latex
\begin{proof}
Conditional on $r$, summing the signs in
\eqref{eq:ising-bond-expansion} gives a bond weight proportional to
\begin{equation}\label{eq:conditional-fk-weight}
 2^{k(\eta)}\prod_{e\in E}p_e^{\eta_e}(1-p_e)^{1-\eta_e},
\end{equation}
where $k(\eta)$ counts all open components, including singletons.
The function $k$ is supermodular: the change on opening an edge is
$-1$ if its endpoints were disconnected and $0$ otherwise, and this
change increases as other bonds are opened.  Hence the weight is
log-supermodular and its law is associated by
Lemma~\ref{lem:association}.  Increasing any nondegenerate parameter
$p_e$ tilts the weight by an increasing likelihood ratio.  Association
therefore implies stochastic increase of the conditional bond law
in each $p_e$, and consequently in $r$.  A zero parameter forces its
edge closed; restriction to the remaining edges, or continuity, handles
these degenerate cases.

For bounded increasing functions $F(r,\eta)$ and $H(r,\eta)$,
conditional association makes
$\E[\operatorname{Cov}(F,H\mid r)]$ nonnegative.  Their conditional
means are increasing in $r$, by stochastic monotonicity of the bond
law and their direct monotonicity in $r$.  Their covariance is
nonnegative by Lemma~\ref{lem:amplitude-association}.  The total
covariance identity proves the claim.
```

**Proposition 4.5 (pinning comparison).** Condition the free system on r_x ≥ 1 − ε, τ_x = 1 for x ∈ B; given (r, η) the sign requirement has probability 2^{−k_B(η)}, so the induced (r, η)-law is tilted by

`build/sections/clusters.tex` lines 275–303:

```latex
\begin{proof}
For $0<\varepsilon<1$, condition the free system on
$r_x\geq1-\varepsilon$ and $\tau_x=1$ for every $x\in B$.
Conditional on $r,\eta$, the probability of the sign requirement is
$2^{-k_B(\eta)}$, where $k_B$ is the number of open components meeting
$B$.  Thus the induced law of $r,\eta$ is tilted by
\[
 T_\varepsilon(r,\eta)
 =\mathbf1_{\{r_x\geq1-\varepsilon\ \text{for all }x\in B\}}
       2^{-k_B(\eta)}.
\]
Opening an edge cannot increase $k_B$; this also holds when it merges
a component meeting $B$ with one not meeting $B$.  Consequently
$T_\varepsilon$ is increasing.  Lemma~\ref{lem:joint-association} gives
\[
 \frac{\E[\mathbf1_A T_\varepsilon]}{\E[T_\varepsilon]}
 \geq\mathbb P(A).
\]

As $\varepsilon\downarrow0$, the conditioned reference measures on
the positive caps converge to point masses at $s_*$.  The finite-volume
Gibbs weight is bounded, positive, and continuous.  The spin-to-bond
kernel is also continuous, since its edge-opening probability is
\[
 1-\exp\big(-2b_{xy}\max\{s_x^1s_y^1,0\}\big).
\]
In particular there is no discontinuity at a zero first coordinate.
Integrating any event of the finitely many bonds and passing to the
limit proves the asserted comparison for exact pinning.
```

**The crossing identity, eq. (4.9) (p. 17) — packet item 5.**

`build/sections/clusters.tex` lines 326–350:

```latex
\begin{proof}
First work in one restricted annulus graph.  Pin both of its layers
to $s_*$, and let $Z(0)$ be its partition function.  Let $Z(\pi)$ be
the partition function with the inner layer changed to $-s_*$ and
the outer layer left at $s_*$.  Then
\begin{equation}\label{eq:crossing-partition-ratio}
 \mathbb P_{\mathrm{pinned}}(\text{crossing})
   =1-\frac{Z(\pi)}{Z(0)}.
\end{equation}
Indeed, in the expansion \eqref{eq:ising-bond-expansion}, every open
component missing the pinned layers has two sign choices.  A component
meeting a pinned layer has one sign choice, provided all prescribed
signs on it agree.  Changing the inner signs therefore removes exactly
the bond configurations joining the two layers and changes the weight
of none of the remaining configurations.  All other factors are
unchanged: the boundary amplitudes are $r=1,q=0$, and the reference
normalization for the unpinned signs is the same.  This argument holds
at every interior amplitude configuration and hence after integration.
If a layer has no vertices, both sides of
\eqref{eq:crossing-partition-ratio} are zero.

By Proposition~\ref{prop:annulus}, the right side of
\eqref{eq:crossing-partition-ratio} is bounded uniformly by a
nonnegative sequence tending to zero.  Enlarge it slightly to obtain
a strictly positive sequence $p_N$ with the same limit.
```

In my words: on the restricted annulus graph with both layers pinned to s_* (that is Z(0) of §3, t = 0) versus inner layer at −s_* = exp(πc)s_* (Z(π) of §3), pinned sites have r = 1, q = 0, so the plane-rotator part and the amplitude weights are identical in both cases; in the bond expansion (4.5) of the Ising part, flipping the inner signs kills exactly the bond configurations with an open path from the inner to the outer layer and leaves every other configuration's weight unchanged. Hence Z(π) = Z(0)·P_pinned(no crossing). Combined with (3.2): P_pinned(crossing) ≤ C_β/log N.

**Product over disjoint annuli (rest of Prop. 4.6 proof, p. 17).**

`build/sections/clusters.tex` lines 352–372:

```latex
For several annuli, pin all their present boundary-layer sites to
$s_*$ simultaneously.  The intersection of their crossing events is
increasing, so Proposition~\ref{prop:pinning} bounds its free
probability by its pinned probability.  In the pinned system, each
unpinned annulus region is separated from every vertex outside that
annulus by the pinned layers: a nearest-neighbor path cannot jump
across either layer.  The unpinned annulus spins are therefore
independent with precisely the required fixed boundary values.
For each resulting annulus marginal, apply
Equation~\eqref{eq:crossing-partition-ratio} to the partition functions
of that restricted annulus graph.
Boundary-only interactions contribute constants, and the independent
sampling of edge bonds preserves independence of the crossing events.
Their probabilities multiply, proving the bound.

This separation uses only edges present in $G$.  Missing vertices or
edges, including those beyond the finite volume, cannot create an
unseparated nearest-neighbor interaction.  Pinning the layer sites
that are present and applying Proposition~\ref{prop:annulus} to the
restricted subgraph proves the assertion also for truncated boxes.
\end{proof}
```

**Proposition 4.7 (coarse paths, pp. 17–18).**

`build/sections/clusters.tex` lines 388–437:

```latex
\begin{proof}
Translate the lattice so that $x=0$.  Choose a finite integer $N\geq2$
so large that, for the sequence in Proposition~\ref{prop:crossings},
\begin{equation}\label{eq:coarse-path-alpha}
 0<p:=p_N<1,
 \qquad \alpha:=4p^{1/25}<1.
\end{equation}
Partition $\Z^2$ into boxes
$z+\{0,\ldots,N-1\}^2$, with coarse vertices $z\in N\Z^2$.
Suppose $0$ and $y$ are joined by an open path and
$d:=\|y\|_\infty>4N$.  Recording the visited boxes and erasing loops
produces a simple nearest-neighbor coarse path from the origin's box
to the box containing $y$.  If its length is $k$, then
\begin{equation}\label{eq:coarse-path-length}
 k\geq d/N-1.
\end{equation}

Each vertex $z$ of this simple coarse path corresponds to a box
visited by the original open path, at some vertex strictly within
sup distance $N$ of $z$.  At least one of the original endpoints is
at distance greater than $2N$ from $z$: otherwise their mutual
distance would be at most $4N$.  A segment of the original path from
the visited vertex to that endpoint thus contains a crossing of
$D_N(z)$: take its first hit on the outer layer and its last preceding
hit on the inner layer. The intervening segment stays in the annulus.  It does not matter whether that segment survived coarse
loop erasure; its edges are still open.

Among the $k+1$ distinct coarse vertices, one of the 25 residue classes
modulo $5N$ in the two coordinates contains at least $(k+1)/25$
vertices.  Their closed outer boxes of sup radius $2N$ are disjoint,
because distinct centers in the class differ by at least $5N$ in
some coordinate.  Choose one maximizing class deterministically for
each candidate coarse path.  Proposition~\ref{prop:crossings} bounds
the probability that this path has all its required crossings by
$p^{(k+1)/25}$.  There are at most $4^k$ coarse paths of length $k$
from the initial box.  The union bound and
\eqref{eq:coarse-path-length} therefore give
\begin{align*}
 \mathbb P(0\leftrightarrow y)
 &\leq p^{1/25}\sum_{k\geq\lceil d/N-1\rceil}\alpha^k\\
 &\leq \frac{p^{1/25}}{\alpha(1-\alpha)}
           \exp\left(-\frac{-\log\alpha}{N}\,d\right).
\end{align*}
Since $d\geq\|y\|_2/\sqrt2$, take
\[
 m_0=\frac{-\log\alpha}{\sqrt2\,N}>0.
\]
Increasing the finite prefactor to cover $d\leq4N$ proves the result.
Every estimate was uniform in the finite subgraph and its strengths.
\end{proof}
```

### 4.4 Packet item 5: correlation inequalities and what each needs

- **Crossing identity:** P_pinned(crossing) = 1 − Z(π)/Z(0), eq. (4.9), for the restricted annulus graph with both layers pinned (Z(0): both at s_*; Z(π): inner at −s_*, outer at s_*). Needs: the Edwards–Sokal expansion (4.5) of the Ising sign part; pinned sites have r = 1, q = 0, so only the sign part sees the flip.
- **Ginibre** (Lemma 4.2): edge energies of the free Ising system and of the free plane-rotator system have nonnegative means and pairwise covariances. Needs nonnegative couplings K, L ≥ 0 (from b ≥ 0 and r, q ≥ 0), free boundary (no field), and the specific reference measures (Ising signs; the circle with its torus doubling and nonnegative Fourier coefficients of cos). It is used only on the two conditional systems, not on the S² model. Consumer: Lemma 4.3.
- **FKG** (Lemma 4.1, finite-product Fortuin–Kasteleyn–Ginibre criterion): needs a positive log-supermodular density (continuous in the interval case) on a finite product of compact intervals or finite chains. Consumers: Lemma 4.3 (amplitude density ρ(r), log-supermodular via Lemma 4.2 and the coupling sign structure); Lemma 4.4 (bond weight (4.8), log-supermodular because the component count k(η) is supermodular and p_e ∈ [0,1)).
- **Campbell–Chayes** [CC]: the paper cites (i) the weighted connection identity [CC, Corollary 4, Equation (7)] = (4.6), needing the Ising-sign embedding of the first coordinate; (ii) amplitude association [CC, Lemma 2] = Lemma 4.3, needing the sign structure ∂K ≥ 0, ∂L ≤ 0, Ginibre for both conditional systems, and the boundary-free hypothesis; (iii) bond association [CC, Corollary 1] → joint version Lemma 4.4, needing Lemma 4.3 plus FK monotonicity in p_e; (iv) boundary comparison, implemented here as Prop. 4.5 by pinning an arbitrary set simultaneously, needing Lemma 4.4 and the monotonicity of the tilt 1{r ≥ 1−ε on B} 2^{−k_B(η)}. The introduction (p. 3) lists "[CC, Lemma 2 and Corollaries 1, 3, and 4]" for amplitude association, bond association, boundary comparison, and connectivity; the assignment of Corollary 3 to boundary comparison is my inference from the order of that list.
- **Edwards–Sokal / Wolff** (not correlation inequalities, but used): bond coupling (4.4)–(4.5), needing K_e ≥ 0 so that p_e ∈ [0,1).

### 4.5 What §4 uses from earlier sections

Prop. 3.1, eq. (3.2), applied to the restricted annulus graph G ∩ D_N(z) (proof of Prop. 4.6, p. 17); the identification of §3's Z(π) with the inner layer at −s_* (exp(πc)s_* = −s_*). The model (1.1). Nothing from §2 directly.

### 4.6 What is specific to the model or dimension in §4

- n = 3 / S²: the decomposition (4.1)–(4.3) into an Ising sign system and a plane-rotator angle system, which is what makes Ginibre available (Lemma 4.2's scope statement, p. 14).
- Ferromagnetism: K, L ≥ 0 and p_e ∈ [0,1) throughout.
- Z²: annuli D_N(z), the separation property of pinned layers, the coarse lattice with 4^k paths and 25 residue classes, and ‖y‖_∞ ≥ ‖y‖_2/√2 (Props. 4.6–4.7). Lemmas 4.2–4.4 and Prop. 4.5 are stated for arbitrary finite graphs; Lemma 4.1 is a statement about finite product spaces.

### 4.7 GAPs in §4

- GAP (minor): the ε ↓ 0 limit in Prop. 4.5 (p. 16) is stated as convergence of conditioned cap measures to point masses plus continuity of the Gibbs weight and of the spin-to-bond kernel. I followed the outline but did not write out the weak-convergence argument.
- GAP (minor): in Lemma 4.1 the continuous-interval case of the induction (measurability and boundedness of the tilts, "bounded in the continuous case by positivity and compactness", p. 14) is accepted at the level stated.

## 5. Reduction from n to 3 components (§4.4, pp. 18–19)

### 5.1 Statement

Lemma 4.8 (p. 18), quoted:

`build/sections/clusters.tex` lines 447–455:

```latex
\begin{lemma}[Three-component conditional law]\label{lem:conditioning}
Let $n>3$, and let $\sigma_x\in S^{n-1}$ have the free Gibbs law
\eqref{eq:model} on a finite graph. Write
$\sigma_x=(\rho_x s_x,v_x)$, with $v_x\in\R^{n-3}$,
$\rho_x=\sqrt{1-|v_x|^2}$, and $s_x\in S^2$ whenever $\rho_x>0$.
Conditional on $(v_x)_{x\in V}$, the directions $(s_x)_{x\in V}$ have
exactly the free three-component Gibbs law with strengths
$b_{xy}\rho_x\rho_y\in[0,\beta]$.
\end{lemma}
```

### 5.2 Mechanism

`build/sections/clusters.tex` lines 457–479:

```latex
\begin{proof}
Represent a uniform sphere point by
$(U,V)/\sqrt{|U|^2+|V|^2}$, where $U\in\R^3$ and
$V\in\R^{n-3}$ have independent standard Gaussian coordinates.
The direction $U/|U|$ is uniform on $S^2$ and is independent of
$(|U|,V)$. Since $v=V/\sqrt{|U|^2+|V|^2}$ depends only on these
latter variables, the three-component direction is independent of
$v$ under the reference law. Applying this representation independently
at all sites proves that, conditional on all $v_x$, the reference
directions remain independent and uniform on $S^2$.

The conditioned interaction is
\[
 \sum_{\{x,y\}\in E}b_{xy}v_x\cdot v_y
 +\sum_{\{x,y\}\in E}b_{xy}\rho_x\rho_y s_x\cdot s_y.
\]
The first sum is constant under the conditional law and cancels from
its normalization. The second gives the asserted couplings. The event
$\rho_x=0$ has reference probability zero and still has Gibbs
probability zero because the finite Gibbs density is positive and
bounded. Alternatively, choosing any independent reference direction
at such a site would decouple it and leave all spin products unchanged.
\end{proof}
```

Proof of Theorem 1.1 (p. 19), quoted:

`build/sections/clusters.tex` lines 481–504:

```latex
\begin{proof}[Proof of Theorem~\ref{thm:main} and Corollary~\ref{cor:limits}]
For $n=3$, Equation~\eqref{eq:spin-connection} and
Proposition~\ref{prop:connection} give
\[
 0\leq\langle s_x^1s_y^1\rangle_{G,b}^{(3)}
 \leq A_0e^{-m_0\|x-y\|_2}.
\]
For $n>3$, apply Lemma~\ref{lem:conditioning}. Conditional on all
$v_x$, the first-coordinate correlation is $\rho_x\rho_y$ times
the corresponding three-component first-coordinate correlation.
The latter is nonnegative and satisfies the same bound, uniformly
in the conditional strength array. Since $0\leq\rho_x\rho_y\leq1$,
averaging gives, for every $n\ge3$,
\[
 0\leq\langle\sigma_x^1\sigma_y^1\rangle_{G,b}^{(n)}
 \leq A_0e^{-m_0\|x-y\|_2}.
\]
The original free spin law is invariant under coordinate permutations,
so
\[
 \langle\sigma_x\cdot\sigma_y\rangle_{G,b}^{(n)}
 =n\langle\sigma_x^1\sigma_y^1\rangle_{G,b}^{(n)}\geq0.
\]
Taking $A=nA_0$ and $m=m_0$ proves Theorem~\ref{thm:main}.
```

In my words: conditioning on the last n − 3 coordinates leaves an O(3) model with strengths b_{xy}ρ_xρ_y ∈ [0,β]; because Prop. 4.7 and (4.6) are uniform over all strength arrays in [0,β]^E, the bound holds conditionally for every realization and then on average, with the extra factor ρ_xρ_y ∈ [0,1]. Coordinate-permutation invariance turns the first-coordinate bound into the bound on σ_x·σ_y with A = nA_0, m = m_0.

### 5.3 Uses, specificity, GAPs

- Uses: (4.6), Prop. 4.7, uniformity in [0,β]^E.
- Specific: Gaussian representation of the uniform measure on S^{n−1} (the S² direction U/|U| is independent of (|U|, V)); the quadratic form σ_x·σ_y splits as v_x·v_y + ρ_xρ_y s_x·s_y; b ≥ 0 keeps the conditioned strengths in [0,β].
- GAP: none.

## 6. What is two-dimensional, what is spin-specific

The paper does not use planarity in the topological sense (no planar duality, dual circuits, or Jordan-curve argument appears in the source). The word "planar" in the source (clusters.tex lines 125, 133, 144) refers to the plane-rotator system of Lemma 4.2. The entries under 6.1 are therefore the steps that use the two-dimensional square lattice.

### 6.1 Steps that use two-dimensionality / the square lattice

- §2 (Lemma 2.1, Prop. 2.2, (2.1)–(2.22)): no use; the statements are for an arbitrary finite graph (rotation.tex line 17).
- (3.8), Lemma 3.2: the lower bound λ_N ≥ ½Σ_{r≤N} 1/r counts 2r+1 lattice frequencies on the face p_1 = r of shell r in Z²; with the normalization r_p^{−3/2} of (3.6) each shell contributes ≳ 1/r.
- (3.9), Lemma 3.2: ‖h‖_2 ≤ C uses O(N²) annulus edges with |h_e| ≤ C/N (the 2D scale invariance of the Dirichlet energy of Θ(·/N)); ‖h − h̃‖_1 ≤ C/N uses O(N²) edges with error CN^{−3}.
- (3.10) via (3.12)–(3.13): "There are O(r) modes on shell r" (p. 11) and the pair count O(rs) give Σ_r Σ_{s≤r} r^{−2} ≤ C(1+log N).
- (3.10) via (3.14)–(3.15): the torus Q_N = (Z/16NZ)², F*F = (16N)²I, and Σ_p r_p^{−2} ≤ C Σ_r r·r^{−2}.
- (3.16)–(3.20), Prop. 3.1: the balance gain λ_N² ≳ (log N)² against cost ≲ log N, giving 1/log N; both are 2D harmonic sums.
- Prop. 4.6 (p. 17): pinned layers of D_N(z) separate the annulus interior from the rest ("a nearest-neighbor path cannot jump across either layer"), giving conditional independence of different annuli. As written the argument uses only that sup-norm layers are vertex cut sets for nearest-neighbor paths.
- Prop. 4.7, (4.10)–(4.11): coarse lattice NZ² with at most 4^k coarse paths of length k, 25 = 5² residue classes mod 5N, α = 4p^{1/25}, and d ≥ ‖y‖_2/√2.
- Theorem 1.1 / Cor. 1.2: G is a finite subgraph of the square lattice only because §3 and Props. 4.6–4.7 need it.

### 6.2 Steps that use that the target is a sphere (and n ≥ 3)

- (2.1), Lemma 2.1, (2.10), (3.4): the changes of variables s_x ↦ exp(−F_x)s_x preserve the uniform probability measure on S² (rotation invariance; S² is an SO(3)-homogeneous space).
- (2.2): |J_e| ≤ β uses |s_x| = 1 and unit operator norm of products of a, b, c.
- (2.11), (2.14): existence of two noncommuting generators with [a,b] = c; the commutator direction c is the twist direction. This is where n ≥ 3 enters §2 (for S¹ there is no such pair).
- (2.12): the fixed spaces of the two half-turn conjugations on so(3) are ℝa and ℝb with zero intersection.
- (3.4), (3.20): Z(t) is 2π-periodic because exp(tc) is; (4.9) uses exp(πc)s_* = −s_*.
- (4.1)–(4.3): uniform measure on S² has first coordinate uniform on [−1,1] independent of azimuth, so that conditional on amplitudes the model splits into an Ising system (4.2, K) and a plane-rotator system (4.2, L).
- Lemma 4.2: Ginibre is proved for those two conditional systems (Ising; circle via the torus homomorphism θ = φ ± ψ), not for the S² model.
- Lemma 4.3: the sign structure ∂K/∂r ≥ 0, ∂L/∂r ≤ 0 comes from r² + q² = 1.
- (4.6), (4.9): spin correlation and boundary flip are expressed through the sign τ of the first coordinate.
- Lemma 4.8: Gaussian representation of the uniform measure on S^{n−1}, independence of U/|U| from (|U|, V), and the splitting σ_x·σ_y = v_x·v_y + ρ_xρ_y s_x·s_y.
- Proof of Thm 1.1 (p. 19): ⟨σ_x·σ_y⟩ = n⟨σ^1_xσ^1_y⟩ by invariance under coordinate permutations.
- §1.1 (p. 2): the exclusion of n = 2 is attributed to Fröhlich–Spencer's polynomial lower bound for the plane rotator.

### 6.3 Steps that use ferromagnetism (b_e ≥ 0)

- §2 and §3: the hypothesis is 0 ≤ b_e ≤ β (rotation.tex line 18; annulus.tex line 18). In the proofs of Lemma 2.1, Prop. 2.2 and Prop. 3.1, b_e enters only through the coefficient bound |J_e| ≤ β in (2.2) (and its consequences (2.5), (2.6), (2.22), (3.17), (3.19)). I did not find a step in §2 or §3 that uses the sign of b_e rather than |b_e| ≤ β.
- (4.4)–(4.5): p_e = 1 − e^{−2K_e} ∈ [0,1) requires K_e = b_e r_x r_y ≥ 0.
- Lemma 4.2: hypothesis "nonnegative couplings" for both Ising (K) and plane rotator (L).
- Lemma 4.3: ∂_x∂_y K_{xy} = b_{xy} ≥ 0 and ∂_x∂_y L_{xy} = b_{xy} r_xr_y/(q_xq_y) ≥ 0; products of first derivatives nonnegative.
- Lemma 4.4, (4.8): random-cluster weight with p_e ∈ [0,1).
- Props. 4.5, 4.6: rely on Lemma 4.4 association.
- (4.6) and Thm 1.1: ⟨s^1_x s^1_y⟩ = E[r_xr_y 1{x↔y}] ≥ 0 gives the lower bound 0 ≤ ⟨σ_x·σ_y⟩ in (1.2).
- Lemma 4.8: conditional strengths b_{xy}ρ_xρ_y ∈ [0,β] stay in the class over which Prop. 4.7 is uniform.

## 7. Consolidated GAP list

1. GAP (§2, (2.10)→(2.14)): passage from the real identity "for F sufficiently small" to the complex-multilinear fourth-derivative identity accepted at the level stated (pp. 5–6).
2. GAP (§2, table p. 7): seven of the eleven table rows not recomputed individually (rows H_{13}X_2X_4, H_{123}X_4, H_{1234} and the 14|23 pairing signs were checked).
3. GAP (§3, (3.9)): Taylor-expansion constants not written out (orders checked).
4. GAP (§4, Prop. 4.5): ε ↓ 0 weak-convergence step followed in outline only.
5. GAP (§4, Lemma 4.1): continuous-interval case of the induction accepted at the level stated.
