# Algebra behind the first corpus entries

These are re-derivations of known Srednicki results, not original physics.
Source statement locations and executable checks are in ../corpus.json.
The supplied handwritten LaTeX notes are incomplete in several later
sections, so the existing book transcriptions are secondary references, not
assumed error-free ground truth.

## 35.1: raise both spinor indices

With epsilon^{12}=1, let
E = [[0,1],[-1,0]]. The definition, with row index dotted, is

    barsigma^(mu,dota,a) = sum_(b,dotb) E[a,b] E[dota,dotb] sigma^(mu)[b,dotb].

Evaluating each of the four matrices gives (I,-sigma1,-sigma2,-sigma3).
The executable check sums the indices explicitly rather than guessing a
matrix order. It separately verifies

    sum_mu eta[mu] sigma^mu[a,dota] sigma^mu[b,dotb]
      = -2 epsilon_lower[a,b] epsilon_lower[dota,dotb].

For 35.2 substitute barsigma into
S_L^(mu,nu)=i(sigma^mu barsigma^nu-sigma^nu barsigma^mu)/4.
The spatial and boost results match (34.9)–(34.10), and all 256 components
of the Lorentz commutator are checked exactly.

## 47: ordered Clifford products and traces

Define h=-eta so {gamma(a),gamma(b)}=2h(a,b). The oriented swap is

    gamma(b) gamma(a) = -gamma(a) gamma(b) + 2h(a,b)  for b>a.

Repeated labels obey gamma(a)^2=h(a,a), with **no sum**.
Normal ordering reduces either the number of inversions or word length.

For an even word, cyclicity and the Clifford relation give

    Tr(gamma(a1)...gamma(a2m))
      = sum_(j=2..2m) (-1)^j h(a1,aj)
          Tr(gamma(a2)...omit gamma(aj)...gamma(a2m)),
    Tr(1)=4.

In particular,

    Tr(abcd)=4[h(a,b)h(c,d)-h(a,c)h(b,d)+h(a,d)h(b,c)].

Two algorithms are tested independently: oriented word reduction and trace
pairings, both against explicit Q(i) chiral matrices for every four-component
word. This is a useful finite representation check; it is not general
tensor/dummy-index canonicalisation.

Gamma5=i gamma0 gamma1 gamma2 gamma3 is only used in four dimensions.
For 47.3, form the 16 matrices 1, four gammas, six i[gamma,gamma]/4,
four gamma*gamma5 and gamma5. Their Hilbert–Schmidt Gram matrix is diagonal
with nonzero entries. Orthogonality proves independence; 16 independent
elements span the 16-dimensional complex vector space of 4x4 matrices.

47.2 is checked **only at d=4**. Do not infer its d=4-epsilon continuation.

## 50.1: completeness and slash sandwiches

Let phi and phi-star initially be independent commuting two-component
variables. Raising a spinor gives (phi1,-phi0).
Srednicki's definitions give

    p_(a,dota)=-phi_a phi-star_dota,
    |p] = (phi,0),       <p| = (0,phi-star),
    |p> = (0,raised phi-star), [p| = (raised phi,0).

The sum of the two outer products is minus the off-diagonal slash matrix:

    -slash(p) = outer(|p>, [p|) + outer(|p], <p|).

The source/tests construct these two outer products directly.
Substituting this decomposition for -slash(k) yields

    [r|(-slash(k))|p> = [rk]<kp>,
    <r|(-slash(k))|p] = <rk>[kp],

and the opposite-helicity sandwiches vanish. The tests construct slash
from p_mu gamma^mu and expand all components, not merely assert these formulas.

The same rank-one blocks show

    <pq>[qr]<rs>[sp] = Tr[(1-gamma5) slash(p)slash(q)slash(r)slash(s)/2],

the 50.4 four-slash identity. The chiral trace is expanded symbolically.

## 50.3: Schouten

For D(a,b)=a0*b1-a1*b0, expand

    D(a,b)D(c,d)+D(a,c)D(d,b)+D(a,d)D(b,c).

Each component monomial occurs twice with opposite sign, so the residual
polynomial is identically zero. A deliberately reversed middle sign leaves
a nonzero polynomial. Both Rust and live Cadabra perform this cancellation.

## 60.1: polarization identities, not just printed Ward checks

Define B_+^mu=<q|gamma^mu|k] and B_-^mu=[q|gamma^mu|k>.
The represented quantity is sqrt(2)*epsilon_+=-B_+/<qk>, and similarly for
minus with denominator [qk]. Exact contraction gives

    -p dot B_+ = <qp>[pk],
    -p dot B_- = [qp]<pk].

Setting p=k or p=q makes a bracket vanish, proving transversality, provided
the reference denominator is nonzero. Choosing q=k violates that assumption.

The pairwise B contractions produce twice the bracket numerators of
(60.38)–(60.40). Dividing by the two denominators and by two (from the
two sqrt(2)-scaled polarizations) recovers the stated formulas.

## 60.2a: make momentum conservation explicit

For every individual j, the same slash sandwich yields

    <ij>[jk] = -p_j dot (<i|gamma|k]),
    [ij]<jk> = -p_j dot ([i|gamma|k>).

Sum this equality over j. Linearity turns the right-hand side into the
same fixed bilinear contracted with sum_j p_j. Momentum conservation then
sets both sums to zero. The executable polynomial check expands a
five-leg instance of this factorization; the argument above works for any n.
No unsupported reduction modulo an implicit momentum-conservation ideal occurs.

## Parke–Taylor and SMGA: different evidence classes

The generic Parke–Taylor ratio checks cyclicity and reflection by exact
cross-multiplication, and little-group weights at exact rational witnesses.
These checks do not derive the formula from a Lagrangian.

The SMGA implementation computes the existing half-collinear (2,2)
sign-function formula and its separately implemented preamplitude recurrence.
At 200 conserved rational fixtures, n=3..10, both agree exactly.
They also agree with both old NumPy paths within 1e-10.
This is finite-case reproduction, **not an all-n proof**, an ordinary
Lorentzian BG calculation, or an extension of the cited result.
