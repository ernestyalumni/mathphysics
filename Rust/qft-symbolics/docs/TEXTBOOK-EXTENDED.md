# Further Srednicki checks — 2026-10-16

Source identifiers were checked against the existing book marker transcription:
50.2a at line9154, 50.4 at9163–9169; the chapter 47 problem list is at8678–8680.
Existing transcriptions contain malformed formulas. This work does not claim
fresh visual certification of those pages. All results here are known results.

## 47.2: formal dimension, rather than only a 4x4 specialization

Let W=gamma(a1)...gamma(ak), and use {gamma^mu,gamma^nu}=-2eta^mu nu.
Commuting the left contracted gamma successively through W gives

    gamma^mu W gamma_mu
       = (-1)^(k+1) d W
         + 2 sum_(j=1..k) (-1)^j
             gamma(a1)...omit gamma(aj)...gamma(ak) gamma(aj).

The first term uses gamma^mu gamma_mu=-d. Each remaining term is the
anticommutator contraction at a particular position. All indices displayed
as a1..ak remain free; this is an explicit specialized contraction operation,
not implicit summation inferred from gamma-label strings.

Using the existing ordered Clifford reducer with h=-eta, this proves

    gamma^mu gamma(a) gamma(b) gamma_mu
      = 4 eta(a,b) + (4-d) gamma(a) gamma(b),

    gamma^mu gamma(a) gamma(b) gamma(c) gamma_mu
      = 2 gamma(c) gamma(b) gamma(a)
        + (d-4) gamma(a) gamma(b) gamma(c).

Tests obtain formal d from the new tensor metric trace, then compare both
sides as normal-word polynomials. Independent explicit matrix sums test all
component assignments of word lengths 0..4 in restricted Clifford subspaces
d=1..4: **498 cases**, including four empty words.

The live Cadabra oracle independently joins gammas and contracts the metric
with an Integer range 1..d. Since Cadabra uses Clifford metric h while the
displayed contraction lowers with eta=-h, the input sandwich has an explicit
minus sign. Both residuals vanish; a wrong metric-sign control does not.

This does not define a d-dimensional gamma5, an epsilon prescription, or a
complete dimensional-regularization scheme.

## 50.4: evaluated chiral trace and epsilon sign

Define

    A=(p.q)(r.s)-(p.r)(q.s)+(p.s)(q.r).

With epsilon^{0123}=+1 and gamma5=i gamma0 gamma1 gamma2 gamma3,

    Tr[(1-gamma5) slash(p)slash(q)slash(r)slash(s)/2]
      = 2A + 2i epsilon^{mu nu rho sigma} p_mu q_nu r_rho s_sigma.

If D is the ordinary determinant of rows of **contravariant** components
(p^0,p^1,p^2,p^3), etc., then lowering the contracted vectors contributes
one minus sign:

    epsilon^{mu nu rho sigma} p_mu q_nu r_rho s_sigma = -D.

Therefore the result is **2A-2iD**, not 2A+2iD. The latter is the
opposite-chirality trace.

New lorentz::chiral_four_slash computes this formula. Independent tests:
- expand actual chiral matrices for sixteen independent vector components;
- compare both chiral traces with the evaluated invariant expressions;
- substitute massless twistor momenta and verify the bracket product
  <pq>[qr]<rs>[sp] equals the evaluated result;
- reject the reversed epsilon sign as a nonzero polynomial;
- anchor determinant orientation at the coordinate basis.

A concrete null-vector witness is

    p=(1, 0, 0, 1), q=(1, 1, 0, 0), r=(1, 0, 1, 0), s=(1, 0, 0,-1).

It has A=2, D=-2 and left trace 4+4i. A rational proper boost preserves this
result; spatial reflection exchanges the two chiral expressions.

## 50.2a: exact angular constructor and witnesses

Write s=sin(theta/2), c=cos(theta/2), x=cos(phi), y=sin(phi), with
s^2+c^2=x^2+y^2=1. Supply a=sqrt(2omega) as an exact positive rational.
Then

    phi = a*(-s*(x-i*y), c),
    p = (a^2/2)*(1, 2scx, 2scy, c^2-s^2).

Direct multiplication gives

    -phi phi^dagger
       = a^2*[[-s^2, sc(x-i*y)], [sc(x+i*y), -c^2]]
       = p_mu sigma^mu.

The last equality uses the two unit-circle constraints. This is the general
algebraic derivation; it is not inferred from sampled answers.

Twistor::from_rational_angles validates signs/ranges and the unit circles,
constructs phi and its true conjugate, and never approximates a square root.
Checked-i128 intermediate overflow remains possible for very large rational
parameters and is documented separately from domain errors.

Forty-eight rational witnesses cover poles, both hemispheres, nonzero complex
phases and three scales. They compare phi/star directly (so phase mistakes
cannot hide behind momentum invariance), independently constructed momenta,
nullness, and p_mu sigma^mu versus the outer product. Off-circle, out-of-range,
negative-scale and huge-invalid-angle controls are rejected.

50.2b's explicit boosted-spinor massless limit is not yet implemented.
