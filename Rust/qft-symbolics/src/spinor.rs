//! Commuting external spinors; no Grassmann fields or reality condition implied.
use crate::{
    polynomial::{Poly, Ratio},
    scalar::Q,
};
pub type Spinor = [Poly; 2];
pub type RealSpinor = [Q; 2];
/// Determinant convention a0*b1-a1*b0 (NOT the SMGA frame convention).
pub fn bracket(a: &Spinor, b: &Spinor) -> Poly {
    a[0].clone() * b[1].clone() - a[1].clone() * b[0].clone()
}
pub fn named(name: &str) -> Spinor {
    [Poly::var(format!("{name}0")), Poly::var(format!("{name}1"))]
}
pub fn schouten(a: &Spinor, b: &Spinor, c: &Spinor, d: &Spinor) -> Poly {
    bracket(a, b) * bracket(c, d) + bracket(a, c) * bracket(d, b) + bracket(a, d) * bracket(b, c)
}
pub fn parke_taylor(legs: &[Spinor], negative: [usize; 2]) -> Result<Ratio, String> {
    if legs.len() < 3 || negative[0] == negative[1] || negative.iter().any(|i| *i >= legs.len()) {
        return Err("need n>=3 and two distinct in-range negative-helicity legs".into());
    }
    let mut denominator = Poly::constant(1);
    for i in 0..legs.len() {
        denominator = denominator * bracket(&legs[i], &legs[(i + 1) % legs.len()]);
    }
    Ratio::new(
        bracket(&legs[negative[0]], &legs[negative[1]]).pow(4),
        denominator,
    )
}
/// [i,j> shift: lambda_i += z lambda_j, tilde_j -= z tilde_i.
/// Acts on independent rational spinors; does not impose Lorentzian conjugacy.
pub fn bcfw_shift(
    lambda: &mut [RealSpinor],
    tilde: &mut [RealSpinor],
    i: usize,
    j: usize,
    z: Q,
) -> Result<(), String> {
    if lambda.len() != tilde.len() || i == j || i >= lambda.len() || j >= lambda.len() {
        return Err("invalid BCFW legs".into());
    }
    for a in 0..2 {
        lambda[i][a] = lambda[i][a] + z * lambda[j][a];
        tilde[j][a] = tilde[j][a] - z * tilde[i][a];
    }
    Ok(())
}
pub fn momentum(lambda: RealSpinor, tilde: RealSpinor) -> [[Q; 2]; 2] {
    std::array::from_fn(|a| std::array::from_fn(|b| lambda[a] * tilde[b]))
}
