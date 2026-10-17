//! Four-dimensional scalar/pseudoscalar trace tools in Srednicki conventions.
//! Arguments contain contravariant components; epsilon^{0123}=+1,
//! epsilon_{0123}=-1, eta=(-,+,+,+), gamma5=i gamma0 gamma1 gamma2 gamma3.
use crate::{
    matrices::ETA,
    polynomial::Poly,
    scalar::C,
    twistor::{dot, FourVector},
};
#[derive(Clone, Copy, Debug)]
pub enum Chirality {
    Left,
    Right,
}
/// epsilon^{mu nu rho sigma} p_mu q_nu r_rho s_sigma.
/// This is MINUS the determinant of rows of contravariant components.
pub fn epsilon_upper_covariant(v: &[FourVector; 4]) -> Poly {
    let mut out = Poly::default();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let indices = [a, b, c, d];
                    if (0..4).any(|i| (i + 1..4).any(|j| indices[i] == indices[j])) {
                        continue;
                    }
                    let inversions = (0..4)
                        .map(|i| (i + 1..4).filter(|j| indices[i] > indices[*j]).count())
                        .sum::<usize>();
                    let epsilon = if inversions.is_multiple_of(2) { 1 } else { -1 };
                    let product = (0..4).fold(Poly::constant(epsilon), |s, i| {
                        s * Poly::constant(ETA[indices[i]]) * v[i][indices[i]].clone()
                    });
                    out = out + product;
                }
            }
        }
    }
    out
}
/// Tr[(1-/+gamma5) slash(p) slash(q) slash(r) slash(s) / 2].
/// Valid for arbitrary commuting four-vectors, not only massless momenta.
pub fn chiral_four_slash(v: &[FourVector; 4], chirality: Chirality) -> Poly {
    let invariant = dot(&v[0], &v[1]) * dot(&v[2], &v[3]) - dot(&v[0], &v[2]) * dot(&v[1], &v[3])
        + dot(&v[0], &v[3]) * dot(&v[1], &v[2]);
    let sign = match chirality {
        Chirality::Left => 2,
        Chirality::Right => -2,
    };
    Poly::constant(2) * invariant
        + Poly::constant(C::from(sign) * C::I) * epsilon_upper_covariant(v)
}
