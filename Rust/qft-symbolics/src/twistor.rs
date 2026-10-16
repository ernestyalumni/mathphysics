//! Srednicki's commuting Weyl twistors, with eta=(-,+,+,+).
//! The conjugate components are algebraically independent variables, so
//! polynomial identities remain valid after imposing Lorentzian conjugacy.
use crate::{
    matrices,
    polynomial::{Poly, Ratio},
    scalar::{C, Q},
    spinor::{bracket, named, Spinor},
};
#[derive(Clone)]
pub struct Twistor {
    pub phi: Spinor,
    pub star: Spinor,
}
pub type FourVector = [Poly; 4];
pub type DiracVector = [Poly; 4];
pub fn raise(v: &Spinor) -> Spinor {
    [v[1].clone(), -v[0].clone()]
}
pub fn square(p: &Twistor, k: &Twistor) -> Poly {
    -bracket(&p.phi, &k.phi)
}
pub fn angle(p: &Twistor, k: &Twistor) -> Poly {
    bracket(&p.star, &k.star)
}
impl Twistor {
    pub fn named(label: &str) -> Self {
        Self {
            phi: named(&format!("{label}p")),
            star: named(&format!("{label}s")),
        }
    }
    pub fn ket_minus(&self) -> DiracVector {
        [
            self.phi[0].clone(),
            self.phi[1].clone(),
            Poly::default(),
            Poly::default(),
        ]
    }
    pub fn ket_plus(&self) -> DiracVector {
        let s = raise(&self.star);
        [Poly::default(), Poly::default(), s[0].clone(), s[1].clone()]
    }
    pub fn bar_minus(&self) -> DiracVector {
        [
            Poly::default(),
            Poly::default(),
            self.star[0].clone(),
            self.star[1].clone(),
        ]
    }
    pub fn bar_plus(&self) -> DiracVector {
        let p = raise(&self.phi);
        [p[0].clone(), p[1].clone(), Poly::default(), Poly::default()]
    }
    /// p_(a,dota)=-phi_a star_dota; invert p_mu sigma^mu.
    pub fn momentum(&self) -> FourVector {
        let p: Vec<Vec<_>> = (0..2)
            .map(|a| {
                (0..2)
                    .map(|b| -self.phi[a].clone() * self.star[b].clone())
                    .collect()
            })
            .collect();
        let half = Poly::constant(Q::new(1, 2));
        [
            -half.clone() * (p[0][0].clone() + p[1][1].clone()),
            half.clone() * (p[0][1].clone() + p[1][0].clone()),
            Poly::constant(C::new(Q::ZERO, Q::new(-1, 2))) * (p[1][0].clone() - p[0][1].clone()),
            half * (p[0][0].clone() - p[1][1].clone()),
        ]
    }
}
/// Slash of a contravariant momentum, p_mu gamma^mu.
pub fn slash(p: &FourVector) -> [[Poly; 4]; 4] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..4).fold(Poly::default(), |s, mu| {
                s + Poly::constant(matrices::ETA[mu])
                    * p[mu].clone()
                    * Poly::constant(matrices::gamma(mu)[i][j])
            })
        })
    })
}
pub fn dot(a: &FourVector, b: &FourVector) -> Poly {
    (0..4).fold(Poly::default(), |s, mu| {
        s + Poly::constant(matrices::ETA[mu]) * a[mu].clone() * b[mu].clone()
    })
}
pub fn sandwich_gamma(bar: &DiracVector, ket: &DiracVector) -> FourVector {
    std::array::from_fn(|mu| {
        let g = matrices::gamma(mu);
        (0..4).fold(Poly::default(), |s, i| {
            s + (0..4).fold(Poly::default(), |t, j| {
                t + bar[i].clone() * Poly::constant(g[i][j]) * ket[j].clone()
            })
        })
    })
}
/// Numerators of sqrt(2)*epsilon_+ and sqrt(2)*epsilon_-.
/// Denominators are respectively <q k> and [q k], kept explicitly.
pub fn polarization_numerator(k: &Twistor, q: &Twistor, plus: bool) -> FourVector {
    if plus {
        sandwich_gamma(&q.bar_minus(), &k.ket_minus()).map(|p| -p)
    } else {
        sandwich_gamma(&q.bar_plus(), &k.ket_plus()).map(|p| -p)
    }
}
pub fn polarization(k: &Twistor, q: &Twistor, plus: bool) -> Result<[Ratio; 4], String> {
    let d = if plus { angle(q, k) } else { square(q, k) };
    if d.is_zero() {
        return Err("singular polarization reference; q=k is not allowed".into());
    }
    let n = polarization_numerator(k, q, plus);
    Ok(n.map(|p| Ratio::new(p, d.clone()).expect("validated denominator")))
}
