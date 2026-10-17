//! Evaluated 50.4, with independent epsilon-sign controls, and 50.2a witnesses.
use qft_symbolics::{
    lorentz::{self as l, Chirality},
    matrices as m,
    polynomial::Poly,
    scalar::{C, Q},
    twistor::{self as t, FourVector, Twistor},
};
use std::collections::BTreeMap;
fn det4(v: &[FourVector; 4]) -> Poly {
    let mut out = Poly::default();
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let p = [a, b, c, d];
                    if (0..4).any(|i| (i + 1..4).any(|j| p[i] == p[j])) {
                        continue;
                    }
                    let inversions = (0..4)
                        .map(|i| (i + 1..4).filter(|j| p[i] > p[*j]).count())
                        .sum::<usize>();
                    let sign = if inversions.is_multiple_of(2) { 1 } else { -1 };
                    let product = (0..4).fold(Poly::constant(sign), |s, i| s * v[i][p[i]].clone());
                    out = out + product;
                }
            }
        }
    }
    out
}
fn pairings(v: &[FourVector; 4]) -> Poly {
    t::dot(&v[0], &v[1]) * t::dot(&v[2], &v[3]) - t::dot(&v[0], &v[2]) * t::dot(&v[1], &v[3])
        + t::dot(&v[0], &v[3]) * t::dot(&v[1], &v[2])
}
fn evaluated(v: &[FourVector; 4], left: bool) -> Poly {
    Poly::constant(2) * pairings(v)
        + Poly::constant(C::from(if left { -2 } else { 2 }) * C::I) * det4(v)
}
fn traces(v: &[FourVector; 4]) -> (Poly, Poly) {
    let mut product: [[Poly; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| Poly::constant(if i == j { 1 } else { 0 }))
    });
    for vector in v {
        let slash = t::slash(vector);
        product = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..4).fold(Poly::default(), |s, k| {
                    s + product[i][k].clone() * slash[k][j].clone()
                })
            })
        });
    }
    (
        product[0][0].clone() + product[1][1].clone(),
        product[2][2].clone() + product[3][3].clone(),
    )
}
#[test]
fn ch50_4_evaluated_trace_for_arbitrary_symbolic_vectors() {
    let g5 = m::scale(
        (0..4).fold(m::identity(), |p, mu| m::mul(p, m::gamma(mu))),
        C::I,
    );
    assert_eq!(
        g5,
        std::array::from_fn(|i| std::array::from_fn(|j| C::from(if i != j {
            0
        } else if i < 2 {
            -1
        } else {
            1
        })))
    );
    let v = std::array::from_fn(|i| std::array::from_fn(|mu| Poly::var(format!("v{i}_{mu}"))));
    let (left, right) = traces(&v);
    assert!((left.clone() - evaluated(&v, true)).is_zero());
    assert_eq!(
        l::chiral_four_slash(&v, Chirality::Left),
        evaluated(&v, true)
    );
    assert_eq!(
        l::chiral_four_slash(&v, Chirality::Right),
        evaluated(&v, false)
    );
    assert_eq!(l::epsilon_upper_covariant(&v), -det4(&v));
    assert!((right - evaluated(&v, false)).is_zero());
    assert!(
        !(left - evaluated(&v, false)).is_zero(),
        "wrong epsilon sign must not pass"
    );
}
#[test]
fn ch50_4_bracket_product_equals_evaluated_invariants() {
    let legs: [Twistor; 4] = std::array::from_fn(|i| Twistor::named(&format!("t{i}")));
    let vectors = legs.clone().map(|t| t.momentum());
    let brackets = t::angle(&legs[0], &legs[1])
        * t::square(&legs[1], &legs[2])
        * t::angle(&legs[2], &legs[3])
        * t::square(&legs[3], &legs[0]);
    assert!((brackets - evaluated(&vectors, true)).is_zero());
}
#[test]
fn epsilon_orientation_covariance_and_nonzero_sign_witness() {
    let basis = std::array::from_fn(|i| {
        std::array::from_fn(|j| Poly::constant(if i == j { 1 } else { 0 }))
    });
    assert_eq!(det4(&basis), Poly::constant(1)); // component determinant, NOT epsilon_lower0123
    let v =
        [[1, 0, 0, 1], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, -1]].map(|p| p.map(Poly::constant));
    assert_eq!(det4(&v), Poly::constant(-2));
    assert_eq!(pairings(&v), Poly::constant(2));
    assert_eq!(traces(&v).0, Poly::constant(C::new(Q::from(4), Q::from(4))));
    // Proper rational boost in x with cosh=5/4 and sinh=3/4.
    let boosted = v.clone().map(|p| {
        [
            Poly::constant(Q::new(5, 4)) * p[0].clone()
                + Poly::constant(Q::new(3, 4)) * p[1].clone(),
            Poly::constant(Q::new(3, 4)) * p[0].clone()
                + Poly::constant(Q::new(5, 4)) * p[1].clone(),
            p[2].clone(),
            p[3].clone(),
        ]
    });
    assert_eq!(evaluated(&v, true), evaluated(&boosted, true));
    let reflected = v.clone().map(|mut p| {
        p[1] = -p[1].clone();
        p
    });
    assert_eq!(evaluated(&v, true), evaluated(&reflected, false));
}
#[test]
fn rational_angular_constructor_rejects_off_circle_and_nonphysical_parameters() {
    let q = Q::from;
    for args in [
        [q(0), q(0), q(1), q(1), q(0)],
        [q(-1), q(0), q(1), q(1), q(0)],
        [q(1), q(1), q(1), q(1), q(0)],
        [q(1), q(0), q(1), q(1), q(1)],
        [q(1), q(-1), q(0), q(1), q(0)],
        [q(1), Q::new(1_i128 << 64, 1), q(0), q(1), q(0)],
    ] {
        assert!(
            Twistor::from_rational_angles(args[0], args[1], args[2], args[3], args[4]).is_err()
        );
    }
}
#[test]
fn ch50_2a_exact_rational_angular_factorization_witnesses() {
    let circles = [
        (Q::ZERO, Q::ONE),
        (Q::ONE, Q::ZERO),
        (Q::new(3, 5), Q::new(4, 5)),
        (Q::new(4, 5), Q::new(3, 5)),
    ];
    let phases = [
        (Q::ONE, Q::ZERO),
        (Q::ZERO, Q::ONE),
        (Q::new(5, 13), Q::new(12, 13)),
        (Q::new(5, 13), Q::new(-12, 13)),
    ];
    let mut witnesses = 0;
    for amplitude in [1, 2, 5] {
        for (sin_half, cos_half) in circles {
            for (cos_phi, sin_phi) in phases {
                assert_eq!(sin_half * sin_half + cos_half * cos_half, Q::ONE);
                assert_eq!(cos_phi * cos_phi + sin_phi * sin_phi, Q::ONE);
                let a = Q::from(amplitude);
                let omega = a * a / Q::from(2);
                let phi = [
                    C::from(-a * sin_half) * C::new(cos_phi, -sin_phi),
                    C::from(a * cos_half),
                ];
                let twistor =
                    Twistor::from_rational_angles(a, sin_half, cos_half, cos_phi, sin_phi).unwrap();
                assert_eq!(twistor.phi, phi.map(Poly::constant));
                assert_eq!(twistor.star, phi.map(|z| Poly::constant(z.conj())));
                let expected = [
                    omega,
                    Q::from(2) * omega * sin_half * cos_half * cos_phi,
                    Q::from(2) * omega * sin_half * cos_half * sin_phi,
                    omega * (cos_half * cos_half - sin_half * sin_half),
                ];
                let actual = twistor.momentum();
                for mu in 0..4 {
                    assert_eq!(
                        actual[mu].eval(&BTreeMap::new()).unwrap(),
                        C::from(expected[mu])
                    );
                }
                assert!(t::dot(&actual, &actual).is_zero());
                let mut sigma_momentum = m::zero();
                for (mu, component) in expected.iter().enumerate() {
                    sigma_momentum = m::add(
                        sigma_momentum,
                        m::scale(m::sigma(mu), C::from(Q::from(m::ETA[mu]) * *component)),
                    );
                }
                let rank_one =
                    std::array::from_fn(|i| std::array::from_fn(|j| -phi[i] * phi[j].conj()));
                assert_eq!(sigma_momentum, rank_one);
                witnesses += 1;
            }
        }
    }
    assert_eq!(witnesses, 48);
}
