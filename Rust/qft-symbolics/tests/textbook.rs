//! Numbered textbook subproblems and equation checks. Not whole chapters.
use qft_symbolics::{
    matrices as m,
    polynomial::Poly,
    scalar::{C, Q},
    twistor::{self as t, Twistor},
};
fn gamma5() -> m::Matrix<4> {
    m::scale(
        (0..4).fold(m::identity(), |p, mu| m::mul(p, m::gamma(mu))),
        C::I,
    )
}
#[test]
fn ch47_1_gamma5_two_gamma_trace() {
    let g5 = gamma5();
    assert_eq!(m::mul(g5, g5), m::identity());
    assert_eq!(m::trace(g5), C::ZERO);
    for mu in 0..4 {
        for nu in 0..4 {
            assert_eq!(
                m::trace(m::mul(g5, m::mul(m::gamma(mu), m::gamma(nu)))),
                C::ZERO
            );
        }
    }
}
#[test]
fn ch47_2_contracted_gamma_equations_four_dimensions_only() {
    // Eq47.20: gamma^mu gamma^a gamma^b gamma_mu =4 eta^ab (d=4).
    // Eq47.21: gamma^mu gamma^a gamma^b gamma^c gamma_mu=2 gamma^c gamma^b gamma^a (d=4).
    for a in 0..4 {
        for b in 0..4 {
            let mut lhs = m::zero();
            for mu in 0..4 {
                lhs = m::add(
                    lhs,
                    m::scale(
                        m::mul(
                            m::gamma(mu),
                            m::mul(m::gamma(a), m::mul(m::gamma(b), m::gamma(mu))),
                        ),
                        C::from(m::ETA[mu]),
                    ),
                );
            }
            let expected = if a == b { 4 * m::ETA[a] } else { 0 };
            assert_eq!(lhs, m::scale(m::identity(), C::from(expected)));
            for c in 0..4 {
                let mut lhs = m::zero();
                for mu in 0..4 {
                    let product = [mu, a, b, c, mu]
                        .iter()
                        .fold(m::identity(), |p, k| m::mul(p, m::gamma(*k)));
                    lhs = m::add(lhs, m::scale(product, C::from(m::ETA[mu])));
                }
                let rhs = m::scale(
                    [c, b, a]
                        .iter()
                        .fold(m::identity(), |p, k| m::mul(p, m::gamma(*k))),
                    C::from(2),
                );
                assert_eq!(lhs, rhs);
            }
        }
    }
}
#[test]
fn ch47_3_sixteen_dirac_basis_matrices_are_orthogonal() {
    let g5 = gamma5();
    let mut basis = vec![m::identity()];
    for mu in 0..4 {
        basis.push(m::gamma(mu));
    }
    for mu in 0..4 {
        for nu in mu + 1..4 {
            basis.push(m::scale(
                m::add(
                    m::mul(m::gamma(mu), m::gamma(nu)),
                    m::scale(m::mul(m::gamma(nu), m::gamma(mu)), C::from(-1)),
                ),
                C::new(Q::ZERO, Q::new(1, 4)),
            ));
        }
    }
    for mu in 0..4 {
        basis.push(m::mul(m::gamma(mu), g5));
    }
    basis.push(g5);
    assert_eq!(basis.len(), 16);
    for i in 0..16 {
        for j in 0..16 {
            let adj = std::array::from_fn(|a| std::array::from_fn(|b| basis[i][b][a].conj()));
            let inner = m::trace(m::mul(adj, basis[j]));
            if i == j {
                assert_eq!(inner.im, Q::ZERO);
                assert_eq!(inner.re.sign(), 1);
            } else {
                assert_eq!(inner, C::ZERO);
            }
        }
    }
}
#[test]
fn ch50_1b_slash_sandwiches_are_polynomial_identities() {
    let (p, k, r) = (
        Twistor::named("p"),
        Twistor::named("k"),
        Twistor::named("r"),
    );
    let mom = k.momentum();
    for plus in [false, true] {
        let bar = if plus { r.bar_plus() } else { r.bar_minus() };
        let ket = if plus { p.ket_plus() } else { p.ket_minus() };
        let lhs = -t::dot(&mom, &t::sandwich_gamma(&bar, &ket));
        let rhs = if plus {
            t::square(&r, &k) * t::angle(&k, &p)
        } else {
            t::angle(&r, &k) * t::square(&k, &p)
        };
        assert!((lhs - rhs).is_zero());
    }
    assert!(t::dot(&mom, &t::sandwich_gamma(&r.bar_minus(), &p.ket_plus())).is_zero());
    assert!(t::dot(&mom, &t::sandwich_gamma(&r.bar_plus(), &p.ket_minus())).is_zero());
    assert!(t::dot(&mom, &mom).is_zero());
    assert!((t::angle(&p, &k) * t::square(&k, &p)
        + Poly::constant(2) * t::dot(&p.momentum(), &mom))
    .is_zero());
}
#[test]
fn ch50_4_chiral_four_slash_trace() {
    let (p, q, r, s) = (
        Twistor::named("p"),
        Twistor::named("q"),
        Twistor::named("r"),
        Twistor::named("s"),
    );
    let mut product: [[Poly; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| Poly::constant(if i == j { 1 } else { 0 }))
    });
    for t in [&p, &q, &r, &s] {
        let slash = t::slash(&t.momentum());
        product = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..4).fold(Poly::default(), |sum, k| {
                    sum + product[i][k].clone() * slash[k][j].clone()
                })
            })
        });
    }
    let trace = product[0][0].clone() + product[1][1].clone();
    let brackets = t::angle(&p, &q) * t::square(&q, &r) * t::angle(&r, &s) * t::square(&s, &p);
    assert!((trace - brackets).is_zero());
}
#[test]
fn ch60_2a_momentum_conservation_factorization() {
    let (i, k) = (Twistor::named("i"), Twistor::named("k"));
    let legs: Vec<_> = (0..5).map(|j| Twistor::named(&format!("j{j}"))).collect();
    let total = std::array::from_fn(|mu| {
        legs.iter()
            .fold(Poly::default(), |s, j| s + j.momentum()[mu].clone())
    });
    let left = legs.iter().fold(Poly::default(), |s, j| {
        s + t::angle(&i, j) * t::square(j, &k)
    });
    let right = legs.iter().fold(Poly::default(), |s, j| {
        s + t::square(&i, j) * t::angle(j, &k)
    });
    // Each sum is minus a fixed bilinear contracted with sum p_j.
    // Therefore both vanish under momentum conservation; no finite-point inference.
    assert!((left + t::dot(&total, &t::sandwich_gamma(&i.bar_minus(), &k.ket_minus()))).is_zero());
    assert!((right + t::dot(&total, &t::sandwich_gamma(&i.bar_plus(), &k.ket_plus()))).is_zero());
}
#[test]
fn ch60_1a_polarization_dot_products_and_transversality() {
    let (p, k, q) = (
        Twistor::named("p"),
        Twistor::named("k"),
        Twistor::named("q"),
    );
    for plus in [false, true] {
        let num = t::polarization_numerator(&k, &q, plus);
        let lhs = t::dot(&p.momentum(), &num);
        let rhs = if plus {
            t::angle(&q, &p) * t::square(&p, &k)
        } else {
            t::square(&q, &p) * t::angle(&p, &k)
        };
        assert!((lhs - rhs).is_zero());
        assert!(t::dot(&k.momentum(), &num).is_zero());
        assert!(t::dot(&q.momentum(), &num).is_zero());
        assert!(t::polarization(&k, &q, plus).is_ok());
        assert!(t::polarization(&k, &k, plus).is_err());
    }
}
#[test]
fn ch60_1b_polarization_pair_products() {
    let (k, q, l, r) = (
        Twistor::named("k"),
        Twistor::named("q"),
        Twistor::named("l"),
        Twistor::named("r"),
    );
    for (plus, other) in [(true, true), (false, false), (true, false)] {
        let lhs = t::dot(
            &t::polarization_numerator(&k, &q, plus),
            &t::polarization_numerator(&l, &r, other),
        );
        let rhs = Poly::constant(2)
            * match (plus, other) {
                (true, true) => t::angle(&q, &r) * t::square(&k, &l),
                (false, false) => t::square(&q, &r) * t::angle(&k, &l),
                _ => t::angle(&q, &l) * t::square(&k, &r),
            };
        assert!((lhs - rhs).is_zero());
    }
}
