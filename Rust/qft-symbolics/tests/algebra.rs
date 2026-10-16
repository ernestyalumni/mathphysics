use qft_symbolics::{
    clifford, matrices as m,
    polynomial::{Poly, Ratio},
    scalar::{C, Q},
    smga,
    spinor::*,
};
use std::collections::BTreeMap;
#[test]
fn rationals_and_gaussian() {
    assert_eq!(Q::new(-6, -8), Q::new(3, 4));
    assert_eq!(Q::new(1, 6) + Q::new(1, 3), Q::new(1, 2));
    assert_eq!(Q::new(i128::MAX, 2) * Q::new(2, i128::MAX), Q::ONE);
    assert_eq!(C::I * C::I, C::from(-1));
    let z = C::new(Q::new(2, 3), Q::new(-4, 5));
    assert_eq!(z / z, C::ONE);
    assert_eq!(z.conj().conj(), z);
}
#[test]
#[should_panic(expected = "overflow")]
fn overflow_fails_closed() {
    let _ = Q::new(i128::MAX, 1) + Q::ONE;
}
#[test]
#[should_panic(expected = "zero")]
fn zero_denominator_fails_closed() {
    let _ = Q::new(1, 0);
}
#[test]
fn polynomial_ring_and_poles() {
    let x = Poly::var("x");
    let y = Poly::var("y");
    assert_eq!(
        (x.clone() + y.clone()) * (x.clone() - y.clone()),
        x.pow(2) - y.pow(2)
    );
    assert!(Ratio::new(x.clone(), Poly::default()).is_err());
    let r = Ratio::new(x.clone(), x).unwrap();
    assert!(r
        .residual(&Ratio::new(Poly::constant(1), Poly::constant(1)).unwrap())
        .is_zero());
    assert!(r.eval(&BTreeMap::from([("x".into(), C::ZERO)])).is_err());
    assert!(r.eval(&BTreeMap::new()).is_err());
}
#[test]
fn schouten_is_polynomial_identity_and_wrong_sign_is_not() {
    let (a, b, c, d) = (named("a"), named("b"), named("c"), named("d"));
    assert!(schouten(&a, &b, &c, &d).is_zero());
    let wrong = bracket(&a, &b) * bracket(&c, &d) - bracket(&a, &c) * bracket(&d, &b)
        + bracket(&a, &d) * bracket(&b, &c);
    assert!(!wrong.is_zero());
    assert!((bracket(&a, &b) + bracket(&b, &a)).is_zero());
    assert!(bracket(&a, &a).is_zero());
}
#[test]
fn parke_taylor_cyclicity_reflection_and_weights() {
    let legs: Vec<_> = (0..5).map(|i| named(&format!("p{i}"))).collect();
    let a = parke_taylor(&legs, [0, 2]).unwrap();
    let mut cyclic = legs[1..].to_vec();
    cyclic.push(legs[0].clone());
    assert!(a
        .residual(&parke_taylor(&cyclic, [4, 1]).unwrap())
        .is_zero());
    let reverse: Vec<_> = legs.iter().rev().cloned().collect();
    let b = parke_taylor(&reverse, [4, 2]).unwrap();
    assert!((a.numerator().clone() * b.denominator().clone()
        + b.numerator().clone() * a.denominator().clone())
    .is_zero());
    // Exact evaluation: scaling one lambda by t gives t^(4-2) or t^(-2).
    let vals: BTreeMap<_, _> = (0..5)
        .flat_map(|i| {
            [
                (format!("p{i}0"), C::ONE),
                (format!("p{i}1"), C::from((i * i + 1) as i64)),
            ]
        })
        .collect();
    let base = a.eval(&vals).unwrap();
    for i in 0..5 {
        let mut scaled = vals.clone();
        for k in 0..2 {
            let key = format!("p{i}{k}");
            scaled.insert(key.clone(), vals[&key] * C::from(3));
        }
        let factor = if i == 0 || i == 2 {
            Q::from(9)
        } else {
            Q::new(1, 9)
        };
        assert_eq!(a.eval(&scaled).unwrap(), base * C::from(factor));
    }
    assert!(parke_taylor(&legs, [1, 1]).is_err());
}
#[test]
fn bcfw_shift_preserves_rank_one_and_total_momentum() {
    let mut l = vec![[Q::from(1), Q::from(2)], [Q::from(3), Q::from(5)]];
    let mut t = vec![[Q::from(7), Q::from(11)], [Q::from(13), Q::from(17)]];
    let total = |l: &[RealSpinor], t: &[RealSpinor]| {
        let a = momentum(l[0], t[0]);
        let b = momentum(l[1], t[1]);
        std::array::from_fn::<_, 2, _>(|i| std::array::from_fn::<_, 2, _>(|j| a[i][j] + b[i][j]))
    };
    let before = total(&l, &t);
    bcfw_shift(&mut l, &mut t, 0, 1, Q::new(3, 7)).unwrap();
    assert_eq!(before, total(&l, &t));
    for i in 0..2 {
        let p = momentum(l[i], t[i]);
        assert_eq!(p[0][0] * p[1][1] - p[0][1] * p[1][0], Q::ZERO);
    }
    assert!(bcfw_shift(&mut l, &mut t, 0, 0, Q::ONE).is_err());
}
#[test]
fn ch35_1_raised_sigma_and_completeness() {
    let eu = m::epsilon_upper();
    let el = m::epsilon_lower();
    assert_eq!(m::mul(el, eu), m::identity());
    for mu in 0..4 {
        let s = m::sigma(mu);
        let raised = std::array::from_fn(|da| {
            std::array::from_fn(|a| {
                (0..2)
                    .flat_map(|b| (0..2).map(move |db| eu[a][b] * eu[da][db] * s[b][db]))
                    .fold(C::ZERO, |a, b| a + b)
            })
        });
        assert_eq!(raised, m::barsigma(mu));
    }
    for a in 0..2 {
        for da in 0..2 {
            for b in 0..2 {
                for db in 0..2 {
                    let lhs = (0..4).fold(C::ZERO, |s, mu| {
                        s + C::from(m::ETA[mu]) * m::sigma(mu)[a][da] * m::sigma(mu)[b][db]
                    });
                    assert_eq!(lhs, C::from(-2) * el[a][b] * el[da][db]);
                }
            }
        }
    }
}
#[test]
fn ch35_2_lorentz_generators() {
    let s = m::left_generator;
    for i in 1..4 {
        assert_eq!(
            s(i, 0),
            m::scale(m::sigma(i), C::new(Q::ZERO, Q::new(1, 2)))
        );
        assert_eq!(s(i, i), m::zero());
    }
    for (i, j, k) in [(1, 2, 3), (2, 3, 1), (3, 1, 2)] {
        assert_eq!(s(i, j), m::scale(m::sigma(k), C::from(Q::new(1, 2))));
    }
    // All 4^4 Lorentz algebra components, including boosts.
    for mu in 0..4 {
        for nu in 0..4 {
            for rho in 0..4 {
                for sig in 0..4 {
                    let comm = m::add(
                        m::mul(s(mu, nu), s(rho, sig)),
                        m::scale(m::mul(s(rho, sig), s(mu, nu)), C::from(-1)),
                    );
                    let mut rhs = m::zero();
                    for (a, b, c, d, sign) in [
                        (mu, rho, nu, sig, 1),
                        (nu, rho, mu, sig, -1),
                        (mu, sig, nu, rho, -1),
                        (nu, sig, mu, rho, 1),
                    ] {
                        if a == b {
                            rhs = m::add(rhs, m::scale(s(c, d), C::I * C::from(sign * m::ETA[a])));
                        }
                    }
                    assert_eq!(comm, rhs);
                }
            }
        }
    }
}
#[test]
fn ch47_symbolic_trace_and_clifford_anticommutator() {
    let tr4 = clifford::trace(&["a", "b", "c", "d"]);
    let expected = Poly::constant(4)
        * (clifford::metric("a", "b") * clifford::metric("c", "d")
            - clifford::metric("a", "c") * clifford::metric("b", "d")
            + clifford::metric("a", "d") * clifford::metric("b", "c"));
    assert_eq!(tr4, expected);
    assert!(clifford::trace(&["a", "b", "c"]).is_zero());
    let ab = clifford::normal_order(&["a", "b"]);
    let ba = clifford::normal_order(&["b", "a"]);
    assert!((ab[&vec!["a".into(), "b".into()]].clone()
        + ba[&vec!["a".into(), "b".into()]].clone())
    .is_zero());
    assert_eq!(
        ba[&Vec::<String>::new()],
        Poly::constant(2) * clifford::metric("a", "b")
    );
    assert_eq!(
        clifford::normal_order(&["a", "a"])[&Vec::<String>::new()],
        clifford::metric("a", "a")
    );
}
#[test]
fn ch47_symbolic_words_match_exact_gamma_matrices() {
    // Every 4-gamma component word, independently checking trace recurrence
    // AND noncommutative normal-ordering against explicit chiral matrices.
    let names = ["a", "b", "c", "d"];
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                for d in 0..4 {
                    let components = [a, b, c, d];
                    let mut values = BTreeMap::new();
                    for i in 0..4 {
                        for j in 0..4 {
                            // metric's unique variable name, evaluated using polynomial display-free API
                            let h = if components[i] == components[j] {
                                -m::ETA[components[i]]
                            } else {
                                0
                            };
                            let key = format!("h(1:{};1:{})", names[i.min(j)], names[i.max(j)]);
                            values.insert(key, C::from(h));
                        }
                    }
                    let direct = components
                        .iter()
                        .fold(m::identity(), |p, mu| m::mul(p, m::gamma(*mu)));
                    assert_eq!(
                        clifford::trace(&names).eval(&values).unwrap(),
                        m::trace(direct)
                    );
                    let reverse = ["d", "c", "b", "a"];
                    let direct_reverse = components
                        .iter()
                        .rev()
                        .fold(m::identity(), |p, mu| m::mul(p, m::gamma(*mu)));
                    let mut reduced = m::zero();
                    for (word, p) in clifford::normal_order(&reverse) {
                        let mat = word.iter().fold(m::identity(), |p, label| {
                            let i = names.iter().position(|s| s == label).unwrap();
                            m::mul(p, m::gamma(components[i]))
                        });
                        reduced = m::add(reduced, m::scale(mat, p.eval(&values).unwrap()));
                    }
                    assert_eq!(reduced, direct_reverse);
                }
            }
        }
    }
}
#[test]
fn ch50_1a_completeness_symbolic_components() {
    // phi and its conjugate are treated as independent variables: stronger
    // polynomial identity, restricting to conjugate values is then allowed.
    let phi = named("phi");
    let star = named("star");
    let raise = |v: &Spinor| [v[1].clone(), -v[0].clone()];
    let up = raise(&phi);
    let star_up = raise(&star);
    let mut slash: Vec<Vec<Poly>> = vec![vec![Poly::default(); 4]; 4];
    for a in 0..2 {
        for da in 0..2 {
            slash[a][da + 2] = -phi[a].clone() * star[da].clone();
            slash[da + 2][a] = -up[a].clone() * star_up[da].clone();
        }
    }
    let uminus = [
        phi[0].clone(),
        phi[1].clone(),
        Poly::default(),
        Poly::default(),
    ];
    let uplus = [
        Poly::default(),
        Poly::default(),
        star_up[0].clone(),
        star_up[1].clone(),
    ];
    let barminus = [
        Poly::default(),
        Poly::default(),
        star[0].clone(),
        star[1].clone(),
    ];
    let barplus = [
        up[0].clone(),
        up[1].clone(),
        Poly::default(),
        Poly::default(),
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert!((slash[i][j].clone()
                + uminus[i].clone() * barminus[j].clone()
                + uplus[i].clone() * barplus[j].clone())
            .is_zero());
        }
    }
}
#[test]
fn smga_partitions_walls_and_domain() {
    assert_eq!(smga::partitions(0, 4, 1).len(), 8);
    assert_eq!(smga::partitions(0, 4, 3).len(), 4);
    let w = vec![[Q::ONE, Q::ONE]; 3];
    assert!(smga::formula_r1(&w).is_err());
    let walls = vec![[Q::from(-2), Q::ZERO], [Q::ONE, Q::ZERO], [Q::ONE, Q::ZERO]];
    assert!(smga::formula_r1(&walls).is_err());
    assert!(smga::minus_first_recursion(&walls).is_err());
}
#[test]
fn smga_exact_bg_vs_formula_n3_to_n10() {
    let (mut zero, mut nonzero) = (0, 0);
    for n in 3..=10 {
        for trial in 0..25 {
            let v = smga::fixture(n, trial);
            let a = smga::formula_r1(&v).unwrap_or_else(|e| {
                panic!("unexpected formula rejection n={n}, trial={trial}: {e}")
            });
            let b = smga::minus_first_recursion(&v).unwrap_or_else(|e| {
                panic!("unexpected recurrence rejection n={n}, trial={trial}: {e}")
            });
            assert_eq!(a, b, "n={n}, trial={trial}");
            if a.is_zero() {
                zero += 1;
            } else {
                nonzero += 1;
            }
        }
    }
    assert_eq!(zero + nonzero, 200);
    assert!(zero > 0 && nonzero > 0);
    eprintln!("SMGA exact fixtures: zero={zero} nonzero={nonzero} excluded_walls=0");
}
