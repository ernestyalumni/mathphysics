//! Independent acceptance checks, each bounded; see docs/SREDNICKI-CHECKS.md.
use qft_symbolics::{clifford as cl, conventions::EpsilonConvention as E};
use qft_symbolics::{
    computation::{Budget, Limits},
    exact::{BigC, BigQ},
    polynomial::{BigPoly as P, Poly, Ratio},
    scalar::{C, Q},
};
use std::collections::BTreeMap;
fn k(n: i64) -> P {
    P::constant(n)
}
fn bracket(a: &[P; 2], c: &[P; 2], b: &mut Budget) -> P {
    E::Library.bracket(a, c, b).unwrap()
}

#[test]
fn ch35_profiles_preserve_components_and_expose_signs() {
    let mut b = Budget::default();
    let v = [P::var("v1"), P::var("v2")];
    for profile in [E::Library, E::MonoNotes] {
        let raised = profile.raise(&v, &mut b).unwrap();
        assert_eq!(profile.lower(&raised, &mut b).unwrap(), v);
        assert_eq!(
            profile.raise(&raised, &mut b).unwrap(),
            [-v[0].clone(), -v[1].clone()]
        );
        for target in [E::Library, E::MonoNotes] {
            let upper = profile.convert_upper(target, &raised, &mut b).unwrap();
            assert_eq!(target.lower(&upper, &mut b).unwrap(), v);
            assert_eq!(profile.convert_lower(target, &v, &mut b).unwrap(), v);
            assert_eq!(
                target.convert_upper(profile, &upper, &mut b).unwrap(),
                raised
            );
        }
    }
    let a = [k(2), k(3)];
    let c = [k(5), k(7)];
    assert_eq!(E::Library.raise(&a, &mut b).unwrap(), [k(3), k(-2)]);
    assert_eq!(E::MonoNotes.raise(&a, &mut b).unwrap(), [k(-3), k(2)]);
    assert_eq!(E::Library.bracket(&a, &c, &mut b).unwrap(), k(-1));
    assert_eq!(E::MonoNotes.bracket(&a, &c, &mut b).unwrap(), k(1));
    assert_ne!(
        E::MonoNotes
            .lower(&E::Library.raise(&a, &mut b).unwrap(), &mut b)
            .unwrap(),
        a
    );
    assert!(E::Library
        .convert_lower(
            E::Library,
            &v,
            &mut Budget::new(Limits {
                max_operations: 0,
                ..Limits::default()
            })
        )
        .is_err());
}
// Oracle: explicit Pauli entries and ordinary loops, never production matrices,
// spinor/twistor helpers, metric construction, or Clifford trace recurrence.
type M = [[C; 4]; 4];
fn sigma(mu: usize) -> [[C; 2]; 2] {
    let z = C::ZERO;
    let o = C::ONE;
    match mu {
        0 => [[o, z], [z, o]],
        1 => [[z, o], [o, z]],
        2 => [[z, -C::I], [C::I, z]],
        3 => [[o, z], [z, -o]],
        _ => unreachable!(),
    }
}
fn identity() -> M {
    std::array::from_fn(|i| std::array::from_fn(|j| C::from(i64::from(i == j))))
}
fn gamma(mu: usize) -> M {
    let s = sigma(mu);
    let mut g = [[C::ZERO; 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            g[i][j + 2] = s[i][j];
            g[i + 2][j] = if mu == 0 { s[i][j] } else { -s[i][j] };
        }
    }
    g
}
fn mul(a: M, c: M) -> M {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..4).fold(C::ZERO, |s, k| s + a[i][k] * c[k][j]))
    })
}
fn word(w: &[usize]) -> M {
    w.iter().fold(identity(), |a, &mu| mul(a, gamma(mu)))
}
fn trace(a: M) -> C {
    (0..4).fold(C::ZERO, |s, i| s + a[i][i])
}
fn metric_values(names: &[&str; 4]) -> BTreeMap<String, BigC> {
    let mut v = BTreeMap::new();
    for i in 0..4 {
        for j in i..4 {
            let (a, c) = if names[i] <= names[j] {
                (names[i], names[j])
            } else {
                (names[j], names[i])
            };
            v.insert(
                format!("h({}:{};{}:{})", a.len(), a, c.len(), c),
                BigC::from(if i != j {
                    0
                } else if i == 0 {
                    1
                } else {
                    -1
                }),
            );
        }
    }
    v
}
#[test]
fn ch35_sigma_contraction_all_components_independent_epsilon() {
    // (35.4), eta=(-,+,+,+); all 16 spinor index combinations.
    let eps = [[0, -1], [1, 0]];
    for a in 0..2 {
        for d in 0..2 {
            for c in 0..2 {
                for e in 0..2 {
                    let lhs = (0..4).fold(P::default(), |s, mu| {
                        s + P::constant(BigC::from(
                            sigma(mu)[a][d]
                                * sigma(mu)[c][e]
                                * C::from(if mu == 0 { -1 } else { 1 }),
                        ))
                    });
                    assert_eq!(lhs, k(-2 * eps[a][c] * eps[d][e]));
                }
            }
        }
    }
    assert_ne!(k(-2 * eps[0][1] * eps[0][1]), k(2));
    // (35.19)-(35.20): two epsilon raises include the index transpose.
    for profile in [E::Library, E::MonoNotes] {
        for mu in 0..4 {
            for d in 0..2 {
                for a in 0..2 {
                    let mut actual = P::default();
                    for c in 0..2 {
                        for e in 0..2 {
                            let epsilon = |i, j| {
                                if i == j {
                                    0
                                } else if i == 0 {
                                    profile.upper12()
                                } else {
                                    -profile.upper12()
                                }
                            };
                            actual = actual
                                + P::constant(BigC::from(sigma(mu)[c][e]))
                                    * k(epsilon(a, c) * epsilon(d, e));
                        }
                    }
                    assert_eq!(
                        actual,
                        P::constant(BigC::from(if mu == 0 {
                            sigma(mu)[d][a]
                        } else {
                            -sigma(mu)[d][a]
                        }))
                    );
                }
            }
        }
    }
}
#[test]
fn ch47_symbolic_four_trace_pairings_and_rename_control() {
    let mut b = Budget::default();
    for names in [["a", "b", "c", "d"], ["a;b", "β", "x:y", "renamed"]] {
        let h = |i, j, b: &mut Budget| cl::try_metric_with::<BigC>(names[i], names[j], b).unwrap();
        let expected = k(4)
            * (h(0, 1, &mut b) * h(2, 3, &mut b) - h(0, 2, &mut b) * h(1, 3, &mut b)
                + h(0, 3, &mut b) * h(1, 2, &mut b));
        assert_eq!(cl::try_trace_exact(&names, &mut b).unwrap(), expected);
        let wrong = k(4)
            * (h(0, 1, &mut b) * h(2, 3, &mut b)
                + h(0, 2, &mut b) * h(1, 3, &mut b)
                + h(0, 3, &mut b) * h(1, 2, &mut b));
        assert!(!expected.try_sub(&wrong, &mut b).unwrap().is_zero());
    }
    assert_ne!(
        cl::try_metric_with::<BigC>("a;b", "c", &mut b).unwrap(),
        cl::try_metric_with::<BigC>("a", "b;c", &mut b).unwrap()
    );
}
#[test]
fn ch47_component_trace_words_through_six_pauli_oracle() {
    let names = ["time", "x;y", "z:2", "β"];
    let values = metric_values(&names);

    // 1+4+16+64+256+1024+4096 = 5461 evaluated words, not symbolic proof.
    for n in 0..=6 {
        for code in 0..4usize.pow(n) {
            let w: Vec<_> = (0..n).map(|j| (code / 4usize.pow(j)) % 4).collect();
            let labels: Vec<_> = w.iter().map(|&i| names[i]).collect();
            let mut b = Budget::new(Limits {
                max_operations: 1_000_000,
                ..Limits::default()
            });
            assert_eq!(
                cl::try_trace_exact(&labels, &mut b)
                    .unwrap()
                    .try_eval(&values, &mut b)
                    .unwrap(),
                BigC::from(trace(word(&w))),
                "{w:?}"
            );
        }
    }
    assert_ne!(trace(word(&[0, 0])), C::from(-4));
}
#[test]
fn ch47_formal_dimension_and_component_sandwiches() {
    let mut b = Budget::default();
    let d = P::var("d");
    let one = cl::try_sandwich_srednicki_with::<BigC>(&["a"], &d, &mut b).unwrap();
    assert_eq!(one[&vec!["a".into()]], d - k(2)); // (47.19), formal dimension.
    let names = ["t", "x", "y", "z"];
    let values = metric_values(&names);
    for n in 0..=3 {
        for code in 0..4usize.pow(n) {
            let w: Vec<_> = (0..n).map(|j| (code / 4usize.pow(j)) % 4).collect();
            let labels: Vec<_> = w.iter().map(|&i| names[i]).collect();
            let normal = cl::try_sandwich_srednicki_with::<BigC>(&labels, &k(4), &mut b).unwrap();
            let mut actual = [[C::ZERO; 4]; 4];
            let mut expected = actual;
            for (slots, p) in normal {
                let indices: Vec<_> = slots
                    .iter()
                    .map(|s| names.iter().position(|x| x == s).unwrap())
                    .collect();
                let m = word(&indices);
                let c = p.try_eval(&values, &mut b).unwrap().try_to_small().unwrap();
                for i in 0..4 {
                    for j in 0..4 {
                        actual[i][j] = actual[i][j] + c * m[i][j];
                    }
                }
            }
            for mu in 0..4 {
                let m = mul(mul(gamma(mu), word(&w)), gamma(mu));
                for i in 0..4 {
                    for j in 0..4 {
                        expected[i][j] =
                            expected[i][j] + C::from(if mu == 0 { -1 } else { 1 }) * m[i][j];
                    }
                }
            }
            assert_eq!(actual, expected, "{w:?}");
        }
    }
    let empty = cl::try_sandwich_srednicki_with::<BigC>(&[], &k(5), &mut b).unwrap();
    assert_eq!(empty[&Vec::<String>::new()], k(-5));
    assert_ne!(empty[&Vec::<String>::new()], k(-4)); // finite oracle cannot realize d=5.
    assert!(cl::try_trace_exact(&["a"; 14], &mut b).is_err());
}
#[test]
fn ch50_schouten_symbolic_both_profiles_and_renamed_components() {
    let mut b = Budget::default();
    for prefix in ["p", "fresh:β;"] {
        let v: [[P; 2]; 4] =
            std::array::from_fn(|i| std::array::from_fn(|j| P::var(format!("{prefix}{i}_{j}"))));
        for profile in [E::Library, E::MonoNotes] {
            let mut br = |i, j| profile.bracket(&v[i], &v[j], &mut b).unwrap();
            let x = br(0, 1) * br(2, 3);
            let y = br(0, 2) * br(3, 1);
            let z = br(0, 3) * br(1, 2);
            assert!((x.clone() + y.clone() + z.clone()).is_zero());
            assert!(!(x - y + z).is_zero());
        }
    }
    let v = [[k(1), k(0)], [k(0), k(1)], [k(2), k(3)], [k(5), k(7)]];
    assert_eq!(bracket(&v[0], &v[1], &mut b), k(1));
    assert_eq!(bracket(&v[2], &v[3], &mut b), k(-1));
}
#[test]
fn ch50_exact_rational_boundary_and_invalid_domains() {
    let mut b = Budget::new(Limits {
        max_operations: 10_000_000,
        ..Limits::default()
    });
    let numerator = "340282366920938463463374607431768211457";
    let huge = BigQ::parse(numerator, "3", &mut b).unwrap();
    assert!(huge.try_to_small().is_err());
    let a = [P::constant(BigC::new(huge.clone(), BigQ::from(0))), k(0)];
    let c = [k(0), k(3)];
    let expected = BigC::new(BigQ::parse(numerator, "1", &mut b).unwrap(), BigQ::from(0));
    assert_eq!(
        bracket(&a, &c, &mut b)
            .try_eval(&BTreeMap::new(), &mut b)
            .unwrap(),
        expected
    );
    assert!(BigQ::parse("1", "0", &mut b).is_err());
    assert!(P::var("missing")
        .try_eval(&BTreeMap::new(), &mut b)
        .is_err());
    let x = Poly::var("x");
    let ratio = Ratio::new(x.clone(), x).unwrap();
    assert!(ratio
        .residual(&Ratio::new(Poly::constant(1), Poly::constant(1)).unwrap())
        .is_zero());
    assert!(ratio
        .eval(&BTreeMap::from([("x".into(), C::ZERO)]))
        .is_err());
    assert_eq!(
        ratio
            .eval(&BTreeMap::from([("x".into(), C::from(Q::new(1, 3)))]))
            .unwrap(),
        C::ONE
    );
    assert!(Ratio::new(Poly::constant(1), Poly::default()).is_err());
}

#[test]
fn ch50_massless_factorization_symbolic_and_rational_pauli_witness() {
    // (50.9), (50.12), (50.15): star slots are independent until evaluated.
    let [a, c, d, e] = ["phi1", "phi2", "star1", "star2"].map(P::var);
    let half = P::constant(BigC::from(Q::new(1, 2)));
    let inv_two_i = P::constant(BigC::from(C::new(Q::ZERO, Q::new(-1, 2))));
    let p = [
        half.clone() * (a.clone() * d.clone() + c.clone() * e.clone()),
        -half.clone() * (a.clone() * e.clone() + c.clone() * d.clone()),
        inv_two_i * (a.clone() * e.clone() - c.clone() * d.clone()),
        half * (c.clone() * e.clone() - a.clone() * d.clone()),
    ];
    assert!((-p[0].clone().pow(2)
        + p[1].clone().pow(2)
        + p[2].clone().pow(2)
        + p[3].clone().pow(2))
    .is_zero());
    let phi = [
        C::new(Q::new(1, 3), Q::new(2, 5)),
        C::new(Q::new(7, 11), Q::new(-3, 13)),
    ];
    let values = BTreeMap::from([
        ("phi1".into(), BigC::from(phi[0])),
        ("phi2".into(), BigC::from(phi[1])),
        ("star1".into(), BigC::from(phi[0].conj())),
        ("star2".into(), BigC::from(phi[1].conj())),
    ]);
    let mut b = Budget::default();
    let evaluated = p.map(|v| v.try_eval(&values, &mut b).unwrap().try_to_small().unwrap());
    let mut contraction = [[C::ZERO; 2]; 2];
    for (mu, value) in evaluated.iter().enumerate() {
        for (i, row) in contraction.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry = *entry + C::from(if mu == 0 { -1 } else { 1 }) * *value * sigma(mu)[i][j];
            }
        }
    }
    let outer = std::array::from_fn(|i| std::array::from_fn(|j| -phi[i] * phi[j].conj()));
    assert_eq!(contraction, outer);
    assert_ne!(
        contraction[0][1], outer[1][0],
        "complex witness detects a transposed or wrong sigma2 sign"
    );
}

#[test]
fn ch50_four_bracket_chiral_trace_rational_component_witness() {
    // (50.37). This is an evaluated witness, not an all-component identity.
    let phi = [
        [C::ONE, C::I],
        [C::from(2), C::ONE],
        [C::ONE + C::I, C::from(2)],
        [C::ONE, C::from(2) * C::I],
    ];
    let star = phi.map(|v| v.map(|z| z.conj()));
    let mut b = Budget::default();
    let a = phi.map(|v| v.map(|z| P::constant(BigC::from(z))));
    let s = star.map(|v| v.map(|z| P::constant(BigC::from(z))));
    let product = bracket(&s[0], &s[1], &mut b)
        * (-bracket(&a[1], &a[2], &mut b))
        * bracket(&s[2], &s[3], &mut b)
        * (-bracket(&a[3], &a[0], &mut b));
    // Independent inversion of the displayed Pauli block (50.15), in C.
    let slash: [M; 4] = std::array::from_fn(|leg| {
        let u = phi[leg];
        let v = star[leg];
        let half = C::from(Q::new(1, 2));
        let p = [
            half * (u[0] * v[0] + u[1] * v[1]),
            -half * (u[0] * v[1] + u[1] * v[0]),
            -C::I * half * (u[0] * v[1] - u[1] * v[0]),
            half * (u[1] * v[1] - u[0] * v[0]),
        ];
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..4).fold(C::ZERO, |sum, mu| {
                    sum + C::from(if mu == 0 { -1 } else { 1 }) * p[mu] * gamma(mu)[i][j]
                })
            })
        })
    });
    let matrix = slash.into_iter().fold(identity(), mul);
    let left = matrix[0][0] + matrix[1][1];
    let right = matrix[2][2] + matrix[3][3];
    assert_eq!(
        product.try_eval(&BTreeMap::new(), &mut b).unwrap(),
        BigC::from(left)
    );
    assert_ne!(
        left, right,
        "opposite chirality must be detected by this witness"
    );
    assert_ne!(left, C::ZERO);
}

#[test]
fn ch50_nonnull_determinant_sign_control() {
    // Explicit (50.15) at p=(2,1,0,0): determinant=3=-eta(p,p).
    // The following OCR prose says eta(p,p); do not inherit its sign blindly.
    let block = [[-2_i64, 1], [1, -2]];
    let det = block[0][0] * block[1][1] - block[0][1] * block[1][0];
    assert_eq!(P::constant(det), k(3));
    assert_ne!(P::constant(det), k(-3));
}
