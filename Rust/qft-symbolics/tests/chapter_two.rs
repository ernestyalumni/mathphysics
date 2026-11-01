//! Independent source-driven acceptance; see docs/CHAPTER-TWO.md for premises.
use qft_symbolics::{
    computation::{Budget, Limits},
    differential::DifferentialOperator as D,
    exact::{BigC, BigQ},
    ordered::OrderedPolynomial as O,
    polynomial::BigPoly as P,
    relativity::{self as rel, Axis, Rapidity, Tolerance, UnitAxis},
};
use std::collections::BTreeMap;
const ETA: [i64; 4] = [-1, 1, 1, 1];
const PAIRS: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
fn budget() -> Budget {
    Budget::new(Limits {
        max_operations: 100_000_000,
        ..Limits::default()
    })
}
fn c(n: i64) -> P {
    P::try_constant(BigC::from(n), &mut budget()).unwrap()
}
fn v(s: &str) -> P {
    P::try_var(s, &mut budget()).unwrap()
}
fn add(a: &P, b: &P) -> P {
    a.try_add(b, &mut budget()).unwrap()
}
fn mul(a: &P, b: &P) -> P {
    a.try_mul(b, &mut budget()).unwrap()
}
fn neg(a: &P) -> P {
    a.try_neg(&mut budget()).unwrap()
}
fn imaginary(n: i64) -> P {
    P::try_constant(BigC::new(0.into(), n.into()), &mut budget()).unwrap()
}
fn og(s: &str) -> O {
    O::generator(s, &mut budget()).unwrap()
}
fn oa(a: &O, b: &O) -> O {
    a.try_add(b, &mut budget()).unwrap()
}
fn os(a: &O, b: &O) -> O {
    a.try_sub(b, &mut budget()).unwrap()
}
fn om(a: &O, b: &O) -> O {
    a.try_mul(b, &mut budget()).unwrap()
}
fn scale(a: &O, p: &P) -> O {
    a.try_scale(p, &mut budget()).unwrap()
}
fn comm(a: &O, b: &O) -> O {
    a.commutator(b, &mut budget()).unwrap()
}
fn m(mu: usize, nu: usize) -> O {
    if mu == nu {
        O::default()
    } else if mu < nu {
        og(&format!("M{mu}{nu}"))
    } else {
        scale(&m(nu, mu), &c(-1))
    }
}
fn p(mu: usize) -> O {
    og(&format!("P{mu}"))
}
fn metric(mu: usize, nu: usize) -> i64 {
    if mu == nu {
        ETA[mu]
    } else {
        0
    }
}
// Source (2.16), expanded independently; no production structure constants.
fn mm_rhs(mu: usize, nu: usize, r: usize, s: usize) -> O {
    let raw = oa(
        &os(
            &scale(&m(nu, s), &c(metric(mu, r))),
            &scale(&m(mu, s), &c(metric(nu, r))),
        ),
        &os(
            &scale(&m(mu, r), &c(metric(nu, s))),
            &scale(&m(nu, r), &c(metric(mu, s))),
        ),
    );
    scale(&raw, &mul(&imaginary(1), &v("hbar")))
}
fn pm_rhs(mu: usize, r: usize, s: usize) -> O {
    scale(
        &os(
            &scale(&p(r), &c(metric(mu, s))),
            &scale(&p(s), &c(metric(mu, r))),
        ),
        &mul(&imaginary(1), &v("hbar")),
    )
}
fn omega(mu: usize, nu: usize) -> P {
    if mu == nu {
        c(0)
    } else if mu < nu {
        v(&format!("w{mu}{nu}"))
    } else {
        neg(&omega(nu, mu))
    }
}
fn coeff(a: &O, s: &str) -> O {
    a.coefficient(s, 1, &mut budget()).unwrap()
}

#[test]
fn problem_2_1_metric_first_order_forces_lowered_antisymmetry() {
    // Compute epsilon coefficient from all sixteen unrestricted mixed entries.
    for (r, _) in ETA.iter().enumerate() {
        for (s, _) in ETA.iter().enumerate() {
            let mut first = c(0);
            for (mu, _) in ETA.iter().enumerate() {
                for (nu, _) in ETA.iter().enumerate() {
                    let wmr = v(&format!("mixed_{mu}_{r}"));
                    let wns = v(&format!("mixed_{nu}_{s}"));
                    let linear = add(
                        &mul(&c(i64::from(mu == r)), &wns),
                        &mul(&wmr, &c(i64::from(nu == s))),
                    );
                    first = add(&first, &mul(&c(metric(mu, nu)), &linear));
                }
            }
            let lowered_sum = add(
                &mul(&c(ETA[r]), &v(&format!("mixed_{r}_{s}"))),
                &mul(&c(ETA[s]), &v(&format!("mixed_{s}_{r}"))),
            );
            assert_eq!(first, lowered_sum);
            assert!(!first.is_zero()); // metric invariance is a real constraint.
            assert!(add(&omega(r, s), &omega(s, r)).is_zero());
        }
    }
    // Mixed boost entries are symmetric; lowered entries are antisymmetric.
    assert!(add(&mul(&c(-1), &v("boost")), &v("boost")).is_zero());
    assert!(!add(&mul(&c(-1), &v("boost")), &neg(&v("boost"))).is_zero());
}

#[test]
fn problem_2_2_conditional_antisymmetric_coefficient_extraction() {
    // Premise: exact group law and metric invariance give contracted (2.13).
    // Any Lambda yields antisymmetric transformed M; Lorentzness needed upstream.
    let lambda: [[P; 4]; 4] =
        std::array::from_fn(|i| std::array::from_fn(|j| v(&format!("Lambda{i}{j}"))));
    let tensor = |mu: usize, nu: usize| {
        let mut out = O::default();
        for (r, _) in ETA.iter().enumerate() {
            for (s, _) in ETA.iter().enumerate() {
                out = oa(&out, &scale(&m(r, s), &mul(&lambda[mu][r], &lambda[nu][s])));
            }
        }
        out
    };
    let mut contraction = O::default();
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            contraction = oa(&contraction, &scale(&tensor(mu, nu), &omega(mu, nu)));
            assert_eq!(oa(&tensor(mu, nu), &tensor(nu, mu)), O::default());
        }
    }
    for (mu, nu) in PAIRS {
        assert_eq!(
            coeff(&contraction, &format!("w{mu}{nu}")),
            scale(&tensor(mu, nu), &c(2))
        );
    }
    // An arbitrary symmetric coefficient is invisible: cannot extract it.
    assert!(add(
        &mul(&omega(0, 1), &v("symmetric")),
        &mul(&omega(1, 0), &v("symmetric"))
    )
    .is_zero());
    assert_ne!(og("A01"), og("A10"));
}

#[test]
fn problem_2_3_linearized_tensor_matches_commutator_coefficients() {
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            let mut delta = O::default();
            for (r, _) in ETA.iter().enumerate() {
                delta = oa(&delta, &scale(&m(r, nu), &mul(&c(ETA[mu]), &omega(mu, r))));
                delta = oa(&delta, &scale(&m(mu, r), &mul(&c(ETA[nu]), &omega(nu, r))));
            }
            for (r, s) in PAIRS {
                // U^-1 M U=M+(i/hbar)w_rs[M,M^rs]. Solve using hbar/i=-i hbar.
                let derived = scale(
                    &coeff(&delta, &format!("w{r}{s}")),
                    &mul(&imaginary(-1), &v("hbar")),
                );
                assert_eq!(derived, mm_rhs(mu, nu, r, s));
            }
        }
    }
    assert_ne!(mm_rhs(0, 1, 0, 2), scale(&mm_rhs(0, 1, 0, 2), &c(-1)));
}
#[test]
fn problem_2_5_linearized_vector_matches_commutator_coefficients() {
    for (mu, _) in ETA.iter().enumerate() {
        let mut delta = O::default();
        for (r, _) in ETA.iter().enumerate() {
            delta = oa(&delta, &scale(&p(r), &mul(&c(ETA[mu]), &omega(mu, r))));
        }
        for (r, s) in PAIRS {
            let derived = scale(
                &coeff(&delta, &format!("w{r}{s}")),
                &mul(&imaginary(-1), &v("hbar")),
            );
            assert_eq!(derived, pm_rhs(mu, r, s));
        }
    }
}
fn coords() -> [String; 4] {
    ["tau", "east", "north", "up"].map(str::to_owned)
}
fn orbital(mu: usize, nu: usize) -> D<4> {
    D::orbital_lorentz(coords(), mu, nu, 1.into(), &mut budget()).unwrap()
}
fn momentum(mu: usize) -> D<4> {
    D::momentum(coords(), mu, 1.into(), &mut budget()).unwrap()
}
fn ds(a: &D<4>, n: i64) -> D<4> {
    a.try_scale(&c(n), &mut budget()).unwrap()
}
fn di(a: &D<4>, n: i64) -> D<4> {
    a.try_scale(&imaginary(n), &mut budget()).unwrap()
}
fn da(a: &D<4>, b: &D<4>) -> D<4> {
    a.try_add(b, &mut budget()).unwrap()
}
fn dc(a: &D<4>, b: &D<4>) -> D<4> {
    a.commutator(b, &mut budget()).unwrap()
}
fn zero() -> D<4> {
    D::from_terms(coords(), BTreeMap::new(), &mut budget()).unwrap()
}
fn j(i: usize) -> D<4> {
    match i {
        0 => orbital(2, 3),
        1 => orbital(3, 1),
        _ => orbital(1, 2),
    }
}
fn k(i: usize) -> D<4> {
    orbital(i + 1, 0)
}
fn eps(i: usize, j: usize, k: usize) -> i64 {
    if i == j || i == k || j == k {
        0
    } else if (i, j, k) == (0, 1, 2) || (i, j, k) == (1, 2, 0) || (i, j, k) == (2, 0, 1) {
        1
    } else {
        -1
    }
}
#[test]
fn problems_2_4_and_2_6_all_j_k_h_p_relations_in_differential_representation() {
    // Differential representation witness: hbar=1 explicitly, K_i=M^{i0}.
    for i in 0..3 {
        assert_eq!(dc(&j(i), &momentum(0)), zero());
        assert_eq!(dc(&k(i), &momentum(0)), di(&momentum(i + 1), 1));
        for l in 0..3 {
            let mut jj = zero();
            let mut jk = zero();
            let mut jp = zero();
            for n in 0..3 {
                jj = da(&jj, &di(&j(n), eps(i, l, n)));
                jk = da(&jk, &di(&k(n), eps(i, l, n)));
                jp = da(&jp, &di(&momentum(n + 1), eps(i, l, n)));
            }
            assert_eq!(dc(&j(i), &j(l)), jj);
            assert_eq!(dc(&j(i), &k(l)), jk);
            assert_eq!(dc(&k(i), &k(l)), ds(&jj, -1));
            assert_eq!(dc(&j(i), &momentum(l + 1)), jp);
            assert_eq!(
                dc(&k(i), &momentum(l + 1)),
                di(&momentum(0), i64::from(i == l))
            );
        }
    }
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            assert_eq!(dc(&momentum(mu), &momentum(nu)), zero());
        }
    }
    assert_ne!(dc(&k(0), &k(1)), di(&j(2), 1));
    assert_ne!(dc(&k(0), &momentum(0)), di(&momentum(1), -1));
}
#[test]
fn problem_2_7_exact_translation_group_law_requires_commuting_mixed_coefficients() {
    // Independent mixed coefficient: T(a)T(b) gives AB; T(a+b) gives (AB+BA)/2.
    let a = og("a_dot_P");
    let b = og("b_dot_P");
    let product = om(&a, &b);
    let half = P::try_constant(
        BigC::new(BigQ::from(qft_symbolics::scalar::Q::new(1, 2)), 0.into()),
        &mut budget(),
    )
    .unwrap();
    let exponential = scale(&oa(&om(&a, &b), &om(&b, &a)), &half);
    assert_eq!(scale(&os(&product, &exponential), &c(2)), comm(&a, &b));
    assert_ne!(product, exponential); // No translation axiom in the free algebra.
                                      // Central bracket counterexample before imposing exact group law: [d_x,x]=I.
    let derivative = D::partial(coords(), 1, &mut budget()).unwrap();
    let x = D::multiplication(coords(), v("east"), &mut budget()).unwrap();
    let identity = D::multiplication(coords(), c(1), &mut budget()).unwrap();
    assert_eq!(dc(&derivative, &x), identity);
    assert_eq!(dc(&identity, &x), zero());
    assert_eq!(dc(&identity, &derivative), zero());
}
#[test]
fn problem_2_8_a_orbital_sign_and_exact_normal_form_with_hbar() {
    let actual = D::orbital_lorentz(coords(), 3, 0, BigQ::from(3), &mut budget()).unwrap();
    // L30=-3i(x3*(-d0)-x0*d3)=+3i(x3*d0+x0*d3).
    let expected = BTreeMap::from([
        ([1, 0, 0, 0], mul(&imaginary(3), &v("up"))),
        ([0, 0, 0, 1], mul(&imaginary(3), &v("tau"))),
    ]);
    assert_eq!(actual.terms(), &expected);
    let field = mul(&v("tau"), &v("up"));
    let image = actual.apply(&field, &mut budget()).unwrap();
    assert_eq!(
        image,
        mul(
            &imaginary(3),
            &add(&mul(&v("up"), &v("up")), &mul(&v("tau"), &v("tau")))
        )
    );
    assert_ne!(image, neg(&image));
}
#[test]
fn problem_2_8_b_ordered_double_action_includes_product_rule() {
    // Assumes coordinate coefficients commute with Hilbert-space generators M.
    let a = orbital(0, 1);
    let b = orbital(1, 2);
    let ab = a.compose(&b, &mut budget()).unwrap();
    let field = mul(&mul(&v("east"), &v("east")), &v("north"));
    assert_eq!(
        ab.apply(&field, &mut budget()).unwrap(),
        a.apply(&b.apply(&field, &mut budget()).unwrap(), &mut budget())
            .unwrap()
    );
    assert_ne!(ab, b.compose(&a, &mut budget()).unwrap());
    // Independent normal form: (-i(t d_e+e d_t))(-i(e d_n-n d_e)).
    let expected = BTreeMap::from([
        ([0, 0, 1, 0], neg(&v("tau"))),
        ([0, 1, 1, 0], neg(&mul(&v("tau"), &v("east")))),
        ([0, 2, 0, 0], mul(&v("tau"), &v("north"))),
        ([1, 0, 1, 0], neg(&mul(&v("east"), &v("east")))),
        ([1, 1, 0, 0], mul(&v("east"), &v("north"))),
    ]);
    assert_eq!(ab.terms(), &expected);
}
#[test]
fn problem_2_8_c_jacobi_in_free_associative_algebra() {
    let a = og("arbitrary_A");
    let b = og("arbitrary_B");
    let c = og("arbitrary_C");
    assert_eq!(
        oa(
            &oa(&comm(&comm(&a, &b), &c), &comm(&comm(&b, &c), &a)),
            &comm(&comm(&c, &a), &b)
        ),
        O::default()
    );
    let abc = om(&om(&a, &b), &c);
    assert_eq!(abc, om(&a, &om(&b, &c)));
    assert_ne!(abc, om(&om(&b, &a), &c));
}
#[test]
fn problem_2_8_d_nested_commutator_preserves_action_order() {
    let phi = og("phi");
    let a = og("M_a");
    let b = og("M_b");
    assert_eq!(
        comm(&phi, &comm(&a, &b)),
        os(&comm(&comm(&phi, &a), &b), &comm(&comm(&phi, &b), &a))
    );
    assert_ne!(
        comm(&phi, &comm(&a, &b)),
        os(&comm(&comm(&phi, &b), &a), &comm(&comm(&phi, &a), &b))
    );
}
#[test]
fn problem_2_8_e_orbital_brackets_all_indices_exact_normal_forms() {
    // hbar=2, explicitly. All 256 ordered index quadruples, including zeros.
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            for (r, _) in ETA.iter().enumerate() {
                for (s, _) in ETA.iter().enumerate() {
                    let l = |a, b| independent_orbital(a, b, 2);
                    let raw = da(
                        &da(
                            &ds(&l(nu, s), metric(mu, r)),
                            &ds(&l(mu, s), -metric(nu, r)),
                        ),
                        &da(
                            &ds(&l(mu, r), metric(nu, s)),
                            &ds(&l(nu, r), -metric(mu, s)),
                        ),
                    );
                    assert_eq!(
                        dc(&l(mu, nu), &l(r, s)),
                        di(&raw, 2),
                        "indices {mu}{nu},{r}{s}"
                    );
                }
            }
        }
    }
}
#[test]
fn problem_2_8_f_central_ambiguity_survives_field_actions() {
    let identity = D::multiplication(coords(), c(7), &mut budget()).unwrap();
    assert_ne!(identity, zero());
    // Scalar field and its derivative are represented as multiplication here.
    let field = mul(&mul(&v("tau"), &v("east")), &v("east"));
    let field_op = D::multiplication(coords(), field.clone(), &mut budget()).unwrap();
    assert_eq!(dc(&identity, &field_op), zero());
    for coordinate in coords() {
        let derivative =
            qft_symbolics::differential::differentiate(&field, &coordinate, &mut budget()).unwrap();
        let derivative_op = D::multiplication(coords(), derivative, &mut budget()).unwrap();
        assert_eq!(dc(&identity, &derivative_op), zero());
    }
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            assert_eq!(dc(&identity, &orbital(mu, nu)), zero());
            assert_eq!(dc(&identity, &momentum(mu)), zero());
        }
    }
    let r = dc(&orbital(0, 1), &orbital(0, 2));
    let shifted = da(&r, &identity);
    for (mu, _) in ETA.iter().enumerate() {
        assert_eq!(dc(&momentum(mu), &r), dc(&momentum(mu), &shifted));
    }
}
#[test]
fn reusable_rapidity_positive_boost_group_and_rejected_domains() {
    let tolerance = Tolerance {
        absolute: 1e-12,
        relative: 1e-12,
    };
    for axis in [Axis::X, Axis::Y, Axis::Z] {
        let r = Rapidity::from_beta(0.6).unwrap();
        let boost = rel::boost_f64(r, UnitAxis::axial(axis)).unwrap();
        assert!(
            rel::metric_residual_f64(&boost, tolerance)
                .unwrap()
                .within_tolerance
        );
        let moving = rel::apply_f64(&boost, &[1., 0., 0., 0.]).unwrap();
        assert!((moving[0] - 1.25).abs() < 1e-12);
        assert!((moving[1] + moving[2] + moving[3] - 0.75).abs() < 1e-12);
        let restored = rel::apply_f64(
            &rel::boost_f64(r.inverse(), UnitAxis::axial(axis)).unwrap(),
            &moving,
        )
        .unwrap();
        for (got, want) in restored.iter().zip([1., 0., 0., 0.]) {
            assert!((got - want).abs() < 1e-12);
        }
    }
    let a = Rapidity::from_beta(0.3).unwrap();
    let b = Rapidity::from_beta(-0.2).unwrap();
    let sum = a.collinear_add(b).unwrap();
    assert!((sum.beta().unwrap() - (0.3 - 0.2) / (1. - 0.06)).abs() < 1e-12);
    let vector = [2., 3., -5., 7.];
    let first = rel::boost_f64(a, UnitAxis::axial(Axis::Z)).unwrap();
    let second = rel::boost_f64(b, UnitAxis::axial(Axis::Z)).unwrap();
    let combined = rel::boost_f64(sum, UnitAxis::axial(Axis::Z)).unwrap();
    let sequential = rel::apply_f64(&first, &rel::apply_f64(&second, &vector).unwrap()).unwrap();
    let direct = rel::apply_f64(&combined, &vector).unwrap();
    for i in 0..4 {
        assert!((sequential[i] - direct[i]).abs() < 1e-12);
    }
    assert!(
        (rel::dot_f64(&vector, &vector).unwrap() - rel::dot_f64(&direct, &direct).unwrap()).abs()
            < 1e-10
    );
    for beta in [-1., 1., f64::INFINITY, f64::NAN] {
        assert!(Rapidity::from_beta(beta).is_err());
    }
    assert!(Rapidity::new(f64::NAN).is_err());
    assert!(UnitAxis::new([0., 0., 0.], 1e-12).is_err());
}

#[test]
fn problems_2_4_and_2_6_general_component_contractions_symbolic_hbar() {
    // Consequence of (2.16)/(2.18), not a claim these input laws were proved.
    let jpairs = [(2, 3), (3, 1), (1, 2)];
    let ih = mul(&imaginary(1), &v("hbar"));
    for i in 0..3 {
        let (a, b) = jpairs[i];
        assert_eq!(pm_rhs(0, a, b), O::default()); // [J,H]=-[H,J]=0.
        assert_eq!(scale(&pm_rhs(0, i + 1, 0), &c(-1)), scale(&p(i + 1), &ih));
        for l in 0..3 {
            let (r, s) = jpairs[l];
            let mut jj = O::default();
            let mut jk = O::default();
            let mut jp = O::default();
            for (n, &(u, w)) in jpairs.iter().enumerate() {
                let factor = mul(&ih, &c(eps(i, l, n)));
                jj = oa(&jj, &scale(&m(u, w), &factor));
                jk = oa(&jk, &scale(&m(n + 1, 0), &factor));
                jp = oa(&jp, &scale(&p(n + 1), &factor));
            }
            assert_eq!(mm_rhs(a, b, r, s), jj);
            assert_eq!(mm_rhs(a, b, l + 1, 0), jk);
            assert_eq!(mm_rhs(i + 1, 0, l + 1, 0), scale(&jj, &c(-1)));
            assert_eq!(scale(&pm_rhs(l + 1, a, b), &c(-1)), jp);
            assert_eq!(
                scale(&pm_rhs(l + 1, i + 1, 0), &c(-1)),
                scale(&p(0), &mul(&ih, &c(i64::from(i == l))))
            );
        }
    }
}

#[test]
fn problem_2_8_a_scalar_covariance_coefficient_derivation() {
    // phi((1-w)x)=phi(x)-w^r_s x^s partial_r phi, a conditional Taylor step.
    let mut delta = O::default();
    for (r, _) in ETA.iter().enumerate() {
        for (s, _) in ETA.iter().enumerate() {
            let coefficient = mul(&mul(&c(-ETA[r]), &omega(r, s)), &v(&format!("x{s}")));
            delta = oa(
                &delta,
                &scale(&og(&format!("partial_{r}_phi")), &coefficient),
            );
        }
    }
    for (mu, nu) in PAIRS {
        let derived = scale(
            &coeff(&delta, &format!("w{mu}{nu}")),
            &mul(&imaginary(-1), &v("hbar")),
        );
        let expected = os(
            &scale(
                &og(&format!("partial_{nu}_phi")),
                &mul(&c(ETA[nu]), &v(&format!("x{mu}"))),
            ),
            &scale(
                &og(&format!("partial_{mu}_phi")),
                &mul(&c(ETA[mu]), &v(&format!("x{nu}"))),
            ),
        );
        assert_eq!(derived, scale(&expected, &mul(&imaginary(-1), &v("hbar"))));
    }
}

fn independent_orbital(mu: usize, nu: usize, h: i64) -> D<4> {
    // Direct sparse map from definition, independent of orbital constructor.
    let mut terms = BTreeMap::new();
    if mu != nu {
        let mut dnu = [0; 4];
        dnu[nu] = 1;
        let mut dmu = [0; 4];
        dmu[mu] = 1;
        terms.insert(dnu, mul(&imaginary(-h * ETA[nu]), &v(&coords()[mu])));
        terms.insert(dmu, mul(&imaginary(h * ETA[mu]), &v(&coords()[nu])));
    }
    D::from_terms(coords(), terms, &mut budget()).unwrap()
}
#[test]
fn orbital_constructor_all_components_against_independent_sparse_maps() {
    for (mu, _) in ETA.iter().enumerate() {
        for (nu, _) in ETA.iter().enumerate() {
            assert_eq!(
                D::orbital_lorentz(coords(), mu, nu, 3.into(), &mut budget()).unwrap(),
                independent_orbital(mu, nu, 3)
            );
        }
    }
}

#[test]
fn first_order_conjugation_is_ordered_coefficient_identity() {
    // hbar=3 here. Independent truncated multiplication, no commutator axiom.
    let identity = O::scalar(c(1), &mut budget()).unwrap();
    let a = og("generic_A");
    let b = og("generic_generator");
    let third_i = P::try_constant(
        BigC::new(0.into(), qft_symbolics::scalar::Q::new(1, 3).into()),
        &mut budget(),
    )
    .unwrap();
    let perturb = scale(&b, &mul(&v("epsilon"), &third_i));
    let forward = oa(&identity, &perturb);
    let inverse = os(&identity, &perturb);
    let expanded = om(&om(&inverse, &a), &forward)
        .truncate(&["epsilon"], 1, &mut budget())
        .unwrap();
    assert_eq!(coeff(&expanded, "epsilon"), scale(&comm(&a, &b), &third_i));
    assert_ne!(coeff(&expanded, "epsilon"), scale(&comm(&b, &a), &third_i));
}

#[test]
fn translation_unitarity_alone_does_not_force_commutativity() {
    use qft_symbolics::{matrices as mat, scalar::C};
    // Exact finite-dimensional unitaries (also Hermitian): X and Z.
    let x = [[C::ZERO, C::ONE], [C::ONE, C::ZERO]];
    let z = [[C::ONE, C::ZERO], [C::ZERO, -C::ONE]];
    assert_eq!(mat::mul(x, x), mat::identity());
    assert_eq!(mat::mul(z, z), mat::identity());
    assert_ne!(mat::mul(x, z), mat::mul(z, x));
    assert_eq!(mat::mul(x, z), mat::scale(mat::mul(z, x), -C::ONE));
    // Anticommutation is a projective phase witness; no continuous group inferred.
}

#[test]
fn reusable_rapidity_exact_finite_witness_and_conditional_symbolic_metric() {
    use qft_symbolics::scalar::Q;
    // u=tanh(eta/2)=1/2 gives c=5/3,s=4/3 by rational identities.
    let boost =
        rel::exact_half_velocity_boost(&BigQ::from(Q::new(1, 2)), Axis::Z, &mut budget()).unwrap();
    let rat = |n, d| BigC::new(BigQ::from(Q::new(n, d)), 0.into());
    let expected = [
        [rat(5, 3), rat(0, 1), rat(0, 1), rat(4, 3)],
        [rat(0, 1), rat(1, 1), rat(0, 1), rat(0, 1)],
        [rat(0, 1), rat(0, 1), rat(1, 1), rat(0, 1)],
        [rat(4, 3), rat(0, 1), rat(0, 1), rat(5, 3)],
    ];
    assert_eq!(boost, expected);
    // Direct exact interval, independently expected -1 (not a production dot oracle).
    let image = rel::apply(
        &boost,
        &[rat(1, 1), rat(0, 1), rat(0, 1), rat(0, 1)],
        &mut budget(),
    )
    .unwrap();
    assert_eq!(image, [rat(5, 3), rat(0, 1), rat(0, 1), rat(4, 3)]);
    assert_eq!(rel::dot(&image, &image, &mut budget()).unwrap(), rat(-1, 1));
    let symbolic =
        rel::symbolic_axial_boost("hyperbolic_c", "hyperbolic_s", Axis::Z, &mut budget()).unwrap();
    let condition = add(
        &add(
            &mul(&v("hyperbolic_c"), &v("hyperbolic_c")),
            &neg(&mul(&v("hyperbolic_s"), &v("hyperbolic_s"))),
        ),
        &c(-1),
    );
    assert_eq!(symbolic.condition, condition);
    assert!(!symbolic.condition.is_zero());
    // c,s independent: the returned condition MUST NOT silently be erased.
    let mut bad = rel::boost_f64(Rapidity::new(0.4).unwrap(), UnitAxis::axial(Axis::Z)).unwrap();
    bad[0][0] += 0.1;
    assert!(
        !rel::metric_residual_f64(
            &bad,
            Tolerance {
                absolute: 1e-12,
                relative: 1e-12
            }
        )
        .unwrap()
        .within_tolerance
    );
}
