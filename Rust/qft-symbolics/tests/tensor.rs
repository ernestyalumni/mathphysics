use qft_symbolics::{
    polynomial::Poly,
    tensor::{
        Context, Dimension, Error, Expr, Factor as F, Index as I, Slot, Space as S, Symmetry, Term,
        Variance as V,
    },
};
fn up(n: &str) -> I {
    I::up(S::Lorentz, n)
}
fn dn(n: &str) -> I {
    I::down(S::Lorentz, n)
}
fn t(n: &str, i: Vec<I>) -> F {
    F::tensor(n, i)
}
fn term(c: i64, f: Vec<F>) -> Term {
    Term::new(Poly::constant(c), f)
}
fn expr(c: &Context, f: Vec<F>) -> Expr {
    Expr::from_terms(c, vec![term(1, f)]).unwrap()
}
fn ctx() -> Context {
    let mut c = Context::default();
    for name in ["A", "B", "C", "D", "V", "W"] {
        c.declare(name, vec![Slot::Lorentz], Symmetry::None)
            .unwrap();
    }
    for (name, sym) in [
        ("S", Symmetry::Symmetric),
        ("K", Symmetry::Antisymmetric),
        ("T", Symmetry::None),
    ] {
        c.declare(name, vec![Slot::Lorentz; 2], sym).unwrap();
    }
    c
}
fn equal(a: &Expr, b: &Expr) {
    assert!(a.residual(b).unwrap().is_zero(), "{a:?} != {b:?}");
}
#[test]
fn free_signature_and_typed_zero_cannot_be_discarded() {
    let c = ctx();
    let a = expr(&c, vec![t("A", vec![dn("mu")])]);
    let b = expr(&c, vec![t("B", vec![dn("mu")])]);
    assert_eq!(a.add(&b).unwrap().signature(), &[dn("mu")]);
    for b in [
        expr(&c, vec![t("B", vec![dn("nu")])]),
        expr(&c, vec![t("B", vec![up("mu")])]),
    ] {
        assert_eq!(a.add(&b), Err(Error::SignatureMismatch));
    }
    let zero = a.residual(&a).unwrap();
    assert!(zero.is_zero());
    assert_eq!(zero.signature(), &[dn("mu")]);
    let scalar_zero = Expr::from_terms(&c, vec![]).unwrap();
    assert_eq!(zero.add(&scalar_zero), Err(Error::SignatureMismatch));
    assert_eq!(
        zero.multiply(&b).unwrap_err(),
        Error::Invalid("same-variance repeated index".into())
    );
}
#[test]
fn malformed_indices_rejected_even_when_coefficient_is_zero() {
    let c = ctx();
    for factors in [
        vec![
            t("A", vec![up("a")]),
            t("B", vec![dn("a")]),
            t("C", vec![dn("a")]),
        ],
        vec![t("A", vec![up("a")]), t("B", vec![up("a")])],
        vec![t("A", vec![dn("a")]), t("B", vec![dn("a")])],
        vec![t("A", vec![I::down(S::Dotted, "a")])],
        vec![F::metric(dn("a"), up("b"))],
        vec![F::delta(up("a"), I::down(S::Undotted, "b"))],
        vec![t("undefined", vec![])],
    ] {
        assert!(Expr::from_terms(&c, vec![term(0, factors)]).is_err());
    }
    let invalid_sum = vec![
        term(0, vec![t("A", vec![dn("a")])]),
        term(0, vec![t("A", vec![dn("b")])]),
    ];
    assert_eq!(
        Expr::from_terms(&c, invalid_sum),
        Err(Error::SignatureMismatch)
    );
}
#[test]
fn same_printed_names_in_distinct_spaces_are_not_contractions() {
    let mut c = ctx();
    c.declare(
        "u",
        vec![Slot::FixedSpinor(S::Undotted, V::Down)],
        Symmetry::None,
    )
    .unwrap();
    c.declare(
        "v",
        vec![Slot::FixedSpinor(S::Dotted, V::Up)],
        Symmetry::None,
    )
    .unwrap();
    let a = expr(
        &c,
        vec![
            t("u", vec![I::down(S::Undotted, "a")]),
            t("v", vec![I::up(S::Dotted, "a")]),
        ],
    );
    assert_eq!(a.signature().len(), 2);
    assert!(!a.is_zero());
    let swapped = t("u", vec![I::down(S::Dotted, "a")]);
    assert!(Expr::from_terms(&c, vec![term(1, vec![swapped])]).is_err());
}
#[test]
fn lorentz_dotted_undotted_contractions_are_independent() {
    let mut c = ctx();
    c.declare(
        "sigma",
        vec![
            Slot::Lorentz,
            Slot::FixedSpinor(S::Undotted, V::Up),
            Slot::FixedSpinor(S::Dotted, V::Up),
        ],
        Symmetry::None,
    )
    .unwrap();
    c.declare(
        "u",
        vec![Slot::FixedSpinor(S::Undotted, V::Down)],
        Symmetry::None,
    )
    .unwrap();
    c.declare(
        "v",
        vec![Slot::FixedSpinor(S::Dotted, V::Down)],
        Symmetry::None,
    )
    .unwrap();
    let a = expr(
        &c,
        vec![
            t(
                "sigma",
                vec![dn("a"), I::up(S::Undotted, "a"), I::up(S::Dotted, "a")],
            ),
            t("V", vec![up("a")]),
            t("u", vec![I::down(S::Undotted, "a")]),
            t("v", vec![I::down(S::Dotted, "a")]),
        ],
    );
    assert!(a.signature().is_empty());
    assert!(!a.is_zero());
}
#[test]
fn symmetric_and_antisymmetric_slot_orbits_preserve_free_names() {
    let c = ctx();
    let s = expr(&c, vec![t("S", vec![dn("a"), dn("b")])]);
    equal(&s, &expr(&c, vec![t("S", vec![dn("b"), dn("a")])]));
    let a = expr(&c, vec![t("K", vec![dn("a"), dn("b")])]);
    let b = expr(&c, vec![t("K", vec![dn("b"), dn("a")])]);
    assert!(a.add(&b).unwrap().is_zero());
    assert!(!a.residual(&b).unwrap().is_zero());
    assert_eq!(
        s.residual(&expr(&c, vec![t("S", vec![dn("c"), dn("d")])])),
        Err(Error::SignatureMismatch)
    );
    // Full variance-bearing slots move under symmetry of metric-related views.
    equal(
        &expr(&c, vec![t("S", vec![up("a"), dn("b")])]),
        &expr(&c, vec![t("S", vec![dn("b"), up("a")])]),
    );
}
#[test]
fn opposite_sign_stabilizers_force_zero() {
    let c = ctx();
    assert!(expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("S", vec![up("a"), up("b")])
        ]
    )
    .is_zero());
    assert!(expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("V", vec![up("a")]),
            t("V", vec![up("b")])
        ]
    )
    .is_zero());
    assert!(!expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("V", vec![up("a")]),
            t("W", vec![up("b")])
        ]
    )
    .is_zero());
}
#[test]
fn dummy_renaming_factor_order_and_sign_must_be_canonicalised_jointly() {
    let c = ctx();
    let a = expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("V", vec![up("a")]),
            t("W", vec![up("b")]),
        ],
    );
    let b = expr(
        &c,
        vec![
            t("W", vec![up("r")]),
            t("K", vec![dn("r"), dn("s")]),
            t("V", vec![up("s")]),
        ],
    );
    assert!(a.add(&b).unwrap().is_zero());
    assert!(!a.residual(&b).unwrap().is_zero());
    // Same letter choices are not required.
    equal(
        &a,
        &expr(
            &c,
            vec![
                t("W", vec![up("z")]),
                t("V", vec![up("x")]),
                t("K", vec![dn("x"), dn("z")]),
            ],
        ),
    );
}
#[test]
fn capture_avoiding_multiplication_of_local_dummies() {
    let c = ctx();
    let p = expr(&c, vec![t("A", vec![up("mu")]), t("B", vec![dn("mu")])]);
    let q = expr(&c, vec![t("C", vec![up("mu")]), t("D", vec![dn("mu")])]);
    let product = p.multiply(&q).unwrap();
    equal(
        &product,
        &expr(
            &c,
            vec![
                t("A", vec![up("mu")]),
                t("B", vec![dn("mu")]),
                t("C", vec![up("nu")]),
                t("D", vec![dn("nu")]),
            ],
        ),
    );
    let free = expr(&c, vec![t("V", vec![dn("mu")])]);
    let mixed = free.multiply(&p).unwrap();
    assert_eq!(mixed.signature(), &[dn("mu")]);
    equal(
        &mixed,
        &expr(
            &c,
            vec![
                t("V", vec![dn("mu")]),
                t("A", vec![up("rho")]),
                t("B", vec![dn("rho")]),
            ],
        ),
    );
    assert!(Expr::from_terms(
        &c,
        vec![term(
            1,
            vec![
                t("A", vec![up("mu")]),
                t("B", vec![dn("mu")]),
                t("C", vec![up("mu")]),
                t("D", vec![dn("mu")])
            ]
        )]
    )
    .is_err());
}
#[test]
fn free_indices_contract_across_multiplication_but_internal_bound_names_do_not() {
    let c = ctx();
    let a = expr(&c, vec![t("A", vec![up("mu")])]);
    let b = expr(&c, vec![t("B", vec![dn("mu")])]);
    equal(
        &a.multiply(&b).unwrap(),
        &expr(&c, vec![t("A", vec![up("k")]), t("B", vec![dn("k")])]),
    );
    // User labels resembling generated names are still distinct from bound IDs.
    let free = expr(&c, vec![t("V", vec![up("_d0")])]);
    assert_eq!(
        free.multiply(&a.multiply(&b).unwrap()).unwrap().signature(),
        &[up("_d0")]
    );
}
#[test]
fn delta_substitution_and_chains_in_all_spaces() {
    let mut c = ctx();
    c.declare(
        "u",
        vec![Slot::FixedSpinor(S::Undotted, V::Up)],
        Symmetry::None,
    )
    .unwrap();
    c.declare(
        "v",
        vec![Slot::FixedSpinor(S::Dotted, V::Up)],
        Symmetry::None,
    )
    .unwrap();
    for (space, symbol) in [(S::Lorentz, "V"), (S::Undotted, "u"), (S::Dotted, "v")] {
        let a = I::up(space, "a");
        let b = I::down(space, "b");
        let upper_b = I::up(space, "b");
        equal(
            &expr(&c, vec![F::delta(a.clone(), b), t(symbol, vec![upper_b])]),
            &expr(&c, vec![t(symbol, vec![a])]),
        );
    }
    equal(
        &expr(
            &c,
            vec![F::delta(up("a"), dn("b")), F::delta(up("b"), dn("c"))],
        ),
        &expr(&c, vec![F::delta(up("a"), dn("c"))]),
    );
}
#[test]
fn traces_give_dimensions_not_metric_signature_sum() {
    let c = ctx();
    for (s, d) in [(S::Lorentz, 4), (S::Undotted, 2), (S::Dotted, 2)] {
        assert_eq!(
            expr(&c, vec![F::delta(I::up(s, "a"), I::down(s, "a"))]).scalar(),
            Some(Poly::constant(d))
        );
    }
    let trace = expr(
        &c,
        vec![F::metric(dn("a"), dn("b")), F::metric(up("a"), up("b"))],
    );
    assert_eq!(trace.scalar(), Some(Poly::constant(4)));
    assert_ne!(trace.scalar(), Some(Poly::constant(2)));
    let formal = c
        .with_lorentz_dimension(Dimension::Symbol("d".into()))
        .unwrap();
    assert_eq!(
        expr(
            &formal,
            vec![F::metric(dn("a"), dn("b")), F::metric(up("a"), up("b"))]
        )
        .scalar(),
        Some(Poly::var("d"))
    );
}
#[test]
fn symmetric_metric_raising_lowering_and_inverse_chains() {
    let c = ctx();
    equal(
        &expr(&c, vec![F::metric(dn("a"), dn("b")), t("V", vec![up("b")])]),
        &expr(&c, vec![t("V", vec![dn("a")])]),
    );
    equal(
        &expr(&c, vec![F::metric(up("a"), up("b")), t("V", vec![dn("b")])]),
        &expr(&c, vec![t("V", vec![up("a")])]),
    );
    equal(
        &expr(
            &c,
            vec![F::metric(dn("a"), dn("b")), F::metric(up("b"), up("c"))],
        ),
        &expr(&c, vec![F::delta(up("c"), dn("a"))]),
    );
    equal(
        &expr(
            &c,
            vec![
                F::metric(dn("a"), dn("b")),
                F::metric(up("b"), up("c")),
                t("V", vec![dn("c")]),
            ],
        ),
        &expr(&c, vec![t("V", vec![dn("a")])]),
    );
}
#[test]
fn spinor_epsilon_is_not_a_symmetric_metric() {
    let mut c = ctx();
    c.declare(
        "eps",
        vec![Slot::FixedSpinor(S::Undotted, V::Down); 2],
        Symmetry::Antisymmetric,
    )
    .unwrap();
    c.declare(
        "psi",
        vec![Slot::FixedSpinor(S::Undotted, V::Up)],
        Symmetry::None,
    )
    .unwrap();
    let a = expr(
        &c,
        vec![
            t(
                "eps",
                vec![I::down(S::Undotted, "a"), I::down(S::Undotted, "b")],
            ),
            t("psi", vec![I::up(S::Undotted, "b")]),
        ],
    );
    assert!(!a.is_zero());
    assert_eq!(a.term_count(), 1);
    assert_eq!(a.signature(), &[I::down(S::Undotted, "a")]);
    assert!(matches!(
        c.require_symmetric_metric(S::Undotted),
        Err(Error::Unsupported(_))
    ));
    assert!(Expr::from_terms(
        &c,
        vec![term(
            1,
            vec![F::metric(
                I::down(S::Undotted, "a"),
                I::down(S::Undotted, "b")
            )]
        )]
    )
    .is_err());
}
#[test]
fn invalid_symmetry_groups_and_contexts_are_not_accepted() {
    let mut c = ctx();
    for generators in [
        vec![(vec![0, 0], 1)],
        vec![(vec![1, 0], 0)],
        vec![(vec![1], 1)],
    ] {
        assert!(c
            .declare(
                "bad",
                vec![Slot::Lorentz; 2],
                Symmetry::Generators(generators)
            )
            .is_err());
    }
    assert!(c
        .declare(
            "bad",
            vec![Slot::Lorentz, Slot::FixedSpinor(S::Undotted, V::Up)],
            Symmetry::Symmetric
        )
        .is_err());
    assert!(c.declare("A", vec![], Symmetry::None).is_err());
    let a = expr(&c, vec![t("A", vec![dn("a")])]);
    let d = c.with_lorentz_dimension(Dimension::Fixed(3)).unwrap();
    assert_eq!(
        a.residual(&expr(&d, vec![t("A", vec![dn("a")])])),
        Err(Error::ContextMismatch)
    );
}
#[test]
fn bounded_orbit_reports_inconclusive_not_partial_representative() {
    let c = ctx().with_budget(2);
    assert_eq!(
        Expr::from_terms(
            &c,
            vec![term(
                1,
                vec![
                    t("K", vec![dn("a"), dn("b")]),
                    t("V", vec![up("a")]),
                    t("W", vec![up("b")])
                ]
            )]
        ),
        Err(Error::Budget)
    );
    let mut c = Context::default().with_budget(1);
    assert_eq!(
        c.declare("large", vec![Slot::Lorentz; 3], Symmetry::Symmetric),
        Err(Error::Budget)
    );
}
#[test]
fn high_symmetry_preparation_is_charged_before_cloning_variants() {
    let mut c = Context::default().with_budget(500_000);
    c.declare("High", vec![Slot::Lorentz; 8], Symmetry::Symmetric)
        .unwrap();
    let c = c.with_budget(2);
    let indices = (0..8).map(|i| dn(&format!("a{i}"))).collect();
    assert_eq!(
        Expr::from_terms(&c, vec![term(1, vec![t("High", indices)])]),
        Err(Error::Budget)
    );
}
#[test]
fn lorentz_dummy_variance_orientation_and_fixed_points() {
    let c = ctx();
    let a = expr(&c, vec![t("A", vec![up("a")]), t("B", vec![dn("a")])]);
    let b = expr(&c, vec![t("A", vec![dn("b")]), t("B", vec![up("b")])]);
    equal(&a, &b);
    // This is not a licence to change a free index's variance.
    assert_eq!(
        expr(&c, vec![t("V", vec![up("a")])]).residual(&expr(&c, vec![t("V", vec![dn("a")])])),
        Err(Error::SignatureMismatch)
    );
    let invariant = expr(
        &c,
        vec![
            t("K", vec![dn("a"), up("b")]),
            t("V", vec![up("a")]),
            t("W", vec![dn("b")]),
        ],
    );
    let raw = invariant
        .monomials()
        .map(|(f, p)| Term::new(p.clone(), f.to_vec()))
        .collect();
    equal(&invariant, &Expr::from_terms(&c, raw).unwrap());
}
#[test]
fn antisymmetric_form_self_contraction_is_not_zero() {
    let c = ctx();
    let norm = expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("K", vec![up("a"), up("b")]),
        ],
    );
    assert!(!norm.is_zero());
    let swapped = expr(
        &c,
        vec![
            t("K", vec![dn("a"), dn("b")]),
            t("K", vec![up("b"), up("a")]),
        ],
    );
    assert!(norm.add(&swapped).unwrap().is_zero());
}
#[test]
fn shared_cases_include_nonzero_controls() {
    let cases = qft_symbolics::tensor_cases::cases().unwrap();
    assert_eq!(cases.len(), 12);
    assert_eq!(cases.iter().filter(|c| !c.expected_zero).count(), 3);
    for c in cases {
        assert_eq!(c.residual.is_zero(), c.expected_zero, "{}", c.id);
    }
}
#[test]
fn general_signed_symmetry_generators() {
    let mut c = ctx();
    // Pair antisymmetries plus pair exchange, but NOT the multi-term Bianchi identity.
    c.declare(
        "R",
        vec![Slot::Lorentz; 4],
        Symmetry::Generators(vec![
            (vec![1, 0, 2, 3], -1),
            (vec![0, 1, 3, 2], -1),
            (vec![2, 3, 0, 1], 1),
        ]),
    )
    .unwrap();
    let a = expr(&c, vec![t("R", vec![dn("a"), dn("b"), dn("c"), dn("d")])]);
    let b = expr(&c, vec![t("R", vec![dn("c"), dn("d"), dn("a"), dn("b")])]);
    equal(&a, &b);
    let reversed = expr(&c, vec![t("R", vec![dn("b"), dn("a"), dn("c"), dn("d")])]);
    assert!(a.add(&reversed).unwrap().is_zero());
}
