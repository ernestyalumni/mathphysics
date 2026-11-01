use qft_symbolics::{
    computation::{Budget, Limits},
    exact::BigC,
    expression::*,
};
fn sym(n: &str) -> Expr {
    Expr::Symbol(n.into())
}
fn can(e: &Expr) -> String {
    e.canonical(&mut Budget::default()).unwrap()
}
fn atom(n: &str, label: IndexName, v: Variance) -> Expr {
    Expr::Tensor(TensorAtom {
        name: n.into(),
        slots: vec![(Space::Lorentz, v)],
        indices: vec![Index {
            name: label,
            space: Space::Lorentz,
            variance: v,
        }],
    })
}
fn contraction(n: &str) -> Expr {
    Expr::Contract {
        binder: n.into(),
        space: Space::Lorentz,
        body: Box::new(Expr::OrderedProduct(vec![
            atom("A", IndexName::Bound(n.into()), Variance::Up),
            atom("B", IndexName::Bound(n.into()), Variance::Down),
        ])),
    }
}
#[test]
fn scalar_commutes_operator_order_survives() {
    assert_eq!(
        can(&Expr::Product(vec![sym("a"), sym("b")])),
        can(&Expr::Product(vec![sym("b"), sym("a")]))
    );
    assert_eq!(
        can(&Expr::Sum(vec![sym("a"), sym("b")])),
        can(&Expr::Sum(vec![sym("b"), sym("a")]))
    );
    let a = Expr::Operator("a".into());
    let b = Expr::Operator("b".into());
    assert_ne!(
        can(&Expr::OrderedProduct(vec![a.clone(), b.clone()])),
        can(&Expr::OrderedProduct(vec![b, a]))
    );
}
#[test]
fn denominators_are_obligations_even_under_zero() {
    let q = Expr::Quotient(
        Box::new(Expr::Constant(BigC::default())),
        Box::new(sym("x")),
    );
    assert_eq!(
        q.nonzero_obligations(&mut Budget::default()).unwrap(),
        vec!["nonzero(qft-expr.v1:S(1:x))"]
    );
    assert!(can(&q).starts_with("qft-expr.v1:Q("));
    let e = Expr::Product(vec![Expr::Constant(BigC::default()), q]);
    assert_eq!(
        e.nonzero_obligations(&mut Budget::default()).unwrap().len(),
        1
    );
    let bad = Expr::Quotient(
        Box::new(sym("x")),
        Box::new(Expr::Constant(BigC::default())),
    );
    assert!(Expr::Product(vec![Expr::Constant(BigC::default()), bad])
        .canonical(&mut Budget::default())
        .is_err());
}
#[test]
fn tensor_validation_and_domain_rejection() {
    let a = atom("A", IndexName::Free("i".into()), Variance::Up);
    assert_eq!(
        a.validate(&mut Budget::default())
            .unwrap()
            .free_indices
            .len(),
        1
    );
    assert!(Expr::Product(vec![a.clone()])
        .validate(&mut Budget::default())
        .is_err());
    assert!(Expr::Sum(vec![a.clone(), sym("a")])
        .validate(&mut Budget::default())
        .is_err());
    assert!(Expr::Sum(vec![
        a.clone(),
        atom("B", IndexName::Free("i".into()), Variance::Down)
    ])
    .validate(&mut Budget::default())
    .is_err());
    assert!(Expr::OrderedProduct(vec![a.clone(), a.clone()])
        .validate(&mut Budget::default())
        .is_err());
    assert!(
        Expr::OrderedProduct(vec![a.clone(), Expr::Operator("X".into())])
            .validate(&mut Budget::default())
            .is_err()
    );
    assert!(
        Expr::OrderedProduct(vec![Expr::Operator("X".into()), a.clone()])
            .validate(&mut Budget::default())
            .is_err()
    );
    assert!(Expr::Quotient(Box::new(sym("x")), Box::new(a))
        .validate(&mut Budget::default())
        .is_err());
    let mut bad = match atom("A", IndexName::Free("i".into()), Variance::Up) {
        Expr::Tensor(t) => t,
        _ => unreachable!(),
    };
    bad.slots[0].0 = Space::Dotted;
    assert!(Expr::Tensor(bad).validate(&mut Budget::default()).is_err());
}
#[test]
fn alpha_hygiene_and_malformed_binders() {
    assert_eq!(can(&contraction("i")), can(&contraction("renamed")));
    assert_eq!(
        contraction("i")
            .validate(&mut Budget::default())
            .unwrap()
            .domain,
        Domain::Scalar
    );
    let dangling = atom("A", IndexName::Bound("i".into()), Variance::Up);
    assert!(dangling.validate(&mut Budget::default()).is_err());
    let once = Expr::Contract {
        binder: "i".into(),
        space: Space::Lorentz,
        body: Box::new(dangling),
    };
    assert!(once.validate(&mut Budget::default()).is_err());
    let shadow = Expr::Contract {
        binder: "i".into(),
        space: Space::Lorentz,
        body: Box::new(contraction("i")),
    };
    assert!(shadow.validate(&mut Budget::default()).is_err());
    let wrong = Expr::Contract {
        binder: "i".into(),
        space: Space::Dotted,
        body: Box::new(match contraction("i") {
            Expr::Contract { body, .. } => *body,
            _ => unreachable!(),
        }),
    };
    assert!(wrong.validate(&mut Budget::default()).is_err());
    // A free name is never captured by a binder with the same spelling.
    let free = Expr::OrderedProduct(vec![
        contraction("i"),
        atom("F", IndexName::Free("i".into()), Variance::Up),
    ]);
    assert_eq!(
        free.validate(&mut Budget::default())
            .unwrap()
            .free_indices
            .len(),
        1
    );
    assert_eq!(
        can(&Expr::Product(vec![contraction("i"), contraction("i")])),
        can(&Expr::Product(vec![contraction("j"), contraction("k")]))
    );
}
#[test]
fn branchwise_binder_counts_and_variance() {
    let Expr::Contract { body, .. } = contraction("i") else {
        unreachable!()
    };
    let e = Expr::Contract {
        binder: "i".into(),
        space: Space::Lorentz,
        body: Box::new(Expr::Sum(vec![*body.clone(), *body])),
    };
    assert!(e.validate(&mut Budget::default()).is_ok());
    let bad = Expr::Contract {
        binder: "i".into(),
        space: Space::Lorentz,
        body: Box::new(Expr::OrderedProduct(vec![
            atom("A", IndexName::Bound("i".into()), Variance::Up),
            atom("B", IndexName::Bound("i".into()), Variance::Up),
        ])),
    };
    assert!(bad.validate(&mut Budget::default()).is_err());
}
#[test]
fn limits_before_derived_allocation() {
    let e = Expr::Product(vec![sym("x"), sym("y")]);
    for limits in [
        Limits {
            max_terms: 2,
            ..Limits::default()
        },
        Limits {
            max_degree: 1,
            ..Limits::default()
        },
        Limits {
            max_input_bytes: 1,
            ..Limits::default()
        },
        Limits {
            max_operations: 1,
            ..Limits::default()
        },
    ] {
        assert!(e.canonical(&mut Budget::new(limits)).is_err());
    }
    let mut deep = sym("x");
    for _ in 0..257 {
        deep = Expr::Product(vec![deep]);
    }
    assert!(deep
        .validate(&mut Budget::new(Limits {
            max_degree: 1000,
            ..Limits::default()
        }))
        .is_err());
}
#[test]
fn names_are_unambiguous_and_serialization_is_nonmutating() {
    let e = sym("a),S(1:b");
    let original = e.clone();
    assert_ne!(can(&e), can(&Expr::Sum(vec![sym("a"), sym("b")])));
    assert_eq!(e, original);
    assert!(sym("").validate(&mut Budget::default()).is_err());
}

#[test]
fn rank_zero_tensor_is_not_an_operator() {
    let t = Expr::Tensor(TensorAtom {
        name: "t".into(),
        slots: vec![],
        indices: vec![],
    });
    for factors in [
        vec![t.clone(), Expr::Operator("O".into())],
        vec![Expr::Operator("O".into()), t],
    ] {
        assert!(Expr::OrderedProduct(factors)
            .validate(&mut Budget::default())
            .is_err());
    }
    assert!(can(&sym("x")).starts_with("qft-expr.v1:"));
}
