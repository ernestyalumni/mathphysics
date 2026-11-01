use qft_symbolics::{
    computation::{Budget, Error, Limits, Resource},
    exact::{BigC, BigQ},
    expression::Expr,
    replay::*,
};
fn c(n: i64) -> Expr {
    Expr::Constant(BigC::from(n))
}
fn s(n: &str) -> Expr {
    Expr::Symbol(n.into())
}
fn sum(v: Vec<Expr>) -> Expr {
    Expr::Sum(v)
}
fn prod(v: Vec<Expr>) -> Expr {
    Expr::Product(v)
}
fn div(a: Expr, b: Expr) -> Expr {
    Expr::Quotient(Box::new(a), Box::new(b))
}
fn cert(a: Expr, z: Expr, r: Rule) -> Certificate {
    Certificate {
        version: CERTIFICATE_VERSION,
        domain: SCALAR_DOMAIN.into(),
        start: a.clone(),
        end: z.clone(),
        steps: vec![ReplayStep {
            before: a,
            after: z,
            rule: r,
        }],
        nonzero_assumptions: vec![],
    }
}
fn check(c: &Certificate) -> qft_symbolics::computation::Result<ReplayReport> {
    replay(c, &mut Budget::default())
}
#[test]
fn independent_binomial_and_wrong_middle_coefficient() {
    let a = prod(vec![sum(vec![s("x"), c(2)]), sum(vec![s("x"), c(2)])]);
    let z = sum(vec![
        prod(vec![s("x"), s("x")]),
        prod(vec![c(4), s("x")]),
        c(4),
    ]);
    assert_eq!(
        check(&cert(a.clone(), z, Rule::PolynomialIdentity))
            .unwrap()
            .steps_checked,
        1
    );
    let wrong = sum(vec![
        prod(vec![s("x"), s("x")]),
        prod(vec![c(3), s("x")]),
        c(4),
    ]);
    assert!(check(&cert(a, wrong, Rule::PolynomialIdentity)).is_err());
}
#[test]
fn gaussian_and_exact_fraction_fixtures() {
    let i = Expr::Constant(BigC::new(0.into(), 1.into()));
    check(&cert(
        prod(vec![i.clone(), i]),
        c(-1),
        Rule::PolynomialIdentity,
    ))
    .unwrap();
    let mut b = Budget::default();
    let half = Expr::Constant(BigC::new(BigQ::parse("1", "2", &mut b).unwrap(), 0.into()));
    let third = Expr::Constant(BigC::new(BigQ::parse("1", "3", &mut b).unwrap(), 0.into()));
    let five_sixths = Expr::Constant(BigC::new(BigQ::parse("5", "6", &mut b).unwrap(), 0.into()));
    check(&cert(
        sum(vec![half, third]),
        five_sixths,
        Rule::PolynomialIdentity,
    ))
    .unwrap();
}
#[test]
fn rational_cancellation_keeps_conditional_domain() {
    let mut q = cert(
        div(prod(vec![s("x"), s("y")]), s("x")),
        s("y"),
        Rule::RationalIdentity,
    );
    assert!(check(&q).is_err());
    q.nonzero_assumptions.push(s("x"));
    let report = check(&q).unwrap();
    assert_eq!(report.conditional_domain.len(), 1);
    assert!(!report.formal_proof);
    assert!(report
        .statement
        .contains("Independently replayed bounded algebra"));
    assert!(report.statement.contains("Not a Lean/formal QFT proof"));
    q.steps[0].rule = Rule::PolynomialIdentity;
    assert!(check(&q).is_err());
}
#[test]
fn rational_crossmultiplication_and_wrong_numerator() {
    let a = sum(vec![div(c(1), s("x")), div(c(1), s("y"))]);
    let z = div(sum(vec![s("x"), s("y")]), prod(vec![s("x"), s("y")]));
    let mut q = cert(a, z, Rule::RationalIdentity);
    q.nonzero_assumptions = vec![s("x"), s("y"), prod(vec![s("x"), s("y")])];
    check(&q).unwrap();
    q.steps[0].after = div(sum(vec![s("x"), c(2)]), prod(vec![s("x"), s("y")]));
    q.end = q.steps[0].after.clone();
    assert!(check(&q).is_err());
}
#[test]
fn nested_and_zero_masked_poles_are_not_deleted() {
    let inner = div(c(1), s("x"));
    let mut q = cert(div(c(1), inner.clone()), s("x"), Rule::RationalIdentity);
    q.nonzero_assumptions = vec![inner];
    assert!(check(&q).is_err());
    q.nonzero_assumptions.push(s("x"));
    check(&q).unwrap();
    let mut masked = cert(
        prod(vec![c(0), div(c(1), s("x"))]),
        c(0),
        Rule::RationalIdentity,
    );
    assert!(check(&masked).is_err());
    masked.nonzero_assumptions.push(s("x"));
    check(&masked).unwrap();
}
#[test]
fn normalized_zero_denominator_and_false_assumptions_rejected() {
    let zero = sum(vec![s("x"), prod(vec![c(-1), s("x")])]);
    let mut q = cert(div(c(0), zero.clone()), c(0), Rule::RationalIdentity);
    q.nonzero_assumptions.push(zero.clone());
    assert!(check(&q).is_err());
    let mut q = cert(c(1), c(1), Rule::PolynomialIdentity);
    q.nonzero_assumptions.push(zero);
    assert!(check(&q).is_err());
    q.nonzero_assumptions = vec![div(c(1), s("unused"))];
    assert!(check(&q).is_err());
}
#[test]
fn strict_chain_and_forged_end() {
    let mut q = cert(sum(vec![s("x"), c(0)]), s("x"), Rule::PolynomialIdentity);
    q.steps.push(ReplayStep {
        before: s("x"),
        after: prod(vec![c(1), s("x")]),
        rule: Rule::PolynomialIdentity,
    });
    q.end = q.steps[1].after.clone();
    check(&q).unwrap();
    q.steps[1].before = sum(vec![s("x"), c(0)]);
    assert!(check(&q).is_err());
    q.steps[1].before = s("x");
    q.end = s("x");
    assert!(check(&q).is_err());
    q.end = c(9);
    assert!(check(&q).is_err());
    q.end = q.steps[1].after.clone();
    q.start = s("x");
    assert!(check(&q).is_err());
}
#[test]
fn empty_and_convention_tampering() {
    let mut q = cert(c(1), c(1), Rule::PolynomialIdentity);
    check(&q).unwrap();
    q.version += 1;
    assert!(check(&q).is_err());
    q.version = CERTIFICATE_VERSION;
    q.domain.push_str("changed");
    assert!(check(&q).is_err());
    q.domain = SCALAR_DOMAIN.into();
    q.steps.clear();
    assert!(check(&q).is_err());
}
#[test]
fn ordered_false_commutation_and_all_unsupported_nodes() {
    use qft_symbolics::expression::{Space, TensorAtom};
    for e in [
        Expr::OrderedProduct(vec![s("a"), s("b")]),
        Expr::Operator("A".into()),
        Expr::Tensor(TensorAtom {
            name: "T".into(),
            slots: vec![],
            indices: vec![],
        }),
        Expr::Contract {
            binder: "j".into(),
            space: Space::Lorentz,
            body: Box::new(c(1)),
        },
    ] {
        assert!(check(&cert(e.clone(), e, Rule::RationalIdentity)).is_err());
    }
    let a = Expr::OrderedProduct(vec![Expr::Operator("A".into()), Expr::Operator("B".into())]);
    let z = Expr::OrderedProduct(vec![Expr::Operator("B".into()), Expr::Operator("A".into())]);
    assert!(check(&cert(a, z, Rule::PolynomialIdentity)).is_err());
    let mut q = cert(c(1), c(1), Rule::PolynomialIdentity);
    q.nonzero_assumptions.push(s(""));
    assert!(check(&q).is_err());
}
#[test]
fn shared_budget_and_aggregate_input_preflight() {
    let q = cert(c(1), c(1), Rule::PolynomialIdentity);
    let mut b = Budget::default();
    check(&q).unwrap();
    replay(&q, &mut b).unwrap();
    let used = b.used();
    let mut b = Budget::new(Limits {
        max_operations: used,
        ..Limits::default()
    });
    replay(&q, &mut b).unwrap();
    assert!(replay(&q, &mut b).is_err());
    let mut b = Budget::new(Limits {
        max_input_bytes: 200,
        ..Limits::default()
    });
    assert!(matches!(
        replay(&q, &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        })
    ));
}
#[test]
fn degree_dense_product_and_coefficient_limits() {
    let a = prod(vec![s("x"); 6]);
    let q = cert(a.clone(), a, Rule::PolynomialIdentity);
    let mut b = Budget::new(Limits {
        max_degree: 4,
        ..Limits::default()
    });
    assert!(matches!(
        replay(&q, &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::Degree,
            ..
        })
    ));
    let factors = (0..7)
        .map(|i| sum(vec![s(&format!("x{i}")), s(&format!("y{i}"))]))
        .collect();
    let a = prod(factors);
    let q = cert(a.clone(), a, Rule::PolynomialIdentity);
    let mut b = Budget::new(Limits {
        max_terms: 100,
        max_operations: 100_000_000,
        ..Limits::default()
    });
    assert!(matches!(
        replay(&q, &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::Terms,
            ..
        })
    ));
    let q = cert(c(1000), c(1000), Rule::PolynomialIdentity);
    let mut b = Budget::new(Limits {
        max_coefficient_bits: 2,
        ..Limits::default()
    });
    assert!(matches!(
        replay(&q, &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
}
#[test]
fn hard_depth_cap_even_under_huge_user_limits() {
    let mut a = c(1);
    for _ in 0..140 {
        a = prod(vec![a]);
    }
    let q = cert(a.clone(), a, Rule::PolynomialIdentity);
    let mut b = Budget::new(Limits {
        max_operations: u64::MAX,
        max_terms: usize::MAX,
        max_degree: usize::MAX,
        max_coefficient_bits: u64::MAX,
        max_input_bytes: usize::MAX,
        max_word_length: usize::MAX,
    });
    assert!(replay(&q, &mut b).is_err());
}
