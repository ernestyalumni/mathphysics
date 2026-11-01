//! Runnable explanations of independently checked algebra; no textbook runtime.
use qft_symbolics::{
    checked, clifford,
    computation::{Budget, Error, Limits},
    conventions::EpsilonConvention,
    exact::BigC,
    expression::{Expr, Index, IndexName, Space, TensorAtom, Variance},
    polynomial::BigPoly,
    replay::{replay, Certificate, ReplayStep, Rule, CERTIFICATE_VERSION, SCALAR_DOMAIN},
    scalar::C,
};
use std::collections::BTreeMap;
fn symbol(s: &str) -> Expr {
    Expr::Symbol(s.into())
}
fn product(v: Vec<Expr>) -> Expr {
    Expr::Product(v)
}
fn certificate(start: Expr, end: Expr, rule: Rule, assumptions: Vec<Expr>) -> Certificate {
    Certificate {
        version: CERTIFICATE_VERSION,
        domain: SCALAR_DOMAIN.into(),
        start: start.clone(),
        end: end.clone(),
        steps: vec![ReplayStep {
            before: start,
            after: end,
            rule,
        }],
        nonzero_assumptions: assumptions,
    }
}
fn contracted(name: &str) -> Expr {
    let atom = |tensor: &str, variance| {
        Expr::Tensor(TensorAtom {
            name: tensor.into(),
            slots: vec![(Space::Lorentz, variance)],
            indices: vec![Index {
                name: IndexName::Bound(name.into()),
                space: Space::Lorentz,
                variance,
            }],
        })
    };
    Expr::Contract {
        binder: name.into(),
        space: Space::Lorentz,
        body: Box::new(Expr::OrderedProduct(vec![
            atom("probe_vector", Variance::Up),
            atom("probe_covector", Variance::Down),
        ])),
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut b = Budget::new(Limits {
        max_operations: 1_000_000_000,
        ..Limits::default()
    });
    println!("1. TWO-GAMMA TRACE: tr(gamma^mu gamma^nu) = -4 eta^mu,nu");
    let symbolic = clifford::try_trace_exact(&["fresh_mu", "fresh_nu"], &mut b)?;
    let metric = clifford::try_metric_with::<BigC>("fresh_mu", "fresh_nu", &mut b)?;
    let key = metric.terms().keys().next().unwrap()[0].clone();
    for mu in 0..4 {
        for nu in 0..4 {
            let eta = if mu == nu {
                if mu == 0 {
                    -1
                } else {
                    1
                }
            } else {
                0
            };
            let gamma_mu = checked::gamma(mu, &mut b)?;
            let gamma_nu = checked::gamma(nu, &mut b)?;
            let actual = checked::trace(&checked::mul(&gamma_mu, &gamma_nu, &mut b)?, &mut b)?;
            let expected = C::from(-4 * eta);
            assert_eq!(actual, expected);
            // Symbolic h=-eta is evaluated separately from direct matrix multiplication.
            let values = BTreeMap::from([(key.clone(), BigC::from(-eta))]);
            assert_eq!(symbolic.try_eval(&values, &mut b)?, BigC::from(expected));
        }
    }
    let time_gamma = checked::gamma(0, &mut b)?;
    let time_trace = checked::trace(&checked::mul(&time_gamma, &time_gamma, &mut b)?, &mut b)?;
    assert_ne!(time_trace, C::from(-4));
    println!("   16 component cases agree: symbolic reduction vs explicit 4x4 matrices.");
    println!("   At mu=nu=0 the answer is +4; a candidate -4 is rejected.");
    println!(
        "   This checks evaluated components; symbolic pairing identities have separate tests.\n"
    );

    println!("2. EPSILON CONVENTIONS: raising then lowering returns the input.");
    let spinor = [
        BigPoly::try_constant(2, &mut b)?,
        BigPoly::try_constant(3, &mut b)?,
    ];
    for profile in [EpsilonConvention::Library, EpsilonConvention::MonoNotes] {
        let upper = profile.raise(&spinor, &mut b)?;
        assert_eq!(profile.lower(&upper, &mut b)?, spinor);
    }
    let upper = EpsilonConvention::Library.raise(&spinor, &mut b)?;
    assert_eq!(
        EpsilonConvention::Library.convert_upper(EpsilonConvention::MonoNotes, &upper, &mut b)?,
        EpsilonConvention::MonoNotes.raise(&spinor, &mut b)?
    );
    assert_ne!(upper, EpsilonConvention::MonoNotes.raise(&spinor, &mut b)?);
    println!("   Both profiles invert correctly; comparing upper components needs the explicit sign conversion.\n");

    println!("3. RATIONAL DOMAINS: x/x is one only where x is nonzero.");
    let x = BigPoly::try_var("domain_probe", &mut b)?;
    let rational = checked::RationalFunction::new(x.clone(), x, &mut b)?;
    assert_eq!(
        rational.try_eval(
            &BTreeMap::from([("domain_probe".into(), BigC::from(5))]),
            &mut b
        )?,
        BigC::from(1)
    );
    assert!(matches!(
        rational.try_eval(
            &BTreeMap::from([("domain_probe".into(), BigC::from(0))]),
            &mut b
        ),
        Err(Error::DivisionByZero)
    ));
    println!("   x=5 returns 1; x=0 returns a typed pole error, never a successful value.\n");

    println!("4. TYPED EXPRESSIONS: operator order and dummy-index scope are retained.");
    let op = |s: &str| Expr::Operator(s.into());
    let ab = Expr::OrderedProduct(vec![op("A"), op("B")]);
    let ba = Expr::OrderedProduct(vec![op("B"), op("A")]);
    assert_ne!(ab.canonical(&mut b)?, ba.canonical(&mut b)?);
    assert!(product(vec![op("A"), op("B")]).validate(&mut b).is_err());
    assert_eq!(
        contracted("mu").canonical(&mut b)?,
        contracted("fresh_dummy").canonical(&mut b)?
    );
    println!("   AB differs from BA; putting operators in a commuting product is rejected.");
    println!("   Renaming a properly bound dummy index leaves canonical structure unchanged.\n");

    println!("5. REPLAY: (a+b)c = ac+bc, checked without the production polynomial normalizer.");
    let (a, c, d) = (
        symbol("heldout_a"),
        symbol("heldout_b"),
        symbol("heldout_c"),
    );
    let left = product(vec![Expr::Sum(vec![a.clone(), c.clone()]), d.clone()]);
    let right = Expr::Sum(vec![product(vec![a, d.clone()]), product(vec![c, d])]);
    let cert = certificate(left, right, Rule::PolynomialIdentity, vec![]);
    let report = replay(&cert, &mut b)?;
    assert_eq!(report.steps_checked, 1);
    assert!(!report.formal_proof);
    let mut forged = cert.clone();
    forged.steps[0].after = symbol("unrelated_result");
    assert!(replay(&forged, &mut b).is_err());
    println!(
        "   The independent sparse-ring checker accepts the step and rejects a tampered result."
    );
    let x = symbol("pole_probe");
    let quotient = Expr::Quotient(Box::new(x.clone()), Box::new(x.clone()));
    let cert = certificate(
        quotient,
        Expr::Constant(BigC::from(1)),
        Rule::RationalIdentity,
        vec![x],
    );
    let report = replay(&cert, &mut b)?;
    assert!(!report.conditional_domain.is_empty());
    assert!(!report.formal_proof);
    let mut missing = cert;
    missing.nonzero_assumptions.clear();
    assert!(replay(&missing, &mut b).is_err());
    println!(
        "   x/x=1 replays only with its explicit nonzero assumption; omitting it is rejected."
    );
    println!("\nAll demonstrated checks met their positive and negative expectations.");
    println!("Evidence: exact computation and bounded independent algebra replay, not a formal QFT proof.");
    Ok(())
}
