use lie_groups::su2_observation::{observe, Chart, Slice};
use qft_symbolics::{
    clifford,
    computation::{Budget, Error, Limits},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
    scalar::Q,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut b = Budget::new(Limits {
        max_operations: 100_000_000,
        ..Limits::default()
    });
    let huge = BigQ::parse(
        "115792089237316195423570985008687907853269984665640564039457584007913129639936",
        "1",
        &mut b,
    )?;
    assert_eq!(BigQ::from(2).checked_pow(256, &mut b)?, huge);
    assert_eq!(
        huge.checked_add(&BigQ::from(1), &mut b)?.to_string(),
        "115792089237316195423570985008687907853269984665640564039457584007913129639937"
    );
    let trace = clifford::try_trace_exact(&["independent_probe_u", "independent_probe_v"], &mut b)?;
    let scaled = trace.try_mul(
        &BigPoly::try_constant(BigC::new(huge, 0.into()), &mut b)?,
        &mut b,
    )?;
    assert_eq!(scaled.term_count(), 1);
    assert_eq!(
        scaled.terms().values().next().unwrap().re.to_string(),
        "463168356949264781694283940034751631413079938662562256157830336031652518559744"
    );
    assert!(matches!(
        Q::new(i128::MAX, 1).checked_add(&Q::ONE, &mut b),
        Err(Error::Overflow { .. })
    ));
    let mut zero = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(clifford::try_trace(&["u", "v"], &mut zero).is_err());
    let o = observe([2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0], std::f64::consts::PI, 0.0)?;
    // Direct Pauli expansion at pi: U=-i n.sigma; no HTTP/JSON or constructor round trip.
    let expected = [
        [[0.0, -1.0 / 3.0], [-2.0 / 3.0, -2.0 / 3.0]],
        [[2.0 / 3.0, -2.0 / 3.0], [0.0, 1.0 / 3.0]],
    ];
    for (row, want) in o.matrix.iter().zip(expected) {
        for (entry, want) in row.iter().zip(want) {
            for (v, want) in entry.iter().zip(want) {
                assert!((v - want).abs() < 2e-12);
            }
        }
    }
    assert!(matches!(o.stereographic, Chart::Finite(_)));
    assert!(matches!(o.slice, Slice::Sphere { .. }));
    assert!(observe([0.0, 0.0, 0.0], 0.0, 0.0).is_err());
    // New capabilities are consumed through public library interfaces.
    use qft_symbolics::{
        checked,
        conventions::EpsilonConvention,
        expression::Expr,
        replay::{self, Certificate, ReplayStep, Rule, CERTIFICATE_VERSION, SCALAR_DOMAIN},
    };
    let exact = BigC::new(
        BigQ::parse("340282366920938463463374607431768211461", "1", &mut b)?,
        0.into(),
    );
    let matrix = [
        [exact.clone(), BigC::from(0)],
        [BigC::from(0), BigC::from(7)],
    ];
    assert_eq!(
        checked::trace(&matrix, &mut b)?.re.to_string(),
        "340282366920938463463374607431768211468"
    );
    let spinor = [
        BigPoly::try_constant(7, &mut b)?,
        BigPoly::try_constant(-11, &mut b)?,
    ];
    let raised = EpsilonConvention::MonoNotes.raise(&spinor, &mut b)?;
    assert_eq!(raised[0], BigPoly::try_constant(11, &mut b)?);
    assert_eq!(raised[1], BigPoly::try_constant(7, &mut b)?);
    let x = Expr::Symbol("external_heldout".into());
    let start = Expr::Sum(vec![x.clone(), x.clone()]);
    let end = Expr::Product(vec![Expr::Constant(BigC::from(2)), x.clone()]);
    assert_ne!(start.canonical(&mut b)?, end.canonical(&mut b)?);
    let cert = Certificate {
        version: CERTIFICATE_VERSION,
        domain: SCALAR_DOMAIN.into(),
        start: start.clone(),
        end: end.clone(),
        steps: vec![ReplayStep {
            before: start,
            after: end,
            rule: Rule::PolynomialIdentity,
        }],
        nonzero_assumptions: vec![],
    };
    let report = replay::replay(&cert, &mut b)?;
    assert_eq!(report.steps_checked, 1);
    assert!(!report.formal_proof);
    let ratio = checked::RationalFunction::new(
        BigPoly::try_constant(3, &mut b)?,
        BigPoly::try_var("external_denominator", &mut b)?,
        &mut b,
    )?;
    assert!(ratio
        .try_eval(
            &std::collections::BTreeMap::from([("external_denominator".into(), BigC::from(0))]),
            &mut b
        )
        .is_err());
    use qft_symbolics::{
        differential::DifferentialOperator as D, ordered::OrderedPolynomial as O, relativity as rel,
    };
    let u = BigQ::parse("-1", "5", &mut b)?;
    let boost = rel::exact_half_velocity_boost(&u, rel::Axis::X, &mut b)?;
    let vector = [BigC::from(3), BigC::from(4), BigC::from(0), BigC::from(0)];
    let image = rel::apply(&boost, &vector, &mut b)?;
    assert_eq!(image[0].re, BigQ::parse("19", "12", &mut b)?);
    assert_eq!(image[1].re, BigQ::parse("37", "12", &mut b)?);
    assert_eq!(rel::dot(&image, &image, &mut b)?, BigC::from(7));
    let a = O::generator("heldout_operator_1", &mut b)?;
    let c = O::generator("heldout_operator_2", &mut b)?;
    assert_ne!(a.try_mul(&c, &mut b)?, c.try_mul(&a, &mut b)?);
    let coords = ["heldout_coordinate".into()];
    let derivative = D::partial(coords.clone(), 0, &mut b)?;
    let coordinate = D::multiplication(
        coords.clone(),
        BigPoly::try_var("heldout_coordinate", &mut b)?,
        &mut b,
    )?;
    let identity = D::multiplication(coords, BigPoly::try_constant(1, &mut b)?, &mut b)?;
    assert_eq!(derivative.commutator(&coordinate, &mut b)?, identity);
    println!("research_consumer_ok: arbitrary coefficients, shared trace, typed errors, independent SU2 matrix, rapidity, ordered and differential operators");
    Ok(())
}
