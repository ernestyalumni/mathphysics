use qft_symbolics::{
    coefficient::Coefficient,
    computation::{Budget, Error, Limits, Resource},
    exact::{BigC, BigQ},
    scalar::{C, Q},
};
fn budget() -> Budget {
    Budget::new(Limits {
        max_operations: 1_000_000_000,
        ..Limits::default()
    })
}
#[test]
fn independent_decimal_fixtures() {
    let mut b = budget();
    let x = BigQ::parse(
        "115792089237316195423570985008687907853269984665640564039457584007913129639936",
        "1",
        &mut b,
    )
    .unwrap();
    assert_eq!(BigQ::from(2).checked_pow(256, &mut b).unwrap(), x);
    assert_eq!(
        x.checked_add(&BigQ::from(1), &mut b).unwrap().to_string(),
        "115792089237316195423570985008687907853269984665640564039457584007913129639937"
    );
    let y = BigQ::parse(
        "1",
        "115792089237316195423570985008687907853269984665640564039457584007913129639936",
        &mut b,
    )
    .unwrap();
    assert_eq!(x.checked_mul(&y, &mut b).unwrap(), BigQ::from(1));
    let z = BigC::new(x.clone(), x.clone());
    assert_eq!(z.checked_div(&z, &mut b).unwrap(), BigC::from(1));
    let product = z.checked_mul(&BigC::new(y.clone(), y), &mut b).unwrap();
    assert_eq!(product, BigC::new(0.into(), 2.into()));
    assert!(x.try_to_small().is_err());
}
#[test]
fn rational_gaussian_fixture() {
    let mut b = budget();
    // Independently: norm = 221/1225, real numerator = 16/105, imaginary = 43/210.
    let z = BigC::new(
        BigQ::parse("1", "2", &mut b).unwrap(),
        BigQ::parse("1", "3", &mut b).unwrap(),
    );
    let w = BigC::new(
        BigQ::parse("2", "5", &mut b).unwrap(),
        BigQ::parse("-1", "7", &mut b).unwrap(),
    );
    let r = z.checked_div(&w, &mut b).unwrap();
    assert_eq!(r.re.to_string(), "560/663");
    assert_eq!(r.im.to_string(), "1505/1326");
    assert_eq!(r.checked_mul(&w, &mut b).unwrap(), z);
    assert_eq!(z.conjugate().conjugate(), z);
    assert!(BigC::default().is_real());
    assert_eq!(BigQ::default(), 0.into());
}
#[test]
fn errors_and_limits() {
    let mut b = budget();
    assert_eq!(BigQ::parse("1", "0", &mut b), Err(Error::ZeroDenominator));
    assert_eq!(BigQ::parse("1x", "2", &mut b), Err(Error::InvalidInteger));
    assert_eq!(BigQ::parse("", "2", &mut b), Err(Error::InvalidInteger));
    assert_eq!(
        BigQ::from(1).checked_div(&0.into(), &mut b),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        BigC::from(1).checked_div(&BigC::default(), &mut b),
        Err(Error::DivisionByZero)
    );
    let mut b = Budget::new(Limits {
        max_input_bytes: 2,
        ..Limits::default()
    });
    assert!(matches!(
        BigQ::parse("123", "1", &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        })
    ));
    let mut b = Budget::new(Limits {
        max_coefficient_bits: 4,
        ..Limits::default()
    });
    assert!(matches!(
        BigQ::parse("99", "1", &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    assert!(matches!(
        BigQ::from(15).checked_mul(&15.into(), &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    let mut b = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(matches!(
        BigC::from(1).checked_neg(&mut b),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        })
    ));
}
#[test]
fn coefficient_and_small_crosschecks() {
    fn generic<T: Coefficient>(a: T, b: T, budget: &mut Budget) -> T {
        a.checked_mul(&b, budget).unwrap()
    }
    let mut b = budget();
    for n in -4..=4 {
        for m in -4..=4 {
            let a = C::new(Q::new(n, 3), Q::new(m, 5));
            let c = C::new(Q::new(m, 7), Q::new(n, 2));
            let x = BigC::from(a);
            let y = BigC::from(c);
            assert_eq!(generic(x.clone(), y.clone(), &mut b), BigC::from(a * c));
            assert_eq!(x.checked_add(&y, &mut b).unwrap(), BigC::from(a + c));
            assert_eq!(x.checked_sub(&y, &mut b).unwrap(), BigC::from(a - c));
            if !c.is_zero() {
                assert_eq!(x.checked_div(&y, &mut b).unwrap(), BigC::from(a / c));
            }
            assert_eq!(x.checked_neg(&mut b).unwrap(), BigC::from(-a));
            assert_eq!(x.checked_pow(3, &mut b).unwrap(), BigC::from(a * a * a));
        }
    }
}

#[test]
fn normalization_and_intermediate_limits() {
    let mut b = budget();
    let q = BigQ::parse("+0006", "-0008", &mut b).unwrap();
    assert_eq!(q.numerator().to_string(), "-3");
    assert_eq!(q.denominator().to_string(), "4");
    assert_eq!(
        BigQ::parse("0", "-9", &mut b)
            .unwrap()
            .denominator()
            .to_string(),
        "1"
    );
    assert_eq!(BigQ::from(0).checked_pow(0, &mut b).unwrap(), 1.into());
    let x = BigQ::parse(
        "100000000000000000000000000000000000000000000000000000000000000000000000000000000",
        "3",
        &mut b,
    )
    .unwrap();
    let y = BigQ::parse(
        "3",
        "100000000000000000000000000000000000000000000000000000000000000000000000000000000",
        &mut b,
    )
    .unwrap();
    assert_eq!(x.checked_mul(&y, &mut b).unwrap(), 1.into());
    let mut tiny = Budget::new(Limits {
        max_coefficient_bits: 4,
        ..Limits::default()
    });
    // Final answer is zero, but cross-products exceed the declared bound.
    assert!(matches!(
        BigQ::from(15).checked_sub(&15.into(), &mut tiny),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    assert!(matches!(
        x.checked_neg(&mut tiny),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    assert_eq!(tiny.used(), 0);
}
