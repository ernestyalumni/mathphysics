use qft_symbolics::{
    checked::*,
    coefficient::Coefficient,
    computation::{Budget, Error, Limits, Resource},
    exact::{BigC, BigQ},
    polynomial::Polynomial,
    scalar::C,
    tensor,
};
use std::collections::BTreeMap;
fn budget() -> Budget {
    Budget::new(Limits {
        max_operations: 100_000_000,
        ..Limits::default()
    })
}
fn p(n: i64, b: &mut Budget) -> Polynomial<BigC> {
    Polynomial::try_constant(BigC::from_i64(n), b).unwrap()
}
#[test]
fn independent_chiral_matrices_and_clifford_form() {
    let mut b = budget();
    // Independent literal chiral representation, eta=(-,+,+,+).
    let z = C::ZERO;
    let o = C::ONE;
    let n = C::from_i64(-1);
    let i = C::I;
    let ni = i.checked_neg(&mut b).unwrap();
    let expected = [
        [[z, z, o, z], [z, z, z, o], [o, z, z, z], [z, o, z, z]],
        [[z, z, z, o], [z, z, o, z], [z, n, z, z], [n, z, z, z]],
        [[z, z, z, ni], [z, z, i, z], [z, i, z, z], [ni, z, z, z]],
        [[z, z, o, z], [z, z, z, n], [n, z, z, z], [z, o, z, z]],
    ];
    for mu in 0..4 {
        assert_eq!(gamma(mu, &mut b).unwrap(), expected[mu]);
        for nu in 0..4 {
            let a = mul(&expected[mu], &expected[nu], &mut b).unwrap();
            let c = mul(&expected[nu], &expected[mu], &mut b).unwrap();
            let factor = if mu != nu {
                0
            } else if mu == 0 {
                2
            } else {
                -2
            };
            assert_eq!(
                add(&a, &c, &mut b).unwrap(),
                scale(
                    &identity::<C, 4>(&mut b).unwrap(),
                    &C::from_i64(factor),
                    &mut b
                )
                .unwrap()
            );
        }
    }
    assert!(gamma(usize::MAX, &mut b).is_err());
    assert!(sigma(4, &mut b).is_err());
}
#[test]
fn large_coefficients_matrix_and_rational() {
    let mut b = budget();
    let h = BigC::new(
        BigQ::parse("340282366920938463463374607431768211456", "1", &mut b).unwrap(),
        BigQ::default(),
    );
    let m = [
        [h.clone(), BigC::from_i64(2)],
        [BigC::from_i64(3), BigC::from_i64(4)],
    ];
    let square = mul(&m, &m, &mut b).unwrap();
    let expected = BigQ::parse(
        "115792089237316195423570985008687907853269984665640564039457584007913129639942",
        "1",
        &mut b,
    )
    .unwrap();
    assert_eq!(square[0][0], BigC::new(expected, BigQ::default()));
    assert_eq!(
        trace(&transpose(&m, &mut b).unwrap(), &mut b).unwrap(),
        h.checked_add(&BigC::from_i64(4), &mut b).unwrap()
    );
    let r = RationalFunction::new(
        Polynomial::try_constant(h.clone(), &mut b).unwrap(),
        p(2, &mut b),
        &mut b,
    )
    .unwrap();
    assert_eq!(
        r.try_eval(&BTreeMap::new(), &mut b).unwrap(),
        h.checked_div(&BigC::from_i64(2), &mut b).unwrap()
    );
}
#[test]
fn poles_and_domains_survive_all_operations() {
    let mut b = budget();
    let x = Polynomial::<BigC>::try_var("u", &mut b).unwrap();
    let r = RationalFunction::new(p(1, &mut b), x.clone(), &mut b).unwrap();
    let one = RationalFunction::new(p(1, &mut b), p(1, &mut b), &mut b).unwrap();
    let zero = RationalFunction::new(p(0, &mut b), p(1, &mut b), &mut b).unwrap();
    assert_eq!(one.try_div(&zero, &mut b), Err(Error::DivisionByZero));
    assert_eq!(
        RationalFunction::new(p(1, &mut b), p(0, &mut b), &mut b),
        Err(Error::ZeroDenominator)
    );
    let at_zero = BTreeMap::from([("u".into(), BigC::from_i64(0))]);
    for v in [
        r.try_div(&r, &mut b).unwrap(),
        one.try_div(&r, &mut b).unwrap(),
        r.try_sub(&r, &mut b).unwrap(),
        r.try_mul(&zero, &mut b).unwrap(),
        r.try_neg(&mut b).unwrap(),
        r.try_add(&one, &mut b).unwrap(),
    ] {
        assert_eq!(v.try_eval(&at_zero, &mut b), Err(Error::DivisionByZero));
    }
    assert!(r
        .try_div(&r, &mut b)
        .unwrap()
        .try_residual(&one, &mut b)
        .unwrap()
        .is_zero());
    assert!(matches!(
        r.try_eval(&BTreeMap::new(), &mut b),
        Err(Error::MissingVariable { .. })
    ));
    let at_two = BTreeMap::from([("u".into(), BigC::from_i64(2))]);
    assert_eq!(
        one.try_div(&r, &mut b)
            .unwrap()
            .try_eval(&at_two, &mut b)
            .unwrap(),
        BigC::from_i64(2)
    );
}
#[test]
fn polynomial_schouten_and_independent_numeric_determinant() {
    let mut b = budget();
    let mut named = |s: &str| {
        [
            Polynomial::<BigC>::try_var(format!("{s}0"), &mut b).unwrap(),
            Polynomial::try_var(format!("{s}1"), &mut b).unwrap(),
        ]
    };
    let a = named("fresh");
    let c = named("other");
    let d = named("third");
    let e = named("fourth");
    assert!(schouten(&a, &c, &d, &e, &mut b).unwrap().is_zero());
    let vals = BTreeMap::from([
        ("fresh0".into(), BigC::from_i64(7)),
        ("fresh1".into(), BigC::from_i64(-2)),
        ("other0".into(), BigC::from_i64(3)),
        ("other1".into(), BigC::from_i64(5)),
    ]);
    assert_eq!(
        bracket(&a, &c, &mut b)
            .unwrap()
            .try_eval(&vals, &mut b)
            .unwrap(),
        BigC::from_i64(41)
    );
}
#[test]
fn bcfw_preserves_total_momentum_and_is_atomic() {
    let mut b = budget();
    let mut l = [[1, 2], [3, 5]].map(|v| v.map(BigC::from_i64));
    let mut t = [[7, 11], [13, 17]].map(|v| v.map(BigC::from_i64));
    let before = (l.clone(), t.clone());
    bcfw_shift(&mut l, &mut t, 0, 1, &BigC::from_i64(2), &mut b).unwrap();
    // Independently integer-computed total momentum components.
    for (a, row) in [[46, 62], [79, 107]].iter().enumerate() {
        for (c, n) in row.iter().enumerate() {
            let s = l[0][a]
                .checked_mul(&t[0][c], &mut b)
                .unwrap()
                .checked_add(&l[1][a].checked_mul(&t[1][c], &mut b).unwrap(), &mut b)
                .unwrap();
            assert_eq!(s, BigC::from_i64(*n));
        }
    }
    let snapshot = (l.clone(), t.clone());
    let mut tiny = Budget::new(Limits {
        max_operations: 10,
        ..Limits::default()
    });
    assert!(bcfw_shift(&mut l, &mut t, 0, 1, &BigC::from_i64(2), &mut tiny).is_err());
    assert_eq!((l.clone(), t.clone()), snapshot);
    assert!(bcfw_shift(&mut l, &mut t, 0, 0, &BigC::from_i64(2), &mut b).is_err());
    assert_ne!((l, t), before);
}
#[test]
fn budgets_and_typed_tensor_errors() {
    let mut b = budget();
    let m = identity::<BigC, 2>(&mut b).unwrap();
    let mut tiny = Budget::new(Limits {
        max_terms: 3,
        ..Limits::default()
    });
    assert!(matches!(
        mul(&m, &m, &mut tiny),
        Err(Error::ResourceLimit {
            resource: Resource::Terms,
            ..
        })
    ));
    let mut tiny = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(identity::<BigC, 2>(&mut tiny).is_err());
    let mut validator = TensorValidator::new(100);
    assert!(matches!(
        validator.require_symmetric_metric(tensor::Space::Undotted, &mut b),
        Err(TensorError::Validation(tensor::Error::Unsupported(_)))
    ));
    validator
        .declare(
            "v",
            vec![tensor::Slot::Lorentz],
            tensor::Symmetry::None,
            &mut b,
        )
        .unwrap();
    assert!(matches!(
        validator.declare("v", vec![], tensor::Symmetry::None, &mut b),
        Err(TensorError::Validation(tensor::Error::Invalid(_)))
    ));
    assert!(matches!(
        validator.declare(
            "bad",
            vec![tensor::Slot::Lorentz; 2],
            tensor::Symmetry::Generators(vec![(vec![0, 0], 1)]),
            &mut b
        ),
        Err(TensorError::Validation(tensor::Error::Invalid(_)))
    ));
    assert!(matches!(
        TensorValidator::new(0).declare("v", vec![], tensor::Symmetry::None, &mut b),
        Err(TensorError::Validation(tensor::Error::Budget))
    ));
}

#[test]
fn coefficient_degree_and_input_limits_propagate() {
    let mut b = budget();
    let h = BigC::new(
        BigQ::parse("340282366920938463463374607431768211456", "1", &mut b).unwrap(),
        BigQ::default(),
    );
    let mut bits = Budget::new(Limits {
        max_coefficient_bits: 128,
        ..Limits::default()
    });
    assert!(matches!(
        transpose(&[[h]], &mut bits),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    let x = Polynomial::<BigC>::try_var("long_variable", &mut b).unwrap();
    let r = RationalFunction::new(p(1, &mut b), x, &mut b).unwrap();
    let mut bytes = Budget::new(Limits {
        max_input_bytes: 2,
        ..Limits::default()
    });
    assert!(matches!(
        r.try_neg(&mut bytes),
        Err(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        })
    ));
    let mut degree = Budget::new(Limits {
        max_degree: 1,
        ..Limits::default()
    });
    assert!(matches!(
        r.try_mul(&r, &mut degree),
        Err(Error::ResourceLimit {
            resource: Resource::Degree,
            ..
        })
    ));
    let mut work = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(matches!(
        r.try_residual(&r, &mut work),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        })
    ));
}

#[test]
fn tensor_factorial_closure_preflight_and_signed_group_bound() {
    let mut validator = TensorValidator::new(400_000);
    let snapshot = validator.context().clone();
    let mut tiny = Budget::new(Limits {
        max_operations: 9,
        max_terms: 1,
        max_input_bytes: 1,
        ..Limits::default()
    });
    assert!(matches!(
        validator.declare(
            "T",
            vec![tensor::Slot::Lorentz; 8],
            tensor::Symmetry::Symmetric,
            &mut tiny
        ),
        Err(TensorError::Computation(Error::ResourceLimit { .. }))
    ));
    assert_eq!(validator.context(), &snapshot);
    assert_eq!(validator.retained_bounds(), (0, 0, 0));
    // Isolate closure-work exhaustion with permissive term/storage limits.
    let mut work = Budget::new(Limits {
        max_operations: 9,
        max_terms: 1_000_000,
        max_input_bytes: 1_000_000_000,
        ..Limits::default()
    });
    assert!(matches!(
        validator.declare(
            "T",
            vec![tensor::Slot::Lorentz; 8],
            tensor::Symmetry::Symmetric,
            &mut work
        ),
        Err(TensorError::Computation(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        }))
    ));
    assert_eq!(validator.context(), &snapshot);
    let mut b = budget();
    validator
        .declare(
            "signed",
            vec![tensor::Slot::Lorentz],
            tensor::Symmetry::Generators(vec![(vec![0], -1)]),
            &mut b,
        )
        .unwrap();
    assert_eq!(validator.retained_bounds().0, 1);
    assert_eq!(validator.retained_bounds().1, 2);
    let mut builtin = TensorValidator::new(100);
    let mut b = budget();
    builtin
        .declare(
            "symmetric",
            vec![tensor::Slot::Lorentz; 2],
            tensor::Symmetry::Symmetric,
            &mut b,
        )
        .unwrap();
    assert_eq!(builtin.retained_bounds().1, 2);
}

#[test]
fn tensor_registry_aggregate_limits_survive_renewed_budgets() {
    let mut validator = TensorValidator::new(1);
    let limits = Limits {
        max_terms: 7,
        ..Limits::default()
    };
    validator
        .declare(
            "A",
            vec![],
            tensor::Symmetry::None,
            &mut Budget::new(limits),
        )
        .unwrap();
    validator
        .declare(
            "B",
            vec![],
            tensor::Symmetry::None,
            &mut Budget::new(limits),
        )
        .unwrap();
    let snapshot = validator.context().clone();
    let counters = validator.retained_bounds();
    assert!(matches!(
        validator.declare(
            "C",
            vec![],
            tensor::Symmetry::None,
            &mut Budget::new(limits)
        ),
        Err(TensorError::Computation(Error::ResourceLimit {
            resource: Resource::Terms,
            ..
        }))
    ));
    assert_eq!(validator.context(), &snapshot);
    assert_eq!(validator.retained_bounds(), counters);
    let mut bytes_validator = TensorValidator::new(1);
    let bytes = Limits {
        max_input_bytes: 4 * 769,
        ..Limits::default()
    };
    bytes_validator
        .declare("A", vec![], tensor::Symmetry::None, &mut Budget::new(bytes))
        .unwrap();
    let snapshot = bytes_validator.context().clone();
    let counters = bytes_validator.retained_bounds();
    assert!(matches!(
        bytes_validator.declare("B", vec![], tensor::Symmetry::None, &mut Budget::new(bytes)),
        Err(TensorError::Computation(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        }))
    ));
    assert_eq!(bytes_validator.context(), &snapshot);
    assert_eq!(bytes_validator.retained_bounds(), counters);
    // Native validation and native orbit errors also leave counters unchanged.
    let registry_snapshot = validator.context().clone();
    let registry_counters = validator.retained_bounds();
    let mut b = budget();
    assert!(matches!(
        validator.declare("A", vec![], tensor::Symmetry::None, &mut b),
        Err(TensorError::Validation(tensor::Error::Invalid(_)))
    ));
    assert_eq!(validator.context(), &registry_snapshot);
    assert_eq!(validator.retained_bounds(), registry_counters);
    let mut exhausted = TensorValidator::new(0);
    let empty = exhausted.context().clone();
    assert!(matches!(
        exhausted.declare("X", vec![], tensor::Symmetry::None, &mut b),
        Err(TensorError::Validation(tensor::Error::Budget))
    ));
    assert_eq!(exhausted.context(), &empty);
    assert_eq!(exhausted.retained_bounds(), (0, 0, 0));
}
