use qft_symbolics::computation::{Budget, Error, Limits, Resource};
#[test]
fn exact_charge_and_failed_charge_are_deterministic() {
    let mut b = Budget::new(Limits {
        max_operations: 3,
        ..Limits::default()
    });
    assert_eq!(b.charge(3), Ok(()));
    assert!(matches!(
        b.charge(1),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            limit: 3,
            required: 4
        })
    ));
    assert_eq!(b.used(), 3);
}
#[test]
fn counter_overflow_never_wraps_or_succeeds() {
    let mut b = Budget::new(Limits {
        max_operations: u64::MAX,
        ..Limits::default()
    });
    b.charge(u64::MAX).unwrap();
    assert!(b.charge(1).is_err());
    assert_eq!(b.used(), u64::MAX);
}
#[test]
fn resources_and_zero_limits_are_explicit() {
    let b = Budget::new(Limits {
        max_terms: 0,
        max_degree: 0,
        max_coefficient_bits: 0,
        max_input_bytes: 0,
        max_word_length: 0,
        ..Limits::default()
    });
    for r in [
        Resource::Terms,
        Resource::Degree,
        Resource::CoefficientBits,
        Resource::InputBytes,
        Resource::WordLength,
    ] {
        assert!(b.check(r, 0).is_ok());
        assert!(matches!(
            b.check(r, 1),
            Err(Error::ResourceLimit {
                limit: 0,
                required: 1,
                ..
            })
        ));
    }
}
