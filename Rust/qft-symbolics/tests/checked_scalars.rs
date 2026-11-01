use qft_symbolics::{
    coefficient::Coefficient,
    computation::{Budget, Error, Limits, Resource},
    scalar::{C, Q},
};
#[test]
fn checked_boundaries() {
    assert_eq!(Q::try_new(1, 0), Err(Error::ZeroDenominator));
    assert!(matches!(
        Q::try_new(i128::MIN, 1),
        Err(Error::Overflow { .. })
    ));
    assert_eq!(Q::try_new(i128::MIN, i128::MIN).unwrap(), Q::ONE);
    assert_eq!(Q::try_new(i128::MIN, -2).unwrap(), Q::new(1i128 << 126, 1));
    assert_eq!(Q::try_new(0, i128::MIN).unwrap(), Q::ZERO);
    assert!(Q::try_new(1, i128::MIN).is_err());
    assert_eq!(
        Q::try_new(i128::MIN, 2).unwrap(),
        Q::new(-(1i128 << 126), 1)
    );
    assert_eq!(Q::try_new(2, i128::MIN).unwrap(), Q::new(-1, 1i128 << 126));
    assert_eq!(Q::try_new(-2, i128::MIN).unwrap(), Q::new(1, 1i128 << 126));
    assert!(matches!(
        Q::try_new(i128::MIN, -1),
        Err(Error::Overflow { .. })
    ));
    assert_eq!(Q::try_new(i128::MIN, 0), Err(Error::ZeroDenominator));

    let mut b = Budget::default();
    assert!(matches!(
        Q::new(i128::MAX, 1).checked_add(&Q::ONE, &mut b),
        Err(Error::Overflow { .. })
    ));
    assert_eq!(
        Q::ONE.checked_div(&Q::ZERO, &mut b),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        C::ONE.checked_div(&C::ZERO, &mut b),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        Q::new(i128::MAX, 2)
            .checked_mul(&Q::new(2, i128::MAX), &mut b)
            .unwrap(),
        Q::ONE
    );
    let mut b = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(matches!(
        C::ONE.checked_add(&C::ONE, &mut b),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        })
    ));
}
#[test]
fn unchanged_legacy_and_checked_agree() {
    fn copy<T: Copy>(v: T) -> T {
        v
    }
    let mut b = Budget::new(Limits {
        max_operations: 100_000_000,
        ..Limits::default()
    });
    assert_eq!(copy(Q::ONE), Q::default() + Q::ONE);
    assert_eq!(copy(C::I), C::default() + C::I);
    for n in -3..=3 {
        for m in -3..=3 {
            let x = Q::new(n, 5);
            let y = Q::new(m, 7);
            assert_eq!(x.checked_add(&y, &mut b).unwrap(), x + y);
            assert_eq!(x.checked_sub(&y, &mut b).unwrap(), x - y);
            assert_eq!(x.checked_mul(&y, &mut b).unwrap(), x * y);
            if !y.is_zero() {
                assert_eq!(x.checked_div(&y, &mut b).unwrap(), x / y);
            }
            assert_eq!(x.checked_neg(&mut b).unwrap(), -x);
            assert_eq!(x.checked_pow(4, &mut b).unwrap(), x.pow(4));
            let z = C::new(x, y);
            let w = C::new(y, x);
            assert_eq!(z.checked_add(&w, &mut b).unwrap(), z + w);
            assert_eq!(z.checked_sub(&w, &mut b).unwrap(), z - w);
            assert_eq!(z.checked_mul(&w, &mut b).unwrap(), z * w);
            if !w.is_zero() {
                assert_eq!(z.checked_div(&w, &mut b).unwrap(), z / w);
            }
            assert_eq!(z.checked_neg(&mut b).unwrap(), -z);
            assert_eq!(z.checked_pow(3, &mut b).unwrap(), z * z * z);
            assert_eq!(Coefficient::conjugate(&z), z.conj());
        }
    }
}
