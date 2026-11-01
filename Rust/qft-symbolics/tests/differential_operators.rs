use qft_symbolics::{
    computation::{Budget, Error, Limits, Resource},
    differential::{differentiate, DifferentialOperator as Op},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
};
use std::collections::BTreeMap;
fn coords<const N: usize>() -> [String; N] {
    std::array::from_fn(|i| format!("x{i}"))
}
fn constant(n: i64, b: &mut Budget) -> BigPoly {
    BigPoly::try_constant(n, b).unwrap()
}
fn monomial(powers: &[usize], b: &mut Budget) -> BigPoly {
    let mut p = constant(1, b);
    for (i, n) in powers.iter().enumerate() {
        let v = BigPoly::try_var(format!("x{i}"), b).unwrap();
        p = p.try_mul(&v.try_pow(*n, b).unwrap(), b).unwrap();
    }
    p
}
fn imag(n: i64, b: &mut Budget) -> BigPoly {
    BigPoly::try_constant(BigC::new(BigQ::from(0), BigQ::from(n)), b).unwrap()
}
#[test]
fn derivative_product_rule_and_mixed_binomials() {
    let b = &mut Budget::default();
    let f = monomial(&[3, 2], b);
    let g = monomial(&[2, 1], b).try_add(&constant(7, b), b).unwrap();
    let lhs = differentiate(&f.try_mul(&g, b).unwrap(), "x0", b).unwrap();
    let rhs = differentiate(&f, "x0", b)
        .unwrap()
        .try_mul(&g, b)
        .unwrap()
        .try_add(
            &f.try_mul(&differentiate(&g, "x0", b).unwrap(), b).unwrap(),
            b,
        )
        .unwrap();
    assert_eq!(lhs, rhs);
    let expected = monomial(&[2, 1], b).try_mul(&constant(6, b), b).unwrap();
    let mixed = differentiate(&differentiate(&f, "x0", b).unwrap(), "x1", b).unwrap();
    assert_eq!(mixed, expected);
    assert_eq!(
        mixed,
        differentiate(&differentiate(&f, "x1", b).unwrap(), "x0", b).unwrap()
    );
    assert!(differentiate(&f, "parameter", b).unwrap().is_zero());
    // D_x^2 D_y composed with multiplication x^2 y: 6 nonzero Leibniz terms.
    let d = Op::<2>::from_terms(coords(), BTreeMap::from([([2, 1], constant(1, b))]), b).unwrap();
    let m = Op::<2>::multiplication(coords(), monomial(&[2, 1], b), b).unwrap();
    let product = d.compose(&m, b).unwrap();
    assert_eq!(product.terms().len(), 6);
    assert_eq!(
        product.terms()[&[1, 0]],
        monomial(&[1, 0], b).try_mul(&constant(4, b), b).unwrap()
    );
    assert_eq!(product.terms()[&[0, 0]], constant(2, b));
    let witness = f.try_add(&g, b).unwrap();
    assert_eq!(
        product.apply(&witness, b).unwrap(),
        d.apply(&m.apply(&witness, b).unwrap(), b).unwrap()
    );
}
#[test]
fn exact_large_coefficients_and_second_order_cancellation() {
    let b = &mut Budget::default();
    let huge = BigQ::parse("123456789012345678901234567890", "7", b).unwrap();
    let c = BigPoly::try_constant(BigC::new(huge.clone(), BigQ::from(0)), b).unwrap();
    let x = monomial(&[1], b);
    assert_eq!(
        differentiate(&x.try_mul(&c, b).unwrap(), "x0", b).unwrap(),
        c
    );
    let d = Op::<1>::partial(coords(), 0, b).unwrap();
    let xd = d.try_scale(&x, b).unwrap();
    assert_eq!(d.commutator(&xd, b).unwrap(), d);
    let dd = d.compose(&d, b).unwrap();
    let m = Op::<1>::multiplication(coords(), x, b).unwrap();
    let two = constant(2, b);
    assert_eq!(dd.commutator(&m, b).unwrap(), d.try_scale(&two, b).unwrap());
    assert!(xd
        .commutator(&d, b)
        .unwrap()
        .terms()
        .keys()
        .all(|a| a[0] <= 1));
}
fn eta(m: usize, n: usize) -> i64 {
    if m != n {
        0
    } else if m == 0 {
        -1
    } else {
        1
    }
}
fn l(m: usize, n: usize, b: &mut Budget) -> Op<4> {
    Op::orbital_lorentz(coords(), m, n, BigQ::from(1), b).unwrap()
}
fn p(m: usize, b: &mut Budget) -> Op<4> {
    Op::momentum(coords(), m, BigQ::from(1), b).unwrap()
}
fn sum(parts: Vec<(i64, Op<4>)>, b: &mut Budget) -> Op<4> {
    let mut out = Op::from_terms(coords(), BTreeMap::new(), b).unwrap();
    for (n, op) in parts {
        let factor = imag(n, b);
        out = out.try_add(&op.try_scale(&factor, b).unwrap(), b).unwrap();
    }
    out
}
#[test]
fn every_lorentz_and_poincare_component() {
    for m in 0..4 {
        for n in 0..4 {
            let b = &mut Budget::default();
            assert_eq!(
                l(m, n, b),
                l(n, m, b).try_scale(&constant(-1, b), b).unwrap()
            );
            for r in 0..4 {
                let b = &mut Budget::default();
                let lhs = l(m, n, b).commutator(&p(r, b), b).unwrap();
                let rhs = sum(vec![(eta(m, r), p(n, b)), (-eta(n, r), p(m, b))], b);
                assert_eq!(lhs, rhs, "[L{m}{n},P{r}]");
                for s in 0..4 {
                    let b = &mut Budget::default();
                    let lhs = l(m, n, b).commutator(&l(r, s, b), b).unwrap();
                    let rhs = sum(
                        vec![
                            (eta(m, r), l(n, s, b)),
                            (-eta(n, r), l(m, s, b)),
                            (eta(n, s), l(m, r, b)),
                            (-eta(m, s), l(n, r, b)),
                        ],
                        b,
                    );
                    assert_eq!(lhs, rhs, "[L{m}{n},L{r}{s}]");
                }
            }
            let b = &mut Budget::default();
            assert!(p(m, b).commutator(&p(n, b), b).unwrap().terms().is_empty());
        }
    }
}
#[test]
fn arbitrary_polynomial_witness_and_generator_signs() {
    let b = &mut Budget::default();
    let f = monomial(&[2, 1, 3, 1], b)
        .try_add(&monomial(&[1, 2, 1, 2], b), b)
        .unwrap()
        .try_add(&constant(11, b), b)
        .unwrap();
    let h = BigQ::parse("7", "3", b).unwrap();
    for mu in 0..4 {
        for nu in 0..4 {
            let a = Op::orbital_lorentz(coords(), mu, nu, h.clone(), b).unwrap();
            let pm = Op::momentum(coords(), mu, h.clone(), b).unwrap();
            let pn = Op::momentum(coords(), nu, h.clone(), b).unwrap();
            let xm = BigPoly::try_var(format!("x{mu}"), b).unwrap();
            let xn = BigPoly::try_var(format!("x{nu}"), b).unwrap();
            let expected = xm
                .try_mul(&pn.apply(&f, b).unwrap(), b)
                .unwrap()
                .try_sub(&xn.try_mul(&pm.apply(&f, b).unwrap(), b).unwrap(), b)
                .unwrap();
            assert_eq!(a.apply(&f, b).unwrap(), expected);
        }
    }
    let b = &mut Budget::default();
    let x0 = BigPoly::try_var("x0", b).unwrap();
    let x1 = BigPoly::try_var("x1", b).unwrap();
    assert_eq!(p(0, b).apply(&x0, b).unwrap(), imag(1, b));
    assert_eq!(p(1, b).apply(&x1, b).unwrap(), imag(-1, b));
}
#[test]
fn invalid_contexts_indices_and_hard_caps() {
    let b = &mut Budget::default();
    assert!(Op::<0>::from_terms([], BTreeMap::new(), b).is_err());
    assert!(Op::<17>::from_terms(coords(), BTreeMap::new(), b).is_err());
    assert!(Op::<2>::partial(["x".into(), "x".into()], 0, b).is_err());
    assert!(Op::<1>::partial(["".into()], 0, b).is_err());
    assert!(Op::<1>::partial(coords(), 1, b).is_err());
    assert!(Op::<4>::momentum(coords(), 4, BigQ::from(1), b).is_err());
    assert!(Op::<4>::orbital_lorentz(coords(), 0, 4, BigQ::from(1), b).is_err());
    assert!(Op::<1>::from_terms(coords(), BTreeMap::from([([33], constant(1, b))]), b).is_err());
    let a = Op::<2>::partial(coords(), 0, b).unwrap();
    let other = Op::<2>::partial(["x1".into(), "x0".into()], 0, b).unwrap();
    assert!(a.try_add(&other, b).is_err());
    assert!(a.try_sub(&other, b).is_err());
    assert!(a.compose(&other, b).is_err());
    assert!(a.commutator(&other, b).is_err());
}
#[test]
fn budgets_validate_zero_paths_and_preflight_growth() {
    let b = &mut Budget::default();
    let x = monomial(&[2], b);
    let a = Op::<1>::partial(coords(), 0, b).unwrap();
    let zero = Op::<1>::from_terms(coords(), BTreeMap::new(), b).unwrap();
    let mut tiny = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(differentiate(&x, "unused", &mut tiny).is_err());
    assert!(a.try_scale(&BigPoly::default(), &mut tiny).is_err());
    let mut tiny = Budget::new(Limits {
        max_degree: 0,
        ..Limits::default()
    });
    assert!(zero.apply(&x, &mut tiny).is_err());
    assert!(a.compose(&zero, &mut tiny).is_err());
    let mut tiny = Budget::new(Limits {
        max_terms: 1,
        ..Limits::default()
    });
    assert!(matches!(
        Op::<1>::partial(coords(), 0, &mut tiny),
        Err(Error::ResourceLimit {
            resource: Resource::Terms,
            ..
        })
    ));
    let mut tiny = Budget::new(Limits {
        max_input_bytes: 1,
        ..Limits::default()
    });
    assert!(Op::<1>::partial(coords(), 0, &mut tiny).is_err());
    let mut tiny = Budget::new(Limits {
        max_word_length: 0,
        ..Limits::default()
    });
    assert!(Op::<1>::partial(coords(), 0, &mut tiny).is_err());
    let mut tiny = Budget::new(Limits {
        max_coefficient_bits: 0,
        ..Limits::default()
    });
    assert!(a.apply(&BigPoly::default(), &mut tiny).is_err());
    let b = &mut Budget::new(Limits {
        max_word_length: 32,
        ..Limits::default()
    });
    let d = Op::<1>::from_terms(coords(), BTreeMap::from([([32], constant(1, b))]), b).unwrap();
    assert!(d.compose(&a, b).is_err());
    let normalized =
        Op::<1>::from_terms(coords(), BTreeMap::from([([1], BigPoly::default())]), b).unwrap();
    assert!(normalized.terms().is_empty());
}

#[test]
fn aggregate_bounds_and_checked_coefficient_growth() {
    let b = &mut Budget::default();
    let c = constant(1, b);
    let terms = BTreeMap::from([([0], c.clone()), ([1], c)]);
    // Two operator nodes and two coefficient monomials exceed three nodes.
    let mut small = Budget::new(Limits {
        max_terms: 3,
        ..Limits::default()
    });
    assert!(matches!(
        Op::<1>::from_terms(coords(), terms, &mut small),
        Err(Error::ResourceLimit {
            resource: Resource::Terms,
            ..
        })
    ));
    let square = monomial(&[2], b);
    let mut small = Budget::new(Limits {
        max_input_bytes: 5,
        ..Limits::default()
    });
    // Four coefficient symbol bytes plus the two coordinate bytes.
    assert!(matches!(
        Op::<1>::multiplication(coords(), square, &mut small),
        Err(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        })
    ));
    let x = monomial(&[1], b);
    let d = Op::<1>::partial(coords(), 0, b).unwrap();
    let m = Op::<1>::multiplication(coords(), x, b).unwrap();
    let before = m.clone();
    // Inputs each fit exactly; the commuted product needs two derivative slots.
    let mut small = Budget::new(Limits {
        max_terms: 3,
        ..Limits::default()
    });
    assert!(d.compose(&m, &mut small).is_err());
    assert_eq!(m, before);
    let high = BigQ::parse("255", "1", b).unwrap();
    let factor = BigPoly::try_constant(BigC::new(high, BigQ::from(0)), b).unwrap();
    let op = Op::<1>::multiplication(coords(), factor.clone(), b).unwrap();
    let mut small = Budget::new(Limits {
        max_coefficient_bits: 8,
        ..Limits::default()
    });
    assert!(matches!(
        op.try_scale(&factor, &mut small),
        Err(Error::ResourceLimit {
            resource: Resource::CoefficientBits,
            ..
        })
    ));
    let mut small = Budget::new(Limits {
        max_operations: 1,
        ..Limits::default()
    });
    assert!(op.try_scale(&factor, &mut small).is_err());
    assert!(small.used() <= 1);
    // A zero hbar must never excuse malformed coordinates or indices.
    assert!(Op::<4>::momentum(
        ["z".into(), "z".into(), "a".into(), "b".into()],
        0,
        BigQ::from(0),
        b
    )
    .is_err());
    assert!(Op::<4>::orbital_lorentz(coords(), 4, 4, BigQ::from(0), b).is_err());
}
