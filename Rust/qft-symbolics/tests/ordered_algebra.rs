use qft_symbolics::{
    computation::{Budget, Error, Limits, Resource},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
};
// This test module also supports a standalone rustc harness before parent exports ordered.
use qft_symbolics::ordered::OrderedPolynomial as O;
use std::collections::BTreeMap;
fn g(s: &str, b: &mut Budget) -> O {
    O::generator(s, b).unwrap()
}
fn k(n: i64, b: &mut Budget) -> BigPoly {
    BigPoly::try_constant(n, b).unwrap()
}
#[test]
fn words_and_commutator_signs() {
    let b = &mut Budget::default();
    let a = g("A", b);
    let c = g("B", b);
    let ab = a.try_mul(&c, b).unwrap();
    let ba = c.try_mul(&a, b).unwrap();
    assert_ne!(ab, ba);
    assert!(ab.terms().contains_key(&vec!["A".into(), "B".into()]));
    let bracket = a.commutator(&c, b).unwrap();
    assert_eq!(bracket, ab.try_sub(&ba, b).unwrap());
    assert_eq!(bracket, c.commutator(&a, b).unwrap().try_neg(b).unwrap());
}
#[test]
fn jacobi_full_expansion() {
    let b = &mut Budget::default();
    let a = g("A", b);
    let c = g("B", b);
    let d = g("C", b);
    let x = a.commutator(&c, b).unwrap().commutator(&d, b).unwrap();
    assert_eq!(x.terms().len(), 4);
    for (word, sign) in [
        (["A", "B", "C"], 1),
        (["B", "A", "C"], -1),
        (["C", "A", "B"], -1),
        (["C", "B", "A"], 1),
    ] {
        assert_eq!(
            x.terms().get(&word.map(String::from).to_vec()),
            Some(&k(sign, b))
        );
    }
    let y = c.commutator(&d, b).unwrap().commutator(&a, b).unwrap();
    let z = d.commutator(&a, b).unwrap().commutator(&c, b).unwrap();
    assert!(x
        .try_add(&y, b)
        .unwrap()
        .try_add(&z, b)
        .unwrap()
        .terms()
        .is_empty());
}
#[test]
fn first_order_inverse_conjugation() {
    let b = &mut Budget::default();
    let a = g("observable", b);
    let j = g("generator", b);
    let eps = BigPoly::try_var("epsilon", b).unwrap();
    let one = O::scalar(k(1, b), b).unwrap();
    let ej = j.try_scale(&eps, b).unwrap();
    let u = one.try_add(&ej, b).unwrap();
    let inv = one.try_sub(&ej, b).unwrap();
    let transformed = inv
        .try_mul(&a, b)
        .unwrap()
        .try_mul(&u, b)
        .unwrap()
        .truncate(&["epsilon"], 1, b)
        .unwrap();
    assert_eq!(transformed.coefficient("epsilon", 0, b).unwrap(), a);
    assert_eq!(
        transformed.coefficient("epsilon", 1, b).unwrap(),
        a.commutator(&j, b).unwrap()
    );
    assert!(transformed
        .coefficient("epsilon", 2, b)
        .unwrap()
        .terms()
        .is_empty());
    assert_eq!(
        inv.try_mul(&u, b)
            .unwrap()
            .truncate(&["epsilon"], 1, b)
            .unwrap(),
        one
    );
}
#[test]
fn selected_total_degree_and_extraction() {
    let b = &mut Budget::default();
    let a = g("O", b);
    let x = BigPoly::try_var("x", b).unwrap();
    let y = BigPoly::try_var("y", b).unwrap();
    let p = x.try_mul(&y, b).unwrap();
    let v = a.try_scale(&p, b).unwrap();
    assert_eq!(v.truncate(&["x"], 1, b).unwrap(), v);
    assert!(v.truncate(&["x", "y"], 1, b).unwrap().terms().is_empty());
    assert_eq!(v.truncate(&["x", "x"], 1, b).unwrap(), v);
    assert_eq!(
        v.coefficient("x", 1, b).unwrap(),
        a.try_scale(&y, b).unwrap()
    );
}
#[test]
fn simultaneous_substitution_preserves_order() {
    let b = &mut Budget::default();
    let a = g("A", b);
    let c = g("B", b);
    let d = g("C", b);
    let map = BTreeMap::from([
        ("A".into(), c.try_add(&d, b).unwrap()),
        ("B".into(), a.clone()),
    ]);
    assert_eq!(
        a.try_mul(&c, b).unwrap().substitute(&map, b).unwrap(),
        c.try_add(&d, b).unwrap().try_mul(&a, b).unwrap()
    );
    assert_eq!(d.substitute(&map, b).unwrap(), d);
}
#[test]
fn invalid_names_are_rejected_even_for_zero() {
    let b = &mut Budget::default();
    let zero = O::default();
    for s in ["", "bad name", "bad\nname", "\0"] {
        assert!(O::generator(s, b).is_err());
        assert!(zero.coefficient(s, 0, b).is_err());
        assert!(zero.truncate(&[s], 0, b).is_err());
        assert!(zero
            .substitute(&BTreeMap::from([(s.into(), zero.clone())]), b)
            .is_err());
        let p = BigPoly::try_var(s, b).unwrap();
        assert!(O::scalar(p, b).is_err());
    }
}
#[test]
fn resource_bounds_validate_zero_operations_and_aggregate_storage() {
    let b = &mut Budget::default();
    let a = g("alpha", b);
    let c = g("bravo", b);
    let mut small = Budget::new(Limits {
        max_input_bytes: 9,
        ..Limits::default()
    });
    assert!(a.try_add(&c, &mut small).is_err());
    assert!(a.try_mul(&c, &mut small).is_err());
    let huge = BigPoly::try_constant(
        BigC::new(
            BigQ::parse("340282366920938463463374607431768211456", "1", b).unwrap(),
            0.into(),
        ),
        b,
    )
    .unwrap();
    let mut small = Budget::new(Limits {
        max_coefficient_bits: 8,
        ..Limits::default()
    });
    assert!(O::default().try_scale(&huge, &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    assert!(O::scalar(BigPoly::default(), &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_word_length: 0,
        ..Limits::default()
    });
    assert!(O::generator("A", &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_terms: 0,
        ..Limits::default()
    });
    assert!(O::generator("A", &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_input_bytes: 4,
        ..Limits::default()
    });
    assert!(a.try_mul(&O::default(), &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_terms: 1,
        ..Limits::default()
    });
    assert!(a.try_add(&c, &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_word_length: 1,
        ..Limits::default()
    });
    assert!(matches!(
        a.try_mul(&c, &mut small),
        Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            ..
        })
    ));
}
#[test]
fn arbitrary_precision_survives_ordered_products() {
    let b = &mut Budget::new(Limits {
        max_operations: 20_000_000,
        ..Limits::default()
    });
    let a = g("A", b);
    let q = BigQ::parse("340282366920938463463374607431768211457", "1", b).unwrap();
    let p = BigPoly::try_constant(BigC::new(q, 0.into()), b).unwrap();
    let x = a.try_scale(&p, b).unwrap().try_mul(&a, b).unwrap();
    assert_eq!(x.terms().get(&vec!["A".into(), "A".into()]), Some(&p));
    assert!(x.try_sub(&x, b).unwrap().terms().is_empty());
}

#[test]
fn aggregate_coefficients_replacements_and_hard_caps() {
    let b = &mut Budget::default();
    let a = g("A", b);
    let c = g("B", b);
    let x = BigPoly::try_var("x", b).unwrap();
    let y = BigPoly::try_var("y", b).unwrap();
    let p = x.try_add(&y, b).unwrap();
    let ap = a.try_scale(&p, b).unwrap();
    let cp = c.try_scale(&p, b).unwrap();
    let mut small = Budget::new(Limits {
        max_terms: 3,
        ..Limits::default()
    });
    assert!(ap.try_add(&cp, &mut small).is_err());
    let mut small = Budget::new(Limits {
        max_coefficient_bits: 1,
        ..Limits::default()
    });
    assert!(O::scalar(p, &mut small).is_err());
    let map = BTreeMap::from([("longkey".into(), g("longword", b))]);
    let mut small = Budget::new(Limits {
        max_input_bytes: 10,
        ..Limits::default()
    });
    assert!(O::default().substitute(&map, &mut small).is_err());
    let mut large = Budget::new(Limits {
        max_word_length: 1024,
        max_operations: 20_000_000,
        ..Limits::default()
    });
    let mut word = g("z", &mut large);
    for _ in 0..8 {
        word = word.try_mul(&word, &mut large).unwrap();
    }
    assert_eq!(word.terms().keys().next().unwrap().len(), 256);
    assert!(matches!(
        word.try_mul(&word, &mut large),
        Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            limit: 256,
            ..
        })
    ));
    let mut large = Budget::new(Limits {
        max_degree: 1024,
        max_operations: 20_000_000,
        ..Limits::default()
    });
    let poly = BigPoly::try_var("z", &mut large)
        .unwrap()
        .try_pow(257, &mut large)
        .unwrap();
    assert!(matches!(
        O::scalar(poly, &mut large),
        Err(Error::ResourceLimit {
            resource: Resource::Degree,
            limit: 256,
            ..
        })
    ));
}
