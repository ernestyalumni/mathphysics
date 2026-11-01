use qft_symbolics::{
    clifford::*,
    computation::{Budget, Error, Limits, Resource},
    exact::BigC,
    matrices,
    polynomial::BigPoly,
};
use std::collections::BTreeMap;
#[test]
fn exact_trace_against_independent_pauli_block_matrices() {
    // Enumerate all component words through length 4, plus longer noncommuting words.
    let names = ["t", "x", "y", "z"];
    let mut words = vec![vec![]];
    for n in 1..=4 {
        for mut code in 0..4usize.pow(n) {
            let mut w = vec![];
            for _ in 0..n {
                w.push(code % 4);
                code /= 4;
            }
            words.push(w);
        }
    }
    words.extend([
        vec![0, 1, 2, 3, 0, 1],
        vec![2, 1, 3, 0, 2, 1, 3, 0],
        vec![0; 12],
    ]);
    let mut values = BTreeMap::new();
    // These are fixed component h=-eta values, independent of the trace recurrence.
    for i in 0..4 {
        for j in i..4 {
            let (a, c) = if names[i] <= names[j] {
                (names[i], names[j])
            } else {
                (names[j], names[i])
            };
            let key = format!("h(1:{a};1:{c})");
            values.insert(
                key,
                BigC::from(if i != j {
                    0
                } else if i == 0 {
                    1
                } else {
                    -1
                }),
            );
        }
    }
    for w in words {
        let labels: Vec<_> = w.iter().map(|&i| names[i]).collect();
        let mut b = Budget::new(Limits {
            max_operations: 100_000_000,
            ..Limits::default()
        });
        let actual = try_trace_exact(&labels, &mut b)
            .unwrap()
            .try_eval(&values, &mut b)
            .unwrap();
        let m = w.iter().fold(matrices::identity::<4>(), |m, &i| {
            matrices::mul(m, matrices::gamma(i))
        });
        assert_eq!(actual, BigC::from(matrices::trace(m)), "{w:?}");
    }
}
#[test]
fn hardcoded_symbolic_pairings_and_normal_order() {
    let mut b = Budget::default();
    let h = |a, c, b: &mut Budget| try_metric_with::<BigC>(a, c, b).unwrap();
    let ab = h("a", "b", &mut b);
    let cd = h("c", "d", &mut b);
    let ac = h("a", "c", &mut b);
    let bd = h("b", "d", &mut b);
    let ad = h("a", "d", &mut b);
    let bc = h("b", "c", &mut b);
    let expected = ab
        .try_mul(&cd, &mut b)
        .unwrap()
        .try_sub(&ac.try_mul(&bd, &mut b).unwrap(), &mut b)
        .unwrap()
        .try_add(&ad.try_mul(&bc, &mut b).unwrap(), &mut b)
        .unwrap()
        .try_mul(&BigPoly::constant(4), &mut b)
        .unwrap();
    assert_eq!(
        try_trace_exact(&["a", "b", "c", "d"], &mut b).unwrap(),
        expected
    );
    let ordered = try_normal_order_with::<BigC>(&["b", "a"], &mut b).unwrap();
    assert_eq!(
        ordered[&vec!["a".into(), "b".into()]],
        BigPoly::constant(-1)
    );
    assert_eq!(
        ordered[&Vec::<String>::new()],
        ab.try_mul(&BigPoly::constant(2), &mut b).unwrap()
    );
    assert_eq!(
        try_sandwich_srednicki_with::<BigC>(&[], &BigPoly::constant(4), &mut b).unwrap()[&vec![]],
        BigPoly::constant(-4)
    );
}
#[test]
fn checked_word_bytes_work_and_unsupported_domains() {
    let generous = Limits {
        max_word_length: 100,
        ..Limits::default()
    };
    assert!(matches!(
        try_trace_exact(&["a"; 14], &mut Budget::new(generous)),
        Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            limit: 12,
            ..
        })
    ));
    assert!(matches!(
        try_normal_order_with::<BigC>(&["a"; 13], &mut Budget::new(generous)),
        Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            ..
        })
    ));
    assert!(matches!(
        try_sandwich_srednicki_with::<BigC>(
            &["a"; 9],
            &BigPoly::constant(4),
            &mut Budget::new(generous)
        ),
        Err(Error::ResourceLimit {
            resource: Resource::WordLength,
            limit: 8,
            ..
        })
    ));
    // Odd words must still validate labels, even though their trace is zero.
    assert!(matches!(
        try_trace_exact(
            &["oversized"],
            &mut Budget::new(Limits {
                max_input_bytes: 2,
                ..Limits::default()
            })
        ),
        Err(Error::ResourceLimit {
            resource: Resource::InputBytes,
            ..
        })
    ));
    assert!(matches!(
        try_trace_exact(
            &[],
            &mut Budget::new(Limits {
                max_operations: 0,
                ..Limits::default()
            })
        ),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        })
    ));
    assert_ne!(metric("a;b", "c"), metric("a", "b;c"));
}
