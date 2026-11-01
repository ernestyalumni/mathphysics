use qft_symbolics::{
    checked,
    computation::{Budget, Limits},
    exact::{BigC, BigQ},
    relativity::*,
};
fn q(n: i64, d: i64) -> BigQ {
    BigQ::new(n.into(), d.into(), &mut Budget::default()).unwrap()
}
// Absolute component tolerance for moderate boosts; diagnostic tolerances are explicit below.
const COMPONENT_TOLERANCE: f64 = 1e-12;
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < COMPONENT_TOLERANCE, "{a} != {b}");
}
fn product(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..4).map(|k| a[i][k] * b[k][j]).sum()))
}
#[test]
fn rapidity_group_and_velocity_policy() {
    let a = Rapidity::from_beta(0.6).unwrap();
    close(a.eta(), 2.0_f64.ln());
    close(a.beta().unwrap(), 0.6);
    close(a.collinear_add(a.inverse()).unwrap().eta(), 0.0);
    close(a.collinear_add(a).unwrap().beta().unwrap(), 15.0 / 17.0);
    assert!(Rapidity::new(30.0).unwrap().beta().is_err());
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(Rapidity::new(x).is_err());
        assert!(Rapidity::from_beta(x).is_err());
    }
    for x in [-1.0, 1.0, 1.1] {
        assert!(Rapidity::from_beta(x).is_err());
    }
    assert!(Rapidity::new(f64::MAX)
        .unwrap()
        .collinear_add(Rapidity::new(f64::MAX).unwrap())
        .is_err());
}
#[test]
fn numeric_sign_inverse_interval_and_composition() {
    let axis = UnitAxis::axial(Axis::X);
    let eta = Rapidity::from_beta(0.6).unwrap();
    let l = boost_f64(eta, axis).unwrap();
    close(l[0][0], 1.25);
    close(l[0][1], 0.75);
    close(l[1][0], 0.75);
    let v = [5.0, 3.0, 2.0, -1.0];
    let w = apply_f64(&l, &v).unwrap();
    close(w[0], 8.5);
    close(w[1], 7.5);
    let opposite = apply_f64(&boost_f64(eta.inverse(), axis).unwrap(), &v).unwrap();
    close(opposite[0], 4.0);
    assert_ne!(opposite, w);
    close(dot_f64(&v, &v).unwrap(), -11.0);
    close(dot_f64(&w, &w).unwrap(), -11.0);
    let back = apply_f64(&boost_f64(eta.inverse(), axis).unwrap(), &w).unwrap();
    for i in 0..4 {
        close(v[i], back[i]);
    }
    let squared = product(&l, &l);
    let doubled = boost_f64(eta.collinear_add(eta).unwrap(), axis).unwrap();
    for i in 0..4 {
        for j in 0..4 {
            close(squared[i][j], doubled[i][j]);
        }
    }
    let tol = Tolerance {
        absolute: 1e-13,
        relative: 1e-13,
    };
    assert!(metric_residual_f64(&l, tol).unwrap().within_tolerance);
    let mut bad = l;
    bad[0][1] = -bad[0][1];
    assert!(!metric_residual_f64(&bad, tol).unwrap().within_tolerance);
}
#[test]
fn general_axis_and_tiny_eta() {
    let n = UnitAxis::new([0.6, 0.8, 0.0], 1e-15).unwrap();
    let l = boost_f64(Rapidity::new(0.7).unwrap(), n).unwrap();
    assert!(
        metric_residual_f64(
            &l,
            Tolerance {
                absolute: 1e-14,
                relative: 1e-14
            }
        )
        .unwrap()
        .within_tolerance
    );
    let tiny = boost_f64(Rapidity::new(1e-10).unwrap(), n).unwrap();
    close(tiny[0][1], 6e-11);
    assert!(tiny[1][2] > 0.0);
    assert!((tiny[1][2] - 2.4e-21).abs() < 1e-35);
}
#[test]
fn numeric_rejection_boundaries() {
    for n in [
        [0.0; 3],
        [2.0, 0.0, 0.0],
        [f64::NAN, 0.0, 0.0],
        [f64::MAX; 3],
    ] {
        assert!(UnitAxis::new(n, 1e-12).is_err());
    }
    for t in [-1.0, 1.0, f64::NAN] {
        assert!(UnitAxis::new([1.0, 0.0, 0.0], t).is_err());
    }
    assert!(boost_f64(Rapidity::new(1000.0).unwrap(), UnitAxis::axial(Axis::X)).is_err());
    let huge = [[f64::MAX; 4]; 4];
    assert!(metric_residual_f64(
        &huge,
        Tolerance {
            absolute: 0.0,
            relative: 0.0
        }
    )
    .is_err());
    assert!(apply_f64(&huge, &[2.0; 4]).is_err());
    assert!(dot_f64(&[f64::MAX; 4], &[2.0; 4]).is_err());
    let id = boost_f64(Rapidity::new(0.0).unwrap(), UnitAxis::axial(Axis::Z)).unwrap();
    assert!(metric_residual_f64(
        &id,
        Tolerance {
            absolute: f64::NAN,
            relative: 0.0
        }
    )
    .is_err());
    assert!(apply_f64(&id, &[f64::NAN; 4]).is_err());
}
#[test]
fn exact_345_and_metric_vector_contracts() {
    let mut b = Budget::default();
    let l = exact_half_velocity_boost(&q(1, 2), Axis::X, &mut b).unwrap();
    assert_eq!(l[0][0], BigC::new(q(5, 3), 0.into()));
    assert_eq!(l[0][1], BigC::new(q(4, 3), 0.into()));
    let v = [5, 3, 2, -1].map(BigC::from);
    let low = lower(&v, &mut b).unwrap();
    assert_eq!(low, [-5, 3, 2, -1].map(BigC::from));
    assert_eq!(raise(&low, &mut b).unwrap(), v);
    let w = apply(&l, &v, &mut b).unwrap();
    assert_eq!(dot(&v, &v, &mut b).unwrap(), BigC::from(-11));
    assert_eq!(dot(&w, &w, &mut b).unwrap(), BigC::from(-11));
    let inv = exact_half_velocity_boost(&q(-1, 2), Axis::X, &mut b).unwrap();
    assert_eq!(
        checked::mul(&l, &inv, &mut b).unwrap(),
        checked::identity::<BigC, 4>(&mut b).unwrap()
    );
    let twice = exact_half_velocity_boost(&q(4, 5), Axis::X, &mut b).unwrap();
    assert_eq!(checked::mul(&l, &l, &mut b).unwrap(), twice);
    let g = metric::<BigC>(&mut b).unwrap();
    let lt = checked::transpose(&l, &mut b).unwrap();
    assert_eq!(
        checked::mul(&checked::mul(&lt, &g, &mut b).unwrap(), &l, &mut b).unwrap(),
        g
    );
    for u in [q(1, 1), q(-1, 1), q(3, 2)] {
        assert!(exact_half_velocity_boost(&u, Axis::Y, &mut b).is_err());
    }
}
#[test]
fn symbolic_condition_is_retained_not_assumed() {
    use std::collections::BTreeMap;
    let mut b = Budget::default();
    let l = symbolic_axial_boost("c", "s", Axis::Z, &mut b).unwrap();
    assert!(!l.condition.is_zero());
    let good = BTreeMap::from([
        ("c".into(), BigC::new(q(5, 3), 0.into())),
        ("s".into(), BigC::new(q(4, 3), 0.into())),
    ]);
    assert!(l.condition.try_eval(&good, &mut b).unwrap().is_zero());
    let bad = BTreeMap::from([("c".into(), BigC::from(2)), ("s".into(), BigC::from(1))]);
    assert_eq!(l.condition.try_eval(&bad, &mut b).unwrap(), BigC::from(2));
    assert_eq!(
        l.matrix[0][3].try_eval(&good, &mut b).unwrap(),
        BigC::new(q(4, 3), 0.into())
    );
    assert!(symbolic_axial_boost("c", "c", Axis::X, &mut b).is_err());
}
#[test]
fn budget_rejections() {
    fn empty() -> Budget {
        Budget::new(Limits {
            max_operations: 0,
            ..Limits::default()
        })
    }
    let v = [1, 2, 3, 4].map(BigC::from);
    assert!(lower(&v, &mut empty()).is_err());
    assert!(raise(&v, &mut empty()).is_err());
    assert!(dot(&v, &v, &mut empty()).is_err());
    assert!(metric::<BigC>(&mut empty()).is_err());
    assert!(exact_half_velocity_boost(&q(1, 2), Axis::X, &mut empty()).is_err());
    assert!(symbolic_axial_boost("c", "s", Axis::X, &mut empty()).is_err());
    let m = checked::identity::<BigC, 4>(&mut Budget::default()).unwrap();
    assert!(apply(&m, &v, &mut empty()).is_err());
    let mut bits = Budget::new(Limits {
        max_coefficient_bits: 1,
        ..Limits::default()
    });
    assert!(lower(&v, &mut bits).is_err());
}

#[test]
fn parent_demo_exact_and_numeric_agree() {
    let mut b = Budget::default();
    let exact = exact_half_velocity_boost(&q(1, 3), Axis::Z, &mut b).unwrap();
    let v = [2, 0, 0, 3].map(BigC::from);
    let w = apply(&exact, &v, &mut b).unwrap();
    assert_eq!(
        w,
        [
            BigC::new(q(19, 4), 0.into()),
            0.into(),
            0.into(),
            BigC::new(q(21, 4), 0.into())
        ]
    );
    assert_eq!(dot(&w, &w, &mut b).unwrap(), BigC::from(5));
    let numeric = boost_f64(Rapidity::from_beta(0.6).unwrap(), UnitAxis::axial(Axis::Z)).unwrap();
    let numeric_w = apply_f64(&numeric, &[2.0, 0.0, 0.0, 3.0]).unwrap();
    for (actual, expected) in numeric_w
        .into_iter()
        .zip([19.0 / 4.0, 0.0, 0.0, 21.0 / 4.0])
    {
        close(actual, expected);
    }
    close(dot_f64(&numeric_w, &numeric_w).unwrap(), 5.0);
}
