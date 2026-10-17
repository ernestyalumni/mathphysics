use lie_groups::{
    default_tolerance, Complex, Half, Matrix, Real, Scalar, So3, So3Algebra, Su2, Su2Algebra,
};
fn v<T: Real>(a: [f64; 3]) -> [T; 3] {
    a.map(T::from_f64)
}
fn close<T: Real>(a: T, b: f64, tol: f64) {
    assert!(
        (a.to_f64() - b).abs() <= tol,
        "{} versus {b} (tol {tol})",
        a.to_f64()
    );
}
fn run<T: Real>() {
    let tol = default_tolerance::<T>().to_f64() * 4.0;
    let angle = T::from_f64(std::f64::consts::FRAC_PI_2);
    let r = So3::from_axis_angle(v::<T>([0., 0., 1.]), angle).unwrap();
    let image = r.apply(&v([1., 0., 0.])).unwrap();
    for (actual, want) in image.into_iter().zip([0., 1., 0.]) {
        close(actual, want, tol);
    }
    let z = So3Algebra::<T>::new(v([0., 0., 0.])).exp().unwrap();
    assert_eq!(z.matrix(), &Matrix::identity());
    for tiny in [1e-3, 1e-6, 1e-10] {
        let theta = T::from_f64(tiny);
        let r = So3Algebra::new([T::zero(), T::zero(), theta])
            .exp()
            .unwrap();
        close(
            r.matrix()[(1, 0)],
            theta.to_f64().sin(),
            tol.max(tiny * 1e-6),
        );
    }
    for sample in 0..32 {
        let x = [
            0.11 * (sample as f64 - 10.0),
            0.07 * (sample as f64 + 1.0),
            -0.03 * (sample as f64 + 3.0),
        ];
        let y = [-x[1], x[2], x[0]];
        let a = So3Algebra::new(v::<T>(x)).exp().unwrap();
        let b = So3Algebra::new(v::<T>(y)).exp().unwrap();
        let c = a.compose(&b).unwrap();
        let backwards = c.inverse().unwrap().compose(&c).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(backwards.matrix()[(i, j)], f64::from(i == j), tol);
            }
        }
        let u = Su2Algebra::new(v::<T>(x)).exp().unwrap();
        let covered = u.to_so3().unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(covered.matrix()[(i, j)], a.matrix()[(i, j)].to_f64(), tol);
            }
        }
        let inverse = u.inverse().compose(&u).unwrap();
        close(inverse.quaternion()[0], 1., tol);
        for q in &inverse.quaternion()[1..] {
            close(*q, 0., tol);
        }
        let antipodal = u.antipode().to_so3().unwrap();
        assert_eq!(covered.matrix(), antipodal.matrix());
        let round = Su2::try_from_matrix(u.matrix(), default_tolerance::<T>()).unwrap();
        for (actual, want) in round.quaternion().iter().zip(u.quaternion()) {
            close(*actual, want.to_f64(), tol);
        }
    }
    let u = Su2Algebra::new(v::<T>([0., 0., std::f64::consts::TAU]))
        .exp()
        .unwrap();
    close(u.quaternion()[0], -1., tol);
    assert!(
        u.quaternion()[0] < T::zero(),
        "SU2 must distinguish -I from +I"
    );
    let full = Su2::exp_pauli_over_i(v::<T>([0., 0., std::f64::consts::FRAC_PI_2])).unwrap();
    close(full.quaternion()[3], 1., tol);
    let (a, b) = (r, So3::from_axis_angle(v([0., 1., 0.]), angle).unwrap());
    let image = a.compose(&b).unwrap().apply(&v([0., 0., 1.])).unwrap();
    for (actual, want) in image.into_iter().zip([0., 1., 0.]) {
        close(actual, want, tol);
    }
    assert!(So3::<T>::try_from_matrix(Matrix::zero(), default_tolerance()).is_err());
    assert!(So3::<T>::try_from_matrix(
        Matrix::from_rows([
            [T::from_f64(-1.), T::zero(), T::zero()],
            [T::zero(), T::one(), T::zero()],
            [T::zero(), T::zero(), T::one()]
        ]),
        default_tolerance()
    )
    .is_err());
    assert!(So3::from_axis_angle(v::<T>([0., 0., 0.]), T::zero()).is_err());
    assert!(Su2::<T>::try_from_unit_quaternion([T::zero(); 4], default_tolerance()).is_err());
    assert!(So3::<T>::try_from_matrix(Matrix::identity(), T::from_f64(-0.1)).is_err());
    assert!(So3::<T>::try_from_matrix(Matrix::identity(), T::one()).is_err());
    let bad = Matrix::from_rows([
        [Complex::real(T::one()), Complex::zero()],
        [Complex::zero(), Complex::real(-T::one())],
    ]);
    assert!(Su2::try_from_matrix(bad, default_tolerance()).is_err());
}
#[test]
fn f32_uses_the_same_generic_api() {
    run::<f32>();
}
#[test]
fn f64_uses_the_same_generic_api() {
    run::<f64>();
}
#[test]
fn binary16_uses_the_same_generic_api() {
    run::<Half>();
}
#[test]
fn approximate_group_inverse_revalidates_instead_of_leaking_an_invalid_value() {
    let m = Matrix::from_rows([
        [0.625_f64, 0.625, 0.625],
        [0.6875, -0.6875, 0.0],
        [0.375, 0.375, -0.75],
    ]);
    let admitted = So3::try_from_matrix(m, 0.12).unwrap();
    assert!(So3::try_from_matrix(m.transpose(), 0.1249).is_err());
    assert!(admitted.inverse().is_err());
}
#[test]
fn nonfinite_values_and_tolerances_fail_closed() {
    let mut bad = Matrix::<f64, 3, 3>::identity();
    bad[(0, 0)] = f64::NAN;
    assert!(So3::try_from_matrix(bad, 1e-10).is_err());
    assert!(So3::<f64>::try_from_matrix(Matrix::identity(), f64::NAN).is_err());
    assert!(So3Algebra::new([f64::INFINITY, 0., 0.]).exp().is_err());
    assert!(Su2Algebra::new([0., f64::NAN, 0.]).exp().is_err());
    assert!(Su2::<f64>::try_from_unit_quaternion([f64::INFINITY, 0., 0., 0.], 1e-10).is_err());
    assert!(So3::<f64>::identity().apply(&[f64::NAN, 0., 0.]).is_err());
}
#[test]
fn generators_square_differently_and_small_angles_have_correct_derivative() {
    let r = So3Algebra::new([0., 0., 1e-8_f64]).exp().unwrap();
    close(r.matrix()[(1, 0)] / 1e-8, 1., 1e-14);
    let u = Su2Algebra::new([0., 0., 1e-8_f64]).exp().unwrap();
    close(u.quaternion()[3] / 1e-8, 0.5, 1e-14);
}
