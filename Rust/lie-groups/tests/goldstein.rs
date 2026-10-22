use lie_groups::{
    goldstein::*, multitrig::MultiTrig, su2, Complex, Error, Half, Matrix, RationalScalar, Scalar,
    So3, So3Algebra, Su2,
};
use qft_symbolics::{polynomial::Poly, scalar::Q};
use std::f64::consts::{FRAC_PI_2, PI};
const CONVENTIONS: [EulerConvention; 3] = [
    EulerConvention::Zxz,
    EulerConvention::Zyz,
    EulerConvention::Zyx,
];
fn close(a: &Matrix<f64, 3, 3>, b: &Matrix<f64, 3, 3>, tol: f64) {
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (a[(i, j)] - b[(i, j)]).abs() < tol,
                "{i},{j}: {} != {}",
                a[(i, j)],
                b[(i, j)]
            );
        }
    }
}
fn rotate(mut v: [f64; 3], axis: usize, angle: f64) -> [f64; 3] {
    let j = (axis + 1) % 3;
    let k = (axis + 2) % 3;
    let (s, c) = angle.sin_cos();
    let a = v[j];
    let b = v[k];
    v[j] = c * a - s * b;
    v[k] = s * a + c * b;
    v
}
#[test]
fn explicit_steps_independently_match_euler_products() {
    for convention in CONVENTIONS {
        for i in 0..129 {
            let angles = [
                (i as f64 - 64.0) * 0.117,
                0.023 + i as f64 * 0.081,
                0.72 - i as f64 * 0.057,
            ];
            let r = from_euler(convention, angles).unwrap();
            let [p, t, s] = angles;
            let steps = match convention {
                EulerConvention::Zxz => [(2, s), (0, t), (2, p)],
                EulerConvention::Zyz => [(2, s), (1, t), (2, p)],
                EulerConvention::Zyx => [(0, s), (1, t), (2, p)],
            };
            for j in 0..3 {
                let mut v = [0.0; 3];
                v[j] = 1.0;
                for (a, x) in steps {
                    v = rotate(v, a, x);
                }
                let actual = r
                    .apply(&std::array::from_fn(|k| if k == j { 1.0 } else { 0.0 }))
                    .unwrap();
                for k in 0..3 {
                    assert!((actual[k] - v[k]).abs() < 3e-15)
                }
            }
            close(
                r.matrix(),
                euler_lift(convention, angles)
                    .unwrap()
                    .to_so3()
                    .unwrap()
                    .matrix(),
                3e-15,
            );
        }
    }
}
#[test]
fn principal_euler_extraction_regular_and_singular() {
    for convention in CONVENTIONS {
        for theta in [
            -PI,
            -FRAC_PI_2,
            -1.0,
            -1e-10,
            0.0,
            1e-10,
            0.8,
            FRAC_PI_2,
            PI - 1e-10,
            PI,
            PI + 1e-10,
        ] {
            for phi in [-2.3, 0.0, 1.2] {
                let r = from_euler(convention, [phi, theta, 0.73]).unwrap();
                let out = extract_euler(&r, convention);
                let r2 = from_euler(convention, out.angles).unwrap();
                close(r.matrix(), r2.matrix(), 1e-13);
                if convention == EulerConvention::Zyx {
                    assert!(out.angles[1].abs() <= FRAC_PI_2 + 1e-15);
                    if out.singular {
                        assert_eq!(out.angles[2], 0.0);
                    }
                } else {
                    assert!((0.0..=PI).contains(&out.angles[1]));
                    if out.singular {
                        assert_eq!(out.angles[2], 0.0);
                    }
                }
            }
        }
    }
    assert!(
        extract_euler(
            &from_euler(EulerConvention::Zxz, [0.5, 0.0, 1.1]).unwrap(),
            EulerConvention::Zxz
        )
        .singular
    );
    assert!(
        extract_euler(
            &from_euler(EulerConvention::Zyx, [0.5, FRAC_PI_2, 1.1]).unwrap(),
            EulerConvention::Zyx
        )
        .singular
    );
    assert!(
        !extract_euler(
            &from_euler(EulerConvention::Zyx, [0.5, FRAC_PI_2 - 1e-8, 1.1]).unwrap(),
            EulerConvention::Zyx
        )
        .singular
    );
}
#[test]
fn axis_angle_survives_zero_pi_and_all_pivots() {
    for axis in [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, -2.0, 3.0],
        [-1.0, -1.0, -1.0],
    ] {
        for angle in [
            0.0,
            1e-18,
            1e-12,
            0.23,
            PI - 1e-9,
            PI,
            PI + 1e-9,
            2.0 * PI,
            6.8,
        ] {
            let r = So3::from_axis_angle(axis, angle).unwrap();
            let aa = axis_angle(&r).unwrap();
            assert!((0.0..=PI).contains(&aa.angle));
            let reconstructed = So3::from_axis_angle(aa.axis, aa.angle).unwrap();
            close(r.matrix(), reconstructed.matrix(), 4e-15);
            close(
                r.matrix(),
                lift_rotation(&r).unwrap().to_so3().unwrap().matrix(),
                4e-15,
            );
            assert!((r.matrix().trace() - (1.0 + 2.0 * aa.angle.cos())).abs() < 3e-15);
            if (0.0..1e-10).contains(&angle) && angle > 0.0 {
                assert!((aa.angle - angle).abs() < angle * 1e-14)
            }
        }
    }
    assert_eq!(
        axis_angle(&So3::identity()).unwrap(),
        AxisAngle {
            axis: [1.0, 0.0, 0.0],
            angle: 0.0
        }
    );
}
#[test]
fn cayley_klein_retains_lifts_and_passive_dictionary() {
    let q = euler_lift(EulerConvention::Zxz, [0.53, 1.2, -0.71]).unwrap();
    let ck = CayleyKlein::from_lift(q);
    let p = ck.parameters();
    assert_eq!(p[2], -p[1].conj());
    assert_eq!(p[3], p[0].conj());
    assert_eq!(ck.matrix(), q.matrix().adjoint());
    let recovered = CayleyKlein::try_from_parameters(p, 1e-14).unwrap();
    assert_eq!(recovered.lift().quaternion(), q.quaternion());
    let opposite = CayleyKlein::from_lift(q.antipode());
    assert_ne!(opposite, ck);
    assert_eq!(opposite.parameters(), p.map(|z| -z));
    close(
        &ck.goldstein_matrix().unwrap(),
        &q.to_so3().unwrap().matrix().transpose(),
        1e-15,
    );
    close(
        &ck.goldstein_matrix().unwrap(),
        &opposite.goldstein_matrix().unwrap(),
        1e-15,
    );
    let x = [0.8, -1.7, 0.3];
    let y = ck.goldstein_matrix().unwrap().mul_vector(&x);
    let k = ck.matrix();
    let sigma = |v: [f64; 3]| {
        let mut result = Matrix::zero();
        for (i, x) in v.into_iter().enumerate() {
            result = &result + &su2::pauli(i).unwrap().scaled(Complex::real(x));
        }
        result
    };
    let lhs = &(&k * &sigma(x)) * &k.adjoint();
    let rhs = sigma(y);
    for i in 0..2 {
        for j in 0..2 {
            assert!((lhs[(i, j)].re - rhs[(i, j)].re).abs() < 2e-15);
            assert!((lhs[(i, j)].im - rhs[(i, j)].im).abs() < 2e-15)
        }
    }
    let mut wrong = p;
    wrong[2] = p[1].conj();
    assert!(CayleyKlein::try_from_parameters(wrong, 1e-14).is_err());
}
#[test]
fn antipodal_lift_tracks_two_pi() {
    let a = euler_lift(EulerConvention::Zxz, [0.2, 0.9, 0.4]).unwrap();
    let b = euler_lift(EulerConvention::Zxz, [0.2 + 2.0 * PI, 0.9, 0.4]).unwrap();
    for i in 0..4 {
        assert!((a.quaternion()[i] + b.quaternion()[i]).abs() < 5e-16)
    }
    close(
        a.to_so3().unwrap().matrix(),
        b.to_so3().unwrap().matrix(),
        2e-15,
    );
}
#[test]
fn rejects_nonfinite_nonrotations_and_bad_parameters() {
    assert!(matches!(
        from_euler(EulerConvention::Zxz, [f64::NAN, 0.0, 0.0]),
        Err(Error::NonFinite)
    ));
    assert!(euler_lift(EulerConvention::Zyx, [0.0, f64::INFINITY, 0.0]).is_err());
    assert!(axis_rotation(3, 1.0, 0.0).is_err());
    assert!(So3::try_from_matrix(
        Matrix::from_rows([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        1e-14
    )
    .is_err());
    let p = [Complex::new(0.0, 0.0); 4];
    assert!(CayleyKlein::try_from_parameters(p, 1e-14).is_err());
    assert!(Su2::try_from_unit_quaternion([2.0, 0.0, 0.0, 0.0], 1e-14).is_err());
}
#[test]
fn forward_scalar_backends_are_generic() {
    for convention in CONVENTIONS {
        let r = from_euler(convention, [0.2_f64, 0.7, -1.2]).unwrap();
        let f = from_euler(convention, [0.2_f32, 0.7, -1.2]).unwrap();
        let h = from_euler(
            convention,
            [
                Half::from_f64(0.2),
                Half::from_f64(0.7),
                Half::from_f64(-1.2),
            ],
        )
        .unwrap();
        close(r.matrix(), &f.matrix().map(|x| f64::from(*x)), 3e-7);
        close(r.matrix(), &h.matrix().map(|x| x.to_f64()), 0.004);
    }
}
type T = MultiTrig<3>;
fn sc() -> ([T; 3], [T; 3]) {
    (
        std::array::from_fn(|i| T::sin(i).unwrap()),
        std::array::from_fn(|i| T::cos(i).unwrap()),
    )
}
#[test]
fn exact_multitrig_ring_domains_and_evaluation() {
    let (s, c) = sc();
    for i in 0..3 {
        assert_eq!(
            s[i].clone() * s[i].clone() + c[i].clone() * c[i].clone(),
            T::one()
        );
    }
    let s0 = Q::new(3, 5);
    let c0 = Q::new(4, 5);
    let ss = [s0, Q::ZERO, Q::ONE];
    let cc = [c0, Q::ONE, Q::ZERO];
    assert_eq!(
        (s[0].clone() * c[2].clone() + c[0].clone() * s[2].clone())
            .at(ss, cc)
            .unwrap(),
        c0
    );
    assert!(s[0].at([Q::ONE; 3], [Q::ONE; 3]).is_err());
    assert!(T::sin(3).is_err());
    assert!(MultiTrig::<9>::cos(0).is_err());
    let a = s[0].clone() + c[1].clone();
    let b = c[0].clone() - s[2].clone();
    let d = c[2].clone() + s[1].clone();
    assert_eq!((a.clone() * b.clone()) * d.clone(), a * (b * d));
}
#[test]
fn exact_three_angle_euler_orthogonality_determinant_ch() {
    let (s, c) = sc();
    for convention in CONVENTIONS {
        let r = euler_matrix(convention, c.clone(), s.clone());
        assert_eq!(&r.transpose() * &r, Matrix::identity());
        assert_eq!(r.determinant(), T::one());
        assert_eq!(cayley_hamilton_residual(&r), Matrix::zero());
        let v = [T::from_i64(2), T::from_i64(-3), T::one()];
        assert_eq!(
            orthogonal_conjugate(&r, &So3Algebra::new(v.clone()).matrix()),
            So3Algebra::new(r.mul_vector(&v)).matrix()
        );
    }
}
#[test]
fn exact_x_y_euler_convention_dictionary() {
    let (s, c) = sc();
    let zxz = euler_matrix(EulerConvention::Zxz, c.clone(), s.clone());
    let zyz = euler_matrix(
        EulerConvention::Zyz,
        [s[0].clone(), c[1].clone(), -s[2].clone()],
        [-c[0].clone(), s[1].clone(), c[2].clone()],
    );
    assert_eq!(zxz, zyz);
}
#[test]
fn exact_half_angles_quaternion_cover_matches_euler() {
    let (s, c) = sc();
    let zero = T::zero();
    for convention in CONVENTIONS {
        let unit = |axis: usize, i: usize| {
            let mut q = [c[i].clone(), zero.clone(), zero.clone(), zero.clone()];
            q[axis + 1] = s[i].clone();
            q
        };
        let (a, b, d) = match convention {
            EulerConvention::Zxz => (unit(2, 0), unit(0, 1), unit(2, 2)),
            EulerConvention::Zyz => (unit(2, 0), unit(1, 1), unit(2, 2)),
            EulerConvention::Zyx => (unit(2, 0), unit(1, 1), unit(0, 2)),
        };
        let q = su2::hamilton(&su2::hamilton(&a, &b), &d);
        assert_eq!(
            q.iter()
                .cloned()
                .map(|x| x.clone() * x)
                .fold(T::zero(), |a, b| a + b),
            T::one()
        );
        let full_c =
            std::array::from_fn(|i| c[i].clone() * c[i].clone() - s[i].clone() * s[i].clone());
        let full_s = std::array::from_fn(|i| T::from_i64(2) * c[i].clone() * s[i].clone());
        assert_eq!(
            su2::quaternion_rotation(&q),
            euler_matrix(convention, full_c, full_s)
        );
        let k = cayley_klein_matrix(&q);
        assert_eq!(&k.adjoint() * &k, Matrix::identity());
        assert_eq!(k.determinant(), Complex::one());
        let p = cayley_klein_parameters(&q);
        assert_eq!(p[2], -p[1].conj());
        assert_eq!(p[3], p[0].conj());
        assert_eq!(k, su2::quaternion_matrix(&q).adjoint());
        assert_ne!(
            su2::quaternion_matrix(&q),
            su2::quaternion_matrix(&q.clone().map(|v| -v))
        );
    }
}
#[test]
fn exact_rodrigues_axis_and_arbitrary_matrix_ch() {
    let (s, c) = sc();
    let n = [T::from_q(Q::new(3, 5)), T::from_q(Q::new(4, 5)), T::zero()];
    let r = rodrigues(n.clone(), c[0].clone(), s[0].clone());
    assert_eq!(r.mul_vector(&n), n);
    assert_eq!(&r.transpose() * &r, Matrix::identity());
    assert_eq!(r.trace(), T::one() + T::from_i64(2) * c[0].clone());
    let a = Matrix::from_rows(std::array::from_fn(|i| {
        std::array::from_fn(|j| Poly::var(format!("a{i}{j}")))
    }));
    assert_eq!(cayley_hamilton_residual(&a), Matrix::zero());
    let wrong = &(&r.transpose() * &r) - &Matrix::identity().scaled(T::from_i64(-1));
    assert_ne!(wrong, Matrix::zero());
}

#[test]
fn reflection_is_not_the_so3_adjoint_vector_dictionary() {
    let b = Matrix::from_rows([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    let w = [2, 3, 5];
    let transformed = orthogonal_conjugate(&b, &So3Algebra::new(w).matrix());
    let naive = So3Algebra::new(b.mul_vector(&w)).matrix();
    assert_ne!(transformed, naive);
    assert_eq!(transformed, naive.scaled(-1));
}
