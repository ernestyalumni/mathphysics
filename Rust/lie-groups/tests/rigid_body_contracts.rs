use lie_groups::{
    goldstein::{euler_lift, EulerConvention},
    rigid_body::{Model, State, StepOptions, Wrench},
    vector, Matrix, So3, Su2,
};
type Mat = Matrix<f64, 3, 3>;
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn diag(a: f64, b: f64, c: f64) -> Mat {
    Matrix::from_rows([[a, 0., 0.], [0., b, 0.], [0., 0., c]])
}
fn offdiag() -> Mat {
    Matrix::from_rows([[2., 0.2, 0.1], [0.2, 3., 0.3], [0.1, 0.3, 4.]])
}
fn close(a: [f64; 3], b: [f64; 3], tol: f64) {
    for i in 0..3 {
        assert!((a[i] - b[i]).abs() <= tol, "{a:?} != {b:?}");
    }
}
fn error(a: &Mat, b: &Mat) -> f64 {
    (0..3)
        .flat_map(|i| (0..3).map(move |j| (a[(i, j)] - b[(i, j)]).abs()))
        .fold(0., f64::max)
}
fn evolve(m: &Model, mut s: State, dt: f64, n: usize, w: Wrench) -> State {
    for _ in 0..n {
        s = m.step(s, dt, w, StepOptions::default()).unwrap().state;
    }
    s
}
fn attitude(s: State) -> Mat {
    *Su2::try_from_unit_quaternion(s.quaternion, 1e-10)
        .unwrap()
        .to_so3()
        .unwrap()
        .matrix()
}
#[test]
fn free_translation_is_exact() {
    let m = Model::new(2., diag(1., 2., 2.5)).unwrap();
    let initial = State {
        position: [0.1, -0.2, 0.3],
        velocity: [0.3, 0.1, -0.2],
        ..State::default()
    };
    let s = evolve(&m, initial, 0.01, 100, Wrench::default());
    close(
        s.position,
        std::array::from_fn(|i| initial.position[i] + initial.velocity[i]),
        1e-14,
    );
    assert_eq!(s.velocity, initial.velocity);
    assert_eq!(s.quaternion, initial.quaternion);
}
#[test]
fn spherical_and_principal_rotation_match_discrete_exact_solution() {
    for sphere in [false, true] {
        let m = Model::new(
            2.,
            if sphere {
                diag(2., 2., 2.)
            } else {
                diag(1., 2., 2.5)
            },
        )
        .unwrap();
        let initial = State {
            omega_body: if sphere {
                [0.2, -0.4, 0.7]
            } else {
                [0., 0., 0.7]
            },
            ..State::default()
        };
        let w = initial.omega_body;
        let n = norm(w);
        let s = evolve(&m, initial, 0.01, 100, Wrench::default());
        close(s.omega_body, w, 2e-14);
        let expected = So3::from_axis_angle(w, 200. * (0.01 * n / 2.).atan()).unwrap();
        assert!(error(&attitude(s), expected.matrix()) < 5e-14);
    }
}
#[test]
fn full_spd_inertia_conserves_energy_and_inertial_momentum_vector() {
    let m = Model::new(2., offdiag()).unwrap();
    let mut s = State {
        omega_body: [0.3, -0.4, 0.8],
        velocity: [0.2, -0.1, 0.3],
        quaternion: *euler_lift(EulerConvention::Zxz, [0.4, 0.7, -0.3])
            .unwrap()
            .quaternion(),
        ..State::default()
    };
    let h = m.momentum_inertial(s).unwrap();
    let energy = m.kinetic_energy(s).unwrap();
    let hnorm = norm(m.momentum_body(s).unwrap());
    for _ in 0..2000 {
        s = m
            .step(
                s,
                0.005,
                Wrench::default(),
                StepOptions {
                    tolerance: 1e-14,
                    max_iterations: 20,
                },
            )
            .unwrap()
            .state;
        close(m.momentum_inertial(s).unwrap(), h, 2e-10);
        assert!((m.kinetic_energy(s).unwrap() - energy).abs() < 2e-10);
        assert!((norm(m.momentum_body(s).unwrap()) - hnorm).abs() < 2e-10);
    }
    close(
        m.solve_inertia(m.inertia().mul_vector(&[0.3, -0.7, 1.1]))
            .unwrap(),
        [0.3, -0.7, 1.1],
        2e-15,
    );
}
#[test]
fn symmetric_top_full_attitude_reference_and_second_order() {
    let m = Model::new(1., diag(2., 2., 3.)).unwrap();
    let s = State {
        omega_body: [0.3, -0.5, 0.8],
        quaternion: *euler_lift(EulerConvention::Zyx, [0.4, -0.2, 0.3])
            .unwrap()
            .quaternion(),
        ..State::default()
    };
    let h = m.momentum_inertial(s).unwrap();
    let a = 0.4_f64;
    let left = So3::from_axis_angle(h, norm(h) / 2.).unwrap();
    let right = So3::from_axis_angle([0., 0., 1.], -a).unwrap();
    let exact = &(left.matrix() * &attitude(s)) * right.matrix();
    let mut previous = 0.;
    for n in [25, 50, 100] {
        let out = evolve(&m, s, 1. / f64::from(n), n as usize, Wrench::default());
        let e = error(&attitude(out), &exact);
        close(
            out.omega_body,
            [
                0.3 * a.cos() + 0.5 * a.sin(),
                0.3 * a.sin() - 0.5 * a.cos(),
                0.8,
            ],
            2e-5,
        );
        if previous > 0. {
            assert!(
                previous / e > 3.8 && previous / e < 4.2,
                "second order ratio {}",
                previous / e
            );
        }
        previous = e;
    }
}
#[test]
fn intermediate_axis_flips_without_external_wrench() {
    let m = Model::new(1., diag(1., 2., 2.5)).unwrap();
    let initial = State {
        omega_body: [0.01, 1., 0.01],
        ..State::default()
    };
    let h = m.momentum_inertial(initial).unwrap();
    let energy = m.kinetic_energy(initial).unwrap();
    let s = evolve(&m, initial, 0.005, 4800, Wrench::default());
    assert!(s.omega_body[1] < -0.99, "unforced intermediate-axis flip");
    close(m.momentum_inertial(s).unwrap(), h, 2e-9);
    assert!((m.kinetic_energy(s).unwrap() - energy).abs() < 2e-9);
}
#[test]
fn body_inertial_force_and_rhs_signs() {
    let m = Model::new(2., diag(2., 2., 2.)).unwrap();
    let mut w = Wrench {
        force_body: [2., -1., 0.5],
        ..Wrench::default()
    };
    let s = evolve(&m, State::default(), 0.01, 100, w);
    close(s.velocity, [1., -0.5, 0.25], 2e-14);
    close(s.position, [0.5, -0.25, 0.125], 2e-14);
    w = Wrench {
        torque_body: [0.1, -0.2, 0.3],
        ..Wrench::default()
    };
    let s = evolve(&m, State::default(), 0.01, 100, w);
    close(s.omega_body, [0.05, -0.1, 0.15], 2e-14);
    let mut tilted = State {
        quaternion: *euler_lift(EulerConvention::Zyx, [std::f64::consts::FRAC_PI_2, 0., 0.])
            .unwrap()
            .quaternion(),
        ..State::default()
    };
    w = Wrench {
        force_body: [2., 0., 0.],
        ..Wrench::default()
    };
    close(m.rhs(tilted, w).unwrap().velocity, [0., 1., 0.], 1e-14);
    tilted.omega_body = [0.2, -0.4, 0.7];
    w = Wrench {
        force_inertial: [2., -1., 0.5],
        ..Wrench::default()
    };
    let s = evolve(&m, tilted, 0.01, 100, w);
    close(s.velocity, [1., -0.5, 0.25], 2e-14);
    close(s.position, [0.5, -0.25, 0.125], 2e-14);
    let m = Model::new(1., diag(1., 2., 2.5)).unwrap();
    let s = State {
        omega_body: [1., 2., 3.],
        ..State::default()
    };
    let d = m.rhs(s, Wrench::default()).unwrap();
    close(d.omega_body, [-3., 2.25, -0.8], 2e-14);
    assert_eq!(d.quaternion, [0., 0.5, 1., 1.5]);
}
#[test]
fn signed_lift_equivariance_and_midpoint_work() {
    let m = Model::new(2., offdiag()).unwrap();
    let mut a = State {
        quaternion: *euler_lift(EulerConvention::Zxz, [0.2, 0.5, -0.4])
            .unwrap()
            .quaternion(),
        omega_body: [0.3, -0.4, 0.8],
        velocity: [0.2, 0.1, -0.3],
        ..State::default()
    };
    let mut b = a;
    b.quaternion = b.quaternion.map(|x| -x);
    let w = Wrench {
        force_body: [0.2, -0.1, 0.3],
        torque_body: [0.03, -0.02, 0.01],
        ..Wrench::default()
    };
    let dt = 0.02;
    let out = m.step(a, dt, w, StepOptions::default()).unwrap().state;
    let vm = std::array::from_fn(|i| (a.velocity[i] + out.velocity[i]) / 2.);
    let wm = std::array::from_fn(|i| (a.omega_body[i] + out.omega_body[i]) / 2.);
    let mut qm = std::array::from_fn(|i| a.quaternion[i] + out.quaternion[i]);
    let n = f64::hypot(f64::hypot(qm[0], qm[1]), f64::hypot(qm[2], qm[3]));
    qm = qm.map(|x| x / n);
    let force = Su2::try_from_unit_quaternion(qm, 1e-10)
        .unwrap()
        .to_so3()
        .unwrap()
        .apply(&w.force_body)
        .unwrap();
    let work = dt * (vector::dot(&vm, &force) + vector::dot(&wm, &w.torque_body));
    assert!((m.kinetic_energy(out).unwrap() - m.kinetic_energy(a).unwrap() - work).abs() < 2e-12);
    a = evolve(&m, a, 0.002, 100, w);
    b = evolve(&m, b, 0.002, 100, w);
    assert_eq!(a.position, b.position);
    assert_eq!(a.velocity, b.velocity);
    assert_eq!(a.omega_body, b.omega_body);
    for i in 0..4 {
        assert!((a.quaternion[i] + b.quaternion[i]).abs() < 1e-15);
    }
}
#[test]
fn invalid_imports_steps_and_newton_failure_are_errors() {
    let j = diag(1., 2., 2.5);
    for mass in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(Model::new(mass, j).is_err());
    }
    for inertia in [diag(1., 1., 3.), diag(0., 1., 1.), diag(-1., 2., 2.)] {
        assert!(Model::new(1., inertia).is_err());
    }
    Model::new(1., diag(1., 1., 2.)).unwrap();
    assert!(Model::new(
        1.,
        Matrix::from_rows([[1., 2., 0.], [2., 1., 0.], [0., 0., 1.]])
    )
    .is_err());
    let rotation = *euler_lift(EulerConvention::Zyx, [0.4, 0.7, -0.3])
        .unwrap()
        .to_so3()
        .unwrap()
        .matrix();
    assert!(Model::new(1., &(&rotation * &diag(1., 1., 3.)) * &rotation.transpose()).is_err());
    Model::new(1., &(&rotation * &diag(1., 1., 2.)) * &rotation.transpose()).unwrap();
    let mut bad = j;
    bad[(0, 1)] = 0.1;
    assert!(Model::new(1., bad).is_err());
    bad = j;
    bad[(0, 0)] = f64::NAN;
    assert!(Model::new(1., bad).is_err());
    let m = Model::new(1., j).unwrap();
    let mut s = State::default();
    for dt in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(m
            .step(s, dt, Wrench::default(), StepOptions::default())
            .is_err());
    }
    s.quaternion = [2., 0., 0., 0.];
    assert!(m.rhs(s, Wrench::default()).is_err());
    assert!(m
        .step(s, 0.01, Wrench::default(), StepOptions::default())
        .is_err());
    s = State::default();
    s.velocity[0] = f64::NAN;
    assert!(m.rhs(s, Wrench::default()).is_err());
    s = State::default();
    let w = Wrench {
        torque_body: [f64::INFINITY, 0., 0.],
        ..Wrench::default()
    };
    assert!(m.step(s, 0.01, w, StepOptions::default()).is_err());
    for options in [
        StepOptions {
            tolerance: 0.,
            max_iterations: 20,
        },
        StepOptions {
            tolerance: 1e-12,
            max_iterations: 0,
        },
        StepOptions {
            tolerance: 1e-3,
            max_iterations: 20,
        },
        StepOptions {
            tolerance: 1e-12,
            max_iterations: 101,
        },
    ] {
        assert!(m.step(s, 0.01, Wrench::default(), options).is_err());
    }
    s.omega_body = [3., 2., 4.];
    assert!(m
        .step(
            s,
            1.2,
            Wrench::default(),
            StepOptions {
                tolerance: 1e-14,
                max_iterations: 1
            }
        )
        .is_err());
    s = State::default();
    s.quaternion[0] += 1e-12;
    assert_eq!(
        m.step(s, 0.01, Wrench::default(), StepOptions::default())
            .unwrap()
            .state
            .quaternion,
        [1., 0., 0., 0.]
    );
}
