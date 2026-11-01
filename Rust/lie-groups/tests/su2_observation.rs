use lie_groups::su2_observation::{observe, Chart, Error, Slice};
use std::f64::consts::PI;
type C = [f64; 2];
fn add(a: C, b: C) -> C {
    [a[0] + b[0], a[1] + b[1]]
}
fn mul(a: C, b: C) -> C {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn conj(a: C) -> C {
    [a[0], -a[1]]
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-12, "{a} != {b}");
}
// Independent Pauli expansion, not quaternion helper or production coordinates.
fn oracle(n: [f64; 3], t: f64) -> [[C; 2]; 2] {
    let pauli = [
        [[[0., 0.], [1., 0.]], [[1., 0.], [0., 0.]]],
        [[[0., 0.], [0., -1.]], [[0., 1.], [0., 0.]]],
        [[[1., 0.], [0., 0.]], [[0., 0.], [-1., 0.]]],
    ];
    std::array::from_fn(|r| {
        std::array::from_fn(|c| {
            let mut z = [if r == c { (t / 2.).cos() } else { 0. }, 0.];
            for k in 0..3 {
                z = add(z, mul([0., -n[k] * (t / 2.).sin()], pauli[k][r][c]));
            }
            z
        })
    })
}
#[test]
fn public_matrix_and_bloch_against_independent_pauli_spinor() {
    let d = 29f64.sqrt();
    for n in [
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [2. / d, -3. / d, 4. / d],
    ] {
        for t in [
            0.,
            PI / 2.,
            PI,
            2. * PI,
            4. * PI,
            -PI / 2.,
            -7.31,
            1e-8,
            1e6,
        ] {
            let o = observe(n, t, 0.25).unwrap();
            let m = oracle(n, t);
            for (actual, expected) in o.matrix.iter().flatten().zip(m.iter().flatten()) {
                near(actual[0], expected[0]);
                near(actual[1], expected[1]);
            }
            // Evolve |+x>=(1,1)/sqrt(2) then evaluate all three Pauli means.
            let psi = m.map(|row| add(row[0], row[1]).map(|x| x / 2f64.sqrt()));
            let cross = mul(conj(psi[0]), psi[1]);
            let bloch = [
                2. * cross[0],
                2. * cross[1],
                mul(conj(psi[0]), psi[0])[0] - mul(conj(psi[1]), psi[1])[0],
            ];
            for (a, b) in o.bloch.into_iter().zip(bloch) {
                near(a, b);
            }
            assert!(o.checks.unitarity_error < 2e-14);
            assert!(o.checks.determinant_error < 2e-14);
        }
    }
}
#[test]
fn sign_periodicity_and_positive_z_orientation() {
    let n = [0., 0., 1.];
    let a = observe(n, 0.73, 0.).unwrap();
    let b = observe(n, 0.73 + 2. * PI, 0.).unwrap();
    for (x, y) in a.coordinates.into_iter().zip(b.coordinates) {
        near(x, -y);
    }
    for (x, y) in a.bloch.into_iter().zip(b.bloch) {
        near(x, y);
    }
    let q = observe(n, PI / 2., 0.).unwrap();
    near(q.bloch[1], 1.);
    near(observe(n, 2. * PI, 0.).unwrap().coordinates[0], -1.);
    near(observe(n, 4. * PI, 0.).unwrap().coordinates[0], 1.);
}
#[test]
fn chart_near_poles_and_subnormal_classification() {
    let n = [0., 0., 1.];
    for t in [
        0.,
        -0.,
        4. * PI,
        -4. * PI,
        8. * PI,
        f64::from_bits(1),
        -f64::from_bits(1),
        1e-310,
    ] {
        let o = observe(n, t, 0.).unwrap();
        assert_eq!(o.stereographic, Chart::Pole);
        assert!(o.coordinates.iter().all(|x| x.is_finite()));
    }
    for t in [
        1e-9,
        -1e-9,
        1e-300,
        -1e-300,
        4. * PI + 1e-8,
        4. * PI - 1e-8,
        2. * PI,
        -2. * PI,
    ] {
        let Chart::Finite(p) = observe(n, t, 0.).unwrap().stereographic else {
            panic!("{t}")
        };
        assert!(p.iter().all(|x| x.is_finite()));
        // Inverse stereographic map is independent of cot and stable here.
        let r = p[2].abs();
        if r < 1e150 {
            let recovered = (r * r - 1.) / (r * r + 1.);
            near(recovered, (t / 2.).cos());
        }
        if t.abs() < 1e-8 {
            assert_eq!(p[2].is_sign_negative(), t.is_sign_positive());
        }
    }
}
#[test]
fn original_axis_domain_without_normalization_and_residuals() {
    for drift in [-9e-13, 9e-13] {
        let axis = [1. + drift, 0., 0.];
        let o = observe(axis, PI, 0.).unwrap();
        assert_eq!(o.axis, axis);
        near(o.coordinates[1], -axis[0]);
        let expected = (axis[0] * axis[0] - 1.).abs();
        assert!((o.checks.unitarity_error - expected).abs() < 1e-15);
        assert!((o.checks.determinant_error - expected).abs() < 1e-15);
        assert!(o.checks.unitarity_error > 1e-12);
    }
    for axis in [
        [0., 0., 0.],
        [1. + 2e-12, 0., 0.],
        [1. - 2e-12, 0., 0.],
        [f64::MAX, f64::MAX, 0.],
    ] {
        assert_eq!(observe(axis, 0., 0.), Err(Error::UnitAxisRequired));
    }
    for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            observe([x, 0., 1.], 0., 0.),
            Err(Error::FiniteNumberRequired)
        );
        assert_eq!(
            observe([0., 0., 1.], x, 0.),
            Err(Error::FiniteNumberRequired)
        );
        assert_eq!(
            observe([0., 0., 1.], 0., x),
            Err(Error::FiniteNumberRequired)
        );
    }
    for t in [-1e6, 1e6] {
        assert!(observe([0., 0., 1.], t, 0.).is_ok());
    }
    for t in [-1e6 - 1., 1e6 + 1.] {
        assert_eq!(observe([0., 0., 1.], t, 0.), Err(Error::TOutOfRange));
    }
}
#[test]
fn slice_topology_and_stable_radius() {
    for c in [-f64::MAX, -1.01, 1.01, f64::MAX] {
        assert_eq!(
            observe([1., 0., 0.], 0., c).unwrap().slice,
            Slice::Empty { coordinate: c }
        );
    }
    for c in [-1., 1.] {
        assert_eq!(
            observe([1., 0., 0.], 0., c).unwrap().slice,
            Slice::Point { coordinate: c }
        );
    }
    for c in [0., 0.3, -0.3, 1. - f64::EPSILON, -1. + f64::EPSILON] {
        let Slice::Sphere { coordinate, radius } = observe([1., 0., 0.], 0., c).unwrap().slice
        else {
            panic!()
        };
        assert_eq!(coordinate, c);
        near(radius * radius + c * c, 1.);
        assert!(radius > 0.);
    }
}

#[test]
fn explicit_minus_i_conversion_matches_existing_algebra_on_unit_axes() {
    use lie_groups::Su2Algebra;
    for axis in [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]] {
        for t in [-4. * PI, -PI, -1e-9, 0., 1e-9, PI, 2. * PI, 4. * PI] {
            let legacy = Su2Algebra::new(axis.map(|n| n * t)).exp().unwrap();
            let o = observe(axis, t, 0.).unwrap();
            let q = legacy.quaternion();
            near(o.coordinates[0], q[0]);
            for i in 0..3 {
                near(o.coordinates[i + 1], -q[i + 1]);
            }
            for i in 0..2 {
                for j in 0..2 {
                    near(o.matrix[i][j][0], legacy.matrix()[(i, j)].re);
                    near(o.matrix[i][j][1], legacy.matrix()[(i, j)].im);
                }
            }
        }
    }
}
