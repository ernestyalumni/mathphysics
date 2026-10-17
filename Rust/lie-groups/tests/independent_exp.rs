use lie_groups::{Complex, Matrix, Scalar, So3Algebra, Su2Algebra};
#[test]
fn closed_exponentials_match_independent_matrix_power_series() {
    for angle in [
        0.0_f64,
        1e-12,
        0.1249,
        0.125,
        0.1251,
        -0.25,
        1.,
        -2.,
        std::f64::consts::PI,
    ] {
        let vector = [0.6 * angle, 0.8 * angle, 0.];
        let algebra = So3Algebra::new(vector);
        let (mut term, mut sum) = (Matrix::identity(), Matrix::identity());
        for k in 1..=80 {
            term = (&term * &algebra.matrix()).scaled(1.0 / (k as f64));
            sum = &sum + &term;
        }
        let closed = algebra.exp().unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((sum[(i, j)] - closed.matrix()[(i, j)]).abs() < 2e-14);
            }
        }
        let algebra = Su2Algebra::new(vector);
        let (mut term, mut sum) = (Matrix::identity(), Matrix::identity());
        for k in 1..=80 {
            term = (&term * &algebra.matrix()).scaled(Complex::real(1.0 / (k as f64)));
            sum = &sum + &term;
        }
        let closed = algebra.exp().unwrap().matrix();
        for i in 0..2 {
            for j in 0..2 {
                let residual = sum[(i, j)] - closed[(i, j)];
                assert!(residual.re.abs() < 2e-14 && residual.im.abs() < 2e-14);
            }
        }
    }
    let x = Complex::new(qft_symbolics::scalar::C::I, qft_symbolics::scalar::C::ONE);
    assert_eq!(x.conj().re, -qft_symbolics::scalar::C::I);
    assert_eq!(Complex::<f64>::one().norm_squared(), 1.);
}
