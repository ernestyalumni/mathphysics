use lie_groups::{
    symbolic::{
        so3_exponential, su2_exponential, Chart, RationalMatrix, SymbolicSo3, SymbolicSu2, Trig,
    },
    Complex, Matrix, RationalScalar, Scalar, So3Algebra, Su2Basis,
};
use qft_symbolics::{
    polynomial::Poly,
    scalar::{C, Q},
};
use std::collections::BTreeMap;
#[test]
fn one_angle_trigonometric_quotient_has_an_exact_normal_form() {
    let (s, c) = (Trig::sin(), Trig::cos());
    assert_eq!(s.clone() * s.clone() + c.clone() * c.clone(), Trig::one());
    assert_ne!(s.clone() * s.clone() - c.clone() * c.clone(), Trig::one());
    let sd = s.clone() * c.clone() * Trig::from_i64(2);
    let cd = c.clone() * c.clone() - s.clone() * s.clone();
    assert_eq!(sd.clone() * sd + cd.clone() * cd, Trig::one());
    assert_eq!(s.at(Q::new(3, 5), Q::new(4, 5)).unwrap(), Q::new(3, 5));
    assert!(s.at(Q::ONE, Q::ONE).is_err());
}
#[test]
fn symbolic_rodrigues_and_both_pauli_exponentials_have_correct_domains() {
    for axis in [
        [Q::ONE, Q::ZERO, Q::ZERO],
        [Q::ZERO, Q::ZERO, Q::ONE],
        [Q::new(3, 5), Q::new(4, 5), Q::ZERO],
    ] {
        let r = so3_exponential(axis).unwrap();
        assert_eq!(&r.transpose() * &r, Matrix::identity());
        assert_eq!(r.determinant(), Trig::one());
        let full = su2_exponential(axis, Su2Basis::PauliOverI).unwrap();
        let half = su2_exponential(axis, Su2Basis::HalfPauliOverI).unwrap();
        assert_eq!(full.argument_scale, Q::ONE);
        assert_eq!(half.argument_scale, Q::new(1, 2));
        assert_eq!(&full.matrix.adjoint() * &full.matrix, Matrix::identity());
        assert_eq!(full.matrix.determinant(), Complex::one());
        assert_eq!(full.matrix, half.matrix); // same formal circle, DIFFERENT argument metadata
    }
    let z = su2_exponential([Q::ZERO, Q::ZERO, Q::ONE], Su2Basis::PauliOverI).unwrap();
    assert_eq!(z.matrix[(0, 0)], Complex::new(Trig::cos(), -Trig::sin()));
    assert_eq!(z.matrix[(1, 1)], Complex::new(Trig::cos(), Trig::sin()));
    assert!(so3_exponential([Q::ONE; 3]).is_err());
    assert!(su2_exponential([Q::ZERO; 3], Su2Basis::PauliOverI).is_err());
}
#[test]
fn symbolic_arbitrary_chart_gives_genuine_groups_on_real_symbols() {
    let parameters = std::array::from_fn(|i| Poly::var(format!("t{i}")));
    let u = SymbolicSu2::stereographic(parameters.clone(), Chart::Plus).unwrap();
    assert!(u.norm_residual().is_zero());
    let m = u.matrix();
    assert!(m
        .multiply(&m.adjoint_real_symbols())
        .equals(&RationalMatrix::<2>::identity()));
    let det = m.numerator().determinant() - m.denominator().pow(2);
    assert!(det.is_zero());
    let r = SymbolicSo3::stereographic(parameters).unwrap();
    let m = r.matrix();
    assert!(m
        .multiply(&m.transpose())
        .equals(&RationalMatrix::<3>::identity()));
    assert!((m.numerator().determinant() - m.denominator().pow(3)).is_zero());
    assert!(u.compose(&u.inverse()).equals(&SymbolicSu2::identity()));
    assert!(!u.equals(&u.antipode()));
    assert!(SymbolicSo3::from_lift(u.clone()).equals(&SymbolicSo3::from_lift(u.antipode())));
}
#[test]
fn symbolic_chart_excludes_complex_parameter_domain_and_has_second_pole_chart() {
    assert!(SymbolicSu2::stereographic(
        [Poly::constant(C::I), Poly::constant(0), Poly::constant(0)],
        Chart::Plus
    )
    .is_err());
    let minus = SymbolicSu2::stereographic(
        [Poly::constant(0), Poly::constant(0), Poly::constant(0)],
        Chart::Minus,
    )
    .unwrap();
    assert!(minus.equals(&SymbolicSu2::identity().antipode()));
    let u = SymbolicSu2::stereographic(
        std::array::from_fn(|i| Poly::var(format!("t{i}"))),
        Chart::Plus,
    )
    .unwrap();
    assert!(u.evaluate(&BTreeMap::new()).is_err());
    let values = BTreeMap::from([
        ("t0".into(), C::I),
        ("t1".into(), C::ZERO),
        ("t2".into(), C::ZERO),
    ]);
    assert!(u.evaluate(&values).is_err());
    assert!(RationalMatrix::new(Matrix::<Poly, 2, 2>::zero(), Poly::constant(0)).is_err());
}
#[test]
fn symbolic_hamilton_norm_and_covering_homomorphism_are_independent_polynomial_identities() {
    let p = std::array::from_fn(|i| Poly::var(format!("p{i}")));
    let q = std::array::from_fn(|i| Poly::var(format!("q{i}")));
    let pq = lie_groups::su2::hamilton(&p, &q);
    let norm = |x: &[Poly; 4]| {
        x.iter()
            .fold(Poly::constant(0), |s, v| s + v.clone() * v.clone())
    };
    assert!((norm(&pq) - norm(&p) * norm(&q)).is_zero());
    let rotation = lie_groups::su2::quaternion_rotation;
    assert_eq!(rotation(&pq), &rotation(&p) * &rotation(&q));
    let matrix = lie_groups::su2::quaternion_matrix;
    assert_eq!(matrix(&pq), &matrix(&p) * &matrix(&q));
}
#[test]
fn rational_chart_exact_witness_matches_active_z_rotation_and_membership() {
    let u = SymbolicSu2::stereographic(
        [
            Poly::constant(0),
            Poly::constant(0),
            Poly::constant(Q::new(1, 2)),
        ],
        Chart::Plus,
    )
    .unwrap();
    let r = u.so3_matrix().evaluate(&BTreeMap::new()).unwrap();
    // q=(3/5,0,0,4/5): sin(theta)=24/25 and cos(theta)=-7/25.
    assert_eq!(r[(0, 0)], C::from(Q::new(-7, 25)));
    assert_eq!(r[(1, 0)], C::from(Q::new(24, 25)));
    assert_eq!(r.determinant(), C::ONE);
}
#[test]
fn half_angle_cover_matches_double_angle_rodrigues_exactly() {
    let (s, c) = (Trig::sin(), Trig::cos());
    let axis = [
        Trig::from_q(Q::new(3, 5)),
        Trig::from_q(Q::new(4, 5)),
        Trig::zero(),
    ];
    let q = [
        c.clone(),
        axis[0].clone() * s.clone(),
        axis[1].clone() * s.clone(),
        Trig::zero(),
    ];
    let projected = lie_groups::su2::quaternion_rotation(&q);
    let k = So3Algebra::new(axis).matrix();
    let expected = &(&Matrix::identity() + &k.scaled(Trig::from_i64(2) * s.clone() * c))
        + &(&k * &k).scaled(Trig::from_i64(2) * s.clone() * s);
    assert_eq!(projected, expected);
}
