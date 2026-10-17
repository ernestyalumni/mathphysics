use lie_groups::{
    levi_civita, su2,
    vector::{cross, dot},
    Complex, Matrix, Scalar, So3Algebra, Su2Algebra, Su2Basis,
};
use qft_symbolics::{polynomial::Poly, scalar::Q};
fn unit<T: Scalar>(i: usize) -> [T; 3] {
    std::array::from_fn(|j| T::from_i64(i64::from(i == j)))
}
#[test]
fn generic_matrices_are_not_mislabelled_rotations() {
    let a = Matrix::from_rows([[12i32, 8, 4], [3, 17, 14], [9, 8, 10]]);
    let b = Matrix::from_rows([[5, 19, 3], [6, 15, 9], [7, 8, 16]]);
    assert_eq!(
        &a * &b,
        Matrix::from_rows([[136, 380, 172], [215, 424, 386], [163, 371, 259]])
    );
    assert_eq!(a[(1, 2)], 14);
    let rectangular = Matrix::from_rows([[1i64, 2, 3], [4, 5, 6]]);
    assert_eq!(
        &rectangular * &rectangular.transpose(),
        Matrix::from_rows([[14, 32], [32, 77]])
    );
    assert_eq!(Matrix::<Q, 3, 3>::identity().determinant(), Q::ONE);
}
#[test]
fn frankel_basis_all_nine_commutators_and_trace_normalization() {
    let expected = [
        [[0, 0, 0], [0, 0, -1], [0, 1, 0]],
        [[0, 0, 1], [0, 0, 0], [-1, 0, 0]],
        [[0, -1, 0], [1, 0, 0], [0, 0, 0]],
    ];
    for (i, want) in expected.iter().enumerate() {
        let e = So3Algebra::<Q>::basis(i).unwrap().matrix();
        assert_eq!(e, Matrix::from_rows(want.map(|r| r.map(Q::from))));
        assert_eq!(e.transpose(), e.scaled(-Q::ONE));
        for j in 0..3 {
            let f = So3Algebra::<Q>::basis(j).unwrap().matrix();
            let lhs = e.commutator(&f);
            let mut rhs = Matrix::<Q, 3, 3>::zero();
            for k in 0..3 {
                rhs = &rhs
                    + &So3Algebra::<Q>::basis(k)
                        .unwrap()
                        .matrix()
                        .scaled(Q::from(i64::from(levi_civita([i, j, k]).unwrap())));
            }
            assert_eq!(lhs, rhs);
            assert_eq!((&e * &f).trace(), Q::from(if i == j { -2 } else { 0 }));
        }
    }
    assert!(So3Algebra::<Q>::basis(3).is_err());
}
#[test]
fn hat_cross_bracket_and_jacobi_are_symbolic_polynomial_identities() {
    let u = std::array::from_fn(|i| Poly::var(format!("u{i}")));
    let v = std::array::from_fn(|i| Poly::var(format!("v{i}")));
    let w = std::array::from_fn(|i| Poly::var(format!("w{i}")));
    let a = So3Algebra::new(u.clone());
    let b = So3Algebra::new(v.clone());
    let c = So3Algebra::new(w);
    assert_eq!(a.matrix().mul_vector(&v), cross(&u, &v));
    assert_eq!(a.matrix().commutator(&b.matrix()), a.bracket(&b).matrix());
    let first = a.bracket(&b.bracket(&c));
    let other = b.bracket(&c.bracket(&a));
    let last = c.bracket(&a.bracket(&b));
    for ((x, y), z) in first
        .coordinates()
        .iter()
        .zip(other.coordinates())
        .zip(last.coordinates())
    {
        assert!((x.clone() + y.clone() + z.clone()).is_zero());
    }
    let expected = Matrix::from_rows(std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            v[i].clone() * u[j].clone()
                - if i == j {
                    dot(&u, &v)
                } else {
                    Poly::constant(0)
                }
        })
    }));
    assert_eq!(&a.matrix() * &b.matrix(), expected);
    assert_eq!(So3Algebra::try_from_matrix(a.matrix()).unwrap(), a);
    let a = Su2Algebra::new(u);
    let b = Su2Algebra::new(v);
    assert_eq!(
        a.matrix().adjoint(),
        a.matrix().scaled(Complex::real(Poly::constant(-1)))
    );
    assert_eq!(a.matrix().commutator(&b.matrix()), a.bracket(&b).matrix());
    assert_eq!(Su2Algebra::try_from_matrix(a.matrix()).unwrap(), a);
}
#[test]
fn skew_matrix_constructor_rejects_symmetric_parts_and_bad_diagonal() {
    assert!(So3Algebra::<Q>::try_from_matrix(Matrix::identity()).is_err());
    let wrong = Matrix::from_rows([
        [Q::ZERO, Q::ONE, Q::ZERO],
        [Q::ONE, Q::ZERO, Q::ZERO],
        [Q::ZERO; 3],
    ]);
    assert!(So3Algebra::try_from_matrix(wrong).is_err());
}
#[test]
fn pauli_and_both_su2_normalizations_are_distinct_and_correct() {
    let i = Complex::new(Q::ZERO, Q::ONE);
    for a in 0..3 {
        let p = su2::pauli::<Q>(a).unwrap();
        assert_eq!(p.adjoint(), p);
        assert_eq!((&p * &p), Matrix::identity());
        assert_eq!(p.determinant(), Complex::real(-Q::ONE));
        let t = su2::generator::<Q>(a, Su2Basis::HalfPauliOverI).unwrap();
        let f = su2::generator::<Q>(a, Su2Basis::PauliOverI).unwrap();
        assert_eq!(f, t.scaled(Complex::real(Q::from(2))));
        assert_eq!(t.adjoint(), t.scaled(Complex::real(-Q::ONE)));
        assert_eq!(t.trace(), Complex::zero());
        assert_ne!(p, t);
        for b in 0..3 {
            let q = su2::pauli::<Q>(b).unwrap();
            let mut target = Matrix::zero();
            for k in 0..3 {
                target = &target
                    + &su2::pauli::<Q>(k).unwrap().scaled(
                        i * Complex::real(Q::from(2 * i64::from(levi_civita([a, b, k]).unwrap()))),
                    );
            }
            assert_eq!(p.commutator(&q), target);
            let tb = su2::generator::<Q>(b, Su2Basis::HalfPauliOverI).unwrap();
            let fb = su2::generator::<Q>(b, Su2Basis::PauliOverI).unwrap();
            let mut target = Matrix::zero();
            for k in 0..3 {
                target = &target
                    + &su2::generator::<Q>(k, Su2Basis::HalfPauliOverI)
                        .unwrap()
                        .scaled(Complex::real(Q::from(i64::from(
                            levi_civita([a, b, k]).unwrap(),
                        ))));
            }
            assert_eq!(t.commutator(&tb), target);
            assert_eq!(f.commutator(&fb), target.scaled(Complex::real(Q::from(4))));
        }
    }
    let a = Su2Algebra::<Q>::new([Q::from(2), Q::from(3), Q::from(-1)]);
    let b = Su2Algebra::<Q>::new([Q::from(-2), Q::from(1), Q::from(4)]);
    assert_eq!(a.matrix().commutator(&b.matrix()), a.bracket(&b).matrix());
    assert_eq!(a.to_so3().bracket(&b.to_so3()), a.bracket(&b).to_so3());
    assert!(su2::pauli::<Q>(3).is_err());
    assert_eq!(Su2Algebra::try_from_matrix(a.matrix()).unwrap(), a);
    assert!(Su2Algebra::try_from_matrix(su2::pauli::<Q>(2).unwrap()).is_err());
    assert!(Su2Algebra::<Q>::basis(3).is_err());
}
#[test]
fn generic_vector_epsilon_agrees_with_cross_product() {
    let a = [2i32, 3, 5];
    let b = [7, 11, 13];
    for i in 0..3 {
        let epsilon_sum: i32 = (0..3)
            .flat_map(|j| (0..3).map(move |k| (j, k)))
            .map(|(j, k)| i32::from(levi_civita([i, j, k]).unwrap()) * a[j] * b[k])
            .sum();
        assert_eq!(cross(&a, &b)[i], epsilon_sum);
    }
    assert_eq!(cross(&unit::<i32>(0), &unit(1)), unit(2));
}
