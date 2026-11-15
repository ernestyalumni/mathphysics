//! Oracles: C = 2 N^2 (SymPy for N = 2, 3 in the Yang-Mills ledger; derivation in the
//! mono-notes, Proposition ym:prop:casimir-sun), and K(X, Y) = 2 N tr(XY).
use lie_groups::special_unitary::{
    commutator, killing_form, negative_real_trace_form, special_unitary_basis,
    special_unitary_casimir_constant, ExactFraction,
};

#[test]
fn casimir_constant_is_twice_the_square_of_the_size() {
    for size in 2..=5usize {
        let expected = ExactFraction::integer(2 * (size * size) as i128);
        assert_eq!(special_unitary_casimir_constant(size), Some(expected), "su({size})");
    }
    assert_eq!(special_unitary_casimir_constant(2), Some(ExactFraction::integer(8)));
    assert_eq!(special_unitary_casimir_constant(3), Some(ExactFraction::integer(18)));
    assert!(special_unitary_casimir_constant(1).is_none());
    assert!(special_unitary_casimir_constant(0).is_none());
}

#[test]
fn killing_form_equals_two_n_times_the_trace_on_all_pairs() {
    for size in 2..=4usize {
        let basis = special_unitary_basis(size).unwrap();
        assert_eq!(basis.len(), size * size - 1);
        for left in &basis {
            for right in &basis {
                // tr(XY) = -B(X, Y) for skew-Hermitian X, Y, so K = -2 N B.
                let expected = -2 * size as i128 * i128::from(negative_real_trace_form(left, right));
                assert_eq!(killing_form(&basis, left, right), ExactFraction::integer(expected));
            }
        }
    }
}

#[test]
fn the_basis_is_orthogonal_traceless_skew_hermitian_and_closed_under_brackets() {
    for size in 2..=4usize {
        let basis = special_unitary_basis(size).unwrap();
        for (left_index, left) in basis.iter().enumerate() {
            let trace: i64 = (0..size).map(|index| left[index][index].imaginary).sum();
            assert_eq!(trace, 0);
            for row in 0..size {
                for column in 0..size {
                    assert_eq!(left[row][column].real, -left[column][row].real);
                    assert_eq!(left[row][column].imaginary, left[column][row].imaginary);
                }
            }
            for (right_index, right) in basis.iter().enumerate() {
                let form = negative_real_trace_form(left, right);
                assert_eq!(form == 0, left_index != right_index);
                let bracket = commutator(left, right);
                for row in 0..size {
                    for column in 0..size {
                        assert_eq!(bracket[row][column].real, -bracket[column][row].real);
                        assert_eq!(bracket[row][column].imaginary, bracket[column][row].imaginary);
                    }
                }
            }
        }
    }
}

#[test]
fn negative_control_a_basis_that_is_not_a_lie_algebra_basis_is_detected() {
    // Dropping one element makes the trace over the "basis" wrong: K is no longer 2 N tr.
    let mut basis = special_unitary_basis(2).unwrap();
    let full = killing_form(&basis, &basis[0].clone(), &basis[0].clone());
    basis.pop();
    let truncated = killing_form(&basis, &basis[0].clone(), &basis[0].clone());
    assert_ne!(full, truncated);
    assert!(ExactFraction::new(1, 0).is_none());
    assert_eq!(ExactFraction::new(-4, -6), ExactFraction::new(2, 3));
}
