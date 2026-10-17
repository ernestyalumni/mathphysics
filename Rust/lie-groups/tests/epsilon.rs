use lie_groups::{levi_civita, IndexError};
#[test]
fn epsilon_is_right_handed_and_alternating_in_all_components() {
    const POSITIVE: i8 = match levi_civita([0, 1, 2]) {
        Ok(s) => s,
        Err(_) => panic!("valid indices"),
    };
    assert_eq!(POSITIVE, 1);
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let expected = match (i, j, k) {
                    (0, 1, 2) | (1, 2, 0) | (2, 0, 1) => 1,
                    (0, 2, 1) | (2, 1, 0) | (1, 0, 2) => -1,
                    _ => 0,
                };
                assert_eq!(levi_civita([i, j, k]).unwrap(), expected);
            }
        }
    }
}
#[test]
fn epsilon_rejects_invalid_indices_instead_of_returning_zero() {
    assert_eq!(
        levi_civita([0, 1, 3]),
        Err(IndexError {
            index: 3,
            dimension: 3
        })
    );
    assert_eq!(
        levi_civita([99, 0, 0]),
        Err(IndexError {
            index: 99,
            dimension: 3
        })
    );
}
#[test]
fn epsilon_contraction_signs_and_dimensions() {
    for i in 0..3 {
        for j in 0..3 {
            for l in 0..3 {
                for m in 0..3 {
                    let lhs: i32 = (0..3)
                        .map(|k| {
                            i32::from(levi_civita([i, j, k]).unwrap())
                                * i32::from(levi_civita([k, l, m]).unwrap())
                        })
                        .sum();
                    let rhs = i32::from(i == l) * i32::from(j == m)
                        - i32::from(i == m) * i32::from(j == l);
                    assert_eq!(lhs, rhs);
                }
            }
        }
    }
    assert_eq!(levi_civita::<0>([]), Ok(1));
    assert_eq!(levi_civita([0]), Ok(1));
    assert_eq!(levi_civita([0, 1, 2, 3]), Ok(1));
    assert_eq!(levi_civita([3, 1, 2, 0]), Ok(-1));
}
