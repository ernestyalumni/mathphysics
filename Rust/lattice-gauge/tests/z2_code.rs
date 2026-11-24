use lattice_gauge::torus::Torus;
use lattice_gauge::z2_code::{brute_force, coset_enumerator, dual_basis, macwilliams, span_basis, Box};
use lattice_gauge::z2_gauge::satisfied_plaquette_histogram;

fn code_enumerator(sides: &[usize], twist: u128) -> (Vec<u64>, usize, usize) {
    let b = Box::new(sides);
    let basis = span_basis(&b.coboundary_columns());
    (coset_enumerator(&basis, twist, b.plaquette_count(), 4), basis.len(), b.links)
}

#[test]
fn matches_independent_brute_force_on_even_boxes() {
    for sides in [vec![2, 2], vec![2, 4], vec![3, 3], vec![2, 6], vec![2, 2, 2]] {
        let b = Box::new(&sides);
        let torus = Torus::new(&sides).unwrap();
        let p = b.plaquette_count();
        assert_eq!(p, torus.plaquette_count());
        let reference = satisfied_plaquette_histogram(&torus, &vec![false; p]).unwrap();
        let (a, r, links) = code_enumerator(&sides, 0);
        for w in 0..=p {
            assert_eq!(a[w] as u128 * (1u128 << (links - r)), reference[p - w], "sides {sides:?} w {w}");
        }
    }
}

#[test]
fn code_walk_matches_link_walk_including_side_one() {
    for sides in [vec![1, 1, 1, 1], vec![2, 1, 1, 1], vec![2, 2, 1, 1], vec![1, 2, 1, 2], vec![2, 1, 2], vec![3, 1, 2]] {
        let b = Box::new(&sides);
        for twist in [0u128, b.twist_mask(0, 1)] {
            let (a, r, links) = code_enumerator(&sides, twist);
            let direct = brute_force(&b, twist, 4);
            for w in 0..a.len() {
                assert_eq!(a[w] * (1u64 << (links - r)), direct[w], "sides {sides:?} twist {twist} w {w}");
            }
        }
    }
}

#[test]
fn all_one_torus_has_only_trivial_plaquettes() {
    let (a, r, _) = code_enumerator(&[1, 1, 1, 1], 0);
    assert_eq!(r, 0);
    assert_eq!(a[0], 1);
}

#[test]
fn macwilliams_identity_holds() {
    for sides in [vec![2, 2, 1, 1], vec![2, 2, 2], vec![2, 2, 2, 1]] {
        let b = Box::new(&sides);
        let n = b.plaquette_count();
        let basis = span_basis(&b.coboundary_columns());
        let a = coset_enumerator(&basis, 0, n, 4);
        let dual = dual_basis(&basis, n);
        assert_eq!(basis.len() + dual.len(), n);
        for &y in &dual {
            for &c in &basis { assert_eq!((y & c).count_ones() % 2, 0); }
        }
        let bw = coset_enumerator(&dual, 0, n, 4);
        let predicted = macwilliams(&a, n).unwrap();
        for v in 0..=n { assert_eq!(predicted[v], bw[v] as i128, "sides {sides:?} v {v}"); }
    }
}
