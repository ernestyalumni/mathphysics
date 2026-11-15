//! Machine check of mono-notes Lemma ym:lem:thooft-class on four-dimensional tori, by two
//! independent methods: the parity argument of the proof, and Gaussian elimination over GF(2).
use lattice_gauge::cochain::{
    is_sum_of_link_coboundaries, link_coboundary_set, odd_cube_count,
    plane_slice_intersection_count, symmetric_difference, thooft_surface,
};
use lattice_gauge::torus::Torus;

fn four_dimensional_tori() -> Vec<Torus> {
    [[2, 2, 2, 2], [3, 2, 2, 3], [3, 3, 3, 3], [4, 2, 3, 2]]
        .iter()
        .map(|sides| Torus::new(sides).unwrap())
        .collect()
}

#[test]
fn the_thooft_surface_creates_no_monopole() {
    for torus in four_dimensional_tori() {
        for first in 0..4 {
            for second in (first + 1)..4 {
                assert_eq!(odd_cube_count(&torus, &thooft_surface(&torus, first, second)), 0);
            }
        }
    }
}

#[test]
fn negative_control_a_single_twisted_plaquette_creates_monopoles() {
    for torus in four_dimensional_tori() {
        let mut set = vec![false; torus.plaquette_count()];
        set[torus.plaquette_index(0, 0, 1)] = true;
        // The plaquette is a face of two cubes in each of the two transverse directions.
        assert_eq!(odd_cube_count(&torus, &set), 4);
        assert_eq!(is_sum_of_link_coboundaries(&torus, &set), Some(false));
    }
}

#[test]
fn parity_argument_of_the_lemma() {
    for torus in four_dimensional_tori() {
        let surface = thooft_surface(&torus, 0, 1);
        for anchor in 0..torus.site_count() {
            assert_eq!(plane_slice_intersection_count(&torus, &surface, 0, 1, anchor), 1);
        }
        for link in 0..torus.link_count() {
            let coboundary = link_coboundary_set(&torus, link);
            for anchor in 0..torus.site_count() {
                let count = plane_slice_intersection_count(&torus, &coboundary, 0, 1, anchor);
                assert!(count == 0 || count == 2, "count {count}");
            }
        }
    }
}

#[test]
fn the_thooft_surface_is_not_a_sum_of_link_coboundaries() {
    for torus in four_dimensional_tori() {
        let empty = vec![false; torus.plaquette_count()];
        assert_eq!(is_sum_of_link_coboundaries(&torus, &empty), Some(true));
        let surface = thooft_surface(&torus, 0, 1);
        assert_eq!(is_sum_of_link_coboundaries(&torus, &surface), Some(false));
        let moved = symmetric_difference(&surface, &link_coboundary_set(&torus, 5)).unwrap();
        assert_eq!(is_sum_of_link_coboundaries(&torus, &moved), Some(false));
        let two_links = symmetric_difference(
            &link_coboundary_set(&torus, 3),
            &link_coboundary_set(&torus, torus.link_count() - 1),
        )
        .unwrap();
        assert_eq!(is_sum_of_link_coboundaries(&torus, &two_links), Some(true));
        // Twists in different planes are different classes.
        let other = thooft_surface(&torus, 2, 3);
        let combined = symmetric_difference(&surface, &other).unwrap();
        assert_eq!(is_sum_of_link_coboundaries(&torus, &combined), Some(false));
        assert_eq!(is_sum_of_link_coboundaries(&torus, &surface[1..]), None);
    }
}

#[test]
fn translating_the_twist_surface_changes_it_by_link_coboundaries() {
    for torus in four_dimensional_tori() {
        let surface = thooft_surface(&torus, 0, 1);
        // The same surface through the sites with coordinates (1, 0) in the directions 0, 1.
        let mut translated = vec![false; torus.plaquette_count()];
        for site in 0..torus.site_count() {
            if torus.site_coordinate(site, 0) == 1 && torus.site_coordinate(site, 1) == 0 {
                translated[torus.plaquette_index(site, 0, 1)] = true;
            }
        }
        let difference = symmetric_difference(&surface, &translated).unwrap();
        assert_eq!(is_sum_of_link_coboundaries(&torus, &difference), Some(true));
    }
}
