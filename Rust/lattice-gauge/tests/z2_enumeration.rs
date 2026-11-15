//! Oracles: the two-dimensional closed form from the character expansion, the total count
//! 2^links, and mono-notes Lemma ym:lem:twist-move (Z[S xor delta b] = Z[S]) checked on the
//! full exact density of states.
use lattice_gauge::cochain::{link_coboundary_set, symmetric_difference, thooft_surface};
use lattice_gauge::torus::Torus;
use lattice_gauge::z2_gauge::{
    binomial, satisfied_plaquette_histogram, two_dimensional_histogram_closed_form,
};

#[test]
fn two_dimensional_enumeration_matches_the_character_expansion() {
    for side in [2u32, 3] {
        let torus = Torus::new(&[side as usize, side as usize]).unwrap();
        let untwisted = vec![false; torus.plaquette_count()];
        let histogram = satisfied_plaquette_histogram(&torus, &untwisted).unwrap();
        assert_eq!(histogram, two_dimensional_histogram_closed_form(side, false).unwrap());
        assert_eq!(histogram.iter().sum::<u128>(), 1u128 << torus.link_count());
        let mut twisted = untwisted.clone();
        twisted[torus.plaquette_index(1, 0, 1)] = true;
        let twisted_histogram = satisfied_plaquette_histogram(&torus, &twisted).unwrap();
        assert_eq!(twisted_histogram, two_dimensional_histogram_closed_form(side, true).unwrap());
        assert_ne!(twisted_histogram, histogram);
        // Two twists cancel in two dimensions: the class depends only on the parity.
        twisted[torus.plaquette_index(0, 0, 1)] = true;
        assert_eq!(satisfied_plaquette_histogram(&torus, &twisted).unwrap(), histogram);
    }
}

#[test]
fn closed_form_values_for_the_two_by_two_torus() {
    // N_j = 2^4 binomial(4, j) (1 + (-1)^(4 - j)): 32, 0, 192, 0, 32. Sum 256 = 2^8.
    assert_eq!(two_dimensional_histogram_closed_form(2, false).unwrap(), vec![32, 0, 192, 0, 32]);
    assert_eq!(two_dimensional_histogram_closed_form(2, true).unwrap(), vec![0, 128, 0, 128, 0]);
    assert_eq!(binomial(9, 4), Some(126));
    assert_eq!(binomial(3, 5), Some(0));
}

#[test]
fn inputs_outside_the_enumerable_range_are_rejected() {
    let torus = Torus::new(&[4, 4]).unwrap();
    assert!(satisfied_plaquette_histogram(&torus, &vec![false; torus.plaquette_count()]).is_none());
    let torus = Torus::new(&[2, 2]).unwrap();
    assert!(satisfied_plaquette_histogram(&torus, &[false; 3]).is_none());
}

#[test]
fn moving_a_twist_by_a_link_coboundary_leaves_the_density_of_states_unchanged() {
    let torus = Torus::new(&[3, 3]).unwrap();
    let mut twist = vec![false; torus.plaquette_count()];
    twist[4] = true;
    let reference = satisfied_plaquette_histogram(&torus, &twist).unwrap();
    for link in [0usize, 7, 17] {
        let moved = symmetric_difference(&twist, &link_coboundary_set(&torus, link)).unwrap();
        assert_ne!(moved, twist);
        assert_eq!(satisfied_plaquette_histogram(&torus, &moved).unwrap(), reference);
    }
}

#[test]
fn three_dimensional_twist_on_the_two_cubed_torus() {
    let torus = Torus::new(&[2, 2, 2]).unwrap();
    assert_eq!(torus.link_count(), 24);
    let untwisted = satisfied_plaquette_histogram(&torus, &vec![false; torus.plaquette_count()]).unwrap();
    let surface = thooft_surface(&torus, 0, 1);
    let twisted = satisfied_plaquette_histogram(&torus, &surface).unwrap();
    let moved_surface = symmetric_difference(&surface, &link_coboundary_set(&torus, 11)).unwrap();
    let moved = satisfied_plaquette_histogram(&torus, &moved_surface).unwrap();
    assert_eq!(untwisted.iter().sum::<u128>(), 1 << 24);
    assert_eq!(twisted.iter().sum::<u128>(), 1 << 24);
    assert_eq!(moved, twisted, "Lemma ym:lem:twist-move");
    assert_ne!(twisted, untwisted, "the twist class is nontrivial");
    // Configurations with every plaquette satisfied are the closed Z2 1-cochains:
    // 2^((V - 1) + d) = 2^10 of them (gauge transformations times flat holonomies).
    assert_eq!(untwisted[torus.plaquette_count()], 1 << 10);
    // With the twist no configuration satisfies every plaquette.
    assert_eq!(twisted[torus.plaquette_count()], 0);
}
