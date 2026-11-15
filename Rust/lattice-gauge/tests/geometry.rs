//! Oracles: cell counts d V, d(d-1)/2 V, binomial(d,3) V; double counting of link-plaquette
//! incidences; delta delta = 0; and the first Betti number d of the d-torus, through
//! rank(delta: links -> plaquettes) = d V - V + 1 - d over GF(2).
use lattice_gauge::cochain::{link_coboundary_rank, link_coboundary_set, odd_cube_count};
use lattice_gauge::torus::Torus;
use std::collections::BTreeSet;

fn tori() -> Vec<Torus> {
    [vec![2, 2, 2, 2], vec![3, 4, 2, 5], vec![3, 3, 3], vec![2, 3], vec![4, 4, 2, 2]]
        .iter()
        .map(|sides| Torus::new(sides).unwrap())
        .collect()
}

#[test]
fn cell_counts_match_the_closed_formulas() {
    let torus = Torus::new(&[2, 2, 2, 2]).unwrap();
    assert_eq!(
        (torus.site_count(), torus.link_count(), torus.plaquette_count(), torus.cube_count()),
        (16, 64, 96, 64)
    );
    let torus = Torus::new(&[3, 4, 2, 5]).unwrap();
    assert_eq!((torus.site_count(), torus.link_count(), torus.plaquette_count()), (120, 480, 720));
    assert_eq!(Torus::new(&[2, 3]).unwrap().cube_count(), 0);
    assert_eq!(torus.cubes().len(), torus.cube_count());
}

#[test]
fn degenerate_tori_are_rejected() {
    assert!(Torus::new(&[]).is_none());
    assert!(Torus::new(&[1, 3]).is_none());
    assert!(Torus::new(&[3, 0]).is_none());
    assert!(Torus::new(&[usize::MAX, 4]).is_none());
}

#[test]
fn indices_round_trip_and_translations_are_periodic() {
    for torus in tori() {
        let d = torus.dimension();
        for site in 0..torus.site_count() {
            assert_eq!(torus.site_index(&torus.site_coordinates(site)), site);
            for direction in 0..d {
                let forward = torus.neighbor(site, direction, true);
                assert_ne!(forward, site);
                assert_eq!(torus.neighbor(forward, direction, false), site);
                let mut walker = site;
                for _ in 0..torus.side_lengths()[direction] {
                    walker = torus.neighbor(walker, direction, true);
                }
                assert_eq!(walker, site);
            }
        }
        for plane in 0..torus.plane_count() {
            let (first, second) = torus.plane_directions(plane);
            assert!(first < second && second < d);
            assert_eq!(torus.plane_index(first, second), plane);
        }
    }
}

#[test]
fn link_plaquette_incidence_is_consistent_in_both_directions() {
    for torus in tori() {
        let d = torus.dimension();
        let mut incidences_from_plaquettes = 0;
        for site in 0..torus.site_count() {
            for first in 0..d {
                for second in (first + 1)..d {
                    let links = torus.plaquette_links(site, first, second);
                    assert_eq!(links.iter().collect::<BTreeSet<_>>().len(), 4);
                    let plaquette = torus.plaquette_index(site, first, second);
                    for link in links {
                        let coboundary = torus.link_coboundary(link);
                        assert_eq!(coboundary.iter().filter(|&&p| p == plaquette).count(), 1);
                        incidences_from_plaquettes += 1;
                    }
                }
            }
        }
        let mut incidences_from_links = 0;
        for link in 0..torus.link_count() {
            let coboundary = torus.link_coboundary(link);
            assert_eq!(coboundary.len(), 2 * (d - 1));
            assert_eq!(coboundary.iter().collect::<BTreeSet<_>>().len(), 2 * (d - 1));
            incidences_from_links += coboundary.len();
        }
        assert_eq!(incidences_from_plaquettes, incidences_from_links);
        assert_eq!(incidences_from_links, 4 * torus.plaquette_count());
    }
}

#[test]
fn the_coboundary_of_a_coboundary_vanishes() {
    for torus in tori() {
        for faces in torus.cubes() {
            assert_eq!(faces.iter().collect::<BTreeSet<_>>().len(), 6);
        }
        for link in 0..torus.link_count() {
            assert_eq!(odd_cube_count(&torus, &link_coboundary_set(&torus, link)), 0);
        }
    }
}

#[test]
fn coboundary_rank_reproduces_the_first_betti_number_of_the_torus() {
    for torus in tori() {
        let (d, v) = (torus.dimension(), torus.site_count());
        assert_eq!(link_coboundary_rank(&torus), d * v - v + 1 - d, "{:?}", torus.side_lengths());
    }
    assert_eq!(link_coboundary_rank(&Torus::new(&[2, 2, 2, 2]).unwrap()), 45);
    assert_ne!(link_coboundary_rank(&Torus::new(&[2, 2, 2, 2]).unwrap()), 64 - 16 + 1);
}
