//! Oracle: the exact SU(2) plaquette average on a 2D periodic torus from the character expansion
//! (character::torus_plaquette_2d), derived independently of the Monte Carlo code.
use lattice_gauge::character::torus_plaquette_2d;
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::{Lattice, Update};

fn run(side: usize, beta: f64, seed: u64, update: Update, therm: usize, meas: usize) -> (f64, f64) {
    let mut l = Lattice::hot(&[side, side], beta, seed, 0).unwrap();
    for _ in 0..therm {
        l.sweep_with(update);
        l.sweep_with(Update::OverRelax);
    }
    let mut p = Vec::with_capacity(meas);
    for _ in 0..meas {
        l.sweep_with(update);
        l.sweep_with(Update::OverRelax);
        p.push(l.plaquette());
    }
    jackknife(&[&p], 20, |m| m[0])
}

#[test]
fn plaquette_matches_exact_torus() {
    for &(side, beta) in &[(4usize, 1.0), (8, 2.0), (8, 4.0), (6, 0.5)] {
        let exact = torus_plaquette_2d(beta, (side * side) as u64);
        let (m, e) = run(side, beta, 11 + side as u64, Update::HeatBath, 200, 4000);
        assert!((m - exact).abs() < 5.0 * e, "L={side} beta={beta}: {m} +- {e} vs exact {exact}");
    }
}

#[test]
fn broken_update_fails() {
    let (side, beta) = (8usize, 2.0);
    let exact = torus_plaquette_2d(beta, 64);
    let (m, e) = run(side, beta, 21, Update::HeatBathBroken, 200, 4000);
    assert!((m - exact).abs() > 8.0 * e, "negative control passed: {m} +- {e} vs {exact}");
}

#[test]
fn independent_of_thread_count() {
    let mut a = Lattice::hot(&[6, 6], 1.5, 7, 0).unwrap();
    let mut b = Lattice::hot(&[6, 6], 1.5, 7, 0).unwrap();
    a.threads = 1;
    b.threads = 5;
    for _ in 0..20 {
        a.update(2);
        b.update(2);
    }
    assert!(a.links == b.links, "configurations differ between 1 and 5 threads");
}

#[test]
fn odd_sides_rejected() {
    assert!(Lattice::cold(&[5, 6], 1.0, 0, 0).is_err());
}
