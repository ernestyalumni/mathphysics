//! Quick look: Q after 0, 5, 10, 20, 40 cooling sweeps on thermalized SU(2) configurations.
//! Usage: topo_check <side> <beta> <configs> <seed>
use lattice_gauge::topology::{charge, cool};
use lattice_gauge::wilson::Lattice;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let n: usize = a[3].parse().unwrap();
    let seed: u64 = a[4].parse().unwrap();
    let mut l = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    for _ in 0..300 { l.update(2); }
    for c in 0..n {
        for _ in 0..20 { l.update(2); }
        let mut k = Lattice::cold(&[side; 4], beta, seed, 1).unwrap();
        k.links = l.links.clone();
        let mut row = vec![format!("{:.3}", charge(&k))];
        let mut done = 0;
        for target in [5, 10, 20, 40] { cool(&mut k, target - done); done = target; row.push(format!("{:.3}", charge(&k))); }
        println!("config {c}: Q after 0/5/10/20/40 cooling sweeps = {}  (plaquette after cooling {:.4})", row.join(" "), k.plaquette());
    }
}
