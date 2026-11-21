//! Plaquette average on one torus over a grid of beta values, annealed upward from beta = 0, for
//! thermodynamic integration: log Z(beta) = N_P * int_0^beta <a0> dbeta' (Haar measure of mass 1).
//! Usage: plaquette_scan <dims, e.g. 4,4,4,8> <beta_max> <intervals> <therm> <sweeps> <seed> <out.tsv> [beta_min, default 0] [group: su2 (default) or u1]
//! Output rows: beta, mean <a0>, jackknife error, tau_int, plaquette count.
use lattice_gauge::stats::{jackknife, tau_int};
use lattice_gauge::u1::U1Lattice;
use lattice_gauge::wilson::Lattice;
use std::io::Write;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let dims: Vec<usize> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let bmax: f64 = a[2].parse().unwrap();
    let n: usize = a[3].parse().unwrap();
    let therm: usize = a[4].parse().unwrap();
    let sweeps: usize = a[5].parse().unwrap();
    let seed: u64 = a[6].parse().unwrap();
    let bmin: f64 = a.get(8).map(|x| x.parse().unwrap()).unwrap_or(0.0);
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[7]).unwrap());
    let u1 = a.get(9).map(|g| g == "u1").unwrap_or(false);
    let mut l = Lattice::hot(&dims, 0.0, seed, 0).unwrap();
    let mut lu = U1Lattice::hot(&dims, 0.0, seed, 0).unwrap();
    let np = l.torus.plaquette_count();
    for i in 0..=n {
        let beta = bmin + (bmax - bmin) * i as f64 / n as f64;
        l.beta = beta;
        lu.beta = beta;
        let mut step = |k: usize| -> f64 {
            if u1 { lu.update(2); lu.plaquette() } else { l.update(2); if k % 1000 == 999 { l.reunitarize(); } l.plaquette() }
        };
        for k in 0..therm { step(k); }
        let p: Vec<f64> = (0..sweeps).map(|k| step(k)).collect();
        let (m, e) = jackknife(&[&p], 50, |v| v[0]);
        writeln!(f, "{}\t{:.12e}\t{:.6e}\t{:.3}\t{}", beta, m, e, tau_int(&p, 6.0), np).unwrap();
        f.flush().unwrap();
    }
}
