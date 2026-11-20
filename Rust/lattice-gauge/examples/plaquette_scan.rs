//! Plaquette average on one torus over a grid of beta values, annealed upward from beta = 0, for
//! thermodynamic integration: log Z(beta) = N_P * int_0^beta <a0> dbeta' (Haar measure of mass 1).
//! Usage: plaquette_scan <dims, e.g. 4,4,4,8> <beta_max> <intervals> <therm> <sweeps> <seed> <out.tsv>
//! Output rows: beta, mean <a0>, jackknife error, tau_int, plaquette count.
use lattice_gauge::stats::{jackknife, tau_int};
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
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[7]).unwrap());
    let mut l = Lattice::hot(&dims, 0.0, seed, 0).unwrap();
    let np = l.torus.plaquette_count();
    for i in 0..=n {
        l.beta = bmax * i as f64 / n as f64;
        for _ in 0..therm { l.update(2); }
        let mut p = Vec::with_capacity(sweeps);
        for k in 0..sweeps {
            l.update(2);
            if k % 1000 == 999 { l.reunitarize(); }
            p.push(l.plaquette());
        }
        let (m, e) = jackknife(&[&p], 50, |v| v[0]);
        writeln!(f, "{}\t{:.12e}\t{:.6e}\t{:.3}\t{}", l.beta, m, e, tau_int(&p, 6.0), np).unwrap();
        f.flush().unwrap();
    }
}
