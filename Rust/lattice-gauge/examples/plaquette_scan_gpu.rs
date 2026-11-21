//! GPU version of plaquette_scan (same output format): updates and plaquette sums on the device.
//! Usage: plaquette_scan_gpu <dims> <beta_max> <intervals> <therm> <sweeps> <seed> <out.tsv> [beta_min, default 0]
use lattice_gauge::cuda::GpuLattice;
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
    let bmin: f64 = a.get(8).map(|x| x.parse().unwrap()).unwrap_or(0.0);
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[7]).unwrap());
    let mut cpu = Lattice::hot(&dims, 0.0, seed, 0).unwrap();
    let np = cpu.torus.plaquette_count();
    let mut g = GpuLattice::from_cpu(&cpu).unwrap();
    for i in 0..=n {
        g.beta = bmin + (bmax - bmin) * i as f64 / n as f64;
        for _ in 0..therm { g.update(2).unwrap(); }
        let mut p = Vec::with_capacity(sweeps);
        for k in 0..sweeps {
            g.update(2).unwrap();
            if k % 2000 == 1999 { cpu.links = g.download().unwrap(); cpu.reunitarize(); g.upload(&cpu.links).unwrap(); }
            p.push(g.plaquette().unwrap());
        }
        let (m, e) = jackknife(&[&p], 50, |v| v[0]);
        writeln!(f, "{}\t{:.12e}\t{:.6e}\t{:.3}\t{}", g.beta, m, e, tau_int(&p, 6.0), np).unwrap();
        f.flush().unwrap();
    }
}
