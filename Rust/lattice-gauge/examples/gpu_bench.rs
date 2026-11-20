//! GPU timing and a statistical check against the CPU path. Usage:
//!   gpu_bench <side> <beta> <therm> <meas> <every> <seed>
//! Plaquettes are measured on the CPU from downloaded configurations.
use lattice_gauge::cuda::GpuLattice;
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;
use std::time::Instant;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let therm: usize = a[3].parse().unwrap();
    let meas: usize = a[4].parse().unwrap();
    let every: usize = a[5].parse().unwrap();
    let seed: u64 = a[6].parse().unwrap();
    let mut cpu = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    let mut gpu = GpuLattice::from_cpu(&cpu).unwrap();
    let t0 = Instant::now();
    for _ in 0..therm { gpu.update(2).unwrap(); }
    let therm_s = t0.elapsed().as_secs_f64();
    let mut p = Vec::with_capacity(meas);
    for _ in 0..meas {
        for _ in 0..every { gpu.update(2).unwrap(); }
        cpu.links = gpu.download().unwrap();
        p.push(cpu.plaquette());
    }
    let ms_per_sweep = therm_s * 1e3 / (therm as f64 * 3.0);
    let (m, e) = if meas >= 20 { jackknife(&[&p], 20, |v| v[0]) } else { (f64::NAN, f64::NAN) };
    println!("gpu side={side} beta={beta} therm={therm} meas={meas} every={every} plaquette={m:.6} err={e:.6} ms_per_sweep={ms_per_sweep:.3} link_updates_per_s={:.3e}",
        (cpu.links.len() as f64) / (ms_per_sweep * 1e-3));
}
