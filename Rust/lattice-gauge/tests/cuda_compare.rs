//! GPU against CPU, configuration by configuration (roadmap T13 acceptance): same start, same
//! streams. Over-relaxation is bitwise identical; the heat bath differs only in the last bits of
//! device log/cos, so links agree to a small drift and no accept decision may differ.
#![cfg(feature = "cuda")]
use lattice_gauge::cuda::GpuLattice;
use lattice_gauge::su2::Quaternion;
use lattice_gauge::wilson::{Lattice, Update};

fn differing(a: &[Quaternion], b: &[Quaternion], tol: f64) -> (usize, f64) {
    let mut n = 0;
    let mut max = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        let d = x.a.iter().zip(y.a.iter()).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        max = max.max(d);
        if d > tol {
            n += 1;
        }
    }
    (n, max)
}

#[test]
fn gpu_matches_cpu_link_by_link() {
    let mut cpu = Lattice::hot(&[6, 4, 4, 8], 2.3, 42, 1).unwrap();
    let mut gpu = GpuLattice::from_cpu(&cpu).unwrap();
    for step in 0..3 {
        cpu.sweep_with(Update::HeatBath);
        gpu.sweep_with(Update::HeatBath).unwrap();
        cpu.sweep_with(Update::OverRelax);
        gpu.sweep_with(Update::OverRelax).unwrap();
        let g = gpu.download().unwrap();
        // Last-bit differences of device log/cos grow chaotically from ~1e-14 per update; a flipped
        // accept decision would show as an O(1) difference. Require no flip and a small drift.
        let (flips, max) = differing(&cpu.links, &g, 1e-8);
        eprintln!("after update {}: {} of {} links differ by more than 1e-8; max difference {:.2e}", step + 1, flips, g.len(), max);
        assert_eq!(flips, 0, "a heat-bath accept decision differs between CPU and GPU");
        assert!(max < 1e-9, "drift {max:.2e} larger than expected after {} updates", step + 1);
    }
}

#[test]
fn overrelaxation_is_bitwise() {
    // Over-relaxation uses only products, sums and one sqrt (correctly rounded on both sides).
    let mut cpu = Lattice::hot(&[4, 4, 4, 4], 2.5, 7, 0).unwrap();
    let mut gpu = GpuLattice::from_cpu(&cpu).unwrap();
    for _ in 0..4 {
        cpu.sweep_with(Update::OverRelax);
        gpu.sweep_with(Update::OverRelax).unwrap();
    }
    assert!(gpu.download().unwrap() == cpu.links, "over-relaxation differs bitwise");
}

#[test]
fn device_plaquette_matches_cpu() {
    let mut cpu = Lattice::hot(&[4, 6, 4, 4], 2.2, 3, 0).unwrap();
    for _ in 0..3 { cpu.update(1); }
    let gpu = GpuLattice::from_cpu(&cpu).unwrap();
    assert!((gpu.plaquette().unwrap() - cpu.plaquette()).abs() < 1e-13);
}
