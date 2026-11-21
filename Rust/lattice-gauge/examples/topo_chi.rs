//! Topological susceptibility from cooled clover charges (Teper hep-th/9812187 Table 15: 15-20 cooling sweeps).
//! Usage: topo_chi <side> <beta> <therm> <configs> <every> <seed> <threads> <out.tsv>
//! Writes per configuration: plaquette, Q after 15 and 20 cooling sweeps. Prints a chi^{1/4} (raw Q_20, and Q_20
//! rescaled by the factor that best maps the charges onto integers, then rounded).
use lattice_gauge::stats::jackknife;
use lattice_gauge::topology::{charge, cool};
use lattice_gauge::wilson::Lattice;
use std::io::Write;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let (therm, configs, every): (usize, usize, usize) = (a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap());
    let seed: u64 = a[6].parse().unwrap();
    let threads: usize = a[7].parse().unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[8]).unwrap());
    let mut l = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    l.threads = threads;
    for i in 0..therm { l.update(2); if i % 100 == 99 { l.reunitarize(); } }
    let mut q20 = Vec::with_capacity(configs);
    for c in 0..configs {
        for _ in 0..every { l.update(2); }
        if c % 20 == 19 { l.reunitarize(); }
        let mut k = Lattice::cold(&[side; 4], beta, seed, 1).unwrap();
        k.threads = threads;
        k.links = l.links.clone();
        cool(&mut k, 15);
        let q15 = charge(&k);
        cool(&mut k, 5);
        let q = charge(&k);
        writeln!(f, "{:.6}\t{:.6}\t{:.6}", l.plaquette(), q15, q).unwrap();
        f.flush().unwrap();
        q20.push(q);
    }
    let v = (side as f64).powi(4);
    let q2: Vec<f64> = q20.iter().map(|q| q * q).collect();
    let (chi_raw, e_raw) = jackknife(&[&q2], 20, |m| (m[0] / v).powf(0.25));
    // rescale factor alpha in [1, 1.3] minimizing the distance of alpha Q to integers (for |Q| > 0.5)
    let alpha = (0..=300).map(|i| 1.0 + i as f64 * 0.001).min_by(|x, y| {
        let d = |al: f64| q20.iter().filter(|q| q.abs() > 0.5).map(|q| (al * q - (al * q).round()).powi(2)).sum::<f64>();
        d(*x).partial_cmp(&d(*y)).unwrap()
    }).unwrap();
    let qi2: Vec<f64> = q20.iter().map(|q| (alpha * q).round().powi(2)).collect();
    let (chi_int, e_int) = jackknife(&[&qi2], 20, |m| (m[0] / v).powf(0.25));
    println!("topo side={side} beta={beta} configs={configs} every={every}: a chi^1/4 raw Q20 = {chi_raw:.4}({:.0}), integer-rounded (alpha {alpha:.3}) = {chi_int:.4}({:.0})", e_raw * 1e4, e_int * 1e4);
}
