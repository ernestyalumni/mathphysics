//! Exact defect, purities and twist free energies of Z2 gauge theory at n = 1, from weight enumerators.
//! Usage: z2_defect <dir with z2_exact outputs for the 2^4 torus: none.txt and optional twist-class files>
//! Boxes with sides in {1, 2} other than 2^4 are enumerated here. Convention: energy sum_p (1 - sigma_p), normalized
//! counting measure, so Z = 2^(-r) sum_w A_w exp(-2 beta w) and Z(0) = 1 (`z2_code`). Prints, against beta:
//!   Delta(1) = 16 log Z(1,1,1,1) - log Z(2,2,2,2), the purities P1..P4 of the telescope (doubling axes 0, 1, 2, 3
//!   in turn), the check log(P1^8 P2^4 P3^2 P4) = -Delta, and F = -log(Z_twist / Z) for each twist class found;
//! then the beta at which Delta(1) = 1/7480 (bisection) against the elementary bound Delta <= 192 beta.
use lattice_gauge::z2_code::{coset_enumerator, span_basis, Box};

struct Enumerator { rank: usize, a: Vec<(usize, u64)> }

impl Enumerator {
    fn of_box(sides: &[usize]) -> Self {
        let b = Box::new(sides);
        let basis = span_basis(&b.coboundary_columns());
        let h = coset_enumerator(&basis, 0, b.plaquette_count(), 8);
        Enumerator { rank: basis.len(), a: h.iter().enumerate().filter(|(_, c)| **c > 0).map(|(w, c)| (w, *c)).collect() }
    }
    fn read(path: &str) -> Option<Self> {
        let text = std::fs::read_to_string(path).ok()?;
        let header = text.lines().next()?;
        let words: Vec<&str> = header.split_whitespace().collect();
        let rank: usize = words[words.iter().position(|w| *w == "rank")? + 1].parse().ok()?;
        let a = text.lines().skip(1).map(|l| { let mut it = l.split('\t'); (it.next().unwrap().parse().unwrap(), it.next().unwrap().parse().unwrap()) }).collect();
        Some(Enumerator { rank, a })
    }
    /// log Z for the untwisted code (Z(0) = 1), accurate at small beta via log1p/expm1.
    fn log_z(&self, beta: f64) -> f64 {
        let norm = (2f64).powi(self.rank as i32);
        if beta < 0.05 {
            let s: f64 = self.a.iter().map(|&(w, c)| c as f64 * (-2.0 * beta * w as f64).exp_m1()).sum::<f64>() / norm;
            s.ln_1p()
        } else {
            self.log_sum(beta) - norm.ln()
        }
    }
    /// log sum_w A_w exp(-2 beta w) (no normalization), by log-sum-exp.
    fn log_sum(&self, beta: f64) -> f64 {
        let terms: Vec<f64> = self.a.iter().map(|&(w, c)| (c as f64).ln() - 2.0 * beta * w as f64).collect();
        let m = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        m + terms.iter().map(|t| (t - m).exp()).sum::<f64>().ln()
    }
}

fn main() {
    let dir = std::env::args().nth(1).expect("usage: z2_defect <dir>");
    let z1111 = Enumerator::of_box(&[1, 1, 1, 1]);
    let z2111 = Enumerator::of_box(&[2, 1, 1, 1]);
    let z2211 = Enumerator::of_box(&[2, 2, 1, 1]);
    let z2221 = Enumerator::of_box(&[2, 2, 2, 1]);
    let z2222 = Enumerator::read(&format!("{dir}/none.txt")).expect("none.txt");
    let mut twists: Vec<(String, Enumerator)> = Vec::new();
    for name in ["e01", "e01-02", "e01-23", "tri", "star", "path", "c4", "paw", "k4-e", "k4"] {
        if let Some(e) = Enumerator::read(&format!("{dir}/{name}.txt")) { twists.push((name.into(), e)); }
    }
    println!("# ranks: 1^4 {} 2x1^3 {} 2^2x1^2 {} 2^3x1 {} 2^4 {}", z1111.rank, z2111.rank, z2211.rank, z2221.rank, z2222.rank);
    let delta = |b: f64| 16.0 * z1111.log_z(b) - z2222.log_z(b);
    print!("beta\tDelta(1)\tDelta/beta\tlogP1\tlogP2\tlogP3\tlogP4\tcheck");
    for (n, _) in &twists { print!("\tF[{n}]"); }
    println!();
    for &b in &[1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0] {
        let (l1, l2, l3, l4, l5) = (z1111.log_z(b), z2111.log_z(b), z2211.log_z(b), z2221.log_z(b), z2222.log_z(b));
        let (p1, p2, p3, p4) = (l2 - 2.0 * l1, l3 - 2.0 * l2, l4 - 2.0 * l3, l5 - 2.0 * l4);
        let d = delta(b);
        let check = 8.0 * p1 + 4.0 * p2 + 2.0 * p3 + p4 + d;
        print!("{b:.1e}\t{d:.6e}\t{:.6}\t{p1:.4e}\t{p2:.4e}\t{p3:.4e}\t{p4:.4e}\t{check:.1e}", d / b);
        for (_, e) in &twists { print!("\t{:.6e}", -(e.log_sum(b) - z2222.log_sum(b))); }
        println!();
    }
    let target = 1.0 / 7480.0;
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..200 { let mid = 0.5 * (lo + hi); if delta(mid) > target { hi = mid } else { lo = mid } }
    println!("# Delta(1) = 1/7480 at beta = {:.10e}; elementary bound 192 beta gives beta = {:.10e}; ratio {:.4}", lo, target / 192.0, lo / (target / 192.0));
    // leading small-beta coefficient: Delta(1) = c beta + O(beta^2), c = 2 <w> on 2^4 (log Z(1^4) = 0 if rank 0)
    let mean_w: f64 = z2222.a.iter().map(|&(w, c)| w as f64 * c as f64).sum::<f64>() / (2f64).powi(z2222.rank as i32);
    println!("# mean frustrated plaquettes on 2^4 at beta = 0: {mean_w}; so Delta(1) = {:.4} beta + O(beta^2) when Z(1^4) = 1", 2.0 * mean_w);
}
