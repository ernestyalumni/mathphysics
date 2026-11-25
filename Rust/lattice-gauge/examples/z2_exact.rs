//! Exact weight enumerator of Z2 lattice gauge theory on a d-dimensional torus (sides >= 1), optionally twisted.
//! Usage: z2_exact <sides, e.g. 2,2,2,2> <cpu|gpu> <threads> [twist planes, e.g. 01 23]
//! Prints the box, the rank r of the coboundary, the backend and the time, then `w A_w` for every w with A_w > 0:
//! the number of codewords (plaquette configurations reachable from link configurations, shifted by the twist)
//! with w plaquettes of value -1. Z(beta) = 2^(-r) sum_w A_w exp(-2 beta w) (energy sum_p (1 - sigma_p), normalized
//! counting measure). See `z2_code`.
use lattice_gauge::z2_code::{coset_enumerator, span_basis, Box};
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 4 {
        eprintln!("usage: z2_exact <sides> <cpu|gpu> <threads> [twist planes]");
        std::process::exit(2);
    }
    let sides: Vec<usize> = a[1].split(',').map(|s| s.parse().unwrap()).collect();
    let threads: usize = a[3].parse().unwrap();
    let b = Box::new(&sides);
    let n = b.plaquette_count();
    let mut basis = span_basis(&b.coboundary_columns());
    // Benchmarking only: Z2_RANK_LIMIT=k enumerates the subcode spanned by the first k basis vectors.
    if let Ok(k) = std::env::var("Z2_RANK_LIMIT") { basis.truncate(k.parse().unwrap()); }
    let mut twist = 0u128;
    for t in &a[4..] {
        let c: Vec<usize> = t.chars().map(|c| c.to_digit(10).unwrap() as usize).collect();
        twist ^= b.twist_mask(c[0], c[1]);
    }
    let t0 = Instant::now();
    let hist = match a[2].as_str() {
        "cpu" => coset_enumerator(&basis, twist, n, threads),
        #[cfg(feature = "cuda")]
        "gpu" => lattice_gauge::z2_code::coset_enumerator_gpu(&basis, twist, n).unwrap(),
        other => panic!("backend {other} not available in this build"),
    };
    let secs = t0.elapsed().as_secs_f64();
    let total: u64 = hist.iter().sum();
    assert_eq!(total, 1u64 << basis.len(), "histogram does not sum to 2^r");
    println!("# sides {} links {} plaquettes {n} rank {} twist [{}] backend {} threads {threads} seconds {secs:.3} codewords/s {:.3e}",
        a[1], b.links, basis.len(), a[4..].join(" "), a[2], total as f64 / secs.max(1e-9));
    for (w, c) in hist.iter().enumerate() {
        if *c > 0 { println!("{w}\t{c}"); }
    }
}
