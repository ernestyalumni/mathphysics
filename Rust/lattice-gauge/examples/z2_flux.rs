//! Exact 't Hooft electric and magnetic flux free energies of Z2 gauge theory on the 2^4 torus.
//! Usage: z2_flux <dir with z2_exact outputs for the 11 twist classes>
//!
//! A twist is a vector n in F2^6 indexed by the planes (01, 02, 03, 12, 13, 23); its class is the isomorphism class of the
//! graph on 4 vertices with those edges, and the 11 class enumerators give Z[n] for all 64 twists. Take axis 0 as time:
//! k = (n01, n02, n03) are the temporal twists and m = (n23, n13, n12) the spatial (magnetic) ones. With z[n] = Z[n]/Z[0],
//!   electric flux e in F2^3 at magnetic flux m:  Z_e(e; m) = 2^-3 sum_k (-1)^(e.k) Z[k, m]   (Fourier transform over k),
//!   F_e(e) = -log( Z_e(e; 0) / Z_e(0; 0) ),     F_m(m) = -log( Z[0, m] / Z[0, 0] ).
//! Confinement in 't Hooft's sense: F_e grows without bound with the box while F_m -> 0; the deconfined (Z2 topologically
//! ordered) phase is the reverse. On a 2^4 box this is a finite-size control, not a phase diagnosis.

fn read(path: &str) -> Vec<(usize, u64)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines().skip(1).map(|l| { let mut it = l.split('\t'); (it.next().unwrap().parse().unwrap(), it.next().unwrap().parse().unwrap()) }).collect()
}

const EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

fn canonical(mask: u8) -> u8 {
    let mut best = u8::MAX;
    let mut p = [0usize, 1, 2, 3];
    loop {
        let mut o = 0u8;
        for (i, &(a, b)) in EDGES.iter().enumerate() {
            if mask >> i & 1 == 1 {
                let e = (p[a].min(p[b]), p[a].max(p[b]));
                o |= 1 << EDGES.iter().position(|&x| x == e).unwrap();
            }
        }
        best = best.min(o);
        // next permutation
        let mut i = 3;
        while i > 0 && p[i - 1] >= p[i] { i -= 1; }
        if i == 0 { break; }
        let mut j = 3;
        while p[j] <= p[i - 1] { j -= 1; }
        p.swap(i - 1, j);
        p[i..].reverse();
    }
    best
}

fn mask_of(planes: &[(usize, usize)]) -> u8 {
    planes.iter().fold(0u8, |m, e| m | 1 << EDGES.iter().position(|x| x == e).unwrap())
}

fn main() {
    let dir = std::env::args().nth(1).expect("usage: z2_flux <dir>");
    let reps: [(&str, Vec<(usize, usize)>); 11] = [
        ("none", vec![]), ("e01", vec![(0, 1)]), ("e01-02", vec![(0, 1), (0, 2)]), ("e01-23", vec![(0, 1), (2, 3)]),
        ("tri", vec![(0, 1), (0, 2), (1, 2)]), ("star", vec![(0, 1), (0, 2), (0, 3)]), ("path", vec![(0, 1), (1, 2), (2, 3)]),
        ("c4", vec![(0, 1), (1, 2), (2, 3), (0, 3)]), ("paw", vec![(0, 1), (0, 2), (1, 2), (2, 3)]),
        ("k4-e", vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)]), ("k4", EDGES.to_vec())];
    let classes: Vec<(u8, Vec<(usize, u64)>)> = reps.iter().map(|(n, e)| (canonical(mask_of(e)), read(&format!("{dir}/{n}.txt")))).collect();
    let mut seen: Vec<u8> = classes.iter().map(|c| c.0).collect();
    seen.sort();
    seen.dedup();
    assert_eq!(seen.len(), 11, "class representatives are not distinct");
    let enumerator = |mask: u8| -> &Vec<(usize, u64)> { let c = canonical(mask); &classes.iter().find(|x| x.0 == c).unwrap().1 };
    let log_sum = |a: &Vec<(usize, u64)>, beta: f64| -> f64 {
        let t: Vec<f64> = a.iter().map(|&(w, c)| (c as f64).ln() - 2.0 * beta * w as f64).collect();
        let m = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        m + t.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
    };
    println!("beta\tF_m(one plane)\tF_e(one unit)\tF_e(e=(1,1,0))\tF_e(e=(1,1,1))");
    for &beta in &[0.05, 0.1, 0.2, 0.3, 0.35, 0.4, 0.42, 0.44, 0.46, 0.5, 0.6, 0.8, 1.0, 1.5] {
        let l0 = log_sum(enumerator(0), beta);
        let z = |mask: u8| (log_sum(enumerator(mask), beta) - l0).exp();
        // temporal planes 01, 02, 03 are bits 0, 1, 2; spatial magnetic twist m = 0
        let ze = |e: [u8; 3]| -> f64 {
            (0u8..8).map(|k| { let sign = if (e[0] & k & 1) ^ (e[1] & (k >> 1) & 1) ^ (e[2] & (k >> 2) & 1) == 1 { -1.0 } else { 1.0 }; sign * z(k) }).sum::<f64>() / 8.0
        };
        let z0 = ze([0, 0, 0]);
        let fe = |e: [u8; 3]| -(ze(e) / z0).ln();
        let fm = -(z(mask_of(&[(1, 2)]))).ln();
        println!("{beta}\t{fm:.6e}\t{:.6e}\t{:.6e}\t{:.6e}", fe([1, 0, 0]), fe([1, 1, 0]), fe([1, 1, 1]));
    }
}
