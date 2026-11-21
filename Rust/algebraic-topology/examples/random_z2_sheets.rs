//! Control for the vortex-sheet homology: sheets from independent random Z_2 links (no dynamics) on L^4.
//! Prints b2 over Z_2 and over Q of the dual surface. Usage: random_z2_sheets <L> <flip probability> <samples> <seed>
use algebraic_topology::surface::{analyze, chain_complex};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let l: usize = a[1].parse().unwrap();
    let p: f64 = a[2].parse().unwrap();
    let samples: usize = a[3].parse().unwrap();
    let mut state: u64 = a[4].parse::<u64>().unwrap() ^ 0x9E3779B97F4A7C15;
    let mut rnd = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; (state >> 11) as f64 / (1u64 << 53) as f64 };
    let sides = [l; 4];
    let mut strides = [1usize; 4];
    for i in 1..4 { strides[i] = strides[i - 1] * sides[i - 1]; }
    let n = l.pow(4);
    let fwd = |s: usize, d: usize| { let c = (s / strides[d]) % l; if c + 1 == l { s - c * strides[d] } else { s + strides[d] } };
    let back = |s: usize, d: usize| { let c = (s / strides[d]) % l; if c == 0 { s + (l - 1) * strides[d] } else { s - strides[d] } };
    for k in 0..samples {
        let z: Vec<i8> = (0..4 * n).map(|_| if rnd() < p { -1 } else { 1 }).collect();
        let mut sq = Vec::new();
        for x in 0..n { for mu in 0..4 { for nu in mu + 1..4 {
            if z[4 * x + mu] * z[4 * fwd(x, mu) + nu] * z[4 * fwd(x, nu) + mu] * z[4 * x + nu] < 0 {
                let comp: Vec<usize> = (0..4).filter(|&d| d != mu && d != nu).collect();
                let (r, s) = (comp[0], comp[1]);
                let b = back(x, r);
                sq.push([x, b, back(b, s), back(x, s)]);
            }
        }}}
        let comps = analyze(&sq);
        let cc = chain_complex(&sq);
        let (b2, bq) = (cc.betti_mod_p(2)[2], cc.betti_mod_p(1_000_003)[2]);
        println!("L={l} p={p} sample {k}: squares {} components {} b2(Z2)={b2} b2(Q)={bq} difference={}", sq.len(), comps.len(), b2 as i64 - bq as i64);
    }
}
