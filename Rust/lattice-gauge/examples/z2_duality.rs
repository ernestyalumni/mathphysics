//! Exact check of Kramers-Wannier duality with all 't Hooft sectors for Z2 gauge theory on the 2^4 torus.
//! Usage: z2_duality <dir with z2_exact outputs for the 11 twist classes>
//!
//! C = im(delta) in F2^96 is the code of untwisted plaquette configurations (rank 45). Its dual code C^perp is the space
//! of mod-2 2-cycles. On the 2^4 torus the dual lattice is again a 2^4 torus, and a 2-cycle of the lattice is a
//! 2-cocycle of the dual lattice; the 2-cocycles are the union of the 2^6 cosets t + im(delta) over the classes t of
//! H^2(T^4; Z2), i.e. the twisted sectors. The MacWilliams identity therefore predicts, for every weight v,
//!   sum_w A_w K_v(w) = |C| * sum_{t in H^2} A^(t)_v,     K_v(w) = sum_j (-1)^j C(w, j) C(96 - w, v - j),
//! where the sum over the 64 twists is the sum over the 11 isomorphism classes of graphs on 4 vertices (edges = planes),
//! each weighted by its number of labelled graphs. Both sides are integers below 2^142 in absolute value; the identity is
//! verified modulo three primes near 2^61, whose product exceeds 2^182, so agreement proves exact equality.

fn read(path: &str) -> Vec<(usize, u64)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
    text.lines().skip(1).map(|l| { let mut it = l.split('\t'); (it.next().unwrap().parse().unwrap(), it.next().unwrap().parse().unwrap()) }).collect()
}

fn main() {
    let dir = std::env::args().nth(1).expect("usage: z2_duality <dir>");
    let n = 96usize;
    let classes: [(&str, u64); 11] = [("none", 1), ("e01", 6), ("e01-02", 12), ("e01-23", 3), ("tri", 4), ("star", 4),
        ("path", 12), ("c4", 3), ("paw", 12), ("k4-e", 6), ("k4", 1)];
    assert_eq!(classes.iter().map(|c| c.1).sum::<u64>(), 64);
    let untwisted = read(&format!("{dir}/none.txt"));
    let size: u64 = untwisted.iter().map(|x| x.1).sum();
    assert_eq!(size, 1u64 << 45);
    let mut sectors = vec![0u128; n + 1];
    for (name, mult) in classes {
        for (v, a) in read(&format!("{dir}/{name}.txt")) { sectors[v] += mult as u128 * a as u128; }
    }
    assert_eq!(sectors.iter().sum::<u128>(), 1u128 << 51, "64 cosets of size 2^45");
    let primes: [u128; 3] = [2305843009213693951, 2305843009213693921, 2305843009213693907]; // 2^61 - 1, -31, -45
    let mut all_ok = true;
    for &p in &primes {
        // binomials mod p by Pascal's triangle
        let mut c = vec![vec![0u128; n + 1]; n + 1];
        for i in 0..=n { c[i][0] = 1; for j in 1..=i { c[i][j] = (c[i - 1][j - 1] + c[i - 1][j]) % p; } }
        let mut bad = 0;
        for v in 0..=n {
            let mut lhs: u128 = 0;
            for &(w, a) in &untwisted {
                let mut k: u128 = 0;
                for j in 0..=v.min(w) {
                    if v - j > n - w { continue; }
                    let t = c[w][j] * c[n - w][v - j] % p;
                    k = if j % 2 == 0 { (k + t) % p } else { (k + p - t) % p };
                }
                lhs = (lhs + (a as u128 % p) * k) % p;
            }
            let rhs = ((size as u128) % p) * (sectors[v] % p) % p;
            if lhs != rhs { bad += 1; }
        }
        println!("prime {p}: {} of {} weights agree", n + 1 - bad, n + 1);
        all_ok &= bad == 0;
    }
    println!("{}", if all_ok { "DUALITY VERIFIED EXACTLY for all 97 weights" } else { "DUALITY FAILS" });
}
