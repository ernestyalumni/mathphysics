//! Exact Z2 lattice gauge theory on small d-dimensional tori as a binary linear code.
//!
//! Links carry sigma = +-1, written additively in F2 (0 for +1). The plaquette variable is the F2 sum of its four
//! links counted with multiplicity, so on a side of length 1, where U(y + e_mu, nu) = U(y, nu), repeated links cancel.
//! The map from link configurations to plaquette configurations is the coboundary delta: F2^links -> F2^plaquettes. Its
//! image C (rank r) is a binary linear code; each codeword is hit by 2^(links - r) link configurations. With the
//! energy sum_p (1 - sigma_p), a plaquette with value 1 costs 2, so with normalized counting measure on the links
//!   Z(beta) = 2^(-r) sum_w A_w exp(-2 beta w),   A_w = #{c in C : |c| = w},
//! and Z(0) = 1. A twist t in F2^plaquettes (plaquettes whose sign is reversed) gives the coset t + C instead.
//! The weight enumerator A_w is the exact answer for every beta; it is computed by a Gray-code walk over C.
//! Duality check (MacWilliams): the enumerator of the dual code C^perp = ker(delta^T) satisfies
//!   B(x, y) = |C|^(-1) A(x + y, x - y),  A(x, y) = sum_w A_w x^(P - w) y^w.

/// Geometry of a d-dimensional torus with sides >= 1, links indexed (site, direction).
pub struct Box {
    pub sides: Vec<usize>,
    pub sites: usize,
    /// For each plaquette (site, mu < nu): the F2 boundary as a set of link indices (multiplicity mod 2 applied).
    pub plaquettes: Vec<Vec<usize>>,
    /// (site, mu, nu) of each plaquette.
    pub plaquette_labels: Vec<(usize, usize, usize)>,
    pub links: usize,
}

impl Box {
    pub fn new(sides: &[usize]) -> Box {
        assert!(!sides.is_empty() && sides.iter().all(|&s| s >= 1));
        let d = sides.len();
        let sites: usize = sides.iter().product();
        let coord = |s: usize, k: usize| -> usize { (s / sides[..k].iter().product::<usize>()) % sides[k] };
        let shift = |s: usize, k: usize| -> usize {
            let stride: usize = sides[..k].iter().product();
            let c = coord(s, k);
            s - c * stride + ((c + 1) % sides[k]) * stride
        };
        let link = |s: usize, k: usize| s * d + k;
        let mut plaquettes = Vec::new();
        let mut labels = Vec::new();
        for s in 0..sites {
            for mu in 0..d {
                for nu in mu + 1..d {
                    let four = [link(s, mu), link(shift(s, mu), nu), link(shift(s, nu), mu), link(s, nu)];
                    let mut set: Vec<usize> = Vec::new();
                    for l in four {
                        if let Some(i) = set.iter().position(|&x| x == l) { set.remove(i); } else { set.push(l); }
                    }
                    set.sort();
                    plaquettes.push(set);
                    labels.push((s, mu, nu));
                }
            }
        }
        Box { sides: sides.to_vec(), sites, plaquettes, plaquette_labels: labels, links: sites * d }
    }

    pub fn plaquette_count(&self) -> usize { self.plaquettes.len() }

    /// Columns of delta: for each link, the set of plaquettes containing it (mod 2), as a 128-bit mask.
    pub fn coboundary_columns(&self) -> Vec<u128> {
        assert!(self.plaquette_count() <= 128, "at most 128 plaquettes");
        let mut cols = vec![0u128; self.links];
        for (p, set) in self.plaquettes.iter().enumerate() {
            for &l in set { cols[l] ^= 1u128 << p; }
        }
        cols
    }

    /// Plaquettes in the plane (mu, nu) at the sites whose mu- and nu-coordinates are 0: the 't Hooft twist surface.
    pub fn twist_mask(&self, mu: usize, nu: usize) -> u128 {
        let (a, b) = (mu.min(nu), mu.max(nu));
        let mut m = 0u128;
        for (p, &(s, x, y)) in self.plaquette_labels.iter().enumerate() {
            if (x, y) == (a, b) {
                let c = |k: usize| (s / self.sides[..k].iter().product::<usize>()) % self.sides[k];
                if c(a) == 0 && c(b) == 0 { m |= 1u128 << p; }
            }
        }
        m
    }
}

/// A basis of the F2 span of `vectors` (Gaussian elimination on 128-bit words).
pub fn span_basis(vectors: &[u128]) -> Vec<u128> {
    let mut basis: Vec<u128> = Vec::new(); // kept in reduced form by pivot = highest set bit
    for &v in vectors {
        let mut x = v;
        for &b in &basis {
            let pivot = 127 - b.leading_zeros();
            if x >> pivot & 1 == 1 { x ^= b; }
        }
        if x != 0 {
            let pivot = 127 - x.leading_zeros();
            for b in basis.iter_mut() { if *b >> pivot & 1 == 1 { *b ^= x; } }
            basis.push(x);
        }
    }
    basis
}

/// A basis of the orthogonal complement {y : y.c = 0 for all c in span(basis)} in F2^n.
pub fn dual_basis(basis: &[u128], n: usize) -> Vec<u128> {
    // Solve the linear system with rows = basis; free variables give the kernel.
    let mut rows: Vec<u128> = basis.to_vec();
    let mut pivots: Vec<usize> = Vec::new();
    let mut r = 0;
    for col in 0..n {
        if let Some(i) = (r..rows.len()).find(|&i| rows[i] >> col & 1 == 1) {
            rows.swap(r, i);
            for j in 0..rows.len() { if j != r && rows[j] >> col & 1 == 1 { rows[j] ^= rows[r]; } }
            pivots.push(col);
            r += 1;
        }
    }
    let mut out = Vec::new();
    for free in (0..n).filter(|c| !pivots.contains(c)) {
        let mut y = 1u128 << free;
        for (i, &pc) in pivots.iter().enumerate() {
            if rows[i] >> free & 1 == 1 { y |= 1u128 << pc; }
        }
        out.push(y);
    }
    out
}

/// Weight enumerator of the coset `offset + span(basis)` over `n` coordinates, on `threads` CPU threads.
pub fn coset_enumerator(basis: &[u128], offset: u128, n: usize, threads: usize) -> Vec<u64> {
    let r = basis.len();
    assert!(r < 64);
    let total: u64 = 1u64 << r;
    let chunks = (threads as u64 * 64).min(total).max(1);
    let next = std::sync::atomic::AtomicU64::new(0);
    let mut result = vec![0u64; n + 1];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let next = &next;
                scope.spawn(move || {
                    let mut hist = vec![0u64; n + 1];
                    loop {
                        let c = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if c >= chunks { break; }
                        let lo = total / chunks * c + (total % chunks).min(c);
                        let hi = total / chunks * (c + 1) + (total % chunks).min(c + 1);
                        walk(basis, offset, lo, hi, &mut hist);
                    }
                    hist
                })
            })
            .collect();
        for h in handles {
            for (a, b) in result.iter_mut().zip(h.join().unwrap()) { *a += b; }
        }
    });
    result
}

/// Gray-code walk over indices [lo, hi): codeword(i) = offset + sum of basis[k] over the set bits k of gray(i).
fn walk(basis: &[u128], offset: u128, lo: u64, hi: u64, hist: &mut [u64]) {
    if lo >= hi { return; }
    let g = lo ^ (lo >> 1);
    let mut v = offset;
    for (k, b) in basis.iter().enumerate() { if g >> k & 1 == 1 { v ^= b; } }
    hist[v.count_ones() as usize] += 1;
    for i in lo + 1..hi {
        v ^= basis[i.trailing_zeros() as usize];
        hist[v.count_ones() as usize] += 1;
    }
}

/// Direct enumeration over all link configurations (reference for small boxes): histogram of the number of
/// plaquettes with value 1 (twisted by `offset`).
pub fn brute_force(boxed: &Box, offset: u128, threads: usize) -> Vec<u64> {
    let cols = boxed.coboundary_columns();
    let unit: Vec<u128> = cols.clone();
    // Walking every link configuration is the coset enumeration with the (dependent) link columns as "basis".
    assert!(unit.len() < 40, "brute force limited to 39 links");
    coset_enumerator(&unit, offset, boxed.plaquette_count(), threads)
}

/// MacWilliams transform: B_v = |C|^(-1) sum_w A_w K_v(w), with the Krawtchouk polynomial
/// K_v(w) = sum_j (-1)^j binom(w, j) binom(n - w, v - j). Exact in i128; None on overflow or a non-integer result.
pub fn macwilliams(a: &[u64], n: usize) -> Option<Vec<i128>> {
    let size: i128 = a.iter().map(|&x| x as i128).sum();
    let binom = |m: usize, k: usize| -> i128 {
        if k > m { return 0; }
        let mut v: i128 = 1;
        for s in 0..k.min(m - k) { v = v * (m - s) as i128 / (s + 1) as i128; }
        v
    };
    let mut out = Vec::with_capacity(n + 1);
    for v in 0..=n {
        let mut s: i128 = 0;
        for (w, &aw) in a.iter().enumerate() {
            if aw == 0 { continue; }
            let mut k: i128 = 0;
            for j in 0..=v.min(w) {
                let term = binom(w, j).checked_mul(binom(n - w, v - j))?;
                k = if j % 2 == 0 { k.checked_add(term)? } else { k.checked_sub(term)? };
            }
            s = s.checked_add((aw as i128).checked_mul(k)?)?;
        }
        if s % size != 0 { return None; }
        out.push(s / size);
    }
    Some(out)
}

#[cfg(feature = "cuda")]
extern "C" {
    fn z2_enum_coset(lo: *const u64, hi: *const u64, r: i32, off_lo: u64, off_hi: u64, n: i32, out: *mut u64) -> i32;
}

/// GPU version of [`coset_enumerator`] (feature `cuda`).
#[cfg(feature = "cuda")]
pub fn coset_enumerator_gpu(basis: &[u128], offset: u128, n: usize) -> Result<Vec<u64>, String> {
    let lo: Vec<u64> = basis.iter().map(|b| *b as u64).collect();
    let hi: Vec<u64> = basis.iter().map(|b| (*b >> 64) as u64).collect();
    let mut out = vec![0u64; n + 1];
    let code = unsafe {
        z2_enum_coset(lo.as_ptr(), hi.as_ptr(), basis.len() as i32, offset as u64, (offset >> 64) as u64, n as i32, out.as_mut_ptr())
    };
    if code == 0 { Ok(out) } else { Err(format!("z2_enum_coset: CUDA error {code}")) }
}
