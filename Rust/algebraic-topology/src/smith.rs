//! Smith normal form (invariant factors) of integer matrices, and rank over Z_p.

/// Nonzero invariant factors d_1 | d_2 | ... of a dense integer matrix (rows x cols, row-major), all positive.
/// Exact in i128 with overflow checks (panics on overflow; intended for small and moderate complexes).
pub fn invariant_factors(rows: usize, cols: usize, entries: &[i64]) -> Vec<i64> {
    let mut a: Vec<Vec<i128>> = (0..rows).map(|r| entries[r * cols..(r + 1) * cols].iter().map(|&x| x as i128).collect()).collect();
    let mut factors = Vec::new();
    let mut t = 0;
    while t < rows.min(cols) {
        // pivot: smallest nonzero |entry| in the remaining block
        let mut best: Option<(usize, usize)> = None;
        for r in t..rows {
            for c in t..cols {
                if a[r][c] != 0 && best.map_or(true, |(br, bc)| a[r][c].abs() < a[br][bc].abs()) {
                    best = Some((r, c));
                }
            }
        }
        let Some((pr, pc)) = best else { break };
        a.swap(t, pr);
        for row in a.iter_mut() {
            row.swap(t, pc);
        }
        loop {
            let p = a[t][t];
            let mut changed = false;
            // clear column below the pivot
            for r in t + 1..rows {
                if a[r][t] != 0 {
                    let q = a[r][t] / p;
                    for c in t..cols {
                        a[r][c] = a[r][c].checked_sub(q.checked_mul(a[t][c]).expect("overflow")).expect("overflow");
                    }
                    if a[r][t] != 0 {
                        changed = true;
                    }
                }
            }
            // clear row right of the pivot
            for c in t + 1..cols {
                if a[t][c] != 0 {
                    let q = a[t][c] / p;
                    for r in t..rows {
                        a[r][c] = a[r][c].checked_sub(q.checked_mul(a[r][t]).expect("overflow")).expect("overflow");
                    }
                    if a[t][c] != 0 {
                        changed = true;
                    }
                }
            }
            if changed {
                // a smaller remainder appeared in the pivot row or column: move it to the pivot and repeat
                let mut best = (t, t);
                for r in t..rows {
                    if a[r][t] != 0 && a[r][t].abs() < a[best.0][best.1].abs() { best = (r, t); }
                }
                for c in t..cols {
                    if a[t][c] != 0 && a[t][c].abs() < a[best.0][best.1].abs() { best = (t, c); }
                }
                a.swap(t, best.0);
                for row in a.iter_mut() { row.swap(t, best.1); }
                continue;
            }
            // divisibility: if some entry of the remaining block is not divisible by the pivot, add its row
            let p = a[t][t];
            let mut fixed = false;
            'outer: for r in t + 1..rows {
                for c in t + 1..cols {
                    if a[r][c] % p != 0 {
                        for k in t..cols {
                            a[t][k] = a[t][k].checked_add(a[r][k]).expect("overflow");
                        }
                        fixed = true;
                        break 'outer;
                    }
                }
            }
            if !fixed {
                break;
            }
        }
        factors.push(a[t][t].abs() as i64);
        t += 1;
    }
    factors
}

/// Rank of a dense integer matrix over Z_p (p prime).
pub fn rank_mod_p(rows: usize, cols: usize, entries: &[i64], p: i64) -> usize {
    let mut a: Vec<Vec<i64>> = (0..rows).map(|r| entries[r * cols..(r + 1) * cols].iter().map(|&x| x.rem_euclid(p)).collect()).collect();
    let inv = |x: i64| -> i64 { // Fermat
        let (mut b, mut e, mut r) = (x % p, p - 2, 1i64);
        while e > 0 { if e & 1 == 1 { r = r * b % p; } b = b * b % p; e >>= 1; }
        r
    };
    let mut rank = 0;
    for c in 0..cols {
        let Some(pr) = (rank..rows).find(|&r| a[r][c] != 0) else { continue };
        a.swap(rank, pr);
        let iv = inv(a[rank][c]);
        for k in c..cols { a[rank][k] = a[rank][k] * iv % p; }
        for r in 0..rows {
            if r != rank && a[r][c] != 0 {
                let f = a[r][c];
                for k in c..cols { a[r][k] = (a[r][k] - f * a[rank][k]).rem_euclid(p); }
            }
        }
        rank += 1;
    }
    rank
}
