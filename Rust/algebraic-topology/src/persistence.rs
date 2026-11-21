//! Persistent homology over Z_2 by the standard column reduction (Edelsbrunner–Letscher–Zomorodian; Zomorodian–Carlsson)
//! with the clearing optimization. Input: a chain complex and a filtration value per cell that is monotone
//! (every face enters no later than its coface). Output: finite pairs (dim, birth, death) and essential classes.
//! `lower_star_cubical` builds the lower-star filtration of a function on the sites of a cubical torus (each cell gets the
//! maximum of its vertex values), for sublevel sets {f <= t}.

use crate::complex::ChainComplex;
use crate::cubical::CubicalTorus;

#[derive(Clone, Debug, PartialEq)]
pub struct Persistence {
    /// (dimension, birth, death) with death > birth or equal (zero-length pairs are kept)
    pub pairs: Vec<(usize, f64, f64)>,
    /// (dimension, birth) of classes that never die
    pub essential: Vec<(usize, f64)>,
}

impl Persistence {
    /// Pairs of the given dimension with persistence (death - birth) above `min_persistence`.
    pub fn bars(&self, dim: usize, min_persistence: f64) -> Vec<(f64, f64)> {
        self.pairs.iter().filter(|p| p.0 == dim && p.2 - p.1 > min_persistence).map(|p| (p.1, p.2)).collect()
    }
}

/// Persistent homology over Z_2. `values[k][j]` is the filtration value of the j-th k-cell.
pub fn persistence(cx: &ChainComplex, values: &[Vec<f64>]) -> Persistence {
    let d = cx.dimension();
    // global order: by value, then dimension, then index (faces precede cofaces at equal value)
    let mut order: Vec<(usize, usize)> = Vec::new();
    for k in 0..=d { for j in 0..cx.cell_counts[k] { order.push((k, j)); } }
    order.sort_by(|a, b| values[a.0][a.1].partial_cmp(&values[b.0][b.1]).unwrap().then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
    let mut pos: Vec<Vec<usize>> = (0..=d).map(|k| vec![0; cx.cell_counts[k]]).collect();
    for (i, &(k, j)) in order.iter().enumerate() { pos[k][j] = i; }
    for k in 1..=d {
        for j in 0..cx.cell_counts[k] {
            for &(f, _) in &cx.boundary[k][j] {
                assert!(pos[k - 1][f] < pos[k][j], "filtration is not monotone");
            }
        }
    }
    let n = order.len();
    let mut columns: Vec<Vec<usize>> = vec![Vec::new(); n]; // reduced boundary columns as sorted positions (Z_2)
    let mut low_owner: Vec<Option<usize>> = vec![None; n]; // pivot row -> column
    let mut paired = vec![false; n];
    let mut pairs = Vec::new();
    // process dimensions from high to low for clearing: a column whose index is already a pivot row is zero
    for k in (1..=d).rev() {
        for i in 0..n {
            let (kk, j) = order[i];
            if kk != k || paired[i] { continue; }
            let mut col: Vec<usize> = cx.boundary[k][j].iter().filter(|&&(_, c)| c.rem_euclid(2) == 1).map(|&(f, _)| pos[k - 1][f]).collect();
            col.sort_unstable();
            loop {
                let Some(&low) = col.last() else { break };
                match low_owner[low] {
                    Some(other) => {
                        // col += columns[other] (symmetric difference)
                        let b = &columns[other];
                        let mut merged = Vec::with_capacity(col.len() + b.len());
                        let (mut x, mut y) = (0, 0);
                        while x < col.len() || y < b.len() {
                            if y == b.len() || (x < col.len() && col[x] < b[y]) { merged.push(col[x]); x += 1; }
                            else if x == col.len() || b[y] < col[x] { merged.push(b[y]); y += 1; }
                            else { x += 1; y += 1; }
                        }
                        col = merged;
                    }
                    None => {
                        low_owner[low] = Some(i);
                        paired[low] = true;
                        paired[i] = true;
                        let (lk, lj) = order[low];
                        pairs.push((lk, values[lk][lj], values[k][j]));
                        break;
                    }
                }
            }
            columns[i] = col;
        }
    }
    let essential = (0..n).filter(|&i| !paired[i]).map(|i| { let (k, j) = order[i]; (k, values[k][j]) }).collect();
    Persistence { pairs, essential }
}

/// Lower-star filtration on a cubical torus from a function on sites: each cell gets the maximum over its vertices.
pub fn lower_star_cubical(torus: &CubicalTorus, f: &[f64]) -> (ChainComplex, Vec<Vec<f64>>) {
    let cx = torus.chain_complex();
    let d = torus.dimension();
    let values = (0..=d).map(|k| (0..cx.cell_counts[k]).map(|idx| {
        let (x, dirs) = torus.cell(k, idx);
        let mut m = f64::NEG_INFINITY;
        for mask in 0u32..(1 << dirs.len()) {
            let mut y = x;
            for (b, &dir) in dirs.iter().enumerate() { if mask & (1 << b) != 0 { y = torus.shift(y, dir); } }
            m = m.max(f[y]);
        }
        m
    }).collect()).collect();
    (cx, values)
}
