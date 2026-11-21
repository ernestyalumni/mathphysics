//! Finite chain complexes C_d -> ... -> C_0 with sparse integer boundary maps, and their homology.

use crate::smith::{invariant_factors, rank_mod_p, rank_sparse_mod_p};

/// boundary[k][j] = boundary of the j-th k-cell as (index of (k-1)-cell, coefficient); boundary[0] is empty.
#[derive(Clone, Debug)]
pub struct ChainComplex {
    pub cell_counts: Vec<usize>,
    pub boundary: Vec<Vec<Vec<(usize, i64)>>>,
}

/// H_k over Z: free rank (Betti number) and torsion coefficients (> 1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Homology {
    pub betti: usize,
    pub torsion: Vec<i64>,
}

impl ChainComplex {
    pub fn dimension(&self) -> usize {
        self.cell_counts.len() - 1
    }

    /// Dense matrix of the boundary map C_k -> C_{k-1} (rows: (k-1)-cells, columns: k-cells).
    pub fn boundary_matrix(&self, k: usize) -> (usize, usize, Vec<i64>) {
        let (rows, cols) = (self.cell_counts[k - 1], self.cell_counts[k]);
        let mut m = vec![0i64; rows * cols];
        for (j, faces) in self.boundary[k].iter().enumerate() {
            for &(i, c) in faces {
                m[i * cols + j] += c;
            }
        }
        (rows, cols, m)
    }

    /// True when every composite boundary of boundary vanishes.
    pub fn is_complex(&self) -> bool {
        (2..=self.dimension()).all(|k| {
            self.boundary[k].iter().all(|faces| {
                let mut acc = std::collections::HashMap::new();
                for &(i, c) in faces {
                    for &(h, e) in &self.boundary[k - 1][i] {
                        *acc.entry(h).or_insert(0i64) += c * e;
                    }
                }
                acc.values().all(|&v| v == 0)
            })
        })
    }

    pub fn euler_characteristic(&self) -> i64 {
        self.cell_counts.iter().enumerate().map(|(k, &n)| if k % 2 == 0 { n as i64 } else { -(n as i64) }).sum()
    }

    /// Integral homology H_0 .. H_d by Smith normal form (dense; for small and moderate complexes).
    pub fn homology_z(&self) -> Vec<Homology> {
        let d = self.dimension();
        let factors: Vec<Vec<i64>> = (0..=d + 1)
            .map(|k| if k == 0 || k > d { vec![] } else { let (r, c, m) = self.boundary_matrix(k); invariant_factors(r, c, &m) })
            .collect();
        (0..=d)
            .map(|k| {
                let rank_k = factors[k].len();
                let rank_k1 = factors[k + 1].len();
                Homology { betti: self.cell_counts[k] - rank_k - rank_k1, torsion: factors[k + 1].iter().cloned().filter(|&x| x > 1).collect() }
            })
            .collect()
    }

    /// Betti numbers over Z_p (p prime).
    pub fn betti_mod_p(&self, p: i64) -> Vec<usize> {
        let d = self.dimension();
        let ranks: Vec<usize> = (0..=d + 1)
            .map(|k| {
                if k == 0 || k > d { return 0; }
                if self.cell_counts[k] * self.cell_counts[k - 1] > 4_000_000 {
                    rank_sparse_mod_p(&self.boundary[k], p)
                } else {
                    let (r, c, m) = self.boundary_matrix(k);
                    rank_mod_p(r, c, &m, p)
                }
            })
            .collect();
        (0..=d).map(|k| self.cell_counts[k] - ranks[k] - ranks[k + 1]).collect()
    }
}
