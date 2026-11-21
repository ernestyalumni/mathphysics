//! Cubical complex of the periodic d-torus with sides n_i >= 2, using the site conventions of
//! `Rust/lattice-gauge` (site = sum c_i stride_i, stride_0 = 1). A k-cell is (base site, k sorted directions);
//! its index is site * C(d,k) + (rank of the direction subset in lexicographic order).
//! Boundary: d[x; s_1 < ... < s_k] = sum_j (-1)^j ([x + e_{s_j}; S minus s_j] - [x; S minus s_j]), j from 0.

use crate::complex::ChainComplex;

pub struct CubicalTorus {
    pub sides: Vec<usize>,
    strides: Vec<usize>,
    pub sites: usize,
    subsets: Vec<Vec<Vec<usize>>>, // subsets[k] = sorted k-subsets of 0..d in lexicographic order
}

fn k_subsets(d: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    fn rec(start: usize, d: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k { out.push(cur.clone()); return; }
        for i in start..d { cur.push(i); rec(i + 1, d, k, cur, out); cur.pop(); }
    }
    rec(0, d, k, &mut cur, &mut out);
    out
}

impl CubicalTorus {
    pub fn new(sides: &[usize]) -> Self {
        assert!(!sides.is_empty() && sides.iter().all(|&n| n >= 2), "sides must be at least 2");
        let mut strides = vec![1; sides.len()];
        for i in 1..sides.len() { strides[i] = strides[i - 1] * sides[i - 1]; }
        let sites = sides.iter().product();
        let d = sides.len();
        CubicalTorus { sides: sides.to_vec(), strides, sites, subsets: (0..=d).map(|k| k_subsets(d, k)).collect() }
    }

    pub fn dimension(&self) -> usize { self.sides.len() }

    pub fn shift(&self, site: usize, dir: usize) -> usize {
        let c = (site / self.strides[dir]) % self.sides[dir];
        if c + 1 == self.sides[dir] { site - c * self.strides[dir] } else { site + self.strides[dir] }
    }

    pub fn coordinate(&self, site: usize, dir: usize) -> usize { (site / self.strides[dir]) % self.sides[dir] }

    pub fn cell_index(&self, site: usize, dirs: &[usize]) -> usize {
        let k = dirs.len();
        let r = self.subsets[k].iter().position(|s| s.as_slice() == dirs).expect("directions must be sorted and distinct");
        site * self.subsets[k].len() + r
    }

    pub fn cell(&self, k: usize, index: usize) -> (usize, &[usize]) {
        let m = self.subsets[k].len();
        (index / m, &self.subsets[k][index % m])
    }

    pub fn chain_complex(&self) -> ChainComplex {
        let d = self.dimension();
        let counts: Vec<usize> = (0..=d).map(|k| self.sites * self.subsets[k].len()).collect();
        let mut boundary = vec![Vec::new()];
        for k in 1..=d {
            let mut bk = Vec::with_capacity(counts[k]);
            for idx in 0..counts[k] {
                let (x, s) = self.cell(k, idx);
                let mut faces = Vec::with_capacity(2 * k);
                for j in 0..k {
                    let rest: Vec<usize> = s.iter().cloned().filter(|&v| v != s[j]).collect();
                    let sign = if j % 2 == 0 { 1 } else { -1 };
                    faces.push((self.cell_index(self.shift(x, s[j]), &rest), sign));
                    faces.push((self.cell_index(x, &rest), -sign));
                }
                bk.push(faces);
            }
            boundary.push(bk);
        }
        ChainComplex { cell_counts: counts, boundary }
    }
}
