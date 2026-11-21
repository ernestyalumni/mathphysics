//! Simplicial complexes from facets; Z_p and integer cochains; the Alexander–Whitney cup product;
//! the Freudenthal triangulation of a periodic d-torus with integer seam cocycles.
//!
//! Simplices are sorted vertex lists; orientation follows the vertex order. Cup product of a p-cochain a and a
//! q-cochain b: (a cup b)[v_0..v_{p+q}] = a[v_0..v_p] b[v_p..v_{p+q}]; delta is the transpose of the boundary.

use crate::complex::ChainComplex;
use std::collections::HashMap;

pub struct SimplicialComplex {
    pub simplices: Vec<Vec<Vec<usize>>>, // by dimension, each sorted
    index: Vec<HashMap<Vec<usize>, usize>>,
}

impl SimplicialComplex {
    pub fn from_facets(facets: &[Vec<usize>]) -> Self {
        let dim = facets.iter().map(|f| f.len() - 1).max().unwrap_or(0);
        let mut index: Vec<HashMap<Vec<usize>, usize>> = vec![HashMap::new(); dim + 1];
        let mut simplices: Vec<Vec<Vec<usize>>> = vec![Vec::new(); dim + 1];
        for f in facets {
            let mut f = f.clone();
            f.sort_unstable();
            let n = f.len();
            for mask in 1u32..(1 << n) {
                let s: Vec<usize> = (0..n).filter(|i| mask & (1 << i) != 0).map(|i| f[i]).collect();
                let k = s.len() - 1;
                if !index[k].contains_key(&s) {
                    index[k].insert(s.clone(), simplices[k].len());
                    simplices[k].push(s);
                }
            }
        }
        SimplicialComplex { simplices, index }
    }

    pub fn dimension(&self) -> usize { self.simplices.len() - 1 }

    pub fn index_of(&self, s: &[usize]) -> Option<usize> {
        self.index.get(s.len() - 1).and_then(|m| m.get(s).copied())
    }

    pub fn chain_complex(&self) -> ChainComplex {
        let d = self.dimension();
        let mut boundary = vec![Vec::new()];
        for k in 1..=d {
            boundary.push(
                self.simplices[k]
                    .iter()
                    .map(|s| (0..=k).map(|i| {
                        let face: Vec<usize> = s.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).collect();
                        (self.index[k - 1][&face], if i % 2 == 0 { 1 } else { -1 })
                    }).collect())
                    .collect(),
            );
        }
        ChainComplex { cell_counts: self.simplices.iter().map(|v| v.len()).collect(), boundary }
    }

    /// Coboundary of an integer k-cochain (values indexed like simplices[k]).
    pub fn coboundary(&self, k: usize, a: &[i64]) -> Vec<i64> {
        self.simplices[k + 1].iter().map(|s| {
            (0..=k + 1).map(|i| {
                let face: Vec<usize> = s.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).collect();
                let sign = if i % 2 == 0 { 1 } else { -1 };
                sign * a[self.index[k][&face]]
            }).sum()
        }).collect()
    }

    /// Alexander–Whitney cup product of integer cochains of degrees p and q.
    pub fn cup(&self, p: usize, a: &[i64], q: usize, b: &[i64]) -> Vec<i64> {
        self.simplices[p + q].iter().map(|s| {
            a[self.index[p][&s[..=p].to_vec()]] * b[self.index[q][&s[p..].to_vec()]]
        }).collect()
    }

    /// Steenrod cup-1 product of integer cochains a (degree p >= 1) and b (degree q >= 1), a (p+q-1)-cochain:
    /// (a cup_1 b)[v_0..v_n] = sum_{i=0}^{p-1} (-1)^{(p-i)(q+1)} a[v_0..v_i, v_{i+q}..v_n] b[v_i..v_{i+q}], n = p+q-1.
    /// It satisfies delta(a cup_1 b) = (-1)^{p+q-1} a cup b + (-1)^{p+q+pq} b cup a + delta a cup_1 b + (-1)^p a cup_1 delta b
    /// (checked in tests/steenrod.rs).
    pub fn cup1(&self, p: usize, a: &[i64], q: usize, b: &[i64]) -> Vec<i64> {
        assert!(p >= 1 && q >= 1);
        let n = p + q - 1;
        self.simplices[n].iter().map(|s| {
            let mut total = 0i64;
            for i in 0..p {
                let mut fa: Vec<usize> = s[..=i].to_vec();
                fa.extend_from_slice(&s[i + q..]);
                let fb: Vec<usize> = s[i..=i + q].to_vec();
                let sign = if ((p - i) * (q + 1)) % 2 == 0 { 1 } else { -1 };
                total += sign * a[self.index[p][&fa]] * b[self.index[q][&fb]];
            }
            total
        }).collect()
    }

    /// Pontryagin square of a mod-2 2-cocycle given by an integer lift w (delta w = 0 mod 2):
    /// P(w) = w cup w + w cup_1 delta w, a 4-cochain whose class mod 4 depends only on the class of w mod 2.
    pub fn pontryagin_square(&self, w: &[i64]) -> Vec<i64> {
        let dw = self.coboundary(2, w);
        let a = self.cup(2, w, 2, w);
        let b = self.cup1(2, w, 3, &dw);
        a.iter().zip(b).map(|(x, y)| x + y).collect()
    }

    /// Basis of H^k(X; Z_p) as cocycle representatives (values in 0..p), by linear algebra mod p:
    /// a basis of ker delta_k completed modulo im delta_{k-1}.
    pub fn cohomology_basis_mod_p(&self, k: usize, p: i64) -> Vec<Vec<i64>> {
        let n = self.simplices[k].len();
        // ker delta_k: solve D x = 0 mod p where D rows = (k+1)-simplices
        let rows: Vec<Vec<i64>> = if k < self.dimension() {
            self.simplices[k + 1].iter().map(|s| {
                let mut r = vec![0i64; n];
                for i in 0..=k + 1 {
                    let face: Vec<usize> = s.iter().enumerate().filter(|&(j, _)| j != i).map(|(_, &v)| v).collect();
                    let sign = if i % 2 == 0 { 1 } else { p - 1 };
                    let c = &mut r[self.index[k][&face]];
                    *c = (*c + sign) % p;
                }
                r
            }).collect()
        } else { vec![] };
        let kernel = nullspace_mod_p(&rows, n, p);
        // image of delta_{k-1}: columns are delta of basis (k-1)-cochains
        let mut image: Vec<Vec<i64>> = Vec::new();
        if k > 0 {
            for j in 0..self.simplices[k - 1].len() {
                let mut e = vec![0i64; self.simplices[k - 1].len()];
                e[j] = 1;
                image.push(self.coboundary(k - 1, &e).iter().map(|v| v.rem_euclid(p)).collect());
            }
        }
        // greedily add kernel vectors independent of (image + chosen)
        let mut basis = Vec::new();
        let mut span = image.clone();
        let mut r0 = rank_rows_mod_p(&span, n, p);
        for v in kernel {
            span.push(v.clone());
            let r1 = rank_rows_mod_p(&span, n, p);
            if r1 > r0 { basis.push(v); r0 = r1; } else { span.pop(); }
        }
        basis
    }
}

fn rank_rows_mod_p(rows: &[Vec<i64>], n: usize, p: i64) -> usize {
    let flat: Vec<i64> = rows.iter().flat_map(|r| r.iter().cloned()).collect();
    crate::smith::rank_mod_p(rows.len(), n, &flat, p)
}

/// Basis of { x in Z_p^n : rows . x = 0 }.
pub fn nullspace_mod_p(rows: &[Vec<i64>], n: usize, p: i64) -> Vec<Vec<i64>> {
    let mut a: Vec<Vec<i64>> = rows.iter().map(|r| r.iter().map(|v| v.rem_euclid(p)).collect()).collect();
    let inv = |x: i64| { let (mut b, mut e, mut r) = (x % p, p - 2, 1i64); while e > 0 { if e & 1 == 1 { r = r * b % p; } b = b * b % p; e >>= 1; } r };
    let mut pivots = Vec::new();
    let mut rank = 0;
    for c in 0..n {
        let Some(pr) = (rank..a.len()).find(|&r| a[r][c] != 0) else { continue };
        a.swap(rank, pr);
        let iv = inv(a[rank][c]);
        for k in 0..n { a[rank][k] = a[rank][k] * iv % p; }
        for r in 0..a.len() {
            if r != rank && a[r][c] != 0 {
                let f = a[r][c];
                for k in 0..n { a[r][k] = (a[r][k] - f * a[rank][k]).rem_euclid(p); }
            }
        }
        pivots.push(c);
        rank += 1;
    }
    let free: Vec<usize> = (0..n).filter(|c| !pivots.contains(c)).collect();
    free.iter().map(|&f| {
        let mut x = vec![0i64; n];
        x[f] = 1;
        for (r, &pc) in pivots.iter().enumerate() { x[pc] = (p - a[r][f]) % p; }
        x
    }).collect()
}

/// Freudenthal triangulation of the d-torus with sides n_i >= 3: each unit cube at base x splits into d! simplices
/// x = w_0 < w_1 < ... < w_d along a permutation pi of the directions, w_{i+1} = w_i + e_{pi(i)}. Returns the complex,
/// the top simplices with their orientation sign sign(pi) (the signed sum is a fundamental cycle), and for each
/// directed edge the seam crossings per direction (for integer seam cocycles).
pub struct TorusTriangulation {
    pub complex: SimplicialComplex,
    pub top: Vec<(Vec<usize>, i64)>, // sorted vertex list and orientation sign relative to the sorted order
    pub sides: Vec<usize>,
    edge_seams: HashMap<(usize, usize), Vec<i64>>, // (u, v) with u < v as vertex ids -> crossing count per direction, oriented u -> v
}

impl TorusTriangulation {
    pub fn new(sides: &[usize]) -> Self {
        assert!(sides.iter().all(|&n| n >= 3), "Freudenthal triangulation of a torus needs sides >= 3");
        let d = sides.len();
        let mut strides = vec![1; d];
        for i in 1..d { strides[i] = strides[i - 1] * sides[i - 1]; }
        let sites: usize = sides.iter().product();
        let coord = |s: usize, i: usize| (s / strides[i]) % sides[i];
        let step = |s: usize, i: usize| if coord(s, i) + 1 == sides[i] { s - coord(s, i) * strides[i] } else { s + strides[i] };
        let mut perms: Vec<Vec<usize>> = vec![vec![]];
        for _ in 0..d {
            let mut next = Vec::new();
            for p in &perms { for i in 0..d { if !p.contains(&i) { let mut q = p.clone(); q.push(i); next.push(q); } } }
            perms = next;
        }
        let perm_sign = |p: &[usize]| { let mut s = 1; for i in 0..p.len() { for j in i + 1..p.len() { if p[i] > p[j] { s = -s; } } } s };
        let mut facets = Vec::new();
        let mut top = Vec::new();
        let mut edge_seams: HashMap<(usize, usize), Vec<i64>> = HashMap::new();
        for x in 0..sites {
            for p in &perms {
                let mut path = vec![x];
                let mut seam = vec![vec![0i64; d]]; // cumulative crossings from x
                for &dir in p {
                    let last = *path.last().unwrap();
                    let mut c = seam.last().unwrap().clone();
                    if coord(last, dir) + 1 == sides[dir] { c[dir] += 1; }
                    path.push(step(last, dir));
                    seam.push(c);
                }
                // record seam crossings for every edge of this simplex (from path order i < j)
                for i in 0..=d { for j in i + 1..=d {
                    let (u, v) = (path[i], path[j]);
                    let diff: Vec<i64> = (0..d).map(|k| seam[j][k] - seam[i][k]).collect();
                    let (key, val) = if u < v { ((u, v), diff) } else { ((v, u), diff.iter().map(|z| -z).collect()) };
                    edge_seams.entry(key).or_insert(val);
                }}
                // orientation: path order vs sorted vertex order, times sign(pi)
                let mut sorted = path.clone();
                sorted.sort_unstable();
                let pos: Vec<usize> = path.iter().map(|v| sorted.iter().position(|w| w == v).unwrap()).collect();
                let s = perm_sign(&pos) * perm_sign(p);
                facets.push(path.clone());
                top.push((sorted, s));
            }
        }
        TorusTriangulation { complex: SimplicialComplex::from_facets(&facets), top, sides: sides.to_vec(), edge_seams }
    }

    /// Integer 1-cocycle x_dir: value on the edge [u < v] = seam crossings in direction dir going from u to v.
    pub fn seam_cocycle(&self, dir: usize) -> Vec<i64> {
        self.complex.simplices[1].iter().map(|e| self.edge_seams[&(e[0], e[1])][dir]).collect()
    }

    /// Evaluate a top-degree integer cochain on the fundamental cycle (signed sum of top simplices).
    pub fn evaluate_top(&self, c: &[i64]) -> i64 {
        let d = self.sides.len();
        self.top.iter().map(|(s, sign)| sign * c[self.complex.index_of(s).unwrap()]).sum::<i64>() * if d == 0 { 1 } else { 1 }
    }
}
