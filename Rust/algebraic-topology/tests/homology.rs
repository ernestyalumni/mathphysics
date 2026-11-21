//! Oracles: homology of tori, the sphere, the real projective plane, the Klein bottle; cup products on RP^2 and T^2.
use algebraic_topology::complex::Homology;
use algebraic_topology::cubical::CubicalTorus;
use algebraic_topology::simplicial::SimplicialComplex;

fn binom(n: usize, k: usize) -> usize { (0..k).fold(1, |a, i| a * (n - i) / (i + 1)) }
fn free(b: usize) -> Homology { Homology { betti: b, torsion: vec![] } }

#[test]
fn cubical_tori() {
    for d in 1..=4 {
        for side in [2usize, 3] {
            let c = CubicalTorus::new(&vec![side; d]).chain_complex();
            assert!(c.is_complex(), "dd != 0 for T^{d}, side {side}");
            assert_eq!(c.euler_characteristic(), 0);
            let h = c.homology_z();
            for k in 0..=d { assert_eq!(h[k], free(binom(d, k)), "H_{k}(T^{d}) side {side}"); }
        }
    }
}

pub fn rp2() -> Vec<Vec<usize>> {
    vec![vec![0,1,2],vec![0,2,3],vec![0,3,4],vec![0,4,5],vec![0,5,1],vec![1,2,4],vec![2,3,5],vec![3,4,1],vec![4,5,2],vec![5,1,3]]
}
pub fn torus7() -> Vec<Vec<usize>> {
    (0..7).flat_map(|i| vec![vec![i, (i + 1) % 7, (i + 3) % 7], vec![i, (i + 2) % 7, (i + 3) % 7]]).collect()
}
/// 3x3 grid with x periodic and y glued with the reflection x -> -x (flip = true: Klein bottle; false: torus).
pub fn grid_surface(flip: bool) -> Vec<Vec<usize>> {
    let n = 3;
    let v = |i: i64, j: i64| -> usize {
        let (mut i, j) = (i, j);
        if j == n && flip { i = -i; }
        (i.rem_euclid(n) + n * j.rem_euclid(n)) as usize
    };
    let mut t = Vec::new();
    for i in 0..n { for j in 0..n {
        t.push(vec![v(i, j), v(i + 1, j), v(i + 1, j + 1)]);
        t.push(vec![v(i, j), v(i + 1, j + 1), v(i, j + 1)]);
    }}
    t
}

#[test]
fn simplicial_surfaces() {
    let s2 = SimplicialComplex::from_facets(&[vec![0,1,2], vec![0,1,3], vec![0,2,3], vec![1,2,3]]).chain_complex();
    assert_eq!(s2.homology_z(), vec![free(1), free(0), free(1)]);
    let p = SimplicialComplex::from_facets(&rp2()).chain_complex();
    assert!(p.is_complex());
    assert_eq!(p.homology_z(), vec![free(1), Homology { betti: 0, torsion: vec![2] }, free(0)]);
    assert_eq!(p.betti_mod_p(2), vec![1, 1, 1]);
    assert_eq!(p.betti_mod_p(3), vec![1, 0, 0]);
    let t = SimplicialComplex::from_facets(&torus7()).chain_complex();
    assert_eq!(t.homology_z(), vec![free(1), free(2), free(1)]);
    let k = SimplicialComplex::from_facets(&grid_surface(true)).chain_complex();
    assert_eq!(k.homology_z(), vec![free(1), Homology { betti: 1, torsion: vec![2] }, free(0)]);
    let t9 = SimplicialComplex::from_facets(&grid_surface(false)).chain_complex();
    assert_eq!(t9.homology_z(), vec![free(1), free(2), free(1)]);
}

/// Evaluate a mod-2 2-cochain on the Z_2 fundamental class (sum of all triangles).
fn eval2(c: &[i64]) -> i64 { c.iter().sum::<i64>().rem_euclid(2) }

#[test]
fn cup_products_mod_2() {
    let p = SimplicialComplex::from_facets(&rp2());
    let h1 = p.cohomology_basis_mod_p(1, 2);
    assert_eq!(h1.len(), 1);
    assert_eq!(eval2(&p.cup(1, &h1[0], 1, &h1[0])), 1, "w cup w must be nonzero on RP^2");
    let t = SimplicialComplex::from_facets(&torus7());
    let h = t.cohomology_basis_mod_p(1, 2);
    assert_eq!(h.len(), 2);
    let (a, b) = (&h[0], &h[1]);
    assert_eq!(eval2(&t.cup(1, a, 1, a)), 0);
    assert_eq!(eval2(&t.cup(1, b, 1, b)), 0);
    assert_eq!(eval2(&t.cup(1, a, 1, b)), 1, "a cup b must be nonzero on T^2");
    // coboundary of a coboundary vanishes
    let e: Vec<i64> = (0..t.simplices[0].len()).map(|i| (i * 7 % 5) as i64).collect();
    let d1 = t.coboundary(0, &e);
    assert!(t.coboundary(1, &d1).iter().all(|&x| x == 0));
}

#[test]
fn sparse_rank_matches_dense() {
    use algebraic_topology::smith::{rank_mod_p, rank_sparse_mod_p};
    let mut s: u64 = 99;
    let mut rnd = || { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); s >> 33 };
    for _ in 0..50 {
        let (rows, cols) = (1 + (rnd() % 30) as usize, 1 + (rnd() % 30) as usize);
        let mut dense = vec![0i64; rows * cols];
        let mut sparse = vec![Vec::new(); cols];
        for c in 0..cols { for r in 0..rows { if rnd() % 4 == 0 { let v = (rnd() % 5) as i64 - 2; dense[r * cols + c] = v; sparse[c].push((r, v)); } } }
        for p in [2i64, 3, 1_000_003] {
            assert_eq!(rank_sparse_mod_p(&sparse, p), rank_mod_p(rows, cols, &dense, p));
        }
    }
    // a big cubical torus through the sparse path: T^4 side 4 has 4096 x 1536 edge-square maps
    let c = CubicalTorus::new(&[4, 4, 4, 4]).chain_complex();
    assert_eq!(c.betti_mod_p(2), vec![1, 4, 6, 4, 1]);
}
