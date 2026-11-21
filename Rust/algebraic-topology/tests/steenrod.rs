//! Oracles for cup-1 and the Pontryagin square: the coboundary identity of cup-1 on random integer cochains, and
//! invariance of int P(w) mod 4 on T^4 under w -> w + delta(lambda) mod 2 (with P(w) = 2 Pf(n) for constant twists).
use algebraic_topology::simplicial::{SimplicialComplex, TorusTriangulation};

struct Lcg(u64);
impl Lcg { fn next(&mut self) -> i64 { self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((self.0 >> 33) % 7) as i64 - 3 } }

fn add(a: &[i64], b: &[i64], s: i64) -> Vec<i64> { a.iter().zip(b).map(|(x, y)| x + s * y).collect() }

#[test]
fn cup1_coboundary_identity() {
    // a 4-simplex and its faces is enough for degrees p, q <= 2 (cochains on all faces)
    let c = SimplicialComplex::from_facets(&[vec![0, 1, 2, 3, 4], vec![2, 3, 4, 5, 6]]);
    let mut r = Lcg(7);
    for (p, q) in [(1usize, 1usize), (1, 2), (2, 1), (2, 2)] {
        for _ in 0..20 {
            let a: Vec<i64> = (0..c.simplices[p].len()).map(|_| r.next()).collect();
            let b: Vec<i64> = (0..c.simplices[q].len()).map(|_| r.next()).collect();
            let n = p + q - 1;
            let lhs = c.coboundary(n, &c.cup1(p, &a, q, &b));
            let s1 = if (p + q - 1) % 2 == 0 { 1 } else { -1 };
            let s2 = if (p + q + p * q) % 2 == 0 { 1 } else { -1 };
            let s3 = if p % 2 == 0 { 1 } else { -1 };
            let mut rhs = add(&vec![0; lhs.len()], &c.cup(p, &a, q, &b), s1);
            rhs = add(&rhs, &c.cup(q, &b, p, &a), s2);
            rhs = add(&rhs, &c.cup1(p + 1, &c.coboundary(p, &a), q, &b), 1);
            rhs = add(&rhs, &c.cup1(p, &a, q + 1, &c.coboundary(q, &b)), s3);
            assert_eq!(lhs, rhs, "coboundary identity fails for p={p} q={q}");
        }
    }
}

#[test]
fn pontryagin_square_is_a_class_invariant() {
    let t = TorusTriangulation::new(&[3, 3, 3, 3]);
    let c = &t.complex;
    let x: Vec<Vec<i64>> = (0..4).map(|mu| t.seam_cocycle(mu)).collect();
    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut r = Lcg(11);
    let mut plain_changed = 0;
    for n in [[1i64, 0, 0, 0, 0, 1], [1, 1, 1, 1, 1, 1], [0, 1, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0]] {
        let mut w = vec![0i64; c.simplices[2].len()];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            if n[k] != 0 { w = add(&w, &c.cup(1, &x[a], 1, &x[b]), n[k]); }
        }
        let p0 = t.evaluate_top(&c.pontryagin_square(&w)).rem_euclid(4);
        let pf = n[0] * n[5] - n[1] * n[4] + n[2] * n[3];
        assert_eq!(p0, (2 * pf).rem_euclid(4), "constant twist n={n:?}");
        for _ in 0..3 {
            // w' = w + delta(lambda) + 2 mu: another integer lift of a cohomologous mod-2 class
            let lambda: Vec<i64> = (0..c.simplices[1].len()).map(|_| r.next()).collect();
            let mu: Vec<i64> = (0..c.simplices[2].len()).map(|_| r.next()).collect();
            let w2 = add(&add(&w, &c.coboundary(1, &lambda), 1), &mu, 2);
            let p1 = t.evaluate_top(&c.pontryagin_square(&w2)).rem_euclid(4);
            assert_eq!(p1, p0, "Pontryagin square changed under w -> w + delta lambda + 2 mu (n={n:?})");
            let plain0 = t.evaluate_top(&c.cup(2, &w, 2, &w)).rem_euclid(4);
            let plain1 = t.evaluate_top(&c.cup(2, &w2, 2, &w2)).rem_euclid(4);
            if plain0 != plain1 { plain_changed += 1; }
        }
    }
    // the cup square alone is not a class invariant for non-cocycle lifts; cup_1 must be doing work
    assert!(plain_changed > 0, "test is vacuous: the plain cup square never changed");
}
