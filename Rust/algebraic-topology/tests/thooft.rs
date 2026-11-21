//! 't Hooft's twisted-torus formula from first principles: on a Freudenthal triangulation of T^4 the seam
//! cocycles x_mu are integer cocycles, x_0 x_1 x_2 x_3 evaluates to +-1 on the fundamental cycle, and for the
//! twist class n = sum_{mu<nu} n_{mu nu} x_mu x_nu, the integral of n cup n is 2 Pf(n) = 2 (n01 n23 - n02 n13 + n03 n12)
//! up to the global orientation sign. For SU(N) twists this gives the fractional instanton number -Pf(n)/N mod 1.
use algebraic_topology::simplicial::TorusTriangulation;

#[test]
fn twisted_torus_cup_square_is_twice_the_pfaffian() {
    let t = TorusTriangulation::new(&[3, 3, 3, 3]);
    let c = &t.complex;
    let x: Vec<Vec<i64>> = (0..4).map(|mu| t.seam_cocycle(mu)).collect();
    for mu in 0..4 { assert!(c.coboundary(1, &x[mu]).iter().all(|&v| v == 0), "x_{mu} is not a cocycle"); }
    let x01 = c.cup(1, &x[0], 1, &x[1]);
    let x012 = c.cup(2, &x01, 1, &x[2]);
    let vol = t.evaluate_top(&c.cup(3, &x012, 1, &x[3]));
    assert!(vol == 1 || vol == -1, "fundamental class normalization {vol}");
    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let cases: [[i64; 6]; 6] = [[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 1, 0], [0, 0, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1], [2, 0, 1, 3, 0, 5], [1, 0, 0, 0, 0, 0]];
    for n in cases {
        let mut cls = vec![0i64; c.simplices[2].len()];
        for (k, &(a, b)) in pairs.iter().enumerate() {
            if n[k] != 0 {
                let xab = c.cup(1, &x[a], 1, &x[b]);
                for (v, w) in cls.iter_mut().zip(xab) { *v += n[k] * w; }
            }
        }
        assert!(c.coboundary(2, &cls).iter().all(|&v| v == 0));
        let sq = t.evaluate_top(&c.cup(2, &cls, 2, &cls)) * vol;
        let pf = n[0] * n[5] - n[1] * n[4] + n[2] * n[3];
        assert_eq!(sq, 2 * pf, "n = {n:?}: integral of n cup n = {sq}, 2 Pf = {}", 2 * pf);
    }
}
