//! Oracles: quaternion group laws; Haar moments on S^3 (E a0^2 = 1/4, E a0^4 = 1/8); the heat-bath
//! marginal for density sqrt(1-x^2) e^{a x}: E x = I_2(a)/I_1(a); over-relaxation invariance.
use lattice_gauge::character::bessel_i_scaled;
use lattice_gauge::philox::Stream;
use lattice_gauge::su2::*;

fn close(a: Quaternion, b: Quaternion, tol: f64) -> bool {
    a.a.iter().zip(b.a.iter()).all(|(x, y)| (x - y).abs() < tol)
}

#[test]
fn group_laws() {
    let mut s = Stream::new(1, 0, 0, 0);
    for _ in 0..1000 {
        let (u, v, w) = (haar(&mut s), haar(&mut s), haar(&mut s));
        assert!((u.norm_sq() - 1.0).abs() < 1e-14);
        assert!(close(u.mul(v).mul(w), u.mul(v.mul(w)), 1e-14));
        assert!(close(u.mul(u.conj()), Quaternion::IDENTITY, 1e-14));
        assert!(close(u.mul(v).conj(), v.conj().mul(u.conj()), 1e-14));
        assert!((u.mul(v).half_trace() - v.mul(u).half_trace()).abs() < 1e-14);
    }
    // Pauli algebra: (i s1)(i s2) = -(i s3)
    let i1 = Quaternion::new(0.0, 1.0, 0.0, 0.0);
    let i2 = Quaternion::new(0.0, 0.0, 1.0, 0.0);
    assert!(close(i1.mul(i2), Quaternion::new(0.0, 0.0, 0.0, -1.0), 1e-15));
}

#[test]
fn haar_moments() {
    let mut s = Stream::new(2, 0, 0, 0);
    let n = 400_000;
    let xs: Vec<f64> = (0..n).map(|_| haar(&mut s).half_trace()).collect();
    let m2 = xs.iter().map(|x| x * x).sum::<f64>() / n as f64;
    let m4 = xs.iter().map(|x| x.powi(4)).sum::<f64>() / n as f64;
    // sd(a0^2) = sqrt(1/8 - 1/16), sd(a0^4) = sqrt(E a0^8 - 1/64) with E a0^8 = 105/(4*6*8*10)... = 7/128
    assert!((m2 - 0.25).abs() < 6.0 * (1.0f64 / 16.0 / n as f64).sqrt(), "E a0^2 = {m2}");
    assert!((m4 - 0.125).abs() < 6.0 * ((7.0 / 128.0 - 1.0 / 64.0) / n as f64).sqrt(), "E a0^4 = {m4}");
}

fn marginal(a: f64, broken: bool) -> (f64, f64) {
    let mut s = Stream::new(3, 0, (a * 1000.0) as u32, 0);
    let staple = haar(&mut s).scale(a); // beta = 1, k = a
    let n = 300_000;
    let xs: Vec<f64> = (0..n)
        .map(|_| {
            let u = if broken { heat_bath_broken(staple, 1.0, &mut s) } else { heat_bath(staple, 1.0, &mut s) };
            u.mul(staple).half_trace() / a
        })
        .collect();
    let m = xs.iter().sum::<f64>() / n as f64;
    let sd = (xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64).sqrt();
    (m, sd / (n as f64).sqrt())
}

#[test]
fn heat_bath_marginal() {
    for &a in &[0.5, 2.0, 5.0] {
        let want = bessel_i_scaled(2, a) / bessel_i_scaled(1, a);
        let (m, se) = marginal(a, false);
        assert!((m - want).abs() < 6.0 * se, "a={a}: {m} vs {want} (se {se})");
    }
}

#[test]
fn broken_heat_bath_is_detected() {
    let a = 2.0;
    let want = bessel_i_scaled(2, a) / bessel_i_scaled(1, a);
    let (m, se) = marginal(a, true);
    assert!((m - want).abs() > 10.0 * se, "negative control not detected: {m} vs {want}");
}

#[test]
fn overrelaxation() {
    let mut s = Stream::new(4, 0, 0, 0);
    for _ in 0..1000 {
        let u = haar(&mut s);
        let w = haar(&mut s).scale(1.7).add(haar(&mut s).scale(0.4));
        let v = overrelax(u, w);
        assert!((v.mul(w).half_trace() - u.mul(w).half_trace()).abs() < 1e-13);
        assert!((v.norm_sq() - 1.0).abs() < 1e-13);
        assert!(close(overrelax(v, w), u, 1e-13));
    }
}
