//! Oracle: AR(1) process x_t = r x_{t-1} + e_t has rho(t) = r^t, tau_int = (1 + r) / (2 (1 - r)).
use lattice_gauge::philox::Stream;
use lattice_gauge::stats::*;

#[test]
fn ar1_tau_and_jackknife() {
    let r = 0.8;
    let mut s = Stream::new(5, 0, 0, 0);
    let n = 400_000;
    let mut x = 0.0;
    let xs: Vec<f64> = (0..n).map(|_| { x = r * x + s.next_gaussian(); x }).collect();
    let tau = tau_int(&xs, 6.0);
    let want = (1.0 + r) / (2.0 * (1.0 - r));
    assert!((tau - want).abs() < 0.05 * want, "tau {tau} vs {want}");
    // jackknife error of the mean ~ sqrt(2 tau var / n), var = 1/(1-r^2)
    let (m, e) = jackknife(&[&xs], 100, |v| v[0]);
    let expect = (2.0 * want / (1.0 - r * r) / n as f64).sqrt();
    assert!(m.abs() < 6.0 * expect);
    assert!((e / expect - 1.0).abs() < 0.3, "jackknife error {e} vs {expect}");
}
