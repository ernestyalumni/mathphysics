//! Oracles: reference digits of I_n (Abramowitz and Stegun table 9.8 / standard tables), the
//! recurrence I_{n-1} - I_{n+1} = (2n/x) I_n, the large-x asymptotic series, and limits of the
//! 2D formulas (small beta: I_2/I_1 -> beta/4; large torus -> infinite volume).
use lattice_gauge::character::*;

fn rel(a: f64, b: f64) -> f64 { ((a - b) / b).abs() }

#[test]
fn bessel_reference_values() {
    let table = [
        (0, 0.5, 1.0634833707413236), (1, 0.5, 0.2578943053908963),
        (0, 1.0, 1.2660658777520082), (1, 1.0, 0.5651591039924851), (2, 1.0, 0.1357476697670383),
        (5, 1.0, 2.714631559569719e-4),
        (0, 2.0, 2.2795853023360673), (1, 2.0, 1.5906368546373291), (2, 2.0, 0.6889484476987382),
        (0, 10.0, 2815.716628466254), (1, 10.0, 2670.988303701255), (2, 10.0, 2281.518967726004),
    ];
    for (n, x, want) in table {
        let got = bessel_i(n, x);
        assert!(rel(got, want) < 1e-13, "I_{n}({x}) = {got}, want {want}");
    }
}

#[test]
fn recurrence_and_asymptotics() {
    for &x in &[0.3, 1.0, 4.0, 17.0, 50.0, 200.0] {
        for n in 1..8u32 {
            // I_{n-1} = I_{n+1} + (2n/x) I_n: both terms positive, so no cancellation in the check.
            let lhs = bessel_i_scaled(n - 1, x);
            let rhs = bessel_i_scaled(n + 1, x) + 2.0 * n as f64 / x * bessel_i_scaled(n, x);
            assert!(rel(lhs, rhs) < 1e-13, "recurrence n={n} x={x}");
        }
    }
    for &x in &[50.0, 200.0] {
        // e^{-x} I_0(x) sqrt(2 pi x) = 1 + 1/(8x) + 9/(128x^2) + 225/(3072x^3) + ...
        let s = bessel_i_scaled(0, x) * (2.0 * std::f64::consts::PI * x).sqrt();
        let series = 1.0 + 1.0 / (8.0 * x) + 9.0 / (128.0 * x * x) + 225.0 / (3072.0 * x * x * x);
        assert!((s - series).abs() < 1e-6, "asymptotic x={x}: {s} vs {series}");
    }
}

#[test]
fn two_dimensional_limits() {
    for &b in &[1e-3, 1e-2] {
        assert!(rel(plaquette_2d_infinite(b), b / 4.0) < 1e-4);
    }
    for &b in &[0.5, 2.0, 5.0] {
        let inf = plaquette_2d_infinite(b);
        assert!((torus_plaquette_2d(b, 4096) - inf).abs() < 1e-12, "large torus beta={b}");
        // On a small torus the higher representations contribute: the result differs from infinite volume.
        assert!((torus_plaquette_2d(b, 4) - inf).abs() > 1e-6, "small torus should differ at beta={b}");
    }
}
