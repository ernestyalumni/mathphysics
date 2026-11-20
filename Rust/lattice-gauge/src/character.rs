//! Modified Bessel functions and the SU(2) character expansion of the Wilson weight.
//!
//! Convention: action S = beta sum_p (1 - a0(U_p)) with a0 = (1/2) tr U, so the plaquette weight is
//! exp(beta a0(U)) up to a constant. Its character expansion is
//!     exp(beta a0(U)) = sum_j c_j(beta) chi_j(U),   c_j = (2j+1) 2 I_{2j+1}(beta) / beta,
//! j = 0, 1/2, 1, ... . On a closed surface of Euler characteristic chi with V plaquettes,
//! integrating every link gives Z = sum_j (2j+1)^chi (c_j/(2j+1))^V; for the torus chi = 0.

/// e^{-x} I_n(x) for integer n >= 0 and 0 <= x <= 600, by the power series
/// sum_k (x/2)^{2k+n} / (k! (k+n)!), whose terms are all positive (no cancellation).
pub fn bessel_i_scaled(n: u32, x: f64) -> f64 {
    assert!(x >= 0.0 && x <= 600.0, "bessel_i_scaled: x out of range");
    if x == 0.0 {
        return if n == 0 { 1.0 } else { 0.0 };
    }
    let h = 0.5 * x;
    // log of the first term, minus x for the scaling
    let mut log_lgamma_n1 = 0.0;
    for m in 2..=n {
        log_lgamma_n1 += (m as f64).ln();
    }
    let mut term = (n as f64 * h.ln() - log_lgamma_n1 - x).exp();
    let mut sum = term;
    let mut k = 0.0f64;
    let q = h * h;
    loop {
        k += 1.0;
        term *= q / (k * (k + n as f64));
        sum += term;
        if term < 1e-18 * sum && k > h {
            break;
        }
        if k > 5000.0 {
            break;
        }
    }
    sum
}

/// I_n(x).
pub fn bessel_i(n: u32, x: f64) -> f64 {
    bessel_i_scaled(n, x) * x.exp()
}

/// Character coefficient c_j(beta) for j = two_j / 2.
pub fn character_coefficient(two_j: u32, beta: f64) -> f64 {
    (two_j + 1) as f64 * 2.0 * bessel_i(two_j + 1, beta) / beta
}

/// Normalized coefficient u_j = c_j / ((2j+1) c_0) = I_{2j+1}(beta) / I_1(beta).
pub fn normalized_coefficient(two_j: u32, beta: f64) -> f64 {
    bessel_i_scaled(two_j + 1, beta) / bessel_i_scaled(1, beta)
}

/// Infinite-volume 2D plaquette average <a0(U_p)> = I_2(beta)/I_1(beta).
pub fn plaquette_2d_infinite(beta: f64) -> f64 {
    normalized_coefficient(1, beta)
}

/// Infinite-volume 2D Wilson loop of area A (in plaquettes): (I_2/I_1)^A.
pub fn wilson_loop_2d_infinite(beta: f64, area: u32) -> f64 {
    plaquette_2d_infinite(beta).powi(area as i32)
}

/// Exact plaquette average on a 2D periodic torus with V plaquettes:
/// Z = sum_n b_n^V with b_n = 2 I_n(beta)/beta, n = 2j+1, and <a0> = (1/V) d log Z / d beta
///   = sum_n w_n (b_n'/b_n) / sum_n w_n, w_n = (b_n/b_1)^V = (I_n/I_1)^V,
/// b_n'/b_n = I_n'/I_n - 1/beta = I_{n+1}/I_n + (n-1)/beta (from I_{n-1} - I_{n+1} = (2n/beta) I_n).
pub fn torus_plaquette_2d(beta: f64, plaquettes: u64) -> f64 {
    let i1 = bessel_i_scaled(1, beta);
    let mut num = 0.0;
    let mut den = 0.0;
    let mut n = 1u32;
    loop {
        let i_n = bessel_i_scaled(n, beta);
        let w = (i_n / i1).powf(plaquettes as f64);
        if i_n == 0.0 || (n > 1 && w < 1e-18 * den) {
            break;
        }
        let ratio = bessel_i_scaled(n + 1, beta) / i_n + (n as f64 - 1.0) / beta;
        num += w * ratio;
        den += w;
        n += 1;
        if n > 2000 {
            break;
        }
    }
    num / den
}
