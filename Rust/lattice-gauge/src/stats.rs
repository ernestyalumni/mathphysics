//! Error analysis for Monte Carlo time series: jackknife over contiguous bins and the integrated
//! autocorrelation time with the automatic window of Madras and Sokal (J. Stat. Phys. 50 (1988) 109).

pub fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

/// Jackknife estimate and error of f(bin means) over `bins` contiguous bins of each series.
/// All series must have the same length; trailing samples that do not fill a bin are dropped.
pub fn jackknife<F: Fn(&[f64]) -> f64>(series: &[&[f64]], bins: usize, f: F) -> (f64, f64) {
    let n = series[0].len();
    let size = n / bins;
    assert!(bins >= 2 && size >= 1, "jackknife: too few samples");
    let used = size * bins;
    let totals: Vec<f64> = series.iter().map(|s| s[..used].iter().sum()).collect();
    let bin_sums: Vec<Vec<f64>> = series.iter().map(|s| (0..bins).map(|b| s[b * size..(b + 1) * size].iter().sum()).collect()).collect();
    let full: Vec<f64> = totals.iter().map(|t| t / used as f64).collect();
    let estimate = f(&full);
    let leave_out: Vec<f64> = (0..bins)
        .map(|b| {
            let means: Vec<f64> = totals.iter().zip(&bin_sums).map(|(t, bs)| (t - bs[b]) / (used - size) as f64).collect();
            f(&means)
        })
        .collect();
    let m = mean(&leave_out);
    let var = leave_out.iter().map(|v| (v - m) * (v - m)).sum::<f64>() * (bins as f64 - 1.0) / bins as f64;
    (estimate, var.sqrt())
}

/// Integrated autocorrelation time tau_int = 1/2 + sum_{t=1}^{W} rho(t), with W the smallest
/// window such that W >= c tau_int(W).
pub fn tau_int(x: &[f64], c: f64) -> f64 {
    let n = x.len();
    let m = mean(x);
    let c0: f64 = x.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / n as f64;
    if c0 == 0.0 {
        return 0.5;
    }
    let mut tau = 0.5;
    for t in 1..n / 2 {
        let ct: f64 = (0..n - t).map(|i| (x[i] - m) * (x[i + t] - m)).sum::<f64>() / (n - t) as f64;
        tau += ct / c0;
        if t as f64 >= c * tau {
            break;
        }
    }
    tau
}
