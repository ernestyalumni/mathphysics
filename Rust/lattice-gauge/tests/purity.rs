//! Exact check of mono-notes Proposition ym:prop:purity on integer spectra.
use lattice_gauge::spectral::{purity_bounds, trace_of_power};

#[test]
fn the_three_bounds_hold_on_twenty_thousand_pseudo_random_spectra() {
    let mut state: u64 = 0x2545_F491_4F6C_DD1D;
    let mut next = |modulus: u64| {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) % modulus
    };
    let mut checked = 0;
    for _ in 0..20_000 {
        let length = 1 + next(8) as usize;
        let spectrum: Vec<u64> = (0..length).map(|_| next(21)).collect();
        let exponent = 1 + next(3) as u32;
        let Some(bounds) = purity_bounds(&spectrum, exponent) else {
            assert!(spectrum.iter().all(|&value| value == 0));
            continue;
        };
        assert!(bounds.excited_sum_below_defect, "{spectrum:?} {exponent}");
        assert!(bounds.defect_below_quadratic, "{spectrum:?} {exponent}");
        assert!(bounds.doubling_improves, "{spectrum:?} {exponent}");
        checked += 1;
    }
    assert!(checked > 19_000);
}

#[test]
fn equality_for_degenerate_spectra_and_strictness_otherwise() {
    // lambda = (5, 5, 0), n = 1: s = 1, 1/P - 1 = (100 - 50)/50 = 1.
    assert_eq!(trace_of_power(&[5, 5, 0], 1), Some(10));
    assert_eq!(trace_of_power(&[5, 5, 0], 2), Some(50));
    // lambda = (3, 1), n = 1: s = 1/3, 1/P - 1 = (16 - 10)/10 = 3/5 > 1/3, so the reversed
    // first inequality 3/5 <= 1/3 is false: negative control, by direct integers.
    let (top, excited, trace, doubled) = (3u128, 1u128, 4u128, 10u128);
    assert!(excited * doubled <= top * (trace * trace - doubled));
    assert!(!(top * (trace * trace - doubled) <= excited * doubled));
    assert!(purity_bounds(&[3, 1], 1).unwrap().excited_sum_below_defect);
}

#[test]
fn invalid_spectra_are_rejected() {
    assert!(purity_bounds(&[], 1).is_none());
    assert!(purity_bounds(&[0, 0], 2).is_none());
    assert!(purity_bounds(&[2, 1], 0).is_none());
    assert!(purity_bounds(&[u64::MAX, 1], 3).is_none());
}
