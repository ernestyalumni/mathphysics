use lie_groups::Half;
#[test]
fn binary16_all_finite_encodings_round_trip_and_nan_stays_nan() {
    assert_eq!(std::mem::size_of::<Half>(), 2);
    for bits in 0..=u16::MAX {
        let h = Half::from_bits(bits);
        let value = h.to_f64();
        if value.is_nan() {
            assert!(Half::from_f64(value).to_f64().is_nan());
        } else {
            assert_eq!(Half::from_f64(value).to_bits(), bits, "bits={bits:04x}");
        }
    }
}
#[test]
fn binary16_known_extremes_rounding_ties_subnormals_and_signed_zero() {
    for (value, bits) in [
        (1., 0x3c00),
        (-2., 0xc000),
        (65504., 0x7bff),
        (65520., 0x7c00),
        (2_f64.powi(-14), 0x0400),
        (2_f64.powi(-24), 1),
        (2_f64.powi(-25), 0),
        (-0., 0x8000),
        (f64::INFINITY, 0x7c00),
    ] {
        assert_eq!(Half::from_f64(value).to_bits(), bits, "{value}");
    }
    assert_eq!(Half::from_f64(1. + 2_f64.powi(-11)).to_bits(), 0x3c00);
    assert_eq!(Half::from_f64(1. + 3. * 2_f64.powi(-11)).to_bits(), 0x3c02);
    // Values just above/below a half-way point must not double-round via f32.
    assert_eq!(
        Half::from_f64(1. + 2_f64.powi(-11) + 2_f64.powi(-40)).to_bits(),
        0x3c01
    );
    assert_eq!(
        Half::from_f64(1. + 2_f64.powi(-11) - 2_f64.powi(-40)).to_bits(),
        0x3c00
    );
    assert_eq!((Half::from_f64(1.) / Half::from_f64(2.)).to_bits(), 0x3800);
    assert!(Half::from_f64(f64::NAN)
        .partial_cmp(&Half::from_f64(0.))
        .is_none());
}
