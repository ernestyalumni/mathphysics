//! Oracles: Random123 known-answer vectors for philox4x32 with 10 rounds (kat_vectors in the
//! Random123 distribution), and uniformity of the f64 stream.
use lattice_gauge::philox::{philox4x32, Stream};

#[test]
fn known_answer_vectors() {
    assert_eq!(philox4x32([0, 0, 0, 0], [0, 0]), [0x6627e8d5, 0xe169c58d, 0xbc57ac4c, 0x9b00dbd8]);
    assert_eq!(philox4x32([0xffffffff; 4], [0xffffffff; 2]), [0x408f276d, 0x41c83b0e, 0xa20bc7c6, 0x6d5451fd]);
    assert_eq!(
        philox4x32([0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344], [0xa4093822, 0x299f31d0]),
        [0xd16cfe09, 0x94fdcceb, 0x5001e420, 0x24126ea1]
    );
}

#[test]
fn uniform_moments() {
    let mut s = Stream::new(12345, 7, 99, 3);
    let n = 1_000_000;
    let xs: Vec<f64> = (0..n).map(|_| s.next_f64_open01()).collect();
    assert!(xs.iter().all(|&x| x > 0.0 && x < 1.0));
    let m = xs.iter().sum::<f64>() / n as f64;
    let v = xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64;
    let se_mean = (1.0f64 / 12.0 / n as f64).sqrt();
    assert!((m - 0.5).abs() < 6.0 * se_mean, "mean {m}");
    assert!((v - 1.0 / 12.0).abs() < 6.0 * (1.0 / 180.0 / n as f64).sqrt(), "variance {v}");
}

#[test]
fn streams_differ_and_repeat() {
    let a: Vec<u32> = { let mut s = Stream::new(1, 2, 3, 4); (0..8).map(|_| s.next_u32()).collect() };
    let b: Vec<u32> = { let mut s = Stream::new(1, 2, 3, 4); (0..8).map(|_| s.next_u32()).collect() };
    let c: Vec<u32> = { let mut s = Stream::new(1, 2, 4, 4); (0..8).map(|_| s.next_u32()).collect() };
    assert_eq!(a, b);
    assert_ne!(a, c);
}
