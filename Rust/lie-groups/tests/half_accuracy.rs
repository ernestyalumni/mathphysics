use lie_groups::{Half, So3Algebra};
#[test]
fn binary16_rotations_agree_with_f64_for_same_represented_inputs() {
    let mut maximum = 0.0_f64;
    for i in 0..128 {
        let v = [
            0.031 * (i as f64 - 64.),
            0.017 * ((i * i % 31) as f64 - 15.),
            0.023 * ((i * 7 % 41) as f64 - 20.),
        ]
        .map(Half::from_f64);
        let half = So3Algebra::new(v).exp().unwrap();
        let reference = So3Algebra::new(v.map(Half::to_f64)).exp().unwrap();
        for a in 0..3 {
            for b in 0..3 {
                maximum = maximum
                    .max((half.matrix()[(a, b)].to_f64() - reference.matrix()[(a, b)]).abs());
            }
        }
    }
    println!("binary16 128-fixture maximum matrix component error vs f64: {maximum}");
    assert!(maximum < 0.01, "binary16 accuracy regression: {maximum}");
}
