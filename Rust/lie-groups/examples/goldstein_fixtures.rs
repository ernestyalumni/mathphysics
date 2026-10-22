//! Trusted-producer NDJSON fixtures for independent cross-language/CAS checks.
use lie_groups::{
    goldstein::{cayley_klein_matrix, euler_lift, euler_matrix, from_euler, EulerConvention},
    multitrig::MultiTrig,
    su2, Matrix, Scalar,
};
const CONVENTIONS: [(&str, EulerConvention); 3] = [
    ("zxz", EulerConvention::Zxz),
    ("zyz", EulerConvention::Zyz),
    ("zyx", EulerConvention::Zyx),
];
fn symbolic_matrix<T: std::fmt::Display, const N: usize>(m: &Matrix<T, N, N>) {
    print!("[");
    for i in 0..N {
        if i > 0 {
            print!(",");
        }
        print!("[");
        for j in 0..N {
            if j > 0 {
                print!(",");
            }
            print!("\"{}\"", m[(i, j)]);
        }
        print!("]");
    }
    print!("]");
}
fn symbolic() {
    type T = MultiTrig<3>;
    let s: [T; 3] = std::array::from_fn(|i| T::sin(i).unwrap());
    let c: [T; 3] = std::array::from_fn(|i| T::cos(i).unwrap());
    for (name, convention) in CONVENTIONS {
        let r = euler_matrix(convention, c.clone(), s.clone());
        print!("{{\"name\":\"{name}_full\",\"convention\":\"{name}\",\"argument\":\"full\",\"matrix\":");
        symbolic_matrix(&r);
        println!("}}");
        let unit = |axis: usize, i: usize| {
            let mut q = [c[i].clone(), T::zero(), T::zero(), T::zero()];
            q[axis + 1] = s[i].clone();
            q
        };
        let (a, b, d) = match convention {
            EulerConvention::Zxz => (unit(2, 0), unit(0, 1), unit(2, 2)),
            EulerConvention::Zyz => (unit(2, 0), unit(1, 1), unit(2, 2)),
            EulerConvention::Zyx => (unit(2, 0), unit(1, 1), unit(0, 2)),
        };
        let q = su2::hamilton(&su2::hamilton(&a, &b), &d);
        let r = su2::quaternion_rotation(&q);
        print!("{{\"name\":\"{name}_half\",\"convention\":\"{name}\",\"argument\":\"half\",\"matrix\":");
        symbolic_matrix(&r);
        print!(",\"quaternion\":[");
        for (i, x) in q.iter().enumerate() {
            if i > 0 {
                print!(",");
            }
            print!("\"{x}\"");
        }
        print!("],\"cayley_klein\":");
        symbolic_matrix(&cayley_klein_matrix(&q));
        println!("}}");
    }
}
fn main() {
    if std::env::args().any(|a| a == "--symbolic") {
        symbolic();
        return;
    }
    for (name, c) in CONVENTIONS {
        for (case, angles) in [
            ("regular", [0.31, 0.73, -1.2]),
            ("identity", [0.0, 0.0, 0.0]),
            ("small", [1e-12, -2e-12, 3e-12]),
            ("pi", [0.0, std::f64::consts::PI, 0.0]),
            ("pole_zero", [0.41, 0.0, -0.7]),
            ("pole_halfpi", [0.41, std::f64::consts::FRAC_PI_2, -0.7]),
            (
                "near_halfpi",
                [0.41, std::f64::consts::FRAC_PI_2 - 1e-10, -0.7],
            ),
            (
                "pole_negative_halfpi",
                [0.41, -std::f64::consts::FRAC_PI_2, -0.7],
            ),
        ] {
            let r = from_euler(c, angles).expect("finite rotation");
            let q = euler_lift(c, angles).expect("unit lift");
            let angle_order = if c == EulerConvention::Zyx {
                "psi,theta,phi"
            } else {
                "phi,theta,psi"
            };
            println!("{{\"name\":\"{name}_{case}\",\"convention\":\"{name}\",\"angles\":{angles:?},\"angle_order\":\"{angle_order}\",\"product_angles\":{angles:?},\"matrix\":{:?},\"quaternion\":{:?}}}",r.matrix().rows(),q.quaternion());
        }
    }
}
