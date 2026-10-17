use lie_groups::{
    su2,
    symbolic::{self, Chart, SymbolicSu2},
    Half, Matrix, Real, So3Algebra, Su2Algebra, Su2Basis,
};
use qft_symbolics::{polynomial::Poly, scalar::Q};
fn numerical<T: Real>() {
    let algebra = So3Algebra::new([
        T::from_f64(0.),
        T::from_f64(0.),
        T::from_f64(std::f64::consts::FRAC_PI_2),
    ]);
    println!(
        "{}: Rz(pi/2) = {:?}",
        std::any::type_name::<T>(),
        algebra.exp().unwrap().matrix()
    );
    let spinor = Su2Algebra::new([
        T::from_f64(0.),
        T::from_f64(0.),
        T::from_f64(std::f64::consts::TAU),
    ])
    .exp()
    .unwrap();
    println!("spinor after 2pi (approximately -I): {:?}", spinor.matrix());
}
fn main() {
    let _: Matrix<f32, 3, 3> = Matrix::identity();
    numerical::<f32>();
    numerical::<f64>();
    numerical::<Half>();
    let e1 = So3Algebra::<Q>::basis(0).unwrap();
    let e2 = So3Algebra::<Q>::basis(1).unwrap();
    println!("[E1,E2] = {:?}", e1.matrix().commutator(&e2.matrix()));
    println!(
        "T3 = {:?}",
        su2::generator::<Q>(2, Su2Basis::HalfPauliOverI).unwrap()
    );
    let axis = [Q::ZERO, Q::ZERO, Q::ONE];
    println!(
        "symbolic exp(theta E3), s=sin(theta), c=cos(theta): {:?}",
        symbolic::so3_exponential(axis).unwrap()
    );
    let u = SymbolicSu2::stereographic(
        std::array::from_fn(|i| Poly::var(format!("t{i}"))),
        Chart::Plus,
    )
    .unwrap();
    println!(
        "symbolic real-chart unit norm residual: {}",
        u.norm_residual()
    );
}
