//! Exact chiral 2x2 and 4x4 matrices over Q(i).
//! eta=(-,+,+,+), gamma={{0,sigma},{barsigma,0}}, {gamma,gamma}=-2eta.
use crate::scalar::{C, Q};
pub type Matrix<const N: usize> = [[C; N]; N];
pub fn zero<const N: usize>() -> Matrix<N> {
    [[C::ZERO; N]; N]
}
pub fn identity<const N: usize>() -> Matrix<N> {
    std::array::from_fn(|a| std::array::from_fn(|b| if a == b { C::ONE } else { C::ZERO }))
}
pub fn add<const N: usize>(a: Matrix<N>, b: Matrix<N>) -> Matrix<N> {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] + b[i][j]))
}
pub fn scale<const N: usize>(a: Matrix<N>, c: C) -> Matrix<N> {
    a.map(|r| r.map(|v| v * c))
}
pub fn mul<const N: usize>(a: Matrix<N>, b: Matrix<N>) -> Matrix<N> {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..N).fold(C::ZERO, |s, k| s + a[i][k] * b[k][j]))
    })
}
pub fn trace<const N: usize>(a: Matrix<N>) -> C {
    (0..N).fold(C::ZERO, |s, i| s + a[i][i])
}
pub const ETA: [i64; 4] = [-1, 1, 1, 1];
pub fn sigma(mu: usize) -> Matrix<2> {
    match mu {
        0 => identity(),
        1 => [[C::ZERO, C::ONE], [C::ONE, C::ZERO]],
        2 => [[C::ZERO, -C::I], [C::I, C::ZERO]],
        3 => [[C::ONE, C::ZERO], [C::ZERO, (-1).into()]],
        _ => panic!("Lorentz component outside 0..4"),
    }
}
pub fn barsigma(mu: usize) -> Matrix<2> {
    scale(sigma(mu), if mu == 0 { C::ONE } else { (-1).into() })
}
pub fn gamma(mu: usize) -> Matrix<4> {
    let (s, b) = (sigma(mu), barsigma(mu));
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i < 2 && j >= 2 {
                s[i][j - 2]
            } else if i >= 2 && j < 2 {
                b[i - 2][j]
            } else {
                C::ZERO
            }
        })
    })
}
pub fn epsilon_upper() -> Matrix<2> {
    [[C::ZERO, C::ONE], [(-1).into(), C::ZERO]]
}
pub fn epsilon_lower() -> Matrix<2> {
    scale(epsilon_upper(), (-1).into())
}
pub fn left_generator(mu: usize, nu: usize) -> Matrix<2> {
    scale(
        add(
            mul(sigma(mu), barsigma(nu)),
            scale(mul(sigma(nu), barsigma(mu)), (-1).into()),
        ),
        C::new(Q::ZERO, Q::new(1, 4)),
    )
}
