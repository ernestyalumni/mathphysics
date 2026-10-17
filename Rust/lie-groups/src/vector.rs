use crate::{Error, Real, Scalar};
pub fn dot<T: Scalar>(a: &[T; 3], b: &[T; 3]) -> T {
    (0..3).fold(T::zero(), |s, i| s + a[i].clone() * b[i].clone())
}
pub fn cross<T: Scalar>(a: &[T; 3], b: &[T; 3]) -> [T; 3] {
    [
        a[1].clone() * b[2].clone() - a[2].clone() * b[1].clone(),
        a[2].clone() * b[0].clone() - a[0].clone() * b[2].clone(),
        a[0].clone() * b[1].clone() - a[1].clone() * b[0].clone(),
    ]
}
pub(crate) fn norm<T: Real>(a: [T; 3]) -> Result<T, Error> {
    if a.iter().any(|v| !v.is_finite()) {
        return Err(Error::NonFinite);
    }
    let n = a[0].hypot(a[1]).hypot(a[2]);
    if !n.is_finite() {
        Err(Error::NonFinite)
    } else {
        Ok(n)
    }
}
pub(crate) fn normalized<T: Real>(a: [T; 3]) -> Result<[T; 3], Error> {
    let n = norm(a)?;
    if n == T::zero() {
        return Err(Error::Domain("axis must be nonzero"));
    }
    Ok(a.map(|x| x / n))
}
