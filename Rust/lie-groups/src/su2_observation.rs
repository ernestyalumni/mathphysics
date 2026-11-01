//! Pure numerical SU(2) observations in the plus-i Pauli convention.
//!
//! Accepted axes retain their supplied values (norm error at most 1e-12).
//! Consequently observations report residuals rather than claiming exact group
//! membership. The existing strict `Su2` constructor and algebra exponential
//! would change this domain or the supplied phase, so neither is forced here.
//! Matrix assembly reuses `su2::quaternion_matrix` with an explicit sign map.

pub type Complex = [f64; 2];
fn mul([a, b]: Complex, [c, d]: Complex) -> Complex {
    [a * c - b * d, a * d + b * c]
}
fn add([a, b]: Complex, [c, d]: Complex) -> Complex {
    [a + c, b + d]
}
fn conj([a, b]: Complex) -> Complex {
    [a, -b]
}
fn magnitude([a, b]: Complex) -> f64 {
    a.hypot(b)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    FiniteNumberRequired,
    UnitAxisRequired,
    TOutOfRange,
}
impl Error {
    pub const fn code(self) -> &'static str {
        match self {
            Self::FiniteNumberRequired => "finite_number_required",
            Self::UnitAxisRequired => "unit_axis_required",
            Self::TOutOfRange => "t_out_of_range",
        }
    }
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.code())
    }
}
impl std::error::Error for Error {}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chart {
    Finite([f64; 3]),
    Pole,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Slice {
    Sphere { coordinate: f64, radius: f64 },
    Point { coordinate: f64 },
    Empty { coordinate: f64 },
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Checks {
    pub unitarity_error: f64,
    pub determinant_error: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Observation {
    pub axis: [f64; 3],
    pub t: f64,
    pub coordinates: [f64; 4],
    pub alpha: Complex,
    pub beta: Complex,
    pub matrix: [[Complex; 2]; 2],
    pub bloch: [f64; 3],
    pub stereographic: Chart,
    pub slice: Slice,
    pub checks: Checks,
}
pub fn validate_axis(axis: [f64; 3]) -> Result<(), Error> {
    if axis.iter().any(|x| !x.is_finite()) {
        return Err(Error::FiniteNumberRequired);
    }
    let norm = axis[0].hypot(axis[1]).hypot(axis[2]);
    if !norm.is_finite() || (norm - 1.0).abs() > 1e-12 {
        return Err(Error::UnitAxisRequired);
    }
    Ok(())
}
pub fn validate_t(t: f64) -> Result<(), Error> {
    if !t.is_finite() {
        return Err(Error::FiniteNumberRequired);
    }
    if t.abs() > 1e6 {
        return Err(Error::TOutOfRange);
    }
    Ok(())
}
/// Stable near the north pole; Pole also means the chart exceeds f64 range.
fn chart(axis: [f64; 3], t: f64) -> Chart {
    let (s, c) = (t / 4.0).sin_cos();
    let cot = c / s;
    let point = axis.map(|n| -n * cot);
    if t % (4.0 * std::f64::consts::PI) == 0.0 || point.iter().any(|x| !x.is_finite()) {
        Chart::Pole
    } else {
        Chart::Finite(point)
    }
}
/// Observe the axis/angle Pauli formula, without normalizing the supplied axis.
/// For an exactly unit axis this is U=exp(-it axis.sigma/2).
pub fn observe(axis: [f64; 3], t: f64, c: f64) -> Result<Observation, Error> {
    validate_axis(axis)?;
    validate_t(t)?;
    if !c.is_finite() {
        return Err(Error::FiniteNumberRequired);
    }
    let (s, a0) = (t / 2.0).sin_cos();
    let [a1, a2, a3] = axis.map(|n| -n * s);
    let alpha = [a0, a3];
    let beta = [a2, a1];
    // Existing helper uses U=q0 I-i q.sigma; negate the plus-i vector.
    let q_matrix = crate::su2::quaternion_matrix(&[a0, -a1, -a2, -a3]);
    let matrix = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let z = q_matrix[(i, j)];
            [z.re, z.im]
        })
    });
    let mut unitarity_error: f64 = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let mut entry = add(
                mul(conj(matrix[0][i]), matrix[0][j]),
                mul(conj(matrix[1][i]), matrix[1][j]),
            );
            entry[0] -= if i == j { 1.0 } else { 0.0 };
            unitarity_error = unitarity_error.max(magnitude(entry));
        }
    }
    let diagonal = mul(matrix[0][0], matrix[1][1]);
    let off_diagonal = mul(matrix[0][1], matrix[1][0]);
    let determinant_error = magnitude([
        diagonal[0] - off_diagonal[0] - 1.0,
        diagonal[1] - off_diagonal[1],
    ]);

    let stereographic = chart(axis, t);
    let slice = if c.abs() > 1.0 {
        Slice::Empty { coordinate: c }
    } else if c.abs() == 1.0 {
        Slice::Point { coordinate: c }
    } else {
        Slice::Sphere {
            coordinate: c,
            radius: ((1.0 - c) * (1.0 + c)).sqrt(),
        }
    };
    // Rodrigues action for U=exp(-it n.sigma/2), on initial +x.
    let bloch = [
        a0 * a0 + a1 * a1 - a2 * a2 - a3 * a3,
        2.0 * (a1 * a2 - a0 * a3),
        2.0 * (a1 * a3 + a0 * a2),
    ];

    Ok(Observation {
        axis,
        t,
        coordinates: [a0, a1, a2, a3],
        alpha,
        beta,
        matrix,
        bloch,
        stereographic,
        slice,
        checks: Checks {
            unitarity_error,
            determinant_error,
        },
    })
}
