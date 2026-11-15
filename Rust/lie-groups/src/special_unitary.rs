//! The Lie algebra su(N) with exact integer arithmetic: Killing form and Casimir constant.
//!
//! Conventions. Elements of su(N) are skew-Hermitian traceless N x N matrices. The invariant
//! inner product is `<X, Y> = -(1/N) Re tr(XY)` (the normalization of the Wilson lattice
//! action `1 - (1/N) Re tr U`). The Killing form is `K(X, Y) = tr_{su(N)}(ad X ad Y)`, a
//! real trace over the real vector space su(N). The Casimir constant `C` is defined by
//! `K = -C <.,.>`. All matrices here have Gaussian-integer entries, so every quantity is an
//! exact integer or an exact fraction of integers.

/// A Gaussian integer `real + i * imaginary`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct GaussianInteger {
    pub real: i64,
    pub imaginary: i64,
}

impl GaussianInteger {
    pub const ZERO: Self = Self { real: 0, imaginary: 0 };
    pub fn new(real: i64, imaginary: i64) -> Self {
        Self { real, imaginary }
    }
    fn add(self, other: Self) -> Self {
        Self::new(self.real + other.real, self.imaginary + other.imaginary)
    }
    fn subtract(self, other: Self) -> Self {
        Self::new(self.real - other.real, self.imaginary - other.imaginary)
    }
    fn multiply(self, other: Self) -> Self {
        Self::new(
            self.real * other.real - self.imaginary * other.imaginary,
            self.real * other.imaginary + self.imaginary * other.real,
        )
    }
}

/// A square matrix of Gaussian integers, row-major.
pub type GaussianMatrix = Vec<Vec<GaussianInteger>>;

/// An exact fraction `numerator / denominator` in lowest terms with positive denominator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExactFraction {
    pub numerator: i128,
    pub denominator: i128,
}

fn greatest_common_divisor(left: i128, right: i128) -> i128 {
    let (mut a, mut b) = (left.abs(), right.abs());
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

impl ExactFraction {
    /// `None` when the denominator is zero.
    pub fn new(numerator: i128, denominator: i128) -> Option<Self> {
        if denominator == 0 {
            return None;
        }
        let divisor = greatest_common_divisor(numerator, denominator).max(1);
        let sign = if denominator < 0 { -1 } else { 1 };
        Some(Self {
            numerator: sign * numerator / divisor,
            denominator: sign * denominator / divisor,
        })
    }
    pub fn integer(value: i128) -> Self {
        Self { numerator: value, denominator: 1 }
    }
    fn add(self, other: Self) -> Self {
        Self::new(
            self.numerator * other.denominator + other.numerator * self.denominator,
            self.denominator * other.denominator,
        )
        .expect("product of nonzero denominators is nonzero")
    }
}

fn zero_matrix(size: usize) -> GaussianMatrix {
    vec![vec![GaussianInteger::ZERO; size]; size]
}

pub fn matrix_product(left: &GaussianMatrix, right: &GaussianMatrix) -> GaussianMatrix {
    let size = left.len();
    let mut product = zero_matrix(size);
    for row in 0..size {
        for column in 0..size {
            let mut entry = GaussianInteger::ZERO;
            for inner in 0..size {
                entry = entry.add(left[row][inner].multiply(right[inner][column]));
            }
            product[row][column] = entry;
        }
    }
    product
}

/// `[X, Y] = XY - YX`.
pub fn commutator(left: &GaussianMatrix, right: &GaussianMatrix) -> GaussianMatrix {
    let forward = matrix_product(left, right);
    let backward = matrix_product(right, left);
    let size = left.len();
    let mut result = zero_matrix(size);
    for row in 0..size {
        for column in 0..size {
            result[row][column] = forward[row][column].subtract(backward[row][column]);
        }
    }
    result
}

/// `-Re tr(XY)`: N times the inner product `<X, Y>`. Positive definite on su(N).
pub fn negative_real_trace_form(left: &GaussianMatrix, right: &GaussianMatrix) -> i64 {
    let size = left.len();
    let mut total = 0;
    for row in 0..size {
        for column in 0..size {
            total += left[row][column].multiply(right[column][row]).real;
        }
    }
    -total
}

/// A basis of su(N), orthogonal for [`negative_real_trace_form`], with `N^2 - 1` elements:
/// `E_ij - E_ji` and `i (E_ij + E_ji)` for `i < j`, and `i * diag(1, ..., 1, -k, 0, ..., 0)`
/// with `k` ones, for `k = 1, ..., N-1`. Returns `None` for `N < 2`.
pub fn special_unitary_basis(size: usize) -> Option<Vec<GaussianMatrix>> {
    if size < 2 {
        return None;
    }
    let mut basis = Vec::with_capacity(size * size - 1);
    for row in 0..size {
        for column in (row + 1)..size {
            let mut antisymmetric = zero_matrix(size);
            antisymmetric[row][column] = GaussianInteger::new(1, 0);
            antisymmetric[column][row] = GaussianInteger::new(-1, 0);
            basis.push(antisymmetric);
            let mut symmetric = zero_matrix(size);
            symmetric[row][column] = GaussianInteger::new(0, 1);
            symmetric[column][row] = GaussianInteger::new(0, 1);
            basis.push(symmetric);
        }
    }
    for count in 1..size {
        let mut diagonal = zero_matrix(size);
        for index in 0..count {
            diagonal[index][index] = GaussianInteger::new(0, 1);
        }
        diagonal[count][count] = GaussianInteger::new(0, -(count as i64));
        basis.push(diagonal);
    }
    Some(basis)
}

/// The Killing form `K(X, Y) = tr_{su(N)}(ad X ad Y)`, computed as a real trace in an
/// orthogonal basis `e_i`: `sum_i B(e_i, [X, [Y, e_i]]) / B(e_i, e_i)` with
/// `B = -Re tr`. No closed formula is used.
pub fn killing_form(
    basis: &[GaussianMatrix],
    left: &GaussianMatrix,
    right: &GaussianMatrix,
) -> ExactFraction {
    let mut total = ExactFraction::integer(0);
    for element in basis {
        let image = commutator(left, &commutator(right, element));
        let numerator = i128::from(negative_real_trace_form(element, &image));
        let norm = i128::from(negative_real_trace_form(element, element));
        total = total.add(ExactFraction::new(numerator, norm).expect("basis elements are nonzero"));
    }
    total
}

/// The Casimir constant `C` of su(N) for `<X, Y> = -(1/N) Re tr(XY)`, computed from the
/// definition `K = -C <.,.>` on every pair of basis elements. Returns `None` if `N < 2` or if
/// the Killing form fails to be proportional to the inner product (which would be a bug).
pub fn special_unitary_casimir_constant(size: usize) -> Option<ExactFraction> {
    let basis = special_unitary_basis(size)?;
    let mut casimir: Option<ExactFraction> = None;
    for (left_index, left) in basis.iter().enumerate() {
        for (right_index, right) in basis.iter().enumerate() {
            let killing = killing_form(&basis, left, right);
            let trace_form = i128::from(negative_real_trace_form(left, right));
            if left_index != right_index {
                if trace_form != 0 || killing.numerator != 0 {
                    return None;
                }
                continue;
            }
            // K = -C * B / N  =>  C = -N K / B.
            let candidate = ExactFraction::new(
                -(size as i128) * killing.numerator,
                killing.denominator * trace_form,
            )?;
            match casimir {
                None => casimir = Some(candidate),
                Some(previous) if previous != candidate => return None,
                Some(_) => {}
            }
        }
    }
    casimir
}
