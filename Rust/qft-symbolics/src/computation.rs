//! Shared fallible computation and deterministic work limits.
//! Counters bound declared algebraic work, not elapsed time or OS allocation.
use std::fmt;
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Resource {
    Operations,
    Terms,
    Degree,
    CoefficientBits,
    InputBytes,
    WordLength,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    ZeroDenominator,
    DivisionByZero,
    InvalidInteger,
    Overflow {
        operation: &'static str,
    },
    ResourceLimit {
        resource: Resource,
        limit: u64,
        required: u64,
    },
    MissingVariable {
        name: String,
    },
    UnsupportedDomain {
        message: &'static str,
    },
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroDenominator => write!(f, "zero rational denominator"),
            Self::DivisionByZero => write!(f, "division by zero"),
            Self::InvalidInteger => write!(f, "invalid integer"),
            Self::Overflow { operation } => write!(f, "arithmetic overflow in {operation}"),
            Self::ResourceLimit {
                resource,
                limit,
                required,
            } => write!(
                f,
                "{resource:?} budget exceeded: required {required}, limit {limit}"
            ),
            Self::MissingVariable { name } => write!(f, "missing variable {name}"),
            Self::UnsupportedDomain { message } => write!(f, "unsupported domain: {message}"),
        }
    }
}
impl std::error::Error for Error {}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    pub max_operations: u64,
    pub max_terms: usize,
    pub max_degree: usize,
    pub max_coefficient_bits: u64,
    pub max_input_bytes: usize,
    pub max_word_length: usize,
}
impl Default for Limits {
    fn default() -> Self {
        Self {
            max_operations: 1_000_000,
            max_terms: 8192,
            max_degree: 64,
            max_coefficient_bits: 65_536,
            max_input_bytes: 1_048_576,
            max_word_length: 12,
        }
    }
}
#[derive(Clone, Debug)]
pub struct Budget {
    limits: Limits,
    used: u64,
}
impl Default for Budget {
    fn default() -> Self {
        Self::new(Limits::default())
    }
}
impl Budget {
    pub fn new(limits: Limits) -> Self {
        Self { limits, used: 0 }
    }
    pub fn limits(&self) -> &Limits {
        &self.limits
    }
    pub fn used(&self) -> u64 {
        self.used
    }
    pub fn check(&self, resource: Resource, required: u64) -> Result<()> {
        let limit = match resource {
            Resource::Operations => self.limits.max_operations,
            Resource::Terms => self.limits.max_terms as u64,
            Resource::Degree => self.limits.max_degree as u64,
            Resource::CoefficientBits => self.limits.max_coefficient_bits,
            Resource::InputBytes => self.limits.max_input_bytes as u64,
            Resource::WordLength => self.limits.max_word_length as u64,
        };
        if required > limit {
            Err(Error::ResourceLimit {
                resource,
                limit,
                required,
            })
        } else {
            Ok(())
        }
    }
    pub fn charge(&mut self, amount: u64) -> Result<()> {
        let required = self.used.checked_add(amount).ok_or(Error::ResourceLimit {
            resource: Resource::Operations,
            limit: self.limits.max_operations,
            required: u64::MAX,
        })?;
        self.check(Resource::Operations, required)?;
        self.used = required;
        Ok(())
    }
}
