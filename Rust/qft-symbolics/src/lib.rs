//! Exact, deliberately bounded QFT calculations; not a Cadabra-compatible CAS.
//! Legacy operators remain bounded; checked APIs return typed errors.
//! See README for domains and conventions; computed evidence is not a formal proof.
pub mod clifford;
pub mod lorentz;
pub mod matrices;
pub mod polynomial;
pub mod scalar;
pub mod smga;
pub mod spinor;
pub mod tensor;
pub mod tensor_cases;
pub mod twistor;

pub mod coefficient;
/// Fallible algebra and explicit deterministic limits.
pub mod computation;
/// Opt-in arbitrary-precision rational and Gaussian-rational coefficients.
pub mod exact;

pub mod checked;
pub mod conventions;
pub mod differential;
pub mod expression;
pub mod ordered;
pub mod relativity;
pub mod replay;
