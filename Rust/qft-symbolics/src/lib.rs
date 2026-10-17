//! Exact, deliberately bounded QFT calculations; not a Cadabra-compatible CAS.
//! Arithmetic overflow and expansion-budget exhaustion panic (never a proof).
//! The CLI catches these as INCONCLUSIVE. See README for domains and conventions.
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
