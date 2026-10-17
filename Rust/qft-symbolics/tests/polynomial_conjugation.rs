use qft_symbolics::{polynomial::Poly, scalar::C};
use std::collections::BTreeMap;
#[test]
fn coefficient_conjugation_preserves_formal_real_symbol_contract() {
    let p = Poly::constant(C::I) * Poly::var("z") + Poly::constant(2);
    assert!(!p.has_real_coefficients());
    assert_eq!(p.conjugate_coefficients().conjugate_coefficients(), p);
    assert!(Poly::var("z").has_real_coefficients());
    let values = BTreeMap::from([("z".into(), C::from(3))]);
    assert_eq!(
        p.conjugate_coefficients().eval(&values).unwrap(),
        p.eval(&values).unwrap().conj()
    );
    let complex_values = BTreeMap::from([("z".into(), C::I)]);
    // It must not be advertised as full conjugation of arbitrary complex variables.
    assert_ne!(
        p.conjugate_coefficients().eval(&complex_values).unwrap(),
        p.eval(&complex_values).unwrap().conj()
    );
}
