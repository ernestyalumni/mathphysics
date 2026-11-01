//! Repository-owned example consuming the public library with fresh inputs.
//! Build artifacts belong in an external CARGO_TARGET_DIR.
use qft_symbolics::{clifford, polynomial::Poly, spinor};
fn main() {
    println!("One shared Rust trace API; different caller-supplied inputs:");
    for labels in [
        vec![],
        vec!["probe_a", "probe_b"],
        vec!["p", "q", "r"],
        vec!["mu", "nu", "rho", "sigma"],
    ] {
        println!("  trace({labels:?}) = {}", clifford::trace(&labels));
    }
    let [a, b, c, d] = ["fresh_a", "fresh_b", "fresh_c", "fresh_d"];
    let expected = Poly::constant(4)
        * (clifford::metric(a, b) * clifford::metric(c, d)
            - clifford::metric(a, c) * clifford::metric(b, d)
            + clifford::metric(a, d) * clifford::metric(b, c));
    assert_eq!(clifford::trace(&[a, b, c, d]), expected);
    let spinors = ["beam_a", "beam_b", "beam_c", "beam_d"].map(spinor::named);
    let residual = spinor::schouten(&spinors[0], &spinors[1], &spinors[2], &spinors[3]);
    assert!(residual.is_zero());
    println!("  fresh-symbol four-gamma identity: checked");
    println!("  Schouten residual: {residual}");
    println!("Domain: formal Clifford h, anticommutator=2h; trace(I)=4; no implicit sums.");
    println!("Exact within the current bounded arithmetic; not a formal proof or general CAS.");
}
