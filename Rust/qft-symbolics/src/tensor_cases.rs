//! Shared named regression cases for the JSON CLI and external Cadabra oracle.
use crate::{
    polynomial::Poly,
    tensor::{Context, Error, Expr, Factor as F, Index as I, Slot, Space, Symmetry, Term},
};
fn u(n: &str) -> I {
    I::up(Space::Lorentz, n)
}
fn d(n: &str) -> I {
    I::down(Space::Lorentz, n)
}
fn t(n: &str, indices: Vec<I>) -> F {
    F::tensor(n, indices)
}
fn term(c: i64, f: Vec<F>) -> Term {
    Term::new(Poly::constant(c), f)
}
pub struct Case {
    pub id: &'static str,
    pub residual: Expr,
    pub expected_zero: bool,
}
pub fn cases() -> Result<Vec<Case>, Error> {
    let mut c = Context::default();
    for name in ["A", "B", "V", "W"] {
        c.declare(name, vec![Slot::Lorentz], Symmetry::None)?;
    }
    c.declare("S", vec![Slot::Lorentz; 2], Symmetry::Symmetric)?;
    c.declare("K", vec![Slot::Lorentz; 2], Symmetry::Antisymmetric)?;
    let raw = vec![
        (
            "symmetric_swap",
            true,
            vec![
                term(1, vec![t("S", vec![d("a"), d("b")])]),
                term(-1, vec![t("S", vec![d("b"), d("a")])]),
            ],
        ),
        (
            "antisymmetric_swap",
            true,
            vec![
                term(1, vec![t("K", vec![d("a"), d("b")])]),
                term(1, vec![t("K", vec![d("b"), d("a")])]),
            ],
        ),
        (
            "antisymmetric_wrong_sign",
            false,
            vec![
                term(1, vec![t("K", vec![d("a"), d("b")])]),
                term(-1, vec![t("K", vec![d("b"), d("a")])]),
            ],
        ),
        (
            "symmetric_antisymmetric_contraction",
            true,
            vec![term(
                1,
                vec![t("K", vec![d("a"), d("b")]), t("S", vec![u("a"), u("b")])],
            )],
        ),
        (
            "antisymmetric_identical_vectors",
            true,
            vec![term(
                1,
                vec![
                    t("K", vec![d("a"), d("b")]),
                    t("V", vec![u("a")]),
                    t("V", vec![u("b")]),
                ],
            )],
        ),
        (
            "antisymmetric_distinct_vectors",
            false,
            vec![term(
                1,
                vec![
                    t("K", vec![d("a"), d("b")]),
                    t("V", vec![u("a")]),
                    t("W", vec![u("b")]),
                ],
            )],
        ),
        (
            "metric_lowering",
            true,
            vec![
                term(1, vec![F::metric(d("a"), d("b")), t("V", vec![u("b")])]),
                term(-1, vec![t("V", vec![d("a")])]),
            ],
        ),
        (
            "inverse_metric",
            true,
            vec![
                term(
                    1,
                    vec![F::metric(d("a"), d("b")), F::metric(u("b"), u("c"))],
                ),
                term(-1, vec![F::delta(u("c"), d("a"))]),
            ],
        ),
        (
            "metric_trace_four",
            true,
            vec![
                term(
                    1,
                    vec![F::metric(d("a"), d("b")), F::metric(u("a"), u("b"))],
                ),
                term(-4, vec![]),
            ],
        ),
        (
            "metric_trace_wrong_two",
            false,
            vec![
                term(
                    1,
                    vec![F::metric(d("a"), d("b")), F::metric(u("a"), u("b"))],
                ),
                term(-2, vec![]),
            ],
        ),
        (
            "dummy_variance_orientation",
            true,
            vec![
                term(1, vec![t("A", vec![u("a")]), t("B", vec![d("a")])]),
                term(-1, vec![t("A", vec![d("b")]), t("B", vec![u("b")])]),
            ],
        ),
        (
            "dummy_rename_exchange_sign",
            true,
            vec![
                term(
                    1,
                    vec![
                        t("K", vec![d("a"), d("b")]),
                        t("V", vec![u("a")]),
                        t("W", vec![u("b")]),
                    ],
                ),
                term(
                    1,
                    vec![
                        t("W", vec![u("c")]),
                        t("K", vec![d("c"), d("d")]),
                        t("V", vec![u("d")]),
                    ],
                ),
            ],
        ),
    ];
    raw.into_iter()
        .map(|(id, expected_zero, terms)| {
            Ok(Case {
                id,
                expected_zero,
                residual: Expr::from_terms(&c, terms)?,
            })
        })
        .collect()
}
