use qft_symbolics::{
    clifford::{self as c, NormalWord},
    matrices as m,
    polynomial::Poly,
    scalar::C,
    tensor::{Context, Dimension, Expr, Factor, Index, Space, Term},
};
use std::collections::BTreeMap;
fn add(a: &mut NormalWord, b: NormalWord, coefficient: Poly) {
    for (w, p) in b {
        let value = a.remove(&w).unwrap_or_default() + p * coefficient.clone();
        if !value.is_zero() {
            a.insert(w, value);
        }
    }
}
fn tensor_dimension() -> Poly {
    let context = Context::default()
        .with_lorentz_dimension(Dimension::Symbol("d".into()))
        .unwrap();
    let dn = |n| Index::down(Space::Lorentz, n);
    let up = |n| Index::up(Space::Lorentz, n);
    Expr::from_terms(
        &context,
        vec![Term::new(
            Poly::constant(1),
            vec![
                Factor::metric(dn("a"), dn("b")),
                Factor::metric(up("a"), up("b")),
            ],
        )],
    )
    .unwrap()
    .scalar()
    .unwrap()
}
#[test]
fn ch47_2_formal_dimension_equations() {
    let d = tensor_dimension();
    assert_eq!(d, Poly::var("d"));
    let lhs = c::sandwich_srednicki(&["a", "b"], &d);
    let mut rhs = NormalWord::from([(vec![], Poly::constant(-4) * c::metric("a", "b"))]); // eta=-h
    add(
        &mut rhs,
        c::normal_order(&["a", "b"]),
        Poly::constant(4) - d.clone(),
    );
    assert_eq!(lhs, rhs);
    let lhs = c::sandwich_srednicki(&["a", "b", "c"], &d);
    let mut rhs = NormalWord::new();
    add(
        &mut rhs,
        c::normal_order(&["c", "b", "a"]),
        Poly::constant(2),
    );
    add(
        &mut rhs,
        c::normal_order(&["a", "b", "c"]),
        d - Poly::constant(4),
    );
    assert_eq!(lhs, rhs);
}
#[test]
fn contracted_zero_and_one_gamma_words() {
    let d = Poly::var("d");
    assert_eq!(
        c::sandwich_srednicki(&[], &d),
        NormalWord::from([(vec![], -d.clone())])
    );
    let mut rhs = NormalWord::new();
    add(&mut rhs, c::normal_order(&["a"]), d - Poly::constant(2));
    assert_eq!(c::sandwich_srednicki(&["a"], &Poly::var("d")), rhs);
}
#[test]
fn explicit_contractions_in_d1_through_d4_independently_match_formal_rule() {
    let names = ["a", "b", "c", "d"];
    let mut cases = 0;
    for dimension in 1usize..=4 {
        for k in 0..=4 {
            for sample in 0..dimension.pow(k as u32) {
                let mut remaining = sample;
                let components: Vec<_> = (0..k)
                    .map(|_| {
                        let value = remaining % dimension;
                        remaining /= dimension;
                        value
                    })
                    .collect();
                let labels = &names[..k];
                let mut values = BTreeMap::new();
                for i in 0..k {
                    for j in i..k {
                        let key = format!("h(1:{};1:{})", names[i], names[j]);
                        let h = if components[i] == components[j] {
                            -m::ETA[components[i]]
                        } else {
                            0
                        };
                        values.insert(key, C::from(h));
                    }
                }
                let word = components
                    .iter()
                    .fold(m::identity(), |p, mu| m::mul(p, m::gamma(*mu)));
                let mut expected = m::zero();
                for mu in 0..dimension {
                    expected = m::add(
                        expected,
                        m::scale(
                            m::mul(m::gamma(mu), m::mul(word, m::gamma(mu))),
                            C::from(m::ETA[mu]),
                        ),
                    );
                }
                let mut actual = m::zero();
                for (normal, p) in c::sandwich_srednicki(labels, &Poly::constant(dimension as i64))
                {
                    let mat = normal.iter().fold(m::identity(), |p, label| {
                        let pos = names.iter().position(|s| s == label).unwrap();
                        m::mul(p, m::gamma(components[pos]))
                    });
                    actual = m::add(actual, m::scale(mat, p.eval(&values).unwrap()));
                }
                assert_eq!(actual, expected, "d={dimension} k={k} sample={sample}");
                cases += 1;
            }
        }
    }
    assert_eq!(cases, 498); // 494 nonempty assignments, plus one empty word per dimension.
}
