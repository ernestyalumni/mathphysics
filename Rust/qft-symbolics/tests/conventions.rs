//! F04 convention boundary gate. Copy into qft-symbolics/tests/ unchanged.
//! Independent fixed witnesses target cross-module interpretation, not CAS completeness.
use qft_symbolics::{
    clifford as c, lorentz as l, matrices as m,
    polynomial::Poly,
    scalar::{C, Q},
    smga, spinor as s,
    twistor::{self as t, Twistor},
};
use std::collections::BTreeMap;
fn vec4(v: [i64; 4]) -> [Poly; 4] {
    v.map(Poly::constant)
}
fn eval(v: &[Poly; 4]) -> [C; 4] {
    v.clone().map(|p| p.eval(&BTreeMap::new()).unwrap())
}
fn dot(a: [C; 4], b: [C; 4]) -> C {
    -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}
fn real_twistor(a: i64, b: i64) -> Twistor {
    Twistor {
        phi: [Poly::constant(a), Poly::constant(b)],
        star: [Poly::constant(a), Poly::constant(b)],
    }
}
#[test]
fn slash_metric_bridge_and_formal_dimension_are_distinct() {
    // Timelike non-null witness: p.p=-3, slash(p)^2=+3 I.
    let p = vec4([2, 1, 0, 0]);
    assert_eq!(t::dot(&p, &p), Poly::constant(-3));
    let slash = t::slash(&p).map(|row| row.map(|x| x.eval(&BTreeMap::new()).unwrap()));
    assert_eq!(m::mul(slash, slash), m::scale(m::identity(), C::from(3)));
    assert_ne!(
        m::mul(slash, slash),
        m::scale(m::identity(), C::from(-3)),
        "h=eta must fail"
    );
    let contracted = (0..4).fold(m::zero(), |a, mu| {
        m::add(
            a,
            m::scale(m::mul(m::gamma(mu), m::gamma(mu)), C::from(m::ETA[mu])),
        )
    });
    assert_eq!(contracted, m::scale(m::identity(), C::from(-4)));
    let formal = c::sandwich_srednicki(&[], &Poly::var("d"));
    let scalar = formal.get(&Vec::<String>::new()).unwrap();
    assert_eq!(
        scalar
            .eval(&BTreeMap::from([("d".into(), C::from(5))]))
            .unwrap(),
        C::from(-5)
    );
    assert_ne!(
        contracted,
        m::scale(m::identity(), C::from(-5)),
        "4x4 gamma matrices cannot realize formal d=5 contractions"
    );
}
#[test]
fn spacetime_orientation_and_chirality_survive_a_parity_control() {
    // Ordered contravariant coordinate basis has determinant +1, but covariant
    // contraction epsilon^abcd p_a q_b r_c s_d is -1 in mostly-plus signature.
    let frame = [
        vec4([1, 0, 0, 0]),
        vec4([0, 1, 0, 0]),
        vec4([0, 0, 1, 0]),
        vec4([0, 0, 0, 1]),
    ];
    assert_eq!(l::epsilon_upper_covariant(&frame), Poly::constant(-1));
    assert_eq!(
        l::chiral_four_slash(&frame, l::Chirality::Left),
        Poly::constant(C::from(-2) * C::I)
    );
    let parity = [
        frame[0].clone(),
        frame[1].clone().map(|x| -x),
        frame[2].clone(),
        frame[3].clone(),
    ];
    assert_eq!(
        l::chiral_four_slash(&parity, l::Chirality::Left),
        l::chiral_four_slash(&frame, l::Chirality::Right)
    );
    assert_ne!(
        l::chiral_four_slash(&frame, l::Chirality::Left),
        l::chiral_four_slash(&frame, l::Chirality::Right)
    );
    let g5 = m::scale(
        (0..4).fold(m::identity(), |a, mu| m::mul(a, m::gamma(mu))),
        C::I,
    );
    assert_eq!(m::mul(g5, g5), m::identity());
    assert_eq!(g5[0][0], C::from(-1)); // left projector selects upper chiral block
    assert_eq!(g5[3][3], C::ONE);
}
#[test]
fn spinor_raising_brackets_and_massless_vector_sign_bridge() {
    let (p, q) = (real_twistor(0, 2), real_twistor(2, 0));
    assert_eq!(
        eval(&p.momentum()),
        [2.into(), 0.into(), 0.into(), 2.into()]
    );
    assert_eq!(
        eval(&q.momentum()),
        [2.into(), 0.into(), 0.into(), (-2).into()]
    );
    assert_eq!(t::raise(&p.phi), [Poly::constant(2), Poly::constant(0)]);
    // Lowering after raising is identity; raising twice is -identity.
    let raised = t::raise(&p.phi);
    let lowered = std::array::from_fn(|i| {
        (0..2).fold(Poly::default(), |a, j| {
            a + Poly::constant(m::epsilon_lower()[i][j]) * raised[j].clone()
        })
    });
    assert_eq!(p.phi, lowered);
    assert_eq!(t::square(&q, &p), -s::bracket(&q.phi, &p.phi));
    let invariant = t::angle(&p, &q) * t::square(&q, &p);
    assert_eq!(invariant, Poly::constant(16));
    assert_eq!(
        invariant,
        Poly::constant(-2) * t::dot(&p.momentum(), &q.momentum())
    );
    assert_ne!(
        invariant,
        Poly::constant(2) * t::dot(&p.momentum(), &q.momentum())
    );
}
#[test]
fn scaled_polarizations_retain_normalization_and_reference_domain() {
    let (k, q) = (real_twistor(0, 2), real_twistor(2, 0));
    let plus = t::polarization(&k, &q, true)
        .unwrap()
        .map(|x| x.eval(&BTreeMap::new()).unwrap());
    let minus = t::polarization(&k, &q, false)
        .unwrap()
        .map(|x| x.eval(&BTreeMap::new()).unwrap());
    assert_eq!(plus, [C::ZERO, C::ONE, -C::I, C::ZERO]);
    assert_eq!(minus, [C::ZERO, -C::ONE, -C::I, C::ZERO]);
    assert_eq!(dot(eval(&k.momentum()), plus), C::ZERO);
    assert_eq!(dot(plus, minus), C::from(-2)); // API vectors = sqrt(2)*epsilon
    assert_ne!(
        dot(plus, minus),
        C::from(-1),
        "unscaled physical normalization must not be assigned to API vectors"
    );
    for helicity in [false, true] {
        assert!(t::polarization(&k, &k, helicity).is_err());
        assert!(
            t::polarization(&k, &real_twistor(0, 4), helicity).is_err(),
            "any collinear reference is singular, not just identical objects"
        );
    }
}
#[test]
fn coefficient_conjugation_is_not_the_lorentzian_spinor_involution() {
    let (p, q) = (Twistor::named("p"), Twistor::named("q"));
    // phi_p=(1,i), phi_q=(2,1); stars are conjugates in this positive-energy witness.
    let values = BTreeMap::from([
        ("pp0".into(), C::ONE),
        ("pp1".into(), C::I),
        ("ps0".into(), C::ONE),
        ("ps1".into(), -C::I),
        ("qp0".into(), C::from(2)),
        ("qp1".into(), C::ONE),
        ("qs0".into(), C::from(2)),
        ("qs1".into(), C::ONE),
    ]);
    let angle = t::angle(&p, &q);
    let reversed_square = t::square(&q, &p).eval(&values).unwrap();
    assert_eq!(angle.eval(&values).unwrap().conj(), reversed_square);
    assert_ne!(
        angle.conjugate_coefficients().eval(&values).unwrap(),
        reversed_square,
        "formal star variables are not swapped by coefficient-only conjugation"
    );
    assert_eq!(
        t::dot(&p.momentum(), &p.momentum()).eval(&values).unwrap(),
        C::ZERO
    );
    // The angular constructor is a restricted positive-energy real slice, not
    // an arbitrary complexified/negative-energy twistor constructor.
    assert!(Twistor::from_rational_angles(Q::ZERO, Q::ZERO, Q::ONE, Q::ONE, Q::ZERO).is_err());
}
#[test]
fn bare_amplitudes_and_split_signature_sign_chambers_are_not_physical_normalizations() {
    let (a, b) = ([Q::ONE, Q::ZERO], [Q::ZERO, Q::ONE]);
    let ap = a.map(Poly::constant);
    let bp = b.map(Poly::constant);
    assert_eq!(s::bracket(&ap, &bp), Poly::constant(1));
    assert_eq!(smga::bracket(a, b), Q::from(-1));
    let legs = [
        ap,
        bp,
        [Poly::constant(1), Poly::constant(1)],
        [Poly::constant(1), Poly::constant(-1)],
    ];
    let ratio = s::parke_taylor(&legs, [0, 1])
        .unwrap()
        .eval(&BTreeMap::new())
        .unwrap();
    assert_eq!(ratio, C::from(Q::new(1, 2)));
    assert_ne!(
        ratio,
        C::I * C::from(Q::new(1, 2)),
        "no universal i prefactor supplied"
    );
    assert!(s::parke_taylor(&legs, [0, 0]).is_err());
    let wall = [[Q::from(-2), Q::ZERO], [Q::ONE, Q::ZERO], [Q::ONE, Q::ZERO]];
    assert!(
        smga::formula_r1(&wall).is_err(),
        "conservation does not remove a split-signature sign-chamber wall"
    );
    assert!(
        smga::formula_r1(&[a, b, a]).is_err(),
        "nonconserved split-frame data must fail"
    );
}
