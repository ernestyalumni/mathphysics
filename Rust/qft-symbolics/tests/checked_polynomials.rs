use qft_symbolics::{
    computation::{Budget, Error, Limits, Resource},
    exact::{BigC, BigQ},
    polynomial::{BigPoly, Poly, Ratio},
};
use std::collections::BTreeMap;
fn budget(l: Limits) -> Budget {
    Budget::new(l)
}
fn resource<T: std::fmt::Debug>(r: Result<T, Error>, want: Resource) {
    assert!(
        matches!(r,Err(Error::ResourceLimit{resource,..}) if resource==want),
        "{r:?}"
    )
}
#[test]
fn hardcoded_binomial_and_public_terms() {
    let mut b = Budget::default();
    let x = BigPoly::try_var("x", &mut b).unwrap();
    let one = BigPoly::try_constant(1, &mut b).unwrap();
    let p = x.try_add(&one, &mut b).unwrap().try_pow(3, &mut b).unwrap();
    let expected = BTreeMap::from([
        (vec![], BigC::from(1)),
        (vec!["x".to_owned()], BigC::from(3)),
        (vec!["x".to_owned(); 2], BigC::from(3)),
        (vec!["x".to_owned(); 3], BigC::from(1)),
    ]);
    assert_eq!(p.terms(), &expected);
    assert_eq!(
        p.try_eval(&BTreeMap::from([("x".into(), 2.into())]), &mut b)
            .unwrap(),
        BigC::from(27)
    );
    assert!(p.try_sub(&p, &mut b).unwrap().is_zero());
    assert_eq!(
        x.try_neg(&mut b).unwrap().terms()[&vec!["x".into()]],
        BigC::from(-1)
    );
}
#[test]
fn arbitrary_precision_coefficient_survives_polynomial() {
    let mut b = Budget::new(Limits {
        max_operations: 10_000_000,
        ..Limits::default()
    });
    let q = BigQ::parse("340282366920938463463374607431768211456", "1", &mut b).unwrap();
    let p = BigPoly::try_constant(BigC::new(q, 0.into()), &mut b).unwrap();
    let squared = p.try_mul(&p, &mut b).unwrap();
    assert_eq!(
        squared.to_string(),
        "115792089237316195423570985008687907853269984665640564039457584007913129639936"
    );
}
#[test]
fn budgets_and_zero_shortcuts() {
    let mut b = Budget::default();
    let x = BigPoly::try_var("long", &mut b).unwrap();
    let zero = BigPoly::default();
    resource(
        zero.try_add(
            &zero,
            &mut budget(Limits {
                max_operations: 0,
                ..Limits::default()
            }),
        ),
        Resource::Operations,
    );
    resource(
        zero.try_neg(&mut budget(Limits {
            max_operations: 0,
            ..Limits::default()
        })),
        Resource::Operations,
    );
    resource(
        zero.try_mul(
            &x,
            &mut budget(Limits {
                max_input_bytes: 1,
                ..Limits::default()
            }),
        ),
        Resource::InputBytes,
    );
    resource(
        x.try_pow(
            0,
            &mut budget(Limits {
                max_degree: 0,
                ..Limits::default()
            }),
        ),
        Resource::Degree,
    );
    resource(
        BigPoly::try_constant(
            2,
            &mut budget(Limits {
                max_coefficient_bits: 1,
                ..Limits::default()
            }),
        ),
        Resource::CoefficientBits,
    );
    resource(
        BigPoly::try_var(
            "x",
            &mut budget(Limits {
                max_terms: 0,
                ..Limits::default()
            }),
        ),
        Resource::Terms,
    );
    resource(
        x.try_add(
            &zero,
            &mut budget(Limits {
                max_operations: 0,
                ..Limits::default()
            }),
        ),
        Resource::Operations,
    );
    resource(x.try_pow(65, &mut Budget::default()), Resource::Degree);
    let y = BigPoly::try_var("y", &mut b).unwrap();
    resource(
        x.try_add(
            &y,
            &mut budget(Limits {
                max_terms: 1,
                ..Limits::default()
            }),
        ),
        Resource::Terms,
    );
    resource(
        x.try_mul(
            &y,
            &mut budget(Limits {
                max_degree: 1,
                ..Limits::default()
            }),
        ),
        Resource::Degree,
    );
    assert!(
        matches!(x.try_eval(&BTreeMap::new(),&mut b),Err(Error::MissingVariable{name}) if name=="long")
    );
    let bad = BigC::new(BigQ::parse("256", "1", &mut b).unwrap(), 0.into());
    resource(
        zero.try_eval(
            &BTreeMap::from([("unused".into(), bad)]),
            &mut budget(Limits {
                max_coefficient_bits: 4,
                ..Limits::default()
            }),
        ),
        Resource::CoefficientBits,
    );
}
#[test]
fn ratio_checked_domain_and_legacy_signatures() {
    let x = Poly::var("x");
    let r = Ratio::new(Poly::constant(1), x).unwrap();
    let mut b = Budget::default();
    assert!(r.try_residual(&r, &mut b).unwrap().is_zero());
    assert_eq!(
        r.try_eval(&BTreeMap::from([("x".into(), 0.into())]), &mut b),
        Err(Error::DivisionByZero)
    );
    assert_eq!(
        r.eval(&BTreeMap::from([("x".into(), 2.into())]))
            .unwrap()
            .to_string(),
        "1/2"
    );
    assert!(Ratio::new(Poly::constant(0), Poly::default()).is_err());
}

#[test]
fn legacy_operators_do_not_inherit_default_global_work_limit() {
    let huge = qft_symbolics::scalar::C::new(
        qft_symbolics::scalar::Q::new(1i128 << 100, 1),
        qft_symbolics::scalar::Q::ZERO,
    );
    let p = (0..8).fold(Poly::default(), |p, i| {
        p + Poly::constant(huge) * Poly::var(format!("x{i}"))
    });
    // Explicitly tiny work budget is independent of coefficient backend costs.
    let mut tight = Budget::new(Limits {
        max_operations: 1,
        ..Limits::default()
    });
    assert!(matches!(
        p.try_mul(&Poly::constant(1), &mut tight),
        Err(Error::ResourceLimit {
            resource: Resource::Operations,
            ..
        })
    ));
    assert_eq!(p.clone() * Poly::constant(1), p);
    // Longer names were legal in the original opaque-string API.
    let name = "x".repeat(1_048_577);
    assert_eq!(Poly::var(name).term_count(), 1);
}

#[test]
fn borrowed_names_and_aggregate_symbol_byte_boundaries() {
    let limited = |bytes| {
        Budget::new(Limits {
            max_input_bytes: bytes,
            ..Limits::default()
        })
    };
    assert!(BigPoly::try_var("é", &mut limited(2)).is_ok());
    resource(BigPoly::try_var("é", &mut limited(1)), Resource::InputBytes);
    let a = BigPoly::var("aa");
    let b = BigPoly::var("bb");
    let sum = a.clone() + b.clone();
    assert_eq!(a.try_add(&b, &mut limited(4)).unwrap(), sum);
    resource(a.try_add(&b, &mut limited(3)), Resource::InputBytes);
    // Each monomial fits; the complete operand does not. Validate before zero.
    let zero = BigPoly::default();
    resource(zero.try_mul(&sum, &mut limited(3)), Resource::InputBytes);
    resource(sum.try_mul(&zero, &mut limited(3)), Resource::InputBytes);
    resource(sum.try_sub(&sum, &mut limited(3)), Resource::InputBytes);
    resource(sum.try_pow(0, &mut limited(3)), Resource::InputBytes);
    assert_eq!(sum.try_add(&zero, &mut limited(4)).unwrap(), sum);
    let product = sum.clone() * sum.clone(); // aa^2 + 2 aa bb + bb^2: 12 stored bytes.
    assert_eq!(sum.try_mul(&sum, &mut limited(12)).unwrap(), product);
    resource(sum.try_mul(&sum, &mut limited(11)), Resource::InputBytes);
    // Existing terms must still merge at capacity, without charging their key twice.
    assert_eq!(
        sum.try_add(&sum, &mut limited(4)).unwrap(),
        sum.clone() + sum.clone()
    );
    let names = BTreeMap::from([("aa".into(), BigC::from(1)), ("bb".into(), BigC::from(1))]);
    resource(zero.try_eval(&names, &mut limited(3)), Resource::InputBytes);
    assert_eq!(
        sum.try_eval(&names, &mut limited(4)).unwrap(),
        BigC::from(2)
    );
    struct LegacyName;
    impl From<LegacyName> for String {
        fn from(_: LegacyName) -> String {
            "opaque".into()
        }
    }
    assert_eq!(Poly::var(LegacyName), Poly::var("opaque"));
}

// Observe large key allocations on this test thread only. Small coefficient and
// harness allocations are irrelevant. This tests sequencing, not a heap cap.
struct ObservedAllocator;
thread_local! {
    static OBSERVE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static LARGEST: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}
fn observe_allocation(size: usize) {
    let _ = OBSERVE.try_with(|active| {
        if active.get() {
            let _ = LARGEST.try_with(|max| max.set(max.get().max(size)));
        }
    });
}
unsafe impl std::alloc::GlobalAlloc for ObservedAllocator {
    unsafe fn alloc(&self, layout: std::alloc::Layout) -> *mut u8 {
        observe_allocation(layout.size());
        std::alloc::GlobalAlloc::alloc(&std::alloc::System, layout)
    }
    unsafe fn alloc_zeroed(&self, layout: std::alloc::Layout) -> *mut u8 {
        observe_allocation(layout.size());
        std::alloc::GlobalAlloc::alloc_zeroed(&std::alloc::System, layout)
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: std::alloc::Layout) {
        std::alloc::GlobalAlloc::dealloc(&std::alloc::System, ptr, layout)
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: std::alloc::Layout, size: usize) -> *mut u8 {
        observe_allocation(size);
        std::alloc::GlobalAlloc::realloc(&std::alloc::System, ptr, layout, size)
    }
}
#[global_allocator]
static ALLOCATOR: ObservedAllocator = ObservedAllocator;
fn measured<T>(f: impl FnOnce() -> T) -> (T, usize) {
    LARGEST.with(|v| v.set(0));
    OBSERVE.with(|v| v.set(true));
    let result = f();
    OBSERVE.with(|v| v.set(false));
    (result, LARGEST.with(|v| v.get()))
}
#[test]
fn known_failures_precede_large_key_allocation() {
    let name = "a".repeat(2_000_000);
    let base = Limits {
        max_input_bytes: 3_000_000,
        max_operations: 10_000,
        ..Limits::default()
    };
    // Positive control proves this thread's allocator observer detects owned keys.
    let (p, allocated) = measured(|| BigPoly::try_var(&name, &mut Budget::new(base)));
    assert!(allocated >= name.len());
    let p = p.unwrap();
    for (limits, want) in [
        (
            Limits {
                max_terms: 0,
                ..base
            },
            Resource::Terms,
        ),
        (
            Limits {
                max_degree: 0,
                ..base
            },
            Resource::Degree,
        ),
        (
            Limits {
                max_operations: 3,
                ..base
            },
            Resource::Operations,
        ),
        (
            Limits {
                max_coefficient_bits: 0,
                ..base
            },
            Resource::CoefficientBits,
        ),
    ] {
        let (result, allocated) = measured(|| {
            qft_symbolics::clifford::try_metric_with::<BigC>(&name, "b", &mut Budget::new(limits))
        });
        resource(result, want);
        assert!(
            allocated < name.len(),
            "allocated {allocated} bytes before {want:?} rejection"
        );
    }
    let mut success = Budget::new(base);
    let negated = p.try_neg(&mut success).unwrap();
    assert_eq!(negated.terms().values().next().unwrap(), &BigC::from(-1));
    // Exhaust immediately after validation, then immediately before insertion;
    // derive the latter from logical work, independently of coefficient costs.
    for work in [3, success.used() - 1] {
        let (result, allocated) = measured(|| {
            p.try_neg(&mut Budget::new(Limits {
                max_operations: work,
                ..base
            }))
        });
        resource(result, Resource::Operations);
        assert!(
            allocated < name.len(),
            "cloned key before exhausted work rejection"
        );
    }
    let (result, allocated) = measured(|| {
        p.try_add(
            &BigPoly::default(),
            &mut Budget::new(Limits {
                max_operations: 4,
                ..base
            }),
        )
    });
    resource(result, Resource::Operations);
    assert!(allocated < name.len());
}
