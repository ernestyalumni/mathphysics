use lie_groups::{
    su2,
    symbolic::{self, Chart, SymbolicSu2},
    Matrix, So3Algebra, Su2Algebra, Su2Basis,
};
use qft_symbolics::{polynomial::Poly, scalar::Q};
use std::fmt::Display;
fn quoted(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}
fn expression(value: &impl Display) -> String {
    quoted(&value.to_string().replace("i)", "*I)"))
}
fn matrix<T: Display, const R: usize, const C: usize>(m: &Matrix<T, R, C>) -> String {
    format!(
        "[{}]",
        m.rows()
            .iter()
            .map(|r| format!(
                "[{}]",
                r.iter().map(expression).collect::<Vec<_>>().join(",")
            ))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn symbolic_json() -> Result<String, String> {
    let basis = |b: Option<Su2Basis>| -> String {
        let entries: Vec<_> = (0..3)
            .map(|i| match b {
                None => matrix(&su2::pauli::<Q>(i).unwrap()),
                Some(b) => matrix(&su2::generator::<Q>(i, b).unwrap()),
            })
            .collect();
        format!("[{}]", entries.join(","))
    };
    let e = (0..3)
        .map(|i| matrix(&So3Algebra::<Q>::basis(i).unwrap().matrix()))
        .collect::<Vec<_>>()
        .join(",");
    let hat =
        matrix(&So3Algebra::new(std::array::from_fn(|i| Poly::var(format!("v{i}")))).matrix());
    let axis = [Q::ZERO, Q::ZERO, Q::ONE];
    let r = matrix(&symbolic::so3_exponential(axis).map_err(|e| e.to_string())?);
    let full = symbolic::su2_exponential(axis, Su2Basis::PauliOverI).map_err(|e| e.to_string())?;
    let half =
        symbolic::su2_exponential(axis, Su2Basis::HalfPauliOverI).map_err(|e| e.to_string())?;
    let u = SymbolicSu2::stereographic(
        std::array::from_fn(|i| Poly::var(format!("t{i}"))),
        Chart::Plus,
    )
    .map_err(|e| e.to_string())?;
    let um = u.matrix();
    let rm = u.so3_matrix();
    let eps: Vec<_> = (0..3)
        .flat_map(|i| {
            (0..3)
                .flat_map(move |j| (0..3).map(move |k| lie_groups::levi_civita([i, j, k]).unwrap()))
        })
        .collect();
    Ok(format!("{{\"status\":\"COMPUTED\",\"epsilon3\":{eps:?},\"E\":[{e}],\"sigma\":{},\"T\":{},\"F\":{},\"hat\":{hat},\"r_exp_z\":{r},\"u_exp_full_z\":{},\"u_exp_half_z\":{},\"full_argument_scale\":{},\"half_argument_scale\":{},\"chart_u\":{{\"numerator\":{},\"denominator\":{}}},\"chart_r\":{{\"numerator\":{},\"denominator\":{}}}}}",
        basis(None),basis(Some(Su2Basis::HalfPauliOverI)),basis(Some(Su2Basis::PauliOverI)),matrix(&full.matrix),matrix(&half.matrix),quoted(&full.argument_scale.to_string()),quoted(&half.argument_scale.to_string()),matrix(um.numerator()),expression(um.denominator()),matrix(rm.numerator()),expression(rm.denominator())))
}
fn fixtures() -> Result<String, String> {
    let mut rows = Vec::new();
    for i in 0..25 {
        let x = i as f64;
        let v = if i == 0 {
            [0.; 3]
        } else {
            [
                0.13 * (x - 12.),
                -0.07 * (x + 2.),
                0.03 * ((i * i % 17 + 1) as f64),
            ]
        };
        let r = So3Algebra::new(v).exp().map_err(|e| e.to_string())?;
        let u = Su2Algebra::new(v).exp().map_err(|e| e.to_string())?;
        rows.push(format!(
            "{{\"id\":{i},\"v\":{v:?},\"rotation\":{:?},\"quaternion\":{:?}}}",
            r.matrix().rows(),
            u.quaternion()
        ));
    }
    Ok(format!(
        "{{\"status\":\"COMPUTED\",\"scalar\":\"f64\",\"cases\":[{}]}}",
        rows.join(",")
    ))
}
fn run(args: &[String]) -> Result<i32, String> {
    match args.first().map(String::as_str).unwrap_or("verify") {
        "symbolic" if args.len() == 1 => {
            println!("{}", symbolic_json()?);
            Ok(0)
        }
        "fixtures" if args.len() == 1 => {
            println!("{}", fixtures()?);
            Ok(0)
        }
        "verify" if args.len() <= 2 && args.get(1).is_none_or(|s| s == "--wrong-sign") => {
            let wrong = args.len() == 2;
            let mut passed = 0;
            for i in 0..3 {
                for j in 0..3 {
                    let mut target = Matrix::<Q, 3, 3>::zero();
                    for k in 0..3 {
                        let coefficient = i64::from(lie_groups::levi_civita([i, j, k]).unwrap())
                            * if wrong { -1 } else { 1 };
                        target = &target
                            + &So3Algebra::<Q>::basis(k)
                                .unwrap()
                                .matrix()
                                .scaled(Q::from(coefficient));
                    }
                    if So3Algebra::<Q>::basis(i)
                        .unwrap()
                        .matrix()
                        .commutator(&So3Algebra::<Q>::basis(j).unwrap().matrix())
                        == target
                    {
                        passed += 1;
                    }
                }
            }
            println!("{{\"status\":{},\"basis_commutators_passed\":{passed},\"required\":9,\"scope\":\"exact Frankel so3 commutators, not arbitrary candidate ingestion\"}}",quoted(if passed==9{"PASS"}else{"FAIL"}));
            Ok(if passed == 9 { 0 } else { 1 })
        }
        _ => Err("commands: verify [--wrong-sign] | symbolic | fixtures".into()),
    }
}
fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let code = match std::panic::catch_unwind(|| run(&args)) {
        Ok(Ok(code)) => code,
        Ok(Err(e)) => {
            println!("{{\"status\":\"INCONCLUSIVE\",\"error\":{}}}", quoted(&e));
            2
        }
        Err(_) => {
            println!("{{\"status\":\"INCONCLUSIVE\",\"error\":\"arithmetic/resource/internal failure\"}}");
            2
        }
    };
    std::process::exit(code);
}
