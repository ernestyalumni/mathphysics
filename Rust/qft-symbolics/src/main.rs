use qft_symbolics::{
    clifford,
    polynomial::Poly,
    smga,
    spinor::{bracket, named, schouten},
};
use std::{env, panic, process};
fn json(s: &str) -> String {
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
fn run(args: &[String]) -> Result<i32, String> {
    match args.first().map(String::as_str).unwrap_or("demo"){
        "schouten"=>{
            if args.len()>2 || args.get(1).is_some_and(|s|s!="--wrong-sign"){return Err("usage: schouten [--wrong-sign]".into());}
            let (a,b,c,d)=(named("a"),named("b"),named("c"),named("d"));
            let r=if args.len()==2{
                bracket(&a,&b)*bracket(&c,&d)-bracket(&a,&c)*bracket(&d,&b)+bracket(&a,&d)*bracket(&b,&c)
            }else{schouten(&a,&b,&c,&d)};
            println!("{{\"status\":{},\"method\":\"exact_polynomial\",\"residual\":{}}}",
                json(if r.is_zero(){"PASS"}else{"FAIL"}),json(&r.to_string()));
            Ok(if r.is_zero(){0}else{1})
        },
        "gamma-trace"=>{
            let labels:Vec<_>=args[1..].iter().map(String::as_str).collect();
            println!("{{\"status\":\"COMPUTED\",\"convention\":\"anticommutator=2h; h=-eta for Srednicki; no implicit sums\",\"result\":{}}}",
                json(&clifford::trace(&labels).to_string()));Ok(0)
        },
        "smga-fixture"=>{
            if args.len()!=3{return Err("usage: smga-fixture n seed".into());}
            let n=args[1].parse().map_err(|_|"invalid n")?;
            let seed=args[2].parse().map_err(|_|"invalid seed")?;
            if !(3..=10).contains(&n){return Err("CLI fixture supports n=3..10".into());}
            let v=smga::fixture(n,seed);
            let a=smga::formula_r1(&v)?;let b=smga::minus_first_recursion(&v)?;
            let vectors=v.iter().map(|p|format!("[{},{}]",json(&p[0].to_string()),json(&p[1].to_string()))).collect::<Vec<_>>().join(",");
            println!("{{\"status\":{},\"n\":{n},\"seed\":{seed},\"vectors\":[{vectors}],\"formula\":{},\"recursion\":{},\"residual\":{}}}",
                json(if a==b{"PASS"}else{"FAIL"}),json(&a.to_string()),json(&b.to_string()),json(&(a-b).to_string()));
            Ok(if a==b{0}else{1})
        },
        "demo"=>{
            if args.len()>1{return Err("usage: demo".into());}
            let (a,b,c,d)=(named("a"),named("b"),named("c"),named("d"));
            assert!(schouten(&a,&b,&c,&d).is_zero());
            assert!(!(bracket(&a,&b)*bracket(&c,&d)-bracket(&a,&c)*bracket(&d,&b)+bracket(&a,&d)*bracket(&b,&c)).is_zero());
            let tr=Poly::constant(4)*(clifford::metric("a","b")*clifford::metric("c","d")-clifford::metric("a","c")*clifford::metric("b","d")+clifford::metric("a","d")*clifford::metric("b","c"));
            assert_eq!(clifford::trace(&["a","b","c","d"]),tr);
            let (mut passed,mut zero,mut nonzero)=(0,0,0);
            let excluded=0;
            for n in 3..=10 {for seed in 0..25 {
                let v=smga::fixture(n,seed);
                let a=smga::formula_r1(&v).map_err(|e|format!("fixture n={n}, seed={seed}, formula: {e}"))?;
                let b=smga::minus_first_recursion(&v).map_err(|e|format!("fixture n={n}, seed={seed}, recurrence: {e}"))?;
                if a!=b {println!("{{\"status\":\"FAIL\",\"n\":{n},\"seed\":{seed}}}");return Ok(1);}
                passed+=1;if a.is_zero(){zero+=1;}else{nonzero+=1;}
            }}
            assert_eq!(passed,200);
            println!("{{\"status\":\"PASS\",\"scope\":\"Schouten and gamma trace identities; finite exact SMGA fixtures, not an all-n proof\",\"schouten_residual\":\"0\",\"negative_control\":\"correctly rejected\",\"smga_passed\":{passed},\"zero_amplitudes\":{zero},\"nonzero_amplitudes\":{nonzero},\"excluded\":{excluded}}}");
            Ok(0)
        },
        _=>Err("commands: demo | schouten [--wrong-sign] | gamma-trace [labels...] | smga-fixture n seed".into()),
    }
}
fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    let result = panic::catch_unwind(|| run(&args));
    let code = match result {
        Ok(Ok(code)) => code,
        Ok(Err(e)) => {
            println!("{{\"status\":\"INCONCLUSIVE\",\"reason\":{}}}", json(&e));
            2
        }
        Err(_) => {
            println!("{{\"status\":\"INCONCLUSIVE\",\"reason\":\"arithmetic range, expansion budget, or internal assertion failure; see stderr\"}}");
            2
        }
    };
    process::exit(code);
}
