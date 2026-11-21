//! ym-audit: does a Lean statement say what an informal statement says?
//!
//! Subcommands:
//!   export  --formal DIR --file F.lean --out F.exports.json   run Lean (Docker, no network) and collect exports
//!   check   --exports F.exports.json --lean NAME --informal NAME [--formalizer] [--trials N]
//!   scan    [--trials N]                                        self-check and hypothesis scan of the informal statements
//!   prepare --formal DIR --run FORMALIZER_RUN --out DIR         benchmark Lean files and item lists
//!   bench   --dir DIR --out REPORT.json [--threads N] [--trials N]
//!   score   --report REPORT.json
mod audit;
mod lean;
mod num;
mod paper;
mod prepare;
mod rng;
mod value;

use serde_json::{json, Value as J};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const IMAGE: &str = "claw-work/lean-agent:lean4.35.0-mathlib";
const SEED: u64 = 20261121;
const TRIALS: usize = 900;

fn flag(args: &[String], name: &str) -> Option<String> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1).cloned())
}
fn need(args: &[String], name: &str) -> String {
    flag(args, name).unwrap_or_else(|| {
        eprintln!("missing {name}");
        std::process::exit(2)
    })
}

/// Run Lean on `file` in Docker and collect `YM_EXPORT` lines and error line numbers.
fn export(formal: &Path, file: &Path) -> Result<J, String> {
    let id = |f: &str| String::from_utf8(Command::new("id").arg(f).output().unwrap().stdout).unwrap().trim().to_string();
    let dir = file.parent().ok_or("file directory")?.canonicalize().map_err(|e| e.to_string())?;
    let name = file.file_name().ok_or("file name")?.to_string_lossy().to_string();
    let output = Command::new("docker")
        .args(["run", "--rm", "--network", "none", "-u", &format!("{}:{}", id("-u"), id("-g"))])
        .args(["-v", &format!("{}:/work:ro", formal.canonicalize().map_err(|e| e.to_string())?.display())])
        .args(["-v", &format!("{}:/candidate:ro", dir.display()), "-w", "/work", "--entrypoint", "bash", IMAGE, "-c"])
        .arg(format!("LEAN_PATH=/work/build:$LEAN_PATH timeout 3600 lean /candidate/{name}"))
        .output()
        .map_err(|e| e.to_string())?;
    let text = String::from_utf8_lossy(&output.stdout).to_string() + &String::from_utf8_lossy(&output.stderr);
    let mut exports = serde_json::Map::new();
    let mut errors = Vec::new();
    let prefix = format!("/candidate/{name}:");
    for line in text.lines() {
        if let Some(p) = line.find("YM_EXPORT ") {
            let j: J = serde_json::from_str(&line[p + 10..]).map_err(|e| format!("export JSON: {e}"))?;
            exports.insert(j["name"].as_str().unwrap_or("").to_string(), j);
        } else if let Some(rest) = line.strip_prefix(&prefix) {
            let mut it = rest.splitn(3, ':');
            let (l, _, msg) = (it.next(), it.next(), it.next().unwrap_or(""));
            if msg.trim_start().starts_with("error") {
                errors.push(json!({"line": l.and_then(|x| x.parse::<u64>().ok()), "message": msg.trim()}));
            }
        }
    }
    Ok(json!({"file": name, "exports": exports, "errors": errors, "returncode": output.status.code()}))
}

fn load(path: &Path) -> J {
    serde_json::from_str(&fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))).unwrap()
}

fn check_one(export: &J, informal: &str, formalizer: bool, trials: usize) -> J {
    let Some(p) = paper::by_name(informal) else { return json!({"verdict": "ERROR", "flags": ["no informal statement"]}) };
    match lean::Statement::from_export(export) {
        Ok(st) => audit::audit(&st, &p, formalizer, trials, SEED),
        Err(e) => json!({"verdict": "ERROR", "flags": ["ERROR"], "first_error": e}),
    }
}

fn bench(dir: &Path, out: &Path, threads: usize, trials: usize) {
    let mut items: Vec<J> = Vec::new();
    for f in ["formalizer-items.json", "reference-items.json", "mutant-items.json"] {
        if let J::Array(a) = load(&dir.join(f)) {
            items.extend(a);
        }
    }
    let mut exports: HashMap<String, J> = HashMap::new();
    for entry in fs::read_dir(dir).unwrap() {
        let p = entry.unwrap().path();
        if p.to_string_lossy().ends_with(".exports.json") {
            let j = load(&p);
            exports.insert(j["file"].as_str().unwrap().to_string(), j);
        }
    }
    let jobs: Vec<(usize, J, Option<J>, Option<String>)> = items
        .iter()
        .enumerate()
        .map(|(k, it)| {
            if let Some(why) = it.get("excluded") {
                return (k, it.clone(), None, why.as_str().map(String::from));
            }
            let file = it["file"].as_str().unwrap_or("");
            let Some(ex) = exports.get(file) else { return (k, it.clone(), None, Some(format!("no export of {file}"))) };
            let (a, b) = (it["lines"][0].as_u64().unwrap_or(0), it["lines"][1].as_u64().unwrap_or(0));
            let err = ex["errors"].as_array().unwrap().iter().find(|e| {
                let l = e["line"].as_u64().unwrap_or(0);
                a > 0 && l >= a && l <= b
            });
            if let Some(e) = err {
                return (k, it.clone(), None, Some(format!("does not elaborate: {}", e["message"].as_str().unwrap_or(""))));
            }
            match ex["exports"].get(it["lean"].as_str().unwrap_or("")) {
                Some(j) => (k, it.clone(), Some(j.clone()), None),
                None => (k, it.clone(), None, Some("does not elaborate (no export)".into())),
            }
        })
        .collect();
    let jobs = std::sync::Arc::new(std::sync::Mutex::new(jobs));
    let results = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
    let mut handles = Vec::new();
    for _ in 0..threads {
        let (jobs, results) = (jobs.clone(), results.clone());
        handles.push(std::thread::spawn(move || loop {
            let job = jobs.lock().unwrap().pop();
            let Some((k, it, ex, why)) = job else { break };
            let started = std::time::Instant::now();
            let report = match (ex, why) {
                (_, Some(why)) => json!({"verdict": "EXCLUDED", "why": why}),
                (Some(ex), None) => {
                    let formalizer = it["kind"] == "formalizer";
                    check_one(&ex, it["informal"].as_str().unwrap_or(""), formalizer, trials)
                }
                _ => unreachable!(),
            };
            eprintln!("{} {} {:?} ({:.0}s)", it["key"], report["verdict"], report["flags"], started.elapsed().as_secs_f64());
            results.lock().unwrap().push((k, json!({"item": it, "report": report})));
        }));
    }
    for h in handles {
        h.join().unwrap();
    }
    let mut res = std::mem::take(&mut *results.lock().unwrap());
    res.sort_by_key(|(k, _)| *k);
    let all: Vec<J> = res.into_iter().map(|(_, j)| j).collect();
    fs::write(out, serde_json::to_string_pretty(&all).unwrap()).unwrap();
}

fn signature(r: &J) -> String {
    let mut f: Vec<String> = r["flags"].as_array().map(|a| a.iter().map(|x| x.as_str().unwrap_or("").to_string()).collect()).unwrap_or_default();
    f.sort();
    format!("{}|weaker={}", f.join(","), r["weaker"])
}

fn score(report: &Path) -> J {
    let rows = load(report);
    let rows = rows.as_array().unwrap();
    // (a) formalizer
    let mut a = Vec::new();
    let (mut misses, mut false_alarms, mut correct, mut excluded) = (0, 0, 0, 0);
    for r in rows.iter().filter(|r| r["item"]["kind"] == "formalizer") {
        let grade = r["item"]["grade"].as_str().unwrap_or("");
        let verdict = r["report"]["verdict"].as_str().unwrap_or("");
        let outcome = match (grade, verdict) {
            (_, "EXCLUDED") => {
                excluded += 1;
                "excluded"
            }
            ("wrong", "PASS") => {
                misses += 1;
                "MISS"
            }
            ("wrong", _) => {
                correct += 1;
                "flagged (correct)"
            }
            (_, "PASS") => {
                correct += 1;
                "pass (correct)"
            }
            _ => {
                false_alarms += 1;
                "FALSE ALARM"
            }
        };
        a.push(json!({"key": r["item"]["key"], "grade": grade, "verdict": verdict, "flags": r["report"]["flags"], "outcome": outcome,
                      "why": r["report"]["why"]}));
    }
    // (b) references
    let refs: Vec<&J> = rows.iter().filter(|r| r["item"]["kind"] == "reference").collect();
    let ref_pass = refs.iter().filter(|r| r["report"]["verdict"] == "PASS").count();
    let mut coverage: HashMap<String, std::collections::BTreeSet<String>> = HashMap::new();
    for r in &refs {
        let e = coverage.entry(r["item"]["informal"].as_str().unwrap().to_string()).or_default();
        for c in r["report"]["covered_informal_conjuncts"].as_array().into_iter().flatten() {
            e.insert(c.as_str().unwrap().to_string());
        }
    }
    let parents: HashMap<String, String> = refs
        .iter()
        .map(|r| (r["item"]["key"].as_str().unwrap().to_string(), signature(&r["report"])))
        .collect();
    // (c) mutants
    let muts: Vec<&J> = rows.iter().filter(|r| r["item"]["kind"] == "mutant" && r["report"]["verdict"] != "EXCLUDED").collect();
    let mut undetected = Vec::new();
    let mut detected = 0;
    for r in &muts {
        let parent = parents.get(r["item"]["parent"].as_str().unwrap()).cloned().unwrap_or_default();
        if signature(&r["report"]) != parent {
            detected += 1;
        } else {
            undetected.push(json!({"key": r["item"]["key"], "mutation": r["item"]["mutation"], "statement": r["item"]["statement"]}));
        }
    }
    let mutants_total = rows.iter().filter(|r| r["item"]["kind"] == "mutant").count();
    json!({
        "a_formalizer": {"scored": correct + misses + false_alarms, "correct": correct, "misses": misses,
                          "false_alarms": false_alarms, "excluded": excluded, "rows": a},
        "b_references": {"pass": ref_pass, "of": refs.len(),
                          "rows": refs.iter().map(|r| json!({"key": r["item"]["key"], "verdict": r["report"]["verdict"], "flags": r["report"]["flags"], "weaker": r["report"]["weaker"]})).collect::<Vec<_>>(),
                          "coverage": coverage.iter().map(|(k, v)| (k.clone(), json!(v))).collect::<serde_json::Map<_, _>>()},
        "c_mutants": {"generated": mutants_total, "elaborated": muts.len(), "detected": detected,
                       "rate": if muts.is_empty() { 0.0 } else { detected as f64 / muts.len() as f64 }, "undetected": undetected},
    })
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(String::as_str).unwrap_or("");
    let trials = flag(&args, "--trials").and_then(|t| t.parse().ok()).unwrap_or(TRIALS);
    match cmd {
        "export" => {
            let j = export(Path::new(&need(&args, "--formal")), Path::new(&need(&args, "--file"))).unwrap();
            eprintln!("{} exports, {} errors", j["exports"].as_object().unwrap().len(), j["errors"].as_array().unwrap().len());
            fs::write(need(&args, "--out"), serde_json::to_string(&j).unwrap()).unwrap();
        }
        "check" => {
            let ex = load(Path::new(&need(&args, "--exports")));
            let name = need(&args, "--lean");
            let j = ex["exports"].get(&name).unwrap_or_else(|| panic!("no export {name}"));
            let r = check_one(j, &need(&args, "--informal"), args.iter().any(|a| a == "--formalizer"), trials);
            println!("{}", serde_json::to_string_pretty(&r).unwrap());
        }
        "scan" => {
            let only = flag(&args, "--only");
            let out: Vec<J> = paper::statements()
                .iter()
                .filter(|p| only.as_deref().map_or(true, |o| o.split(',').any(|x| x == p.name)))
                .map(|p| audit::scan(p, trials * 4, SEED))
                .collect();
            println!("{}", serde_json::to_string_pretty(&out).unwrap());
        }
        "prepare" => {
            let (formal, out) = (PathBuf::from(need(&args, "--formal")), PathBuf::from(need(&args, "--out")));
            prepare::formalizer(Path::new(&need(&args, "--run")), &out).unwrap();
            prepare::references(&out).unwrap();
            // hypothesis names come from the reference exports
            let refs = out.join("References.exports.json");
            if !refs.exists() {
                let j = export(&formal, &out.join("References.lean")).unwrap();
                fs::write(&refs, serde_json::to_string(&j).unwrap()).unwrap();
            }
            let ex = load(&refs);
            let hyps = |thm: &str| -> Vec<String> {
                ex["exports"]
                    .get(format!("YMFormal.{thm}"))
                    .and_then(|j| lean::Statement::from_export(j).ok())
                    .map(|st| st.hyps().map(|b| b.name.clone()).collect())
                    .unwrap_or_default()
            };
            prepare::mutant_files(&formal, &out, &hyps).unwrap();
        }
        "bench" => {
            let threads = flag(&args, "--threads").and_then(|t| t.parse().ok()).unwrap_or(4);
            bench(Path::new(&need(&args, "--dir")), Path::new(&need(&args, "--out")), threads, trials);
        }
        "score" => println!("{}", serde_json::to_string_pretty(&score(Path::new(&need(&args, "--report")))).unwrap()),
        _ => {
            eprintln!("usage: ym-audit export|check|scan|prepare|bench|score ...");
            std::process::exit(2);
        }
    }
}
