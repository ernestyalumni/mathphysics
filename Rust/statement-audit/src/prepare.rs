//! Lean files for the benchmark: formalizer candidates, reference theorems and mutants.
use serde_json::{json, Value as J};
use std::fs;
use std::path::Path;

pub const REFERENCES: &[(&str, &str, &str)] = &[
    ("Purity", "purity_bounds", "purity_bounds"),
    ("Purity", "exp_neg_gap_le", "gap_from_purity"),
    ("Purity", "gap_from_purity", "gap_from_purity"),
    ("Defect", "defect_telescope", "defect_telescope"),
    ("Defect", "defect_nonneg_and_last", "defect_telescope"),
    ("Defect", "gap_from_defect", "gap_from_defect"),
    ("ConditionalGap", "conditional_gap", "conditional_gap"),
    ("Casimir", "trace_sandwich", "trace_ad_ad"),
    ("Casimir", "trace_adjoint_comp", "trace_ad_ad"),
];

fn header() -> String {
    "import Mathlib\nimport Aesop\nimport Audit.Export\n\n".to_string()
}

/// The 28 formalizer answers, with `λ` renamed to `lam`, each in its own namespace.
pub fn formalizer(run: &Path, out: &Path) -> Result<(), String> {
    let compile = fs::read_to_string(run.join("compile.jsonl")).map_err(|e| e.to_string())?;
    let grades: J = serde_json::from_str(&fs::read_to_string(run.join("grades.json")).map_err(|e| e.to_string())?)
        .map_err(|e| e.to_string())?;
    let mut text = header();
    let mut items = Vec::new();
    for line in compile.lines() {
        let row: J = serde_json::from_str(line).map_err(|e| e.to_string())?;
        let name = row["name"].as_str().unwrap_or("");
        let attempt = row["attempt"].as_u64().unwrap_or(0);
        let key = format!("{name}#{attempt}");
        let grade = grades[&key]["grade"].as_str().unwrap_or("ungraded");
        let Some(code) = row["code"].as_str() else {
            items.push(json!({"key": key, "informal": name, "kind": "formalizer", "grade": grade, "excluded": "no code"}));
            continue;
        };
        let code = code.replace('λ', "lam");
        let body: Vec<&str> = code.lines().filter(|l| !l.trim_start().starts_with("import ")).collect();
        let theorem = body
            .iter()
            .find_map(|l| l.trim_start().strip_prefix("theorem ").map(|r| r.split_whitespace().next().unwrap_or("").to_string()))
            .unwrap_or_default();
        let ns = format!("F_{name}_{attempt}");
        let start = text.lines().count() + 1;
        text.push_str(&format!("section\nnamespace {ns}\n{}\nend {ns}\nend\n", body.join("\n")));
        let end = text.lines().count();
        let lean = format!("{ns}.{theorem}");
        text.push_str(&format!("#ym_export {lean}\n\n"));
        items.push(json!({"key": key, "informal": name, "kind": "formalizer", "grade": grade,
                          "lean": lean, "file": "Formalizer.lean", "lines": [start, end]}));
    }
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    fs::write(out.join("Formalizer.lean"), text).map_err(|e| e.to_string())?;
    fs::write(out.join("formalizer-items.json"), serde_json::to_string_pretty(&items).unwrap()).map_err(|e| e.to_string())?;
    Ok(())
}

pub fn references(out: &Path) -> Result<(), String> {
    let mut text = "import Mathlib\nimport Audit.Export\nimport YMFormal.Purity\nimport YMFormal.Defect\nimport YMFormal.ConditionalGap\nimport YMFormal.Casimir\n\n".to_string();
    let mut items = Vec::new();
    for (_, thm, informal) in REFERENCES {
        text.push_str(&format!("#ym_export YMFormal.{thm}\n"));
        items.push(json!({"key": thm, "informal": informal, "kind": "reference", "lean": format!("YMFormal.{thm}"),
                          "file": "References.lean", "lines": [0, 0]}));
    }
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    fs::write(out.join("References.lean"), text).map_err(|e| e.to_string())?;
    fs::write(out.join("reference-items.json"), serde_json::to_string_pretty(&items).unwrap()).map_err(|e| e.to_string())?;
    Ok(())
}

/// The statement text of `theorem name` in `src`: from just after the name to just before `:=`.
fn statement_span(src: &str, name: &str) -> Option<(usize, usize)> {
    let marker = format!("theorem {name} ");
    let at = src.find(&marker)? + marker.len();
    let rest = &src[at..];
    let end = rest.find(" := by").or_else(|| rest.find(" :=\n"))?;
    Some((at, at + end))
}

fn is_ident(c: char) -> bool {
    c.is_alphanumeric() || c == '_' || c == '\'' || ('₀'..='₉').contains(&c) || c == '.'
}

/// Single-edit mutants of a statement text.
pub fn mutants(stmt: &str, hyps: &[String]) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    let mut push = |kind: &str, s: String| {
        if s != stmt && !out.iter().any(|(_, t)| *t == s) {
            out.push((kind.to_string(), s));
        }
    };
    let chars: Vec<(usize, char)> = stmt.char_indices().collect();
    for (k, &(i, c)) in chars.iter().enumerate() {
        let next = chars.get(k + 1).map(|x| x.1).unwrap_or(' ');
        let prev = if k > 0 { chars[k - 1].1 } else { ' ' };
        let rest = &stmt[i + c.len_utf8()..];
        match c {
            '≤' => push("le→lt", format!("{}<{}", &stmt[..i], rest)),
            '<' if !matches!(next, '|' | ';' | '-' | '=' | '>') => push("lt→le", format!("{}≤{}", &stmt[..i], rest)),
            '-' if prev == ' ' && next == ' ' => push("minus→plus", format!("{}+{}", &stmt[..i], rest)),
            _ => {}
        }
        if c.is_ascii_digit() && !is_ident(prev) {
            let len = stmt[i..].chars().take_while(|d| d.is_ascii_digit()).count();
            let after = stmt[i + len..].chars().next().unwrap_or(' ');
            if !is_ident(after) {
                let n: u64 = stmt[i..i + len].parse().unwrap_or(0);
                push("numeral", format!("{}{}{}", &stmt[..i], n + 1, &stmt[i + len..]));
            }
        }
    }
    for (from, to) in [("j + 1", "j"), ("k + 1", "k"), ("lam 0", "lam 1"), ("lam 1", "lam 2"), ("exp (-", "exp (")] {
        let mut start = 0;
        while let Some(p) = stmt[start..].find(from) {
            let at = start + p;
            push("index/sign", format!("{}{}{}", &stmt[..at], to, &stmt[at + from.len()..]));
            start = at + from.len();
        }
    }
    for h in hyps {
        let open = format!("({h} :");
        if let Some(a) = stmt.find(&open) {
            let mut depth = 0;
            for (j, c) in stmt[a..].char_indices() {
                match c {
                    '(' => depth += 1,
                    ')' => {
                        depth -= 1;
                        if depth == 0 {
                            let b = a + j + 1;
                            let s = format!("{}{}", stmt[..a].trim_end(), &stmt[b..]);
                            push("drop hypothesis", s);
                            break;
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    out
}

/// Mutant files: a copy of each module with the mutants appended inside `namespace YMFormal`.
pub fn mutant_files(formal: &Path, out: &Path, hyps: &dyn Fn(&str) -> Vec<String>) -> Result<(), String> {
    fs::create_dir_all(out).map_err(|e| e.to_string())?;
    let mut items = Vec::new();
    for module in ["Purity", "Defect", "ConditionalGap", "Casimir"] {
        let src = fs::read_to_string(formal.join(format!("YMFormal/{module}.lean"))).map_err(|e| e.to_string())?;
        let end_at = src.rfind("end YMFormal").ok_or("module without `end YMFormal`")?;
        let mut text = src[..end_at].to_string();
        text = format!("import Audit.Export\n{text}");
        let mut exports = String::new();
        for (m, thm, informal) in REFERENCES {
            if *m != module {
                continue;
            }
            let (a, b) = statement_span(&src, thm).ok_or(format!("statement of {thm}"))?;
            let stmt = &src[a..b];
            for (k, (kind, s)) in mutants(stmt, &hyps(thm)).into_iter().enumerate() {
                let name = format!("mut_{thm}_{k}");
                text.push_str(&format!("\ntheorem {name} {s} := by sorry\n"));
                let at = text.rfind(&format!("theorem {name} ")).unwrap();
                let start = text[..at].matches('\n').count() + 1;
                let end = text.trim_end_matches('\n').matches('\n').count() + 1;
                exports.push_str(&format!("#ym_export YMFormal.{name}\n"));
                items.push(json!({"key": name, "informal": informal, "kind": "mutant", "parent": thm, "mutation": kind,
                                  "lean": format!("YMFormal.{name}"), "file": format!("Mutants{module}.lean"),
                                  "lines": [start, end], "statement": s}));
            }
        }
        text.push_str("\nend YMFormal\n\n");
        text.push_str(&exports);
        fs::write(out.join(format!("Mutants{module}.lean")), text).map_err(|e| e.to_string())?;
    }
    fs::write(out.join("mutant-items.json"), serde_json::to_string_pretty(&items).unwrap()).map_err(|e| e.to_string())?;
    Ok(())
}
