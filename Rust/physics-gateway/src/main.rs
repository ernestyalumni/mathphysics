mod su2;

use serde_json::{json, Value};
use std::{
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{Read, Seek, SeekFrom},
    path::{Component, Path, PathBuf},
    sync::Arc,
    time::{Instant, UNIX_EPOCH},
};
use tiny_http::{Header, Request, Response, Server, StatusCode};
const BODY: u64 = 1_048_576;
const TEXT: u64 = 2_097_152;
type Err = (u16, String);
type Result<T> = std::result::Result<T, Err>;
struct Config {
    static_dir: PathBuf,
    roots: HashMap<String, PathBuf>,
    port: u16,
}
fn err(code: u16, s: impl ToString) -> Err {
    (code, s.to_string())
}
fn header(k: &str, v: &str) -> Header {
    Header::from_bytes(k, v).unwrap()
}
fn send(req: Request, status: u16, value: Value) {
    let _ = req.respond(
        Response::from_string(value.to_string())
            .with_status_code(StatusCode(status))
            .with_header(header("Content-Type", "application/json"))
            .with_header(header("X-Content-Type-Options", "nosniff")),
    );
}
fn confined(base: &Path, rel: &str) -> Result<PathBuf> {
    if rel.contains('\\')
        || rel.contains('\0')
        || Path::new(rel).components().any(|c| {
            matches!(
                c,
                Component::ParentDir | Component::RootDir | Component::Prefix(_)
            )
        })
    {
        return Err(err(403, "unsafe_path"));
    }
    let p = base.join(rel).canonicalize().map_err(|e| {
        err(
            if e.kind() == std::io::ErrorKind::NotFound {
                404
            } else {
                403
            },
            e,
        )
    })?;
    if !p.starts_with(base) {
        return Err(err(403, "outside_root"));
    }
    Ok(p)
}
fn kind(p: &Path) -> &'static str {
    if p.is_dir() {
        return "directory";
    }
    extension_kind(p)
}
fn extension_kind(p: &Path) -> &'static str {
    match p
        .extension()
        .and_then(|x| x.to_str())
        .unwrap_or("")
        .to_ascii_lowercase()
        .as_str()
    {
        "pdf" => "pdf",
        "tex" => "tex",
        "md" | "mmd" => "markdown",
        "txt" | "json" | "jsonl" | "yaml" | "yml" | "csv" | "prompt" | "rs" | "toml" => "text",
        _ => "file",
    }
}
fn modified(m: &fs::Metadata) -> u64 {
    m.modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_secs())
        .unwrap_or(0)
}
fn skipped(s: &str) -> bool {
    matches!(
        s,
        ".git"
            | "node_modules"
            | "target"
            | "build"
            | ".cache"
            | "__pycache__"
            | "dist"
            | ".next"
            | ".venv"
    )
}
fn files(base: &Path, rel: &str, q: &str) -> Result<Value> {
    let start = confined(base, rel)?;
    if !start.is_dir() {
        return Err(err(400, "not_directory"));
    }
    let clock = Instant::now();
    let mut stack = vec![(start, rel.to_string())];
    let mut seen = HashSet::new();
    let mut entries = vec![];
    let mut visited = 0;
    let mut truncated = false;
    'scan: while let Some((dir, display)) = stack.pop() {
        if !seen.insert(dir.clone()) {
            continue;
        }
        let iter = fs::read_dir(dir).map_err(|e| err(403, e))?;
        for item in iter {
            if visited >= 10000 || entries.len() >= 200 || clock.elapsed().as_millis() >= 500 {
                truncated = true;
                break 'scan;
            }
            visited += 1;
            let item = match item {
                Ok(x) => x,
                Err(_) => continue,
            };
            let name = item.file_name().to_string_lossy().to_string();
            let path = if display.is_empty() {
                name.clone()
            } else {
                format!("{display}/{name}")
            };
            let actual = match confined(base, &path) {
                Ok(x) => x,
                Err(_) => continue,
            };
            let meta = match fs::metadata(&actual) {
                Ok(x) => x,
                Err(_) => continue,
            };
            if !q.is_empty() && meta.is_dir() && !skipped(&name) {
                stack.push((actual.clone(), path.clone()));
            }
            if q.is_empty() || name.to_lowercase().contains(&q.to_lowercase()) {
                entries.push(json!({"name":name,"path":path,"kind":kind(&actual),"size":meta.len(),"modified_secs":modified(&meta)}));
            }
        }
        if q.is_empty() {
            break;
        }
    }
    entries.sort_by(|a, b| {
        let ad = a["kind"] == "directory";
        let bd = b["kind"] == "directory";
        bd.cmp(&ad)
            .then_with(|| a["name"].as_str().cmp(&b["name"].as_str()))
    });
    Ok(json!({"path":rel,"entries":entries,"truncated":truncated,"visited":visited}))
}
fn range(s: &str, len: u64) -> Result<(u64, u64)> {
    let s = s
        .strip_prefix("bytes=")
        .ok_or_else(|| err(416, "invalid_range"))?;
    let (a, b) = s.split_once('-').ok_or_else(|| err(416, "invalid_range"))?;
    if len == 0 {
        return Err(err(416, "empty_file"));
    }
    if a.is_empty() {
        let n: u64 = b.parse().map_err(|_| err(416, "invalid_range"))?;
        if n == 0 {
            return Err(err(416, "invalid_range"));
        }
        return Ok((len.saturating_sub(n), len - 1));
    }
    let start: u64 = a.parse().map_err(|_| err(416, "invalid_range"))?;
    let end = if b.is_empty() {
        len - 1
    } else {
        b.parse::<u64>()
            .map_err(|_| err(416, "invalid_range"))?
            .min(len - 1)
    };
    if start >= len || end < start {
        return Err(err(416, "invalid_range"));
    }
    Ok((start, end))
}
fn stream(req: Request, p: PathBuf, pdf: bool) -> Result<()> {
    let mut f = match File::open(&p) {
        Ok(f) => f,
        Err(e) => {
            send(
                req,
                403,
                json!({"error":{"code":403,"message":e.to_string()}}),
            );
            return Ok(());
        }
    };
    let len = match f.metadata() {
        Ok(m) => m.len(),
        Err(e) => {
            send(
                req,
                403,
                json!({"error":{"code":403,"message":e.to_string()}}),
            );
            return Ok(());
        }
    };
    let r = req
        .headers()
        .iter()
        .find(|h| h.field.equiv("Range"))
        .map(|h| h.value.as_str());
    let (start, end, status) = if pdf {
        if let Some(r) = r {
            match range(r, len) {
                Ok((a, b)) => (a, b, 206),
                Err(_) => {
                    let _ = req.respond(
                        Response::from_string(
                            json!({"error":{"code":416,"message":"invalid_range"}}).to_string(),
                        )
                        .with_status_code(416)
                        .with_header(header("Content-Range", &format!("bytes */{len}")))
                        .with_header(header("Content-Type", "application/json")),
                    );
                    return Ok(());
                }
            }
        } else {
            (0, len.saturating_sub(1), 200)
        }
    } else {
        (0, len.saturating_sub(1), 200)
    };
    let size = if len == 0 { 0 } else { end - start + 1 };
    if let Err(e) = f.seek(SeekFrom::Start(start)) {
        send(
            req,
            500,
            json!({"error":{"code":500,"message":e.to_string()}}),
        );
        return Ok(());
    }
    let mime = if pdf {
        "application/pdf".into()
    } else {
        mime_guess::from_path(&p)
            .first_or_octet_stream()
            .to_string()
    };
    let mut headers = vec![
        header("Content-Type", &mime),
        header("X-Content-Type-Options", "nosniff"),
        header("X-Frame-Options", "SAMEORIGIN"),
    ];
    if pdf {
        headers.push(header("Content-Disposition", "inline"));
        headers.push(header("Accept-Ranges", "bytes"));
        if status == 206 {
            headers.push(header(
                "Content-Range",
                &format!("bytes {start}-{end}/{len}"),
            ));
        }
    }
    let response = Response::new(
        StatusCode(status),
        headers,
        f.take(size),
        Some(size as usize),
        None,
    );
    let _ = req.respond(response);
    Ok(())
}
fn decode(s: &str, plus: bool) -> Result<String> {
    let b = s.as_bytes();
    for i in 0..b.len() {
        if b[i] == b'%'
            && (i + 2 >= b.len() || !b[i + 1].is_ascii_hexdigit() || !b[i + 2].is_ascii_hexdigit())
        {
            return Err(err(400, "invalid_encoding"));
        }
    }
    let v = if plus {
        s.replace('+', " ")
    } else {
        s.to_string()
    };
    percent_encoding::percent_decode_str(&v)
        .decode_utf8()
        .map(|x| x.into_owned())
        .map_err(|_| err(400, "invalid_encoding"))
}
fn query_params(s: &str) -> Result<HashMap<String, String>> {
    let mut result = HashMap::new();
    for item in s.split('&').filter(|x| !x.is_empty()) {
        let (k, v) = item.split_once('=').unwrap_or((item, ""));
        let k = decode(k, true)?;
        if result.insert(k, decode(v, true)?).is_some() {
            return Err(err(400, "duplicate_parameter"));
        }
    }
    Ok(result)
}
fn handle(mut req: Request, c: &Config) {
    let hosts: Vec<_> = req
        .headers()
        .iter()
        .filter(|h| h.field.equiv("Host"))
        .collect();
    let valid = [
        format!("127.0.0.1:{}", c.port),
        format!("localhost:{}", c.port),
    ];
    if hosts.len() != 1
        || !valid
            .iter()
            .any(|v| v.eq_ignore_ascii_case(hosts[0].value.as_str()))
    {
        send(
            req,
            403,
            json!({"error":{"code":403,"message":"denied_host"}}),
        );
        return;
    }
    for h in req.headers() {
        if h.field.equiv("Origin")
            && h.value.as_str() != format!("http://{}", hosts[0].value.as_str())
        {
            send(
                req,
                403,
                json!({"error":{"code":403,"message":"denied_origin"}}),
            );
            return;
        }
    }
    if req.body_length().is_some_and(|n| n > BODY as usize) {
        send(
            req,
            413,
            json!({"error":{"code":413,"message":"body_too_large"}}),
        );
        return;
    }
    let mut body = vec![];
    if req
        .as_reader()
        .take(BODY + 1)
        .read_to_end(&mut body)
        .is_err()
    {
        send(
            req,
            400,
            json!({"error":{"code":400,"message":"body_read"}}),
        );
        return;
    }
    if body.len() > BODY as usize {
        send(
            req,
            413,
            json!({"error":{"code":413,"message":"body_too_large"}}),
        );
        return;
    }
    let raw = req.url().to_string();
    let method = req.method().as_str().to_string();
    // Decode raw paths before URL normalization so ../ cannot disappear.
    let (raw_path, query) = raw.split_once('?').unwrap_or((&raw, ""));
    let decoded = match decode(raw_path, false) {
        Ok(v) => v,
        Err((s, m)) => {
            send(req, s, json!({"error":{"code":s,"message":m}}));
            return;
        }
    };
    let params = match query_params(query) {
        Ok(v) => v,
        Err((s, m)) => {
            send(req, s, json!({"error":{"code":s,"message":m}}));
            return;
        }
    };
    let rel = params.get("path").map(String::as_str).unwrap_or("");
    let result = (|| -> Result<Option<Value>> {
        if method != "GET"
            && !(method == "POST" && matches!(decoded.as_str(), "/api/tool" | "/api/su2"))
            && !(method == "HEAD" && (decoded == "/api/pdf" || !decoded.starts_with("/api/")))
        {
            return Err(err(405, "method_not_allowed"));
        }
        match decoded.as_str() {
            "/api/health" => Ok(Some(
                json!({"status":"ok","version":env!("CARGO_PKG_VERSION"),"read_only":true}),
            )),
            "/api/roots" => {
                let mut roots: Vec<_> = c
                    .roots
                    .iter()
                    .map(|(id, p)| json!({"id":id,"label":id,"path":p}))
                    .collect();
                roots.sort_by_key(|v| v["id"].as_str().unwrap().to_string());
                Ok(Some(json!({"roots":roots})))
            }
            "/api/harnesses" => {
                use std::os::unix::fs::PermissionsExt;
                let list = [
                    ("claude-code", "ClaudeCode", "claude"),
                    ("codex", "Codex", "codex"),
                    ("grok-build", "GrokBuild", "grok"),
                    ("hermes-agent", "HermesAgent", "hermes"),
                    ("pi", "Pi", "pi"),
                ];
                let hs:Vec<_>=list.iter().map(|(id,name,exe)|{let available=std::env::var_os("PATH").map(|p|std::env::split_paths(&p).any(|d|fs::metadata(d.join(exe)).is_ok_and(|m|m.is_file()&&m.permissions().mode()&0o111!=0))).unwrap_or(false);json!({"id":id,"name":name,"available":available,"status":if available{"not_connected"}else{"not_found"},"detail":format!("PATH executable detection only: {exe}; no process started")})}).collect();
                Ok(Some(json!({"harnesses":hs})))
            }
            "/api/su2" => {
                if method != "POST" {
                    return Err(err(405, "method_not_allowed"));
                }
                Ok(Some(su2::handle_bytes(&body).map_err(|e| err(400, e))?))
            }
            "/api/tool" => {
                serde_json::from_slice::<Value>(&body).map_err(|e| err(400, e))?;
                Ok(Some(physics_tools::handle_bytes(&body).response))
            }
            "/api/files" | "/api/text" | "/api/pdf" => {
                let id = params.get("root").ok_or_else(|| err(400, "missing_root"))?;
                let base = c.roots.get(id).ok_or_else(|| err(404, "unknown_root"))?;
                if decoded == "/api/files" {
                    let mut v =
                        files(base, rel, params.get("q").map(String::as_str).unwrap_or(""))?;
                    v["root"] = json!(id);
                    return Ok(Some(v));
                }
                let p = confined(base, rel)?;
                let k = extension_kind(Path::new(rel));
                if decoded == "/api/pdf" {
                    if k != "pdf" || !p.is_file() {
                        return Err(err(415, "pdf_required"));
                    }
                    return Ok(None);
                }
                if !matches!(k, "tex" | "markdown" | "text") {
                    return Err(err(415, "unsupported_text_extension"));
                }
                let mut f = File::open(&p).map_err(|e| err(403, e))?;
                let meta = f.metadata().map_err(|e| err(403, e))?;
                if meta.len() > TEXT {
                    return Err(err(413, "text_too_large"));
                }
                let mut bytes = vec![];
                f.by_ref()
                    .take(TEXT + 1)
                    .read_to_end(&mut bytes)
                    .map_err(|e| err(500, e))?;
                if bytes.len() > TEXT as usize {
                    return Err(err(413, "text_too_large"));
                }
                let text = String::from_utf8(bytes).map_err(|_| err(422, "invalid_utf8"))?;
                let paired = if k == "tex" {
                    let candidate = Path::new(rel)
                        .with_extension("pdf")
                        .to_string_lossy()
                        .to_string();
                    if confined(base, &candidate).is_ok_and(|p| p.is_file()) {
                        json!({"path":candidate,"kind":"pdf"})
                    } else {
                        Value::Null
                    }
                } else {
                    Value::Null
                };
                Ok(Some(
                    json!({"root":id,"path":rel,"text":text,"kind":k,"paired":paired,"modified_secs":modified(&meta)}),
                ))
            }
            p if p.starts_with("/api/") => Err(err(404, "unknown_endpoint")),
            _ => {
                let p = confined(
                    &c.static_dir,
                    if decoded == "/" {
                        "index.html"
                    } else {
                        decoded.trim_start_matches('/')
                    },
                )?;
                if !p.is_file() {
                    return Err(err(404, "not_file"));
                }
                File::open(p).map_err(|e| err(403, e))?;
                Ok(None)
            }
        }
    })();
    match result {
        Ok(Some(v)) => send(req, 200, v),
        Err((s, m)) => send(req, s, json!({"error":{"code":s,"message":m}})),
        Ok(None) => {
            let pdf = decoded == "/api/pdf";
            let p = if pdf {
                confined(&c.roots[&params["root"]], rel)
            } else {
                confined(
                    &c.static_dir,
                    if decoded == "/" {
                        "index.html"
                    } else {
                        decoded.trim_start_matches('/')
                    },
                )
            };
            match p {
                Ok(p) => {
                    let _ = stream(req, p, pdf);
                }
                Err((s, m)) => send(req, s, json!({"error":{"code":s,"message":m}})),
            }
        }
    }
}
fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1)
    }
}
fn run() -> std::result::Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let mut args = std::env::args().skip(1);
    let mut static_dir = None;
    let mut roots = HashMap::new();
    let mut port = 8876;
    while let Some(a) = args.next() {
        let v = args.next().ok_or("missing argument value")?;
        match a.as_str() {
            "--static-dir" => {
                if !Path::new(&v).is_absolute() {
                    return Err("static directory must be absolute".into());
                }
                static_dir = Some(PathBuf::from(v).canonicalize()?);
            }
            "--port" => port = v.parse()?,
            "--root" => {
                let (id, p) = v.split_once('=').ok_or("root must be id=ABS")?;
                if id.is_empty()
                    || !id
                        .bytes()
                        .all(|x| x.is_ascii_alphanumeric() || x == b'_' || x == b'-')
                    || !Path::new(p).is_absolute()
                {
                    return Err("invalid root".into());
                }
                let p = PathBuf::from(p).canonicalize()?;
                if !p.is_dir() || roots.insert(id.to_string(), p).is_some() {
                    return Err("duplicate root or not directory".into());
                }
            }
            _ => return Err("unknown argument".into()),
        }
    }
    let static_dir = static_dir.ok_or("--static-dir required")?;
    if !static_dir.is_dir() {
        return Err("static directory required".into());
    }
    let server = Arc::new(Server::http(("127.0.0.1", port))?);
    let port = server
        .server_addr()
        .to_ip()
        .ok_or("missing address")?
        .port();
    let c = Arc::new(Config {
        static_dir,
        roots,
        port,
    });
    println!("http://127.0.0.1:{port}");
    let mut workers = vec![];
    for _ in 0..4 {
        let s = server.clone();
        let c = c.clone();
        workers.push(std::thread::spawn(move || {
            for req in s.incoming_requests() {
                handle(req, &c)
            }
        }));
    }
    for w in workers {
        let _ = w.join();
    }
    Ok(())
}
