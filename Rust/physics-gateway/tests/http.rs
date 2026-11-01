use serde_json::{json, Value};
use std::{
    fs,
    io::{BufRead, BufReader, Read, Write},
    net::TcpStream,
    path::PathBuf,
    process::{Child, Command, Stdio},
    time::Duration,
};
struct Fixture {
    child: Child,
    dir: PathBuf,
    port: u16,
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
        let _ = fs::remove_dir_all(&self.dir);
    }
}
impl Fixture {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "gateway-http-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(dir.join("static")).unwrap();
        fs::create_dir_all(dir.join("root/sub")).unwrap();
        fs::create_dir_all(dir.join("root/.git")).unwrap();
        fs::write(dir.join("static/index.html"), "app").unwrap();
        fs::write(dir.join("root/sub/note.tex"), "\\alpha+1").unwrap();
        fs::write(dir.join("root/sub/note.pdf"), b"%PDF-123456789").unwrap();
        fs::write(dir.join("root/.hidden.md"), "hidden").unwrap();
        fs::write(dir.join("root/bad.txt"), [255]).unwrap();
        fs::write(dir.join("root/page.html"), "<script>bad</script>").unwrap();
        fs::write(dir.join("root/.git/secret.tex"), "excluded").unwrap();
        fs::write(dir.join("outside.txt"), "outside").unwrap();
        std::os::unix::fs::symlink(dir.join("outside.txt"), dir.join("root/escape.txt")).unwrap();
        std::os::unix::fs::symlink(dir.join("root/sub"), dir.join("root/alias")).unwrap();
        std::os::unix::fs::symlink(dir.join("outside.txt"), dir.join("static/escape")).unwrap();
        let mut child = Command::new(env!("CARGO_BIN_EXE_physics-gateway"))
            .args([
                "--static-dir",
                dir.join("static").to_str().unwrap(),
                "--root",
                &format!("docs={}", dir.join("root").display()),
                "--port",
                "0",
            ])
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let mut line = String::new();
        BufReader::new(child.stdout.take().unwrap())
            .read_line(&mut line)
            .unwrap();
        let port = line.trim().rsplit(':').next().unwrap().parse().unwrap();
        Self { child, dir, port }
    }
    fn raw(&self, method: &str, path: &str, headers: &str, body: &[u8]) -> (u16, String, Vec<u8>) {
        let mut stream = TcpStream::connect(("127.0.0.1", self.port)).unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();
        write!(stream,"{method} {path} HTTP/1.1\r\nHost: 127.0.0.1:{}\r\nConnection: close\r\n{headers}Content-Length: {}\r\n\r\n",self.port,body.len()).unwrap();
        stream.write_all(body).unwrap();
        let mut bytes = vec![];
        stream.read_to_end(&mut bytes).unwrap();
        let i = bytes.windows(4).position(|w| w == b"\r\n\r\n").unwrap();
        let h = String::from_utf8_lossy(&bytes[..i]).to_string();
        let status = h.split_whitespace().nth(1).unwrap().parse().unwrap();
        (status, h, bytes[i + 4..].to_vec())
    }
    fn get(&self, path: &str) -> (u16, Value) {
        let (s, _, b) = self.raw("GET", path, "", &[]);
        (s, serde_json::from_slice(&b).unwrap())
    }
}
#[test]
fn filesystem_and_security() {
    let f = Fixture::new();
    assert_eq!(f.get("/api/health").1["status"], "ok");
    assert_eq!(f.get("/api/roots").1["roots"][0]["id"], "docs");
    let (_, v) = f.get("/api/files?root=docs");
    assert!(v["entries"]
        .as_array()
        .unwrap()
        .iter()
        .any(|e| e["name"] == ".hidden.md"));
    assert_eq!(v["entries"][0]["kind"], "directory");
    let (_, v) = f.get("/api/files?root=docs&q=note");
    assert!(v["entries"]
        .as_array()
        .unwrap()
        .iter()
        .any(|e| e["path"] == "alias/note.tex" || e["path"] == "sub/note.tex"));
    assert_eq!(
        f.get("/api/files?root=docs&q=secret").1["entries"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    let (s, v) = f.get("/api/text?root=docs&path=alias/note.tex");
    assert_eq!(s, 200);
    assert_eq!(v["path"], "alias/note.tex");
    assert_eq!(v["paired"]["path"], "alias/note.pdf");
    assert_eq!(v["text"], "\\alpha+1");
    for (p, s) in [
        ("/api/files?root=missing", 404),
        ("/api/text?root=docs&path=%2E%2E/outside.txt", 403),
        ("/api/text?root=docs&path=escape.txt", 403),
        ("/api/text?root=docs&path=page.html", 415),
        ("/api/text?root=docs&path=bad.txt", 422),
        ("/api/text?root=docs&path=/etc/passwd", 403),
        ("/%2e%2e/outside.txt", 403),
        ("/escape", 403),
        ("/api/text?root=docs&path=%zz", 400),
    ] {
        assert_eq!(f.get(p).0, s, "{p}");
    }
    assert_eq!(
        f.raw(
            "GET",
            "/api/health",
            "Origin: https://evil.invalid\r\n",
            &[]
        )
        .0,
        403
    );
    assert_eq!(
        f.raw("GET", "/api/health", "Host: evil.invalid\r\n", &[]).0,
        403
    );
    let (s, _, b) = f.raw("GET", "/", "", &[]);
    assert_eq!(s, 200);
    assert_eq!(b, b"app");
    let (_, v) = f.get("/api/harnesses");
    assert_eq!(v["harnesses"].as_array().unwrap().len(), 5);
    for h in v["harnesses"].as_array().unwrap() {
        assert!(h["status"] == "not_connected" || h["status"] == "not_found");
    }
}
#[test]
fn pdf_ranges() {
    let f = Fixture::new();
    let p = "/api/pdf?root=docs&path=sub/note.pdf";
    let (s, h, b) = f.raw("GET", p, "Range: bytes=2-5\r\n", &[]);
    assert_eq!(s, 206);
    assert_eq!(b, b"DF-1");
    assert!(h.to_lowercase().contains("content-range: bytes 2-5/14"));
    let (s, _, b) = f.raw("GET", p, "Range: bytes=-3\r\n", &[]);
    assert_eq!(s, 206);
    assert_eq!(b, b"789");
    assert_eq!(f.raw("GET", p, "Range: bytes=99-100\r\n", &[]).0, 416);
    assert_eq!(f.raw("GET", p, "Range: bytes=2-1\r\n", &[]).0, 416);
    assert_eq!(f.raw("GET", p, "Range: bytes=0-1,3-4\r\n", &[]).0, 416);
    let (s, h, b) = f.raw("HEAD", p, "", &[]);
    assert_eq!(s, 200);
    assert!(h.to_lowercase().contains("content-length: 14"));
    assert!(b.is_empty());
}
#[test]
fn limits_and_real_tool() {
    let f = Fixture::new();
    let file = fs::File::create(f.dir.join("root/large.txt")).unwrap();
    file.set_len(2_097_153).unwrap();
    assert_eq!(f.get("/api/text?root=docs&path=large.txt").0, 413);
    // Headers alone must reject an oversized request before waiting for its payload.
    let mut s = TcpStream::connect(("127.0.0.1", f.port)).unwrap();
    s.set_read_timeout(Some(Duration::from_secs(3))).unwrap();
    write!(s,"POST /api/tool HTTP/1.1\r\nHost: 127.0.0.1:{}\r\nContent-Length: 1048577\r\nConnection: close\r\n\r\n",f.port).unwrap();
    let mut b = [0; 1024];
    let n = s.read(&mut b).unwrap();
    assert!(String::from_utf8_lossy(&b[..n]).contains("413"));
    let body=json!({"command":"gamma_trace","labels":["x","y"],"domain":"formal_symmetric_h","budget":1}).to_string();
    let (s, _, b) = f.raw(
        "POST",
        "/api/tool",
        "Content-Type: application/json\r\n",
        body.as_bytes(),
    );
    assert_eq!(s, 200);
    let v: Value = serde_json::from_slice(&b).unwrap();
    assert_eq!(v["status"], "computed");
    assert_eq!(v["formal_proof"], false);
    assert_eq!(v["evidence_kind"], "exact_algebra");
    assert_eq!(v, physics_tools::handle_bytes(body.as_bytes()).response);
    let (s, _, b) = f.raw("POST", "/api/tool", "", br#"{"command":"nope"}"#);
    assert_eq!(s, 200);
    assert_ne!(
        serde_json::from_slice::<Value>(&b).unwrap()["status"],
        "computed"
    );
    assert_eq!(f.raw("POST", "/api/tool", "", b"{").0, 400);
    fs::create_dir(f.dir.join("root/many")).unwrap();
    for i in 0..205 {
        fs::write(f.dir.join(format!("root/many/{i}.txt")), "").unwrap();
    }
    let (_, v) = f.get("/api/files?root=docs&path=many");
    assert_eq!(v["entries"].as_array().unwrap().len(), 200);
    assert_eq!(v["truncated"], true);
}

#[test]
fn encoded_names_pairs_and_pdf_whitelist() {
    let f = Fixture::new();
    fs::write(f.dir.join("root/A B.TEX"), "source").unwrap();
    fs::write(f.dir.join("root/A B.pdf"), "%PDF-pair").unwrap();
    let (status, v) = f.get("/api/text?root=docs&path=A%20B.TEX");
    assert_eq!(status, 200);
    assert_eq!(v["paired"]["path"], "A B.pdf");
    fs::write(f.dir.join("root/UPPER.PDF"), b"%PDF-UPPER").unwrap();
    assert_eq!(
        f.raw("GET", "/api/pdf?root=docs&path=UPPER.PDF", "", &[]).0,
        200
    );
    std::os::unix::fs::symlink(f.dir.join("root/UPPER.PDF"), f.dir.join("root/alias.txt")).unwrap();
    assert_eq!(f.get("/api/pdf?root=docs&path=alias.txt").0, 415);
    std::os::unix::fs::symlink(f.dir.join("outside.txt"), f.dir.join("root/escape.pdf")).unwrap();
    fs::write(f.dir.join("root/escape.tex"), "source").unwrap();
    assert_eq!(
        f.get("/api/text?root=docs&path=escape.tex").1["paired"],
        Value::Null
    );
    assert_eq!(f.get("/api/pdf?root=docs&path=escape.pdf").0, 403);
    assert_eq!(f.get("/api/files?root=docs&root=other").0, 400);
    let origin = format!("Origin: http://127.0.0.1:{}\r\n", f.port);
    assert_eq!(f.raw("GET", "/api/health", &origin, &[]).0, 200);
    let (status, headers, body) = f.raw(
        "HEAD",
        "/api/pdf?root=docs&path=UPPER.PDF",
        "Range: bytes=5-\r\n",
        &[],
    );
    assert_eq!(status, 206);
    assert!(headers
        .to_lowercase()
        .contains("content-range: bytes 5-9/10"));
    assert!(body.is_empty());
    fs::write(f.dir.join("root/empty.pdf"), []).unwrap();
    assert_eq!(
        f.raw(
            "GET",
            "/api/pdf?root=docs&path=empty.pdf",
            "Range: bytes=0-\r\n",
            &[]
        )
        .0,
        416
    );
}

#[test]
fn parsed_source_extensions_and_same_origin_frames() {
    let f = Fixture::new();
    let source = "# Parsed α\n<script>alert('source only')</script>\n\\alpha";
    for (extension, kind) in [
        ("mmd", "markdown"),
        ("jsonl", "text"),
        ("yaml", "text"),
        ("yml", "text"),
        ("csv", "text"),
    ] {
        fs::write(f.dir.join(format!("root/parsed.{extension}")), source).unwrap();
        let (status, headers, body) = f.raw(
            "GET",
            &format!("/api/text?root=docs&path=parsed.{extension}"),
            "",
            &[],
        );
        assert_eq!(status, 200);
        assert!(headers
            .to_lowercase()
            .contains("content-type: application/json"));
        let value: Value = serde_json::from_slice(&body).unwrap();
        assert_eq!(value["text"], source);
        assert_eq!(value["kind"], kind);
        assert_eq!(value["paired"], Value::Null);
    }
    fs::write(f.dir.join("root/invalid.mmd"), [255]).unwrap();
    assert_eq!(f.get("/api/text?root=docs&path=invalid.mmd").0, 422);
    assert_eq!(f.get("/api/text?root=docs&path=page.html").0, 415);
    assert_eq!(
        f.get("/api/text?root=docs&path=sub/note.tex").1["paired"]["path"],
        "sub/note.pdf"
    );
    for path in ["/", "/api/pdf?root=docs&path=sub/note.pdf"] {
        let (status, headers, _) = f.raw("GET", path, "", &[]);
        assert_eq!(status, 200);
        assert!(headers
            .to_lowercase()
            .contains("x-frame-options: sameorigin"));
        assert!(!headers.to_lowercase().contains("x-frame-options: deny"));
    }
}

// Independent oracle: explicit Pauli entries and ordinary complex matrix products.
type C = [f64; 2];
type M = [[C; 2]; 2];
fn c_add(a: C, b: C) -> C {
    [a[0] + b[0], a[1] + b[1]]
}
fn c_mul(a: C, b: C) -> C {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn adjoint(a: M) -> M {
    std::array::from_fn(|i| std::array::from_fn(|j| [a[j][i][0], -a[j][i][1]]))
}
fn product(a: M, b: M) -> M {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| c_add(c_mul(a[i][0], b[0][j]), c_mul(a[i][1], b[1][j])))
    })
}
fn identity() -> M {
    [[[1., 0.], [0., 0.]], [[0., 0.], [1., 0.]]]
}
fn pauli() -> [M; 3] {
    [
        [[[0., 0.], [1., 0.]], [[1., 0.], [0., 0.]]],
        [[[0., 0.], [0., -1.]], [[0., 1.], [0., 0.]]],
        [[[1., 0.], [0., 0.]], [[0., 0.], [-1., 0.]]],
    ]
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11, "{a} != {b}");
}
fn matrix_close(a: M, b: M) {
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                close(a[i][j][k], b[i][j][k]);
            }
        }
    }
}
fn oracle(axis: [f64; 3], t: f64) -> M {
    let mut out = identity();
    for (i, row) in out.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = c_mul(*entry, [(t / 2.).cos(), 0.]);
            for (k, n) in axis.iter().enumerate() {
                *entry = c_add(*entry, c_mul([0., -n * (t / 2.).sin()], pauli()[k][i][j]));
            }
        }
    }
    out
}
fn lab(f: &Fixture, axis: [f64; 3], t: f64, slice: f64) -> Value {
    let body = json!({"axis":axis,"t":t,"slice":slice}).to_string();
    let (status, headers, bytes) = f.raw(
        "POST",
        "/api/su2",
        "Content-Type: application/json\r\n",
        body.as_bytes(),
    );
    assert_eq!(status, 200, "{}", String::from_utf8_lossy(&bytes));
    assert!(headers
        .to_lowercase()
        .contains("content-type: application/json"));
    let v: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(v["convention"], "pauli-plus-i-v1");
    assert_eq!(v["evidence_kind"], "numerical");
    assert_eq!(v["formal_proof"], false);
    assert_eq!(v["axis"], json!(axis));
    close(v["t"].as_f64().unwrap(), t);
    v
}
fn response_matrix(v: &Value) -> M {
    serde_json::from_value(v["matrix"].clone()).unwrap()
}

#[test]
fn su2_direct_pauli_matrix_unitarity_determinant_and_bloch() {
    let f = Fixture::new();
    let axes = [
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [1. / 3., 2. / 3., -2. / 3.],
        [-0.6, 0., 0.8],
    ];
    let pi = std::f64::consts::PI;
    for axis in axes {
        for t in [0., pi / 2., pi, 2. * pi, 4. * pi, -pi, 1e6, -1e6, 0.713] {
            let v = lab(&f, axis, t, 0.25);
            let m = response_matrix(&v);
            matrix_close(m, oracle(axis, t));
            matrix_close(product(adjoint(m), m), identity());
            let d = c_add(
                c_mul(m[0][0], m[1][1]),
                c_mul([-1., 0.], c_mul(m[0][1], m[1][0])),
            );
            close(d[0], 1.);
            close(d[1], 0.);
            // Transform +x density matrix, then trace against all three Paulis.
            let rho = [[[0.5, 0.], [0.5, 0.]], [[0.5, 0.], [0.5, 0.]]];
            let evolved = product(product(m, rho), adjoint(m));
            for (k, sigma) in pauli().into_iter().enumerate() {
                let p = product(evolved, sigma);
                let tr = c_add(p[0][0], p[1][1]);
                close(v["bloch"][k].as_f64().unwrap(), tr[0]);
                close(tr[1], 0.);
            }
            let coordinates: [f64; 4] = serde_json::from_value(v["coordinates"].clone()).unwrap();
            close(coordinates.iter().map(|x| x * x).sum(), 1.);
            assert_eq!(v["alpha"], json!([coordinates[0], coordinates[3]]));
            assert_eq!(v["beta"], json!([coordinates[2], coordinates[1]]));
            assert!(v["checks"]["unitarity_error"].as_f64().unwrap() < 3e-12);
            assert!(v["checks"]["determinant_error"].as_f64().unwrap() < 3e-12);
        }
    }
    let z = lab(&f, [0., 0., 1.], pi / 2., 0.);
    close(z["bloch"][1].as_f64().unwrap(), 1.);
    let u0 = lab(&f, [0., 0., 1.], 0., 0.);
    let u2 = lab(&f, [0., 0., 1.], 2. * pi, 0.);
    let u4 = lab(&f, [0., 0., 1.], 4. * pi, 0.);
    close(u0["coordinates"][0].as_f64().unwrap(), 1.);
    close(u2["coordinates"][0].as_f64().unwrap(), -1.);
    close(u4["coordinates"][0].as_f64().unwrap(), 1.);
    // Order/sign oracle: (-i sigma_x)(-i sigma_y)=-i sigma_z.
    let x = response_matrix(&lab(&f, [1., 0., 0.], pi, 0.));
    let y = response_matrix(&lab(&f, [0., 1., 0.], pi, 0.));
    let z = response_matrix(&lab(&f, [0., 0., 1.], pi, 0.));
    matrix_close(product(x, y), z);
    matrix_close(product(y, x), oracle([0., 0., -1.], pi));
    assert_eq!(
        lab(&f, [1. / 3., 2. / 3., -2. / 3.], 0.713, 0.25),
        lab(&f, [1. / 3., 2. / 3., -2. / 3.], 0.713, 0.25)
    );
}

#[test]
fn su2_chart_and_slice_edges() {
    let f = Fixture::new();
    let pi = std::f64::consts::PI;
    for t in [0., 4. * pi, -4. * pi] {
        let v = lab(&f, [0., 0., 1.], t, 0.);
        assert_eq!(v["stereographic"], json!({"status":"pole","point":null}));
    }
    let south = lab(&f, [0., 0., 1.], 2. * pi, 0.);
    assert_eq!(south["stereographic"]["status"], "finite");
    close(south["stereographic"]["point"][2].as_f64().unwrap(), 0.);
    for t in [pi, 0.71, -0.71, 1e-10] {
        let v = lab(&f, [0., 0., 1.], t, 0.);
        assert_eq!(v["stereographic"]["status"], "finite");
        let p: [f64; 3] = serde_json::from_value(v["stereographic"]["point"].clone()).unwrap();
        assert!((p[2] / (-(t / 4.).cos() / (t / 4.).sin()) - 1.).abs() < 1e-12);
        if t.abs() > 1e-5 {
            close(
                p[2],
                v["coordinates"][3].as_f64().unwrap()
                    / (1. - v["coordinates"][0].as_f64().unwrap()),
            );
        }
    }
    for c in [
        -1.2,
        -1.,
        -1. + f64::EPSILON,
        0.,
        0.6,
        1. - f64::EPSILON,
        1.,
        1.2,
        f64::MAX,
    ] {
        let v = lab(&f, [1., 0., 0.], 0.5, c);
        close(v["slice"]["coordinate"].as_f64().unwrap(), c);
        if c.abs() > 1. {
            assert_eq!(v["slice"]["status"], "empty");
            assert!(v["slice"]["radius"].is_null());
        } else if c.abs() == 1. {
            assert_eq!(v["slice"]["status"], "point");
            close(v["slice"]["radius"].as_f64().unwrap(), 0.);
        } else {
            assert_eq!(v["slice"]["status"], "sphere");
            let r = v["slice"]["radius"].as_f64().unwrap();
            close(r * r + c * c, 1.);
            assert!(r > 0.);
        }
    }
}

#[test]
fn su2_invalid_requests_and_inherited_http_guards() {
    let f = Fixture::new();
    for body in [
        "{",
        "null",
        "[]",
        "{}",
        r#"{"axis":[0,0,1],"t":0}"#,
        r#"{"axis":[0,0,1],"t":0,"slice":0,"extra":1}"#,
        r#"{"axis":[0,0],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,1,0],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,0],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,2],"t":0,"slice":0}"#,
        r#"{"axis":[1,1,1],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,1.0000001],"t":0,"slice":0}"#,
        r#"{"axis":[1e308,1e308,1e308],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,NaN],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,1e999],"t":0,"slice":0}"#,
        r#"{"axis":[0,0,"1"],"t":0,"slice":0}"#,
        r#"{"axis":null,"t":0,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":null,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":true,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":1000000.1,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":-1000000.1,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":Infinity,"slice":0}"#,
        r#"{"axis":[0,0,1],"t":0,"slice":1e999}"#,
        r#"{"axis":[0,0,1],"t":0,"slice":[]}"#,
    ] {
        let (status, _, bytes) = f.raw("POST", "/api/su2", "", body.as_bytes());
        assert_eq!(status, 400, "{body}");
        let v: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(v["error"]["code"], 400);
    }
    let body = br#"{"axis":[0,0,1],"t":0,"slice":0}"#;
    for method in ["GET", "HEAD", "PUT", "DELETE"] {
        assert_eq!(f.raw(method, "/api/su2", "", body).0, 405);
    }
    for headers in [
        "Origin: https://evil.invalid\r\n",
        "Origin: null\r\n",
        "Host: evil.invalid\r\n",
    ] {
        assert_eq!(f.raw("POST", "/api/su2", headers, body).0, 403);
    }
    assert_eq!(
        f.raw(
            "POST",
            "/api/su2",
            &format!("Origin: http://127.0.0.1:{}\r\n", f.port),
            body
        )
        .0,
        200
    );
    let mut stream = TcpStream::connect(("127.0.0.1", f.port)).unwrap();
    stream
        .set_read_timeout(Some(Duration::from_secs(5)))
        .unwrap();
    write!(stream,"POST /api/su2 HTTP/1.1\r\nHost: 127.0.0.1:{}\r\nContent-Length: 1048577\r\nConnection: close\r\n\r\n",f.port).unwrap();
    let mut response = String::new();
    stream.read_to_string(&mut response).unwrap();
    assert!(response.starts_with("HTTP/1.1 413"));
    // Chunked requests have no declared length; exercise the limited body reader.
    let mut stream = TcpStream::connect(("127.0.0.1", f.port)).unwrap();
    stream
        .set_read_timeout(Some(Duration::from_secs(5)))
        .unwrap();
    write!(stream,"POST /api/su2 HTTP/1.1\r\nHost: 127.0.0.1:{}\r\nTransfer-Encoding: chunked\r\nConnection: close\r\n\r\n100001\r\n",f.port).unwrap();
    stream.write_all(&vec![b' '; 1_048_577]).unwrap();
    stream.write_all(b"\r\n0\r\n\r\n").unwrap();
    let mut response = String::new();
    stream.read_to_string(&mut response).unwrap();
    assert!(response.starts_with("HTTP/1.1 413"));
}
