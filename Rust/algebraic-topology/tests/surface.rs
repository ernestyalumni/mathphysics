//! Oracles: surface of a cube (sphere, chi 2, orientable); square torus and Klein bottle (chi 0; orientable / not).
use algebraic_topology::surface::analyze;

#[test]
fn cube_torus_klein() {
    // cube vertices: bits (x, y, z)
    let v = |x: usize, y: usize, z: usize| x + 2 * y + 4 * z;
    let cube = vec![
        [v(0,0,0), v(1,0,0), v(1,1,0), v(0,1,0)], [v(0,0,1), v(0,1,1), v(1,1,1), v(1,0,1)],
        [v(0,0,0), v(0,0,1), v(1,0,1), v(1,0,0)], [v(0,1,0), v(1,1,0), v(1,1,1), v(0,1,1)],
        [v(0,0,0), v(0,1,0), v(0,1,1), v(0,0,1)], [v(1,0,0), v(1,0,1), v(1,1,1), v(1,1,0)],
    ];
    let c = &analyze(&cube)[0];
    assert_eq!((c.euler, c.closed, c.orientable), (2, true, Some(true)));
    let n = 4;
    let grid = |flip: bool| -> Vec<[usize; 4]> {
        let id = |i: i64, j: i64| -> usize {
            let mut i = i;
            if j == n && flip { i = -i - 1; }
            (i.rem_euclid(n) + n * j.rem_euclid(n)) as usize
        };
        let mut sq = Vec::new();
        for i in 0..n { for j in 0..n { sq.push([id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]); } }
        sq
    };
    let t = &analyze(&grid(false))[0];
    assert_eq!((t.euler, t.closed, t.orientable, t.squares), (0, true, Some(true), 16));
    let k = &analyze(&grid(true))[0];
    assert_eq!((k.euler, k.closed, k.orientable), (0, true, Some(false)));
}

#[test]
fn homological_orientability_certificate() {
    use algebraic_topology::surface::chain_complex;
    let n = 4i64;
    let grid = |flip: bool| -> Vec<[usize; 4]> {
        let id = |i: i64, j: i64| -> usize { let mut i = i; if j == n && flip { i = -i - 1; } (i.rem_euclid(n) + n * j.rem_euclid(n)) as usize };
        let mut sq = Vec::new();
        for i in 0..n { for j in 0..n { sq.push([id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)]); } }
        sq
    };
    let t = chain_complex(&grid(false));
    assert!(t.is_complex());
    assert_eq!(t.betti_mod_p(2), vec![1, 2, 1]);
    assert_eq!(t.betti_mod_p(1_000_003), vec![1, 2, 1]);
    let k = chain_complex(&grid(true));
    assert_eq!(k.betti_mod_p(2), vec![1, 2, 1]);
    assert_eq!(k.betti_mod_p(1_000_003), vec![1, 1, 0], "Klein bottle: b2 over Q must vanish");
}
