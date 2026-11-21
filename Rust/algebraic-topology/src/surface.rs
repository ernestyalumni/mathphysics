//! Topology of surfaces made of squares (for example P-vortex sheets on the dual lattice): edge-connected
//! components, Euler characteristic chi = V - E + F, closedness (every edge in an even number of squares, a Z_2
//! cycle), branching edges (in 4 or more squares, which occur for Z_2 surfaces in 4D), and orientability of each
//! component by propagating square orientations across edges shared by exactly two squares.
//! Squares are given as 4 vertex ids in cyclic order; distinct squares must not share a vertex set.

use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq)]
pub struct ComponentTopology {
    pub squares: usize,
    pub vertices: usize,
    pub edges: usize,
    pub euler: i64,
    pub closed: bool,
    pub branching_edges: usize,
    /// Some(true/false) when every edge has exactly two squares; None when branching makes orientation ambiguous.
    pub orientable: Option<bool>,
}

fn find(p: &mut Vec<usize>, x: usize) -> usize {
    let mut r = x;
    while p[r] != r { r = p[r]; }
    let mut y = x;
    while p[y] != r { let n = p[y]; p[y] = r; y = n; }
    r
}

/// Topology of every edge-connected component, largest first.
pub fn analyze(squares: &[[usize; 4]]) -> Vec<ComponentTopology> {
    // edge -> list of (square, direction +1 if square traverses (min -> max), -1 otherwise)
    let mut edges: HashMap<(usize, usize), Vec<(usize, i8)>> = HashMap::new();
    for (s, q) in squares.iter().enumerate() {
        for i in 0..4 {
            let (a, b) = (q[i], q[(i + 1) % 4]);
            let (key, dir) = if a < b { ((a, b), 1) } else { ((b, a), -1) };
            edges.entry(key).or_default().push((s, dir));
        }
    }
    let mut parent: Vec<usize> = (0..squares.len()).collect();
    for list in edges.values() {
        for w in list.windows(2) {
            let (a, b) = (find(&mut parent, w[0].0), find(&mut parent, w[1].0));
            if a != b { parent[a] = b; }
        }
    }
    let mut comp: HashMap<usize, Vec<usize>> = HashMap::new();
    for s in 0..squares.len() { let r = find(&mut parent, s); comp.entry(r).or_default().push(s); }
    let mut out: Vec<ComponentTopology> = comp.values().map(|members| {
        let set: std::collections::HashSet<usize> = members.iter().cloned().collect();
        let mut verts = std::collections::HashSet::new();
        let mut my_edges = Vec::new();
        for &s in members { for &v in &squares[s] { verts.insert(v); } }
        for (k, list) in &edges { if set.contains(&list[0].0) { my_edges.push((*k, list)); } }
        let closed = my_edges.iter().all(|(_, l)| l.len() % 2 == 0);
        let branching = my_edges.iter().filter(|(_, l)| l.len() > 2).count();
        let orientable = if branching > 0 || !closed { None } else {
            // flip[s] in {+1, -1}; consistent iff across each edge the two squares traverse it oppositely after flips
            let mut flip: HashMap<usize, i8> = HashMap::new();
            let mut adj: HashMap<usize, Vec<(usize, i8)>> = HashMap::new(); // neighbour, required relative flip
            for (_, l) in &my_edges {
                let ((s1, d1), (s2, d2)) = (l[0], l[1]);
                let rel = -d1 * d2; // flip[s2] must equal rel * flip[s1] for opposite traversal: d1 f1 = -d2 f2
                adj.entry(s1).or_default().push((s2, rel));
                adj.entry(s2).or_default().push((s1, rel));
            }
            let mut ok = true;
            let start = members[0];
            flip.insert(start, 1);
            let mut stack = vec![start];
            while let Some(s) = stack.pop() {
                let fs = flip[&s];
                for &(t, rel) in adj.get(&s).map(|v| v.as_slice()).unwrap_or(&[]) {
                    let want = rel * fs;
                    match flip.get(&t) {
                        Some(&ft) => if ft != want { ok = false; },
                        None => { flip.insert(t, want); stack.push(t); }
                    }
                }
            }
            Some(ok)
        };
        ComponentTopology {
            squares: members.len(), vertices: verts.len(), edges: my_edges.len(),
            euler: verts.len() as i64 - my_edges.len() as i64 + members.len() as i64,
            closed, branching_edges: branching, orientable,
        }
    }).collect();
    out.sort_by(|a, b| b.squares.cmp(&a.squares));
    out
}
