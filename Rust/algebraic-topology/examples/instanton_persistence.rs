//! Instanton lumps by persistent homology of the cooled charge density q(x) (dump from lattice-gauge dump_charge).
//! Sublevel sets of -q: every H_0 class is a separate peak of positive charge; its persistence is the peak height
//! above the saddle that merges it with a higher peak. Likewise +q for anti-instantons. N_I and N_A count classes
//! with persistence above tau (the essential class counts when its peak exceeds tau). Consistency: N_I - N_A = Q.
//! Peak-height persistence is scale dependent (an instanton of size rho peaks at ~6/(pi^2 rho^4)), so the program also
//! segments q into steepest-ascent basins (each site climbs to a local maximum of q, or of -q) and counts basins whose
//! integrated charge exceeds 1/2 in magnitude: a scale-invariant lump count.
//! Usage: instanton_persistence <dump> [tau]
use algebraic_topology::cubical::CubicalTorus;
use algebraic_topology::persistence::{lower_star_cubical, persistence};
fn count(t: &CubicalTorus, f: &[f64], tau: f64) -> (usize, Vec<f64>) {
    let (cx, vals) = lower_star_cubical(t, f);
    let p = persistence(&cx, &vals);
    let mut pers: Vec<f64> = p.pairs.iter().filter(|q| q.0 == 0).map(|q| q.2 - q.1).collect();
    let global_min = p.essential.iter().find(|e| e.0 == 0).unwrap().1;
    let background = 0.0; // q is ~0 away from lumps after cooling
    pers.push(background - global_min);
    pers.sort_by(|a, b| b.partial_cmp(a).unwrap());
    (pers.iter().filter(|&&x| x > tau).count(), pers.into_iter().take(4).collect())
}
/// Basins of steepest ascent of g on the torus (2d neighbours); returns the integrated value of `w` per basin.
fn basin_charges(t: &CubicalTorus, g: &[f64], w: &[f64]) -> Vec<f64> {
    let d = t.dimension();
    let up = |x: usize| -> usize {
        let mut best = x;
        for dir in 0..d {
            let f = t.shift(x, dir);
            let mut b = x; // backward neighbour: shift sides-1 times
            for _ in 0..t.sides[dir] - 1 { b = t.shift(b, dir); }
            for y in [f, b] { if g[y] > g[best] { best = y; } }
        }
        best
    };
    let mut root = vec![usize::MAX; g.len()];
    let mut sums = std::collections::HashMap::new();
    for x in 0..g.len() {
        let mut y = x;
        let mut path = vec![];
        while root[y] == usize::MAX { let z = up(y); path.push(y); if z == y { root[y] = y; break; } y = z; }
        let r = root[y];
        for p in path { root[p] = r; }
        *sums.entry(r).or_insert(0.0) += w[x];
    }
    sums.into_values().collect()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let text = std::fs::read_to_string(&a[1]).unwrap();
    let tau: f64 = a.get(2).map(|x| x.parse().unwrap()).unwrap_or(0.002);
    let mut lines = text.lines();
    let sides: Vec<usize> = lines.next().unwrap().trim_start_matches("# sides ").split(',').map(|x| x.parse().unwrap()).collect();
    let t = CubicalTorus::new(&sides);
    let mut agree = 0;
    let mut basin_agree = 0;
    let mut n = 0;
    while let (Some(ql), Some(fl)) = (lines.next(), lines.next()) {
        let q_total: f64 = ql.split('\t').nth(1).unwrap().parse().unwrap();
        let q: Vec<f64> = fl.split('\t').map(|x| x.parse().unwrap()).collect();
        let neg: Vec<f64> = q.iter().map(|x| -x).collect();
        let (ni, top_i) = count(&t, &neg, tau);
        let (na, top_a) = count(&t, &q, tau);
        let qi = (q_total * 1.15).round() as i64;
        let pos: Vec<f64> = q.iter().map(|x| x.max(0.0)).collect();
        let negw: Vec<f64> = q.iter().map(|x| (-x).max(0.0)).collect();
        let bi = basin_charges(&t, &q, &pos).iter().filter(|&&c| c * 1.15 > 0.5).count() as i64;
        let ba = basin_charges(&t, &neg, &negw).iter().filter(|&&c| c * 1.15 > 0.5).count() as i64;
        if bi - ba == qi { basin_agree += 1; }
        n += 1;
        if ni as i64 - na as i64 == qi { agree += 1; }
        println!("Q = {q_total:+.3} (rounded {qi:+}): basins I/A = {bi}/{ba}; N_I = {ni}, N_A = {na}, N_I - N_A = {:+}; top instanton persistences {:?}; top anti {:?}",
            ni as i64 - na as i64, top_i.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>(), top_a.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>());
    }
    println!("tau = {tau}: N_I - N_A equals the rounded charge in {agree}/{n} configurations; basin count difference equals it in {basin_agree}/{n}");
}
