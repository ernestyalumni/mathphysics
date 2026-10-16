//! Exact port of the existing half-collinear (2,2) SMGA calculation.
//! This is NOT ordinary Lorentzian Yang--Mills Berends--Giele recursion.
//! Boundaries of sign chambers return errors, rather than a guessed Theta(0).
use crate::{scalar::Q, spinor::RealSpinor};
pub fn bracket(a: RealSpinor, b: RealSpinor) -> Q {
    a[1] * b[0] - a[0] * b[1]
}
fn sum(v: &[RealSpinor]) -> RealSpinor {
    v.iter()
        .fold([Q::ZERO; 2], |a, b| [a[0] + b[0], a[1] + b[1]])
}
fn validate(v: &[RealSpinor]) -> Result<(), String> {
    if !(3..=12).contains(&v.len()) {
        return Err("SMGA supports 3..12 legs".into());
    }
    if sum(v) != [Q::ZERO; 2] {
        return Err("half-collinear momentum not conserved".into());
    }
    Ok(())
}
fn nonzero_sign(q: Q) -> Result<i8, String> {
    if q.is_zero() {
        Err("on a sign-chamber wall".into())
    } else {
        Ok(q.sign())
    }
}
pub fn formula_r1(v: &[RealSpinor]) -> Result<Q, String> {
    validate(v)?;
    if v[0][0].sign() != -1 || v[1..].iter().any(|v| v[0].sign() != 1) {
        return Err("R1 requires omega0<0 and all other omegas>0".into());
    }
    let mut a = Q::ONE;
    for m in 1..v.len() - 1 {
        let s = nonzero_sign(bracket(v[m], v[m + 1]))?;
        let t = nonzero_sign(bracket(v[0], sum(&v[1..=m])))?;
        a = a * Q::new((s + t) as i128, 2);
    }
    Ok(a)
}
/// Pair (V,Vbar); evaluate every boundary condition even if a factor vanishes.
pub fn vertex(v: &[RealSpinor]) -> Result<(Q, Q), String> {
    let (mut direct, mut barred) = (Q::ONE, Q::ONE);
    for k in 1..v.len() {
        let adjacent = nonzero_sign(bracket(v[k - 1], v[k]))?;
        let split = nonzero_sign(bracket(sum(&v[..k]), sum(&v[k..])))?;
        let s = Q::from(adjacent as i64);
        if adjacent * split < 0 {
            direct = direct * s;
            barred = Q::ZERO;
        } else {
            barred = barred * s;
            direct = Q::ZERO;
        }
    }
    Ok((direct, barred))
}
/// Consecutive-block partitions of [start,end), represented by intervals.
pub fn partitions(start: usize, end: usize, min_parts: usize) -> Vec<Vec<(usize, usize)>> {
    assert!(end > start && end - start <= 12, "partition domain");
    let mut out = vec![];
    for mask in 0usize..(1usize << (end - start - 1)) {
        if (mask.count_ones() as usize) + 1 < min_parts {
            continue;
        }
        let mut blocks = vec![];
        let mut begin = start;
        for gap in start + 1..end {
            if mask & (1 << (gap - start - 1)) != 0 {
                blocks.push((begin, gap));
                begin = gap;
            }
        }
        blocks.push((begin, end));
        out.push(blocks);
    }
    out
}
/// Last leg is minus helicity. Dynamic programming over contiguous intervals
/// implements the preamplitudes without recursive Python calls or label maps.
pub fn recursion(minus_last: &[RealSpinor]) -> Result<Q, String> {
    validate(minus_last)?;
    let n = minus_last.len() - 1;
    let mut current = vec![vec![Q::ZERO; n + 1]; n];
    for len in 1..=n {
        for start in 0..=n - len {
            let end = start + len;
            if len == 1 {
                current[start][end] = Q::ONE;
                continue;
            }
            if len == 2 {
                continue;
            }
            let mut value = Q::ZERO;
            for blocks in partitions(start, end, 3) {
                let product = blocks.iter().fold(Q::ONE, |a, (s, e)| a * current[*s][*e]);
                if product.is_zero() {
                    continue;
                }
                let vectors: Vec<_> = blocks
                    .iter()
                    .map(|(s, e)| sum(&minus_last[*s..*e]))
                    .collect();
                value = value - vertex(&vectors)?.0 * product;
            }
            current[start][end] = value;
        }
    }
    let mut result = Q::ZERO;
    for blocks in partitions(0, n, 1) {
        let product = blocks.iter().fold(Q::ONE, |a, (s, e)| a * current[*s][*e]);
        if product.is_zero() {
            continue;
        }
        let vectors: Vec<_> = blocks
            .iter()
            .map(|(s, e)| sum(&minus_last[*s..*e]))
            .collect();
        let (v, vbar) = vertex(&vectors)?;
        result = result - (v - vbar) * product;
    }
    Ok(result)
}
pub fn minus_first_recursion(v: &[RealSpinor]) -> Result<Q, String> {
    validate(v)?;
    let mut rotated = v[1..].to_vec();
    rotated.push(v[0]);
    recursion(&rotated)
}
/// Deterministic exact test data, not a physical phase-space sampler.
/// seed=0 gives ordered positive slopes (nonzero examples); other seeds mix them.
pub fn fixture(n: usize, mut seed: u64) -> Vec<RealSpinor> {
    assert!((3..=12).contains(&n));
    let ordered = seed == 0;
    let mut v = vec![[Q::ZERO; 2]];
    for i in 1..n {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let w = if ordered {
            1
        } else {
            1 + ((seed >> 32) % 17) as i64
        };
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let wz = if ordered {
            (i * i) as i64
        } else {
            ((seed >> 32) % 1009) as i64 - 504
        };
        v.push([Q::from(w), Q::from(wz)]);
    }
    let total = sum(&v);
    v[0] = [-total[0], -total[1]];
    v
}
