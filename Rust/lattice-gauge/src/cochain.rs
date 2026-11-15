//! Sets of plaquettes modulo 2 (Z2 2-cochains): twist surfaces and their classes.
//!
//! A set of plaquettes is a `Vec<bool>` indexed by plaquette index. The coboundary `delta b`
//! of a link `b` is the set of plaquettes containing `b`. Twisted partition functions depend
//! on a set only through its class modulo sums of link coboundaries
//! (mono-notes, Lemma `ym:lem:twist-move`).
use crate::torus::Torus;

/// Symmetric difference of two sets of plaquettes. `None` if the lengths differ.
pub fn symmetric_difference(left: &[bool], right: &[bool]) -> Option<Vec<bool>> {
    if left.len() != right.len() {
        return None;
    }
    Some(left.iter().zip(right).map(|(&a, &b)| a != b).collect())
}

/// The set `delta b` for a link `b`.
pub fn link_coboundary_set(torus: &Torus, link: usize) -> Vec<bool> {
    let mut set = vec![false; torus.plaquette_count()];
    for plaquette in torus.link_coboundary(link) {
        set[plaquette] = !set[plaquette];
    }
    set
}

/// The 't Hooft twist surface: all plaquettes in the plane `(first, second)` at sites whose
/// coordinates in the directions `first` and `second` are both zero. Requires
/// `first < second < dimension`.
pub fn thooft_surface(torus: &Torus, first: usize, second: usize) -> Vec<bool> {
    let mut set = vec![false; torus.plaquette_count()];
    for site in 0..torus.site_count() {
        if torus.site_coordinate(site, first) == 0 && torus.site_coordinate(site, second) == 0 {
            set[torus.plaquette_index(site, first, second)] = true;
        }
    }
    set
}

/// Number of 3-cubes with an odd number of faces in the set (Z2 monopoles of the twist).
pub fn odd_cube_count(torus: &Torus, set: &[bool]) -> usize {
    torus
        .cubes()
        .iter()
        .filter(|faces| faces.iter().filter(|&&face| set[face]).count() % 2 == 1)
        .count()
}

/// Number of plaquettes of the set lying in the two-dimensional slice through `anchor_site`
/// spanned by the directions `first < second`.
pub fn plane_slice_intersection_count(
    torus: &Torus,
    set: &[bool],
    first: usize,
    second: usize,
    anchor_site: usize,
) -> usize {
    (0..torus.site_count())
        .filter(|&site| {
            (0..torus.dimension()).all(|direction| {
                direction == first
                    || direction == second
                    || torus.site_coordinate(site, direction)
                        == torus.site_coordinate(anchor_site, direction)
            })
        })
        .filter(|&site| set[torus.plaquette_index(site, first, second)])
        .count()
}

fn pack(set: &[bool]) -> Vec<u64> {
    let mut words = vec![0u64; set.len().div_ceil(64)];
    for (index, _) in set.iter().enumerate().filter(|(_, &member)| member) {
        words[index / 64] |= 1 << (index % 64);
    }
    words
}

fn lowest_set_bit(words: &[u64]) -> Option<usize> {
    words
        .iter()
        .enumerate()
        .find(|(_, &word)| word != 0)
        .map(|(index, &word)| index * 64 + word.trailing_zeros() as usize)
}

/// An echelon basis over GF(2) of the span of all link coboundaries, indexed by pivot bit.
fn link_coboundary_basis(torus: &Torus) -> Vec<Option<Vec<u64>>> {
    let mut basis: Vec<Option<Vec<u64>>> = vec![None; torus.plaquette_count()];
    for link in 0..torus.link_count() {
        let mut vector = pack(&link_coboundary_set(torus, link));
        while let Some(pivot) = lowest_set_bit(&vector) {
            match &basis[pivot] {
                Some(row) => vector.iter_mut().zip(row).for_each(|(word, &other)| *word ^= other),
                None => {
                    basis[pivot] = Some(vector);
                    break;
                }
            }
        }
    }
    basis
}

/// Rank over GF(2) of the coboundary map from links to plaquettes. For a torus it equals
/// `d V - V + 1 - d`: the closed Z2 1-cochains have dimension `(V - 1) + d`.
pub fn link_coboundary_rank(torus: &Torus) -> usize {
    link_coboundary_basis(torus).iter().flatten().count()
}

/// Whether the set is a symmetric difference of link coboundaries, decided by Gaussian
/// elimination over GF(2). `None` if the set has the wrong length.
pub fn is_sum_of_link_coboundaries(torus: &Torus, set: &[bool]) -> Option<bool> {
    if set.len() != torus.plaquette_count() {
        return None;
    }
    let basis = link_coboundary_basis(torus);
    let mut vector = pack(set);
    while let Some(pivot) = lowest_set_bit(&vector) {
        match &basis[pivot] {
            Some(row) => vector.iter_mut().zip(row).for_each(|(word, &other)| *word ^= other),
            None => return Some(false),
        }
    }
    Some(true)
}
