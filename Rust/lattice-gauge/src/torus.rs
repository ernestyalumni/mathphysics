//! Sites, links, plaquettes and cubes of the periodic torus `(Z/n_0) x ... x (Z/n_{d-1})`.
//!
//! Index conventions (stable; other modules and future GPU code rely on them):
//! - site index: `sum_i c_i * stride_i` with `stride_0 = 1`, `stride_{i+1} = stride_i * n_i`;
//! - link index: `site * d + direction`;
//! - plane index of `first < second`: lexicographic, `(0,1), (0,2), ..., (1,2), ...`;
//! - plaquette index: `site * plane_count + plane_index`.

/// A periodic torus with every side length at least 2.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Torus {
    side_lengths: Vec<usize>,
    strides: Vec<usize>,
    site_count: usize,
}

impl Torus {
    /// `None` if there are no sides, a side is shorter than 2 (a side of 1 makes a link a
    /// loop and breaks the plaquette conventions), or the number of sites overflows.
    pub fn new(side_lengths: &[usize]) -> Option<Self> {
        if side_lengths.is_empty() || side_lengths.iter().any(|&side| side < 2) {
            return None;
        }
        let mut strides = Vec::with_capacity(side_lengths.len());
        let mut site_count: usize = 1;
        for &side in side_lengths {
            strides.push(site_count);
            site_count = site_count.checked_mul(side)?;
        }
        site_count.checked_mul(side_lengths.len() * side_lengths.len())?;
        Some(Self {
            side_lengths: side_lengths.to_vec(),
            strides,
            site_count,
        })
    }
    pub fn dimension(&self) -> usize {
        self.side_lengths.len()
    }
    pub fn side_lengths(&self) -> &[usize] {
        &self.side_lengths
    }
    pub fn site_count(&self) -> usize {
        self.site_count
    }
    pub fn plane_count(&self) -> usize {
        self.dimension() * (self.dimension() - 1) / 2
    }
    pub fn link_count(&self) -> usize {
        self.site_count * self.dimension()
    }
    pub fn plaquette_count(&self) -> usize {
        self.site_count * self.plane_count()
    }
    /// Number of 3-cubes: `binomial(d, 3)` per site.
    pub fn cube_count(&self) -> usize {
        let d = self.dimension();
        if d < 3 {
            return 0;
        }
        self.site_count * (d * (d - 1) * (d - 2) / 6)
    }
    /// Coordinates are reduced modulo the side lengths.
    pub fn site_index(&self, coordinates: &[usize]) -> usize {
        coordinates
            .iter()
            .zip(&self.side_lengths)
            .zip(&self.strides)
            .map(|((&coordinate, &side), &stride)| (coordinate % side) * stride)
            .sum()
    }
    pub fn site_coordinate(&self, site: usize, direction: usize) -> usize {
        (site / self.strides[direction]) % self.side_lengths[direction]
    }
    pub fn site_coordinates(&self, site: usize) -> Vec<usize> {
        (0..self.dimension())
            .map(|direction| self.site_coordinate(site, direction))
            .collect()
    }
    /// The site one step along `direction`, forward or backward, with periodic wrap.
    pub fn neighbor(&self, site: usize, direction: usize, forward: bool) -> usize {
        let side = self.side_lengths[direction];
        let stride = self.strides[direction];
        let coordinate = self.site_coordinate(site, direction);
        if forward {
            if coordinate + 1 == side {
                site - (side - 1) * stride
            } else {
                site + stride
            }
        } else if coordinate == 0 {
            site + (side - 1) * stride
        } else {
            site - stride
        }
    }
    pub fn link_index(&self, site: usize, direction: usize) -> usize {
        site * self.dimension() + direction
    }
    /// Requires `first < second < dimension`.
    pub fn plane_index(&self, first: usize, second: usize) -> usize {
        debug_assert!(first < second && second < self.dimension());
        let d = self.dimension();
        first * (2 * d - first - 1) / 2 + (second - first - 1)
    }
    /// Inverse of [`Torus::plane_index`].
    pub fn plane_directions(&self, plane: usize) -> (usize, usize) {
        let d = self.dimension();
        let mut remaining = plane;
        for first in 0..d {
            let width = d - first - 1;
            if remaining < width {
                return (first, first + 1 + remaining);
            }
            remaining -= width;
        }
        panic!("plane index out of range");
    }
    pub fn plaquette_index(&self, site: usize, first: usize, second: usize) -> usize {
        site * self.plane_count() + self.plane_index(first, second)
    }
    /// The four links of the plaquette, in the order of the plaquette variable
    /// `U(y,mu) U(y+e_mu,nu) U(y+e_nu,mu)^{-1} U(y,nu)^{-1}` with `mu = first`, `nu = second`.
    pub fn plaquette_links(&self, site: usize, first: usize, second: usize) -> [usize; 4] {
        [
            self.link_index(site, first),
            self.link_index(self.neighbor(site, first, true), second),
            self.link_index(self.neighbor(site, second, true), first),
            self.link_index(site, second),
        ]
    }
    /// The plaquettes that contain the link: `2 (d - 1)` of them.
    pub fn link_coboundary(&self, link: usize) -> Vec<usize> {
        let d = self.dimension();
        let (site, direction) = (link / d, link % d);
        let mut plaquettes = Vec::with_capacity(2 * (d - 1));
        for other in (0..d).filter(|&other| other != direction) {
            let (first, second) = (direction.min(other), direction.max(other));
            plaquettes.push(self.plaquette_index(site, first, second));
            plaquettes.push(self.plaquette_index(self.neighbor(site, other, false), first, second));
        }
        plaquettes
    }
    /// The six faces of the 3-cube at `site` spanned by `first < second < third`.
    pub fn cube_plaquettes(&self, site: usize, first: usize, second: usize, third: usize) -> [usize; 6] {
        [
            self.plaquette_index(site, first, second),
            self.plaquette_index(self.neighbor(site, third, true), first, second),
            self.plaquette_index(site, first, third),
            self.plaquette_index(self.neighbor(site, second, true), first, third),
            self.plaquette_index(site, second, third),
            self.plaquette_index(self.neighbor(site, first, true), second, third),
        ]
    }
    /// All 3-cubes, each as its six faces.
    pub fn cubes(&self) -> Vec<[usize; 6]> {
        let d = self.dimension();
        let mut cubes = Vec::with_capacity(self.cube_count());
        for site in 0..self.site_count {
            for first in 0..d {
                for second in (first + 1)..d {
                    for third in (second + 1)..d {
                        cubes.push(self.cube_plaquettes(site, first, second, third));
                    }
                }
            }
        }
        cubes
    }
}
