//! GPU path for SU(2) Wilson updates (feature `cuda`): a thin wrapper over the C ABI of
//! cuda/su2_update.cu. The GPU lattice mirrors `wilson::Lattice` (same link layout, update classes,
//! Philox streams and sweep counter), so a run can be compared with the CPU path link by link.

use crate::su2::Quaternion;
use crate::wilson::{Lattice, Update};
use std::os::raw::{c_int, c_void};

extern "C" {
    fn lg_cuda_init(n_sites: c_int, d: c_int, nb: *const c_int, class_links: *const c_int, offsets: *const c_int, n_classes: c_int, out: *mut *mut c_void) -> c_int;
    fn lg_cuda_upload(h: *mut c_void, links: *const f64) -> c_int;
    fn lg_cuda_download(h: *mut c_void, links: *mut f64) -> c_int;
    fn lg_cuda_sweep(h: *mut c_void, kind: c_int, beta: f64, seed: u64, sweep: u32, replica: u32) -> c_int;
    fn lg_cuda_free(h: *mut c_void);
}

pub struct GpuLattice {
    handle: *mut c_void,
    n_links: usize,
    pub beta: f64,
    pub seed: u64,
    pub replica: u32,
    pub sweep: u32,
}

fn check(code: c_int, what: &str) -> Result<(), String> {
    if code == 0 { Ok(()) } else { Err(format!("{what}: CUDA error {code}")) }
}

impl GpuLattice {
    /// Copy geometry, update classes, links and counters from a CPU lattice.
    pub fn from_cpu(l: &Lattice) -> Result<Self, String> {
        let t = &l.torus;
        let d = t.dimension();
        let n = t.site_count();
        let mut nb = Vec::with_capacity(n * d * 2);
        for site in 0..n {
            for dir in 0..d {
                nb.push(t.neighbor(site, dir, false) as c_int);
                nb.push(t.neighbor(site, dir, true) as c_int);
            }
        }
        let classes = l.update_classes();
        let mut offsets = vec![0 as c_int];
        let mut flat = Vec::new();
        for c in classes {
            flat.extend(c.iter().map(|&x| x as c_int));
            offsets.push(flat.len() as c_int);
        }
        let mut handle = std::ptr::null_mut();
        check(unsafe { lg_cuda_init(n as c_int, d as c_int, nb.as_ptr(), flat.as_ptr(), offsets.as_ptr(), classes.len() as c_int, &mut handle) }, "init")?;
        let g = GpuLattice { handle, n_links: l.links.len(), beta: l.beta, seed: l.seed, replica: l.replica, sweep: l.sweep };
        g.upload(&l.links)?;
        Ok(g)
    }

    pub fn upload(&self, links: &[Quaternion]) -> Result<(), String> {
        assert_eq!(links.len(), self.n_links);
        check(unsafe { lg_cuda_upload(self.handle, links.as_ptr() as *const f64) }, "upload")
    }

    pub fn download(&self) -> Result<Vec<Quaternion>, String> {
        let mut links = vec![Quaternion::ZERO; self.n_links];
        check(unsafe { lg_cuda_download(self.handle, links.as_mut_ptr() as *mut f64) }, "download")?;
        Ok(links)
    }

    /// One sweep; increments the sweep counter, as `Lattice::sweep_with`.
    pub fn sweep_with(&mut self, update: Update) -> Result<(), String> {
        let kind = match update {
            Update::HeatBath => 0,
            Update::OverRelax => 1,
            Update::HeatBathBroken => return Err("the broken update exists only on the CPU".into()),
        };
        check(unsafe { lg_cuda_sweep(self.handle, kind, self.beta, self.seed, self.sweep, self.replica) }, "sweep")?;
        self.sweep = self.sweep.wrapping_add(1);
        Ok(())
    }

    pub fn update(&mut self, overrelax: usize) -> Result<(), String> {
        self.sweep_with(Update::HeatBath)?;
        for _ in 0..overrelax {
            self.sweep_with(Update::OverRelax)?;
        }
        Ok(())
    }
}

impl Drop for GpuLattice {
    fn drop(&mut self) {
        unsafe { lg_cuda_free(self.handle) }
    }
}
