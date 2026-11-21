// CUDA C++ path for SU(2) Wilson updates (roadmap T13). C ABI only; Rust calls it through src/cuda.rs.
//
// Mirrors src/philox.rs, src/su2.rs and Lattice::staple/sweep_with in src/wilson.rs operation for
// operation: the same Philox4x32-10 stream keyed by (seed, sweep, link, replica), the same staple
// summation order, the same quaternion product formula. Compile with -fmad=false so that products
// round as on the CPU. Device log/cos may differ from glibc in the last bit, which can flip a rare
// accept decision; tests therefore count differing links rather than require bitwise identity.

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdlib.h>

struct Q { double a0, a1, a2, a3; };

__device__ __forceinline__ Q qmul(Q x, Q y) {
    Q r;
    r.a0 = x.a0 * y.a0 - x.a1 * y.a1 - x.a2 * y.a2 - x.a3 * y.a3;
    r.a1 = x.a0 * y.a1 + y.a0 * x.a1 - (x.a2 * y.a3 - x.a3 * y.a2);
    r.a2 = x.a0 * y.a2 + y.a0 * x.a2 - (x.a3 * y.a1 - x.a1 * y.a3);
    r.a3 = x.a0 * y.a3 + y.a0 * x.a3 - (x.a1 * y.a2 - x.a2 * y.a1);
    return r;
}
__device__ __forceinline__ Q qconj(Q x) { return Q{x.a0, -x.a1, -x.a2, -x.a3}; }
__device__ __forceinline__ Q qadd(Q x, Q y) { return Q{x.a0 + y.a0, x.a1 + y.a1, x.a2 + y.a2, x.a3 + y.a3}; }
__device__ __forceinline__ Q qscale(Q x, double s) { return Q{x.a0 * s, x.a1 * s, x.a2 * s, x.a3 * s}; }
__device__ __forceinline__ double qnorm_sq(Q x) {
    // same left-to-right order as the Rust iterator sum
    double s = 0.0;
    s += x.a0 * x.a0; s += x.a1 * x.a1; s += x.a2 * x.a2; s += x.a3 * x.a3;
    return s;
}

// ---- Philox4x32-10 -------------------------------------------------------------------------------
struct Stream { uint32_t key[2]; uint32_t ctr[4]; uint32_t buf[4]; int used; };

__device__ void philox(const uint32_t c_in[4], const uint32_t k_in[2], uint32_t out[4]) {
    uint32_t c0 = c_in[0], c1 = c_in[1], c2 = c_in[2], c3 = c_in[3];
    uint32_t k0 = k_in[0], k1 = k_in[1];
    for (int round = 0; round < 10; ++round) {
        if (round > 0) { k0 += 0x9E3779B9u; k1 += 0xBB67AE85u; }
        uint32_t hi0 = __umulhi(0xD2511F53u, c0), lo0 = 0xD2511F53u * c0;
        uint32_t hi1 = __umulhi(0xCD9E8D57u, c2), lo1 = 0xCD9E8D57u * c2;
        uint32_t n0 = hi1 ^ c1 ^ k0, n2 = hi0 ^ c3 ^ k1;
        c0 = n0; c1 = lo1; c2 = n2; c3 = lo0;
    }
    out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
}

__device__ void stream_init(Stream &s, uint64_t seed, uint32_t sweep, uint32_t link, uint32_t replica) {
    s.key[0] = (uint32_t)seed; s.key[1] = (uint32_t)(seed >> 32);
    s.ctr[0] = 0; s.ctr[1] = sweep; s.ctr[2] = link; s.ctr[3] = replica;
    s.used = 4;
}
__device__ uint32_t next_u32(Stream &s) {
    if (s.used == 4) { philox(s.ctr, s.key, s.buf); s.ctr[0] += 1; s.used = 0; }
    return s.buf[s.used++];
}
__device__ double next_open01(Stream &s) {
    uint64_t hi = next_u32(s), lo = next_u32(s);
    uint64_t bits = ((hi << 32) | lo) >> 11;
    return ((double)bits + 0.5) * (1.0 / 9007199254740992.0);
}
__device__ double next_gaussian(Stream &s) {
    double r1 = next_open01(s), r2 = next_open01(s);
    return sqrt(-2.0 * log(r1)) * cos(2.0 * 3.141592653589793 * r2);
}

// ---- SU(2) sampling ------------------------------------------------------------------------------
__device__ Q haar(Stream &s) {
    for (;;) {
        Q q{next_gaussian(s), next_gaussian(s), next_gaussian(s), next_gaussian(s)};
        double n = qnorm_sq(q);
        if (n > 1e-300) return qscale(q, 1.0 / sqrt(n));
    }
}

__device__ Q heat_bath(Q staple, double beta, Stream &s) {
    double k = sqrt(qnorm_sq(staple));
    double a = beta * k;
    if (a < 1e-12) return haar(s);
    Q v = qscale(staple, 1.0 / k);
    double x0;
    for (;;) {
        double r1 = next_open01(s), r2 = next_open01(s), r3 = next_open01(s);
        double c = cos(2.0 * 3.141592653589793 * r2);
        double lambda2 = -(log(r1) + c * c * log(r3)) / (2.0 * a);
        double r4 = next_open01(s);
        if (r4 * r4 <= 1.0 - lambda2) { x0 = fmax(1.0 - 2.0 * lambda2, -1.0); break; }
    }
    double r = sqrt(fmax(1.0 - x0 * x0, 0.0));
    double cos_t = 2.0 * next_open01(s) - 1.0;
    double sin_t = sqrt(fmax(1.0 - cos_t * cos_t, 0.0));
    double phi = 2.0 * 3.141592653589793 * next_open01(s);
    Q x{x0, r * sin_t * cos(phi), r * sin_t * sin(phi), r * cos_t};
    return qmul(x, qconj(v));
}

__device__ Q overrelax(Q u, Q staple) {
    double n = qnorm_sq(staple);
    if (n < 1e-300) return u;
    Q vd = qconj(qscale(staple, 1.0 / sqrt(n)));
    return qmul(qmul(vd, qconj(u)), vd);
}

// ---- lattice ---------------------------------------------------------------------------------------
// nb[(site * d + dir) * 2 + forward] = neighbour site; links[site * d + dir].
__global__ void update_class(Q *links, const int *nb, const int *class_links, int count, int d, int kind,
                             double beta, uint64_t seed, uint32_t sweep, uint32_t replica) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= count) return;
    int link = class_links[i];
    int site = link / d, mu = link % d;
    int x_mu = nb[(site * d + mu) * 2 + 1];
    Q w{0.0, 0.0, 0.0, 0.0};
    for (int nu = 0; nu < d; ++nu) {
        if (nu == mu) continue;
        int x_nu = nb[(site * d + nu) * 2 + 1];
        int x_mnu = nb[(site * d + nu) * 2 + 0];
        int x_mu_mnu = nb[(x_mu * d + nu) * 2 + 0];
        Q up = qmul(qmul(links[x_mu * d + nu], qconj(links[x_nu * d + mu])), qconj(links[site * d + nu]));
        Q down = qmul(qmul(qconj(links[x_mu_mnu * d + nu]), qconj(links[x_mnu * d + mu])), links[x_mnu * d + nu]);
        w = qadd(qadd(w, up), down);
    }
    if (kind == 0) {
        Stream s; stream_init(s, seed, sweep, (uint32_t)link, replica);
        links[link] = heat_bath(w, beta, s);
    } else {
        links[link] = overrelax(links[link], w);
    }
}

// Sum of a0 over all plaquettes; one thread per site, double atomics (sm_60+). Summation order varies,
// so the result can differ from the CPU in the last bits; it never feeds back into the Markov chain.
__global__ void plaquette_sum(const Q *links, const int *nb, int n_sites, int d, double *out) {
    int site = blockIdx.x * blockDim.x + threadIdx.x;
    if (site >= n_sites) return;
    double s = 0.0;
    for (int mu = 0; mu < d; ++mu) {
        int x_mu = nb[(site * d + mu) * 2 + 1];
        for (int nu = mu + 1; nu < d; ++nu) {
            int x_nu = nb[(site * d + nu) * 2 + 1];
            Q p = qmul(qmul(qmul(links[site * d + mu], links[x_mu * d + nu]), qconj(links[x_nu * d + mu])), qconj(links[site * d + nu]));
            s += p.a0;
        }
    }
    atomicAdd(out, s);
}

struct Handle {
    int n_links, d, n_classes;
    Q *links; int *nb; int *class_links;
    int *offsets; // host copy, n_classes + 1
};

extern "C" {

// Returns 0 on success, a CUDA error code otherwise.
int lg_cuda_init(int n_sites, int d, const int *nb, const int *class_links, const int *offsets, int n_classes,
                 void **out) {
    Handle *h = (Handle *)calloc(1, sizeof(Handle));
    h->n_links = n_sites * d; h->d = d; h->n_classes = n_classes;
    h->offsets = (int *)malloc(sizeof(int) * (n_classes + 1));
    for (int c = 0; c <= n_classes; ++c) h->offsets[c] = offsets[c];
    cudaError_t e;
    if ((e = cudaMalloc(&h->links, sizeof(Q) * h->n_links))) return (int)e;
    if ((e = cudaMalloc(&h->nb, sizeof(int) * n_sites * d * 2))) return (int)e;
    if ((e = cudaMalloc(&h->class_links, sizeof(int) * offsets[n_classes]))) return (int)e;
    if ((e = cudaMemcpy(h->nb, nb, sizeof(int) * n_sites * d * 2, cudaMemcpyHostToDevice))) return (int)e;
    if ((e = cudaMemcpy(h->class_links, class_links, sizeof(int) * offsets[n_classes], cudaMemcpyHostToDevice))) return (int)e;
    *out = h;
    return 0;
}

int lg_cuda_upload(void *hp, const double *links) {
    Handle *h = (Handle *)hp;
    return (int)cudaMemcpy(h->links, links, sizeof(Q) * h->n_links, cudaMemcpyHostToDevice);
}

int lg_cuda_download(void *hp, double *links) {
    Handle *h = (Handle *)hp;
    return (int)cudaMemcpy(links, h->links, sizeof(Q) * h->n_links, cudaMemcpyDeviceToHost);
}

// kind 0 = heat bath, 1 = over-relaxation. Classes are processed in order, as on the CPU.
int lg_cuda_sweep(void *hp, int kind, double beta, uint64_t seed, uint32_t sweep, uint32_t replica) {
    Handle *h = (Handle *)hp;
    for (int c = 0; c < h->n_classes; ++c) {
        int count = h->offsets[c + 1] - h->offsets[c];
        int threads = 128, blocks = (count + threads - 1) / threads;
        update_class<<<blocks, threads>>>(h->links, h->nb, h->class_links + h->offsets[c], count, h->d, kind, beta,
                                          seed, sweep, replica);
    }
    return (int)cudaDeviceSynchronize();
}

// Average of a0 over all plaquettes.
int lg_cuda_plaquette(void *hp, double *result) {
    Handle *h = (Handle *)hp;
    int n_sites = h->n_links / h->d;
    double *dev; cudaError_t e;
    if ((e = cudaMalloc(&dev, sizeof(double)))) return (int)e;
    cudaMemset(dev, 0, sizeof(double));
    plaquette_sum<<<(n_sites + 127) / 128, 128>>>(h->links, h->nb, n_sites, h->d, dev);
    double sum = 0.0;
    e = cudaMemcpy(&sum, dev, sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(dev);
    *result = sum / ((double)n_sites * h->d * (h->d - 1) / 2);
    return (int)e;
}

void lg_cuda_free(void *hp) {
    Handle *h = (Handle *)hp;
    if (!h) return;
    cudaFree(h->links); cudaFree(h->nb); cudaFree(h->class_links);
    free(h->offsets); free(h);
}

}  // extern "C"
