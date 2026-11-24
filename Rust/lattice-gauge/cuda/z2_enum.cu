// Weight enumerator of a coset offset + span(basis) in F2^n, n <= 128, by Gray-code walks on the GPU.
// Codeword(i) = offset + sum of basis[k] over the set bits k of gray(i) = i ^ (i >> 1); consecutive indices differ
// by basis[ctz(i)]. Each thread walks a contiguous range with a private histogram in local memory (no contention),
// then adds it to the global histogram. Launches cover at most 2^34 indices each, so no kernel runs long.
#include <cstdint>
#include <cuda_runtime.h>

__constant__ uint64_t c_lo[64];
__constant__ uint64_t c_hi[64];

__global__ void z2_walk(int r, uint64_t off_lo, uint64_t off_hi, uint64_t start, uint64_t end, uint64_t per_thread,
                        int n, unsigned long long *hist) {
    unsigned int cnt[129];
    for (int w = 0; w <= n; ++w) cnt[w] = 0;
    uint64_t tid = (uint64_t)blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t lo = start + tid * per_thread;
    if (lo < end) {
        uint64_t hi = lo + per_thread < end ? lo + per_thread : end;
        uint64_t g = lo ^ (lo >> 1);
        uint64_t vlo = off_lo, vhi = off_hi;
        for (int k = 0; k < r; ++k)
            if ((g >> k) & 1ull) { vlo ^= c_lo[k]; vhi ^= c_hi[k]; }
        cnt[__popcll(vlo) + __popcll(vhi)] += 1u;
        for (uint64_t i = lo + 1; i < hi; ++i) {
            int k = __ffsll((long long)i) - 1;
            vlo ^= c_lo[k];
            vhi ^= c_hi[k];
            cnt[__popcll(vlo) + __popcll(vhi)] += 1u;
        }
    }
    for (int w = 0; w <= n; ++w)
        if (cnt[w]) atomicAdd(&hist[w], (unsigned long long)cnt[w]);
}

extern "C" {

// Returns 0 on success, a CUDA error code otherwise. out has n + 1 entries.
int z2_enum_coset(const uint64_t *lo, const uint64_t *hi, int r, uint64_t off_lo, uint64_t off_hi, int n, uint64_t *out) {
    if (r < 0 || r > 62 || n < 0 || n > 128) return -1;
    cudaError_t e;
    if ((e = cudaMemcpyToSymbol(c_lo, lo, sizeof(uint64_t) * (r > 0 ? r : 1))) != cudaSuccess) return (int)e;
    if ((e = cudaMemcpyToSymbol(c_hi, hi, sizeof(uint64_t) * (r > 0 ? r : 1))) != cudaSuccess) return (int)e;
    unsigned long long *d_hist = nullptr;
    if ((e = cudaMalloc(&d_hist, sizeof(unsigned long long) * (n + 1))) != cudaSuccess) return (int)e;
    cudaMemset(d_hist, 0, sizeof(unsigned long long) * (n + 1));
    const uint64_t total = 1ull << r;
    const uint64_t per_thread = 1ull << 14;
    const uint64_t batch = 1ull << 34;
    const int threads = 256;
    for (uint64_t start = 0; start < total; start += batch) {
        uint64_t end = start + batch < total ? start + batch : total;
        uint64_t work = (end - start + per_thread - 1) / per_thread;
        unsigned int blocks = (unsigned int)((work + threads - 1) / threads);
        z2_walk<<<blocks, threads>>>(r, off_lo, off_hi, start, end, per_thread, n, d_hist);
        if ((e = cudaDeviceSynchronize()) != cudaSuccess) { cudaFree(d_hist); return (int)e; }
    }
    e = cudaMemcpy(out, d_hist, sizeof(unsigned long long) * (n + 1), cudaMemcpyDeviceToHost);
    cudaFree(d_hist);
    return e == cudaSuccess ? 0 : (int)e;
}

}  // extern "C"
