#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <type_traits>
#include <vector>
#include "ducc0/fft/fft.h"
#include "ducc0/fft/fftnd_impl.h"
namespace fftconv {

// M_E is POSIX, not standard C++ (absent on MSVC-style compilers without
// _USE_MATH_DEFINES), so use our own constant for Euler's number.
constexpr double euler_e = 2.718281828459045235360287471352662498;

// Real lower branch W_{-1}(x) for x in [-1/e, 0), via Halley iteration.
// Matches scipy.special.lambertw(x, k=-1).real, used by _calc_oa_lens.
inline double lambertw_wm1(double x) {
    // W_{-1}(-1/e) = -1 exactly; both real branches meet here. Guard the
    // boundary where the Halley denominator degenerates to 0/0 (w -> -1).
    if (x <= -1.0 / euler_e) return -1.0;
    // Initial guess from the asymptotic expansion near 0^- (Corless et al.).
    const double L1 = std::log(-x);
    const double L2 = std::log(-L1);
    double w = L1 - L2 + L2 / L1;
    for (int i = 0; i < 20; ++i) {
        const double ew = std::exp(w);
        const double f = w * ew - x;
        // Halley step
        const double denom = ew * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0);
        const double dw = f / denom;
        w -= dw;
        if (std::fabs(dw) <= 1e-16 * (1.0 + std::fabs(w))) break;
    }
    return w;
}

struct OaLens { size_t block_size, overlap, in1_step, in2_step; bool fallback; };

// Port of _calc_oa_lens(s1, s2) from _signaltools.py.
inline OaLens calc_oa_lens(size_t s1, size_t s2) {
    OaLens f{s1 + s2 - 1, 0, s1, s2, true};
    if (s1 == s2 || s1 == 1 || s2 == 1) return f;
    bool swapped = false;
    size_t a = s1, b = s2;
    if (b > a) { std::swap(a, b); swapped = true; }
    if (b >= a / 2.0) return f;                       // s2 >= s1/2
    const double overlap = double(b) - 1.0;
    const double opt = -overlap *
        lambertw_wm1(-1.0 / (2.0 * euler_e * overlap));
    size_t block_size = ducc0::good_size_complex(size_t(std::ceil(opt)));
    if (block_size >= a) return f;                    // only one block
    OaLens r;
    r.fallback = false;
    r.block_size = block_size;
    r.overlap = size_t(overlap);
    if (!swapped) { r.in1_step = block_size - b + 1; r.in2_step = b; }
    else          { r.in1_step = b; r.in2_step = block_size - b + 1; }
    return r;
}

template<typename T> struct is_complex : std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : std::true_type {};
template<typename T> inline constexpr bool is_complex_v = is_complex<T>::value;

template<typename T> struct real_type { using type = T; };
template<typename T> struct real_type<std::complex<T>> { using type = T; };
template<typename T> using real_t = typename real_type<T>::type;

enum class ConvMode { Full, Same, Valid };

// N-D single-axis overlap-add convolution. in1 is (outer, n, inner), in2 is
// (outer, m, inner), out is (outer, out_n, inner), all C-contiguous. Convolves
// the middle axis, batching over outer*inner. Real via r2c/c2r, complex via c2c
// (if constexpr). Returns false (decline) on the calc_oa_lens fallback.
template<typename T>
bool oaconvolve_axis(const T* in1, size_t outer, size_t n, size_t inner,
                     const T* in2, size_t m, T* out, ConvMode mode, size_t nthreads) {
    OaLens L = calc_oa_lens(n, m);
    if (L.fallback) return false;

    // Orient so `sig` has the longer conv length. Keep original n,m for slicing.
    const T *sig = in1, *ker = in2;
    size_t ns = n, nk = m;
    if (nk > ns) { std::swap(sig, ker); std::swap(ns, nk); }

    const size_t B  = L.block_size;
    const size_t S  = B - (nk - 1);
    using R  = real_t<T>;
    using Cx = std::complex<R>;
    constexpr bool cplx = is_complex_v<T>;
    const size_t Bc   = cplx ? B : (B / 2 + 1);
    const size_t nblk = (ns + S - 1) / S;
    const R fct       = R(1) / R(B);
    const size_t full = ns + nk - 1;

    // Kernel spectrum: pad ker (outer, nk, inner) -> (outer, B, inner), FFT over axis 1.
    std::vector<T> kbuf(outer * B * inner, T(0));
    for (size_t o = 0; o < outer; ++o)
        std::copy(ker + o * nk * inner, ker + (o + 1) * nk * inner,
                  kbuf.begin() + o * B * inner);
    std::vector<Cx> kspec(outer * Bc * inner);
    {
        ducc0::cfmav<T>  kin(kbuf.data(), ducc0::fmav_info::shape_t{outer, B, inner});
        ducc0::vfmav<Cx> kout(kspec.data(), ducc0::fmav_info::shape_t{outer, Bc, inner});
        if constexpr (cplx) ducc0::c2c(kin, kout, ducc0::fmav_info::shape_t{1}, true, R(1), nthreads);
        else                ducc0::r2c(kin, kout, /*axis=*/1, true, R(1), nthreads);
    }

    // Blocked signal (outer, nblk, B, inner), zero-padded.
    std::vector<T> blocks(outer * nblk * B * inner, T(0));
    for (size_t o = 0; o < outer; ++o)
        for (size_t blk = 0; blk < nblk; ++blk) {
            size_t start = blk * S;
            size_t len = std::min(S, ns > start ? ns - start : 0);
            if (len)
                std::copy(sig + (o * ns + start) * inner,
                          sig + (o * ns + start + len) * inner,
                          blocks.begin() + ((o * nblk + blk) * B) * inner);
        }

    // Batched forward FFT over the B axis (axis 2 of the 4-D layout).
    std::vector<Cx> spec(outer * nblk * Bc * inner);
    {
        ducc0::cfmav<T>  bin(blocks.data(), ducc0::fmav_info::shape_t{outer, nblk, B, inner});
        ducc0::vfmav<Cx> bout(spec.data(), ducc0::fmav_info::shape_t{outer, nblk, Bc, inner});
        if constexpr (cplx) ducc0::c2c(bin, bout, ducc0::fmav_info::shape_t{2}, true, R(1), nthreads);
        else                ducc0::r2c(bin, bout, /*axis=*/2, true, R(1), nthreads);
    }
    // Multiply by kernel spectrum, broadcast over the nblk axis.
    for (size_t o = 0; o < outer; ++o)
        for (size_t blk = 0; blk < nblk; ++blk) {
            Cx* sp = spec.data() + ((o * nblk + blk) * Bc) * inner;
            const Cx* kp = kspec.data() + (o * Bc) * inner;
            for (size_t t = 0; t < Bc * inner; ++t) sp[t] *= kp[t];
        }
    // Batched inverse FFT over the B axis, scaled 1/B.
    {
        ducc0::cfmav<Cx> sin(spec.data(), ducc0::fmav_info::shape_t{outer, nblk, Bc, inner});
        ducc0::vfmav<T>  sout(blocks.data(), ducc0::fmav_info::shape_t{outer, nblk, B, inner});
        if constexpr (cplx) ducc0::c2c(sin, sout, ducc0::fmav_info::shape_t{2}, false, fct, nthreads);
        else                ducc0::c2r(sin, sout, /*axis=*/2, false, fct, nthreads);
    }

    // Overlap-add into acc (outer, full, inner).
    std::vector<T> acc(outer * full * inner, T(0));
    for (size_t o = 0; o < outer; ++o)
        for (size_t blk = 0; blk < nblk; ++blk) {
            size_t start = blk * S;
            size_t cnt = (start < full ? std::min(B, full - start) : 0) * inner;
            T* dst = acc.data() + (o * full + start) * inner;
            const T* src = blocks.data() + ((o * nblk + blk) * B) * inner;
            for (size_t t = 0; t < cnt; ++t) dst[t] += src[t];
        }

    // Mode slice along the middle axis into out. Offsets use the ORIGINAL n,m
    // (so the internal sig/ker swap is invisible), identical to the 1-D core.
    size_t offset, out_n;
    if (mode == ConvMode::Full)      { offset = 0;           out_n = n + m - 1; }
    else if (mode == ConvMode::Same) { offset = (m - 1) / 2; out_n = n; }
    else { size_t big = std::max(n, m), sm = std::min(n, m); offset = sm - 1; out_n = big - sm + 1; }
    for (size_t o = 0; o < outer; ++o)
        std::copy(acc.begin() + (o * full + offset) * inner,
                  acc.begin() + (o * full + offset + out_n) * inner,
                  out + o * out_n * inner);
    return true;
}

// Full 1-D overlap-add convolution. out must have length s1+s2-1 for Full.
// Returns false (decline) when calc_oa_lens says fallback; caller handles it.
template<typename T>
bool oaconvolve_1d(const T* in1, size_t s1, const T* in2, size_t s2,
                   T* out, ConvMode mode, size_t nthreads) {
    OaLens L = calc_oa_lens(s1, s2);
    if (L.fallback) return false;  // decline: caller routes to fftconvolve

    // Orient so `sig` is the longer, `ker` the shorter (OA is symmetric).
    const T *sig = in1, *ker = in2;
    size_t ns = s1, nk = s2;
    if (nk > ns) { std::swap(sig, ker); std::swap(ns, nk); }

    const size_t B  = L.block_size;              // FFT length
    const size_t S  = B - (nk - 1);              // step (valid samples/block)
    using R  = real_t<T>;              // float or double
    using Cx = std::complex<R>;        // spectrum element type
    constexpr bool cplx = is_complex_v<T>;
    const size_t Bc = cplx ? B : (B / 2 + 1);   // c2c: full length; r2c: half+1
    const size_t nblk = (ns + S - 1) / S;
    const R fct = R(1) / R(B);

    // Kernel spectrum (once).
    std::vector<T> kbuf(B, T(0));
    std::copy(ker, ker + nk, kbuf.begin());
    std::vector<Cx> kspec(Bc);
    {
        ducc0::cfmav<T>  kin(kbuf.data(), ducc0::fmav_info::shape_t{B});
        ducc0::vfmav<Cx> kout(kspec.data(), ducc0::fmav_info::shape_t{Bc});
        if constexpr (cplx) ducc0::c2c(kin, kout, ducc0::fmav_info::shape_t{0}, /*forward=*/true, R(1), nthreads);
        else                ducc0::r2c(kin, kout, /*axis=*/0, /*forward=*/true, R(1), nthreads);
    }

    // Blocked signal buffer (nblk x B), zero-padded.
    std::vector<T> blocks(nblk * B, T(0));
    for (size_t i = 0; i < nblk; ++i) {
        size_t start = i * S;
        size_t len = std::min(S, ns > start ? ns - start : 0);
        if (len) std::copy(sig + start, sig + start + len, blocks.begin() + i * B);
    }

    // Batched forward transform over the block-length axis (axis 1).
    std::vector<Cx> spec(nblk * Bc);
    {
        ducc0::cfmav<T>  bin(blocks.data(), ducc0::fmav_info::shape_t{nblk, B});
        ducc0::vfmav<Cx> bout(spec.data(), ducc0::fmav_info::shape_t{nblk, Bc});
        if constexpr (cplx) ducc0::c2c(bin, bout, ducc0::fmav_info::shape_t{1}, true, R(1), nthreads);
        else                ducc0::r2c(bin, bout, /*axis=*/1, true, R(1), nthreads);
    }
    // Multiply each block spectrum by the kernel spectrum.
    // For the complex path this is a full-length (Bc == B) elementwise product;
    // the real path uses the half spectrum (Bc == B/2+1).
    for (size_t i = 0; i < nblk; ++i)
        for (size_t j = 0; j < Bc; ++j)
            spec[i * Bc + j] *= kspec[j];
    // Batched inverse transform, scaled by 1/B.
    {
        ducc0::cfmav<Cx> sin(spec.data(), ducc0::fmav_info::shape_t{nblk, Bc});
        ducc0::vfmav<T>  sout(blocks.data(), ducc0::fmav_info::shape_t{nblk, B});
        if constexpr (cplx) ducc0::c2c(sin, sout, ducc0::fmav_info::shape_t{1}, false, fct, nthreads);
        else                ducc0::c2r(sin, sout, /*axis=*/1, false, fct, nthreads);
    }

    // Overlap-add into a full-length accumulator.
    const size_t full = ns + nk - 1;
    std::vector<T> acc(full, T(0));
    for (size_t i = 0; i < nblk; ++i) {
        size_t start = i * S;
        for (size_t j = 0; j < B && start + j < full; ++j)
            acc[start + j] += blocks[i * B + j];
    }

    // Mode slicing. `acc` holds the length-(ns+nk-1) full convolution.
    // Note: for Same, the output length is s1 (the ORIGINAL first input),
    // centered in the full result. Because we may have swapped sig/ker,
    // recompute offsets from the original s1/s2.
    if (mode == ConvMode::Full) {
        std::copy(acc.begin(), acc.end(), out);
    } else if (mode == ConvMode::Same) {
        size_t offset = (s2 - 1) / 2;          // center wrt in1
        std::copy(acc.begin() + offset, acc.begin() + offset + s1, out);
    } else {  // Valid: length max-min+1, starts at min-1
        size_t big = std::max(s1, s2), small = std::min(s1, s2);
        size_t nvalid = big - small + 1;
        std::copy(acc.begin() + (small - 1),
                  acc.begin() + (small - 1) + nvalid, out);
    }
    return true;
}

}  // namespace fftconv
