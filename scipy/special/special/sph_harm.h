#pragma once

#include "error.h"
#include "legendre.h"
#include "mdspan.h"
#include "specfun.h"

extern "C" double cephes_poch(double x, double m);

namespace special {

inline std::complex<double> sph_harm(long m, long n, double theta, double phi) {
    double x, prefactor;
    std::complex<double> val;
    int mp;
    x = std::cos(phi);

    if (std::abs(m) > n) {
        set_error("sph_harm", SF_ERROR_ARG, "m should not be greater than n");
        return NAN;
    }
    if (n < 0) {
        set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return NAN;
    }

    if (m < 0) {
        mp = -m;
        prefactor = std::pow(-1, mp) * cephes_poch(n + mp + 1, -2 * mp);
    } else {
        mp = m;
    }

    val = pmv(mp, n, x);
    if (m < 0) {
        val *= prefactor;
    }

    val *= std::sqrt((2 * n + 1) / 4.0 / M_PI);
    val *= std::sqrt(cephes_poch(n + m + 1, -2 * m));
    val *= std::exp(std::complex<double>(0, m * theta));

    return val;
}

inline std::complex<double> sph_harm(double m, double n, double theta, double phi) {
    return sph_harm(static_cast<long>(m), static_cast<long>(n), theta, phi);
}

inline std::complex<float> sph_harm(long m, long n, float theta, float phi) {
    return static_cast<std::complex<float>>(sph_harm(m, n, static_cast<double>(theta), static_cast<double>(phi)));
}

inline std::complex<float> sph_harm(float m, float n, float theta, float phi) {
    return sph_harm(static_cast<long>(m), static_cast<long>(n), theta, phi);
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    using out_submat = std::mdspan<T, typename OutMat::extents_type, typename OutMat::layout_type>;

    long m = y.extent(0) - 1;
    long n = y.extent(1) - 1;

    std::array<ptrdiff_t, 2> suby_strides{2 * y.stride(0), 2 * y.stride(1)};

    out_submat y_real(reinterpret_cast<T *>(y.data_handle()), {y.extents(), suby_strides});
    lpmn_no_deriv(std::cos(phi), y_real);

    out_submat y_imag(reinterpret_cast<T *>(y.data_handle()) + 1, {y.extents(), suby_strides});
    for (long i = 0; i <= m; ++i) {
        for (long j = 0; j <= n; ++j) {
            y_imag(i, j) = 0;
        }
    }

    for (long i = 0; i <= m; ++i) {
        for (long j = i; j <= n; ++j) {
            y(i, j) *= std::sqrt((2 * j + 1) * cephes_poch(j + i + 1, -2 * i) / (4 * M_PI)) *
                       std::exp(std::complex<double>(0, i * theta));
        }
    }
}

} // namespace special
