#pragma once

#include "error.h"
#include "legendre.h"
#include "mdspan.h"
#include "specfun.h"

extern "C" double cephes_poch(double x, double m);

namespace special {

inline std::complex<double> sph_harm(long m, long n, double theta, double phi) {
    double prefactor;
    std::complex<double> val;
    int mp;

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

    val = pmv(mp, n, std::cos(phi));
    if (m < 0) {
        val *= prefactor;
    }

    val *= std::sqrt((2 * n + 1) / 4.0 / M_PI);
    val *= std::sqrt(cephes_poch(n + m + 1, -2 * m));
    val *= std::exp(std::complex<double>(0, m * theta));

    return val;
}

inline std::complex<float> sph_harm(long m, long n, float theta, float phi) {
    return static_cast<std::complex<float>>(sph_harm(m, n, static_cast<double>(theta), static_cast<double>(phi)));
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    const long m = (y.extent(0) - 1) / 2;
    const long n = y.extent(1) - 1;

    OutMat y_pos = std::submdspan(y, std::make_tuple(0, m + 1), std::full_extent);
    lpmn(std::cos(phi), y_pos);

    for (long j = 0; j <= n; ++j) {
        y(0, j) *= std::sqrt((2 * j + 1) / (4 * M_PI));
        for (long i = 1; i <= j; ++i) {
            y(i, j) *= static_cast<T>(std::sqrt((2 * j + 1) * cephes_poch(j + i + 1, -2 * i) / (4 * M_PI))) *
                       std::exp(std::complex<T>(0, i * theta));
            y(y.extent(0) - i, j) = static_cast<T>(std::pow(-1, i)) * std::conj(y(i, j));
        }
    }
}

} // namespace special
