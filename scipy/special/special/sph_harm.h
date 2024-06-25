#pragma once

#include "error.h"
#include "legendre.h"
#include "specfun.h"
#include "third_party/kokkos/mdspan.hpp"

#include "cephes/poch.h"

namespace special {

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (n < 0) {
        set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return NAN;
    }

    long m_abs = std::abs(m);
    if (m_abs > n) {
        return 0;
    }

    std::complex<T> val = pmv(m_abs, n, std::cos(phi));
    if (m < 0) {
        val *= std::pow(-1, m_abs) * cephes::poch(n + m_abs + 1, -2 * m_abs);
    }

    val *= std::sqrt((2 * n + 1) * cephes::poch(n + m + 1, -2 * m) / (4 * M_PI));
    val *= std::exp(std::complex(static_cast<T>(0), m * theta));

    return val;
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    long m = (y.extent(0) - 1) / 2;
    long n = y.extent(1) - 1;

    OutMat y_pos = std::submdspan(y, std::make_tuple(0, m + 1), std::full_extent);
    sph_legendre_all(phi, y_pos);

    for (long j = 0; j <= n; ++j) {
        for (long i = 1; i <= j; ++i) {
            y(i, j) *= std::exp(std::complex(static_cast<T>(0), i * theta));
            y(2 * m + 1 - i, j) = static_cast<T>(std::pow(-1, i)) * std::conj(y(i, j));
        }
        for (long i = j + 1; i <= m; ++i) {
            y(2 * m + 1 - i, j) = 0;
        }
    }
}

} // namespace special
