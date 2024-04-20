#pragma once

#include "error.h"
#include "legendre.h"

namespace special {

template <typename T>
std::complex<T> sph_harm(long m, long n, T theta, T phi) {
    if (n < 0) {
        set_error("sph_harm", SF_ERROR_ARG, "n should not be negative");
        return std::numeric_limits<T>::quiet_NaN();
    }

    long m_abs = std::abs(m);
    if (m_abs > n) {
        return 0;
    }

    std::complex<T> y = sph_legendre(m_abs, n, phi) * std::exp(std::complex(T(0), m * theta));
    if (m < 0) {
        y *= std::pow(-1, m_abs);
    }

    return y;
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    long m = (y.extent(0) - 1) / 2;
    long n = y.extent(1) - 1;

    sph_legendre(0, n, phi, [y](long i, long j, T phi, T value) { y(i, j) = value; });

    for (long i = 1; i <= m; ++i) {
        sph_legendre(i, n, phi, [theta, y](long i, long j, T phi, T value) {
            y(i, j) = value * std::exp(std::complex(T(0), i * theta));
            y(y.extent(0) - i, j) = T(std::pow(-1, i)) * std::conj(y(i, j));
        });
    }
}

} // namespace special
