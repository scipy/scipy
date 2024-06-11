#pragma once

#include "legendre.h"

namespace special {

template <typename T>
std::complex<T> sph_harm(unsigned int n, int m, T theta, T phi) {
    std::complex<T> y = sph_legendre_p(n, std::abs(m), phi) * std::exp(std::complex(T(0), m * theta));
    if (m < 0) {
        y *= std::pow(-1, m);
    }

    return y;
}

template <typename T, typename OutMat>
void sph_harm_all(T theta, T phi, OutMat y) {
    unsigned int m = (y.extent(0) - 1) / 2;
    unsigned int n = y.extent(1) - 1;

    sph_legendre_p(n, 0, phi, [y](unsigned int j, unsigned int i, T phi, T value) { y(i, j) = value; });

    for (unsigned int i = 1; i <= m; ++i) {
        sph_legendre_p(n, i, phi, [theta, y](unsigned int j, unsigned int i, T phi, T value) {
            y(i, j) = value * std::exp(std::complex(T(0), i * theta));
            y(y.extent(0) - i, j) = T(std::pow(-1, i)) * std::conj(y(i, j));
        });
    }
}

} // namespace special
