/*

Implementation of spherical Bessel functions and modified spherical Bessel
functions of the first and second kinds, as well as their derivatives.

Author: Tadeusz Pudlik

Distributed under the same license as SciPy.

I attempt to correctly handle the edge cases (0 and infinity), but this is
tricky: the values of the functions often depend on the direction in which
the limit is taken. At zero, I follow the convention of numpy (1.9.2),
which treats zero differently depending on its type:

    >>> np.cos(0)/0
    inf
    >>> np.cos(0+0j)/(0+0j)
    inf + nan*j

So, real zero is assumed to be "positive zero", while complex zero has an
unspecified sign and produces nans.  Similarly, complex infinity is taken to
represent the "point at infinity", an ambiguity which for some functions
makes `nan` the correct return value.

Translated to C++ by SciPy developers in 2024.

*/

#pragma once

#include "amos.h"
#include "error.h"

namespace special {

template <typename T>
T sph_bessel_j(long n, T x) {
    if (std::isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_jn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if ((x == std::numeric_limits<T>::infinity()) || (x == -std::numeric_limits<T>::infinity())) {
        return 0;
    }

    if (x == 0) {
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    if ((n > 0) && (n >= x)) {
        return std::sqrt(M_PI_2 / x) * cyl_bessel_j(n + 1 / static_cast<T>(2), x);
    }

    T s0 = std::sin(x) / x;
    if (n == 0) {
        return s0;
    }

    T s1 = (s0 - std::cos(x)) / x;
    if (n == 1) {
        return s1;
    }

    T sn;
    for (int i = 0; i < n - 1; ++i) {
        sn = (2 * i + 3) * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if (std::isinf(sn)) {
            // Overflow occurred already : terminate recurrence.
            return sn;
        }
    }

    return sn;
}

template <typename T>
std::complex<T> sph_bessel_j(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_jn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == std::numeric_limits<T>::infinity() || std::real(z) == -std::numeric_limits<T>::infinity()) {
        // https://dlmf.nist.gov/10.52.E3
        if (std::imag(z) == 0) {
            return 0;
        }

        return std::complex<T>(1, 1) * std::numeric_limits<T>::infinity();
    }

    if ((std::real(z) == 0) && (std::imag(z) == 0)) {
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    std::complex<T> out = std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_j(n + 1 / static_cast<T>(2), z);
    if (std::imag(z) == 0) {
        return std::real(out); // Small imaginary part is spurious
    }

    return out;
}

template <typename T>
T sph_bessel_j_jac(long n, T z) {
    if (n == 0) {
        return -sph_bessel_j(1, z);
    }

    if (z == static_cast<T>(0)) {
        // DLMF 10.51.2 doesn't work, so use 10.51.1 to get the exact value
        if (n == 1) {
            return static_cast<T>(1) / static_cast<T>(3);
        }

        return 0;
    }

    // DLMF 10.51.2
    return sph_bessel_j(n - 1, z) - static_cast<T>(n + 1) * sph_bessel_j(n, z) / z;
}

template <typename T>
T sph_bessel_y(long n, T x) {
    T s0, s1, sn;
    int idx;

    if (isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_yn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (x < 0) {
        return std::pow(-1, n + 1) * sph_bessel_y(n, -x);
    }

    if (x == std::numeric_limits<T>::infinity() || x == -std::numeric_limits<T>::infinity()) {
        return 0;
    }

    if (x == 0) {
        return -std::numeric_limits<T>::infinity();
    }

    s0 = -cos(x) / x;
    if (n == 0) {
        return s0;
    }

    s1 = (s0 - sin(x)) / x;
    if (n == 1) {
        return s1;
    }

    for (idx = 0; idx < n - 1; ++idx) {
        sn = (2 * idx + 3) * s1 / x - s0;
        s0 = s1;
        s1 = sn;
        if (isinf(sn)) {
            // Overflow occurred already: terminate recurrence.
            return sn;
        }
    }

    return sn;
}

inline float sph_bessel_y(long n, float x) { return sph_bessel_y(n, static_cast<double>(x)); }

template <typename T>
std::complex<T> sph_bessel_y(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_yn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == 0 && std::imag(z) == 0) {
        // https://dlmf.nist.gov/10.52.E2
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::real(z) == std::numeric_limits<T>::infinity() || std::real(z) == -std::numeric_limits<T>::infinity()) {
        // https://dlmf.nist.gov/10.52.E3
        if (std::imag(z) == 0) {
            return 0;
        }

        return std::complex<T>(1, 1) * std::numeric_limits<T>::infinity();
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_y(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_y_jac(long n, T x) {
    if (n == 0) {
        return -sph_bessel_y(1, x);
    }

    return sph_bessel_y(n - 1, x) - static_cast<T>(n + 1) * sph_bessel_y(n, x) / x;
}

template <typename T>
T sph_bessel_i(long n, T x) {
    if (isnan(x)) {
        return x;
    }

    if (n < 0) {
        set_error("spherical_in", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (x == 0) {
        // https://dlmf.nist.gov/10.52.E1
        if (n == 0) {
            return 1;
        }
        return 0;
    }

    if (isinf(x)) {
        // https://dlmf.nist.gov/10.49.E8
        if (x == -std::numeric_limits<T>::infinity()) {
            return std::pow(-1, n) * std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::infinity();
    }

    return sqrt(static_cast<T>(M_PI_2) / x) * cyl_bessel_i(n + 1 / static_cast<T>(2), x);
}

template <typename T>
std::complex<T> sph_bessel_i(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_in", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(z) == 0) {
        // https://dlmf.nist.gov/10.52.E1
        if (n == 0) {
            return 1;
        }

        return 0;
    }

    if (std::isinf(std::real(z)) || std::isinf(std::imag(z))) {
        // https://dlmf.nist.gov/10.52.E5
        if (std::imag(z) == 0) {
            if (std::real(z) == -std::numeric_limits<T>::infinity()) {
                return std::pow(-1, n) * std::numeric_limits<T>::infinity();
            }

            return std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::quiet_NaN();
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_i(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_i_jac(long n, T z) {
    if (n == 0) {
        return sph_bessel_i(1, z);
    }

    if (z == static_cast<T>(0)) {
        if (n == 1) {
            return 1./3.;
        }
        else {
            return 0;
        }
    }

    return sph_bessel_i(n - 1, z) - static_cast<T>(n + 1) * sph_bessel_i(n, z) / z;
}

template <typename T>
T sph_bessel_k(long n, T z) {
    if (isnan(z)) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_kn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (z == 0) {
        return std::numeric_limits<T>::infinity();
    }

    if (isinf(z)) {
        // https://dlmf.nist.gov/10.52.E6
        if (z == std::numeric_limits<T>::infinity()) {
            return 0;
        }

        return -std::numeric_limits<T>::infinity();
    }

    return std::sqrt(M_PI_2 / z) * cyl_bessel_k(n + 1 / static_cast<T>(2), z);
}

template <typename T>
std::complex<T> sph_bessel_k(long n, std::complex<T> z) {
    if (std::isnan(std::real(z)) || std::isnan(std::imag(z))) {
        return z;
    }

    if (n < 0) {
        set_error("spherical_kn", SF_ERROR_DOMAIN, nullptr);
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::abs(z) == 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }

    if (std::isinf(std::real(z)) || std::isinf(std::imag(z))) {
        // https://dlmf.nist.gov/10.52.E6
        if (std::imag(z) == 0) {
            if (std::real(z) == std::numeric_limits<T>::infinity()) {
                return 0;
            }

            return -std::numeric_limits<T>::infinity();
        }

        return std::numeric_limits<T>::quiet_NaN();
    }

    return std::sqrt(static_cast<T>(M_PI_2) / z) * cyl_bessel_k(n + 1 / static_cast<T>(2), z);
}

template <typename T>
T sph_bessel_k_jac(long n, T x) {
    if (n == 0) {
        return -sph_bessel_k(1, x);
    }

    return -sph_bessel_k(n - 1, x) - static_cast<T>(n + 1) * sph_bessel_k(n, x) / x;
}

} // namespace special
