#pragma once

#include "cephes/dd_real.h"
#include "trig.h"

namespace xsf {

inline double log1p(double x) { return cephes::log1p(x); }

inline float log1p(float x) { return log1p(static_cast<double>(x)); }

inline std::complex<double> clog1p_ddouble(double zr, double zi) {
    double x, y;

    cephes::detail::double_double r(zr);
    cephes::detail::double_double i(zi);
    cephes::detail::double_double two(2.0);

    cephes::detail::double_double rsqr = r * r;
    cephes::detail::double_double isqr = i * i;
    cephes::detail::double_double rtwo = two * r;
    cephes::detail::double_double absm1 = rsqr + isqr;
    absm1 = absm1 + rtwo;

    x = 0.5 * log1p(static_cast<double>(absm1));
    y = atan2(zi, zr + 1.0);
    return std::complex<double>{x, y};
}

// log(z + 1) = log(x + 1 + 1j*y)
//             = log(sqrt((x+1)**2 + y**2)) + 1j*atan2(y, x+1)
//
// Using atan2(y, x+1) for the imaginary part is always okay.  The real part
// needs to be calculated more carefully.  For |z| large, the naive formula
// log(z + 1) can be used.  When |z| is small, rewrite as
//
// log(sqrt((x+1)**2 + y**2)) = 0.5*log(x**2 + 2*x +1 + y**2)
//       = 0.5 * log1p(x**2 + y**2 + 2*x)
//       = 0.5 * log1p(hypot(x,y) * (hypot(x, y) + 2*x/hypot(x,y)))
//
// This expression suffers from cancellation when x < 0 and
// y = +/-sqrt(2*fabs(x)). To get around this cancellation problem, we use
// double-double precision when necessary.
inline std::complex<double> log1p(std::complex<double> z) {
    double x, y, az, azi;

    if (!std::isfinite(std::real(z)) || !std::isfinite(std::imag(z))) {
        z = z + 1.0;
        return std::log(z);
    }

    double zr = z.real();
    double zi = z.imag();

    if (zi == 0.0 && zr >= -1.0) {
        return log1p(zr);
    }

    az = std::abs(z);
    if (az < 0.707) {
        azi = std::fabs(zi);
        if (zr < 0 && std::abs(-zr - azi * azi / 2) / (-zr) < 0.5) {
            return clog1p_ddouble(zr, zi);
        } else {
            x = 0.5 * log1p(az * (az + 2 * zr / az));
            y = atan2(zi, zr + 1.0);
            return std::complex<double>(x, y);
        }
    }

    z = z + 1.0;
    return std::log(z);
}

inline std::complex<float> log1p(std::complex<float> z) {
    return static_cast<std::complex<float>>(log1p(static_cast<std::complex<double>>(z)));
}

inline double log1pmx(double x) { return cephes::log1pmx(x); }

inline float log1pmx(float x) { return log1pmx(static_cast<double>(x)); }

template <typename T>
T xlogy(T x, T y) {
    if (x == 0 && !std::isnan(y)) {
        return 0;
    }

    return x * std::log(y);
}

template <typename T>
std::complex<T> xlogy(std::complex<T> x, std::complex<T> y) {
    if (x == T(0) && !std::isnan(std::real(y)) && !std::isnan(std::imag(y))) {
        return 0;
    }

    return x * std::log(y);
}

template <typename T>
T xlog1py(T x, T y) {
    if (x == 0 && !std::isnan(y)) {
        return 0;
    }

    return x * log1p(y);
}

template <typename T>
std::complex<T> xlog1py(std::complex<T> x, std::complex<T> y) {
    if (x == T(0) && !std::isnan(std::real(y)) && !std::isnan(std::imag(y))) {
        return 0;
    }

    return x * log1p(y);
}

} // namespace xsf
