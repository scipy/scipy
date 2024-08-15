/* Translated into C++ by SciPy developers in 2024.
 *
 * Original author: Takuma Yoshimura, 2024.
 */

 /*
  * Implement sin(pi * x), cos(pi * x), tan(pi * x) and cot(pi * x) for real x. 
  * Return exactly 0 for integer or half-integer input values.
  * sinpi and cospi conform to the IEEE754-2008 specification.
  * tanpi and cotpi return sinpi/cospi and cospi/sinpi.
  */

#pragma once

#include "../config.h"

namespace xsf {
namespace cephes {

    /* Compute sin(pi * x). */
    template <typename T>
    XSF_HOST_DEVICE T sinpi(T x) {
        T s = std::copysign(1.0, x);
        x = std::abs(x);

        T r = x - std::ldexp(std::floor(std::ldexp(x, -1)), 1);
        s = r <= 1.0 ? s : -s;
        r = r - std::floor(r);
        r = r <= 0.5 ? r : (1.0 - r);

        T y = std::copysign(std::sin(M_PI * r), s);
 
        return y;
    }

    /* Compute cos(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T cospi(T x) {
        x = std::abs(x);

        T r = x - std::ldexp(std::floor(std::ldexp(x, -1)), 1);
        r = r <= 1.0 ? (0.5 - r) : (r - 1.5);

        T y = std::sin(M_PI * r);

        return y;
    }

    /* Compute tan(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T tanpi(T x) {
        T s = std::copysign(1.0, x);
        x = std::abs(x);

        T r = x - std::ldexp(std::floor(std::ldexp(x, -1)), 1);
        s = r <= 1.0 ? s : -s;

        T rs = r - std::floor(r);
        rs = rs <= 0.5 ? rs : (1.0 - rs);
        T rc = r <= 1.0 ? (0.5 - r) : (r - 1.5);

        T y = std::copysign(std::sin(M_PI * rs), s) / std::sin(M_PI * rc);

        return y;
    }

    /* Compute cot(pi * x) */
    template <typename T>
    XSF_HOST_DEVICE T cotpi(T x) {
        T s = std::copysign(1.0, x);
        x = std::abs(x);

        T r = x - std::ldexp(std::floor(std::ldexp(x, -1)), 1);
        s = r <= 1.0 ? s : -s;

        T rs = r - std::floor(r);
        rs = rs <= 0.5 ? rs : (1.0 - rs);
        T rc = r <= 1.0 ? (0.5 - r) : (r - 1.5);

        T y = std::sin(M_PI * rc) / std::copysign(std::sin(M_PI * rs), s);

        return y;
    }
} // namespace cephes
} // namespace xsf
