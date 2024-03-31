#pragma once

#include "error.h"
#include "specfun/specfun.h"

#define SPECFUN_ZCONVINF(func, z)                                                                                      \
    do {                                                                                                               \
        if ((double) (z).real() == (double) 1.0e300) {                                                                 \
            set_error(func, SF_ERROR_OVERFLOW, NULL);                                                                  \
            (z).real(std::numeric_limits<double>::infinity());                                                         \
        }                                                                                                              \
        if ((double) (z).real() == (double) -1.0e300) {                                                                \
            set_error(func, SF_ERROR_OVERFLOW, NULL);                                                                  \
            (z).real(-std::numeric_limits<double>::infinity());                                                        \
        }                                                                                                              \
    } while (0)
#define SPECFUN_CONVINF(func, x)                                                                                       \
    do {                                                                                                               \
        if ((double) (x) == (double) 1.0e300) {                                                                        \
            set_error(func, SF_ERROR_OVERFLOW, NULL);                                                                  \
            (x) = std::numeric_limits<double>::infinity();                                                             \
        }                                                                                                              \
        if ((double) (x) == (double) -1.0e300) {                                                                       \
            set_error(func, SF_ERROR_OVERFLOW, NULL);                                                                  \
            (x) = -std::numeric_limits<double>::infinity();                                                            \
        }                                                                                                              \
    } while (0)

namespace special {

inline std::complex<double> chyp2f1(double a, double b, double c, std::complex<double> z) {
    int l0 = ((c == floor(c)) && (c < 0));
    int l1 = ((fabs(1 - z.real()) < 1e-15) && (z.imag() == 0) && (c - a - b <= 0));
    if (l0 || l1) {
        set_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
        return std::numeric_limits<double>::infinity();
    }

    int isfer = 0;
    std::complex<double> outz = specfun::hygfz(a, b, c, z, &isfer);
    if (isfer == 3) {
        set_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
        outz.real(std::numeric_limits<double>::infinity());
        outz.imag(0.0);
    } else if (isfer == 5) {
        set_error("chyp2f1", SF_ERROR_LOSS, NULL);
    } else if (isfer != 0) {
        set_error("chyp2f1", static_cast<sf_error_t>(isfer), NULL);
        outz.real(std::numeric_limits<double>::quiet_NaN());
        outz.imag(std::numeric_limits<double>::quiet_NaN());
    }
    return outz;
}

inline std::complex<double> chyp1f1(double a, double b, std::complex<double> z) {
    std::complex<double> outz = specfun::cchg(a, b, z);
    if (outz.real() == 1e300) {
        set_error("chyp1f1", SF_ERROR_OVERFLOW, NULL);
        outz.real(std::numeric_limits<double>::infinity());
    }

    return outz;
}

inline double hypu(double a, double b, double x) {
    double out;
    int md; /* method code --- not returned */
    int isfer = 0;

    out = specfun::chgu(x, a, b, &md, &isfer);
    if (out == 1e300) {
        set_error("hypU", SF_ERROR_OVERFLOW, NULL);
        out = std::numeric_limits<double>::infinity();
    }
    if (isfer == 6) {
        set_error("hypU", SF_ERROR_NO_RESULT, NULL);
        out = std::numeric_limits<double>::quiet_NaN();
    } else if (isfer != 0) {
        set_error("hypU", static_cast<sf_error_t>(isfer), NULL);
        out = std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

inline double hyp1f1(double a, double b, double x) {
    double outy;

    outy = specfun::chgm(x, a, b);
    return outy;
}

template <typename T>
T exp1(T x) {
    T out = specfun::e1xb(x);
    SPECFUN_CONVINF("exp1", out);
    return out;
}

template <typename T>
std::complex<T> exp1(std::complex<T> z) {
    std::complex<T> outz = specfun::e1z(z);
    SPECFUN_ZCONVINF("exp1", outz);
    return outz;
}

template <typename T>
T expi(T x) {
    T out = specfun::eix(x);
    SPECFUN_CONVINF("expi", out);
    return out;
}

template <typename T>
std::complex<T> expi(std::complex<T> z) {
    std::complex<T> outz = specfun::eixz(z);
    SPECFUN_ZCONVINF("cexpi", outz);
    return outz;
}

inline std::complex<double> cerf(std::complex<double> z) { return specfun::cerror(z); }

inline double pmv(double m, double v, double x) {
    int int_m;
    double out;

    if (m != floor(m)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    out = specfun::lpmv(x, int_m, v);
    SPECFUN_CONVINF("pmv", out);
    return out;
}

} // namespace special
