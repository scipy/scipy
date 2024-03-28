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

template <typename T>
T sem_cva(T m, T q);

template <typename T>
void sem(T m, T q, T x, T *csf, T *csd);

inline std::complex<double> chyp2f1(double a, double b, double c, std::complex<double> z) {
    std::complex<double> outz;
    std::complex<double> z99;
    std::complex<double> outz99;
    int l1, l0, isfer = 0;

    l0 = ((c == floor(c)) && (c < 0));
    l1 = ((fabs(1 - z.real()) < 1e-15) && (z.imag() == 0) && (c - a - b <= 0));
    if (l0 || l1) {
        set_error("chyp2f1", SF_ERROR_OVERFLOW, NULL);
        outz.real(std::numeric_limits<double>::infinity());
        outz.imag(0.0);
        return outz;
    }
    z99 = std::complex<double>(z.real(), z.imag());
    outz99 = specfun::hygfz(a, b, c, z99, &isfer);
    outz.real(outz99.real());
    outz.imag(outz99.imag());
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
    std::complex<double> outz;
    std::complex<double> outz99;
    std::complex<double> z99(z.real(), z.imag());

    outz99 = specfun::cchg(a, b, z99);
    outz.real(outz99.real());
    outz.imag(outz99.imag());

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
void itairy(T x, T *apt, T *bpt, T *ant, T *bnt) {
    T tmp;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itairy(x, apt, bpt, ant, bnt);
    if (flag) { /* negative limit -- switch signs and roles */
        tmp = *apt;
        *apt = -*ant;
        *ant = -tmp;
        tmp = *bpt;
        *bpt = -*bnt;
        *bnt = -tmp;
    }
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

inline std::complex<double> cerf(std::complex<double> z) {
    std::complex<double> outz;
    std::complex<double> z99(z.real(), z.imag());

    std::complex<double> outz99 = specfun::cerror(z99);
    outz.real(outz99.real());
    outz.imag(outz99.imag());
    return outz;
}


/* Integrals of bessel functions */

/* int(j0(t),t=0..x) */
/* int(y0(t),t=0..x) */

template <typename T>
void it1j0y0(T x, T *j0int, T *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itjya(x, j0int, y0int);
    if (flag) {
        *j0int = -(*j0int);
        *y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* int((1-j0(t))/t,t=0..x) */
/* int(y0(t)/t,t=x..inf) */

template <typename T>
void it2j0y0(T x, T *j0int, T *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittjya(x, j0int, y0int);
    if (flag) {
        *y0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* Integrals of modified bessel functions */

template <typename T>
void it1i0k0(T x, T *i0int, T *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itika(x, i0int, k0int);
    if (flag) {
        *i0int = -(*i0int);
        *k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

template <typename T>
void it2i0k0(T x, T *i0int, T *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittika(x, i0int, k0int);
    if (flag) {
        *k0int = std::numeric_limits<T>::quiet_NaN(); /* domain error */
    }
}

/* Fresnel integrals of complex numbers */

inline void cfresnl(std::complex<double> z, std::complex<double> *zfs, std::complex<double> *zfc) {
    std::complex<double> zfd;

    specfun::cfs(z, zfs, &zfd);
    specfun::cfc(z, zfc, &zfd);
}


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
