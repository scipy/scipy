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

template <>
inline void itairy(float xf, float *aptf, float *bptf, float *antf, float *bntf) {
    double x = xf;
    double apt;
    double bpt;
    double ant;
    double bnt;
    itairy(x, &apt, &bpt, &ant, &bnt);

    *aptf = apt;
    *bptf = bpt;
    *antf = ant;
    *bntf = bnt;
}

template <typename T>
T exp1(T x) {
    T out = specfun::e1xb(x);
    SPECFUN_CONVINF("exp1", out);
    return out;
}

inline std::complex<double> cexp1(std::complex<double> z) {
    std::complex<double> outz;
    std::complex<double> z99(z.real(), z.imag());

    std::complex<double> outz99 = specfun::e1z(z99);
    outz.real(outz99.real());
    outz.imag(outz99.imag());
    SPECFUN_ZCONVINF("cexp1", outz);
    return outz;
}

template <typename T>
T expi(T x) {
    T out = specfun::eix(x);
    SPECFUN_CONVINF("expi", out);
    return out;
}

inline std::complex<double> cexpi(std::complex<double> z) {
    std::complex<double> outz;
    std::complex<double> z99(z.real(), z.imag());

    std::complex<double> outz99 = specfun::eixz(z99);
    outz.real(outz99.real());
    outz.imag(outz99.imag());
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

template <typename T>
T itstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = specfun::itsh0(x);
    SPECFUN_CONVINF("itstruve0", out);
    return out;
}

template <typename T>
T it2struve0(T x) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }

    T out = specfun::itth0(x);
    SPECFUN_CONVINF("it2struve0", out);
    if (flag) {
        out = M_PI - out;
    }
    return out;
}

template <typename T>
T itmodstruve0(T x) {
    if (x < 0) {
        x = -x;
    }

    T out = specfun::itsl0(x);
    SPECFUN_CONVINF("itmodstruve0", out);
    return out;
}

template <typename T>
T ber(T x) {
    std::complex<T> Be;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("ber", Be);
    return Be.real();
}

template <>
inline float ber(float xf) {
    double x = xf;
    return ber(x);
}

template <typename T>
T bei(T x) {
    std::complex<T> Be;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("bei", Be);
    return Be.imag();
}

template <>
inline float bei(float xf) {
    double x = xf;
    return bei(x);
}

template <typename T>
inline T ker(T x) {
    std::complex<T> Ke;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("ker", Ke);
    return Ke.real();
}

template <>
inline float ker(float xf) {
    double x = xf;
    return ker(x);
}

template <typename T>
T kei(T x) {
    std::complex<T> Ke;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("kei", Ke);
    return Ke.imag();
}

template <>
float kei(float x) {
    return kei(static_cast<double>(x));
}

template <typename T>
T berp(T x) {
    std::complex<T> Bep;
    T ber, bei, ger, gei, der, dei, her, hei;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Bep.real(der);
    Bep.imag(dei);
    SPECFUN_ZCONVINF("berp", Bep);
    if (flag) {
        return -Bep.real();
    }
    return Bep.real();
}

template <>
inline float berp(float xf) {
    double x = xf;
    return berp(x);
}

template <typename T>
T beip(T x) {
    std::complex<T> Bep;
    T ber, bei, ger, gei, der, dei, her, hei;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Bep.real(der);
    Bep.imag(dei);
    SPECFUN_ZCONVINF("beip", Bep);
    if (flag) {
        return -Bep.imag();
    }
    return Bep.imag();
}

template <>
float beip(float xf) {
    double x = xf;
    return beip(x);
}

template <typename T>
T kerp(T x) {
    std::complex<T> Kep;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("kerp", Kep);
    return Kep.real();
}

template <>
inline float kerp(float xf) {
    double x = xf;
    return kerp(x);
}

template <typename T>
T keip(T x) {
    std::complex<T> Kep;
    T ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<T>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("keip", Kep);
    return Kep.imag();
}

template <>
inline float keip(float xf) {
    double x = xf;
    return keip(x);
}

template <typename T>
void kelvin(T x, std::complex<T> *Be, std::complex<T> *Ke, std::complex<T> *Bep, std::complex<T> *Kep) {
    int flag = 0;
    T ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        x = -x;
        flag = 1;
    }

    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be->real(ber);
    Be->imag(bei);
    Ke->real(ger);
    Ke->imag(gei);
    Bep->real(der);
    Bep->imag(dei);
    Kep->real(her);
    Kep->imag(hei);

    SPECFUN_ZCONVINF("klvna", *Be);
    SPECFUN_ZCONVINF("klvna", *Ke);
    SPECFUN_ZCONVINF("klvna", *Bep);
    SPECFUN_ZCONVINF("klvna", *Kep);
    if (flag) {
        Bep->real(-Bep->real());
        Bep->imag(-Bep->imag());
        Ke->real(std::numeric_limits<T>::quiet_NaN());
        Ke->imag(std::numeric_limits<T>::quiet_NaN());
        Kep->real(std::numeric_limits<T>::quiet_NaN());
        Kep->imag(std::numeric_limits<T>::quiet_NaN());
    }
}

template <>
inline void kelvin(float xf, std::complex<float> *Bef, std::complex<float> *Kef, std::complex<float> *Bepf,
                   std::complex<float> *Kepf) {
    double x = xf;
    std::complex<double> Be;
    std::complex<double> Ke;
    std::complex<double> Bep;
    std::complex<double> Kep;
    kelvin(x, &Be, &Ke, &Bep, &Kep);

    *Bef = Be;
    *Kef = Ke;
    *Bepf = Bep;
    *Kepf = Kep;
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

template <>
inline void it1j0y0(float xf, float *j0intf, float *y0intf) {
    double x = xf;
    double j0int;
    double y0int;
    it1j0y0(x, &j0int, &y0int);

    *j0intf = j0int;
    *y0intf = y0int;
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

template <>
inline void it2j0y0(float xf, float *j0intf, float *y0intf) {
    double x = xf;
    double j0int, y0int;
    it2j0y0(x, &j0int, &y0int);

    *j0intf = j0int;
    *y0intf = y0int;
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

template <>
inline void it1i0k0(float xf, float *i0intf, float *k0intf) {
    double x = xf;
    double i0int;
    double k0int;
    it1i0k0(x, &i0int, &k0int);

    *i0intf = i0int;
    *k0intf = k0int;
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

template <>
inline void it2i0k0(float xf, float *i0intf, float *k0intf) {
    double x = xf;
    double i0int, k0int;
    it2i0k0(x, &i0int, &k0int);

    *i0intf = i0int;
    *k0intf = k0int;
}

/* Fresnel integrals of complex numbers */

inline void cfresnl(std::complex<double> z, std::complex<double> *zfs, std::complex<double> *zfc) {
    std::complex<double> zfd;

    specfun::cfs(z, zfs, &zfd);
    specfun::cfc(z, zfc, &zfd);
}

/* Mathieu functions */
/* Characteristic values */
template <typename T>
T cem_cva(T m, T q) {
    int int_m, kd = 1;

    if ((m < 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E26 */
        if (int_m % 2 == 0) {
            return cem_cva(m, -q);
        } else {
            return sem_cva(m, -q);
        }
    }

    if (int_m % 2) {
        kd = 2;
    }
    return specfun::cva2(kd, int_m, q);
}

template <>
inline float cem_cva(float mf, float qf) {
    double m = mf;
    double q = qf;

    return cem_cva(m, q);
}

template <typename T>
T sem_cva(T m, T q) {
    int int_m, kd = 4;

    if ((m <= 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<T>::quiet_NaN();
    }
    int_m = (int) m;
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E26 */
        if (int_m % 2 == 0) {
            return sem_cva(m, -q);
        } else {
            return cem_cva(m, -q);
        }
    }
    if (int_m % 2) {
        kd = 3;
    }
    return specfun::cva2(kd, int_m, q);
}

template <>
inline float sem_cva(float mf, float qf) {
    double m = mf;
    double q = qf;

    return sem_cva(m, q);
}

/* Mathieu functions */
template <typename T>
inline void cem(T m, T q, T x, T *csf, T *csd) {
    int int_m, kf = 1, sgn;
    T f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        *csf = std::numeric_limits<T>::quiet_NaN();
        *csd = std::numeric_limits<T>::quiet_NaN();
        set_error("cem", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        if (q < 0) {
            /* https://dlmf.nist.gov/28.2#E34 */
            if (int_m % 2 == 0) {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                cem(m, -q, 90 - x, &f, &d);
                *csf = sgn * f;
                *csd = -sgn * d;

            } else {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                sem(m, -q, 90 - x, &f, &d);
                *csf = sgn * f;
                *csd = -sgn * d;
            }
        } else {
            specfun::mtu0(kf, int_m, q, x, csf, csd);
        }
    }
}

template <>
inline void cem(float mf, float qf, float xf, float *csff, float *csdf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double csf;
    double csd;
    cem(m, q, x, &csf, &csd);

    *csff = csf;
    *csdf = csd;
}

template <typename T>
void sem(T m, T q, T x, T *csf, T *csd) {
    int int_m, kf = 2, sgn;
    T f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        *csf = std::numeric_limits<T>::quiet_NaN();
        *csd = std::numeric_limits<T>::quiet_NaN();
        set_error("sem", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        if (int_m == 0) {
            *csf = 0;
            *csd = 0;
        } else if (q < 0) {
            /* https://dlmf.nist.gov/28.2#E34 */
            if (int_m % 2 == 0) {
                sgn = ((int_m / 2) % 2 == 0) ? -1 : 1;
                sem(m, -q, 90 - x, &f, &d);
                *csf = sgn * f;
                *csd = -sgn * d;
            } else {
                sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
                cem(m, -q, 90 - x, &f, &d);
                *csf = sgn * f;
                *csd = -sgn * d;
            }
        } else {
            specfun::mtu0(kf, int_m, q, x, csf, csd);
        }
    }
}

template <>
inline void sem(float mf, float qf, float xf, float *csff, float *csdf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double csf;
    double csd;
    sem(m, q, x, &csf, &csd);

    *csff = csf;
    *csdf = csd;
}

template <typename T>
void mcm1(T m, T q, T x, T *f1r, T *d1r) {
    int int_m, kf = 1, kc = 1;
    T f2r = 0.0, d2r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        *f1r = std::numeric_limits<T>::quiet_NaN();
        *d1r = std::numeric_limits<T>::quiet_NaN();
        set_error("mcm1", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
    }
}

template <>
inline void mcm1(float mf, float qf, float xf, float *f1rf, float *d1rf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double f1r;
    double d1r;
    mcm1(m, q, x, &f1r, &d1r);

    *f1rf = f1r;
    *d1rf = d1r;
}

template <typename T>
void msm1(T m, T q, T x, T *f1r, T *d1r) {
    int int_m, kf = 2, kc = 1;
    T f2r = 0.0, d2r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        *f1r = std::numeric_limits<T>::quiet_NaN();
        *d1r = std::numeric_limits<T>::quiet_NaN();
        set_error("msm1", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
    }
}

template <>
inline void msm1(float mf, float qf, float xf, float *f1rf, float *d1rf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double f1r;
    double d1r;
    msm1(m, q, x, &f1r, &d1r);

    *f1rf = f1r;
    *d1rf = d1r;
}

template <typename T>
void mcm2(T m, T q, T x, T *f2r, T *d2r) {
    int int_m, kf = 1, kc = 2;
    T f1r = 0.0, d1r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        *f2r = std::numeric_limits<T>::quiet_NaN();
        *d2r = std::numeric_limits<T>::quiet_NaN();
        set_error("mcm2", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
    }
}

template <>
inline void mcm2(float mf, float qf, float xf, float *f2rf, float *d2rf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double f2r;
    double d2r;
    mcm2(m, q, x, &f2r, &d2r);

    *f2rf = f2r;
    *d2rf = d2r;
}

template <typename T>
inline void msm2(T m, T q, T x, T *f2r, T *d2r) {
    int int_m, kf = 2, kc = 2;
    T f1r = 0.0, d1r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        *f2r = std::numeric_limits<T>::quiet_NaN();
        *d2r = std::numeric_limits<T>::quiet_NaN();
        set_error("msm2", SF_ERROR_DOMAIN, NULL);
    } else {
        int_m = (int) m;
        specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
    }
}

template <>
inline void msm2(float mf, float qf, float xf, float *f2rf, float *d2rf) {
    double m = mf;
    double q = qf;
    double x = xf;
    double f2r;
    double d2r;
    msm2(m, q, x, &f2r, &d2r);

    *f2rf = f2r;
    *d2rf = d2r;
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

/*
 * If x > 0 return w1f and w1d. Otherwise set x = abs(x) and return
 * w2f and -w2d.
 */
inline void pbwa(double a, double x, double *wf, double *wd) {
    int flag = 0;
    double w1f = 0.0, w1d = 0.0, w2f = 0.0, w2d = 0.0;

    if (x < -5 || x > 5 || a < -5 || a > 5) {
        /*
         * The Zhang and Jin implementation only uses Taylor series;
         * return NaN outside of the range which they are accurate.
         */
        *wf = std::numeric_limits<double>::quiet_NaN();
        *wd = std::numeric_limits<double>::quiet_NaN();
        set_error("pbwa", SF_ERROR_LOSS, NULL);
    } else {
        if (x < 0) {
            x = -x;
            flag = 1;
        }
        specfun::pbwa(a, x, &w1f, &w1d, &w2f, &w2d);
        if (flag) {
            *wf = w2f;
            *wd = -w2d;
        } else {
            *wf = w1f;
            *wd = w1d;
        }
    }
}

inline void pbdv(double v, double x, double *pdf, double *pdd) {
    double *dv;
    double *dp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pdf = std::numeric_limits<double>::quiet_NaN();
        *pdd = std::numeric_limits<double>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        dv = (double *) malloc(sizeof(double) * 2 * num);
        if (dv == NULL) {
            set_error("pbdv", SF_ERROR_OTHER, "memory allocation error");
            *pdf = std::numeric_limits<double>::quiet_NaN();
            *pdd = std::numeric_limits<double>::quiet_NaN();
        } else {
            dp = dv + num;
            specfun::pbdv(x, v, dv, dp, pdf, pdd);
            free(dv);
        }
    }
}

inline void pbvv(double v, double x, double *pvf, double *pvd) {
    double *vv;
    double *vp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pvf = std::numeric_limits<double>::quiet_NaN();
        *pvd = std::numeric_limits<double>::quiet_NaN();
    } else {
        /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
        num = std::abs((int) v) + 2;
        vv = (double *) malloc(sizeof(double) * 2 * num);
        if (vv == NULL) {
            set_error("pbvv", SF_ERROR_OTHER, "memory allocation error");
            *pvf = std::numeric_limits<double>::quiet_NaN();
            *pvd = std::numeric_limits<double>::quiet_NaN();
        } else {
            vp = vv + num;
            specfun::pbvv(x, v, vv, vp, pvf, pvd);
            free(vv);
        }
    }
}

inline double prolate_segv(double m, double n, double c) {
    int kd = 1;
    int int_m, int_n;
    double cv = 0.0, *eg;

    if ((m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("prolate_segv", SF_ERROR_OTHER, "memory allocation error");
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    return cv;
}

inline double oblate_segv(double m, double n, double c) {
    int kd = -1;
    int int_m, int_n;
    double cv = 0.0, *eg;

    if ((m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("oblate_segv", SF_ERROR_OTHER, "memory allocation error");
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    free(eg);
    return cv;
}

double prolate_aswfa_nocv(double m, double n, double c, double x, double *s1d) {
    int kd = 1;
    int int_m, int_n;
    double cv = 0.0, s1f, *eg;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("prolate_aswfa_nocv", SF_ERROR_DOMAIN, NULL);
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("prolate_aswfa_nocv", SF_ERROR_OTHER, "memory allocation error");
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::aswfa(x, int_m, int_n, c, kd, cv, &s1f, s1d);
    free(eg);
    return s1f;
}

double oblate_aswfa_nocv(double m, double n, double c, double x, double *s1d) {
    int kd = -1;
    int int_m, int_n;
    double cv = 0.0, s1f = 0.0, *eg;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("oblate_aswfa_nocv", SF_ERROR_DOMAIN, NULL);
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("oblate_aswfa_nocv", SF_ERROR_OTHER, "memory allocation error");
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::aswfa(x, int_m, int_n, c, kd, cv, &s1f, s1d);
    free(eg);
    return s1f;
}

inline void prolate_aswfa(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_aswfa", SF_ERROR_DOMAIN, NULL);
        *s1f = std::numeric_limits<double>::quiet_NaN();
        *s1d = std::numeric_limits<double>::quiet_NaN();
    } else {
        specfun::aswfa(x, static_cast<int>(m), static_cast<int>(n), c, 1, cv, s1f, s1d);
    }
}

inline void oblate_aswfa(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_aswfa", SF_ERROR_DOMAIN, NULL);
        *s1f = std::numeric_limits<double>::quiet_NaN();
        *s1d = std::numeric_limits<double>::quiet_NaN();
    } else {
        specfun::aswfa(x, static_cast<int>(m), static_cast<int>(n), c, -1, cv, s1f, s1d);
    }
}

inline double prolate_radial1_nocv(double m, double n, double c, double x, double *r1d) {
    int kf = 1, kd = 1;
    double r2f = 0.0, r2d = 0.0, r1f = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("prolate_radial1_nocv", SF_ERROR_DOMAIN, NULL);
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("prolate_radial1_nocv", SF_ERROR_OTHER, "memory allocation error");
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, r1d, &r2f, &r2d);
    free(eg);
    return r1f;
}

inline double prolate_radial2_nocv(double m, double n, double c, double x, double *r2d) {
    int kf = 2, kd = 1;
    double r1f = 0.0, r1d = 0.0, r2f = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("prolate_radial2_nocv", SF_ERROR_DOMAIN, NULL);
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("prolate_radial2_nocv", SF_ERROR_OTHER, "memory allocation error");
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, r2d);
    free(eg);
    return r2f;
}

inline void prolate_radial1(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    int kf = 1;
    double r2f = 0.0, r2d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_radial1", SF_ERROR_DOMAIN, NULL);
        *r1f = std::numeric_limits<double>::quiet_NaN();
        *r1d = std::numeric_limits<double>::quiet_NaN();
    } else {
        int_m = (int) m;
        int_n = (int) n;
        specfun::rswfp(int_m, int_n, c, x, cv, kf, r1f, r1d, &r2f, &r2d);
    }
}

inline void prolate_radial2(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    int kf = 2;
    double r1f = 0.0, r1d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_radial2", SF_ERROR_DOMAIN, NULL);
        *r2f = std::numeric_limits<double>::quiet_NaN();
        *r2d = std::numeric_limits<double>::quiet_NaN();
    } else {
        int_m = (int) m;
        int_n = (int) n;
        specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, r2f, r2d);
    }
}

inline double oblate_radial1_nocv(double m, double n, double c, double x, double *r1d) {
    int kf = 1, kd = -1;
    double r2f = 0.0, r2d = 0.0, r1f = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("oblate_radial1_nocv", SF_ERROR_DOMAIN, NULL);
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("oblate_radial1_nocv", SF_ERROR_OTHER, "memory allocation error");
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::rswfo(int_m, int_n, c, x, cv, kf, &r1f, r1d, &r2f, &r2d);
    free(eg);
    return r1f;
}

inline double oblate_radial2_nocv(double m, double n, double c, double x, double *r2d) {
    int kf = 2, kd = -1;
    double r1f = 0.0, r1d = 0.0, r2f = 0.0, cv = 0.0, *eg;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n)) || ((n - m) > 198)) {
        set_error("oblate_radial2_nocv", SF_ERROR_DOMAIN, NULL);
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    int_m = (int) m;
    int_n = (int) n;
    eg = (double *) malloc(sizeof(double) * (n - m + 2));
    if (eg == NULL) {
        set_error("oblate_radial2_nocv", SF_ERROR_OTHER, "memory allocation error");
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::segv(int_m, int_n, c, kd, &cv, eg);
    specfun::rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, &r2f, r2d);
    free(eg);
    return r2f;
}

inline void oblate_radial1(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    int kf = 1;
    double r2f = 0.0, r2d = 0.0;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_radial1", SF_ERROR_DOMAIN, NULL);
        *r1f = std::numeric_limits<double>::quiet_NaN();
        *r1d = std::numeric_limits<double>::quiet_NaN();
    } else {
        specfun::rswfo(static_cast<int>(m), static_cast<int>(n), c, x, cv, kf, r1f, r1d, &r2f, &r2d);
    }
}

inline void oblate_radial2(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    int kf = 2;
    double r1f = 0.0, r1d = 0.0;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_radial2", SF_ERROR_DOMAIN, NULL);
        *r2f = std::numeric_limits<double>::quiet_NaN();
        *r2d = std::numeric_limits<double>::quiet_NaN();
    } else {
        specfun::rswfo(static_cast<int>(m), static_cast<int>(n), c, x, cv, kf, &r1f, &r1d, r2f, r2d);
    }
}

template <typename T>
void modified_fresnel_plus(T x, std::complex<T> *Fplus, std::complex<T> *Kplus) {
    int ks = 0;
    T fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

    specfun::ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
    Fplus->real(fr);
    Fplus->imag(fi);
    Kplus->real(gr);
    Kplus->imag(gi);
}

template <>
inline void modified_fresnel_plus(float xf, std::complex<float> *Fplusf, std::complex<float> *Kplusf) {
    double x = xf;
    std::complex<double> Fplus = *Fplusf;
    std::complex<double> Kplus = *Kplusf;
    modified_fresnel_plus(x, &Fplus, &Kplus);

    *Fplusf = Fplus;
    *Kplusf = Kplus;
}

template <typename T>
void modified_fresnel_minus(T x, std::complex<T> *Fminus, std::complex<T> *Kminus) {
    int ks = 1;
    T fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

    specfun::ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
    Fminus->real(fr);
    Fminus->imag(fi);
    Kminus->real(gr);
    Kminus->imag(gi);
}

template <>
inline void modified_fresnel_minus(float xf, std::complex<float> *Fminusf, std::complex<float> *Kminusf) {
    double x = xf;
    std::complex<double> Fminus = *Fminusf;
    std::complex<double> Kminus = *Kminusf;
    modified_fresnel_minus(x, &Fminus, &Kminus);

    *Fminusf = Fminus;
    *Kminusf = Kminus;
}

} // namespace special
