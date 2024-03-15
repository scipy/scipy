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

inline int cem(double m, double q, double x, double *csf, double *csd);
inline double sem_cva(double m, double q);
inline int sem(double m, double q, double x, double *csf, double *csd);

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

inline int itairy(double x, double *apt, double *bpt, double *ant, double *bnt) {
    double tmp;
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
    return 0;
}

inline double exp1(double x) {
    double out;

    out = specfun::e1xb(x);
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

inline double expi(double x) {
    double out;

    out = specfun::eix(x);
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

inline double itstruve0(double x) {
    double out;

    if (x < 0)
        x = -x;
    out = specfun::itsh0(x);
    SPECFUN_CONVINF("itstruve0", out);
    return out;
}

inline double it2struve0(double x) {
    double out;
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    out = specfun::itth0(x);
    SPECFUN_CONVINF("it2struve0", out);
    if (flag) {
        out = M_PI - out;
    }
    return out;
}

inline double itmodstruve0(double x) {
    double out;

    if (x < 0) {
        x = -x;
    }
    out = specfun::itsl0(x);
    SPECFUN_CONVINF("itmodstruve0", out);
    return out;
}

inline double ber(double x) {
    std::complex<double> Be;
    double ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("ber", Be);
    return Be.real();
}

inline double bei(double x) {
    std::complex<double> Be;
    double ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        x = -x;
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Be.real(ber);
    Be.imag(bei);
    SPECFUN_ZCONVINF("bei", Be);
    return Be.imag();
}

inline double ker(double x) {
    std::complex<double> Ke;
    double ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("ker", Ke);
    return Ke.real();
}

inline double kei(double x) {
    std::complex<double> Ke;
    double ber, bei, ger, gei, der, dei, her, hei;
    if (x < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Ke.real(ger);
    Ke.imag(gei);
    SPECFUN_ZCONVINF("kei", Ke);
    return Ke.imag();
}

inline double berp(double x) {
    std::complex<double> Bep;
    double ber, bei, ger, gei, der, dei, her, hei;
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

inline double beip(double x) {
    std::complex<double> Bep;
    double ber, bei, ger, gei, der, dei, her, hei;
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

inline double kerp(double x) {
    std::complex<double> Kep;
    double ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("kerp", Kep);
    return Kep.real();
}

inline double keip(double x) {
    std::complex<double> Kep;
    double ber, bei, ger, gei, der, dei, her, hei;

    if (x < 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    specfun::klvna(x, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
    Kep.real(her);
    Kep.imag(hei);
    SPECFUN_ZCONVINF("keip", Kep);
    return Kep.imag();
}

inline int kelvin(double x, std::complex<double> *Be, std::complex<double> *Ke, std::complex<double> *Bep,
                  std::complex<double> *Kep) {
    int flag = 0;
    double ber, bei, ger, gei, der, dei, her, hei;
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
        Ke->real(std::numeric_limits<double>::quiet_NaN());
        Ke->imag(std::numeric_limits<double>::quiet_NaN());
        Kep->real(std::numeric_limits<double>::quiet_NaN());
        Kep->imag(std::numeric_limits<double>::quiet_NaN());
    }
    return 0;
}

/* Integrals of bessel functions */

/* int(j0(t),t=0..x) */
/* int(y0(t),t=0..x) */

inline int it1j0y0(double x, double *j0int, double *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itjya(x, j0int, y0int);
    if (flag) {
        *j0int = -(*j0int);
        *y0int = std::numeric_limits<double>::quiet_NaN(); /* domain error */
    }
    return 0;
}

/* int((1-j0(t))/t,t=0..x) */
/* int(y0(t)/t,t=x..inf) */

inline int it2j0y0(double x, double *j0int, double *y0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittjya(x, j0int, y0int);
    if (flag) {
        *y0int = std::numeric_limits<double>::quiet_NaN(); /* domain error */
    }
    return 0;
}

/* Integrals of modified bessel functions */

inline int it1i0k0(double x, double *i0int, double *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::itika(x, i0int, k0int);
    if (flag) {
        *i0int = -(*i0int);
        *k0int = std::numeric_limits<double>::quiet_NaN(); /* domain error */
    }
    return 0;
}

inline int it2i0k0(double x, double *i0int, double *k0int) {
    int flag = 0;

    if (x < 0) {
        x = -x;
        flag = 1;
    }
    specfun::ittika(x, i0int, k0int);
    if (flag) {
        *k0int = std::numeric_limits<double>::quiet_NaN(); /* domain error */
    }
    return 0;
}

/* Fresnel integrals of complex numbers */

inline int cfresnl(std::complex<double> z, std::complex<double> *zfs, std::complex<double> *zfc) {
    std::complex<double> z99(z.real(), z.imag());
    std::complex<double> zfs99, zfc99, zfd;

    specfun::cfs(z99, &zfs99, &zfd);
    specfun::cfc(z99, &zfc99, &zfd);

    zfs->real(zfs99.real());
    zfs->imag(zfs99.imag());
    zfc->real(zfc99.real());
    zfc->imag(zfc99.imag());

    return 0;
}

/* Mathieu functions */
/* Characteristic values */
inline double cem_cva(double m, double q) {
    int int_m, kd = 1;
    double out;

    if ((m < 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
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
    out = specfun::cva2(kd, int_m, q);
    return out;
}

inline double sem_cva(double m, double q) {
    int int_m, kd = 4;
    double out;

    if ((m <= 0) || (m != floor(m))) {
        set_error("cem_cva", SF_ERROR_DOMAIN, NULL);
        return std::numeric_limits<double>::quiet_NaN();
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
    out = specfun::cva2(kd, int_m, q);
    return out;
}

/* Mathieu functions */
inline int cem(double m, double q, double x, double *csf, double *csd) {
    int int_m, kf = 1, sgn;
    double f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        *csf = std::numeric_limits<double>::quiet_NaN();
        *csd = std::numeric_limits<double>::quiet_NaN();
        set_error("cem", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E34 */
        if (int_m % 2 == 0) {
            sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
            cem(m, -q, 90 - x, &f, &d);
            *csf = sgn * f;
            *csd = -sgn * d;
            return 0;
        } else {
            sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
            sem(m, -q, 90 - x, &f, &d);
            *csf = sgn * f;
            *csd = -sgn * d;
            return 0;
        }
    }
    specfun::mtu0(kf, int_m, q, x, csf, csd);
    return 0;
}

inline int sem(double m, double q, double x, double *csf, double *csd) {
    int int_m, kf = 2, sgn;
    double f = 0.0, d = 0.0;
    if ((m < 0) || (m != floor(m))) {
        *csf = std::numeric_limits<double>::quiet_NaN();
        *csd = std::numeric_limits<double>::quiet_NaN();
        set_error("sem", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    if (int_m == 0) {
        *csf = 0;
        *csd = 0;
        return 0;
    }
    if (q < 0) {
        /* https://dlmf.nist.gov/28.2#E34 */
        if (int_m % 2 == 0) {
            sgn = ((int_m / 2) % 2 == 0) ? -1 : 1;
            sem(m, -q, 90 - x, &f, &d);
            *csf = sgn * f;
            *csd = -sgn * d;
            return 0;
        } else {
            sgn = ((int_m / 2) % 2 == 0) ? 1 : -1;
            cem(m, -q, 90 - x, &f, &d);
            *csf = sgn * f;
            *csd = -sgn * d;
            return 0;
        }
    }
    specfun::mtu0(kf, int_m, q, x, csf, csd);
    return 0;
}

inline int mcm1(double m, double q, double x, double *f1r, double *d1r) {
    int int_m, kf = 1, kc = 1;
    double f2r = 0.0, d2r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        *f1r = std::numeric_limits<double>::quiet_NaN();
        *d1r = std::numeric_limits<double>::quiet_NaN();
        set_error("mcm1", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    specfun::mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
    return 0;
}

inline int msm1(double m, double q, double x, double *f1r, double *d1r) {
    int int_m, kf = 2, kc = 1;
    double f2r = 0.0, d2r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        *f1r = std::numeric_limits<double>::quiet_NaN();
        *d1r = std::numeric_limits<double>::quiet_NaN();
        set_error("msm1", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    specfun::mtu12(kf, kc, int_m, q, x, f1r, d1r, &f2r, &d2r);
    return 0;
}

inline int mcm2(double m, double q, double x, double *f2r, double *d2r) {
    int int_m, kf = 1, kc = 2;
    double f1r = 0.0, d1r = 0.0;

    if ((m < 0) || (m != floor(m)) || (q < 0)) {
        *f2r = std::numeric_limits<double>::quiet_NaN();
        *d2r = std::numeric_limits<double>::quiet_NaN();
        set_error("mcm2", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
    return 0;
}

inline int msm2(double m, double q, double x, double *f2r, double *d2r) {
    int int_m, kf = 2, kc = 2;
    double f1r = 0.0, d1r = 0.0;

    if ((m < 1) || (m != floor(m)) || (q < 0)) {
        *f2r = std::numeric_limits<double>::quiet_NaN();
        *d2r = std::numeric_limits<double>::quiet_NaN();
        set_error("msm2", SF_ERROR_DOMAIN, NULL);
        return -1;
    }
    int_m = (int) m;
    specfun::mtu12(kf, kc, int_m, q, x, &f1r, &d1r, f2r, d2r);
    return 0;
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
inline int pbwa(double a, double x, double *wf, double *wd) {
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
        return 0;
    }

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
    return 0;
}

inline int pbdv(double v, double x, double *pdf, double *pdd) {
    double *dv;
    double *dp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pdf = std::numeric_limits<double>::quiet_NaN();
        *pdd = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    /* NB. Indexing of DV/DP in specfun.f:PBDV starts from 0, hence +2 */
    num = std::abs((int) v) + 2;
    dv = (double *) malloc(sizeof(double) * 2 * num);
    if (dv == NULL) {
        set_error("pbdv", SF_ERROR_OTHER, "memory allocation error");
        *pdf = std::numeric_limits<double>::quiet_NaN();
        *pdd = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    dp = dv + num;
    specfun::pbdv(x, v, dv, dp, pdf, pdd);
    free(dv);
    return 0;
}

inline int pbvv(double v, double x, double *pvf, double *pvd) {
    double *vv;
    double *vp;
    int num;

    if (isnan(v) || isnan(x)) {
        *pvf = std::numeric_limits<double>::quiet_NaN();
        *pvd = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    /* NB. Indexing of DV/DP in specfun.f:PBVV starts from 0, hence +2 */
    num = std::abs((int) v) + 2;
    vv = (double *) malloc(sizeof(double) * 2 * num);
    if (vv == NULL) {
        set_error("pbvv", SF_ERROR_OTHER, "memory allocation error");
        *pvf = std::numeric_limits<double>::quiet_NaN();
        *pvd = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    vp = vv + num;
    specfun::pbvv(x, v, vv, vp, pvf, pvd);
    free(vv);
    return 0;
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

inline int prolate_aswfa(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    int kd = 1;
    int int_m, int_n;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_aswfa", SF_ERROR_DOMAIN, NULL);
        *s1f = std::numeric_limits<double>::quiet_NaN();
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::aswfa(x, int_m, int_n, c, kd, cv, s1f, s1d);
    return 0;
}

inline int oblate_aswfa(double m, double n, double c, double cv, double x, double *s1f, double *s1d) {
    int kd = -1;
    int int_m, int_n;

    if ((x >= 1) || (x <= -1) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_aswfa", SF_ERROR_DOMAIN, NULL);
        *s1f = std::numeric_limits<double>::quiet_NaN();
        *s1d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::aswfa(x, int_m, int_n, c, kd, cv, s1f, s1d);
    return 0;
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

inline int prolate_radial1(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    int kf = 1;
    double r2f = 0.0, r2d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_radial1", SF_ERROR_DOMAIN, NULL);
        *r1f = std::numeric_limits<double>::quiet_NaN();
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::rswfp(int_m, int_n, c, x, cv, kf, r1f, r1d, &r2f, &r2d);
    return 0;
}

inline int prolate_radial2(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    int kf = 2;
    double r1f = 0.0, r1d = 0.0;
    int int_m, int_n;

    if ((x <= 1.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("prolate_radial2", SF_ERROR_DOMAIN, NULL);
        *r2f = std::numeric_limits<double>::quiet_NaN();
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::rswfp(int_m, int_n, c, x, cv, kf, &r1f, &r1d, r2f, r2d);
    return 0;
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

inline int oblate_radial1(double m, double n, double c, double cv, double x, double *r1f, double *r1d) {
    int kf = 1;
    double r2f = 0.0, r2d = 0.0;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_radial1", SF_ERROR_DOMAIN, NULL);
        *r1f = std::numeric_limits<double>::quiet_NaN();
        *r1d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::rswfo(int_m, int_n, c, x, cv, kf, r1f, r1d, &r2f, &r2d);
    return 0;
}

inline int oblate_radial2(double m, double n, double c, double cv, double x, double *r2f, double *r2d) {
    int kf = 2;
    double r1f = 0.0, r1d = 0.0;
    int int_m, int_n;

    if ((x < 0.0) || (m < 0) || (n < m) || (m != floor(m)) || (n != floor(n))) {
        set_error("oblate_radial2", SF_ERROR_DOMAIN, NULL);
        *r2f = std::numeric_limits<double>::quiet_NaN();
        *r2d = std::numeric_limits<double>::quiet_NaN();
        return 0;
    }
    int_m = (int) m;
    int_n = (int) n;
    specfun::rswfo(int_m, int_n, c, x, cv, kf, &r1f, &r1d, r2f, r2d);
    return 0;
}

inline int modified_fresnel_plus(double x, std::complex<double> *Fplus, std::complex<double> *Kplus) {
    int ks = 0;
    double fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

    specfun::ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
    Fplus->real(fr);
    Fplus->imag(fi);
    Kplus->real(gr);
    Kplus->imag(gi);
    return 0;
}

inline int modified_fresnel_minus(double x, std::complex<double> *Fminus, std::complex<double> *Kminus) {
    int ks = 1;
    double fr = 0.0, gr = 0.0, fi = 0.0, gi = 0.0, fa = 0.0, ga = 0.0, fm = 0.0, gm = 0.0;

    specfun::ffk(ks, x, &fr, &fi, &fm, &fa, &gr, &gi, &gm, &ga);
    Fminus->real(fr);
    Fminus->imag(fi);
    Kminus->real(gr);
    Kminus->imag(gi);
    return 0;
}

} // namespace special
