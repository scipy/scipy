#pragma once

#include "amos.h"
#include "cephes/airy.h"
#include "config.h"
#include "error.h"

inline int cephes_airy(float xf, float *aif, float *aipf, float *bif, float *bipf) {
    double ai;
    double aip;
    double bi;
    double bip;
    int res = xsf::cephes::airy(xf, &ai, &aip, &bi, &bip);

    *aif = ai;
    *aipf = aip;
    *bif = bi;
    *bipf = bip;
    return res;
}

namespace xsf {
namespace detail {

    template <typename T>
    void itairy(T x, T &apt, T &bpt, T &ant, T &bnt) {

        // ======================================================
        // Purpose: Compute the integrals of Airy fnctions with
        //          respect to t from 0 and x ( x ≥ 0 )
        // Input  : x   --- Upper limit of the integral
        // Output : APT --- Integration of Ai(t) from 0 and x
        //          BPT --- Integration of Bi(t) from 0 and x
        //          ANT --- Integration of Ai(-t) from 0 and x
        //          BNT --- Integration of Bi(-t) from 0 and x
        // ======================================================

        int k, l;
        T fx, gx, r, su1, su2, su3, su4, su5, su6, xp6, xe, xr1, xr2;

        const T pi = 3.141592653589793;
        const T c1 = 0.355028053887817;
        const T c2 = 0.258819403792807;
        const T sr3 = 1.732050807568877;
        const T q0 = 1.0 / 3.0;
        const T q1 = 2.0 / 3.0;
        const T q2 = 1.4142135623730951;
        const T eps = 1e-5;
        static const T a[16] = {0.569444444444444,     0.891300154320988,     0.226624344493027e+01,
                                0.798950124766861e+01, 0.360688546785343e+02, 0.198670292131169e+03,
                                0.129223456582211e+04, 0.969483869669600e+04, 0.824184704952483e+05,
                                0.783031092490225e+06, 0.822210493622814e+07, 0.945557399360556e+08,
                                0.118195595640730e+10, 0.159564653040121e+11, 0.231369166433050e+12,
                                0.358622522796969e+13};

        if (x == 0.0) {
            apt = 0.0;
            bpt = 0.0;
            ant = 0.0;
            bnt = 0.0;
        } else {
            if (std::abs(x) <= 9.25) {
                for (l = 0; l < 2; l++) {
                    x = x * std::pow(-1, l);
                    fx = x;
                    r = x;
                    for (k = 1; k < 41; k++) {
                        r = r * (3.0 * k - 2.0) / (3.0 * k + 1.0) * x / (3.0 * k) * x / (3.0 * k - 1.0) * x;
                        fx += r;
                        if (std::abs(r) < std::abs(fx) * eps) {
                            break;
                        }
                    }
                    gx = 0.5 * x * x;
                    r = gx;
                    for (k = 1; k < 41; k++) {
                        r = r * (3.0 * k - 1.0) / (3.0 * k + 2.0) * x / (3.0 * k) * x / (3.0 * k + 1.0) * x;
                        gx += r;
                        if (std::abs(r) < std::abs(gx) * eps) {
                            break;
                        }
                    }
                    ant = c1 * fx - c2 * gx;
                    bnt = sr3 * (c1 * fx + c2 * gx);
                    if (l == 0) {
                        apt = ant;
                        bpt = bnt;
                    } else {
                        ant = -ant;
                        bnt = -bnt;
                        x = -x;
                    }
                }
            } else {
                xe = x * std::sqrt(x) / 1.5;
                xp6 = 1.0 / std::sqrt(6.0 * pi * xe);
                su1 = 1.0;
                r = 1.0;
                xr1 = 1.0 / xe;
                for (k = 1; k < 17; k++) {
                    r = -r * xr1;
                    su1 += a[k - 1] * r;
                }
                su2 = 1.0;
                r = 1.0;
                for (k = 1; k < 17; k++) {
                    r = r * xr1;
                    su2 += a[k - 1] * r;
                }
                apt = q0 - std::exp(-xe) * xp6 * su1;
                bpt = 2.0 * std::exp(xe) * xp6 * su2;
                su3 = 1.0;
                r = 1.0;
                xr2 = 1.0 / (xe * xe);
                for (k = 1; k < 9; k++) {
                    r = -r * xr2;
                    su3 += a[2 * k - 1] * r;
                }
                su4 = a[0] * xr1;
                r = xr1;
                for (k = 1; k < 8; k++) {
                    r = -r * xr2;
                    su4 += a[2 * k] * r;
                }
                su5 = su3 + su4;
                su6 = su3 - su4;
                ant = q1 - q2 * xp6 * (su5 * std::cos(xe) - su6 * std::sin(xe));
                bnt = q2 * xp6 * (su5 * std::sin(xe) + su6 * std::cos(xe));
            }
        }
    }

} // namespace detail

inline void airyb(double x, double *ai, double *bi, double *ad, double *bd) {

    // =======================================================
    // Purpose: Compute Airy functions and their derivatives
    // Input:   x  --- Argument of Airy function
    // Output:  AI --- Ai(x)
    //          BI --- Bi(x)
    //          AD --- Ai'(x)
    //          BD --- Bi'(x)
    // =======================================================

    int k, km, km2, kmax;
    double ck[52], dk[52];
    double xa, xq, xm, fx, r, gx, df, dg, sai, sad, sbi, sbd, xp1, xr2;
    double xf, rp, xar, xe, xr1, xcs, xss, ssa, sda, ssb, sdb;

    const double eps = 1.0e-15;
    const double pi = 3.141592653589793;
    const double c1 = 0.355028053887817;
    const double c2 = 0.258819403792807;
    const double sr3 = 1.732050807568877;

    km2 = 0;
    xa = fabs(x);
    xq = sqrt(xa);
    xm = 8.0;
    if (x > 0.0)
        xm = 5.0;

    if (x == 0.0) {
        *ai = c1;
        *bi = sr3 * c1;
        *ad = -c2;
        *bd = sr3 * c2;
        return;
    }

    if (xa <= xm) {
        fx = 1.0;
        r = 1.0;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k - 1.0) * x;
            fx += r;
            if (fabs(r) < fabs(fx) * eps)
                break;
        }

        gx = x;
        r = x;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k + 1.0) * x;
            gx += r;
            if (fabs(r) < fabs(gx) * eps)
                break;
        }

        *ai = c1 * fx - c2 * gx;
        *bi = sr3 * (c1 * fx + c2 * gx);

        df = 0.5 * x * x;
        r = df;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k + 2.0) * x;
            df += r;
            if (fabs(r) < fabs(df) * eps)
                break;
        }

        dg = 1.0;
        r = 1.0;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k - 2.0) * x;
            dg += r;
            if (fabs(r) < fabs(dg) * eps)
                break;
        }

        *ad = c1 * df - c2 * dg;
        *bd = sr3 * (c1 * df + c2 * dg);
    } else {
        km = (int) (24.5 - xa);
        if (xa < 6.0)
            km = 14;
        if (xa > 15.0)
            km = 10;

        if (x > 0.0) {
            kmax = km;
        } else {
            // Choose cutoffs so that the remainder term in asymptotic
            // expansion is epsilon size. The X<0 branch needs to be fast
            // in order to make AIRYZO efficient
            if (xa > 70.0)
                km = 3;
            if (xa > 500.0)
                km = 2;
            if (xa > 1000.0)
                km = 1;

            km2 = km;
            if (xa > 150.0)
                km2 = 1;
            if (xa > 3000.0)
                km2 = 0;
            kmax = 2 * km + 1;
        }
        xe = xa * xq / 1.5;
        xr1 = 1.0 / xe;
        xar = 1.0 / xq;
        xf = sqrt(xar);
        rp = 0.5641895835477563;
        r = 1.0;
        for (k = 1; k <= kmax; k++) {
            r = r * (6.0 * k - 1.0) / 216.0 * (6.0 * k - 3.0) / k * (6.0 * k - 5.0) / (2.0 * k - 1.0);
            ck[k - 1] = r;
            dk[k - 1] = -(6.0 * k + 1.0) / (6.0 * k - 1.0) * r;
        }

        if (x > 0.0) {
            sai = 1.0;
            sad = 1.0;
            r = 1.0;
            for (k = 1; k <= km; k++) {
                r *= -xr1;
                sai += ck[k - 1] * r;
                sad += dk[k - 1] * r;
            }
            sbi = 1.0;
            sbd = 1.0;
            r = 1.0;
            for (k = 1; k <= km; k++) {
                r *= xr1;
                sbi += ck[k - 1] * r;
                sbd += dk[k - 1] * r;
            }
            xp1 = exp(-xe);
            *ai = 0.5 * rp * xf * xp1 * sai;
            *bi = rp * xf / xp1 * sbi;
            *ad = -0.5 * rp / xf * xp1 * sad;
            *bd = rp / xf / xp1 * sbd;
        } else {
            xcs = cos(xe + pi / 4.0);
            xss = sin(xe + pi / 4.0);
            ssa = 1.0;
            sda = 1.0;
            r = 1.0;
            xr2 = 1.0 / (xe * xe);
            for (k = 1; k <= km; k++) {
                r *= -xr2;
                ssa += ck[2 * k - 1] * r;
                sda += dk[2 * k - 1] * r;
            }
            ssb = ck[0] * xr1;
            sdb = dk[0] * xr1;
            r = xr1;
            for (k = 1; k <= km2; k++) {
                r *= -xr2;
                ssb += ck[2 * k] * r;
                sdb += dk[2 * k] * r;
            }

            *ai = rp * xf * (xss * ssa - xcs * ssb);
            *bi = rp * xf * (xcs * ssa + xss * ssb);
            *ad = -rp / xf * (xcs * sda + xss * sdb);
            *bd = rp / xf * (xss * sda - xcs * sdb);
        }
    }
    return;
}

inline void airyzo(int nt, int kf, double *xa, double *xb, double *xc, double *xd) {

    // ========================================================
    // Purpose: Compute the first NT zeros of Airy functions
    //          Ai(x) and Ai'(x), a and a', and the associated
    //          values of Ai(a') and Ai'(a); and the first NT
    //          zeros of Airy functions Bi(x) and Bi'(x), b and
    //          b', and the associated values of Bi(b') and
    //          Bi'(b)
    // Input :  NT    --- Total number of zeros
    //          KF    --- Function code
    //                    KF=1 for Ai(x) and Ai'(x)
    //                    KF=2 for Bi(x) and Bi'(x)
    // Output:  XA(m) --- a, the m-th zero of Ai(x) or
    //                    b, the m-th zero of Bi(x)
    //          XB(m) --- a', the m-th zero of Ai'(x) or
    //                    b', the m-th zero of Bi'(x)
    //          XC(m) --- Ai(a') or Bi(b')
    //          XD(m) --- Ai'(a) or Bi'(b)
    //                    ( m --- Serial number of zeros )
    // Routine called: AIRYB for computing Airy functions and
    //                 their derivatives
    // =======================================================

    const double pi = 3.141592653589793;
    int i;
    double rt = 0.0, rt0, u = 0.0, u1 = 0.0, x, ai, bi, ad, bd, err;

    for (i = 1; i <= nt; ++i) {
        rt0 = 0.0;
        if (kf == 1) {
            u = 3.0 * pi * (4.0 * i - 1) / 8.0;
            u1 = 1 / (u * u);
        } else if (kf == 2) {
            if (i == 1) {
                rt0 = -1.17371;
            } else {
                u = 3.0 * pi * (4.0 * i - 3.0) / 8.0;
                u1 = 1 / (u * u);
            }
        }

        if (rt0 == 0) {
            // DLMF 9.9.18
            rt0 = -pow(u * u, 1.0 / 3.0) *
                  (1.0 +
                   u1 * (5.0 / 48.0 + u1 * (-5.0 / 36.0 + u1 * (77125.0 / 82944.0 + u1 * (-108056875.0 / 6967296.0)))));
        }

        while (1) {
            x = rt0;
            airyb(x, &ai, &bi, &ad, &bd);

            if (kf == 1) {
                rt = rt0 - ai / ad;
            } else if (kf == 2) {
                rt = rt0 - bi / bd;
            }

            err = fabs((rt - rt0) / rt);
            if (err <= 1.0e-12) {
                break;
            } else {
                rt0 = rt;
            }
        }

        xa[i - 1] = rt;
        if (err > 1.0e-14) {
            airyb(rt, &ai, &bi, &ad, &bd);
        }

        if (kf == 1) {
            xd[i - 1] = ad;
        } else if (kf == 2) {
            xd[i - 1] = bd;
        }
    }

    for (i = 1; i <= nt; ++i) {
        rt0 = 0.0;

        if (kf == 1) {
            if (i == 1) {
                rt0 = -1.01879;
            } else {
                u = 3.0 * pi * (4.0 * i - 3.0) / 8.0;
                u1 = 1 / (u * u);
            }
        } else if (kf == 2) {
            if (i == 1) {
                rt0 = -2.29444;
            } else {
                u = 3.0 * pi * (4.0 * i - 1.0) / 8.0;
                u1 = 1 / (u * u);
            }
        }

        if (rt0 == 0) {
            // DLMF 9.9.19
            rt0 = -pow(u * u, 1.0 / 3.0) *
                  (1.0 + u1 * (-7.0 / 48.0 +
                               u1 * (35.0 / 288.0 + u1 * (-181223.0 / 207360.0 + u1 * (18683371.0 / 1244160.0)))));
        }

        while (1) {
            x = rt0;
            airyb(x, &ai, &bi, &ad, &bd);

            if (kf == 1) {
                rt = rt0 - ad / (ai * x);
            } else if (kf == 2) {
                rt = rt0 - bd / (bi * x);
            }

            err = fabs((rt - rt0) / rt);
            if (err <= 1.0e-12) {
                break;
            } else {
                rt0 = rt;
            }
        }
        xb[i - 1] = rt;

        if (err > 1.0e-14) {
            airyb(rt, &ai, &bi, &ad, &bd);
        }

        if (kf == 1) {
            xc[i - 1] = ai;
        } else if (kf == 2) {
            xc[i - 1] = bi;
        }
    }
    return;
}

template <typename T>
void airy(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int ierr = 0;
    int kode = 1;
    int nz;

    ai = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), ai);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), bi);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), aip);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airy:", ierr_to_sferr(nz, ierr), bip);
}

template <typename T>
void airye(std::complex<T> z, std::complex<T> &ai, std::complex<T> &aip, std::complex<T> &bi, std::complex<T> &bip) {
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;

    ai = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), ai);

    nz = 0;
    bi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), bi);

    id = 1;
    aip = amos::airy(z, id, kode, &nz, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), aip);

    nz = 0;
    bip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), bip);
}

template <typename T>
void airye(T z, T &ai, T &aip, T &bi, T &bip) {
    int id = 0;
    int kode = 2; /* Exponential scaling */
    int nz, ierr;
    std::complex<T> cai, caip, cbi, cbip;

    cai.real(NAN);
    cai.imag(NAN);
    cbi.real(NAN);
    cbi.imag(NAN);
    caip.real(NAN);
    caip.imag(NAN);
    cbip.real(NAN);
    cbip.real(NAN);

    if (z < 0) {
        ai = NAN;
    } else {
        cai = amos::airy(z, id, kode, &nz, &ierr);
        set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cai);
        ai = std::real(cai);
    }

    nz = 0;
    cbi = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cbi);
    bi = std::real(cbi);

    id = 1;
    if (z < 0) {
        aip = NAN;
    } else {
        caip = amos::airy(z, id, kode, &nz, &ierr);
        set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), caip);
        aip = std::real(caip);
    }

    nz = 0;
    cbip = amos::biry(z, id, kode, &ierr);
    set_error_and_nan("airye:", ierr_to_sferr(nz, ierr), cbip);
    bip = std::real(cbip);
}

template <typename T>
void airy(T x, T &ai, T &aip, T &bi, T &bip) {
    /* For small arguments, use Cephes as it's slightly faster.
     * For large arguments, use AMOS as it's more accurate.
     */
    if (x < -10 || x > 10) {
        std::complex<T> zai, zaip, zbi, zbip;
        airy(std::complex(x), zai, zaip, zbi, zbip);
        ai = std::real(zai);
        aip = std::real(zaip);
        bi = std::real(zbi);
        bip = std::real(zbip);
    } else {
        cephes::airy(x, &ai, &aip, &bi, &bip);
    }
}

template <typename T>
void itairy(T x, T &apt, T &bpt, T &ant, T &bnt) {
    bool x_signbit = std::signbit(x);
    if (x_signbit) {
        x = -x;
    }

    detail::itairy(x, apt, bpt, ant, bnt);
    if (x_signbit) { /* negative limit -- switch signs and roles */
        T tmp = apt;
        apt = -ant;
        ant = -tmp;

        tmp = bpt;
        bpt = -bnt;
        bnt = -tmp;
    }
}

} // namespace xsf
