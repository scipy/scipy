/*
 *
 * This file accompanied with the header file specfun.h is a partial
 * C translation of the Fortran code by Zhang and Jin following
 * original description:
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *       COMPUTATION OF SPECIAL FUNCTIONS
 *
 *          Shanjie Zhang and Jianming Jin
 *
 *       Copyrighted but permission granted to use code in programs.
 *       Buy their book:
 *
 *          Shanjie Zhang, Jianming Jin,
 *          Computation of Special Functions,
 *          Wiley, 1996,
 *          ISBN: 0-471-11963-6,
 *          LC: QA351.C45.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *       Scipy changes:
 *       - Compiled into a single source file and changed REAL To DBLE throughout.
 *       - Changed according to ERRATA.
 *       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
 *       - Made functions return sf_error codes in ISFER variables instead
 *         of printing warnings. The codes are
 *         - SF_ERROR_OK        = 0: no error
 *         - SF_ERROR_SINGULAR  = 1: singularity encountered
 *         - SF_ERROR_UNDERFLOW = 2: floating point underflow
 *         - SF_ERROR_OVERFLOW  = 3: floating point overflow
 *         - SF_ERROR_SLOW      = 4: too many iterations required
 *         - SF_ERROR_LOSS      = 5: loss of precision
 *         - SF_ERROR_NO_RESULT = 6: no result obtained
 *         - SF_ERROR_DOMAIN    = 7: out of domain
 *         - SF_ERROR_ARG       = 8: invalid input parameter
 *         - SF_ERROR_OTHER     = 9: unclassified error
 *       - Improved initial guesses for roots in JYZO.
 *
 *
 */

/*
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "specfun.h"

#ifndef CMPLX
#define CMPLX(x, y) ((double complex)((double)(x) + I * (double)(y)))
#endif /* CMPLX */

static void specfun_airyb(double, double *, double *, double *, double *);
static void specfun_bjndd(double, int, double *, double *, double *);
static void specfun_cbk(int, int, double, double, double, double *, double *);
static void specfun_cerf(double complex, double complex *, double complex *);
static double specfun_chgubi(double, double, double, int *);
static double specfun_chguit(double, double, double, int *);
static double specfun_chgul(double, double, double, int *);
static double specfun_chgus(double, double, double, int *);
static double complex specfun_cpdla(int, double complex);
static double complex specfun_cpdsa(int, double complex);
static double specfun_cv0(double, double, double);
static double specfun_cvf(int, int, double, double, int);
static double specfun_cvql(int, int, double);
static double specfun_cvqm(int, double);
static double specfun_dvla(double, double);
static double specfun_dvsa(double, double);
static double specfun_gaih(double);
static double specfun_gam0(double);
static double specfun_gamma2(double);
static void specfun_gmn(int, int, double, double, double *, double *, double *);
static void specfun_jynb(int, double, int *, double *, double *, double *, double *);
static void specfun_jynbh(int, int, double, int *, double *, double *);
static void specfun_jyndd(int, double, double *, double *, double *, double *, double *, double *);
static void specfun_kmn(int, int, double, double, int, double *, double *, double *, double *);
static double specfun_lpmv0(double, int, double);
static void specfun_lpmns(int, int, double, double *, double *);
static void specfun_lqmns(int, int, double, double *, double *);
static int specfun_msta1(double, int);
static int specfun_msta2(double, int, int);
static double specfun_psi_spec(double);
static void specfun_qstar(int, int, double, double, double *, double *, double *);
static double specfun_refine(int, int, double, double);
static void specfun_rmn1(int, int, double, double, int, double *, double *, double *);
static void specfun_rmn2l(int, int, double, double, int, double *, double *, double *, int *);
static void specfun_rmn2so(int, int, double, double, double, int, double *, double *, double *);
static void specfun_rmn2sp(int, int, double, double, double, int, double *, double *, double *);
static void specfun_sckb(int, int, double, double *, double *);
static void specfun_sphj(double, int, int *, double *, double *);
static void specfun_sphy(double, int, int *, double *, double *);
static double specfun_vvla(double, double);
static double specfun_vvsa(double, double);


void specfun_airyb(double x, double* ai, double* bi, double* ad, double* bd) {

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
    if (x > 0.0) xm = 5.0;

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
            if (fabs(r) < fabs(fx) * eps) break;
        }

        gx = x;
        r = x;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k + 1.0) * x;
            gx += r;
            if (fabs(r) < fabs(gx) * eps) break;
        }

        *ai = c1 * fx - c2 * gx;
        *bi = sr3 * (c1 * fx + c2 * gx);

        df = 0.5 * x * x;
        r = df;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k + 2.0) * x;
            df += r;
            if (fabs(r) < fabs(df) * eps) break;
        }

        dg = 1.0;
        r = 1.0;
        for (k = 1; k <= 40; k++) {
            r = r * x / (3.0 * k) * x / (3.0 * k - 2.0) * x;
            dg += r;
            if (fabs(r) < fabs(dg) * eps) break;
        }

        *ad = c1 * df - c2 * dg;
        *bd = sr3 * (c1 * df + c2 * dg);
    } else {
        km = (int)(24.5 - xa);
        if (xa < 6.0) km = 14;
        if (xa > 15.0) km = 10;

        if (x > 0.0) {
            kmax = km;
        } else {
            // Choose cutoffs so that the remainder term in asymptotic
            // expansion is epsilon size. The X<0 branch needs to be fast
            // in order to make AIRYZO efficient
            if (xa > 70.0) km = 3;
            if (xa > 500.0) km = 2;
            if (xa > 1000.0) km = 1;

            km2 = km;
            if (xa > 150.0) km2 = 1;
            if (xa > 3000.0) km2 = 0;
            kmax = 2 * km + 1;
        }
        xe = xa * xq / 1.5;
        xr1 = 1.0 / xe;
        xar = 1.0 / xq;
        xf = sqrt(xar);
        rp = 0.5641895835477563;
        r = 1.0;
        for (k = 1; k <= kmax; k++) {
            r = r * (6.0 * k - 1.0) / 216.0 * (6.0 * k - 3.0)
                / k * (6.0 * k - 5.0) / (2.0 * k - 1.0);
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


void specfun_airyzo(int nt, int kf, double *xa, double *xb, double *xc, double *xd) {

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
            rt0 = -pow(u * u, 1.0 / 3.0) * (
                    1.0
                    + u1 * (5.0 / 48.0
                    + u1 * (-5.0 / 36.0
                    + u1 * (77125.0 / 82944.0
                    + u1 * (-108056875.0 / 6967296.0)))));
        }

        while (1) {
            x = rt0;
            specfun_airyb(x, &ai, &bi, &ad, &bd);

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
            specfun_airyb(rt, &ai, &bi, &ad, &bd);
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
            rt0 = -pow(u * u, 1.0 / 3.0) * (
                    1.0
                    + u1 * (-7.0 / 48.0
                    + u1 * (35.0 / 288.0
                    + u1 * (-181223.0 / 207360.0
                    + u1 * (18683371.0 / 1244160.0)))));
        }

        while (1) {
            x = rt0;
            specfun_airyb(x, &ai, &bi, &ad, &bd);

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
            specfun_airyb(rt, &ai, &bi, &ad, &bd);
        }

        if (kf == 1) {
            xc[i - 1] = ai;
        } else if (kf == 2) {
            xc[i - 1] = bi;
        }
    }
    return;
}


void specfun_aswfa(double x, int m, int n, double c, int kd, double cv, double *s1f, double *s1d) {

    // ===========================================================
    // Purpose: Compute the prolate and oblate spheroidal angular
    //          functions of the first kind and their derivatives
    // Input :  m  --- Mode parameter,  m = 0,1,2,...
    //          n  --- Mode parameter,  n = m,m+1,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument of angular function, |x| < 1.0
    //          KD --- Function code
    //                 KD=1 for prolate;  KD=-1 for oblate
    //          cv --- Characteristic value
    // Output:  S1F --- Angular function of the first kind
    //          S1D --- Derivative of the angular function of
    //                  the first kind
    // Routine called:
    //          SCKB for computing expansion coefficients ck
    // ===========================================================

    int ip, k, nm, nm2;
    double a0, d0, d1, r, su1, su2, x0, x1;
    double *ck = calloc(200, sizeof(double));
    double *df = calloc(200, sizeof(double));
    const double eps = 1e-14;
    x0 = x;
    x = fabs(x);
    ip = ((n-m) % 2 == 0 ? 0 : 1);
    nm = 40 + (int)((n-m)/2 + c);
    nm2 = nm/2 - 2;
    specfun_sdmn(m, n, c, cv, kd, df);
    specfun_sckb(m, n, c, df, ck);
    x1 = 1.0 - x*x;
    if ((m == 0) && (x1 == 0.0)) {
        a0 = 1.0;
    } else {
        a0 = pow(x1, 0.5*m);
    }
    su1 = ck[0];
    for (k = 1; k <= nm2; k++) {
        r = ck[k]*pow(x1, k);
        su1 += r;
        if ((k >= 10) && (fabs(r/su1) < eps)) { break; }
    }
    *s1f = a0*pow(x, ip)*su1;
    if (x == 1.0) {
        if (m == 0) {
            *s1d = ip*ck[0] - 2.0*ck[1];
        } else if (m == 1) {
            *s1d = -1e100;
        } else if (m == 2) {
            *s1d = -2.0*ck[0];
        } else if (m >= 3) {
            *s1d = 0.0;
        }
    } else {
        d0 = ip - m/x1*pow(x, ip+1.0);
        d1 = -2.0*a0*pow(x, ip+1.0);
        su2 = ck[1];
        for (k = 2; k <= nm2; k++) {
            r = k*ck[k]*pow(x1, (k-1.0));
            su2 += r;
            if ((k >= 10) && (fabs(r/su2) < eps)) { break; }
        }
        *s1d = d0*a0*su1 + d1*su2;
    }
    if ((x0 < 0.0) && (ip == 0)) { *s1d = -*s1d; }
    if ((x0 < 0.0) && (ip == 1)) { *s1f = -*s1f; }
    x = x0;
    free(ck); free(df);
    return;
}


void specfun_bernob(int n, double *bn) {

    // ======================================
    // Purpose: Compute Bernoulli number Bn
    // Input :  n >= 3 --- Serial number
    // Output:  BN(n)  --- Bn
    // ======================================

    int k, m;
    double r1, r2, s;
    const double tpi = 6.283185307179586;

    bn[0] = 1.0;
    bn[1] = -0.5;
    bn[2] = 1.0 / 6.0;
    r1 = pow(2.0 / tpi, 2);
    for ( m = 4; m < (n+1); m += 2) {
        r1 = -r1 * (m-1)*m/(tpi*tpi);
        r2 = 1.0;
        for (k = 2; k < 10001; k++) {
            s = pow(1.0/k, m);
            r2 += s;
            if (s < 1e-15) { break; }
        }
        bn[m] = r1*r2;
    }
    return;
}


void specfun_bjndd(double x, int n, double *bj, double *dj, double *fj) {

    // =====================================================
    // Purpose: Compute Bessel functions Jn(x) and their
    //          first and second derivatives ( 0 <= n <= 100)
    // Input:   x ---  Argument of Jn(x)  ( x ≥ 0 )
    //          n ---  Order of Jn(x)
    // Output:  BJ(n+1) ---  Jn(x)
    //          DJ(n+1) ---  Jn'(x)
    //          FJ(n+1) ---  Jn"(x)
    // =====================================================

    int k, m, mt;
    double bs = 0.0, f = 0.0, f0 = 0.0, f1 = 1e-35;

    for (m = 1; m < 901; m++) {
        mt = (int)(0.5*log10(6.28*m)-m*log10(1.36*fabs(x)/m));
        if (mt > 20) { break; }
    }
    if (m == 901) { m -= 1; }
    for (k = m; k > -1; k--) {
        f = 2.0*(k+1.0)*f1/x - f0;
        if (k <= n) { bj[k] = f; }
        if (k % 2 == 0) { bs += 2.0*f; }
        f0 = f1;
        f1 = f;
    }
    for (k = 0; k < (n+1); k++) {
        bj[k] /= (bs - f);
    }
    dj[0] = -bj[1];
    fj[0] = -bj[0] - dj[0]/x;
    for (k = 1; k < (n+1); k++) {
        dj[k] = bj[k-1] - k*bj[k]/x;
        fj[k] = (k*k/(x*x)-1.0)*bj[k] - dj[k]/x;
    }
    return;
}


void specfun_cbk(int m, int n, double c, double cv, double qt, double *ck, double *bk) {
    const double eps = 1.0e-14;

    int i, i1, ip, j, k, n2, nm;
    double r1, s1, sw, t;

    ip = ((n - m) % 2 == 0 ? 0 : 1);
    nm = 25 + (int)(0.5 * (n - m) + c);
    double *u = calloc(200, sizeof(double));
    double *v = calloc(200, sizeof(double));
    double *w = calloc(200, sizeof(double));

    u[0] = 0.0;
    n2 = nm - 2;

    for (j = 1; j < n2 + 1; j++) {
        u[j] = c * c;
    }
    for (j = 1; j < n2 + 1; j++) {
        v[j - 1] = (2.0 * j - 1.0 - ip) * (2.0 * (j - m) - ip) + m * (m - 1.0) - cv;
    }
    for (j = 1; j < nm; j++)
        w[j - 1] = (2.0 * j - ip) * (2.0 * j + 1.0 - ip);

    if (ip == 0) {
        sw = 0.0;
        for (k = 0; k < n2; k++) {
            s1 = 0.0;
            i1 = k - m + 1;

            for (i = i1; i < nm + 1; i++) {
                if (i < 0) { continue; }
                r1 = 1.0;
                for (j = 1; j <= k; j++) {
                    r1 =  r1*(i + m - j) / (1.0*j);
                }
                s1 += ck[i] * (2.0 * i + m) * r1;
                if (fabs(s1 - sw) < fabs(s1) * eps) { break; }
                sw = s1;
            }
            bk[k] = qt * s1;
        }
    } else if (ip == 1) {
        sw = 0.0;
        for (k = 0; k < n2; k++) {
            s1 = 0.0;
            i1 = k - m + 1;

            for (int i = i1; i < nm + 1; i++) {
                if (i < 0) { continue; }
                r1 = 1.0;
                for (j = 1; j <= k; j++) {
                    r1 = r1* (i + m - j) / (1.0*j);
                }
                if (i > 0) {
                    s1 += ck[i - 1] * (2.0 * i + m - 1) * r1;
                }
                s1 -= ck[i] * (2.0 * i + m) * r1;

                if (fabs(s1 - sw) < fabs(s1) * eps) { break; }
                sw = s1;
            }
            bk[k] = qt * s1;
        }
    }

    w[0] /= v[0];
    bk[0] /= v[0];

    for (k = 2; k <= n2; k++) {
        t = v[k - 1] - w[k - 2] * u[k - 1];
        w[k - 1] /= t;
        bk[k - 1] = (bk[k - 1] - bk[k - 2] * u[k - 1]) / t;
    }

    for (k = n2 - 1; k >= 1; k--) {
        bk[k - 1] -= w[k - 1] * bk[k];
    }
    free(u); free(v); free(w);
    return;
}


void specfun_cerf(double complex z, double complex *cer, double complex *cder) {

    // ==========================================================
    // Purpose: Compute complex Error function erf(z) & erf'(z)
    // Input:   z   --- Complex argument of erf(z)
    //          x   --- Real part of z
    //          y   --- Imaginary part of z
    // Output:  CER --- erf(z)
    //          CDER --- erf'(z)
    // ==========================================================

    int k;
    double c0, cs, er0, er, er1, ei1, er2, ei2, err, eri, r, ss, w, w1, w2;
    const double eps = 1.0e-12;
    const double pi = 3.141592653589793;

    double x = creal(z);
    double y = cimag(z);
    double x2 = x * x;

    if (x <= 3.5) {
        er = 1.0;
        r = 1.0;
        w = 0.0;

        for (k = 1; k <= 100; k++) {
            r = r * x2 / (k + 0.5);
            er += r;
            if (fabs(er - w) <= eps * fabs(er))
                break;
            w = er;
        }

        c0 = 2.0 / sqrt(pi) * x * exp(-x2);
        er0 = c0 * er;
        *cer = CMPLX(er0, 0.0);
    } else {
        er = 1.0;
        r = 1.0;

        for (k = 1; k <= 12; k++) {
            r = -r * (k - 0.5) / x2;
            er += r;
        }

        c0 = exp(-x2) / (x * sqrt(pi));
        er0 = 1.0 - c0 * er;
        *cer = CMPLX(er0, 0.0);
    }

    if (y == 0.0) {
        err = creal(*cer);
        eri = 0.0;
        *cer = CMPLX(err, eri);
    } else {
        cs = cos(2.0 * x * y);
        ss = sin(2.0 * x * y);
        er1 = exp(-x2) * (1.0 - cs) / (2.0 * pi * x);
        ei1 = exp(-x2) * ss / (2.0 * pi * x);
        er2 = 0.0;
        w1 = 0.0;

        for (int n = 1; n <= 100; n++) {
            er2 += exp(-0.25 * n * n) / (n * n + 4.0 * x2) * (2.0 * x - 2.0 * x * cosh(n * y) * cs + n * sinh(n * y) * ss);
            if (fabs((er2 - w1) / er2) < eps)
                break;
            w1 = er2;
        }

        c0 = 2.0 * exp(-x2) / pi;
        err = creal(*cer) + er1 + c0 * er2;
        ei2 = 0.0;
        w2 = 0.0;

        for (int n = 1; n <= 100; n++) {
            ei2 += exp(-0.25 * n * n) / (n * n + 4.0 * x2) * (2.0 * x * cosh(n * y) * ss + n * sinh(n * y) * cs);
            if (fabs((ei2 - w2) / ei2) < eps)
                break;
            w2 = ei2;
        }
        *cer = CMPLX(err, ei1 + c0 * ei2);
    }
    *cder = 2.0 / sqrt(pi) * cexp(-z*z);

}


double complex specfun_cerror(double complex z) {

    // ====================================================
    // Purpose: Compute error function erf(z) for a complex
    //          argument (z=x+iy)
    // Input :  z   --- Complex argument
    // Output:  CER --- erf(z)
    // ====================================================

    int k;
    double complex cer, cl, cr, cs, z1;
    double complex c0 = cexp(-z*z);
    const double sqpi = 1.7724538509055160273;
    z1 = z;
    if (creal(z) < 0.0) { z1 = -z; }
    // Cutoff radius R = 4.36; determined by balancing rounding error
    // and asymptotic expansion error, see below.
    //
    // The resulting maximum global accuracy expected is around 1e-8
    //
    if (cabs(z) <= 4.36) {
        // Rounding error in the Taylor expansion is roughly
        // ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2))
        cs = z1;
        cr = z1;
        for (k = 1; k < 121; k++) {
            cr = cr*(z1*z1) / (k+0.5);
            cs += cr;
            if (cabs(cr/cs) < 1e-15) { break; }
        }
        cer = 2.0*c0*cs/sqpi;
    } else {
        cl = 1.0 / z1;
        cr = cl;
        // Asymptotic series; maximum K must be at most ~ R^2.
        //
        // The maximum accuracy obtainable from this expansion is roughly
        //
        // ~ Gamma(2R**2 + 2) / (
        //          (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2))
        for (k = 1; k < 21; k++) {
            cr = -cr*(k-0.5) / (z1*z1);
            cl += cr;
            if (cabs(cr/cl) < 1e-15) { break; }
        }
        cer = 1.0 - c0*cl/sqpi;
    }
    if (creal(z) < 0.0) { cer = -cer; }
    return cer;
}


void specfun_cerzo(int nt, double complex *zo) {

    // ===============================================================
    // Purpose : Evaluate the complex zeros of error function erf(z)
    //           using the modified Newton's iteration method
    // Input :   NT --- Total number of zeros
    // Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
    // Routine called: CERF for computing erf(z) and erf'(z)
    // ===============================================================

    int i, j, nr, it = 0;
    double pu, pv, px, py, w0;
    double complex z, zf, zd, zp, zw, zq, zfd, zgd;
    double w = 0.0;
    const double pi = 3.141592653589793;

    for (nr = 1; nr <= nt; nr++) {
        pu = sqrt(pi * (4.0 * nr - 0.5));
        pv = pi * sqrt(2.0 * nr - 0.25);
        px = 0.5 * pu - 0.5 * log(pv) / pu;
        py = 0.5 * pu + 0.5 * log(pv) / pu;
        z = px + I * py;
        it = 0;

        do {
            it++;
            specfun_cerf(z, &zf, &zd);
            zp = 1.0;

            for (i = 1; i < nr; i++) {
                zp *= (z - zo[i - 1]);
            }
            zfd = zf / zp;
            zq = 0.0;
            for (i = 1; i < nr; i++) {
                zw = 1.0;
                for (j = 1; j < nr; j++) {
                    if (j == i) continue;
                    zw *= (z - zo[j - 1]);
                }
                zq += zw;
            }
            zgd = (zd - zq * zfd) / zp;
            z -= zfd / zgd;
            w0 = w;
            w = cabs(z);
        } while ((it <= 50) && (fabs((w - w0) / w) > 1.0e-11));
        zo[nr - 1] = z;
    }
    return;
}


void specfun_cfc(double complex z, double complex *zf, double complex *zd) {

    // =========================================================
    // Purpose: Compute complex Fresnel integral C(z) and C'(z)
    // Input :  z --- Argument of C(z)
    // Output:  ZF --- C(z)
    //          ZD --- C'(z)
    // =========================================================

    int k, m;
    double wa0, wa;
    double complex c, cr, cf, cf0, cf1, cg, d;
    const double eps = 1.0e-14;
    const double pi = 3.141592653589793;

    double w0 = cabs(z);
    double complex zp = 0.5 * pi * z * z;
    double complex zp2 = zp * zp;
    double complex z0 = 0.0;

    if (z == z0) {
        c = z0;
    } else if (w0 <= 2.5) {
        cr = z;
        c = cr;
        wa0 = 0.0;
        for (k = 1; k <= 80; k++) {
            cr = -0.5*cr*(4.0*k - 3.0)/k/(2.0*k - 1.0)/(4.0*k + 1.0)*zp2;
            c += cr;
            wa = cabs(c);
            if ((fabs((wa - wa0) / wa) < eps) && (k > 10)) {
                *zf = c;
                *zd = ccos(0.5*pi*z*z);
                return;
            }
            wa0 = wa;
        }
    } else if ((w0 > 2.5) && (w0 < 4.5)) {
        m = 85;
        c = z0;
        cf1 = z0;
        cf0 = 1.0e-100;
        for (k = m; k >= 0; k--) {
            cf = (2.0*k + 3.0)*cf0/zp - cf1;
            if (k % 2 == 0) { c += cf; }
            cf1 = cf0;
            cf0 = cf;
        }
        c *= 2.0/(pi*z)*csin(zp)/cf;
    } else {
        // See comment at CFS(), use C(z) = iC(-iz)
        if ((cimag(z) > -creal(z)) && (cimag(z) <= creal(z))) {
            // right quadrant
            d = 0.5;
        } else if ((cimag(z) > creal(z)) && (cimag(z) >= -creal(z))) {
            // upper quadrant
            d = 0.5*I;
        } else if ((cimag(z) < -creal(z)) && (cimag(z) >= creal(z))) {
            // left quadrant
            d = -0.5;
        } else {
            d = -0.5*I;
        }
        cr = 1.0;
        cf = 1.0;
        for (k = 1; k <= 20; k++) {
            cr = -0.25*cr*(4.0*k - 1.0)*(4.0*k - 3.0)/zp2;
            cf += cr;
        }
        cr = 1.0/(pi*z*z);
        cg = cr;
        for (k = 1; k <= 12; k++) {
            cr = -0.25*cr*(4.0*k + 1.0)*(4.0*k - 1.0)/zp2;
            cg += cr;
        }
        c = d + (cf*csin(zp) - cg*ccos(zp))/(pi*z);
    }
    *zf = c;
    *zd = ccos(0.5*pi*z*z);
    return;
}


void specfun_cfs(double complex z, double complex *zf, double complex *zd) {

    // =========================================================
    // Purpose: Compute complex Fresnel Integral S(z) and S'(z)
    // Input :  z  --- Argument of S(z)
    // Output:  ZF --- S(z)
    //          ZD --- S'(z)
    // =========================================================

    int k, m;
    double wb0, wb;
    double complex s, cr, cf, cf0, cf1, cg, d;
    const double eps = 1.0e-14;
    const double pi = 3.141592653589793;

    double w0 = cabs(z);
    double complex zp = 0.5 * pi * z * z;
    double complex zp2 = zp * zp;
    double complex z0 = 0.0;

    if (z == z0) {
        s = z0;
    } else if (w0 <= 2.5) {
        s = z * zp / 3.0;
        cr = s;
        wb0 = 0.0;
        for (k = 1; k <= 80; k++) {
            cr = -0.5 * cr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * zp2;
            s += cr;
            wb = cabs(s);
            if ((fabs(wb - wb0) < eps) && (k > 10)) {
                    *zf = s;
                    *zd = csin(0.5*pi*z*z);
                    return;
            }
            wb0 = wb;
        }
    } else if ((w0 > 2.5) && (w0 < 4.5)) {
        m = 85;
        s = z0;
        cf1 = z0;
        cf0 = 1.0e-100;
        for (k = m; k >= 0; k--) {
            cf = (2.0*k + 3.0)*cf0/zp - cf1;
            if (k % 2 == 1) { s += cf; }
            cf1 = cf0;
            cf0 = cf;
        }
        s = 2.0/(pi * z)*csin(zp)/cf*s;
    } else {
        // Auxiliary functions f(z) and g(z) can be computed using an
        // asymptotic expansion in the right quadrant |arg(z)| <= pi/4, not pi/2
        // as sometimes suggested. Use the symmetry S(z) = -iS(-iz).
        // Interestingly, most of the expansion code is the same across
        // the quadrants. (The forth power in Z is the equalizer here.)
        // Only one constant has to be adapted.
        if ((cimag(z) > -creal(z)) && (cimag(z) <= creal(z))) {
            // right quadrant
            d = 0.5;
        } else if ((cimag(z) > creal(z)) && (cimag(z) >= -creal(z))) {
            // upper quadrant
            d = -0.5*I;
        } else if ((cimag(z) < -creal(z)) && (cimag(z) >= creal(z))) {
            // left quadrant
            d = -0.5;
        } else {
            d = 0.5*I;
        }
        cr = 1.0;
        cf = 1.0;
        for (k = 1; k <= 20; k++) {
            cr = -0.25*cr*(4.0*k - 1.0)*(4.0*k - 3.0)/zp2;
            cf += cr;
        }
        cr = 1.0;
        cg = 1.0;
        for (k = 1; k <= 12; k++) {
            cr = -0.25*cr*(4.0*k + 1.0)*(4.0*k - 1.0)/zp2;
            cg += cr;
        }
        cg = cg/(pi*z*z);
        s = d - (cf*ccos(zp) + cg*csin(zp))/(pi*z);
    }
    *zf = s;
    *zd = csin(0.5*pi*z*z);
    return;
}


double complex specfun_cchg(double a, double b, double complex z) {

    // ===================================================
    // Purpose: Compute confluent hypergeometric function
    //          M(a,b,z) with real parameters a, b and a
    //          complex argument z
    // Input :  a --- Parameter
    //          b --- Parameter
    //          z --- Complex argument
    // Output:  CHG --- M(a,b,z)
    // Routine called: CGAMA for computing complex ln[Г(x)]
    // ===================================================

    int i, j, k, la, m, n, nl, ns;
    double a0, a1, phi, x0, x, y;
    double complex cfac, cg1, cg2, cg3, chg, chg1, chg2, chw, cr, cr1, cr2, cs1,\
                   cs2, crg, cy0, cy1, z0;
    const double pi = 3.141592653589793;
    const double complex ci = CMPLX(0.0, 1.0);
    a0 = a;
    a1 = a;
    z0 = z;
    cy0 = 0.0;
    cy1 = 0.0;
    if ((b == 0.0) || (b == -(int)fabs(b))) { return 1e300; }
    if ((a == 0.0) || (z == 0.0)) { return 1.0; }
    if (a == -1.0) { return 1.0 - z/b; }
    if (a == b) { return cexp(z); }
    if (a - b == 1.0) { return (1.0 + z/b)*cexp(z); }
    if ((a == 1.0) && (b == 2.0)) { return (cexp(z)-1.0) / z; }
    if ((a == (int)a) && (a < 0.0)) {
        m = (int)(-a);
        cr = 1.0;
        chg = 1.0;
        for (k = 1; k < (m+1); k++) {
            cr = cr * (a+k-1.0)/k/(b+k-1.0)*z;
            chg += cr;
        }
    } else {
        x0 = creal(z);
        if (x0 < 0.0) {
            a = b-a;
            a0 = a;
            z = -z;
        }
        nl = 0;
        la = 0;
        if (a >= 2.0) {
            nl = 1;
            la = (int)a;
            a -= la + 1;
        }
        ns = 0;
        for (n = 0; n < (nl+1); n++) {
            if (a0 >= 2.0) { a += 1.0; }
            if ((cabs(z) < 20.0+fabs(b)) || (a < 0.0)) {
                chg = 1.0;
                chw = 0.0;
                crg = 1.0;
                for (j = 1; j < 501; j++) {
                    crg = crg * (a+j-1.0)/(j*(b+j-1.0))*z;
                    if (cabs((chg-chw)/chg) < 1e-15) { break; }
                    chw = chg;
                }
            } else {
                y = 0.0;
                cg1 = specfun_cgama(a, 0);
                cg2 = specfun_cgama(b, 0);
                cg3 = specfun_cgama(b-a, 0);
                cs1 = 1.0;
                cs2 = 1.0;
                cr1 = 1.0;
                cr2 = 1.0;
                for (i = 1; i <= 8; i++) {
                    cr1 = -cr1 * (a+i-1.0)*(a-b+i)/(z*i);
                    cr2 = cr2 * (b-a+i-1.0)*(i-a)/(z*i);
                    cs1 += cr1;
                    cs2 += cr2;
                }
                x = creal(z);
                y = cimag(z);
                if ((x == 0.0) && (y >= 0.0)) {
                    phi = 0.5*pi;
                } else if ((x == 0.0) && (y <= 0.0)) {
                    phi = -0.5*pi;
                } else {
                    phi = atan(y/x);
                }
                if ((phi > -0.5*pi) && (phi < 1.5*pi)) { ns = 1; }
                if ((phi > -1.5*pi) && (phi <= -0.5*pi)) { ns = -1; }
                cfac = cexp(ns*ci*pi*a);
                if (y == 0.0) { cfac = cos(pi*a); }
                chg1 = cexp(cg2-cg3)*cpow(z, -a)*cfac*cs1;
                chg2 = cexp(cg2-cg1+z)*cpow(z, a-b)*cs2;
                chg = chg1 + chg2;
            }
            if (n == 0) { cy0 = chg; }
            if (n == 1) { cy1 = chg; }
        }
        if (a0 >= 2.0) {
            for (i = 1; i < la; i++) {
                chg = ((2.0*a-b+z)*cy1 + (b-a)*cy0)/a;
                cy0 = cy1;
                cy1 = chg;
                a += 1.0;
            }
        }
        if (x0 < 0.0) { chg *= cexp(-z); }
    }
    a = a1;
    z = z0;
    return chg;
}


double complex specfun_cgama(double complex z, int kf) {

    // =========================================================
    // Purpose: Compute the gamma function Г(z) or ln[Г(z)]
    //          for a complex argument
    // Input :  z  --- Complex argument
    //          kf --- Function code
    //                 kf=0 for ln[Г(z)]
    //                 kf=1 for Г(z)
    // Output:  g  --- ln[Г(z)] or Г(z)
    // ========================================================

    double complex g, z1;
    double az0, az1, gi, gi1, gr, gr1, t, th, th1, th2, sr, si, x0, xx, yy;
    int j, k, na;
    const double pi = 3.141592653589793;
    static const double a[10] = {
        8.333333333333333e-02, -2.777777777777778e-03,
        7.936507936507937e-04, -5.952380952380952e-04,
        8.417508417508418e-04, -1.917526917526918e-03,
        6.410256410256410e-03, -2.955065359477124e-02,
        1.796443723688307e-01, -1.392432216905900e+00
    };
    xx = creal(z);
    yy = cimag(z);
    if ((yy == 0.0) && (xx <= 0.0) && (xx == (int)xx)) {
        return 1e300;
    } else if (xx < 0.0) {
        z1 = z;
        z = -z;
        xx = -xx;
        yy = -yy;
    } else {
        z1 = CMPLX(xx, 0.0);
    }
    x0 = xx;
    na = 0;
    if (xx <= 7.0) {
        na = (int)(7 - xx);
        x0 = xx + na;
    }
    az0 = cabs(CMPLX(x0, yy));
    th = atan(yy / x0);
    gr = (x0 - 0.5)*log(az0) - th*yy - x0 + 0.5*log(2.0*pi);
    gi = th*(x0 - 0.5) + yy*log(az0) - yy;
    for (k = 1; k < 11; k++) {
        t = pow(az0, 1-2*k);
        gr += a[k - 1]*t*cos((2.0*k - 1.0)*th);
        gi += -a[k - 1]*t*sin((2.0*k - 1.0)*th);
    }
    if (xx <= 7.0) {
        gr1 = 0.0;
        gi1 = 0.0;
        for (j = 0; j < na; j++) {
            gr1 += 0.5*log(pow(xx + j, 2) + yy*yy);
            gi1 += atan(yy/(xx + j));
        }
        gr -= gr1;
        gi -= gi1;
    }
    if (creal(z1) < 0.0) {
        az0 = cabs(z);
        th1 = atan(yy/xx);
        sr = -sin(pi*xx)*cosh(pi*yy);
        si = -cos(pi*xx)*sinh(pi*yy);
        az1 = cabs(CMPLX(sr, si));
        th2 = atan(si/sr);
        if (sr < 0.0) {
            th2 += pi;
        }
        gr = log(pi/(az0*az1)) - gr;
        gi = - th1 - th2 - gi;
        z = z1;
    }
    if (kf == 1) {
        g = exp(gr)*CMPLX(cos(gi), sin(gi));
    } else {
        g = CMPLX(gr, gi);
    }
    return g;
}


double specfun_chgm(double x, double a, double b) {

    // ===================================================
    // Purpose: Compute confluent hypergeometric function
    //          M(a,b,x)
    // Input  : a  --- Parameter
    //          b  --- Parameter ( b <> 0,-1,-2,... )
    //          x  --- Argument
    // Output:  HG --- M(a,b,x)
    // Routine called: CGAMA for computing complex ln[Г(x)]
    // ===================================================

    int i, j, la, n, nl;
    double a0 = a, a1 = a, x0 = x, y0, y1, hg1, hg2, r1, r2, rg, xg, sum1, sum2;
    double complex cta, ctb, ctba;
    const double pi = 3.141592653589793;
    double hg = 0.0;

    //  DLMF 13.2.39
    if (x < 0.0) {
        a = b - a;
        a0 = a;
        x = fabs(x);
    }
    nl = 0;
    la = 0;
    if (a >= 2.0) {
        // preparing terms for DLMF 13.3.1
        nl = 1;
        la = (int)a;
        a -= la+1;
    }
    y0 = 0.0;
    y1 = 0.0;
    for (n = 0; n < (nl + 1); n++) {
        if (a0 >= 2.0) { a += 1.0; }
        if ((x <= 30.0 + fabs(b)) || (a < 0.0)) {
            hg = 1.0;
            rg = 1.0;
            for (j = 1; j < 501; j++) {
                rg = rg * (a + j - 1.0) / (j*(b + j - 1.0))*x;
                hg += rg;
                if ((hg != 0.0) && (fabs(rg/hg) < 1e-15)) {
                    // DLMF 13.2.39 (cf. above)
                    if (x0 < 0.0) { hg *= exp(x0); }
                    break;
                }
            }
        } else {
            // DLMF 13.7.2 & 13.2.4, SUM2 corresponds to first sum
            cta = specfun_cgama(a, 0);
            ctb = specfun_cgama(b, 0);
            xg = b-a;
            ctba = specfun_cgama(xg, 0);
            sum1 = 1.0;
            sum2 = 1.0;
            r1 = 1.0;
            r2 = 1.0;
            for (i = 1; i < 9; i++) {
                r1 = -r1*(a+i-1.0)*(a-b+i)/(x*i);
                r2 = -r2*(b-a+i-1.0)*(a-i)/(x*i);
                sum1 += r1;
                sum2 += r2;
            }
            if (x0 >= 0.0) {
                hg1 = creal(cexp(ctb-ctba))*pow(x, -a)*cos(pi*a)*sum1;
                hg2 = creal(cexp(ctb-cta+x))*pow(x, a-b)*sum2;
            } else {
            // DLMF 13.2.39 (cf. above)
                 hg1 = creal(cexp(ctb-ctba+x0))*pow(x, -a)*cos(pi*a)*sum1;
                 hg2 = creal(cexp(ctb-cta))*pow(x, a-b)*sum2;
            }
            hg = hg1 + hg2;
        }
        /* 25 */
        if (n == 0) { y0 = hg; }
        if (n == 1) { y1 = hg; }
    }
    if (a0 >= 2.0) {
        // DLMF 13.3.1
        for (i = 1; i < la; i++) {
            hg = ((2.0*a - b + x)*y1 + (b - a)*y0) / a;
            y0 = y1;
            y1 = hg;
            a += 1.0;
        }
    }
    a = a1;
    x = x0;
    return hg;
}


double specfun_chgu(double x, double a, double b, int *md, int *isfer) {

    // =======================================================
    // Purpose: Compute the confluent hypergeometric function
    //          U(a,b,x)
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument  ( x > 0 )
    // Output:  HU --- U(a,b,x)
    //          MD --- Method code
    //          ISFER --- Error flag
    // Routines called:
    //      (1) CHGUS for small x ( MD=1 )
    //      (2) CHGUL for large x ( MD=2 )
    //      (3) CHGUBI for integer b ( MD=3 )
    //      (4) CHGUIT for numerical integration ( MD=4 )
    // =======================================================

    int il1, il2, il3, bl1, bl2, bl3, bn, id1 = 0, id;
    double aa, hu = 0.0, hu1;

    aa = a - b + 1.0;
    *isfer = 0;
    il1 = (a == (int)a) && (a <= 0.0);
    il2 = (aa == (int)aa) && (aa <= 0.0);
    il3 = fabs(a*(a-b+1.0))/x <= 2.0;
    bl1 = (x <= 5.0) || (x <= 10.0 && a <= 2.0);
    bl2 = (x > 5.0) && (x <= 12.5) && ((a >= 1.0) && (b >= a+4.0));
    bl3 = (x > 12.5) && (a >= 5.0) && (b >= a + 5.0);
    bn = (b == (int)b) && (b != 0.0);

    id = -100;
    hu1 = 0.0;
    if (b != (int)b) {
        hu = specfun_chgus(x, a, b, &id1);
        *md = 1;
        if (id1 >= 9) { return hu; }
        hu1 = hu;
    }
    if (il1 || il2 || il3) {
        hu = specfun_chgul(x, a, b, &id);
        *md = 2;
        if (id >= 9) { return hu; }
        if (id1 > id) {
            *md = 1;
            id = id1;
            hu = hu1;
        }
    }
    if (a >= 1.0) {
        if (bn && (bl1 || bl2 || bl3)) {
            hu = specfun_chgubi(x, a, b, &id);
            *md = 3;
        } else {
            hu = specfun_chguit(x, a, b, &id);
            *md = 4;
        }
    } else {
        if (b <= a) {
            a -= b - 1.0;
            b = 2.0 - b;
            hu = specfun_chguit(x, a, b, &id);
            hu *= pow(x, 1.0 - b);
            *md = 4;
        } else if (bn && (~il1)) {
            hu = specfun_chgubi(x, a, b, &id);
            *md = 3;
        }
    }
    if (id < 6) { *isfer = 6; }
    return hu;
}


double specfun_chgubi(double x, double a, double b, int *id) {

    // ======================================================
    // Purpose: Compute confluent hypergeometric function
    //          U(a,b,x) with integer b ( b = ±1,±2,... )
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // Routines called:
    //      (1) GAMMA2 for computing gamma function Г(x)
    //      (2) PSI_SPEC for computing psi function
    // ======================================================

    int id1, id2, j, k, m, n;
    double a0, a1, a2, da1, da2, db1, db2, ga, ga1, h0, hm1, hm2, hm3,\
           hmax, hmin, hu, hu1, hu2, hw, ps, r, rn, rn1, s0, s1, s2,\
           sa, sb, ua, ub;
    const double el = 0.5772156649015329;

    *id = -100;
    n = (int)fabs(b-1);
    rn1 = 1.0;
    rn = 1.0;
    for (j = 1; j <= n; j++) {
        rn *= j;
        if (j == n-1) {
            rn1 = rn;
        }
    }
    ps = specfun_psi_spec(a);
    ga = specfun_gamma2(a);
    if (b > 0.0) {
        a0 = a;
        a1 = a - n;
        a2 = a1;
        ga1 = specfun_gamma2(a1);
        ua = pow(-1, n-1) / (rn * ga1);
        ub = rn1 / ga * pow(x, -n);
    } else {
        a0 = a + n;
        a1 = a0;
        a2 = a;
        ga1 = specfun_gamma2(a1);
        ua = pow(-1, n-1) / (rn * ga) * pow(x, n);
        ub = rn1 / ga1;
    }
    hm1 = 1.0;
    r = 1.0;
    hmax = 0.0;
    hmin = 1e300;
    h0 = 0.0;
    for (k = 1; k <= 150; k++) {
        r = r * (a0 + k - 1) * x / ((n + k) * k);
        hm1 += r;
        hu1 = fabs(hm1);

        if (hu1 > hmax) {
            hmax = hu1;
        }
        if (hu1 < hmin) {
            hmin = hu1;
        }
        if (fabs(hm1 - h0) < fabs(hm1) * 1.0e-15) { break; }
        h0 = hm1;
    }

    da1 = log10(hmax);
    da2 = 0;

    if (hmin != 0) {
        da2 = log10(hmin);
    }

    *id = 15 - (int)fabs(da1 - da2);
    hm1 *= log(x);
    s0 = 0;

    for (m = 1; m <= n; m++) {
        if (b >= 0) {
            s0 -= 1.0 / m;
        }
        if (b < 0) {
            s0 += (1.0 - a) / (m * (a + m - 1));
        }
    }

    hm2 = ps + 2 * el + s0;
    r = 1;
    hmax = 0;
    hmin = 1.0e+300;

    for (k = 1; k <= 150; k++) {
        s1 = 0;
        s2 = 0;

        if (b > 0) {
            for (m = 1; m <= k; m++) {
                s1 -= (m + 2 * a - 2) / (m * (m + a - 1));
            }
            for (m = 1; m <= n; m++) {
                s2 += 1.0 / (k + m);
            }
        } else {
            for (m = 1; m <= k + n; m++) {
                s1 += (1.0 - a) / (m * (m + a - 1));
            }
            for (m = 1; m <= k; m++) {
                s2 += 1.0 / m;
            }
        }

        hw = 2 * el + ps + s1 - s2;
        r = r * (a0 + k - 1) * x / ((n + k) * k);
        hm2 += r * hw;
        hu2 = fabs(hm2);

        if (hu2 > hmax) {
            hmax = hu2;
        }

        if (hu2 < hmin) {
            hmin = hu2;
        }

        if (fabs((hm2 - h0) / hm2) < 1.0e-15) {
            break;
        }

        h0 = hm2;
    }

    db1 = log10(hmax);
    db2 = 0.0;
    if (hmin != 0.0) { db2 = log10(hmin); }
    id1 = 15 - (int)fabs(db1 - db2);
    if (id1 < *id) { *id = id1; }
    hm3 = 1.0;
    if (n == 0) { hm3 = 0.0; }
    r = 1.0;
    for (k = 1; k < n; k++) {
        r = r * (a2 + k - 1.0) / ((k - n)*k)*x;
        hm3 += r;
    }
    sa = ua*(hm1 + hm2);
    sb = ub*hm3;
    hu = sa + sb;
    id2 = 0;
    if (sa != 0.0) { id1 = (int)(log10(fabs(sa))); }
    if (hu != 0.0) { id2 = (int)(log10(fabs(hu))); }
    if (sa*sb < 0.0) { *id -= abs(id1-id2); }
    return hu;
}


double specfun_chguit(double x, double a, double b, int *id) {

    // ======================================================
    // Purpose: Compute hypergeometric function U(a,b,x) by
    //          using Gaussian-Legendre integration (n=60)
    // Input  : a  --- Parameter ( a > 0 )
    //          b  --- Parameter
    //          x  --- Argument ( x > 0 )
    // Output:  HU --- U(a,b,z)
    //          ID --- Estimated number of significant digits
    // Routine called: GAMMA2 for computing Г(x)
    // ======================================================

    int k, j, m;
    double a1, b1, c, d, f1, f2, g, ga, hu, hu0, hu1, hu2, s, t1, t2, t3, t4;
    static const double t[30] = {
        0.259597723012478e-01, 0.778093339495366e-01, 0.129449135396945e+00, 0.180739964873425e+00,
        0.231543551376029e+00, 0.281722937423262e+00, 0.331142848268448e+00, 0.379670056576798e+00,
        0.427173741583078e+00, 0.473525841761707e+00, 0.518601400058570e+00, 0.562278900753945e+00,
        0.604440597048510e+00, 0.644972828489477e+00, 0.683766327381356e+00, 0.720716513355730e+00,
        0.755723775306586e+00, 0.788693739932264e+00, 0.819537526162146e+00, 0.848171984785930e+00,
        0.874519922646898e+00, 0.898510310810046e+00, 0.920078476177628e+00, 0.939166276116423e+00,
        0.955722255839996e+00, 0.969701788765053e+00, 0.981067201752598e+00, 0.989787895222222e+00,
        0.995840525118838e+00, 0.999210123227436e+00
    };
    static const double w[30] = {
        0.519078776312206e-01, 0.517679431749102e-01, 0.514884515009810e-01, 0.510701560698557e-01,
        0.505141845325094e-01, 0.498220356905502e-01, 0.489955754557568e-01, 0.480370318199712e-01,
        0.469489888489122e-01, 0.457343797161145e-01, 0.443964787957872e-01, 0.429388928359356e-01,
        0.413655512355848e-01, 0.396806954523808e-01, 0.378888675692434e-01, 0.359948980510845e-01,
        0.340038927249464e-01, 0.319212190192963e-01, 0.297524915007890e-01, 0.275035567499248e-01,
        0.251804776215213e-01, 0.227895169439978e-01, 0.203371207294572e-01, 0.178299010142074e-01,
        0.152746185967848e-01, 0.126781664768159e-01, 0.100475571822880e-01, 0.738993116334531e-02,
        0.471272992695363e-02, 0.202681196887362e-02
    };
    *id = 9;
    // DLMF 13.4.4, integration up to C=12/X
    a1 = a - 1.0;
    b1 = b - a - 1.0;
    c = 12.0 / x;
    hu0 = 0.0;
    for (m = 10; m <= 100; m += 5) {
        hu1 = 0.0;
        g=0.5 * c / m;
        d=g;
        for (j = 1; j < (m + 1); j++) {
            s = 0.0;
                for (k = 1; k <= 30; k++) {
                    t1 = d + g * t[k-1];
                    t2 = d - g * t[k-1];
                    f1 = exp(-x*t1) * pow(t1, a1) * pow(1.0 + t1, b1);
                    f2 = exp(-x*t2) * pow(t2, a1) * pow(1.0 + t2, b1);
                    s += w[k-1]*(f1 + f2);
                }
            hu1 += s * g;
            d += 2.0 * g;
        }
        if (fabs(1.0 - hu0/hu1) < 1.0e-9) { break; }
        hu0 = hu1;
    }
    ga = specfun_gamma2(a);
    hu1 /= ga;
    // DLMF 13.4.4 with substitution t=C/(1-u)
    // integration u from 0 to 1, i.e. t from C=12/X to infinity
    for (m = 2; m <= 10; m += 2) {
        hu2 = 0.0;
        g = 0.5 / m;
        d = g;
        for (j = 1; j <= m; j++) {
            s = 0.0;
            for (k = 1; k <= 30; k++) {
                t1 = d + g * t[k-1];
                t2 = d - g * t[k-1];
                t3 = c / (1.0 - t1);
                t4 = c / (1.0 - t2);
                f1 = t3*t3 / c * exp(-x*t3)*pow(t3, a1)*pow(1.0 + t3, b1);
                f2 = t4*t4 / c * exp(-x*t4)*pow(t4, a1)*pow(1.0 + t4, b1);
                s += w[k-1]*(f1 + f2);
            }
            hu2 += s*g;
            d += 2.0*g;
        }
        if (fabs(1.0 - hu0/hu2) < 1.0e-9) { break; }
        hu0 = hu2;
    }
    ga = specfun_gamma2(a);
    hu2 /= ga;
    hu = hu1 + hu2;
    return hu;
}


double specfun_chgul(double x, double a, double b, int *id) {

    // =======================================================
    // Purpose: Compute the confluent hypergeometric function
    //          U(a,b,x) for large argument x
    // Input  : a  --- Parameter
    //          b  --- Parameter
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // =======================================================

    int il1, il2, k, nm;
    double aa, hu, r, r0 = 0.0, ra = 0.0;

    *id = -100;
    aa = a - b + 1.0;
    il1 = (a == (int)a) && (a <= 0.0);
    il2 = (aa == (int)aa) && (aa <= 0.0);
    nm = 0;
    if (il1) { nm = (int)fabs(a); }
    if (il2) { nm = (int)fabs(aa); }
    // IL1: DLMF 13.2.7 with k=-s-a
    // IL2: DLMF 13.2.8
    if (il1 || il2) {
        hu = 1.0;
        r = 1.0;
        for (k = 1; k <= nm; k++) {
            r = -r*(a + k - 1.0)*(a - b + k) / (k*x);
            hu += r;
        }
        hu *= pow(x, -a);
        *id = 10;
    } else {
        // DLMF 13.7.3
        hu = 1.0;
        r = 1.0;
        for (k = 1; k <= 25; k++) {
            r = -r*(a + k - 1.0)*(a - b + k) / (k*x);
            ra = fabs(r);
            if (((k > 5) && (ra >= r0)) || (ra < 1e-15)) { break; }
            r0 = ra;
            hu += r;
        }
        *id = (int)fabs(log10(ra));
        hu *= pow(x, -a);
    }
    return hu;
}


double specfun_chgus(double x, double a, double b, int *id) {

    // ======================================================
    // Purpose: Compute confluent hypergeometric function
    //          U(a,b,x) for small argument x
    // Input  : a  --- Parameter
    //          b  --- Parameter ( b <> 0,-1,-2,...)
    //          x  --- Argument
    // Output:  HU --- U(a,b,x)
    //          ID --- Estimated number of significant digits
    // Routine called: GAMMA2 for computing gamma function
    // ======================================================

    // DLMF 13.2.42 with prefactors rewritten according to
    // DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2
    int j;
    double d1, d2, ga, gb, gab, gb2, h0, hmax, hmin, hu, hu0, hua, r1, r2;
    const double pi = 3.141592653589793;

    *id = 100;
    ga = specfun_gamma2(a);
    gb = specfun_gamma2(b);
    gab = specfun_gamma2(1.0 + a - b);
    gb2 = specfun_gamma2(2.0 - b);
    hu0 = pi / sin(pi*b);
    r1 = hu0 / (gab*gb);
    r2 = hu0*pow(x, 1.0 - b) / (ga*gb2);
    hu = r1 - r2;
    hmax = 0.0;
    hmin = 1e300;
    h0 = 0.0;
    for (j = 1; j < 151; j++) {
        r1 = r1*(a + j - 1.0) / (j*(b + j - 1.0))*x;
        r2 = r2*(a - b + j) / (j*(1.0 - b + j))*x;
        hu += r1 - r2;
        hua = fabs(hu);
        if (hua > hmax) { hmax = hua; }
        if (hua < hmin) { hmin = hua; }
        if (fabs(hu - h0) < fabs(hu)*1e-15) { break; }
        h0 = hu;
    }
    d1 = log10(hmax);
    d2 = 0.0;
    if (hmin != 0.0) { d2 = log10(hmin); }
    *id = 15 - (d1 - d2 < 0 ? d2 - d1 : d1 - d2);
    return hu;
}


void specfun_clpmn(double complex z, int m, int n, int ntype, double complex *cpm, double complex *cpd) {

    // =========================================================
    // Purpose: Compute the associated Legendre functions Pmn(z)
    //          and their derivatives Pmn'(z) for a complex
    //          argument
    // Input :  x     --- Real part of z
    //          y     --- Imaginary part of z
    //          m     --- Order of Pmn(z),  m = 0,1,2,...,n
    //          n     --- Degree of Pmn(z), n = 0,1,2,...,N
    //          mm    --- Physical dimension of CPM and CPD
    //          ntype --- type of cut, either 2 or 3
    // Output:  CPM(m,n) --- Pmn(z)
    //          CPD(m,n) --- Pmn'(z)
    //
    // SciPy mod: C translation uses a contiguous memory block
    // =========================================================

    int i, j, ls;
    double complex zq, zs;
    double x = creal(z);
    double y = cimag(z);

    for (i = 0; i < (m+1)*(n+1); i++) {
            cpm[i] = 0.0;
            cpd[i] = 0.0;
    }
    cpm[0] = 1.0;
    if (n == 0) {
        return;
    }
    if ((fabs(x) == 1.0) && (y == 0.0)) {
        for (i = 1; i <= n; i++) {
            cpm[i] = pow(x, i);
            cpd[i] = 0.5*i*(i+1)*pow(x, i+1);
        }
        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                if (i == 1) {
                    cpd[i*(n+1) + j] = INFINITY;
                } else if (i == 2) {
                    cpd[i*(n+1) + j] = -0.25*(j+2)*(j+1)*j*(j-1)*pow(x, j+1);
                }
            }
        }
        return;
    }
    if (ntype == 2) {
        // sqrt(1 - z**2) with branch cut on |x|>1
        zs = (1.0 - z*z);
        zq = -csqrt(zs);
        ls = -1;
    } else {
        // sqrt(z**2 - 1) with branch cut between [-1, 1]
        zs = (z*z - 1.0);
        zq = csqrt(zs);
        if (x < 0.) { zq = -zq; }
        ls = 1;
    }
    for (i = 1; i <= m; i++) {
        // DLMF 14.7.15
        cpm[i*(n + 2)] = (2.*i - 1.)*zq*cpm[(i-1)*(n+2)];
    }
    for (i = 0; i <= (m > n-1 ? n-1 : m); i++) {
        // DLMF 14.10.7
        cpm[i*(n + 2) + 1] = (2.*i + 1)*z*cpm[i*(n + 2)];
    }
    for (i = 0; i <= m; i++) {
        for (j = i+2; j <= n; j++) {
            // DLMF 14.10.3
            cpm[i*(n + 1) + j] = ((2.*j - 1)*z*cpm[i*(n + 1) + j-1] - (i+j-1)*cpm[i*(n+1) + j-2])/(j-i);
        }
    }
    cpd[0] = 0.0;
    for (j = 1; j <= n; j++) {
        // DLMF 14.10.5
        cpd[j] = ls*j*(z*cpm[j] - cpm[j-1])/zs;
    }
    for (i = 1; i <= m; i++) {
        for (j = i; j <= n; j++) {
            // derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
            // derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
            cpd[i * (n + 1) + j] = ls*(-i*z*cpm[i * (n + 1) + j]/zs +
                                   (j+i)*(j-i+1.0)/zq*cpm[(i - 1)*(n + 1) + j]);
        }
    }
    return;
}


void specfun_clpn(int n, double complex z, double complex *cpn, double complex *cpd) {

    // ==================================================
    // Purpose: Compute Legendre polynomials Pn(z) and
    //          their derivatives Pn'(z) for a complex
    //          argument
    // Input :  x --- Real part of z
    //          y --- Imaginary part of z
    //          n --- Degree of Pn(z), n = 0,1,2,...
    // Output:  CPN(n) --- Pn(z)
    //          CPD(n) --- Pn'(z)
    // ==================================================

    int k;
    double complex cp0, cp1, cpf;

    cpn[0] = 1.0;
    cpn[1] = z;
    cpd[0] = 0.0;
    cpd[1] = 1.0;
    cp0 = 1.0;
    cp1 = z;
    for (k = 2; k <= n; k++) {
        cpf = (2.0 * k -1.0) / k * z * cp1 - (k - 1.0) / k * cp0;
        cpn[k] = cpf;
        if (z == 1.0) {
            cpd[k] = 0.5 * pow(creal(z) , k+1) * k * (k + 1.0);
        } else {
            cpd[k] = k * (cp1 - z * cpf) / (1.0 - z * z);
        }
        cp0 = cp1;
        cp1 = cpf;
    }
    return;
}


void specfun_clqmn(double complex z, int m, int n, double complex *cqm, double complex *cqd) {

    // =======================================================
    // Purpose: Compute the associated Legendre functions of
    //          the second kind, Qmn(z) and Qmn'(z), for a
    //          complex argument
    // Input :  x  --- Real part of z
    //          y  --- Imaginary part of z
    //          m  --- Order of Qmn(z)  ( m = 0,1,2,… )
    //          n  --- Degree of Qmn(z) ( n = 0,1,2,… )
    //          mm --- Physical dimension of CQM and CQD
    // Output:  CQM(m,n) --- Qmn(z)
    //          CQD(m,n) --- Qmn'(z)
    // =======================================================

    int i, j, k, km, ls;
    double xc;
    double complex cq0, cq1, cq10, cqf0 = 0.0, cqf, cqf1, cqf2, zq, zs;
    double x = creal(z);
    double y = cimag(z);
    if ((fabs(x) == 1.0) && (y == 0.0)) {
        for (i = 0; i < (m + 1) * (n + 1); i++) {
            cqm[i] = 1e300;
            cqd[i] = 1e300;
        }
        return;
    }
    xc = cabs(z);
    ls = 0;
    if ((cimag(z) == 0.0) || (xc < 1.0)) {
        ls = 1;
    }
    if (xc > 1.0) {
        ls = -1;
    }
    zs = ls*(1.0 - z*z);
    zq = csqrt(zs);
    cq0 = 0.5*clog(ls*(1.0 + z)/(1.0 - z));

    if (xc < 1.0001) {
        cqm[0] = cq0;
        cqm[1] = z*cq0 - 1.0;
        cqm[n + 1] = -1.0 / zq;
        cqm[n + 2] = -zq*(cq0 + z / (1.0 - z*z));
        for (i = 0; i <= 1; i++) {
            for (j = 2; j <= n; j++) {
                cqm[i * (n + 1) + j] = ((2.0*j-1.0)*z*cqm[i * (n + 1) + j - 1]
                                       -(j+i-1.0)*cqm[i * (n + 1) + j - 2])/(j-i);
            }
        }
        for (i = 2; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                cqm[i * (n + 1) + j] = -2.0*(i-1.0)*z/zq*cqm[(i - 1) * (n + 1) + j]
                                       - ls*(j+i-1.0)*(j-i+2.0)*cqm[(i - 2) * (n + 1) + j];
            }
        }
    } else {
        if (xc > 1.1) {
            km = 40 + m + n;
        } else {
            km = (40 + m + n)*((int)(-1.0 - 1.8*log(xc - 1.)));
        }
        cqf2 = 0.0;
        cqf1 = 1.0;
        for (k = km; k >= 0; k--) {
            cqf0 = ((2*k+3.0)*z*cqf1 - (k+2.0)*cqf2) / (k+1.0);
            if (k <= n) {
                cqm[k] = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }
        /* 25 */
        for (k = 0; k <= n; k++) {
            cqm[k] *= cq0 / cqf0;
        }
        /* 30 */
        cqf2 = 0.0;
        cqf1 = 1.0;
        for (k = km; k >= 0; k--) {
            cqf0 = ((2*k+3.0)*z*cqf1 - (k+1.0)*cqf2) / (k+2.0);
            if (k <= n) {
                cqm[n + 1 + k] = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }
        /* 35 */
        cq10 = -1.0 / zq;
        for (k = 0; k <= n; k++) {
            cqm[n + 1 + k] *= cq10 / cqf0;
        }
        for (j = 0; j <= n; j++) {
            cq0 = cqm[j];
            cq1 = cqm[n + 1 + j];
            for (i = 0; i <= (m-2); i++) {
                cqf = -2.0*(i+1)*z/zq*cq1 + (j-i)*(j+i+1.0)*cq0;
                cqm[(i + 2)*(n + 1) + j] = cqf;
                cq0 = cq1;
                cq1 = cqf;
            }
        }
        cqd[0] = ls / zs;
        for (j = 1; j <= n; j++) {
            cqd[j] = ls*j*(cqm[j-1] - z*cqm[j])/zs;
        }
        /* 50 */
        for (i = 1; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                cqd[i*(n + 1) + j] = ls*i*z/zs*cqm[i*(n + 1) + j] + (i+j)*(j-i+1.0)/zq*cqm[(i - 1)*(n + 1) + j];
            }
        }
        return;
    }
}


void specfun_clqn(int n, double complex z, double complex *cqn, double complex *cqd) {

    // ==================================================
    // Purpose: Compute the Legendre functions Qn(z) and
    //          their derivatives Qn'(z) for a complex
    //          argument
    // Input :  x --- Real part of z
    //          y --- Imaginary part of z
    //          n --- Degree of Qn(z), n = 0,1,2,...
    // Output:  CQN(n) --- Qn(z)
    //          CQD(n) --- Qn'(z)
    // ==================================================

    int k, km, ls;
    double complex cq0, cq1, cqf0 = 0.0, cqf1, cqf2;

    if (z == 1.0) {
        for (int k = 0; k <= n; ++k) {
            cqn[k] = 1e300;
            cqd[k] = 1e300;
        }
        return;
    }
    ls = ((cabs(z) > 1.0) ? -1 : 1);

    cq0 = 0.5 * clog(ls * (1.0 + z) / (1.0 - z));
    cq1 = z * cq0 - 1.0;

    cqn[0] = cq0;
    cqn[1] = cq1;

    if (cabs(z) < 1.0001) {
        cqf0 = cq0;
        cqf1 = cq1;
        for (k = 2; k <= n; k++) {
            cqf2 = ((2.0 * k - 1.0) * z * cqf1 - (k - 1.0) * cqf0) / k;
            cqn[k] = cqf2;
            cqf0 = cqf1;
            cqf1 = cqf2;
        }
    } else {
        if (cabs(z) > 1.1) {
            km = 40 + n;
        } else {
            km = (int)((40 + n) * floor(-1.0 - 1.8 * log(cabs(z - 1.0))));
        }

        cqf2 = 0.0;
        cqf1 = 1.0;
        for (int k = km; k >= 0; k--) {
            cqf0 = ((2 * k + 3.0) * z * cqf1 - (k + 2.0) * cqf2) / (k + 1.0);
            if (k <= n) {
                cqn[k] = cqf0;
            }
            cqf2 = cqf1;
            cqf1 = cqf0;
        }
        for (int k = 0; k <= n; ++k) {
            cqn[k] *= cq0 / cqf0;
        }
    }
    cqd[0] = (cqn[1] - z * cqn[0]) / (z * z - 1.0);

    for (int k = 1; k <= n; ++k) {
        cqd[k] = (k * z * cqn[k] - k * cqn[k - 1]) / (z * z - 1.0);
    }
    return;
}


void specfun_cpbdn(int n, double complex z, double complex *cpb, double complex *cpd) {

    // ==================================================
    // Purpose: Compute the parabolic cylinder functions
    //           Dn(z) and Dn'(z) for a complex argument
    // Input:   z --- Complex argument of Dn(z)
    //          n --- Order of Dn(z)  ( n=0,±1,±2,… )
    // Output:  CPB(|n|) --- Dn(z)
    //          CPD(|n|) --- Dn'(z)
    // Routines called:
    //      (1) CPDSA for computing Dn(z) for a small |z|
    //      (2) CPDLA for computing Dn(z) for a large |z|
    // ==================================================

    int n0, n1, nm1;
    double a0, x;
    double complex ca0, cf, cf0, cf1, cfa, cfb, cs0, z1;
    const double pi = 3.141592653589793;

    x = creal(z);
    a0 = cabs(z);
    ca0 = cexp(-0.25 * z * conj(z));
    n0 = 0;

    if (n >= 0) {
        cf0 = ca0;
        cf1 = z * ca0;

        cpb[0] = cf0;
        cpb[1] = cf1;

        for (int k = 2; k <= n; ++k) {
            cf = z * cf1 - (k - 1.0) * cf0;
            cpb[k] = cf;
            cf0 = cf1;
            cf1 = cf;
        }
    } else {
        n0 = -n;

        if (x <= 0.0 || a0 == 0.0) {
            cf0 = ca0;
            cpb[0] = cf0;

            z1 = -z;

            if (a0 <= 7.0) {
                cpb[1] = specfun_cpdsa(-1, z1);
            } else {
                cpb[1] = specfun_cpdla(-1, z1);
            }

            cf1 = csqrt(2.0 * pi) / ca0 - cpb[1];
            cpb[1] = cf1;

            for (int k = 2; k < n0; ++k) {
                cf = (-z * cf1 + cf0) / (k - 1.0);
                cpb[k] = cf;
                cf0 = cf1;
                cf1 = cf;
            }
        } else if (a0 <= 3.0) {
            cpb[n0] = specfun_cpdsa(-n0, z);
            n1 = n0 + 1;
            cpb[n1] = specfun_cpdsa(-n1, z);

            nm1 = n0 - 1;
            for (int k = nm1; k >= 0; --k) {
                cf = z * cpb[n0] + (k + 1.0) * cpb[n1];
                cpb[k] = cf;
                cpb[n1] = cpb[n0];
                cpb[n0] = cf;
            }
        } else {
            int m = 100 + abs(n);
            cfa = 0.0;
            cfb = 1.0e-30 + 0.0 * I;

            for (int k = m; k >= 0; --k) {
                cf = z * cfb + (k + 1.0) * cfa;

                if (k <= n0) {
                    cpb[k] = cf;
                }

                cfa = cfb;
                cfb = cf;
            }

            cs0 = ca0 / cfb;

            for (int k = 0; k <= n0; ++k) {
                cpb[k] = cs0 * cpb[k];
            }
        }
    }

    cpd[0] = -0.5 * z * cpb[0];

    if (n >= 0) {
        for (int k = 1; k <= n; ++k) {
            cpd[k] = -0.5 * z * cpb[k] + k * cpb[k - 1];
        }
    } else {
        for (int k = 1; k < n0; ++k) {
            cpd[k] = 0.5 * z * cpb[k] - cpb[k - 1];
        }
    }
}


double complex specfun_cpdla(int n, double complex z) {

    // ===========================================================
    // Purpose: Compute complex parabolic cylinder function Dn(z)
    //          for large argument
    // Input:   z   --- Complex argument of Dn(z)
    //          n   --- Order of Dn(z) (n = 0,±1,±2,…)
    // Output:  CDN --- Dn(z)
    // ===========================================================

    int k;
    double complex cb0, cr, cdn;

    cb0 = cpow(z, n)*cexp(-0.25*z*z);
    cr = 1.0;
    cdn = 1.0;
    for (k = 1; k <= 16; k++) {
        cr = - 0.5 * cr * (2.0 * k - n - 1.0) * (2.0 * k - n - 2.0) / (k * z * z);
        cdn += cr;
        if (cabs(cr) < cabs(cdn) * 1e-12) { break; }
    }
    return cdn * cb0;
}


double complex specfun_cpdsa(int n, double complex z) {

    // ===========================================================
    // Purpose: Compute complex parabolic cylinder function Dn(z)
    //          for small argument
    // Input:   z   --- complex argument of D(z)
    //          n   --- Order of D(z) (n = 0,-1,-2,...)
    // Output:  CDN --- Dn(z)
    // Routine called: GAIH for computing Г(x), x=n/2 (n=1,2,...)
    // ===========================================================

    int m;
    double va0, pd, vm, vt, xn;
    double complex ca0, cb0, cdn, cr, cdw, g0, g1, ga0, gm;
    const double eps = 1.0e-15;
    const double pi = 3.141592653589793;
    const double sq2 = sqrt(2.0);

    ca0 = cexp(-0.25 * z * z);
    va0 = 0.5 * (1.0 - n);
    if (n == 0.0) {
        cdn = ca0;
    } else {
        if (cabs(z) == 0.0) {
            if ((va0 <= 0.0) && (va0 == (int)va0)) {
                cdn = 0.0;
            } else {
                ga0 = specfun_gaih(va0);
                pd = sqrt(pi) / (pow(2.0, -0.5 * n) * ga0);
                cdn = pd;
            }
        } else {
            xn = -n;
            g1 = specfun_gaih(xn);
            cb0 = pow(2.0, -0.5 * n - 1.0) * ca0 / g1;
            vt = -0.5 * n;
            g0 = specfun_gaih(vt);
            cdn = g0;
            cr = CMPLX(1.0, 0.0);

            for (m = 1; m <= 250; m++) {
                vm = 0.5 * (m - n);
                gm = specfun_gaih(vm);
                cr = -cr*sq2 * z / m;
                cdw = gm * cr;
                cdn += cdw;
                if (cabs(cdw) < cabs(cdn) * eps) {
                    break;
                }
            }
            cdn *= cb0;
        }
    }
    return cdn;
}


double specfun_cv0(double kd, double m, double q) {

    // =====================================================
    // Purpose: Compute the initial characteristic value of
    //          Mathieu functions for m ≤ 12  or q ≤ 300 or
    //          q ≥ m*m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Characteristic value
    // Routines called:
    //       (1) CVQM for computing initial characteristic
    //           value for q ≤ 3*m
    //       (2) CVQL for computing initial characteristic
    //           value for q ≥ m*m
    // ====================================================

    double a0 = 0.0, q2 = q*q;

    if (m == 0) {
        if (q <= 1.0) {
            a0 = (((0.0036392 * q2 - 0.0125868) * q2 + 0.0546875) * q2 - 0.5) * q2;
        } else if (q <= 10.0) {
            a0 = ((3.999267e-3 * q - 9.638957e-2) * q - 0.88297) * q + 0.5542818;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 1) {
        if ((q <= 1.0) && (kd == 2)) {
            a0 = (((-6.51e-4 * q - 0.015625) * q - 0.125) * q + 1.0) * q + 1.0;
        } else if (q <= 1.0 && kd == 3) {
            a0 = (((-6.51e-4 * q + 0.015625) * q - 0.125) * q - 1.0) * q + 1.0;
        } else if (q <= 10.0 && kd == 2) {
            a0 = (((-4.94603e-4 * q + 1.92917e-2) * q - 0.3089229) * q + 1.33372) * q + 0.811752;
        } else if (q <= 10.0 && kd == 3) {
            a0 = ((1.971096e-3 * q - 5.482465e-2) * q - 1.152218) * q + 1.10427;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 2) {
        if (q <= 1.0 && kd == 1) {
            a0 = (((-0.0036391 * q2 + 0.0125888) * q2 - 0.0551939) * q2 + 0.416667) * q2 + 4.0;
        } else if (q <= 1.0 && kd == 4) {
            a0 = (0.0003617 * q2 - 0.0833333) * q2 + 4.0;
        } else if (q <= 15.0 && kd == 1) {
            a0 = (((3.200972e-4 * q - 8.667445e-3) * q - 1.829032e-4) * q + 0.9919999) * q + 3.3290504;
        } else if (q <= 10.0 && kd == 4) {
            a0 = ((2.38446e-3 * q - 0.08725329) * q - 4.732542e-3) * q + 4.00909;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 3) {
        if (q <= 1.0 && kd == 2) {
            a0 = (((6.348e-4 * q + 0.015625) * q + 0.0625) * q2 + 9.0);
        } else if (q <= 1.0 && kd == 3) {
            a0 = (((6.348e-4 * q - 0.015625) * q + 0.0625) * q2 + 9.0);
        } else if (q <= 20.0 && kd == 2) {
            a0 = (((3.035731e-4 * q - 1.453021e-2) * q + 0.19069602) * q - 0.1039356) * q + 8.9449274;
        } else if (q <= 15.0 && kd == 3) {
            a0 = ((9.369364e-5 * q - 0.03569325) * q + 0.2689874) * q + 8.771735;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 4) {
        if (q <= 1.0 && kd == 1) {
            a0 = ((-2.1e-6 * q2 + 5.012e-4) * q2 + 0.0333333) * q2 + 16.0;
        } else if (q <= 1.0 && kd == 4) {
            a0 = ((3.7e-6 * q2 - 3.669e-4) * q2 + 0.0333333) * q2 + 16.0;
        } else if (q <= 25.0 && kd == 1) {
            a0 = (((1.076676e-4 * q - 7.9684875e-3) * q + 0.17344854) * q - 0.5924058) * q + 16.620847;
        } else if (q <= 20.0 && kd == 4) {
            a0 = ((-7.08719e-4 * q + 3.8216144e-3) * q + 0.1907493) * q + 15.744;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 5) {
        if (q <= 1.0 && kd == 2) {
            a0 = ((6.8e-6 * q + 1.42e-5) * q2 + 0.0208333) * q2 + 25.0;
        } else if (q <= 1.0 && kd == 3) {
            a0 = ((-6.8e-6 * q + 1.42e-5) * q2 + 0.0208333) * q2 + 25.0;
        } else if (q <= 35.0 && kd == 2) {
            a0 = (((2.238231e-5 * q - 2.983416e-3) * q + 0.10706975) * q - 0.600205) * q + 25.93515;
        } else if (q <= 25.0 && kd == 3) {
            a0 = ((-7.425364e-4 * q + 2.18225e-2) * q + 4.16399e-2) * q + 24.897;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 6) {
        if (q <= 1.0) {
            a0 = (4e-6 * q2 + 0.0142857) * q2 + 36.0;
        } else if (q <= 40.0 && kd == 1) {
            a0 = (((-1.66846e-5 * q + 4.80263e-4) * q + 2.53998e-2) * q - 0.181233) * q + 36.423;
        } else if (q <= 35.0 && kd == 4) {
            a0 = ((-4.57146e-4 * q + 2.16609e-2) * q - 2.349616e-2) * q + 35.99251;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m == 7) {
        if (q <= 10.0) {
            a0 = specfun_cvqm(m, q);
        } else if (q <= 50.0 && kd == 2) {
            a0 = (((-1.411114e-5 * q + 9.730514e-4) * q - 3.097887e-3) * q + 3.533597e-2) * q + 49.0547;
        } else if (q <= 40.0 && kd == 3) {
            a0 = ((-3.043872e-4 * q + 2.05511e-2) * q - 9.16292e-2) * q + 49.19035;
        } else {
            a0 = specfun_cvql(kd, m, q);
        }
    } else if (m >= 8) {
        if (q <= 3*m) {
            a0 = specfun_cvqm(m, q);
        } else if (q > m * m) {
            a0 = specfun_cvql(kd, m, q);
        } else {
            if (m == 8 && kd == 1) {
                a0 = (((8.634308e-6 * q - 2.100289e-3) * q + 0.169072) * q - 4.64336) * q + 109.4211;
            } else if (m == 8 && kd == 4) {
                a0 = ((-6.7842e-5 * q + 2.2057e-3) * q + 0.48296) * q + 56.59;
            } else if (m == 9 && kd == 2) {
                a0 = (((2.906435e-6 * q - 1.019893e-3) * q + 0.1101965) * q - 3.821851) * q + 127.6098;
            } else if (m == 9 && kd == 3) {
                a0 = ((-9.577289e-5 * q + 0.01043839) * q + 0.06588934) * q + 78.0198;
            } else if (m == 10 && kd == 1) {
                a0 = (((5.44927e-7 * q - 3.926119e-4) * q + 0.0612099) * q - 2.600805) * q + 138.1923;
            } else if (m == 10 && kd == 4) {
                a0 = ((-7.660143e-5 * q + 0.01132506) * q - 0.09746023) * q + 99.29494;
            } else if (m == 11 && kd == 2) {
                a0 = (((-5.67615e-7 * q + 7.152722e-6) * q + 0.01920291) * q - 1.081583) * q + 140.88;
            } else if (m == 11 && kd == 3) {
                a0 = ((-6.310551e-5 * q + 0.0119247) * q - 0.2681195) * q + 123.667;
            } else if (m == 12 && kd == 1) {
                a0 = (((-2.38351e-7 * q - 2.90139e-5) * q + 0.02023088) * q - 1.289) * q + 171.2723;
            } else if (m == 12 && kd == 4) {
                a0 = (((3.08902e-7 * q - 1.577869e-4) * q + 0.0247911) * q - 1.05454) * q + 161.471;
            }
        }
    }
    return a0;
}


double specfun_cva2(int kd, int m, double q) {

    // ======================================================
    // Purpose: Calculate a specific characteristic value of
    //          Mathieu functions
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    //          KD --- Case code
    //                 KD=1 for cem(x,q)  ( m = 0,2,4,...)
    //                 KD=2 for cem(x,q)  ( m = 1,3,5,...)
    //                 KD=3 for sem(x,q)  ( m = 1,3,5,...)
    //                 KD=4 for sem(x,q)  ( m = 2,4,6,...)
    // Output:  A  --- Characteristic value
    // Routines called:
    //       (1) REFINE for finding accurate characteristic
    //           value using an iteration method
    //       (2) CV0 for finding initial characteristic
    //           values using polynomial approximation
    //       (3) CVQM for computing initial characteristic
    //           values for q ≤ 3*m
    //       (3) CVQL for computing initial characteristic
    //           values for q ≥ m*m
    // ======================================================

    int ndiv, nn, i;
    double a = 0.0, delta, q1, q2, qq, a1, a2;

    if ((m <= 12) || (q <= 3.0 * m) || (q > m * m)) {
        a = specfun_cv0(kd, m, q);
        if ((q != 0.0) && (m != 2)) { a = specfun_refine(kd, m, q, a); }
        if ((q > 2.0e-3) && (m == 2)) { a = specfun_refine(kd, m, q, a); }
    } else {
        ndiv = 10;
        delta = (m - 3.0) * m / ndiv;

        if ((q - 3.0 * m) <= (m * m - q)) {
            nn = (int)((q - 3.0 * m) / delta) + 1;
            delta = (q - 3.0 * m) / nn;
            q1 = 2.0 * m;
            a1 = specfun_cvqm(m, q1);
            q2 = 3.0 * m;
            a2 = specfun_cvqm(m, q2);
            qq = 3.0 * m;
            for (i = 1; i <= nn; i++) {
                qq = qq + delta;
                a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
                a = specfun_refine(kd, m, qq, a);
                q1 = q2;
                q2 = qq;
                a1 = a2;
                a2 = a;
            }
        } else {
            nn = (int)((m * m - q) / delta) + 1;
            delta = (m * m - q) / nn;
            q1 = m * (m - 1.0);
            a1 = specfun_cvql(kd, m, q1);
            q2 = m * m;
            a2 = specfun_cvql(kd, m, q2);
            qq = m * m;
            for (i = 1; i <= nn; ++i) {
                qq = qq - delta;
                a = (a1 * q2 - a2 * q1 + (a2 - a1) * qq) / (q2 - q1);
                a = specfun_refine(kd, m, qq, a);
                q1 = q2;
                q2 = qq;
                a1 = a2;
                a2 = a;
            }
        }
    }
    return a;
}


double specfun_cvf(int kd, int m, double q, double a, int mj) {

    // ======================================================
    // Purpose: Compute the value of F for characteristic
    //          equation of Mathieu functions
    // Input :  m --- Order of Mathieu functions
    //          q --- Parameter of Mathieu functions
    //          A --- Characteristic value
    // Output:  F --- Value of F for characteristic equation
    // ======================================================

    int j, ic = m / 2, l = 0, l0 = 0, j0 = 2;
    int jf = ic;
    double t0 = 0.0, t1 = 0.0, t2 = 0.0, b = a, f;

    if (kd == 1) {
        j0 = 3;
        l0 = 2;
    }
    if ((kd == 2) || (kd == 3)) { l = 1; }
    if (kd == 4) { jf = ic-1; }

    for (j = mj; j >= ic+1; j--) {
        t1 = -q*q/(pow(2.0*j + l, 2) - b + t1);
    }

    if (m <= 2) {
        if ((kd == 1) && (m == 0)) { t1 += t1; }
        if ((kd == 1) && (m == 2)) { t1 = -2.0*q*q/(4.0-b+t1) - 4.0; }
        if ((kd == 2) && (m == 1)) { t1 += q; }
        if ((kd == 3) && (m == 1)) { t1 -= q; }
    } else {
        if (kd == 1) { t0 = 4.0 - b + 2.0*q*q / b; }
        if (kd == 2) { t0 = 1.0 - b + q; }
        if (kd == 3) { t0 = 1.0 - b - q; }
        if (kd == 4) { t0 = 4.0 - b; }
        t2 = -q*q / t0;
        for (j = j0; j <= jf; j++) {
            t2 = -q*q/(pow(2.0*j -l-l0, 2.0) - b + t2);
        }
    }
    f = pow(2.0*ic+l, 2) + t1 + t2 - b;
    return f;
}


double specfun_cvql(int kd, int m, double q) {

    // ========================================================
    // Purpose: Compute the characteristic value of Mathieu
    //          functions  for q ≥ 3m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Initial characteristic value
    // ========================================================

    double a0, w, w2, w3, w4, w6, d1, d2, d3, d4, c1, p1, p2, cv1, cv2;

    w = 0.0;
    if ((kd == 1) || (kd == 2)) { w=2.0*m + 1.0; }
    if ((kd == 3) || (kd == 4)) { w=2.0*m - 1.0; }
    w2 = w*w;
    w3 = w*w2;
    w4 = w2*w2;
    w6 = w2*w4;
    d1 = 5.0+34.0/w2+9.0/w4;
    d2 = (33.0 + 410.0/w2 + 405.0/w4)/w;
    d3 = (63.0 + 1260.0/w2 + 2943.0/w4 + 486.0/w6)/w2;
    d4 = (527.0 + 15617.0/w2 + 69001.0/w4 + 41607.0/w6)/w3;
    c1 = 128.0;
    p2 = q/w4;
    p1 = sqrt(p2);
    cv1 = -2.0*q+2.0*w*sqrt(q) - (w2+1.0)/8.0;
    cv2 = (w+3.0/w) + d1/(32.0*p1) + d2/(8.0*c1*p2);
    cv2 = cv2 + d3/(64.0*c1*p1*p2)+d4/(16.0*c1*c1*p2*p2);
    a0 = cv1 - cv2/(c1*p1);
    return a0;
}


double specfun_cvqm(int m, double q) {

    // =====================================================
    // Purpose: Compute the characteristic value of Mathieu
    //          functions for q ≤ m*m
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    // Output:  A0 --- Initial characteristic value
    // =====================================================

    double hm1, hm3, hm5, a0;

    hm1= 0.5*q/(m*m-1.0);
    hm3=.25*pow(hm1, 3)/(m*m - 4.0);
    hm5 = hm1*hm3*q/((m*m - 1.0)*(m*m - 9.0));
    a0 = m*m + q*(hm1+(5.0*m*m + 7.0)*hm3 + (9.0*pow(m, 4) + 58.0*m*m + 29.0)*hm5);
    return a0;
}


void specfun_cy01(int kf, double complex z, double complex *zf, double complex *zd) {

    // ===========================================================
    // Purpose: Compute complex Bessel functions Y0(z), Y1(z)
    //          and their derivatives
    // Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
    //          KF --- Function choice code
    //              KF=0 for ZF=Y0(z) and ZD=Y0'(z)
    //              KF=1 for ZF=Y1(z) and ZD=Y1'(z)
    //              KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
    // Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
    //          ZD --- Y0'(z) or Y1'(z) or Y1''(z)
    // ===========================================================

    int k, k0;
    double a0, w0, w1;
    double complex cr, cp, cp0, cq0, cu, cp1, cq1, cbj0, cbj1,\
                   cby0, cby1, cdy0, cdy1, cs, ct1, ct2, z1, z2;

    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    const double rp2 = 2.0 / pi;
    const double complex ci = CMPLX(0.0, 1.0);

    static const double a[12] = {-0.703125e-01,        0.112152099609375, -0.5725014209747314,
                                  0.6074042001273483, -0.1100171402692467, 0.3038090510922384,
                                 -0.1188384262567832, 0.6252951493434797, -0.4259392165047669,
                                  0.3646840080706556, -0.3833534661393944, 0.4854014686852901};

    static const double b[12] = { 0.732421875e-01, -0.2271080017089844, 0.1727727502584457,
                                 -0.2438052969955606, 0.5513358961220206, -0.1825775547429318,
                                  0.8328593040162893, -0.5006958953198893, 0.3836255180230433,
                                 -0.3649010818849833, 0.4218971570284096, -0.5827244631566907};

    static const double a1[12] = { 0.1171875,          -0.144195556640625,   0.6765925884246826,
                                  -0.6883914268109947,  0.1215978918765359, -0.3302272294480852,
                                   0.1276412726461746, -0.6656367718817688,  0.4502786003050393,
                                  -0.3833857520742790,  0.4011838599133198, -0.5060568503314727};

    static const double b1[12] = {-0.1025390625,        0.2775764465332031, -0.1993531733751297,
                                   0.2724882731126854, -0.6038440767050702,  0.1971837591223663,
                                  -0.8902978767070678,  0.5310411010968522, -0.4043620325107754,
                                   0.3827011346598605, -0.4406481417852278,  0.6065091351222699};

    a0 = cabs(z);
    z1 = z;
    z2 = z * z;
    if (a0 == 0.0) {
        cbj0 = CMPLX(1.0, 0.0);
        cbj1 = CMPLX(0.0, 0.0);
        cby0 = CMPLX(-1e300, 0.0);
        cby1 = CMPLX(-1e300, 0.0);
        cdy0 = CMPLX( 1e300, 0.0);
        cdy1 = CMPLX( 1e300, 0.0);
        if (kf == 0) {
            *zf = cby0;
            *zd = cdy0;
        } else if (kf == 1) {
            *zf = cby1;
            *zd = cdy1;
        } else if (kf == 2) {
            *zf = cdy1;
            *zd = -cdy1 / z - (1.0 - 1.0 / (z * z)) * cby1;
        }
        return;
    }

    if (creal(z) < 0.0) {
        z1 = -z;
    }

    if (a0 <= 12.0) {
        cbj0 = CMPLX(1.0, 0.0);
        cr = CMPLX(1.0, 0.0);
        for (k = 1; k <= 40; k++) {
            cr = -0.25 * cr * z2 / (k * k);
            cbj0 += cr;
            if (cabs(cr) < cabs(cbj0) * 1.0e-15) break;
        }

        cbj1 = CMPLX(1.0, 0.0);
        cr = CMPLX(1.0, 0.0);
        for (k = 1; k <= 40; k++) {
            cr = -0.25 * cr * z2 / (k * (k + 1.0));
            cbj1 += cr;
            if (cabs(cr) < cabs(cbj1) * 1.0e-15) break;
        }

        cbj1 *= 0.5 * z1;
        w0 = 0.0;
        cr = CMPLX(1.0, 0.0);
        cs = CMPLX(0.0, 0.0);
        for (k = 1; k <= 40; k++) {
            w0 += 1.0 / k;
            cr = -0.25 * cr / (k * k) * z2;
            cp = cr * w0;
            cs += cp;
            if (cabs(cp) < cabs(cs) * 1.0e-15) break;
        }

        cby0 = rp2 * (clog(z1 / 2.0) + el) * cbj0 - rp2 * cs;
        w1 = 0.0;
        cr = 1.0;
        cs = 1.0;
        for (k = 1; k <= 40; k++) {
            w1 += 1.0 / k;
            cr = -0.25 * cr / (k * (k + 1)) * z2;
            cp = cr * (2.0 * w1 + 1.0 / (k + 1.0));
            cs += cp;
            if (cabs(cp) < cabs(cs) * 1.0e-15) break;
        }
        cby1 = rp2 * ((clog(z1 / 2.0) + el) * cbj1 - 1.0 / z1 - 0.25 * z1 * cs);
    } else {
        k0 = 12;
        if (a0 >= 35.0) k0 = 10;
        if (a0 >= 50.0) k0 = 8;

        ct1 = z1 - 0.25 * pi;
        cp0 = 1.0;
        for (k = 1; k <= k0; k++) {
            cp0 += a[k - 1] * pow(z1, -2 * k);
        }
        cq0 = -0.125 / z1;
        for (k = 1; k <= k0; k++)
            cq0 += b[k - 1] * pow(z1, -2 * k - 1);

        cu = csqrt(rp2 / z1);
        cbj0 = cu * (cp0 * cos(ct1) - cq0 * sin(ct1));
        cby0 = cu * (cp0 * sin(ct1) + cq0 * cos(ct1));

        ct2 = z1 - 0.75 * pi;
        cp1 = 1.0;
        for (k = 1; k <= k0; k++)
            cp1 += a1[k - 1] * pow(z1, -2 * k);

        cq1 = 0.375 / z1;
        for (k = 1; k <= k0; k++) {
            cq1 = cq1 + b1[k - 1] * pow(z1, -2 * k - 1);
        }
        cbj1 = cu * (cp1 * cos(ct2) - cq1 * sin(ct2));
        cby1 = cu * (cp1 * sin(ct2) + cq1 * cos(ct2));
    }

    if (creal(z) < 0.0) {
        if (cimag(z) < 0.0) cby0 = cby0 - 2.0 * ci * cbj0;
        if (cimag(z) > 0.0) cby0 = cby0 + 2.0 * ci * cbj0;
        if (cimag(z) < 0.0) cby1 = -(cby1 - 2.0 * ci * cbj1);
        if (cimag(z) > 0.0) cby1 = -(cby1 + 2.0 * ci * cbj1);
        cbj1 = -cbj1;
    }

    cdy0 = -cby1;
    cdy1 = cby0 - 1.0 / z * cby1;

    if (kf == 0) {
        *zf = cby0;
        *zd = cdy0;
    } else if (kf == 1) {
        *zf = cby1;
        *zd = cdy1;
    } else if (kf == 2) {
        *zf = cdy1;
        *zd = -cdy1 / z - (1.0 - 1.0 / (z * z)) * cby1;
    }
    return;
}


void specfun_cyzo(int nt, int kf, int kc, double complex *zo, double complex *zv) {

    // ===========================================================
    // Purpose : Compute the complex zeros of Y0(z), Y1(z) and
    //           Y1'(z), and their associated values at the zeros
    //           using the modified Newton's iteration method
    // Input:    NT --- Total number of zeros/roots
    //           KF --- Function choice code
    //                  KF=0 for  Y0(z) & Y1(z0)
    //                  KF=1 for  Y1(z) & Y0(z1)
    //                  KF=2 for  Y1'(z) & Y1(z1')
    //           KC --- Choice code
    //                  KC=0 for complex roots
    //                  KC=1 for real roots
    // Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
    //           ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
    //                     at the L-th zero
    // Routine called: CY01 for computing Y0(z) and Y1(z), and
    //                 their derivatives
    // ===========================================================

    int i, it, j, nr;
    double x, h, w, y, w0;
    double complex z, zf, zd, zfd, zgd, zp, zq, zw;

    x = 0.0;
    y = 0.0;
    h = 0.0;

    if (kc == 0) {
        x = -2.4;
        y = 0.54;
        h = 3.14;
    } else if (kc == 1) {
        x = 0.89;
        y = 0.0;
        h = -3.14;
    }

    if (kf == 1) {
        x = -0.503;
    }

    if (kf == 2) {
        x = 0.577;
    }
    z = CMPLX(x, y);
    w = 0.0;
    for (nr = 1; nr <= nt; nr++) {
        if (nr > 1) {
            z = zo[nr - 2] - h;
        }
        it = 0;
        do {
            it += 1;
            specfun_cy01(kf, z, &zf, &zd);
            zp = 1.0;
            for (i = 1; i < nr; i++) {
                zp *= (z - zo[i - 1]);
            }
            zfd = zf / zp;
            zq = 0.0;
            for (i = 1; i < nr; i++) {
                zw = 1.0;
                for (j = 1; j < nr; j++) {
                    if (j == i) { continue; }
                    zw *= (z - zo[j - 1]);
                }
                zq += zw;
            }
            zgd = (zd - zq * zfd) / zp;
            z -= zfd / zgd;
            w0 = w;
            w = cabs(z);
        } while ((it <= 50) && (fabs((w - w0) / w) > 1.0e-12));

        zo[nr - 1] = z;
    }

    for (i = 1; i <= nt; i++) {
        z = zo[i - 1];
        if ((kf == 0) || (kf == 2)) {
            specfun_cy01(1, z, &zf, &zd);
            zv[i - 1] = zf;
        } else if (kf == 1) {
            specfun_cy01(0, z, &zf, &zd);
            zv[i - 1] = zf;
        }
    }
    return;
}


double specfun_dvla(double x, double va) {

    // ====================================================
    // Purpose: Compute parabolic cylinder functions Dv(x)
    //          for large argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PD --- Dv(x)
    // Routines called:
    //       (1) VVLA for computing Vv(x) for large |x|
    //       (2) GAMMA2 for computing Г(x)
    // ====================================================

    int k;
    double ep, a0, gl, pd, r, vl, x1;

    const double pi=3.141592653589793;
    const double eps=1.0e-12;
    ep = exp(-.25*x*x);
    a0 = pow(fabs(x), va)*ep;
    r=1.0;
    pd=1.0;
    for (k = 1; k <= 16; k++) {
        r = -0.5*r*(2.0*k-va-1.0)*(2.0*k-va-2.0)/(k*x*x);
        pd += r;
        if (fabs(r/pd) < eps) { break; }
    }
    pd *= a0;
    if (x < 0.0) {
        x1 = -x;
        vl = specfun_vvla(x1, va);
        gl = specfun_gamma2(-va);
        pd = pi*vl/gl + cos(pi*va)*pd;
    }
    return pd;
}


double specfun_dvsa(double x, double va) {

    // ===================================================
    // Purpose: Compute parabolic cylinder function Dv(x)
    //          for small argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PD --- Dv(x)
    // Routine called: GAMMA2 for computing Г(x)
    // ===================================================

    int m;
    double ep, a0, va0, pd, ga0, g1, g0, r, vm, gm, r1, vt;

    const double pi = 3.141592653589793;
    const double eps = 1.0e-15;
    const double sq2 = sqrt(2);

    ep = exp(-0.25*x*x);
    va0 = 0.5*(1.0 - va);
    if (va == 0.0) {
        pd = ep;
    } else {
        if (x == 0.0) {
            if ((va0 <= 0.0) && (va0 == (int)va0)) {
                pd = 0.0;
            } else {
                ga0 = specfun_gamma2(va0);
                pd = sqrt(pi)/(pow(2.0, -0.5*va)*ga0);
            }
        } else {
            g1 = specfun_gamma2(-va);
            a0 = pow(2.0, -0.5*va - 1.0)*ep/g1;
            vt = -0.5*va;
            g0 = specfun_gamma2(vt);
            pd = g0;
            r = 1.0;
            for (m = 1; m <= 250; m++) {
                vm = 0.5*(m-va);
                gm = specfun_gamma2(vm);
                r = -r*sq2*x/m;
                r1 = gm*r;
                pd += r1;
                if (fabs(r1) < fabs(pd)*eps) { break; }
            }
            pd *= a0;
        }
    }
    return pd;
}


double specfun_e1xb(double x) {

    // ============================================
    // Purpose: Compute exponential integral E1(x)
    // Input :  x  --- Argument of E1(x)
    // Output:  E1 --- E1(x)  ( x > 0 )
    // ============================================

    int k, m;
    double e1, r, t, t0;
    const double ga = 0.5772156649015328;

    if (x == 0.0) {
        e1 = 1e300;
    } else if (x <= 1.0) {
        e1 = 1.0;
        r = 1.0;
        for (k = 1; k < 26; k++) {
            r = -r*k*x/pow(k+1.0, 2);
            e1 += r;
            if (fabs(r) <= fabs(e1)*1e-15) { break; }
        }
        e1 = -ga - log(x) + x*e1;
    } else {
        m = 20 + (int)(80.0/x);
        t0 = 0.0;
        for (k = m; k > 0; k--) {
            t0 = k / (1.0 + k / (x+t0));
        }
        t = 1.0 / (x + t0);
        e1 = exp(-x)*t;
    }
    return e1;
}


double complex specfun_e1z(double complex z) {

    // ====================================================
    // Purpose: Compute complex exponential integral E1(z)
    // Input :  z   --- Argument of E1(z)
    // Output:  CE1 --- E1(z)
    // ====================================================

    const double pi = 3.141592653589793;
    const double el = 0.5772156649015328;
    int k;
    double complex ce1, cr, zc, zd, zdc;
    double x = creal(z);
    double a0 = cabs(z);
    // Continued fraction converges slowly near negative real axis,
    // so use power series in a wedge around it until radius 40.0
    double xt = -2.0*fabs(cimag(z));

    if (a0 == 0.0) { return 1e300; }
    if ((a0 < 5.0) || ((x < xt) && (a0 < 40.0))) {
        // Power series
        ce1 = 1.0;
        cr = 1.0;
        for (k = 1; k < 501; k++) {
            cr = -cr*z*(k / pow(k + 1.0, 2));
            ce1 += cr;
            if (cabs(cr) < cabs(ce1)*1e-15) { break; }
        }
        if ((x <= 0.0) && (cimag(z) == 0.0)) {
            //Careful on the branch cut -- use the sign of the imaginary part
            // to get the right sign on the factor if pi.
            ce1 = -el - clog(-z) + z*ce1 - copysign(pi, cimag(z))*CMPLX(0.0, 1.0);
        } else {
            ce1 = -el - clog(z) + z*ce1;
        }
    } else {
        // Continued fraction https://dlmf.nist.gov/6.9
        //                  1     1     1     2     2     3     3
        // E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
        //                Z +   1 +   Z +   1 +   Z +   1 +   Z +
        zc = 0.0;
        zd = 1 / z;
        zdc = zd;
        zc += zdc;
        for (k = 1; k < 501; k++) {
            zd = 1.0 / (zd*k + 1.0);
            zdc *= (1.0*zd - 1.0);
            zc += zdc;

            zd = 1.0 / (zd*k + z);
            zdc *= (z*zd - 1.0);
            zc += zdc;
            if ((cabs(zdc) <= cabs(zc)*1e-15) && (k > 20)) { break; }
        }
        ce1 = cexp(-z)*zc;
        if ((x <= 0.0) && (cimag(z) == 0.0)) {
            ce1 -= pi*CMPLX(0.0, 1.0);
        }
    }
    return ce1;
}


double specfun_eix(double x) {

    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================

    const double ga = 0.5772156649015328;
    double ei, r;

    if (x == 0.0) {
        ei = -1.0e+300;
    } else if (x < 0) {
        ei = -specfun_e1xb(-x);
    } else if (fabs(x) <= 40.0) {
        // Power series around x=0
        ei = 1.0;
        r = 1.0;

        for (int k = 1; k <= 100; k++) {
            r = r * k * x / ((k + 1.0) * (k + 1.0));
            ei += r;
            if (fabs(r / ei) <= 1.0e-15) { break; }
        }
        ei = ga + log(x) + x * ei;
    } else {
        // Asymptotic expansion (the series is not convergent)
        ei = 1.0;
        r = 1.0;
        for (int k = 1; k <= 20; k++) {
            r = r * k / x;
            ei += r;
        }
        ei = exp(x) / x * ei;
    }
    return ei;
}


double complex specfun_eixz(double complex z) {

    // ============================================
    // Purpose: Compute exponential integral Ei(x)
    // Input :  x  --- Complex argument of Ei(x)
    // Output:  EI --- Ei(x)
    // ============================================

    double complex cei;
    const double pi = 3.141592653589793;
    cei = - specfun_e1z(-z);
    if (cimag(z) > 0.0) {
        cei += CMPLX(0.0, pi);
    } else if (cimag(z) < 0.0 ) {
        cei -= CMPLX(0.0, pi);
    } else {
        if (creal(z) > 0.0) {
            cei += CMPLX(0.0, copysign(pi, cimag(z)));
        }
    }
    return cei;
}


void specfun_eulerb(int n, double *en) {

    // ======================================
    // Purpose: Compute Euler number En
    // Input :  n --- Serial number
    // Output:  EN(n) --- En
    // ======================================

    int k, m, isgn;
    double r1, r2, s;
    const double hpi = 2.0 / 3.141592653589793;
    en[0] = 1.0;
    en[2] = -1.0;
    r1 = -4.0*pow(hpi, 3);
    for (m = 4; m <= n; m += 2) {
        r1 = -r1 * (m-1) * m * hpi * hpi;
        r2 = 1.0;
        isgn = 1;
        for (k = 3; k <= 1000; k += 2) {
            isgn = -isgn;
            s = pow(1.0 / k, m + 1);
            r2 += isgn * s;
            if (s < 1e-15) { break; }
        }
        en[m] = r1*r2;
    }
    return;
}


void specfun_fcoef(int kd, int m, double q, double a, double *fc) {

    // =====================================================
    // Purpose: Compute expansion coefficients for Mathieu
    //          functions and modified Mathieu functions
    // Input :  m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions
    //          KD --- Case code
    //                 KD=1 for cem(x,q)  ( m = 0,2,4,...)
    //                 KD=2 for cem(x,q)  ( m = 1,3,5,...)
    //                 KD=3 for sem(x,q)  ( m = 1,3,5,...)
    //                 KD=4 for sem(x,q)  ( m = 2,4,6,...)
    //          A  --- Characteristic value of Mathieu
    //                 functions for given m and q
    // Output:  FC(k) --- Expansion coefficients of Mathieu
    //                 functions ( k= 1,2,...,KM )
    //                 FC(1),FC(2),FC(3),... correspond to
    //                 A0,A2,A4,... for KD=1 case, A1,A3,
    //                 A5,... for KD=2 case, B1,B3,B5,...
    //                 for KD=3 case and B2,B4,B6,... for
    //                 KD=4 case
    // =====================================================

    int i, k, j, jm = 0, km, kb;
    double f1, fnan, qm, s, f, u, v, f2, f3, sp, ss, s0;

    for (i = 0; i < 251; ++i) { fc[i] = 0.0; }

    if (fabs(q) <= 1.0e-7) {
        // Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28)
        if (kd == 1) {
            jm = m / 2 + 1;
        } else if ((kd == 2) || (kd == 3)) {
            jm = (m - 1) / 2 + 1;
        } else if (kd == 4) {
            jm = m / 2;
        }

        if (jm + 1 > 251) {
            fnan = NAN;
            for (i = 0; i < 251; ++i) {
                fc[i] = fnan;
            }
            return;
        }
        // Proceed using the simplest expansion
        if (kd == 1 || kd == 2) {
            if (m == 0) {
                fc[0] = 1.0 / sqrt(2.0);
                fc[1] = -q / (2.0 * sqrt(2.0));
            } else if (m == 1) {
                fc[0] = 1.0;
                fc[1] = -q / 8.0;
            } else if (m == 2) {
                fc[0] = q / 4.0;
                fc[1] = 1.0;
                fc[2] = -q / 12.0;
            } else {
                fc[jm - 1] = 1.0;
                fc[jm] = -q / (4.0 * (m + 1));
                fc[jm - 2] = q / (4.0 * (m - 1));
            }
        } else if (kd == 3 || kd == 4) {
            if (m == 1) {
                fc[0] = 1.0;
                fc[1] = -q / 8.0;
            } else if (m == 2) {
                fc[0] = 1.0;
                fc[1] = -q / 12.0;
            } else {
                fc[jm - 1] = 1.0;
                fc[jm] = -q / (4.0 * (m + 1));
                fc[jm - 2] = q / (4.0 * (m - 1));
            }
        }
        return;
    } else if (q <= 1.0) {
        qm = 7.5 + 56.1 * sqrt(q) - 134.7 * q + 90.7 * sqrt(q) * q;
    } else {
        qm = 17.0 + 3.1 * sqrt(q) - 0.126 * q + 0.0037 * sqrt(q) * q;
    }

    km = (int)(qm + 0.5 * m);

    if (km > 251) {
        // Overflow, generate NaNs
        for (i = 0; i < 251; ++i) {
            fc[i] = NAN;
        }
        return;
    }

    kb = 0;
    s = 0.0;
    f = 1.0e-100;
    u = 0.0;
    fc[km - 1] = 0.0;
    f2 = 0.0;

    if (kd == 1) {
        for (k = km; k >= 3; k--) {
            v = u;
            u = f;
            f = (a - 4.0 * k * k) * u / q - v;

            if (fabs(f) < fabs(fc[k])) {
                kb = k;
                fc[0] = 1.0e-100;
                sp = 0.0;
                f3 = fc[k];
                fc[1] = a / q * fc[0];
                fc[2] = (a - 4.0) * fc[1] / q - 2.0 * fc[0];
                u = fc[1];
                f1 = fc[2];

                for (i = 3; i <= kb; i++) {
                    v = u;
                    u = f1;
                    f1 = (a - 4.0 * (i - 1.0) * (i - 1.0)) * u / q - v;
                    fc[i] = f1;

                    if (i == kb) { f2 = f1; }
                    if (i != kb) { sp += f1*f1; }
                }

                sp += 2.0*fc[0]*fc[0] + fc[1]*fc[1] + fc[2]*fc[2];
                ss = s + sp * (f3 / f2) * (f3 / f2);
                s0 = sqrt(1.0 / ss);

                for (j = 1; j <= km; j++) {
                    if (j <= kb + 1) {
                        fc[j - 1] = s0 * fc[j-1] * f3 / f2;
                    } else {
                        fc[j - 1] *= s0;
                    }
                }
                if (fc[0] < 0.0) { for (j = 0; j < km; j++) { fc[j] = -fc[j]; } }
                return;
            } else {
                fc[k - 1] = f;
                s += f*f;
            }
        }

        fc[1] = q * fc[2] / (a - 4.0 - 2.0 * q * q / a);
        fc[0] = q / a * fc[1];
        s += 2.0 * fc[0] * fc[0] + fc[1] * fc[1];
        s0 = sqrt(1.0 / s);

        for (k = 1; k <= km; k++) {
            fc[k - 1] *= s0;
        }
    } else if ((kd == 2) || (kd == 3)) {
        for (k = km; k >= 3; k--) {
            v = u;
            u = f;
            f = (a - (2.0 * k - 1) * (2.0 * k - 1)) * u / q - v;

            if (fabs(f) >= fabs(fc[k - 1])) {
                fc[k - 2] = f;
                s += f * f;
            } else {
                kb = k;
                f3 = fc[k - 1];
                goto L45;
            }
        }

        fc[0] = q / (a - 1.0 - pow(-1, kd) * q) * fc[1];
        s += fc[0] * fc[0];
        s0 = sqrt(1.0 / s);

        for (k = 1; k <= km; k++) {
            fc[k - 1] *= s0;
        }
        if (fc[0] < 0.0) { for (j = 0; j < km; j++) { fc[j] = -fc[j]; } }
        return;
L45:
        fc[0] = 1.0e-100;
        fc[1] = (a - 1.0 - pow(-1, kd) * q) / q * fc[0];
        sp = 0.0;
        u = fc[0];
        f1 = fc[1];

        for (i = 2; i <= kb - 1; i++) {
            v = u;
            u = f1;
            f1 = (a - (2.0 * i - 1) * (2.0 * i - 1)) * u / q - v;

            if (i != kb - 1) {
                fc[i] = f1;
                sp += f1 * f1;
            } else {
                f2 = f1;
            }
        }

        sp += fc[0] * fc[0] + fc[1] * fc[1];
        ss = s + sp * (f3 / f2) * (f3 / f2);
        s0 = sqrt(1.0 / ss);

        for (j = 1; j <= km; j++) {
            if (j < kb) {
                fc[j - 1] *= s0 * f3 / f2;
            }

            if (j >= kb) {
                fc[j - 1] *= s0;
            }
        }

    } else if (kd == 4) {
        for (k = km; k >= 3; k--) {
            v = u;
            u = f;
            f = (a - 4.0 * k * k) * u / q - v;

            if (fabs(f) >= fabs(fc[k])) {
                fc[k - 2] = f;
                s += f*f;
            } else {
                kb = k;
                f3 = fc[k - 1];
                goto L70;
            }
        }

        fc[0] = q / (a - 4.0) * fc[1];
        s += fc[0] * fc[0];
        s0 = sqrt(1.0 / s);

        for (k = 1; k <= km; k++) {
            fc[k - 1] *= s0;
        }
        if (fc[0] < 0.0) { for (j = 0; j < km; j++) { fc[j] = -fc[j]; } }
        return;
L70:
        fc[0] = 1.0e-100;
        fc[1] = (a - 4.0) / q * fc[0];
        sp = 0.0;
        u = fc[0];
        f1 = fc[1];

        for (i = 2; i <= kb - 1; i++) {
            v = u;
            u = f1;
            f1 = (a - 4.0 * i * i) * u / q - v;

            if (i != kb - 1) {
                fc[i] = f1;
                sp = sp + f1 * f1;
            } else {
                f2 = f1;
            }
        }

        sp += fc[0] * fc[0] + fc[1] * fc[1];
        ss = s + sp * (f3 / f2) * (f3 / f2);
        s0 = sqrt(1.0 / ss);

        for (j = 1; j <= km; j++) {
            if (j < kb) {
                fc[j - 1] *= s0 * f3 / f2;
            } else {
                fc[j - 1] *= s0;
            }
        }
    }
    if (fc[0] < 0.0) { for (j = 0; j < km; j++) { fc[j] = -fc[j]; } }
    return;
}


void specfun_fcszo(int kf, int nt, double complex *zo) {

    // ===============================================================
    // Purpose: Compute the complex zeros of Fresnel integral C(z)
    //          or S(z) using modified Newton's iteration method
    // Input :  KF  --- Function code
    //                  KF=1 for C(z) or KF=2 for S(z)
    //          NT  --- Total number of zeros
    // Output:  ZO(L) --- L-th zero of C(z) or S(z)
    // Routines called:
    //      (1) CFC for computing Fresnel integral C(z)
    //      (2) CFS for computing Fresnel integral S(z)
    // ==============================================================


    int i, j, it, nr;
    double psq, px, py, w, w0;
    double complex z, zp, zf, zd, zfd, zgd, zq, zw;
    const double pi = 3.141592653589793;
    psq = 0.0;
    w = 0.0;

    for (nr = 1; nr <= nt; ++nr) {
        if (kf == 1)
            psq = sqrt(4.0 * nr - 1.0);
        if (kf == 2)
            psq = 2.0 * sqrt(nr);

        px = psq - log(pi * psq) / (pi * pi * psq * psq * psq);
        py = log(pi * psq) / (pi * psq);
        z = CMPLX(px, py);

        if (kf == 2) {
            if (nr == 2) { z = CMPLX(2.8334, 0.2443); }
            if (nr == 3) { z = CMPLX(3.4674, 0.2185); }
            if (nr == 4) { z = CMPLX(4.0025, 0.2008); }
        }

        it = 0;
        do {
            it++;
            if (kf == 1) { specfun_cfc(z, &zf, &zd); }
            if (kf == 2) { specfun_cfs(z, &zf, &zd); }

            zp = 1.0;
            for (i = 1; i < nr; i++)
                zp *= (z - zo[i - 1]);

            zfd = zf / zp;
            zq = 0.0;
            for (i = 1; i < nr ; i++) {
                zw = 1.0;
                for (j = 1; j < nr; j++) {
                    if (j == i) { continue; }
                    zw *= (z - zo[j - 1]);
                }
                zq += zw;
            }
            zgd = (zd - zq * zfd) / zp;
            z -= zfd / zgd;
            w0 = w;
            w = cabs(z);
        } while ((it <= 50) && (fabs((w - w0) / w) > 1.0e-12));
        zo[nr - 1] = z;
    }
    return;
}


void specfun_ffk(int ks, double x, double *fr, double *fi, double *fm, double *fa,
         double *gr, double *gi, double *gm, double *ga) {

    // =======================================================
    // Purpose: Compute modified Fresnel integrals F±(x)
    //          and K±(x)
    // Input :  x   --- Argument of F±(x) and K±(x)
    //          KS  --- Sign code
    //                  KS=0 for calculating F+(x) and K+(x)
    //                  KS=1 for calculating F_(x) and K_(x)
    // Output:  FR  --- Re[F±(x)]
    //          FI  --- Im[F±(x)]
    //          FM  --- |F±(x)|
    //          FA  --- Arg[F±(x)]  (Degs.)
    //          GR  --- Re[K±(x)]
    //          GI  --- Im[K±(x)]
    //          GM  --- |K±(x)|
    //          GA  --- Arg[K±(x)]  (Degs.)
    // ======================================================

    const double srd = 57.29577951308233;
    const double eps = 1.0e-15;
    const double pi = 3.141592653589793;
    const double pp2 = 1.2533141373155;
    const double p2p = 0.7978845608028654;

    double fi0, c1, s1, cs, ss, xa, x2, x4, xc, xf, xf0, xf1, xg, xp, xq, xq2, xr, xs, xsu, xw;

    xa = fabs(x);
    x2 = x * x;
    x4 = x2 * x2;

    if (x == 0.0) {
        *fr = 0.5 * sqrt(0.5 * pi);
        *fi = pow(-1, ks) * (*fr);
        *fm = sqrt(0.25 * pi);
        *fa = pow(-1, ks) * 45.0;
        *gr = 0.5;
        *gi = 0.0;
        *gm = 0.5;
        *ga = 0.0;
    } else {
        if (xa <= 2.5) {
            xr = p2p * xa;
            c1 = xr;

            for (int k = 1; k <= 50; ++k) {
                xr = -0.5 * xr * (4.0 * k - 3.0) / k / (2.0 * k - 1.0) / (4.0 * k + 1.0) * x4;
                c1 += xr;
                if (fabs(xr / c1) < eps)
                    break;
            }

            s1 = p2p * xa * xa * xa / 3.0;
            xr = s1;

            for (int k = 1; k <= 50; ++k) {
                xr = -0.5 * xr * (4.0 * k - 1.0) / k / (2.0 * k + 1.0) / (4.0 * k + 3.0) * x4;
                s1 += xr;
                if (fabs(xr / s1) < eps)
                    break;
            }

            *fr = pp2 * (0.5 - c1);
            fi0 = pp2 * (0.5 - s1);
            *fi = pow(-1, ks) * fi0;
        } else if (xa < 5.5) {
            int m = (int)(42 + 1.75 * x2);
            xsu = 0.0;
            xc = 0.0;
            xs = 0.0;
            xf1 = 0.0;
            xf0 = 1.0e-100;

            for (int k = m; k >= 0; --k) {
                xf = (2.0 * k + 3.0) * xf0 / x2 - xf1;
                if (k % 2 == 0) {
                    xc += xf;
                } else {
                    xs += xf;
                }
                xsu += (2.0 * k + 1.0) * xf * xf;
                xf1 = xf0;
                xf0 = xf;
            }

            xq = sqrt(xsu);
            xw = p2p * xa / xq;
            c1 = xc * xw;
            s1 = xs * xw;
        } else {
            xr = 1.0;
            xf = 1.0;

            for (int k = 1; k <= 12; ++k) {
                xr = -0.25 * xr * (4.0 * k - 1.0) * (4.0 * k - 3.0) / x4;
                xf += xr;
            }

            xr = 1.0 / (2.0 * xa * xa);
            xg = xr;

            for (int k = 1; k <= 12; ++k) {
                xr = -0.25 * xr * (4.0 * k + 1.0) * (4.0 * k - 1.0) / x4;
                xg += xr;
            }

            c1 = 0.5 + (xf * sin(x2) - xg * cos(x2)) / sqrt(2.0 * pi) / xa;
            s1 = 0.5 - (xf * cos(x2) + xg * sin(x2)) / sqrt(2.0 * pi) / xa;
        }

        *fr = pp2 * (0.5 - c1);
        fi0 = pp2 * (0.5 - s1);
        *fi = pow(-1, ks) * fi0;
        *fm = cabs(CMPLX(*fr, *fi));

        if (*fr >= 0.0) {
            *fa = srd * atan((*fi) / (*fr));
        } else if (*fi > 0.0) {
            *fa = srd * (atan((*fi) / (*fr)) + pi);
        } else if (*fi < 0.0) {
            *fa = srd * (atan((*fi) / (*fr)) - pi);
        }

        xp = x2 + pi / 4.0;
        cs = cos(xp);
        ss = sin(xp);
        xq2 = 1.0 / sqrt(pi);

        *gr = xq2 * ((*fr) * cs + fi0 * ss);
        *gi = pow(-1, ks) * xq2 * (fi0 * cs - (*fr) * ss);
        *gm = sqrt((*gr) * (*gr) + (*gi) * (*gi));

        if (*gr >= 0.0) {
            *ga = srd * atan((*gi) / (*gr));
        } else if (*gi > 0.0) {
            *ga = srd * (atan((*gi) / (*gr)) + pi);
        } else if (*gi < 0.0) {
            *ga = srd * (atan((*gi) / (*gr)) - pi);
        }

        if (x < 0.0) {
            *fr = pp2 - (*fr);
            *fi = pow(-1, ks) * pp2 - (*fi);
            *fm = sqrt((*fr) * (*fr) + (*fi) * (*fi));
            *fa = srd * atan((*fi) / (*fr));
            *gr = cos(x2) - (*gr);
            *gi = -pow(-1, ks) * sin(x2) - (*gi);
            *gm = sqrt((*gr) * (*gr) + (*gi) * (*gi));
            *ga = srd * atan((*gi) / (*gr));
        }
    }
}


double specfun_gaih(double x) {

    // =====================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x), x = n/2, n=1,2,…
    // Output:  GA --- Г(x)
    // =====================================================

    int k, m;
    const double pi = 3.141592653589793;
    double ga = 0.0;

    if ((x == (int)x) && (x > 0.0)) {
        ga = 1.0;
        m = (int)(x - 1.0);
        for (k = 2; k < (m+1); k++) {
            ga *= k;
        }
    } else if (((x+0.5) == (int)(x+0.5)) && (x > 0.0)) {
        m = (int)x;
        ga = sqrt(pi);
        for (k = 1; k < (m+1); k++) {
            ga *= 0.5*(2.0*k - 1.0);
        }
    } else {
        ga = NAN;
    }
    return ga;
}


double specfun_gam0(double x) {

    // ================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
    // Output:  GA --- Г(x)
    // ================================================
    double gr;
    static const double g[25] = {
         1.0e0,
         0.5772156649015329e0, -0.6558780715202538e0, -0.420026350340952e-1, 0.1665386113822915e0,
        -0.421977345555443e-1, -0.96219715278770e-2, 0.72189432466630e-2, -0.11651675918591e-2,
        -0.2152416741149e-3, 0.1280502823882e-3, -0.201348547807e-4, -0.12504934821e-5,
         0.11330272320e-5, -0.2056338417e-6, 0.61160950e-8, 0.50020075e-8, -0.11812746e-8,
         0.1043427e-9, 0.77823e-11, -0.36968e-11, 0.51e-12, -0.206e-13, -0.54e-14, 0.14e-14
    };
    gr = g[24];
    for (int k = 23; k >= 0; k--) {
        gr = gr*x + g[k];
    }
    return 1.0 / (gr * x);
}


double specfun_gamma2(double x) {

    // ==================================================
    // Purpose: Compute gamma function Г(x)
    // Input :  x  --- Argument of Г(x)
    //                 ( x is not equal to 0,-1,-2,…)
    // Output:  GA --- Г(x)
    // ==================================================

    double ga, gr, r, z;
    int k, m;
    const double pi = 3.141592653589793;
    static const double g[26] = {
        1.0000000000000000e+00,  0.5772156649015329e+00, -0.6558780715202538e+00, -0.4200263503409520e-01,
        0.1665386113822915e+00, -0.4219773455554430e-01, -0.9621971527877000e-02,  0.7218943246663000e-02,
       -0.1165167591859100e-02, -0.2152416741149000e-03,  0.1280502823882000e-03, -0.2013485478070000e-04,
       -0.1250493482100000e-05,  0.1133027232000000e-05, -0.2056338417000000e-06,  0.6116095000000000e-08,
        0.5002007500000000e-08, -0.1181274600000000e-08,  0.1043427000000000e-09,  0.7782300000000000e-11,
       -0.3696800000000000e-11,  0.5100000000000000e-12, -0.2060000000000000e-13, -0.5400000000000000e-14,
        0.1400000000000000e-14,  0.1000000000000000e-15
    };
    if (x == (int)x) {
        if (x > 0.0) {
            ga = 1.0;
            m = (int)(x - 1);
            for  (k = 2; k < (m+1); k++) {
                ga *= k;
            }
        } else {
            ga = 1e300;
        }
    } else {
        r = 1.0;
        if (fabs(x) > 1.0) {
            z = fabs(x);
            m = (int)z;
            for  (k = 1; k < (m+1); k++) {
                r *= (z-k);
            }
            z -= m;
        } else {
            z = x;
        }
        gr = g[25];
        for ( k = 25; k > 0; k--) {
            gr *= z;
            gr += g[k-1];
        }
        ga = 1.0 / (gr*z);
        if (fabs(x) > 1.0) {
            ga *= r;
            if (x < 0.0) {
                ga = -pi / (x*ga*sin(pi*x));
            }
        }
    }
    return ga;
}


void specfun_gmn(int m, int n, double c, double x, double *bk, double *gf, double *gd) {

    // ===========================================================
    // Purpose: Compute gmn(-ic,ix) and its derivative for oblate
    //          radial functions with a small argument
    // ===========================================================

    int ip, k, nm;
    double xm, gf0, gw, gd0, gd1;
    const double eps = 1.0e-14;
    ip = ((n - m) % 2 == 0 ? 0 : 1);
    nm = 25 + (int)(0.5 * (n - m) + c);
    xm = pow(1.0 + x * x, -0.5 * m);
    gf0 = 0.0;
    gw = 0.0;

    for (k = 1; k <= nm; k++) {
        gf0 += bk[k - 1] * pow(x, 2.0 * k - 2.0);
        if ((fabs((gf0 - gw) / gf0) < eps) && (k >= 10)) { break; }
        gw = gf0;
    }

    *gf = xm * gf0 * pow(x, 1 - ip);
    gd1 = -m * x / (1.0 + x * x) * (*gf);
    gd0 = 0.0;

    for (k = 1; k < nm; ++k) {
        if (ip == 0) {
            gd0 += (2.0 * k - 1.0) * bk[k - 1] * pow(x, 2.0 * k - 2.0);
        } else {
            gd0 += 2.0 * k * bk[k - 1] * pow(x, 2.0 * k - 1.0);
        }
        if ((fabs((gd0 - gw) / gd0) < eps) && (k >= 10)) { break; }
        gw = gd0;
    }
    *gd = gd1 + xm * gd0;
}


double complex specfun_hygfz(double a, double b, double c, double complex z, int *isfer) {

    // ======================================================
    // Purpose: Compute the hypergeometric function for a
    //          complex argument, F(a,b,c,z)
    // Input :  a --- Parameter
    //          b --- Parameter
    //          c --- Parameter,  c <> 0,-1,-2,...
    //          z --- Complex argument
    // Output:  ZHF --- F(a,b,c,z)
    //          ISFER --- Error flag
    // Routines called:
    //      (1) GAMMA2 for computing gamma function
    //      (2) PSI_SPEC for computing psi function
    // ======================================================

    int L0 = 0, L1 = 0, L2 = 0, L3 = 0, L4 = 0, L5 = 0, L6 = 0;
    int j, k=1, m, mab, mcab, nca, ncb, nm;
    double a0, aa, bb, ca, cb, g0, g1, g2, g3, ga, gab, gam, gabc, gb, gba, gbm, gc, gcab,\
           gca, gcb, gm, pa, pac, pb, pca, rk1, rk2, rm, sp0, sm, sp, sq, sj1, sj2, w0, ws;
    double complex z00, z1, zc0, zc1, zf0, zf1, zhf = 0.0, zp, zr, zp0, zr0, zr1, zw = 0.0;
    double x = creal(z);
    double y = cimag(z);
    double eps = 1e-15;
    double pi = 3.141592653589793;
    double el = 0.5772156649015329;
    *isfer = 0;

    if ((c == (int)c) && (c < 0.0)) { L0 = 1; }
    if ((fabs(1 - x) < eps) && (y == 0.0) && (c-a-b <= 0.0)) { L1 = 1; }
    if ((cabs(z+1) < eps) && (fabs(c-a+b - 1.0) < eps)) { L2 = 1; }
    if ((a == (int)a) && (a < 0.0)) { L3 = 1; }
    if ((b == (int)b) && (b < 0.0)) { L4 = 1; }
    if (((c-a) == (int)(c-a)) && (c-a <= 0.0)) { L5 = 1; }
    if (((c-b) == (int)(c-b)) && (c-b <= 0.0)) { L6 = 1; }
    aa = a;
    bb = b;
    a0 = cabs(z);
    if (a0 > 0.95) { eps = 1e-8; }
    if (L0 || L1) {
        *isfer = 3;
        return 0.0;
    }

    if ((a0 == 0.0) || (a == 0.0) || (b == 0.0)) {
        zhf = 1.0;
    } else if ((z == 1.0) && (c-a-b > 0.0)) {
        gc = specfun_gamma2(c);
        gcab = specfun_gamma2(c-a-b);
        gca = specfun_gamma2(c-a);
        gcb = specfun_gamma2(c-b);
        zhf = gc*gcab/(gca*gcb);
    } else if (L2) {
        g0 = sqrt(pi)*pow(2.0, -a);
        g1 = specfun_gamma2(c);
        g2 = specfun_gamma2(1.0 + 0.5*a - b);
        g3 = specfun_gamma2(0.5 + 0.5*a);
        zhf = g0*g1/(g2*g3);
    } else if (L3 || L4) {
        if (L3) { nm = (int)fabs(a); }
        if (L4) { nm = (int)fabs(b); }
        zhf = 1.0;
        zr = 1.0;
        for (k = 1; k < (nm+1); k++) {
            zr = zr*(c-a+k-1.0)*(c-b+k-1.0)/(k*(c+k-1.0))*z;
            zhf += zr;
        }
    } else if (L5 || L6) {
        if (L5) { nm = (int)fabs(c-a); }
        if (L6) { nm = (int)fabs(c-b); }
        zhf = 1.0;
        zr = 1.0;
        for (k = 1; k < (nm+1); k++) {
            zr = zr*(c-a+k-1.0)*(c-b+k-1.0)/(k*(c+k-1.0))*z;
            zhf += zr;
        }
        zhf *= cpow(1.0-z, c-a-b);
    } else if (a0 <= 1.0) {
        if (x < 0.0) {
            z1 = z / (z - 1.0);
            if ((c > a) && (b < a) && (b > 0.0)) {
                a = aa;
                b = bb;
            }
            zc0 = 1.0 / cpow(1.0 - z, a);
            zhf = 1.0;
            zr0 = 1.0;
            zw = 0.0;
            for (k = 1; k <501; k++) {
                zr0 = zr0*(a+k-1.0)*(c-b+k-1.0)/(k*(c+k-1.0))*z1;
                zhf += zr0;
                if (cabs(zhf-zw) < cabs(zhf)*eps) { break; }
                zw = zhf;
            }
            zhf *= zc0;
        } else if (a0 >= 0.9) {
            gm = 0.0;
            mcab = (int)(c-a-b + eps*copysign(1.0, c-a-b));
            if (fabs(c-a-b-mcab) < eps) {
                m = (int)(c-a-b);
                ga = specfun_gamma2(a);
                gb = specfun_gamma2(b);
                gc = specfun_gamma2(c);
                gam = specfun_gamma2(a+m);
                gbm = specfun_gamma2(b+m);
                pa = specfun_psi_spec(a);
                pb = specfun_psi_spec(b);
                if (m != 0) { gm = 1.0; }
                for (j = 1; j < abs(m); j++) {
                    gm *= j;
                }
                rm = 1.0;
                for (j = 1; j < abs(m)+1; j++) {
                    rm *= j;
                }
                zf0 = 1.0;
                zr0 = 1.0;
                zr1 = 1.0;
                sp0 = 0.0;
                sp = 0.0;
                if (m >= 0) {
                    zc0 = gm*gc/(gam*gbm);
                    zc1 = -gc*cpow(z-1.0, m)/(ga*gb*rm);
                    for (k = 1; k < m; k++) {
                        zr0 = zr0*(a+k-1.0)*(b+k-1.0)/(k*(k-m))*(1.0-z);
                        zf0 += zr0;
                    }
                    for (k = 1; k < (m+1); k++) {
                        sp0 += 1.0/(a+k-1.0) + 1.0/(b+k-1.0) - 1.0/k;
                    }
                    zf1 = pa + pb + sp0 + 2.0*el + clog(1.0 - z);
                    zw = 0.0;
                    for (k = 1; k <501; k++) {
                        sp += (1.0-a)/(k*(a+k-1.0)) + (1.0-b)/(k*(b+k-1.0));
                        sm = 0.0;
                        for (j = 1; j < (m+1); j++) {
                            sm += (1.0-a)/((j+k)*(a+j+k-1.0)) + 1.0/(b+j+k-1.0);
                        }
                        zp = pa + pb + 2.0*el + sp + sm + clog(1.0 - z);
                        zr1 = zr1*(a+m+k-1.0)*(b+m+k-1.0) / (k*(m+k))*(1.0-z);
                        zf1 += zr1*zp;
                        if (cabs(zf1-zw) < cabs(zf1)*eps) { break; }
                        zw = zf1;
                    }
                    zhf = zf0*zc0 + zf1*zc1;
                } else if (m < 0) {
                    m = -m;
                    zc0 = gm*gc/(ga*gb*cpow(1.0 - z, m));
                    zc1 = -(pow(-1.0, m))*gc/(gam*gbm*rm);
                    for (k = 1; k < m; k++) {
                        zr0 = zr0*(a-m+k-1.0)*(b-m+k-1.0)/(k*(k-m))*(1.0-z);
                        zf0 += zr0;
                    }
                    for (k = 1; k < (m+1); k++) {
                        sp0 += 1.0 / k;
                    }
                    zf1 = pa + pb -sp0 + 2.0*el + clog(1.0 - z);
                    zw = 0.0;
                    for (k = 1; k <501; k++) {
                        sp += (1.0-a)/(k*(a+k-1.0)) + (1.0-b)/(k*(b+k-1.0));
                        sm = 0.0;
                        for (j = 1; j < (m+1); j++) {
                            sm += 1.0/(j+k);
                        }
                        zp = pa + pb+2.0*el + sp - sm + clog(1.0 -z );
                        zr1 = zr1*(a+k-1.0)*(b+k-1.0)/(k*(m+k))*(1.0-z);
                        zf1 += zr1*zp;
                        if (cabs(zf1-zw) < cabs(zf1)*eps) { break; }
                        zw = zf1;
                    }
                    zhf = zf0*zc0 + zf1*zc1;
                }
            } else {
                ga = specfun_gamma2(a);
                gb = specfun_gamma2(b);
                gc = specfun_gamma2(c);
                gca = specfun_gamma2(c-a);
                gcb = specfun_gamma2(c-b);
                gcab = specfun_gamma2(c-a-b);
                gabc = specfun_gamma2(a+b-c);
                zc0 = gc*gcab/(gca*gcb);
                zc1 = gc*gabc/(ga*gb)*cpow(1.0-z, c-a-b);
                zhf = 0.0;
                zr0 = zc0;
                zr1 = zc1;
                zw = 0.0;
                for (k = 1; k < 501; k++) {
                    zr0 = zr0*(a+k-1.0)*(b+k-1.0)/(k*(a+b-c+k))*(1.0-z);
                    zr1 = zr1*(c-a+k-1.0)*(c-b+k-1.0)/(k*(c-a-b+k))*(1.0-z);
                    zhf += zr0+zr1;
                    if (cabs(zhf-zw) < cabs(zhf)*eps) { break; }
                    zw = zhf;
                }
                zhf += zc0 + zc1;
            }
        } else {
            z00 = 1.0;
            if ((c-a < a) && (c-b < b)) {
                z00 = cpow(1.0 - z, c-a-b);
                a = c-a;
                b = c-b;
            }
            zhf = 1.0;
            zr = 1.0;
            zw = 0.0;
            for (k = 1; k < 1501; k++) {
                zr = zr*(a+k-1.0)*(b+k-1.0)/(k*(c+k-1.0))*z;
                zhf += zr;
                if (cabs(zhf-zw) < cabs(zhf)*eps) { break; }
                zw = zhf;
            }
            zhf *= z00;
        }
    } else if (a0 > 1.0) {
        mab = (int)(a - b + eps*copysign(1.0, a - b));
        if ((fabs(a-b-mab) < eps) && (a0 <= 1.1)) { b += eps; }
        if (fabs(a-b-mab) > eps) {
            ga = specfun_gamma2(a);
            gb = specfun_gamma2(b);
            gc = specfun_gamma2(c);
            gab = specfun_gamma2(a-b);
            gba = specfun_gamma2(b-a);
            gca = specfun_gamma2(c-a);
            gcb = specfun_gamma2(c-b);
            zc0 = gc*gba/(gca*gb*cpow(-z, a));
            zc1 = gc*gab/(gcb*ga*cpow(-z, b));
            zr0 = zc0;
            zr1 = zc1;
            zhf = 0.0;
            for (k = 1; k < 501; k++) {
                zr0 = zr0*(a+k-1.0)*(a-c+k)/((a-b+k)*k*z);
                zr1 = zr1*(b+k-1.0)*(b-c+k)/((b-a+k)*k*z);
                zhf += zr0+zr1;
                if (cabs(zhf-zw) < cabs(zhf)*eps) { break; }
                zw = zhf;
            }
            zhf += zc0 + zc1;
        } else {
            if (a-b < 0.0) {
                a = bb;
                b = aa;
            }
            ca = c - a;
            cb = c - b;
            nca = (int)(ca + eps*copysign(1.0, ca));
            ncb = (int)(cb + eps*copysign(1.0, cb));
            if ((fabs(ca-nca) < eps) || (fabs(cb-ncb) < eps)) { c += eps; }
            ga = specfun_gamma2(a);
            gc = specfun_gamma2(c);
            gcb = specfun_gamma2(c-b);
            pa = specfun_psi_spec(a);
            pca = specfun_psi_spec(c-a);
            pac = specfun_psi_spec(a-c);
            mab = (int)(a-b+eps);
            zc0 = gc / (ga*cpow(-z, b));
            gm = specfun_gamma2(a-b);
            zf0 = gm/gcb*zc0;
            zr = zc0;
            for (k = 1; k < mab; k++) {
                zr = zr*(b+k-1.0)/(k*z);
                g0 = specfun_gamma2(a-b-k);
                zf0 += zr*g0/specfun_gamma2(c-b-k);
            }
            if (mab == 0) { zf0 = 0.0; }
            zc1 = gc/(ga*gcb*cpow(-z, a));
            sp = -2.0*el - pa- pca;
            for (j = 1; j < (mab+1); j++) {
                sp += 1.0 / j;
            }
            zp0 = sp + clog(-z);
            sq = 1.0;
            for (j = 1; j < (mab+1); j++) {
                sq = sq * (b+j-1.0)*(b-c+j)/j;
            }
            zf1 = (sq*zp0)*zc1;
            zr = zc1;
            rk1 = 1.0;
            sj1 = 0.0;
            w0 = 0.0;
            for (k = 1; k < 10001; k++) {
                zr /= z;
                rk1 = rk1*(b+k-1.0)*(b-c+k)/(k*k);
                rk2 = rk1;
                for (j = k+1; j <= (k+mab); j++) {
                    rk2 = rk2 * (b+j-1.0)*(b-c+j)/j;
                }
                sj1 += (a-1.0)/(k*(a+k-1.0)) + (a-c-1.0)/(k*(a-c+k-1.0));
                sj2 = sj1;
                for (j = k+1; j <= (k+mab); j++) {
                    sj2 += 1.0 / j;
                }
                zp= -2.0*el -pa - pac + sj2 - 1.0/(k+a-c) - pi/tan(pi*(k+a-c)) + clog(-z);
                zf1 += rk2*zr*zp;
                ws = cabs(zf1);
                if (fabs((ws-w0)/ws) < eps) { break; }
                w0 = ws;
            }
            zhf = zf0 + zf1;
        }
    }
    a = aa;
    b = bb;
    if (k > 150) { *isfer = 5; }
    return zhf;
}


void specfun_itairy(double x, double *apt, double *bpt, double *ant, double *bnt) {

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
    double fx, gx, r, su1, su2, su3, su4, su5, su6, xp6, xe, xr1, xr2;

    const double pi = 3.141592653589793;
    const double c1 = 0.355028053887817;
    const double c2 = 0.258819403792807;
    const double sr3 = 1.732050807568877;
    const double q0 = 1.0 / 3.0;
    const double q1 = 2.0 / 3.0;
    const double q2 = 1.4142135623730951;
    const double eps = 1e-5;
    static const double a[16] = {
        0.569444444444444    , 0.891300154320988    , 0.226624344493027e+01, 0.798950124766861e+01,
        0.360688546785343e+02, 0.198670292131169e+03, 0.129223456582211e+04, 0.969483869669600e+04,
        0.824184704952483e+05, 0.783031092490225e+06, 0.822210493622814e+07, 0.945557399360556e+08,
        0.118195595640730e+10, 0.159564653040121e+11, 0.231369166433050e+12, 0.358622522796969e+13
    };

    if (x == 0.0) {
        *apt = 0.0;
        *bpt = 0.0;
        *ant = 0.0;
        *bnt = 0.0;
    } else {
        if (fabs(x) <= 9.25) {
            for (l = 0; l < 2; l++) {
                x = x * pow(-1, l);
                fx = x;
                r = x;
                for (k = 1; k < 41; k++) {
                    r = r*(3.0*k - 2.0)/(3.0*k + 1.0)*x/(3.0*k)*x/(3.0*k - 1.0)*x;
                    fx += r;
                    if (fabs(r) < fabs(fx)*eps) { break; }
                }
                gx = 0.5*x*x;
                r = gx;
                for (k = 1; k < 41; k++) {
                    r = r*(3.0*k - 1.0)/(3.0*k + 2.0)*x/(3.0*k)*x/(3.0*k + 1.0)*x;
                    gx += r;
                    if (fabs(r) < fabs(gx)*eps) { break; }
                }
                *ant = c1*fx - c2*gx;
                *bnt = sr3*(c1*fx + c2*gx);
                if (l == 0) {
                    *apt = *ant;
                    *bpt = *bnt;
                } else {
                    *ant = -*ant;
                    *bnt = -*bnt;
                    x= -x;
                }
            }
        } else {
            xe = x*sqrt(x)/1.5;
            xp6 = 1.0/sqrt(6.0*pi*xe);
            su1 = 1.0;
            r = 1.0;
            xr1 = 1.0/xe;
            for (k = 1; k < 17; k++) {
                r = -r * xr1;
                su1 += a[k-1]*r;
            }
            su2 = 1.0;
            r = 1.0;
            for (k = 1; k < 17; k++) {
                r = r * xr1;
                su2 += a[k-1]*r;
            }
            *apt = q0 - exp(-xe)*xp6*su1;
            *bpt = 2.0*exp(xe)*xp6*su2;
            su3 = 1.0;
            r = 1.0;
            xr2 = 1.0/(xe*xe);
            for (k = 1; k < 9; k++) {
                r = -r * xr2;
                su3 += a[2*k - 1]*r;
            }
            su4 = a[0]*xr1;
            r = xr1;
            for (k = 1; k < 8; k++) {
                r = -r * xr2;
                su4 += a[2*k]*r;
            }
            su5 = su3 + su4;
            su6 = su3 - su4;
            *ant = q1 - q2 * xp6 * (su5*cos(xe) - su6*sin(xe));
            *bnt = q2 * xp6 * (su5*sin(xe) + su6*cos(xe));
        }
    }
    return;
}


void specfun_itika(double x, double *ti, double *tk) {

    // =======================================================
    // Purpose: Integrate modified Bessel functions I0(t) and
    //          K0(t) with respect to t from 0 to x
    // Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
    // Output:  TI --- Integration of I0(t) from 0 to x
    //          TK --- Integration of K0(t) from 0 to x
    // =======================================================

    int k;
    double rc1, rc2, e0, b1, b2, rs, r, tw, x2;
    static const double a[10] = {
        0.625, 1.0078125, 2.5927734375, 9.1868591308594,
        4.1567974090576e+1, 2.2919635891914e+2,1.491504060477e+3,
        1.1192354495579e+4, 9.515939374212e+4,9.0412425769041e+5
    };
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;

    if (x == 0.0) {
        *ti = 0.0;
        *tk = 0.0;
        return;
    } else if (x < 20.0) {
        x2 = x*x;
        *ti = 1.0;
        r = 1.0;
        for (k = 1; k <= 50; k++) {
            r = 0.25 * r *(2 * k - 1.0) / (2 * k + 1.0) / (k * k) * x2;
            *ti += r;
            if (fabs(r/(*ti)) < 1.0e-12) { break; }
        }
        *ti *= x;
    } else {
        x2 = 0.0;
        *ti = 1.0;
        r = 1.0;
        for (k = 1; k <= 10; k++) {
            r = r/x;
            *ti += a[k-1]*r;
        }
        rc1 = 1.0 / sqrt(2.0*pi*x);
        *ti = rc1 * exp(x) * (*ti);
    }

    if (x < 12.0) {
        e0 = el + log(x / 2.0);
        b1 = 1.0 - e0;
        b2 = 0.0;
        rs = 0.0;
        r = 1.0;
        tw = 0.0;
        for (k = 1; k <= 50; k++) {
            r = 0.25*r*(2*k - 1.0) / (2*k + 1.0) / (k*k)*x2;
            b1 += r * (1.0 / (2 * k + 1) - e0);
            rs += 1.0 / k;
            b2 += r * rs;
            *tk = b1 + b2;
            if (fabs((*tk - tw) / (*tk)) < 1.0e-12) { break; }
            tw = *tk;
        }
        *tk *= x;
    } else {
        *tk = 1.0;
        r = 1.0;
        for (k = 1; k <= 10; k++) {
            r = -r / x;
            *tk = *tk + a[k-1]*r;
        }
        rc2 = sqrt(pi/(2.0*x));
        *tk = pi/2.0 - rc2*(*tk)*exp(-x);
    }
    return;
}


void specfun_itjya(double x, double *tj, double *ty) {


    int k;
    double a[18], a0, a1, af, bf, bg, r, r2, rc, rs, ty1, ty2, x2, xp;
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    const double eps = 1.0e-12;

    if (x == 0.0) {
        *tj = 0.0;
        *ty = 0.0;
    } else if (x <= 20.0) {
        x2 = x*x;
        *tj = x;
        r = x;
        for (k = 1; k < 61; k++) {
            r = -0.25*r*(2.0*k - 1.0) / (2.0*k + 1.0) / (k*k)*x2;
            *tj += r;
            if (fabs(r) < fabs(*tj)*eps) { break; }
        }
        ty1 = (el + log(x/2.0)) * (*tj);
        rs = 0.0;
        ty2 = 1.0;
        r = 1.0;
        for (k = 1; k < 61; k++) {
            r = -0.25*r*(2.0*k - 1.0) / (2.0*k + 1.0) / (k*k) * x2;
            rs += 1.0 / k;
            r2 = r * (rs + 1.0 / (2.0*k + 1.0));
            ty2 += r2;
            if (fabs(r2) < fabs(ty2)*eps) { break; }
        }
        *ty = (ty1 - x*ty2)*2.0/pi;
    } else {
        a0 = 1.0;
        a1 = 5.0 / 8.0;
        a[0] = a1;

        for (int k = 1; k <= 16; k++) {
            af = ((1.5*(k+0.5)*(k+5.0/6.0)*a1-0.5*(k+0.5)*(k+0.5)*(k-0.5)*a0))/(k+1.0);
            a[k] = af;
            a0 = a1;
            a1 = af;
        }
        bf = 1.0;
        r = 1.0;
        for (int k = 1; k <= 8; k++) {
            r = -r / (x*x);
            bf += a[2*k-1] * r;
        }
        bg = a[0] / x;
        r = 1.0 / x;
        for (int k = 1; k <= 8; k++) {
            r = -1.0/(x*x);
            bg += a[2*k] * r;
        }
        xp = x + 0.25 * pi;
        rc = sqrt(2.0 / (pi * x));
        *tj = 1.0 - rc*(bf*cos(xp) + bg*sin(xp));
        *ty = rc*(bg*cos(xp) - bf*sin(xp));
    }
    return;
}


double specfun_itsh0(double x) {

    // ===================================================
    // Purpose: Evaluate the integral of Struve function
    //          H0(t) with respect to t from 0 and x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TH0 --- Integration of H0(t) from 0 and x
    // ===================================================

    int k;
    double a[25], a0, a1, af, bf, bg, r, rd, s, s0, th0, ty, xp;
    const double pi = 3.141592653589793;
    const double el = 0.57721566490153;

    r = 1.0;
    if (x <= 30.0) {
        s = 0.5;
        for (k = 1; k < 101; k++) {
            rd = 1.0;
            if (k == 1) { rd = 0.5; }
            r = -r*rd*k/(k+1.0)*pow(x/(2.0*k+1.0), 2);
            s += r;
            if (fabs(r) < fabs(s)*1.0e-12) { break; }
        }
        th0 = 2.0 / pi * x * x * s;
    } else {
        s = 1.0;
        for (k = 1; k < 13; k++) {
            r = -r*k/(k+1.0)*pow((2.0*k+1.0)/x, 2);
            s += r;
            if (fabs(r) < fabs(s)*1.0e-12) { break; }
        }
        s0 = s / (pi*x*x) + 2.0 / pi *(log(2.0*x) + el);
        a0 = 1.0;
        a1 = 5.0 / 8.0;
        a[0] = a1;
        for (k = 1; k < 21; k++) {
            af=((1.5*(k+0.5)*(k+5.0/6.0)*a1-0.5*(k+0.5)*(k+0.5)*(k-0.5)*a0))/(k+1.0);
            a[k] = af;
            a0 = a1;
            a1 = af;
        }
        bf = 1.0;
        r = 1.0;
        for (k = 1; k < 11; k++) {
            r = -r / (x*x);
            bf += a[2*k - 1]*r;
        }
        bg = a[0]*x;
        r = 1.0 / x;
        for (k = 1; k < 10; k++) {
            r = -r / (x*x);
            bg += a[2*k]*r;
        }
        xp = x + 0.25*pi;
        ty = sqrt(2.0/(pi*x))*(bg*cos(xp) - bf*sin(xp));
        th0 = ty + s0;
    }
    return th0;
}


double specfun_itsl0(double x) {

    // ===========================================================
    // Purpose: Evaluate the integral of modified Struve function
    //          L0(t) with respect to t from 0 to x
    // Input :  x   --- Upper limit  ( x ≥ 0 )
    // Output:  TL0 --- Integration of L0(t) from 0 to x
    // ===========================================================

    int k;
    double a[18], a0, a1, af, r, rd, s, s0, ti, tl0;
    const double pi = 3.141592653589793;
    const double el = 0.57721566490153;
    r = 1.0;
    if (x <= 20.0) {
        s=0.5;
        for (k = 1; k < 101; k++) {
            rd = 1.0;
            if (k == 1) { rd=0.5; }
            r = r*rd*k/(k+1.0)*pow(x/(2.0*k+1.0), 2);
            s += r;
            if (fabs(r/s) < 1.0e-12) { break; }
        }
        tl0=2.0/pi*x*x*s;
    } else {
        s = 1.0;
        for (k = 1; k < 11; k++) {
            r = r*k/(k + 1.0)*pow((2.0*k + 1.0)/x, 2);
            s += r;
            if (fabs(r/s) < 1.0e-12) { break; }
        }
        s0 = -s/(pi*x*x)+2.0/pi*(log(2.0*x)+el);
        a0 = 1.0;
        a1 = 5.0/8.0;
        a[0] = a1;
        for (k = 1; k < 11; k++) {
            af=((1.5*(k+.50)*(k+5.0/6.0)*a1-0.5*pow(k+0.5, 2)*(k-0.5)*a0))/(k+1.0);
            a[k] = af;
            a0 = a1;
            a1 = af;
        }
        ti = 1.0;
        r = 1.0;
        for (k = 1; k < 11; k++) {
            r = r / x;
            ti += a[k-1]*r;
        }
        tl0 = ti/sqrt(2*pi*x)*exp(x) + s0;
    }
    return tl0;
}


double specfun_itth0(double x) {

    // ===========================================================
    // Purpose: Evaluate the integral H0(t)/t with respect to t
    //          from x to infinity
    // Input :  x   --- Lower limit  ( x ≥ 0 )
    // Output:  TTH --- Integration of H0(t)/t from x to infinity
    // ===========================================================

    int k;
    double f0, g0, r, s, t, tth, tty, xt;
    const double pi = 3.141592653589793;
    s=1.0;
    r=1.0;
    if (x < 24.5) {
        for (k = 1; k < 61; k++) {
            r = -r*x*x*(2.0*k - 1.0)/pow(2.0*k + 1.0, 3);
            s += r;
            if (fabs(r) < fabs(s)*1.0e-12) { break; }
        }
        tth = pi/2.0 - 2.0/pi*x*s;
    } else {
        for (k = 1; k < 11; k++) {
            r = -r*pow(2.0*k-1.0, 3)/((2.0*k+1.0)*x*x);
            s += r;
            if (fabs(r) < fabs(s)*1.0e-12) { break; }
        }
        tth = 2.0/(pi*x)*s;
        t = 8.0 / x;
        xt = x + 0.25*pi;
        f0 = (((((0.18118e-2*t-0.91909e-2)*t+0.017033)*t-0.9394e-3)*t-0.051445)*t-0.11e-5)*t+0.7978846;
        g0 = (((((-0.23731e-2*t+0.59842e-2)*t+0.24437e-2)*t-0.0233178)*t+0.595e-4)*t+0.1620695)*t;
        tty = (f0*sin(xt) - g0*cos(xt))/(sqrt(x)*x);
        tth = tth + tty;
    }
    return tth;
}


void specfun_ittika(double x, double *tti, double *ttk) {

    // =========================================================
    // Purpose: Integrate [I0(t)-1]/t with respect to t from 0
    //          to x, and K0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
    //          TTK --- Integration of K0(t)/t from x to ∞
    // =========================================================

    int k;
    double b1, e0, r, r2, rc, rs;
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    static const double c[8] = {
        1.625, 4.1328125, 1.45380859375, 6.553353881835, 3.6066157150269,
        2.3448727161884, 1.7588273098916,1.4950639538279
    };

    if (x == 0.0) {
        *tti = 0.0;
        *ttk = 1.0e+300;
        return;
    }
    if (x < 40.0) {
        *tti = 1.0;
        r = 1.0;
        for (int k = 2; k <= 50; k++) {
            r = 0.25*r*(k - 1.0) / (k*k*k)*x*x;
            *tti += r;
            if (fabs(r / (*tti)) < 1.0e-12) { break; }
        }
        *tti *= 0.125*x*x;
    } else {
        *tti = 1.0;
        r = 1.0;
        for (int k = 1; k <= 8; k++) {
            r = r / x;
            *tti += c[k-1] * r;
        }
        rc = x * sqrt(2.0 * pi * x);
        *tti = (*tti) * exp(x) / rc;
    }
    if (x <= 12.0) {
        e0 = (0.5 * log(x / 2.0) + el) * log(x / 2.0)
                + pi * pi / 24.0 + 0.5 * el * el;
        b1 = 1.5 - (el + log(x / 2.0));
        rs = 1.0;
        r = 1.0;
        for (k = 2; k <= 50; k++) {
            r = 0.25*r*(k - 1.0) / (k*k*k)*x*x;
            rs += 1.0 / k;
            r2 = r * (rs + 1.0 / (2.0*k) - (el + log(x / 2.0)));
            b1 += r2;
            if (fabs(r2 / b1) < 1.0e-12) { break; }
        }
        *ttk = e0 - 0.125*x*x*b1;
    } else {
        *ttk = 1.0;
        r = 1.0;
        for (k = 1; k <= 8; k++) {
            r = -r / x;
            *ttk += c[k-1] * r;
        }
        rc = x*sqrt(2.0 / (pi * x));
        *ttk = (*ttk) * exp(-x) / rc;
    }
    return;
}


void specfun_ittjya(double x, double *ttj, double *tty) {

    // =========================================================
    // Purpose: Integrate [1-J0(t)]/t with respect to t from 0
    //          to x, and Y0(t)/t with respect to t from x to ∞
    // Input :  x   --- Variable in the limits  ( x ≥ 0 )
    // Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
    //          TTY --- Integration of Y0(t)/t from x to ∞
    // =========================================================

    int k, l;
    double a0, bj0, bj1, by0, by1, e0, b1, g0, g1, px, qx, r, r0, r1, r2, rs, t, vt, xk;
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;

    if (x == 0.0) {
        *ttj = 0.0;
        *tty = -1.0e+300;
    } else if (x <= 20.0) {
        *ttj = 1.0;
        r = 1.0;
        for (k = 2; k <= 100; k++) {
            r = -0.25 * r * (k - 1.0) / (k*k*k)*x*x;
            *ttj += r;
            if (fabs(r) < fabs(*ttj) * 1.0e-12) { break; }
        }
        *ttj *= 0.125*x*x;
        e0 = 0.5*(pi*pi/6.0 - el*el) - (0.5*log(x / 2.0) + el)*log(x/2.0);
        b1 = el + log(x / 2.0) - 1.5;
        rs = 1.0;
        r = -1.0;
        for (k = 2; k <= 100; k++) {
            r = -0.25 * r * (k - 1.0) / (k*k*k)*x*x;
            rs = rs + 1.0 / k;
            r2 = r * (rs + 1.0/(2.0 * k) - (el + log(x/2.0)));
            b1 = b1 + r2;
            if (fabs(r2) < fabs(b1) * 1.0e-12) { break; }
        }
        *tty = 2.0/pi*(e0 + 0.125 * x * x * b1);
    } else {
        a0 = sqrt(2.0 / (pi * x));
        bj0 = 0.0;
        by0 = 0.0;
        bj1 = 0.0;
        for (l = 0; l <= 1; l++) {
            vt = 4.0 * l * l;
            px = 1.0;
            r = 1.0;
            for (k = 1; k <= 14; k++) {
                r = -0.0078125*r*(vt-pow((4.0*k-3.0), 2))/(x*k)*(vt-pow((4.0*k-1.0), 2))/((2.0*k-1.0)*x);
                px += r;
                if (fabs(r) < fabs(px) * 1.0e-12) { break; }
            }
            qx = 1.0;
            r = 1.0;
            for (k = 1; k <= 14; k++) {
                r = -0.0078125*r*(vt-pow((4.0*k-1.0), 2))/(x*k)*(vt-pow((4.0*k+1.0), 2))/((2.0*k+1.0)*x);
                qx += r;
                if (fabs(r) < fabs(qx) * 1.0e-12) { break; }
            }
            qx = 0.125 * (vt - 1.0) / x * qx;
            xk = x - (0.25 + 0.5 * l)*pi;
            bj1 = a0*(px*cos(xk) - qx*sin(xk));
            by1 = a0*(px*sin(xk) + qx*cos(xk));
            if (l == 0) {
                bj0 = bj1;
                by0 = by1;
            }
        }
        t = 2.0 / x;
        g0 = 1.0;
        r0 = 1.0;
        for (k = 1; k <= 10; k++) {
            r0 = -r0*k*k*t*t;
            g0 += r0;
        }
        g1 = 1.0;
        r1 = 1.0;
        for (k = 1; k <= 10; k++) {
            r1 = -r1 * k * (k + 1.0) * t * t;
            g1 += r1;
        }
        *ttj = 2.0*g1*bj0 / (x*x) - g0*bj1 / x + el + log(x / 2.0);
        *tty = 2.0*g1*by0 / (x*x) - g0*by1 / x;
    }
    return;
}


void specfun_jdzo(int nt, double *zo, int *n, int *m, int *p) {

    // ===========================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x) and
    //          Jn'(x), and arrange them in the order of their
    //          magnitudes
    // Input :  NT    --- Number of total zeros ( NT ≤ 1200 )
    // Output:  ZO(L) --- Value of the L-th zero of Jn(x)
    //                    and Jn'(x)
    //          N(L)  --- n, order of Jn(x) or Jn'(x) associated
    //                    with the L-th zero
    //          M(L)  --- m, serial number of the zeros of Jn(x)
    //                    or Jn'(x) associated with the L-th zero
    //                    ( L is the serial number of all the
    //                      zeros of Jn(x) and Jn'(x) )
    //          P(L)  --- 0 (TM) or 1 (TE), a code for designating the
    //                    zeros of Jn(x)  or Jn'(x).
    //                    In the waveguide applications, the zeros
    //                    of Jn(x) correspond to TM modes and
    //                    those of Jn'(x) correspond to TE modes
    // Routine called:    BJNDD for computing Jn(x), Jn'(x) and
    //                    Jn''(x)
    // =============================================================

    int i, j, k, L, L0, L1, L2, mm, nm;
    double x, x0, x1, x2, xm;

    int* p1 = calloc(70, sizeof(int));
    // Compared to specfun.f we use a single array instead of separate
    // three arrays and use pointer arithmetic to access. Their usage
    // is pretty much one-shot hence does not complicate the code.

    // Note: ZO and ZOC arrays are 0-indexed in specfun.f

    // m1, n1, zoc -> 70 + 70 + 71
    double* mnzoc = calloc(211, sizeof(double));

    // bj, dj, fj -> 101 + 101 + 101
    double* bdfj = calloc(303, sizeof(double));

    x = 0;

    if (nt < 600) {
        xm = -1.0 + 2.248485*sqrt(nt) - 0.0159382*nt + 3.208775e-4*pow(nt, 1.5);
        nm = (int)(14.5 + 0.05875*nt);
        mm = (int)(0.02*nt) + 6;
    } else {
        xm = 5.0 + 1.445389*sqrt(nt) + 0.01889876*nt - 2.147763e-4*pow(nt, 1.5);
        nm = (int)(27.8 + 0.0327*nt);
        mm = (int)(0.01088*nt) + 10;
    }

    L0 = 0;
    /* 45 LOOP */
    for (i = 1; i < (nm+1); i++) {
        x1 = 0.407658 + 0.4795504*sqrt(i-1.0) + 0.983618*(i-1);
        x2 = 1.99535  + 0.8333883*sqrt(i-1.0) + 0.984584*(i-1);
        L1 = 0;

        /* 30 LOOP */
        for (j = 1; j < (mm+1); j++) {
            if ((i != 1) || (j != 1)) {
                x = x1;
                do
                {
                    specfun_bjndd(x, i, &bdfj[0], &bdfj[101], &bdfj[202]);
                    x0 = x;
                    x -= bdfj[100+i]/bdfj[201+i];
                    if (x1 > xm) { goto L20; }
                } while (fabs(x-x0) > 1e-10);
            }
            /* 15 */
            L1 += 1;
            mnzoc[69 + L1] = i-1;                    /* N[L1] */
            mnzoc[L1-1] = j;                         /* M[L1] */
            if (i == 1) { mnzoc[L1 - 1] = j-1; }
            p1[L1-1] = 1;
            mnzoc[140+L1] = x;                       /* ZOC[L1] */
            if (i <= 15) {
                x1 = x + 3.057 + 0.0122*(i-1) + (1.555 + 0.41575*(i-1))/pow(j+1, 2.0);
            } else {
                x1 = x + 2.918 + 0.01924*(i-1) + (6.26 + 0.13205*(i-1))/pow(j+1, 2.0);
            }
L20:
            x = x2;
            do {
                specfun_bjndd(x, i, &bdfj[0], &bdfj[101], &bdfj[202]);
                x0 = x;
                x -= bdfj[i-1]/bdfj[100+i];
                if (x > xm) { goto L30; }  /* Need to "continue;" twice hence goto is simpler */
            } while (fabs(x-x0) > 1e-10);
            L1 += 1;
            mnzoc[69 + L1] = i-1;
            mnzoc[L1-1] = j;
            p1[L1-1] = 0;
            mnzoc[140+L1] = x;
            if (i <= 15) {
                x2 = x + 3.11 + 0.0138*(i-1) + (0.04832 + 0.2804*(i-1))/pow(j+1, 2);
            } else {
                x2 = x + 3.001 + 0.0105*(i-1) + (11.52 + 0.48525*(i-1))/pow(j+3, 2);
            }
L30:
        ; /* Do nothing line to silence compiler */
        }
        L = L0 + L1;
        L2 = L;
        do {
            if (L0 == 0) {
                for (k = 1; k < (L+1); k++) {
                    p[k-1] = p1[k-1];
                    m[k-1] = mnzoc[k-1];   /* m[k-1] = mnzoc[k-1] */
                    n[k-1] = mnzoc[69+k];  /* n[k-1] = mnzoc[70 + (k-1)] */
                    zo[k] = mnzoc[140+k];
                }
                L1 = 0;
            } else if (L0 != 0) {
                if (zo[L0] >= mnzoc[140+L1]) {
                    p[L0+L1-1] = p[L0-1];
                    m[L0+L1-1] = m[L0-1];
                    n[L0+L1-1] = n[L0-1];
                    zo[L0+L1] = zo[L0];
                    L0 -= 1;
                } else {
                    p[L0+L1-1] = p1[L1-1];
                    m[L0+L1-1] = mnzoc[L1-1];
                    n[L0+L1-1] = mnzoc[69+L1];
                    zo[L0+L1] = mnzoc[140+L1];
                    L1 -= 1;
                }
            }
        } while (L1 != 0);
        /* 45 */
        L0 = L2;
    }
    free(bdfj);
    free(mnzoc);
    free(p1);
    return;
}


void specfun_jynb(int n, double x, int *nm, double *bj, double *dj, double *by, double *dy) {

    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x) and
    //          their derivatives
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Order of Jn(x) and Yn(x)
    // Output:  BJ(n) --- Jn(x)
    //          DJ(n) --- Jn'(x)
    //          BY(n) --- Yn(x)
    //          DY(n) --- Yn'(x)
    //          NM --- Highest order computed
    // Routines called:
    //          JYNBH to calculate the Jn and Yn
    // =====================================================

    int k;
    specfun_jynbh(n, 0, x, nm, bj, by);
    // Compute derivatives by differentiation formulas
    if (x < 1.0e-100) {
        for (k = 0; k <= n; k++) {
            dj[k] = 0.0;
            dy[k] = 1.0e+300;
        }
        dj[1] = 0.5;
    } else {
        dj[0] = -bj[1];
        for (k = 1; k <= *nm; k++) {
            dj[k] = bj[k - 1] - k / x * bj[k];
        }

        dy[0] = -by[1];
        for (k = 1; k <= *nm; k++) {
            dy[k] = by[k - 1] - k * by[k] / x;
        }
    }
    return;
}


void specfun_jynbh(int n, int nmin, double x, int *nm, double *bj, double *by) {

    // =====================================================
    // Purpose: Compute Bessel functions Jn(x), Yn(x)
    // Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
    //          n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
    //          nmin -- Lowest order computed  ( nmin ≥ 0 )
    // Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
    //          BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 to calculate the starting
    //          point for backward recurrence
    // =====================================================

    int k, m, ky;
    double pi = 3.141592653589793;
    double r2p = 0.63661977236758;
    double bs, s0, su, sv, f2, f1, f;
    double bj0, bj1, ec, by0, by1, bjk, byk;
    double p0, q0, cu, t1, p1, q1, t2;

    double a[4] = { -0.0703125, 0.112152099609375, -0.5725014209747314, 0.6074042001273483e+01 };
    double b[4] = { 0.0732421875, -0.2271080017089844, 0.1727727502584457e+01, -0.2438052969955606e+02 };
    double a1[4] = { 0.1171875, -0.144195556640625, 0.6765925884246826, -0.6883914268109947e+01 };
    double b1[4] = { -0.1025390625, 0.2775764465332031, -0.1993531733751297e+01, 0.2724882731126854e+02 };

    *nm = n;
    if (x < 1.0e-100) {
        for (k = nmin; k <= n; k++) {
            bj[k - nmin] = 0.0;
            by[k - nmin] = -1.0e+300;
        }

        if (nmin == 0) { bj[0] = 1.0; }
        return;
    }

    if ((x <= 300.0) || (n > (int)(0.9 * x))) {
        // Backward recurrence for Jn
        if (n == 0) {
            *nm = 1;
        }
        m = specfun_msta1(x, 200);
        if (m < *nm) {
            *nm = m;
        } else {
            m = specfun_msta2(x, *nm, 15);
        }
        bs = 0.0;
        su = 0.0;
        sv = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;
        f = 0.0;

        for (k = m; k >= 0; k--) {
            f = 2.0*(k + 1.0)/x*f1 - f2;
            if ((k <= *nm) && (k >= nmin)) {
                bj[k - nmin] = f;
            }
            if (k == 2 * (int)(k / 2) && k != 0) {
                bs += 2.0 * f;
                su += pow(-1, k / 2) * f / k;
            } else if (k > 1) {
                sv += pow(-1, k / 2) * k / (k * k - 1.0) * f;
            }
            f2 = f1;
            f1 = f;
        }
        s0 = bs + f;

        for (k = nmin; k <= *nm; k++) {
            bj[k - nmin] /= s0;
        }
        // Estimates for Yn at start of recurrence
        bj0 = f1 / s0;
        bj1 = f2 / s0;
        ec = log(x / 2.0) + 0.5772156649015329;
        by0 = r2p * (ec * bj0 - 4.0*su/s0);
        by1 = r2p * ((ec - 1.0)*bj1 - bj0/x - 4.0*sv/s0);

        if (0 >= nmin) { by[0 - nmin] = by0; }
        if (1 >= nmin) { by[1 - nmin] = by1; }
        ky = 2;
    } else {
        // Hankel expansion
        t1 = x - 0.25*pi;
        p0 = 1.0;
        q0 = -0.125/x;

        for (k = 1; k <= 4; k++) {
            p0 += a[k - 1] * pow(x,-2*k);
            q0 += b[k - 1] * pow(x, -2*k - 1);
        }

        cu = sqrt(r2p / x);
        bj0 = cu * (p0*cos(t1) - q0*sin(t1));
        by0 = cu * (p0*sin(t1) + q0*cos(t1));

        if (0 >= nmin) {
            bj[0 - nmin] = bj0;
            by[0 - nmin] = by0;
        }

        t2 = x - 0.75*pi;
        p1 = 1.0;
        q1 = 0.375/x;

        for (k = 1; k <= 4; k++) {
            p1 += a1[k - 1] * pow(x, -2*k);
            q1 += b1[k - 1] * pow(x, -2*k - 1);
        }

        bj1 = cu * (p1*cos(t2) - q1*sin(t2));
        by1 = cu * (p1*sin(t2) + q1*cos(t2));

        if (1 >= nmin) {
            bj[1 - nmin] = bj1;
            by[1 - nmin] = by1;
        }

        for (k = 2; k <= *nm; k++) {
            bjk = 2.0*(k - 1.0)/x*bj1 - bj0;
            if (k >= nmin) { bj[k - nmin] = bjk; }
            bj0 = bj1;
            bj1 = bjk;
        }
        ky = 2;
    }
    // Forward recurrence for Yn
    for (k = ky; k <= *nm; k++) {
        byk = 2.0 * (k - 1.0) * by1 / x - by0;

        if (k >= nmin)
            by[k - nmin] = byk;

        by0 = by1;
        by1 = byk;
    }
}


void specfun_jyndd(int n, double x, double *bjn, double *djn, double *fjn, double *byn, double *dyn, double *fyn) {

    // ===========================================================
    // purpose: compute bessel functions jn(x) and yn(x), and
    //          their first and second derivatives
    // input:   x   ---  argument of jn(x) and yn(x) ( x > 0 )
    //          n   ---  order of jn(x) and yn(x)
    // output:  bjn ---  jn(x)
    //          djn ---  jn'(x)
    //          fjn ---  jn"(x)
    //          byn ---  yn(x)
    //          dyn ---  yn'(x)
    //          fyn ---  yn"(x)
    // routines called:
    //          jynbh to compute jn and yn
    // ===========================================================

    int nm = 0;
    double bj[2], by[2];

    specfun_jynbh(n+1, n, x, &nm, bj, by);
    // compute derivatives by differentiation formulas
    *bjn = bj[0];
    *byn = by[0];
    *djn = -bj[1] + n*bj[0]/x;
    *dyn = -by[1] + n*by[0]/x;
    *fjn = (n*n/(x*x) - 1.0)*(*bjn) - (*djn)/x;
    *fyn = (n*n/(x*x) - 1.0)*(*byn) - (*dyn)/x;
    return;
}


void specfun_jyzo(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1) {

    // ======================================================
    // Purpose: Compute the zeros of Bessel functions Jn(x),
    //          Yn(x), and their derivatives
    // Input :  n  --- Order of Bessel functions  (n >= 0)
    //          NT --- Number of zeros (roots)
    // Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
    //          RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
    //          RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
    //          RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
    // Routine called: JYNDD for computing Jn(x), Yn(x), and
    //                 their first and second derivatives
    // ======================================================

    /*
     * SciPy Note:
     * See GH-18859 for additional changes done by SciPy for
     * better initial condition selection in Newton iteration
     */

    int L;
    double b, h, x, x0, bjn, djn, fjn, byn, dyn, fyn;
    const double pi = 3.141592653589793;
    // -- Newton method for j_{N,L}
    // initial guess for j_{N,1}
    if (n == 0) {
        x = 2.4;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 1.85576*pow(n, 0.33333) + 1.03315/ pow(n, 0.33333);
    }
    // iterate
    L = 0;
L10:
    x0 = x;
    specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= bjn/djn;
    if (fabs(x - x0) > 1e-11) { goto L10; }

    L += 1;
    rj0[L - 1] = x;
    // initial guess for j_{N,L+1}
    if (L == 1) {
        if (n == 0) {
            x = 5.52;
        } else {
            // Expansion from https://dlmf.nist.gov/10.21#E32 and
            // coefficients from Olver 1951
            x= n + 3.24460 * pow(n, 0.33333) + 3.15824 / pow(n, 0.33333);
        }
    } else {
        // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
        x = rj0[L - 1] + (rj0[L - 1] - rj0[L - 2]);
    }
    if (L <= (n + 10)) {
        specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
        h = atan(fabs(djn) / sqrt(fabs(fjn * bjn)));
        b = -djn / (bjn * atan(h));
        x -= (h - pi/2) / b;
    }

    if (L < nt) { goto L10; }

    // -- Newton method for j_{N,L+1}'
    if (n == 0) {
        x = 3.8317;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 0.80861 * pow(n, 0.33333) + 0.07249 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L15:
    x0 = x;
    specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= djn/fjn;
    if (fabs(x-x0) > 1e-11) goto L15;
    L += 1;
    rj1[L - 1] = x;
    if (L < nt) {
        // https://dlmf.nist.gov/10.21#E20
        x = rj1[L - 1] + (rj0[L] - rj0[L - 1]);
        goto L15;
    }

    // -- Newton method for y_{N,L}
    // initial guess for y_{N,1}
    if (n == 0) {
           x = 0.89357697;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 0.93158 * pow(n, 0.33333) + 0.26035 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L20:
    x0 = x;
    specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= byn/dyn;
    if (fabs(x - x0) > 1.0e-11) goto L20;
    L += 1;
    ry0[L - 1] = x;
    // initial guess for y_{N,L+1}
    if (L == 1) {
        if (n == 0) {
            x = 3.957678419314858;
        } else {
            // Expansion from https://dlmf.nist.gov/10.21#E33 and
            // coefficients from Olver 1951
            x = n + 2.59626 * pow(n, 0.33333) + 2.022183 / pow(n, 0.33333);
        }
    } else {
        // growth of roots is approximately linear (https://dlmf.nist.gov/10.21#E19)
        x = ry0[L - 1] + (ry0[L - 1] - ry0[L - 2]);
    }
    if (L <= n+10) {
        specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
        h = atan(fabs(dyn) / sqrt(fabs(fyn * byn)));
        b = -dyn / (byn * tan(h));
        x -= (h - pi/2) / b;
    }

    if (L < nt) goto L20;

    // -- Newton method for y_{N,L+1}'
    if (n == 0) {
        x = 2.67257;
    } else {
        // https://dlmf.nist.gov/10.21#E40
        x = n + 1.8211 * pow(n, 0.33333) + 0.94001 / pow(n, 0.33333);
    }
    // iterate
    L=0;
L25:
    x0 = x;
    specfun_jyndd(n, x, &bjn, &djn, &fjn, &byn, &dyn, &fyn);
    x -= dyn/fyn;
    if (fabs(x-x0) > 1.0e-11) goto L25;
    L += 1;
    ry1[L - 1] = x;
    if (L < nt) {
        // https://dlmf.nist.gov/10.21#E20
        x=ry1[L - 1] + (ry0[L] - ry0[L - 1]);
        goto L25;
    }
    return;
}


void specfun_klvna(double x, double *ber, double *bei, double *ger, double *gei,
                   double *der, double *dei, double *her, double *hei) {

    // ======================================================
    // Purpose: Compute Kelvin functions ber x, bei x, ker x
    //          and kei x, and their derivatives  ( x > 0 )
    // Input :  x   --- Argument of Kelvin functions
    // Output:  BER --- ber x
    //          BEI --- bei x
    //          GER --- ker x
    //          GEI --- kei x
    //          DER --- ber'x
    //          DEI --- bei'x
    //          HER --- ker'x
    //          HEI --- kei'x
    // ================================================

    int k, km, m;
    double gs, r, x2, x4, pp1, pn1, qp1, qn1, r1, pp0, pn0, qp0, qn0, r0,\
           fac, xt, cs, ss, xd, xe1, xe2, xc1, xc2, cp0, cn0, sp0, sn0, rc, rs;
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    const double eps = 1.0e-15;

    if (x == 0.0) {
        *ber = 1.0;
        *bei = 0.0;
        *ger = 1.0e+300;
        *gei = -0.25 * pi;
        *der = 0.0;
        *dei = 0.0;
        *her = -1.0e+300;
        *hei = 0.0;
        return;
    }

    x2 = 0.25 * x * x;
    x4 = x2 * x2;

    if (fabs(x) < 10.0) {
        *ber = 1.0;
        r = 1.0;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / pow((2.0 * m - 1.0), 2) * x4;
            *ber += r;
            if (fabs(r) < fabs(*ber)*eps) { break; }
        }

        *bei = x2;
        r = x2;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / pow((2.0 * m + 1.0), 2) * x4;
            *bei += r;
            if (fabs(r) < fabs(*bei)*eps) { break; }
        }

        *ger = -(log(x / 2.0) + el) * (*ber) + 0.25*pi*(*bei);

        r = 1.0;
        gs = 0.0;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / pow((2.0 * m - 1.0), 2) * x4;
            gs = gs + 1.0 / (2.0 * m - 1.0) + 1.0 / (2.0 * m);
            *ger += r * gs;
            if (fabs(r * gs) < fabs(*ger)*eps) { break; }
        }

        *gei = x2 - (log(x / 2.0) + el) * (*bei) - 0.25 * pi * (*ber);

        r = x2;
        gs = 1.0;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / pow((2.0 * m + 1.0), 2) * x4;
            gs = gs + 1.0 / (2.0 * m) + 1.0 / (2.0 * m + 1.0);
            *gei += r * gs;
            if (fabs(r * gs) < fabs(*gei)*eps) { break; }
        }

        *der = -0.25 * x * x2;
        r = *der;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / m / (m + 1.0) / pow((2.0 * m + 1.0), 2) * x4;
            *der += r;
            if (fabs(r) < fabs(*der)*eps) { break; }
        }

        *dei = 0.5 * x;
        r = *dei;
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / (2.0 * m - 1.0) / (2.0 * m + 1.0) * x4;
            *dei += r;
            if (fabs(r) < fabs(*dei)*eps) { break; }
        }

        r = -0.25 * x * x2;
        gs = 1.5;
        *her = 1.5 * r - (*ber) / x - (log(x / 2.0) + el)*(*der) + 0.25*pi*(*dei);
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / m / (m + 1.0) / pow((2.0 * m + 1.0), 2) * x4;
            gs = gs + 1.0 / (2.0 * m + 1.0) + 1.0 / (2 * m + 2.0);
            *her += r * gs;
            if (fabs(r * gs) < fabs(*her)*eps) { break; }
        }

        r = 0.5 * x;
        gs = 1.0;
        *hei = 0.5 * x - (*bei) / x - (log(x / 2.0) + el) * (*dei) - 0.25 * pi * (*der);
        for (m = 1; m <= 60; m++) {
            r = -0.25 * r / (m * m) / (2 * m - 1.0) / (2 * m + 1.0) * x4;
            gs = gs + 1.0 / (2.0 * m) + 1.0 / (2 * m + 1.0);
            *hei += r * gs;
            if (fabs(r * gs) < fabs(*hei)*eps) { return; }
        }
    } else {
        pp0 = 1.0;
        pn0 = 1.0;
        qp0 = 0.0;
        qn0 = 0.0;
        r0 = 1.0;
        km = 18;
        if (fabs(x) >= 40.0) km = 10;
        fac = 1.0;
        for (k = 1; k <= km; k++) {
            fac = -fac;
            xt = 0.25 * k * pi - trunc(0.125 * k) * 2.0 * pi;
            cs = cos(xt);
            ss = sin(xt);
            r0 = 0.125 * r0 * pow((2.0 * k - 1.0), 2) / k / x;
            rc = r0 * cs;
            rs = r0 * ss;
            pp0 += rc;
            pn0 += fac * rc;
            qp0 += rs;
            qn0 += fac * rs;
        }

        xd = x / sqrt(2.0);
        xe1 = exp(xd);
        xe2 = exp(-xd);
        xc1 = 1.0 / sqrt(2.0 * pi * x);
        xc2 = sqrt(0.5 * pi / x);
        cp0 = cos(xd + 0.125 * pi);
        cn0 = cos(xd - 0.125 * pi);
        sp0 = sin(xd + 0.125 * pi);
        sn0 = sin(xd - 0.125 * pi);

        *ger = xc2 * xe2 * (pn0 * cp0 - qn0 * sp0);
        *gei = xc2 * xe2 * (-pn0 * sp0 - qn0 * cp0);
        *ber = xc1 * xe1 * (pp0 * cn0 + qp0 * sn0) - (*gei) / pi;
        *bei = xc1 * xe1 * (pp0 * sn0 - qp0 * cn0) + (*ger) / pi;

        pp1 = 1.0;
        pn1 = 1.0;
        qp1 = 0.0;
        qn1 = 0.0;
        r1 = 1.0;
        fac = 1.0;
        for (int k = 1; k <= km; k++) {
            fac = -fac;
            xt = 0.25*k*pi - (int)(0.125 * k)*2.0*pi;
            cs = cos(xt);
            ss = sin(xt);
            r1 = 0.125*r1*(4.0 - pow(2.0*k - 1.0, 2))/(k*x);
            rc = r1*cs;
            rs = r1*ss;
            pp1 += fac*rc;
            pn1 += rc;
            qp1 += fac*rs;
            qn1 += rs;
        }
        *her = xc2 * xe2 * (-pn1*cn0 + qn1*sn0);
        *hei = xc2 * xe2 * (pn1*sn0 + qn1*cn0);
        *der = xc1 * xe1 * (pp1*cp0 + qp1*sp0) - (*hei)/pi;
        *dei = xc1 * xe1 * (pp1*sp0 - qp1*cp0) + (*her)/pi;
    }
    return;
}


void specfun_klvnzo(int nt, int kd, double *zo) {

    // ====================================================
    // Purpose: Compute the zeros of Kelvin functions
    // Input :  NT  --- Total number of zeros
    //          KD  --- Function code
    //          KD=1 to 8 for ber x, bei x, ker x, kei x,
    //                    ber'x, bei'x, ker'x and kei'x,
    //                    respectively.
    // Output:  ZO(M) --- the M-th zero of Kelvin function
    //                    for code KD
    // Routine called:
    //          KLVNA for computing Kelvin functions and
    //          their derivatives
    // ====================================================

    double ber, bei, ger, gei, der, dei, her, hei;
    double rt0[9] = {0.0, 2.84891, 5.02622, 1.71854, 3.91467,
                          6.03871, 3.77268, 2.66584, 4.93181};
    double rt = rt0[kd];

    for (int m = 1; m <= nt; m++) {
        while (1) {
            specfun_klvna(rt, &ber, &bei, &ger, &gei, &der, &dei, &her, &hei);
            if (kd == 1) {
                rt -= ber / der;
            } else if (kd == 2) {
                rt -= bei / dei;
            } else if (kd == 3) {
                rt -= ger / her;
            } else if (kd == 4) {
                rt -= gei / hei;
            } else if (kd == 5) {
                rt -= der / (-bei - der / rt);
            } else if (kd == 6) {
                rt -= dei / (ber - dei / rt);
            } else if (kd == 7) {
                rt -= her / (-gei - her / rt);
            } else {
                rt -= hei / (ger - hei / rt);
            }

            if (fabs(rt - rt0[kd]) <= 5e-10) {
                break;
            } else {
                rt0[kd] = rt;
            }
        }
        zo[m - 1] = rt;
        rt += 4.44;
    }
}


void specfun_kmn(int m, int n, double c, double cv, int kd, double *df, double *dn, double *ck1, double *ck2) {

    // ===================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions
    //          and joining factors
    // ===================================================

    int nm, nn, ip, k, i, l, j;
    double cs, gk0, gk1, gk2, gk3, t, r, dnp, su0, sw, r1, r2, r3, sa0, r4, r5, g0, sb0;
    nm = 25 + (int)(0.5 * (n - m) + c);
    nn = nm + m;
    double *u =  malloc((nn + 4) * sizeof(double));
    double *v =  malloc((nn + 4) * sizeof(double));
    double *w =  malloc((nn + 4) * sizeof(double));
    double *tp = malloc((nn + 4) * sizeof(double));
    double *rk = malloc((nn + 4) * sizeof(double));

    const double eps = 1.0e-14;

    cs = c * c * kd;
    *ck1 = 0.0;
    *ck2 = 0.0;

    ip = ((n - m) % 2 == 0 ? 0 : 1);
    k = 0;

    for (i = 1; i <= nn + 3; i++) {
        k = (ip == 0 ? -2 * (i - 1) : -(2 * i - 3));
        gk0 = 2.0 * m + k;
        gk1 = (m + k) * (m + k + 1.0);
        gk2 = 2.0 * (m + k) - 1.0;
        gk3 = 2.0 * (m + k) + 3.0;

        u[i - 1] = gk0 * (gk0 - 1.0) * cs / (gk2 * (gk2 + 2.0));
        v[i - 1] = gk1 - cv + (2.0 * (gk1 - m * m) - 1.0) * cs / (gk2 * gk3);
        w[i - 1] = (k + 1.0) * (k + 2.0) * cs / ((gk2 + 2.0) * gk3);
    }

    for (k = 1; k <= m; k++) {
        t = v[m];
        for (l = 0; l <= m - k - 1; l++)
            t = v[m - l - 1] - w[m - l] * u[m - l - 1] / t;

        rk[k - 1] = -u[k - 1] / t;
    }

    r = 1.0;
    for (k = 1; k <= m; k++) {
        r = r * rk[k - 1];
        dn[k - 1] = df[0] * r;
    }

    tp[nn - 1] = v[nn];

    for (k = nn - 1; k >= m + 1; k--) {
        tp[k - 1] = v[k] - w[k + 1] * u[k] / tp[k];

        if (k > m + 1)
            rk[k - 1] = -u[k - 1] / tp[k - 1];
    }

    dnp = (m == 0 ? df[0] : dn[m - 1]);
    dn[m] = pow(-1, ip) * dnp * cs / ((2.0 * m - 1.0) * (2.0 * m + 1.0 - 4.0 * ip) * tp[m]);

    for (k = m + 2; k <= nn; k++)
        dn[k - 1] *= rk[k - 1];

    r1 = 1.0;
    for (j = 1; j <= (n + m + ip) / 2; j++) {
        r1 = r1*(j + 0.5 * (n + m + ip));
    }
    r = 1.0;
    for (j = 1; j <= 2 * m + ip; ++j){
        r *= j;
    }
    su0 = r * df[0];
    sw = 0.0;

    for (k = 2; k <= nm; ++k) {
        r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        su0 = su0 + r * df[k - 1];
        if (k > (n - m) / 2 && fabs((su0 - sw) / su0) < eps) { break; }
        sw = su0;
    }

    if (kd != 1) {
        r2 = 1.0;

        for (j = 1; j <= m; ++j)
            r2 = 2.0 * c * r2 * j;

        r3 = 1.0;

        for (j = 1; j <= (n - m - ip) / 2; ++j)
            r3 = r3 * j;

        sa0 = (2.0 * (m + ip) + 1.0) * r1 / (pow(2.0, n) * pow(c, ip) * r2 * r3 * df[0]);
        *ck1 = sa0 * su0;

        if (kd == -1) {
            free(u); free(v); free(w); free(tp); free(rk);
            return;
        }
    }

    r4 = 1.0;
    for (j = 1; j <= (n - m - ip) / 2; ++j) {
        r4 *= 4.0 * j;
    }
    r5 = 1.0;
    for (j = 1; j <= m; ++j)
        r5 = r5 * (j + m) / c;

    g0 = dn[m - 1];

    if (m == 0)
        g0 = df[0];

    sb0 = (ip + 1.0) * pow(c, ip + 1) / (2.0 * ip * (m - 2.0) + 1.0) / (2.0 * m - 1.0);
    *ck2 = pow(-1, ip) * sb0 * r4 * r5 * g0 / r1 * su0;

    free(u); free(v); free(w); free(tp); free(rk);
    return;
}


void specfun_lamn(int n, double x, int *nm, double *bl, double *dl) {

    // =========================================================
    // Purpose: Compute lambda functions and their derivatives
    // Input:   x --- Argument of lambda function
    //          n --- Order of lambda function
    // Output:  BL(n) --- Lambda function of order n
    //          DL(n) --- Derivative of lambda function
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the start
    //          point for backward recurrence
    // =========================================================

    int i, k, m;
    double bk, r, uk, bs, f, f0, f1, bg, r0, x2;

    *nm = n;
    if (fabs(x) < 1e-100) {
        for (k = 0; k <= n; k++) {
            bl[k] = 0.0;
            dl[k] = 0.0;
        }
        bl[0] = 1.0;
        dl[1] = 0.5;
        return;
    }
    if (x <= 12.0) {
        x2 = x * x;
        for (k = 0; k <= n; k++) {
            bk = 1.0;
            r = 1.0;
            for (i = 1; i <= 50; i++) {
                r = -0.25 * r * x2 / (i * (i + k));
                bk += r;

                if (fabs(r) < fabs(bk) * 1.0e-15) { break; }
            }
            bl[k] = bk;
            if (k >= 1) {
                dl[k - 1] = -0.5 * x / k * bk;
            }
        }
        uk = 1.0;
        r = 1.0;
        for (i = 1; i <= 50; i++) {
            r = -0.25 * r * x2 / (i * (i + n + 1.0));
            uk += r;

            if (fabs(r) < fabs(uk) * 1.0e-15) { break; }
        }
        dl[n] = -0.5 * x / (n + 1.0) * uk;
        return;
    }
    if (n == 0) {
        *nm = 1;
    }
    m = specfun_msta1(x, 200);
    if (m < *nm) {
        *nm = m;
    } else {
        m = specfun_msta2(x, *nm, 15);
    }
    bs = 0.0;
    f = 0.0;
    f0 = 0.0;
    f1 = 1e-100;
    for (k = m; k >= 0; k--) {
        f = 2.0 * (k + 1.0) * f1 / x - f0;
        if (k <= *nm) {
            bl[k] = f;
        }
        if (k % 2 == 0) {
            bs += 2.0 * f;
        }
        f0 = f1;
        f1 = f;
    }
    bg = bs - f;
    for (k = 0; k <= *nm; k++) {
        bl[k] /= bg;
    }
    r0 = 1.0;
    for (k = 1; k <= *nm; k++) {
        r0 = 2.0 * r0 * k / x;
        bl[k] *= r0;
    }
    dl[0] = -0.5 * x * bl[1];
    for (k = 1; k <= *nm; k++) {
        dl[k] = 2.0 * k / x * (bl[k - 1] - bl[k]);
    }
    return;
}


void specfun_lamv(double v, double x, double *vm, double *vl, double *dl) {

    // =========================================================
    // Purpose: Compute lambda function with arbitrary order v,
    //          and their derivative
    // Input :  x --- Argument of lambda function
    //          v --- Order of lambda function
    // Output:  VL(n) --- Lambda function of order n+v0
    //          DL(n) --- Derivative of lambda function
    //          VM --- Highest order computed
    // Routines called:
    //      (1) MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    //      (2) GAM0 for computing gamma function (|x| ≤ 1)
    // =========================================================

    int i, n, k, j, k0, m;
    double cs, ga, fac, r0, f0, f1, f2, f, xk, vv;
    double x2, v0, vk, bk, r, uk, qx, px, rp, a0, ck, sk, bjv0, bjv1;
    const double pi = 3.141592653589793;
    const double rp2 = 0.63661977236758;

    x = fabs(x);
    x2 = x * x;
    n = (int)v;
    v0 = v - n;
    *vm = v;

    if (x <= 12.0) {
        for (k = 0; k <= n; k++) {
            vk = v0 + k;
            bk = 1.0;
            r = 1.0;

            for (i = 1; i <= 50; i++) {
                r = -0.25 * r * x2 / (i * (i + vk));
                bk = bk + r;

                if (fabs(r) < fabs(bk) * 1.0e-15)
                    break;
            }
            vl[k] = bk;

            uk = 1.0;
            r = 1.0;

            for (i = 1; i <= 50; i++) {
                r = -0.25 * r * x2 / (i * (i + vk + 1.0));
                uk = uk + r;

                if (fabs(r) < fabs(uk) * 1.0e-15)
                    break;
            }
            dl[k] = -0.5 * x / (vk + 1.0) * uk;
        }
        return;
    }

    k0 = (x >= 50.0) ? 8 : ((x >= 35.0) ? 10 : 11);
    bjv0 = 0.0;
    bjv1 = 0.0;

    for (j = 0; j <= 1; j++) {
        vv = 4.0 * (j + v0) * (j + v0);
        px = 1.0;
        rp = 1.0;
        for (k = 1; k <= k0; k++) {
            rp = -0.78125e-2 * rp * (vv - pow((4.0 * k - 3.0), 2.0)) * (vv - pow((4.0 * k - 1.0), 2.0)) / (k * (2.0 * k - 1.0) * x2);
            px += rp;
        }

        qx = 1.0;
        rp = 1.0;
        for (k = 1; k <= k0; k++) {
            rp = -0.78125e-2 * rp * (vv - pow((4.0 * k - 1.0), 2.0)) * (vv - pow((4.0 * k + 1.0), 2.0)) / (k * (2.0 * k + 1.0) * x2);
            qx += rp;
        }

        qx = 0.125 * (vv - 1.0) * qx / x;
        xk = x - (0.5 * (j + v0) + 0.25) * pi;
        a0 = sqrt(rp2 / x);
        ck = cos(xk);
        sk = sin(xk);

        if (j == 0) bjv0 = a0 * (px * ck - qx * sk);
        if (j == 1) bjv1 = a0 * (px * ck - qx * sk);
    }

    if (v0 == 0.0) {
        ga = 1.0;
    } else {
        ga = specfun_gam0(v0);
        ga *= v0;
    }

    fac = pow(2.0 / x, v0) * ga;
    vl[0] = bjv0;
    dl[0] = -bjv1 + v0 / x * bjv0;
    vl[1] = bjv1;
    dl[1] = bjv0 - (1.0 + v0) / x * bjv1;
    r0 = 2.0 * (1.0 + v0) / x;

    if (n <= 1) {
        vl[0] *= fac;
        dl[0] = fac * dl[0] - v0 / x * vl[0];
        vl[1] *= fac * r0;
        dl[1] = fac * r0 * dl[1] - (1.0 + v0) / x * vl[1];
        return;
    }

    if (n >= 2 && n <= (int)(0.9 * x)) {
        f0 = bjv0;
        f1 = bjv1;

        for (k = 2; k <= n; k++) {
            f = 2.0 * (k + v0 - 1.0) / x * f1 - f0;
            f0 = f1;
            f1 = f;
            vl[k] = f;
        }
    } else if (n >= 2) {
        m = specfun_msta1(x, 200);
        if (m < n) {
            n = m;
        } else {
            m = specfun_msta2(x, n, 15);
        }

        f = 0.0;
        f2 = 0.0;
        f1 = 1.0e-100;

        for (k = m; k >= 0; k--) {
            f = 2.0 * (v0 + k + 1.0) / x * f1 - f2;
            if (k <= n) vl[k] = f;
            f2 = f1;
            f1 = f;
        }

        cs = 0.0;
        if (fabs(bjv0) > fabs(bjv1)) {
            cs = bjv0 / f;
        } else {
            cs = bjv1 / f2;
        }

        for (k = 0; k <= n; k++) {
            vl[k] *= cs;
        }
    }

    vl[0] *= fac;
    for (j = 1; j <= n; j++) {
        vl[j] *= fac * r0;
        dl[j - 1] = -0.5 * x / (j + v0) * vl[j];
        r0 = 2.0 * (j + v0 + 1) / x * r0;
    }

    dl[n] = 2.0 * (v0 + n) * (vl[n - 1] - vl[n]) / x;
    *vm = n + v0;
    return;
}


void specfun_lpmn(int m, int n, double x, double *pm, double *pd) {

    // =====================================================
    // Purpose: Compute the associated Legendre functions
    //          Pmn(x) and their derivatives Pmn'(x) for
    //          real argument
    // Input :  x  --- Argument of Pmn(x)
    //          m  --- Order of Pmn(x),  m = 0,1,2,...,n
    //          n  --- Degree of Pmn(x), n = 0,1,2,...,N
    //          mm --- Physical dimension of PM and PD
    // Output:  PM(m,n) --- Pmn(x)
    //          PD(m,n) --- Pmn'(x)
    // =====================================================

    int i, j, ls;
    double xq, xs;


    for (i = 0; i < (m + 1)*(n + 1); i++) {
            pm[i] = 0.0;
            pd[i] = 0.0;
    }

    pm[0] = 1.0;

    if (n == 0) {
        return;
    }

    if (fabs(x) == 1.0) {
        for (i = 1; i <= n; i++) {
            pm[i] = pow(x, i);
            pd[i] = 0.5 * i * (i + 1.0) * pow(x, i + 1);
        }

        for (i = 1; i <= m; i++) {
            for (j = 1; j <= n; j++) {
                if (i == 1) {
                    pd[n + 1 + j] = INFINITY;
                } else if (i == 2) {
                    pd[2*n + 2 + j] = -0.25 * (j + 2) * (j + 1) * j * (j - 1) * pow(x, j + 1);
                }
            }
        }
        return;
    }

    ls = (fabs(x) > 1.0 ? -1 : 1);
    xq = sqrt(ls * (1.0 - x * x));
    // Ensure connection to the complex-valued function for |x| > 1
    if (x < -1.0) {
        xq = -xq;
    }
    xs = ls * (1.0 - x * x);
    /* 30 */
    for (i = 1; i <= m; ++i) {
        pm[i * (n + 2)] = -ls * (2.0 * i - 1.0) * xq * pm[(i - 1)*(n + 2)];
    }
    /* 35 */
    for (i = 0; i <= (m > (n-1) ? n - 1: m); i++) {
        pm[i * (n + 2) + 1] = (2.0*i+1.0)*x*pm[i * (n + 2)];
    }
    /* 40 */
    for (i = 0; i <= m; i++) {
        for (j = i + 2; j <= n; j++) {
            pm[i * (n + 1) + j] = ((2.0 * j - 1.0) * x * pm[i * (n + 1) + j - 1]
                                  - (i + j - 1.0) * pm[i * (n + 1) + j - 2]) / (j - i);
        }
    }

    pd[0] = 0.0;
    /* 45 */
    for (j = 1; j <= n; j++) {
        pd[j] = ls * j * (pm[j-1] - x * pm[j]) / xs;
    }
    /* 50 */
    for (i = 1; i <= m; i++) {
        for (j = i; j <= n; j++) {
            pd[i * (n + 1) + j] = ls * i * x * pm[i * (n + 1) + j] / xs + (j + i) * (j - i + 1.0)
                                  / xq * pm[(i - 1) * (n + 1) + j];
        }
    }
    return;
}


void specfun_lpmns(int m, int n, double x, double* pm, double* pd) {

    // ========================================================
    // Purpose: Compute associated Legendre functions Pmn(x)
    //          and Pmn'(x) for a given order
    // Input :  x --- Argument of Pmn(x)
    //          m --- Order of Pmn(x),  m = 0,1,2,...,n
    //          n --- Degree of Pmn(x), n = 0,1,2,...,N
    // Output:  PM(n) --- Pmn(x)
    //          PD(n) --- Pmn'(x)
    // ========================================================

    int k;
    double coef, x0, pm0, pm1, pm2, pmk;
    for (k = 0; k <= n; k++) {
        pm[k] = 0.0;
        pd[k] = 0.0;
    }

    if (fabs(x) == 1.0) {
        for (k = 0; k <= n; k++) {
            if (m == 0) {
                pm[k] = 1.0;
                pd[k] = 0.5 * k * (k + 1.0);
                if (x < 0.0) {
                    pm[k] *= ((k % 2) == 0 ? 1 : -1 );
                    pd[k] *= (((k + 1) % 2) == 0 ? 1 : -1 );
                }
            } else if (m == 1) {
                pd[k] = 1e300;
            } else if (m == 2) {
                pd[k] = -0.25 * (k + 2.0) * (k + 1.0) * k * (k - 1.0);
                if (x < 0.0)
                    pd[k] *= (((k + 1) % 2) == 0 ? 1 : -1 );
            }
        }
        return;
    }

    x0 = fabs(1.0 - x * x);
    pm0 = 1.0;
    pmk = pm0;
    for (k = 1; k <= m; k++) {
        pmk = (2.0 * k - 1.0) * sqrt(x0) * pm0;
        pm0 = pmk;
    }
    pm1 = (2.0 * m + 1.0) * x * pm0;
    pm[m] = pmk;
    pm[m + 1] = pm1;
    for (k = m + 2; k <= n; k++) {
        pm2 = ((2.0 * k - 1.0) * x * pm1 - (k + m - 1.0) * pmk) / (k - m);
        pm[k] = pm2;
        pmk = pm1;
        pm1 = pm2;
    }

    pd[0] = ((1.0 - m) * pm[1] - x * pm[0]) / (x * x - 1.0);
    for (k = 1; k <= n; k++) {
        pd[k] = (k * x * pm[k] - (k + m) * pm[k - 1]) / (x * x - 1.0);
    }
    coef = ((m % 2) == 0 ? 1 : -1 );
    for (k = 1; k <= n; k++) {
        pm[k] *= coef;
        pd[k] *= coef;
    }
    return;
}


void specfun_lpn(int n, double x, double *pn, double *pd) {

    // ===============================================
    // Purpose: Compute Legendre polynomials Pn(x)
    //          and their derivatives Pn'(x)
    // Input :  x --- Argument of Pn(x)
    //          n --- Degree of Pn(x) ( n = 0,1,...)
    // Output:  PN(n) --- Pn(x)
    //          PD(n) --- Pn'(x)
    // ===============================================

    int k;
    double p0, p1, pf;
    pn[0] = 1.0;
    pn[1] = x;
    pd[0] = 0.0;
    pd[1] = 1.0;
    p0 = 1.0;
    p1 = x;
    for (k = 2; k <= n; k++) {
        pf = (2.0*k - 1.0)/k*x*p1 - (k - 1.0)/k*p0;
        pn[k] = pf;
        if (fabs(x) == 1.0) {
            pd[k] = 0.5*pow(x, k+1)*k*(k+1);
        } else {
            pd[k] = k*(p1 - x*pf)/(1.0 - x*x);
        }
        p0 = p1;
        p1 = pf;
    }
    return;
}


double specfun_lpmv(double x, int m, double v) {

    // =======================================================
    // Purpose: Compute the associated Legendre function
    //          Pmv(x) with an integer order and an arbitrary
    //          degree v, using recursion for large degrees
    // Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
    //          m   --- Order of Pmv(x)
    //          v   --- Degree of Pmv(x)
    // Output:  PMV --- Pmv(x)
    // Routine called:  LPMV0
    // =======================================================

    int mx, neg_m, nv, j;
    double vx, pmv, v0, p0, p1, g1, g2;
    if ((x == -1.0) && (v != (int)v)) {
        if (m == 0) {
            pmv = -1e300;
        } else {
            pmv = 1e300;
        }
        return pmv;
    }
    vx = v;
    mx = m;
    // DLMF 14.9.5
    if (v < 0) { vx = -vx -1.0; }
    neg_m = 0;
    if (m < 0) {
        if (((vx+m+1) > 0) || (vx != (int)vx)) {
            neg_m = 1;
            mx = -m;
        } else {
            // We don't handle cases where DLMF 14.9.3 doesn't help
            return NAN;
        }
    }
    nv = (int)vx;
    v0 = vx - nv;
    if ((nv > 2) && (nv > mx)) {
        // Up-recursion on degree, AMS 8.5.3 / DLMF 14.10.3
        p0 = specfun_lpmv0(v0+mx, mx, x);
        p1 = specfun_lpmv0(v0+mx+1, mx, x);
        pmv = p1;
        for (j = mx+2; j <= nv; j++) {
            pmv = ((2*(v0+j)-1)*x*p1 - (v0+j-1+mx)*p0) / (v0+j-mx);
            p0 = p1;
            p1 = pmv;
        }
    } else {
        pmv = specfun_lpmv0(vx, mx, x);
    }
    if ((neg_m != 0) && (fabs(pmv) < 1.e300)) {
        // DLMF 14.9.3
        g1 = specfun_gamma2(vx-mx+1);
        g2 =  specfun_gamma2(vx+mx+1);
        pmv = pmv*g1/g2 * pow(-1, mx);
    }
    return pmv;
}


double specfun_lpmv0(double v, int m, double x) {

    // =======================================================
    // Purpose: Compute the associated Legendre function
    //          Pmv(x) with an integer order and an arbitrary
    //          nonnegative degree v
    // Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
    //          m   --- Order of Pmv(x)
    //          v   --- Degree of Pmv(x)
    // Output:  PMV --- Pmv(x)
    // Routine called:  PSI_SPEC for computing Psi function
    // =======================================================

    int j, k, nv;
    double c0, v0, vs, pa, pss, pv0, pmv, r, r0, r1, r2, s, s0, s1, s2, qr, rg, xq;

    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    const double eps = 1e-14;

    nv = (int)v;
    v0 = v - nv;

    if (x == -1.0 && v != nv) {
        if (m == 0)
            return -1.0e+300;
        if (m != 0)
            return 1.0e+300;
    }

    c0 = 1.0;
    if (m != 0) {
        rg = v * (v + m);
        for (j = 1; j <= m - 1; j++) {
            rg *= (v * v - j * j);
        }
        xq = sqrt(1.0 - x*x);
        r0 = 1.0;
        for (j = 1; j <= m; j++) {
            r0 = 0.5*r0*xq/j;
        }
        c0 = r0*rg;
    }

    if (v0 == 0.0) {
        // DLMF 14.3.4, 14.7.17, 15.2.4
        pmv = 1.0;
        r = 1.0;
        for (k = 1; k <= nv - m; k++) {
            r = 0.5 * r * (-nv + m + k - 1.0) * (nv + m + k) / (k * (k + m)) * (1.0 + x);
            pmv += r;
        }
        return pow(-1, nv)*c0*pmv;
    } else {
        if (x >= -0.35) {
            // DLMF 14.3.4, 15.2.1
            pmv = 1.0;
            r = 1.0;
            for (k = 1; k <= 100; k++) {
                r = 0.5 * r * (-v + m + k - 1.0) * (v + m + k) / (k * (m + k)) * (1.0 - x);
                pmv += r;
                if (k > 12 && fabs(r / pmv) < eps) { break; }
            }
            return pow(-1, m)*c0*pmv;
        } else {
            // DLMF 14.3.5, 15.8.10
            vs = sin(v * pi) / pi;
            pv0 = 0.0;
            if (m != 0) {
                qr = sqrt((1.0 - x) / (1.0 + x));
                r2 = 1.0;
                for (j = 1; j <= m; j++) {
                    r2 *= qr * j;
                }
                s0 = 1.0;
                r1 = 1.0;
                for (k = 1; k <= m - 1; k++) {
                    r1 = 0.5 * r1 * (-v + k - 1) * (v + k) / (k * (k - m)) * (1.0 + x);
                    s0 += r1;
                }
                pv0 = -vs * r2 / m * s0;
            }

            pa = 2.0 * (specfun_psi_spec(v) + el) + pi / tan(pi * v) + 1.0 / v;
            s1 = 0.0;
            for (j = 1; j <= m; j++) {
                s1 += (j * j + v * v) / (j * (j * j - v * v));
            }
            pmv = pa + s1 - 1.0 / (m - v) + log(0.5 * (1.0 + x));
            r = 1.0;
            for (k = 1; k <= 100; k++) {
                r = 0.5 * r * (-v + m + k - 1.0) * (v + m + k) / (k * (k + m)) * (1.0 + x);
                s = 0.0;
                for (j = 1; j <= m; j++)
                    s += ((k + j) * (k + j) + v * v) / ((k + j) * ((k + j) * (k + j) - v * v));

                s2 = 0.0;
                for (j = 1; j <= k; j++)
                    s2 = s2 + 1.0 / (j * (j * j - v * v));

                pss = pa + s + 2.0 * v * v * s2 - 1.0 / (m + k - v) + log(0.5 * (1.0 + x));
                r2 = pss * r;
                pmv += r2;
                if (fabs(r2 / pmv) < eps) { break; }
            }
            return pv0 + pmv * vs * c0;
        }
    }
}


void specfun_lqmn(double x, int m, int n, double *qm, double *qd) {

    // ==========================================================
    // Purpose: Compute the associated Legendre functions of the
    //          second kind, Qmn(x) and Qmn'(x)
    // Input :  x  --- Argument of Qmn(x)
    //          m  --- Order of Qmn(x)  ( m = 0,1,2,… )
    //          n  --- Degree of Qmn(x) ( n = 0,1,2,… )
    //          mm --- Physical dimension of QM and QD
    // Output:  QM(m,n) --- Qmn(x)
    //          QD(m,n) --- Qmn'(x)
    // ==========================================================

    double q0, q1, q10, qf, qf0, qf1, qf2, xs, xq;
    int i, j, k, km, ls;

    if (fabs(x) == 1.0) {
        for (i = 0; i < (m + 1)*(n + 1); i++) {
            qm[i] = 1e300;
            qd[i] = 1e300;
        }
        return;
    }
    ls = 1;
    if (fabs(x) > 1.0) { ls = -1; }
    xs = ls*(1.0 - x*x);
    xq = sqrt(xs);
    q0 = 0.5*log(fabs((x + 1.0)/(x - 1.0)));
    if (fabs(x) < 1.0001) {
        qm[0] = q0;
        qm[1] = x*q0 - 1.0;
        qm[n + 1] = -1.0 / xq;
        qm[n + 2] = -ls*xq*(q0 + x/(1. - x*x));
        for (i = 0; i <= 1; i++) {
            for (j = 2; j <= n; j++) {
                qm[i * (n + 1) + j] = ((2.0*j - 1.)*x*qm[i * (n + 1) + j - 1]
                                      - (j+i-1)*qm[i * (n + 1) + j - 2]) / (j-i);
            }
        }
        /* 15 */
        for (i = 2; i <= m; i++) {
            for (j = 0; j <= n; j++) {
                qm[i * (n + 1) + j] = -2.0*(i-1.0)*x/xq*qm[(i - 1) * (n + 1) + j]
                                      - ls*(j+i-1.0)*(j-i+2.0)*qm[(i - 2) * (n + 1) + j];
            }
        }
    /* 20 */
    } else {
        if (fabs(x) > 1.1) {
            km = 40 + m + n;
        } else {
            km = (40 + m + n)*((int)(-1. - 1.8*log(x - 1.)));
        }
        qf2 = 0.0;
        qf1 = 1.0;
        qf0 = 0.0;
        for (k = km; k >= 0; k--) {
            qf0=((2.0*k + 3.0)*x*qf1 - (k + 2.0)*qf2)/(k + 1.0);
            if (k <= n) {
                qm[k] = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }
        /* 25 */
        for (k = 0; k <= n; k++) {
            qm[k] *= q0/qf0;
        }
        /* 30 */
        qf2 = 0.0;
        qf1 = 1.0;
        for (k = km; k >= 0; k--) {
            qf0 = ((2.0*k + 3.0)*x*qf1 - (k + 1.0)*qf2)/(k + 2.0);
            if (k <= n) {
                qm[n + 1 + k] = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }
        /* 35 */
        q10 = -1.0 / xq;
        for (k = 0; k <= n; k++) {
            qm[n + 1 + k] *= q10/qf0;
        }
        /* 40 */
        for (j = 0; j <= n; j++) {
            q0 = qm[j];
            q1 = qm[n + 1 + j];
            for (i = 0; i <= (m-2); i++) {
                qf = -2.*(i+1.)*x/xq*q1 + (j-i)*(j+i+1.)*q0;
                qm[(i + 2) * (n + 1) + j] = qf;
                q0 = q1;
                q1 = qf;
            }
        }
    }
    /* 45 */
    qd[0] = ls / xs;
    for (j = 1; j <= n; j++) {
         qd[j] = ls*j*(qm[j-1] - x*qm[j])/xs;
    }
    /* 50 */
    for (i = 1; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            qd[i * (n + 1) + j] = ls*i*x/xs*qm[i * (n + 1) + j]
                                  + (i+j)*(j-i+1.)/xq*qm[(i - 1) * (n + 1) + j];
        }
    }
    return;
}


void specfun_lqnb(int n, double x, double* qn, double* qd) {

    // ====================================================
    // Purpose: Compute Legendre functions Qn(x) & Qn'(x)
    // Input :  x  --- Argument of Qn(x)
    //          n  --- Degree of Qn(x)  ( n = 0,1,2,…)
    // Output:  QN(n) --- Qn(x)
    //          QD(n) --- Qn'(x)
    // ====================================================

    int k, j, l, nl;
    double x2, q0, q1, qf, qc1, qc2, qr, qf0, qf1, qf2;
    const double eps = 1.0e-14;

    if (fabs(x) == 1.0) {
        for (k = 0; k <= n; k++) {
            qn[k] = 1.0e300;
            qd[k] = 1.0e300;
        }
        return;
    }

    if (x <= 1.021) {
        x2 = fabs((1.0 + x) / (1.0 - x));
        q0 = 0.5 * log(x2);
        q1 = x * q0 - 1.0;
        qn[0] = q0;
        qn[1] = q1;
        qd[0] = 1.0 / (1.0 - x * x);
        qd[1] = qn[0] + x * qd[0];

        for (k = 2; k <= n; k++) {
            qf = ((2.0 * k - 1.0) * x * q1 - (k - 1.0) * q0) / k;
            qn[k] = qf;
            qd[k] = (qn[k - 1] - x * qf) * k / (1.0 - x * x);
            q0 = q1;
            q1 = qf;
        }
    } else {
        qc1 = 0.0;
        qc2 = 1.0 / x;

        for (j = 1; j <= n; j++) {
            qc2 *= j / ((2.0 * j + 1.0) * x);
            if (j == n - 1) qc1 = qc2;
        }

        for (l = 0; l <= 1; l++) {
            nl = n + l;
            qf = 1.0;
            qr = 1.0;

            for (k = 1; k <= 500; k++) {
                qr = qr * (0.5 * nl + k - 1.0) * (0.5 * (nl - 1) + k) / ((nl + k - 0.5) * k * x * x);
                qf += qr;
                if (fabs(qr / qf) < eps) break;
            }

            if (l == 0) {
                qn[n - 1] = qf * qc1;
            } else {
                qn[n] = qf * qc2;
            }
        }

        qf2 = qn[n];
        qf1 = qn[n - 1];

        for (k = n; k >= 2; k--) {
            qf0 = ((2 * k - 1.0) * x * qf1 - k * qf2) / (k - 1.0);
            qn[k - 2] = qf0;
            qf2 = qf1;
            qf1 = qf0;
        }

        qd[0] = 1.0 / (1.0 - x * x);

        for (k = 1; k <= n; k++) {
            qd[k] = k * (qn[k - 1] - x * qn[k]) / (1.0 - x * x);
        }
    }
    return;
}


void specfun_lqmns(int m, int n, double x, double *qm, double *qd) {

    // ========================================================
    // Purpose: Compute associated Legendre functions Qmn(x)
    //          and Qmn'(x) for a given order
    // Input :  x --- Argument of Qmn(x)
    //          m --- Order of Qmn(x),  m = 0,1,2,...
    //          n --- Degree of Qmn(x), n = 0,1,2,...
    // Output:  QM(n) --- Qmn(x)
    //          QD(n) --- Qmn'(x)
    // ========================================================

    int l, ls, k, km;
    double xq, q0, q00, q10, q01, q11, qf0, qf1, qm0, qm1, qg0, qg1, qh0, qh1,\
           qh2, qmk, q0l, q1l, qf2, val;

    val = 0.0;
    if (fabs(x) == 1.0) { val = 1e300; }
    for (k = 0; k <= n; k++) {
        qm[k] = val;
        qd[k] = val;
    }

    if (fabs(x) == 1.0) {
        return;
    }
    ls = (fabs(x) > 1.0 ? -1 : 1);

    xq = sqrt(ls*(1.0 - x*x));
    q0 = 0.5 * log(fabs((x + 1.0) / (x - 1.0)));
    q00 = q0;
    q10 = -1.0 / xq;
    q01 = x*q0 - 1.0;
    q11 = -ls*xq*(q0 + x / (1.0 - x*x));
    qf0 = q00;
    qf1 = q10;
    qm0 = 0.0;
    qm1 = 0.0;

    for (k = 2; k <= m; k++) {
        qm0 = -2.0 * (k-1.0) / xq * x * qf1 - ls * (k-1.0) * (2.0 - k) * qf0;
        qf0 = qf1;
        qf1 = qm0;
    }

    if (m == 0) {
        qm0 = q00;
    }
    if (m == 1) {
        qm0 = q10;
    }

    qm[0] = qm0;

    if (fabs(x) < 1.0001) {
        if ((m == 0) && (n > 0)) {
            qf0 = q00;
            qf1 = q01;
            for (k = 2; k <= n; k++) {
                qf2 = ((2.0 * k - 1.0) * x * qf1 - (k - 1.0) * qf0) / k;
                qm[k] = qf2;
                qf0 = qf1;
                qf1 = qf2;
            }
        }

        qg0 = q01;
        qg1 = q11;
        for (k = 2; k <= m; k++) {
            qm1 = -2.0 * (k - 1.0) / xq * x * qg1 - ls * k * (3.0 - k) * qg0;
            qg0 = qg1;
            qg1 = qm1;
        }

        if (m == 0) {
            qm1 = q01;
        }
        if (m == 1) {
            qm1 = q11;
        }

        qm[1] = qm1;

        if ((m == 1) && (n > 1)) {
            qh0 = q10;
            qh1 = q11;
            for (k = 2; k <= n; k++) {
                qh2 = ((2.0 * k - 1.0) * x * qh1 - k * qh0) / (k - 1.0);
                qm[k] = qh2;
                qh0 = qh1;
                qh1 = qh2;
            }
        } else if (m >= 2) {
            qg0 = q00;
            qg1 = q01;
            qh0 = q10;
            qh1 = q11;
            qmk = 0.0;
            for (l = 2; l <= n; l++) {
                q0l = ((2.0 * l - 1.0) * x * qg1 - (l - 1.0) * qg0) / l;
                q1l = ((2.0 * l - 1.0) * x * qh1 - l * qh0) / (l - 1.0);
                qf0 = q0l;
                qf1 = q1l;
                for (k = 2; k <= m; k++) {
                    qmk = -2.0 * (k - 1.0) / xq * x * qf1 - ls * (k + l - 1.0) * (l + 2.0 - k) * qf0;
                    qf0 = qf1;
                    qf1 = qmk;
                }
                qm[l] = qmk;
                qg0 = qg1;
                qg1 = q0l;
                qh0 = qh1;
                qh1 = q1l;
            }
        }
    } else {
        if (fabs(x) > 1.1) {
            km = 40 + m + n;
        }
        else {
            km = (40 + m + n) * (int)(-1.0 - 1.8 * log(x - 1.0));
        }
        qf2 = 0.0;
        qf1 = 1.0;
        for (k = km; k >= 0; k--) {
            qf0 = ((2.0 * k + 3.0) * x * qf1 - (k + 2.0 - m) * qf2) / (k + m + 1.0);
            if (k <= n) {
                qm[k] = qf0;
            }
            qf2 = qf1;
            qf1 = qf0;
        }
        for (k = 0; k <= n; k++) {
            qm[k] = qm[k] * qm0 / qf0;
        }
    }

    if (fabs(x) < 1.0) {
        for (k = 0; k <= n; k++) {
            qm[k] = pow(-1, m) * qm[k];
        }
    }

    qd[0] = ((1.0 - m) * qm[1] - x * qm[0]) / (x*x - 1.0);
    for (k = 1; k <= n; k++) {
        qd[k] = (k * x * qm[k] - (k + m) * qm[k-1]) / (x*x - 1.0);
    }
    return;
}


int specfun_msta1(double x, int mp) {

    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that the magnitude of
    //          Jn(x) at that point is about 10^(-MP)
    // Input :  x     --- Argument of Jn(x)
    //          MP    --- Value of magnitude
    // Output:  MSTA1 --- Starting point
    // ===================================================

    int it, nn, n0, n1;
    double a0, f, f0, f1;

    a0 = fabs(x);
    n0 = (int)(1.1*a0) + 1;
    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0)- mp;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - mp;
    for (it = 1; it <= 20; it++) {
        nn = n1 - (n1 - n0) / (1.0 - f0/f1);
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - mp;
        if (abs(nn-n1) < 1) { break; }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}


int specfun_msta2(double x, int n, int mp) {

    // ===================================================
    // Purpose: Determine the starting point for backward
    //          recurrence such that all Jn(x) has MP
    //          significant digits
    // Input :  x  --- Argument of Jn(x)
    //          n  --- Order of Jn(x)
    //          MP --- Significant digit
    // Output:  MSTA2 --- Starting point
    // ===================================================

    int it, n0, n1, nn;
    double a0, hmp, ejn, obj, f, f0, f1;

    a0 = fabs(x);
    hmp = 0.5*mp;
    ejn = 0.5*log10(6.28*n) - n*log10(1.36*a0/n);
    if (ejn <= hmp ) {
        obj = mp;
        n0 = (int)(1.1*a0) + 1;
    } else {
        obj = hmp + ejn;
        n0 = n;
    }
    f0 = 0.5*log10(6.28*n0) - n0*log10(1.36*a0/n0) - obj;
    n1 = n0 + 5;
    f1 = 0.5*log10(6.28*n1) - n1*log10(1.36*a0/n1) - obj;
    for (it = 1; it <= 20; it++) {
        nn = n1 - (n1 - n0) / (1.0 - f0/f1);
        f = 0.5*log10(6.28*nn) - nn*log10(1.36*a0/nn) - obj;
        if (abs(nn-n1) < 1) { break; }
        n0 = n1;
        f0 = f1;
        n1 = nn;
        f1 = f;
    }
    return nn;
}


void specfun_mtu0(int kf, int m, double q, double x, double *csf, double *csd) {

    // ===============================================================
    // Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
    //          and their derivatives ( q ≥ 0 )
    // Input :  KF  --- Function code
    //                  KF=1 for computing cem(x,q) and cem'(x,q)
    //                  KF=2 for computing sem(x,q) and sem'(x,q)
    //          m   --- Order of Mathieu functions
    //          q   --- Parameter of Mathieu functions
    //          x   --- Argument of Mathieu functions (in degrees)
    // Output:  CSF --- cem(x,q) or sem(x,q)
    //          CSD --- cem'x,q) or sem'x,q)
    // Routines called:
    //      (1) CVA2 for computing the characteristic values
    //      (2) FCOEF for computing the expansion coefficients
    // ===============================================================

    int kd = 0, km = 0, ic, k;
    double a, qm, xr;
    const double eps = 1.0e-14;
    const double rd = 1.74532925199433e-2;

    if (kf == 1 && m == 2 * (int)(m / 2)) { kd = 1; }
    if (kf == 1 && m != 2 * (int)(m / 2)) { kd = 2; }
    if (kf == 2 && m != 2 * (int)(m / 2)) { kd = 3; }
    if (kf == 2 && m == 2 * (int)(m / 2)) { kd = 4; }

    a = specfun_cva2(kd, m, q);

    if (q <= 1.0) {
        qm = 7.5 + 56.1 * sqrt(q) - 134.7 * q + 90.7 * sqrt(q) * q;
    } else {
        qm = 17.0 + 3.1 * sqrt(q) - 0.126 * q + 0.0037 * sqrt(q) * q;
    }

    km = (int)(qm + 0.5 * m);

    if (km > 251) {
        *csf = NAN;
        *csd = NAN;
        return;
    }

    double *fg = calloc(251, sizeof(double));
    specfun_fcoef(kd, m, q, a, fg);

    ic = (int)(m / 2) + 1;
    xr = x * rd;
    *csf = 0.0;
    for (k = 1; k <= km; k++) {
        if (kd == 1) {
            *csf += fg[k - 1] * cos((2*k - 2) * xr);
        } else if (kd == 2) {
            *csf += fg[k - 1] * cos((2*k - 1) * xr);
        } else if (kd == 3) {
            *csf += fg[k - 1] * sin((2*k - 1) * xr);
        } else if (kd == 4) {
            *csf += fg[k - 1] * sin(2*k*xr);
        }
        if ((k >= ic) && (fabs(fg[k]) < fabs(*csf) * eps)) {
            break;
        }
    }

    *csd = 0.0;
    for (k = 1; k <= km; k++) {
        if (kd == 1) {
            *csd -= (2*k - 2) * fg[k - 1] * sin((2*k - 2) * xr);
        } else if (kd == 2) {
            *csd -= (2*k - 1) * fg[k - 1] * sin((2*k - 1) * xr);
        } else if (kd == 3) {
            *csd += (2*k - 1) * fg[k - 1] * cos((2*k - 1) * xr);
        } else if (kd == 4) {
            *csd += 2.0 * k * fg[k - 1] * cos(2*k*xr);
        }
        if ((k >= ic) && (fabs(fg[k - 1]) < fabs(*csd) * eps)) {
            break;
        }
    }
    free(fg);
    return;
}


void specfun_mtu12(int kf, int kc, int m, double q, double x, double *f1r, double *d1r, double *f2r, double *d2r) {

    // ==============================================================
    // Purpose: Compute modified Mathieu functions of the first and
    //          second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
    //          and their derivatives
    // Input:   KF --- Function code
    //                 KF=1 for computing Mcm(x,q)
    //                 KF=2 for computing Msm(x,q)
    //          KC --- Function Code
    //                 KC=1 for computing the first kind
    //                 KC=2 for computing the second kind
    //                      or Msm(2)(x,q) and Msm(2)'(x,q)
    //                 KC=3 for computing both the first
    //                      and second kinds
    //          m  --- Order of Mathieu functions
    //          q  --- Parameter of Mathieu functions ( q ≥ 0 )
    //          x  --- Argument of Mathieu functions
    // Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
    //          D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
    //          F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
    //          D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
    // Routines called:
    //      (1) CVA2 for computing the characteristic values
    //      (2) FCOEF for computing expansion coefficients
    //      (3) JYNB for computing Jn(x), Yn(x) and their
    //          derivatives
    // ==============================================================

    double eps = 1.0e-14;
    double a, qm, c1, c2, u1, u2, w1, w2;
    int kd, km, ic, k, nm = 0;

    if ((kf == 1) && (m % 2 == 0)) { kd = 1; }
    if ((kf == 1) && (m % 2 != 0)) { kd = 2; }
    if ((kf == 2) && (m % 2 != 0)) { kd = 3; }
    if ((kf == 2) && (m % 2 == 0)) { kd = 4; }

    a = specfun_cva2(kd, m, q);

    if (q <= 1.0) {
        qm = 7.5 + 56.1 * sqrt(q) - 134.7 * q + 90.7 * sqrt(q) * q;
    } else {
        qm = 17.0 + 3.1 * sqrt(q) - 0.126 * q + 0.0037 * sqrt(q) * q;
    }

    km = (int)(qm + 0.5 * m);
    if (km >= 251) {
        *f1r = NAN;
        *d1r = NAN;
        *f2r = NAN;
        *d2r = NAN;
        return;
    }

    // allocate memory after a possible NAN return
    double *fg = calloc(251, sizeof(double));
    double *bj1 = calloc(252, sizeof(double));
    double *dj1 = calloc(252, sizeof(double));
    double *bj2 = calloc(252, sizeof(double));
    double *dj2 = calloc(252, sizeof(double));
    double *by1 = calloc(252, sizeof(double));
    double *dy1 = calloc(252, sizeof(double));
    double *by2 = calloc(252, sizeof(double));
    double *dy2 = calloc(252, sizeof(double));

    specfun_fcoef(kd, m, q, a, fg);
    ic = (int)(m / 2) + 1;
    if (kd == 4) { ic = m / 2; }

    c1 = exp(-x);
    c2 = exp(x);
    u1 = sqrt(q) * c1;
    u2 = sqrt(q) * c2;
    specfun_jynb(km+1, u1, &nm, bj1, dj1, by1, dy1);
    specfun_jynb(km+1, u2, &nm, bj2, dj2, by2, dy2);
    w1 = 0.0;
    w2 = 0.0;

    if (kc != 2) {
        *f1r = 0.0;
        for (k = 1; k <= km; k++) {
            if (kd == 1) {
                *f1r += pow(-1, ic + k) * fg[k - 1] * bj1[k - 1] * bj2[k - 1];
            } else if (kd == 2 || kd == 3) {
                *f1r += pow(-1, ic + k) * fg[k - 1] * (bj1[k - 1] * bj2[k] + pow(-1, kd) * bj1[k] * bj2[k - 1]);
            } else {
                *f1r += pow(-1, ic + k) * fg[k - 1] * (bj1[k - 1] * bj2[k + 1] - bj1[k + 1] * bj2[k - 1]);
            }

            if (k >= 5 && fabs(*f1r - w1) < fabs(*f1r) * eps) { break; }
            w1 = *f1r;
        }

        *f1r /= fg[0];

        *d1r = 0.0;
        for (k = 1; k <= km; k++) {
            if (kd == 1) {
                *d1r += pow(-1, ic + k) * fg[k - 1] * (c2 * bj1[k - 1] * dj2[k - 1] - c1 * dj1[k - 1] * bj2[k - 1]);
            } else if (kd == 2 || kd == 3) {
                *d1r += pow(-1, ic + k) * fg[k - 1] * (c2 * (bj1[k - 1] * dj2[k] + pow(-1, kd) * bj1[k] * dj2[k - 1])
                        - c1 * (dj1[k - 1] * bj2[k] + pow(-1, kd) * dj1[k] * bj2[k - 1]));
            } else {
                *d1r += pow(-1, ic + k) * fg[k - 1] * (c2 * (bj1[k - 1] * dj2[k + 1] - bj1[k + 1] * dj2[k - 1])
                        - c1 * (dj1[k - 1] * bj2[k + 1] - dj1[k + 1] * bj2[k - 1]));
            }

            if (k >= 5 && fabs(*d1r - w2) < fabs(*d1r) * eps) { break; }
            w2 = *d1r;
        }
        *d1r *= sqrt(q) / fg[0];
        if (kc == 1) {
            free(fg);
            free(bj1);free(dj1);free(bj2);free(dj2);
            free(by1);free(dy1);free(by2);free(dy2);
            return;
        }
    }

    *f2r = 0.0;
    for (k = 1; k <= km; k++) {
        if (kd == 1) {
            *f2r += pow(-1, ic + k) * fg[k - 1] * bj1[k - 1] * by2[k - 1];
        } else if (kd == 2 || kd == 3) {
            *f2r += pow(-1, ic + k) * fg[k - 1] * (bj1[k - 1] * by2[k] + pow(-1, kd) * bj1[k] * by2[k - 1]);
        } else {
            *f2r += pow(-1, ic + k) * fg[k - 1] * (bj1[k - 1] * by2[k + 1] - bj1[k + 1] * by2[k - 1]);
        }

        if (k >= 5 && fabs(*f2r - w1) < fabs(*f2r) * eps) { break; }
        w1 = *f2r;
    }
    *f2r /= fg[0];

    *d2r = 0.0;
    for (k = 1; k <= km; k++) {
        if (kd == 1) {
            *d2r += pow(-1, ic + k) * fg[k - 1] * (c2 * bj1[k - 1] * dy2[k - 1] - c1 * dj1[k - 1] * by2[k - 1]);
        } else if (kd == 2 || kd == 3) {
            *d2r += pow(-1, ic + k) * fg[k - 1] * (c2 * (bj1[k - 1] * dy2[k] + pow(-1, kd) * bj1[k] * dy2[k - 1])
                    - c1 * (dj1[k - 1] * by2[k] + pow(-1, kd) * dj1[k] * by2[k - 1]));
        } else {
            *d2r += pow(-1, ic + k) * fg[k - 1] * (c2 * (bj1[k - 1] * dy2[k + 1 ] - bj1[k + 1] * dy2[k - 1])
                    - c1 * (dj1[k - 1] * by2[k + 1] - dj1[k + 1] * by2[k - 1]));
        }

        if (k >= 5 && fabs(*d2r - w2) < fabs(*d2r) * eps) { break; }
        w2 = *d2r;
    }
    *d2r = *d2r * sqrt(q) / fg[0];

    free(fg);
    free(bj1);free(dj1);free(bj2);free(dj2);
    free(by1);free(dy1);free(by2);free(dy2);
    return;
}


void specfun_pbdv(double x, double v, double *dv, double *dp, double *pdf, double *pdd) {

    // ====================================================
    // Purpose: Compute parabolic cylinder functions Dv(x)
    //          and their derivatives
    // Input:   x --- Argument of Dv(x)
    //          v --- Order of Dv(x)
    // Output:  DV(na) --- Dn+v0(x)
    //          DP(na) --- Dn+v0'(x)
    //          ( na = |n|, v0 = v-n, |v0| < 1,
    //            n = 0,±1,±2,… )
    //          PDF --- Dv(x)
    //          PDD --- Dv'(x)
    // Routines called:
    //       (1) DVSA for computing Dv(x) for small |x|
    //       (2) DVLA for computing Dv(x) for large |x|
    // ====================================================

    int ja, k, l, m, nk, nv, na;
    double xa, vh, ep, f, f0, f1, v0, v1, v2, pd, pd0, pd1, s0;

    xa = fabs(x);
    vh = v;
    v += copysign(1.0, v);
    nv = (int)v;
    v0 = v - nv;
    na = abs(nv);
    ep = exp(-0.25*x*x);
    ja = 0;
    if (na >= 1) { ja = 1; }
    if (v >= 0.0) {
        if (v0 == 0.0) {
            pd0 = ep;
            pd1 = x*ep;
        } else {
            for (l = 0; l <= ja; l++) {
                v1 = v0 + l;
                if (xa <= 5.8) {
                    pd1 = specfun_dvsa(x, v1);
                } else {
                    pd1 = specfun_dvla(x, v1);
                }
                if (l == 0) {
                    pd0 = pd1;
                }
            }
        }
        dv[0] = pd0;
        dv[1] = pd1;
        for (k = 2; k <= na; k++) {
            *pdf = x*pd1 - (k+v0-1.0)*pd0;
            dv[k] = *pdf;
            pd0 = pd1;
            pd1 = *pdf;
        }
    } else {
        if (x <= 0.0) {
            if (xa <= 5.8)
            {
                pd0 = specfun_dvsa(x, v0);
                v1 = v0 - 1.0;
                pd1 = specfun_dvsa(x, v1);
            } else {
                pd0 = specfun_dvla(x, v0);
                v1 = v0 - 1.0;
                pd1 = specfun_dvla(x, v1);
            }
            dv[0] = pd0;
            dv[1] = pd1;
            for (k = 2; k <= na; k++) {
                pd = (-x*pd1+pd0)/(k-1.0-v0);
                dv[k] = pd;
                pd0 = pd1;
                pd1 = pd;
            }
        } else if (x <= 2.0) {
            v2 = nv + v0;
            if (nv == 0) { v2 -= 1.0; }
            nk = (int)(-v2);
            f1 = specfun_dvsa(x, v2);
            v1 = v2 + 1.0;
            f0 = specfun_dvsa(x, v1);
            dv[nk] = f1;
            dv[nk-1] = f0;
            for (k = nk-2; k >= 0; k--) {
                f = x*f0 + (k-v0+1.0)*f1;
                dv[k] = f;
                f1 = f0;
                f0 = f;
            }
        } else {
            if (xa <= 5.8)
            {
                pd0 = specfun_dvsa(x, v0);
            } else {
                pd0 = specfun_dvla(x, v0);
            }
            dv[0] = pd0;
            m = 100 + na;
            f1 = 0.0;
            f0 = 1e-30;
            f = 0.0;
            for (k = m; k >= 0; k--) {
                f = x*f0 + (k-v0+1.0)*f1;
                if (k <= na) { dv[k] = f; }
                f1 = f0;
                f0 = f;
            }
            s0 = pd0/f;
            for (k = 0; k <= na; k++) { dv[k] *= s0; }
        }
    }
    for (k = 0; k < na; k++) {
        v1 = fabs(v0) + k;
        if (v >= 0.0) {
            dp[k] = 0.5*x*dv[k] - dv[k+1];
        } else {
            dp[k] = -0.5*x*dv[k] - v1*dv[k+1];
        }
    }
    *pdf = dv[na-1];
    *pdd = dp[na-1];
    v = vh;
    return;
}


void specfun_pbvv(double x, double v, double *vv, double *vp, double *pvf, double *pvd) {

    // ===================================================
    // Purpose: Compute parabolic cylinder functions Vv(x)
    //          and their derivatives
    // Input:   x --- Argument of Vv(x)
    //          v --- Order of Vv(x)
    // Output:  VV(na) --- Vv(x)
    //          VP(na) --- Vv'(x)
    //          ( na = |n|, v = n+v0, |v0| < 1
    //            n = 0,±1,±2,… )
    //          PVF --- Vv(x)
    //          PVD --- Vv'(x)
    // Routines called:
    //       (1) VVSA for computing Vv(x) for small |x|
    //       (2) VVLA for computing Vv(x) for large |x|
    // ===================================================

    int ja, k, kv, l, m, na, nv;
    double f, f0, f1, pv0, q2p, qe, s0, v0, v1, v2, vh, xa;

    const double pi = 3.141592653589793;

    xa = fabs(x);
    vh = v;
    v += copysign(1.0, v);
    nv = (int)v;
    v0 = v - nv;
    na = abs(nv);
    qe = exp(0.25*x*x);
    q2p = sqrt(2.0/pi);
    ja = 0;
    if (na >= 1) { ja = 1; }
    f = 0.0;
    if (v <= 0.0) {
        if (v0 == 0.0) {
            if (xa <= 7.5) {
                pv0 = specfun_vvsa(x, v0);
            } else {
                pv0 = specfun_vvla(x, v0);
            }
            f0 = q2p*qe;
            f1 = x*f0;
            vv[0] = pv0;
            vv[1] = f0;
            vv[2] = f1;
        } else {
            for (l = 0; l <= ja; l++) {
                v1 = v0-l;
                if (xa <= 7.5) {
                    f1 = specfun_vvsa(x, v1);
                } else {
                    f1 = specfun_vvla(x, v1);
                }
                if (l == 0) { f0 = f1; }
            }
            vv[0] = f0;
            vv[1] = f1;
        }
        kv = 2;
        if (v0 == 0.0) { kv=3; }
        for (k = kv; k <= na; k++) {
            f= x*f1 + (k-v0-2.0)*f0;
            vv[k] = f;
            f0 = f1;
            f1 = f;
        }
    } else {
        if ((x >= 0.0) && (x <= 7.5)) {
            v2 = v;
            if (v2 < 1.0) { v2 = v2+1.0; }
            f1 = specfun_vvsa(x, v2);
            v1 = v2 - 1.0;
            kv = (int)v2;
            f0 = specfun_vvsa(x, v1);
            vv[kv] = f1;
            vv[kv - 1] = f0;
            for (k = kv-2; k >= 0; k--) {
                f = x*f0 - (k+v0+2.0)*f1;
                if (k <= na) { vv[k] = f; }
                f1 = f0;
                f0 = f;
            }
        } else if (x > 7.5) {
            pv0 = specfun_vvla(x, v0);
            m = 100 + abs(na);
            vv[1] = pv0;
            f1 = 0.0;
            f0 = 1.0e-40;
            for (k = m; k >= 0; k--) {
                f= x*f0 - (k+v0+2.0)*f1;
                if (k <= na) { vv[k] = f; }
                f1 = f0;
                f0 = f;
            }
            s0 = pv0/f;
            for (k = 0; k <= na; k++) {
                vv[k] *= s0;
            }
        } else {
            if (xa <= 7.5) {
                f0 = specfun_vvsa(x, v0);
                v1 = v0 + 1.0;
                f1 = specfun_vvsa(x, v1);
            } else {
                f0 = specfun_vvla(x, v0);
                v1 = v0 + 1.0;
                f1 = specfun_vvla(x, v1);
            }
            vv[0] = f0;
            vv[1] = f1;
            for (k = 2; k <= na; k++) {
                f = (x*f1 - f0)/(k + v0);
                vv[k] = f;
                f0 = f1;
                f1 = f;
            }
        }
    }
    for (k = 0; k < na; k++) {
        v1 = v0 + k;
        if (v >= 0.0) {
            vp[k] = 0.5*x*vv[k] - (v1+1.0)*vv[k+1];
        } else {
            vp[k] = -0.5*x*vv[k] + vv[k+1];
        }
    }
    *pvf = vv[na-1];
    *pvd = vp[na-1];
    v = vh;
    return;
}


void specfun_pbwa(double a, double x, double *w1f, double *w1d, double *w2f, double *w2d) {

    // ======================================================
    // Purpose: Compute parabolic cylinder functions W(a,±x)
    //          and their derivatives
    // Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 )
    //          x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )
    // Output : W1F --- W(a,x)
    //          W1D --- W'(a,x)
    //          W2F --- W(a,-x)
    //          W2D --- W'(a,-x)
    // Routine called:
    //         CGAMA for computing complex gamma function
    // ======================================================

    int k, L1, L2;
    double d[80], d1, d2, dl, f1, f2, g1, g2, h[100], h0, h1, hl, r, r1,\
           y1d, y2d, y1f, y2f;
    double complex ug, vg;
    const double eps = 1e-15;
    const double p0 = 0.59460355750136;

    if (a == 0.0) {
        g1 = 3.625609908222;
        g2 = 1.225416702465;
    } else {
        ug = specfun_cgama(CMPLX(0.25, 0.5*a), 1);
        g1 = cabs(ug);
        vg = specfun_cgama(CMPLX(0.75, 0.5*a), 1);
        g2 = cabs(vg);
    }
    f1 = sqrt(g1/g2);
    f2 = sqrt(2.0*g2/g1);
    h0 = 1.0;
    h1 = a;
    h[0] = a;
    for (L1 = 2; L1 <= 100; L1++) {
        hl = a*h1 - 0.25*(2*L1 - 2.0)*(2*L1 - 3.0)*h0;
        h[L1 - 1] = hl;
        h0 = h1;
        h1 = hl;
    }
    y1f = 1.0;
    r = 1.0;
    for (k = 1; k <= 100; k++) {
        r = 0.5*r*x*x/(k*(2.0*k - 1.0));
        r1 = h[k - 1]*r;
        y1f += r1;
        if ((fabs(r1) <= eps*fabs(y1f)) && (k > 30)) { break; }
    }
    y1d = a;
    r = 1.0;
    for (k = 1; k < 100; k++) {
        r = 0.5*r*x*x/(k*(2.0*k + 1.0));
        r1 = h[k]*r;
        y1d += r1;
        if ((fabs(r1) <= eps*fabs(y1d)) && (k > 30)) { break; }
    }
    y1d *= x;
    d1 = 1.0;
    d2 = a;
    d[0] = 1.0;
    d[1] = a;
    for (L2 = 3; L2 <= 80; L2++) {
        dl = a*d2 - 0.25*((2*L2-1) - 2.0)*((2*L2-1) - 3.0)*d1;
        d[L2 - 1] = dl;
        d1 = d2;
        d2 = dl;
    }
    y2f = 1.0;
    r = 1.0;
    for (k = 1; k < 80; k++) {
        r = 0.5*r*x*x/(k*(2.0*k + 1.0));
        r1 = d[k]*r;
        y2f += r1;
        if ((fabs(r1) <= eps*fabs(y2f)) && (k > 30)) { break; }
    }
    y2f *= x;
    y2d = 1.0;
    r = 1.0;
    for (k = 1; k < 80; k++) {
        r = 0.5*r*x*x/(k*(2.0*k - 1.0));
        r1 = d[k]*r;
        y2d += r1;
        if ((fabs(r1) <= eps*fabs(y2d)) && (k > 30)) { break; }
    }
    *w1f = p0*(f1*y1f - f2*y2f);
    *w2f = p0*(f1*y1f + f2*y2f);
    *w1d = p0*(f1*y1d - f2*y2d);
    *w2d = p0*(f1*y1d + f2*y2d);
    return;
}


double specfun_psi_spec(double x) {

    // ======================================
    // Purpose: Compute Psi function
    // Input :  x  --- Argument of psi(x)
    // Output:  PS --- psi(x)
    // ======================================

    int k, n;
    double ps, s = 0.0, x2, xa = fabs(x);
    const double pi = 3.141592653589793;
    const double el = 0.5772156649015329;
    static const double a[8] = {
        -0.8333333333333e-01,
         0.83333333333333333e-02,
        -0.39682539682539683e-02,
         0.41666666666666667e-02,
        -0.75757575757575758e-02,
         0.21092796092796093e-01,
        -0.83333333333333333e-01,
         0.4432598039215686
    };

    if ((x == (int)x) && (x <= 0.0)) {
        return 1e300;
    } else if (xa == (int)xa) {
        n = (int)xa;
        for (k = 1; k < n; k++) {
            s += 1.0 / k;
        }
        ps = -el + s;
    } else if ((xa + 0.5) == (int)(xa + 0.5)) {
        n = (int)(xa - 0.5);
        for (k = 1; k < (n+1); k++) {
            s += 1.0 / (2.0*k - 1.0);
        }
        ps = -el + 2.0*s - 1.386294361119891;  /* 2*log(2) */
    } else {
        if (xa < 10.0) {
            n = (10.0 - (int)xa);
            for (k = 0; k < n; k++) {
                s += 1.0 / (xa + k);
            }
            xa += n;
        }
        x2 = 1.0 / (xa*xa);
        ps = log(xa) - 0.5 / xa;
        ps += x2*(((((((a[7]*x2+a[6])*x2+a[5])*x2+a[4])*x2+a[3])*x2+a[2])*x2+a[1])*x2+a[0]);
        ps -= s;
    }
    if (x < 0.0) {
        ps -= pi*cos(pi*x)/sin(pi*x) + 1.0 / x;
    }
    return ps;
}


void specfun_qstar(int m, int n, double c, double ck1, double *ck, double *qs, double *qt) {
    int ip, i, l, k;
    double r, s, sk, qs0;
    double *ap = malloc(200*sizeof(double));
    ip = ((n - m) == 2 * ((n - m) / 2) ? 0 : 1);
    r = 1.0 / pow(ck[0], 2);
    ap[0] = r;

    for (i = 1; i <= m; i++) {
        s = 0.0;
        for (l = 1; l <= i; l++) {
            sk = 0.0;
            for (k = 0; k <= l; k++)
                sk += ck[k] * ck[l - k];
            s += sk * ap[i - l];
        }
        ap[i] = -r * s;
    }
    qs0 = ap[m - 1];

    for (l = 1; l < m; l++) {
        r = 1.0;
        for (k = 1; k <= l; ++k) {
            r = r * (2.0 * k + ip) * (2.0 * k - 1.0 + ip) / pow(2.0 * k, 2);
        }
        qs0 += ap[m - l] * r;
    }
    *qs = pow(-1, ip) * (ck1) * (ck1 * qs0) / c;
    *qt = -2.0 / (ck1) * (*qs);
    free(ap);
    return;
}


void specfun_rctj(int n, double x, int *nm, double *rj, double *dj) {

    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the first
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of jn(x)  ( n = 0,1,2,... )
    // Output:  RJ(n) --- x·jn(x)
    //          DJ(n) --- [x·jn(x)]'
    //          NM --- Highest order computed
    // Routines called:
    //          MSTA1 and MSTA2 for computing the starting
    //          point for backward recurrence
    // ========================================================

    int k, m;
    double cs, f, f0, f1, rj0, rj1;

    *nm = n;
    if (fabs(x) < 1.0e-100) {
        for (int k = 0; k <= n; k++) {
            rj[k] = 0.0;
            dj[k] = 0.0;
        }
        dj[0] = 1.0;
        return;
    }
    rj[0] = sin(x);
    rj[1] = rj[0] / x - cos(x);
    rj0 = rj[0];
    rj1 = rj[1];
    cs = 0.0;
    f = 0.0;

    if (n >= 2) {
        m = specfun_msta1(x, 200);
        if (m < n) {
            *nm = m;
        } else {
            m = specfun_msta2(x, n, 15);
        }

        f0 = 0.0;
        f1 = 1.0e-100;

        for (k = m; k >= 0; k--) {
            f = (2.0 * k + 3.0) * f1 / x - f0;
            if (k <= *nm) { rj[k] = f; }
            f0 = f1;
            f1 = f;
        }
        cs = (fabs(rj0) > fabs(rj1) ?  rj0 / f : rj1 / f0);
        for (k = 0; k <= *nm; k++) { rj[k] = cs * rj[k]; }
    }
    dj[0] = cos(x);
    for (int k = 1; k <= *nm; k++) {
        dj[k] = -k * rj[k] / x + rj[k - 1];
    }
}


void specfun_rcty(int n, double x, int *nm, double *ry, double *dy) {

    // ========================================================
    // Purpose: Compute Riccati-Bessel functions of the second
    //          kind and their derivatives
    // Input:   x --- Argument of Riccati-Bessel function
    //          n --- Order of yn(x)
    // Output:  RY(n) --- x·yn(x)
    //          DY(n) --- [x·yn(x)]'
    //          NM --- Highest order computed
    // ========================================================

    int k;
    double rf0, rf1, rf2;
    *nm = n;
    if (x < 1.0e-60) {
        for (k = 0; k <= n; k++) {
            ry[k] = -1.0e+300;
            dy[k] = 1.0e+300;
        }
        ry[0] = -1.0;
        dy[0] = 0.0;
        return;
    }

    ry[0] = -cos(x);
    ry[1] = ry[0] / x - sin(x);
    rf0 = ry[0];
    rf1 = ry[1];

    for (k = 2; k <= n; k++) {
        rf2 = (2.0 * k - 1.0) * rf1 / x - rf0;
        if (fabs(rf2) > 1.0e+300) { break; }
        ry[k] = rf2;
        rf0 = rf1;
        rf1 = rf2;
    }

    *nm = k - 1;
    dy[0] = sin(x);
    for (k = 1; k <= *nm; k++) {
        dy[k] = -k * ry[k] / x + ry[k - 1];
    }
    return;
}


double specfun_refine(int kd, int m, double q, double a) {

    // =====================================================
    // Purpose: calculate the accurate characteristic value
    //          by the secant method
    // Input :  m --- Order of Mathieu functions
    //          q --- Parameter of Mathieu functions
    //          A --- Initial characteristic value
    // Output:  A --- Refineed characteristic value
    // Routine called:  CVF for computing the value of F for
    //                  characteristic equation
    // ========================================================

    int it, mj;
    double x1, x0, x, f, f0, f1;
    const double eps = 1e-14;

    mj = 10 + m;
    x0 = a;
    f0 = specfun_cvf(kd, m, q, x0, mj);
    x1 = 1.002*a;
    f1 = specfun_cvf(kd, m, q, x1, mj);
    for (it = 1; it <= 100; it++) {
        mj += 1;
        x = x1 - (x1-x0)/(1.0 - f0/f1);
        f = specfun_cvf(kd, m, q, x, mj);
        if ((fabs(1.0 - x1/x) < eps) || (f == 0.0)) { break; }
        x0 = x1;
        f0 = f1;
        x1 = x;
        f1 = f;
    }
    return x;
}


void specfun_rmn1(int m, int n, double c, double x, int kd, double *df, double *r1f, double *r1d) {

    // =======================================================
    // Purpose: Compute prolate and oblate spheroidal radial
    //          functions of the first kind for given m, n,
    //          c and x
    // Routines called:
    //      (1) SCKB for computing expansion coefficients c2k
    //      (2) SPHJ for computing the spherical Bessel
    //          functions of the first kind
    // =======================================================

    double a0, b0, cx, r, r0, r1, r2, r3, reg, sa0, suc, sud, sum, sw, sw1;
    int ip, j, k, l, lg, nm, nm1, nm2, np;

    double *ck = calloc(200, sizeof(double));
    double *dj = calloc(252, sizeof(double));
    double *sj = calloc(252, sizeof(double));
    const double eps = 1.0e-14;

    nm1 = (int)((n - m) / 2);
    ip = (n - m == 2 * nm1 ? 0 : 1);
    nm = 25 + nm1 + (int)c;
    reg = (m + nm > 80 ? 1.0e-200 : 1.0);
    r0 = reg;

    for (j = 1; j <= 2 * m + ip; ++j) {
        r0 *= j;
    }
    r = r0;
    suc = r * df[0];
    sw = 0.0;

    for (k = 2; k <= nm; k++) {
        r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        suc += r * df[k - 1];
        if ((k > nm1) && (fabs(suc - sw) < fabs(suc) * eps)) { break; }
        sw = suc;
    }

    if (x == 0.0) {
        specfun_sckb(m, n, c, df, ck);

        sum = 0.0;
        sw1 = 0.0;
        for (j = 1; j <= nm; ++j) {
            sum += ck[j - 1];
            if (fabs(sum - sw1) < fabs(sum) * eps) { break; }
            sw1 = sum;
        }

        r1 = 1.0;
        for (j = 1; j <= (n + m + ip) / 2; j++) {
            r1 = r1 * (j + 0.5 * (n + m + ip));
        }

        r2 = 1.0;
        for (j = 1; j <= m; j++) {
            r2 *= 2.0 * c * j;
        }

        r3 = 1.0;
        for (j = 1; j <= (n - m - ip) / 2; j++) {
            r3 *= j;
        }
        sa0 = (2.0 * (m + ip) + 1.0) * r1 / (pow(2.0, n) * pow(c, ip) * r2 * r3);

        if (ip == 0) {
            *r1f = sum / (sa0 * suc) * df[0] * reg;
            *r1d = 0.0;
        } else if (ip == 1) {
            *r1f = 0.0;
            *r1d = sum / (sa0 * suc) * df[0] * reg;
        }
        free(ck);free(dj);free(sj);
        return;
    }

    cx = c * x;
    nm2 = 2 * nm + m;
    specfun_sphj(nm2, cx, &nm2, sj, dj);

    a0 = pow(1.0 - kd / (x * x), 0.5 * m) / suc;
    *r1f = 0.0;
    sw = 0.0;
    lg = 0;

    for (k = 1; k <= nm; ++k) {
        l = 2 * k + m - n - 2 + ip;
        lg = (l % 4 == 0 ? 1 : -1);

        if (k == 1) {
            r = r0;
        } else {
            r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        *r1f += lg * r * df[k - 1] * sj[np];

        if ((k > nm1) && (fabs(*r1f - sw) < fabs(*r1f) * eps)) { break; }
        sw = *r1f;
    }

    *r1f *= a0;
    b0 = kd * m / pow(x, 3.0) / (1.0 - kd / (x * x)) * (*r1f);

    sud = 0.0;
    sw = 0.0;

    for (k = 1; k <= nm; k++) {
        l = 2 * k + m - n - 2 + ip;
        lg = (l % 4 == 0 ? 1 : -1);

        if (k == 1) {
            r = r0;
        } else {
            r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        sud = sud + lg * r * df[k - 1] * dj[np];

        if ((k > nm1) && (fabs(sud - sw) < fabs(sud) * eps)) { break; }
        sw = sud;
    }
    *r1d = b0 + a0 * c * sud;
    free(ck);free(dj);free(sj);
    return;
}


void specfun_rmn2l(int m, int n, double c, double x, int Kd, double *Df, double *R2f, double *R2d, int *Id) {

    // ========================================================
    // Purpose: Compute prolate and oblate spheroidal radial
    //          functions of the second kind for given m, n,
    //          c and a large cx
    // Routine called:
    //          SPHY for computing the spherical Bessel
    //          functions of the second kind
    // ========================================================


    int ip, nm1, nm, nm2, np, j, k, l, lg, id1, id2;
    double a0, b0, cx, reg, r0, r, suc, sud, sw, eps1, eps2;
    const double eps = 1.0e-14;
    double *sy = calloc(252, sizeof(double));
    double *dy = calloc(252, sizeof(double));

    ip = 1;
    nm1 = (int)((n - m) / 2);
    if (n - m == 2 * nm1) {
        ip = 0;
    }
    nm = 25 + nm1 + (int)c;
    reg = 1.0;
    if (m + nm > 80) {
        reg = 1.0e-200;
    }
    nm2 = 2 * nm + m;
    cx = c * x;
    specfun_sphy(cx, nm2, &nm2, sy, dy);
    r0 = reg;

    for (j = 1; j <= 2 * m + ip; ++j) {
        r0 *= j;
    }

    r = r0;
    suc = r * Df[0];
    sw = 0.0;

    for (k = 2; k <= nm; k++) {
        r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        suc += r * Df[k - 1];
        if ((k > nm1) && (fabs(suc - sw) < fabs(suc) * eps)) { break; }
        sw = suc;
    }

    a0 = pow(1.0 - Kd / (x * x), 0.5 * m) / suc;
    *R2f = 0.0;
    eps1 = 0.0;
    np = 0;

    for (k = 1; k <= nm; k++) {
        l = 2 * k + m - n - 2 + ip;
        lg = (l % 4 == 0 ? 1 : -1);

        if (k == 1) {
            r = r0;
        } else {
            r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        *R2f += lg * r * (Df[k - 1] * sy[np]);
        eps1 = fabs(*R2f - sw);
        if (k > nm1 && eps1 < fabs(*R2f) * eps) { break; }
        sw = *R2f;
    }

    id1 = (int)(log10(eps1 / fabs(*R2f) + eps));
    *R2f *= a0;

    if (np >= nm2) {
        *Id = 10;
        free(sy);
        free(dy);
        return;
    }

    b0 = Kd * m / pow(x, 3) / (1.0 - Kd / (x * x)) * (*R2f);
    sud = 0.0;
    eps2 = 0.0;

    for (k = 1; k <= nm; k++) {
        l = 2 * k + m - n - 2 + ip;
        lg = (l % 4 == 0 ? 1 : -1);

        if (k == 1) {
            r = r0;
        } else {
            r = r * (m + k - 1.0) * (m + k + ip - 1.5) / (k - 1.0) / (k + ip - 1.5);
        }

        np = m + 2 * k - 2 + ip;
        sud += lg * r * (Df[k - 1] * dy[np]);
        eps2 = fabs(sud - sw);
        if ((k > nm1) && (eps2 < fabs(sud) * eps)) { break; }
        sw = sud;
    }

    *R2d = b0 + a0 * c * sud;
    id2 = (int)log10(eps2 / fabs(sud) + eps);
    *Id = (id1 > id2) ? id1 : id2;
    free(sy);
    free(dy);
    return;
}


void specfun_rmn2so(int m, int n, double c, double x, double cv, int kd, double *df, double *r2f, double *r2d) {

    // =============================================================
    // Purpose: Compute oblate radial functions of the second kind
    //          with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
    // Routines called:
    //      (1) SCKB for computing the expansion coefficients c2k
    //      (2) KMN for computing the joining factors
    //      (3) QSTAR for computing the factor defined in (15.7.3)
    //      (4) CBK for computing the the expansion coefficient
    //          defined in (15.7.6)
    //      (5) GMN for computing the function defined in (15.7.4)
    //      (6) RMN1 for computing the radial function of the first
    //          kind
    // =============================================================

    int nm, ip, j;
    double ck1, ck2, r1f, r1d, qs, qt, sum, sw, gf, gd, h0;
    const double eps = 1.0e-14;
    const double pi = 3.141592653589793;

    if (fabs(df[0]) <= 1.0e-280) {
        *r2f = 1.0e+300;
        *r2d = 1.0e+300;
        return;
    }
    double *bk = calloc(200, sizeof(double));
    double *ck = calloc(200, sizeof(double));
    double *dn = calloc(200, sizeof(double));

    nm = 25 + (int)((n - m) / 2 + c);
    ip = (n - m) % 2;
    specfun_sckb(m, n, c, df, ck);
    specfun_kmn(m, n, c, cv, kd, df, dn, &ck1, &ck2);
    specfun_qstar(m, n, c, ck1, ck, &qs, &qt);
    specfun_cbk(m, n, c, cv, qt, ck, bk);

    if (x == 0.0) {
        sum = 0.0;
        sw = 0.0;

        for (j = 0; j < nm; ++j) {
            sum += ck[j];
            if (fabs(sum - sw) < fabs(sum) * eps) { break; }
            sw = sum;
        }

        if (ip == 0) {
            r1f = sum / ck1;
            *r2f = -0.5 * pi * qs * r1f;
            *r2d = qs * r1f + bk[0];
        } else {
            r1d = sum / ck1;
            *r2f = bk[0];
            *r2d = -0.5 * pi * qs * r1d;
        }
    } else {
        specfun_gmn(m, n, c, x, bk, &gf, &gd);
        specfun_rmn1(m, n, c, x, kd, df, &r1f, &r1d);
        h0 = atan(x) - 0.5 * pi;
        *r2f = qs * r1f * h0 + gf;
        *r2d = qs * (r1d * h0 + r1f / (1.0 + x * x)) + gd;
    }
    free(bk); free(ck); free(dn);
    return;
}


void specfun_rmn2sp(int m, int n, double c, double x, double cv, int kd, double *df, double *r2f, double *r2d) {

    // ======================================================
    // Purpose: Compute prolate spheroidal radial function
    //          of the second kind with a small argument
    // Routines called:
    //      (1) LPMNS for computing the associated Legendre
    //          functions of the first kind
    //      (2) LQMNS for computing the associated Legendre
    //          functions of the second kind
    //      (3) KMN for computing expansion coefficients
    //          and joining factors
    // ======================================================

    int k, j, j1, j2, l1, ki, nm3, sum, sdm;
    double ip, nm1, nm, nm2, su0, sw, sd0, su1, sd1, sd2, ga, r1, r2, r3,\
           sf, gb, spl, gc, sd, r4, spd1, spd2, su2, ck1, ck2;

    double *pm = malloc(252*sizeof(double));
    double *pd = malloc(252*sizeof(double));
    double *qm = malloc(252*sizeof(double));
    double *qd = malloc(252*sizeof(double));
    double *dn = malloc(201*sizeof(double));
    const double eps = 1.0e-14;

    nm1 = (n - m) / 2;
    nm = 25.0 + nm1 + c;
    nm2 = 2 * nm + m;
    ip = (n - m) % 2;

    specfun_kmn(m, n, c, cv, kd, df, dn, &ck1, &ck2);
    specfun_lpmns(m, nm2, x, pm, pd);
    specfun_lqmns(m, nm2, x, qm, qd);

    su0 = 0.0;
    sw = 0.0;
    for (k = 1; k <= nm; k++) {
        j = 2 * k - 2 + m + ip;
        su0 += df[k - 1] * qm[j - 1];
        if ((k > nm1) && (fabs(su0 - sw) < fabs(su0) * eps)) { break; }
        sw = su0;
    }

    sd0 = 0.0;
    for (k = 1; k <= nm; k++) {
        j = 2 * k - 2 + m + ip;
        sd0 += df[k - 1] * qd[j - 1];
        if (k > nm1 && fabs(sd0 - sw) < fabs(sd0) * eps) { break; }
        sw = sd0;
    }

    su1 = 0.0;
    sd1 = 0.0;
    for (k = 1; k <= m; k++) {
        j = m - 2 * k + ip;
        if (j < 0) {
            j = -j - 1;
        }
        su1 += dn[k - 1] * qm[j - 1];
        sd1 += dn[k - 1] * qd[j - 1];
    }

    ga = pow((x - 1.0) / (x + 1.0), 0.5 * m);
    for (k = 1; k <= m; k++) {
        j = m - 2 * k + ip;
        if (j >= 0) { continue; }
        if (j < 0) { j = -j - 1; }
        r1 = 1.0;
        for (j1 = 0; j1 < j; j1++) {
            r1 *= (m + j1);
        }
        r2 = 1.0;
        for (j2 = 1; j2 <= (m - j - 2); j2++) {
            r2 *= j2;
        }
        r3 = 1.0;
        sf = 1.0;
        for (l1 = 1; l1 <= j; l1++) {
            r3 = 0.5 * r3 * (-j + l1 - 1.0) * (j + l1) / ((m + l1) * l1) * (1.0 - x);
            sf += r3;
        }
        if (m - j >= 2) {
            gb = (m - j - 1.0) * r2;
        }
        if (m - j <= 1) {
            gb = 1.0;
        }
        spl = r1 * ga * gb * sf;
        su1 += pow(-1, (j + m)) * dn[k-1] * spl;
        spd1 = m / (x * x - 1.0) * spl;
        gc = 0.5 * j * (j + 1.0) / (m + 1.0);
        sd = 1.0;
        r4 = 1.0;
        for (l1 = 1; l1 <= j - 1; l1++) {
            r4 = 0.5 * r4 * (-j + l1) * (j + l1 + 1.0) / ((m + l1 + 1.0) * l1) * (1.0 - x);
            sd += r4;
        }
        spd2 = r1 * ga * gb * gc * sd;
        sd1 += pow(-1, (j + m)) * dn[k - 1] * (spd1 + spd2);
    }
    su2 = 0.0;
    ki = (2 * m + 1 + ip) / 2;
    nm3 = nm + ki;
    for (k = ki; k <= nm3; k++) {
        j = 2 * k - 1 - m - ip;
        su2 += dn[k - 1] * pm[j - 1];
        if ((j > m) && (fabs(su2 - sw) < fabs(su2) * eps)) { break; }
        sw = su2;
    }
    sd2 = 0.0;
    for (k = ki; k < nm3; k++) {
        j = 2 * k - 1 - m - ip;
        sd2 += dn[k - 1] * pd[j - 1];
        if (j > m && fabs(sd2 - sw) < fabs(sd2) * eps) { break; }
        sw = sd2;
    }
    sum = su0 + su1 + su2;
    sdm = sd0 + sd1 + sd2;
    *r2f = sum / ck2;
    *r2d = sdm / ck2;
    free(pm); free(pd); free(qm); free(qd); free(dn);
    return;
}


void specfun_rswfp(int m, int n, double c, double x, double cv, int kf, double *r1f, double *r1d, double *r2f, double *r2d) {

    // ==============================================================
    // Purpose: Compute prolate spheriodal radial functions of the
    //          first and second kinds, and their derivatives
    // Input :  m  --- Mode parameter, m = 0,1,2,...
    //          n  --- Mode parameter, n = m,m+1,m+2,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument of radial function ( x > 1.0 )
    //          cv --- Characteristic value
    //          KF --- Function code
    //                 KF=1 for the first kind
    //                 KF=2 for the second kind
    //                 KF=3 for both the first and second kinds
    // Output:  R1F --- Radial function of the first kind
    //          R1D --- Derivative of the radial function of
    //                  the first kind
    //          R2F --- Radial function of the second kind
    //          R2D --- Derivative of the radial function of
    //                  the second kind
    // Routines called:
    //      (1) SDMN for computing expansion coefficients dk
    //      (2) RMN1 for computing prolate and oblate radial
    //          functions of the first kind
    //      (3) RMN2L for computing prolate and oblate radial
    //          functions of the second kind for a large argument
    //      (4) RMN2SP for computing the prolate radial function
    //          of the second kind for a small argument
    // ==============================================================

    double *df = malloc(200*sizeof(double));
    int id, kd = 1;

    specfun_sdmn(m, n, c, cv, kd, df);

    if (kf != 2) {
        specfun_rmn1(m, n, c, x, kd, df, r1f, r1d);
    }
    if (kf > 1) {
        specfun_rmn2l(m, n, c, x, kd, df, r2f, r2d, &id);
        if (id > -8) {
            specfun_rmn2sp(m, n, c, x, cv, kd, df, r2f, r2d);
        }
    }
    free(df);
    return;
}


void specfun_rswfo(int m, int n, double c, double x, double cv, int kf, double *r1f, double *r1d, double *r2f, double *r2d) {

    // ==========================================================
    // Purpose: Compute oblate radial functions of the first
    //          and second kinds, and their derivatives
    // Input :  m  --- Mode parameter,  m = 0,1,2,...
    //          n  --- Mode parameter,  n = m,m+1,m+2,...
    //          c  --- Spheroidal parameter
    //          x  --- Argument (x ≥ 0)
    //          cv --- Characteristic value
    //          KF --- Function code
    //                 KF=1 for the first kind
    //                 KF=2 for the second kind
    //                 KF=3 for both the first and second kinds
    // Output:  R1F --- Radial function of the first kind
    //          R1D --- Derivative of the radial function of
    //                  the first kind
    //          R2F --- Radial function of the second kind
    //          R2D --- Derivative of the radial function of
    //                  the second kind
    // Routines called:
    //      (1) SDMN for computing expansion coefficients dk
    //      (2) RMN1 for computing prolate or oblate radial
    //          function of the first kind
    //      (3) RMN2L for computing prolate or oblate radial
    //          function of the second kind for a large argument
    //      (4) RMN2SO for computing oblate radial functions of
    //          the second kind for a small argument
    // ==========================================================

    double *df = malloc(200*sizeof(double));
    int id, kd = -1;

    specfun_sdmn(m, n, c, cv, kd, df);

    if (kf != 2) {
        specfun_rmn1(m, n, c, x, kd, df, r1f, r1d);
    }
    if (kf > 1) {
        id = 10;
        if (x > 1e-8) {
            specfun_rmn2l(m, n, c, x, kd, df, r2f, r2d, &id);
        }
        if (id > -1) {
            specfun_rmn2so(m, n, c, x, cv, kd, df, r2f, r2d);
        }
    }
    free(df);
    return;
}


void specfun_sckb(int m, int n, double c, double *df, double *ck) {

    // ======================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          DF(k) --- Expansion coefficients dk
    // Output:  CK(k) --- Expansion coefficients ck;
    //                    CK(1), CK(2), ... correspond to
    //                    c0, c2, ...
    // ======================================================

    int i, ip, i1, i2, k, nm;
    double reg, fac, sw, r, d1, d2, d3, sum, r1;

    if (c <= 1.0e-10) {
        c = 1.0e-10;
    }
    nm = 25 + (int)(0.5 * (n - m) + c);
    ip = (n - m) % 2;
    reg = ((m + nm) > 80 ? 1.0e-200 : 1.0);
    fac = -pow(0.5, m);
    sw = 0.0;

    for (k = 0; k < nm; k++) {
        fac = -fac;
        i1 = 2 * k + ip + 1;
        r = reg;

        for (i = i1; i <= i1 + 2 * m - 1; i++) {
            r *= i;
        }
        i2 = k + m + ip;
        for (i = i2; i <= i2 + k - 1; i++) {
            r *= (i + 0.5);
        }
        sum = r * df[k];
        for (i = k + 1; i <= nm; i++) {
            d1 = 2.0 * i + ip;
            d2 = 2.0 * m + d1;
            d3 = i + m + ip - 0.5;
            r = r * d2 * (d2 - 1.0) * i * (d3 + k) / (d1 * (d1 - 1.0) * (i - k) * d3);
            sum += r * df[i];
            if (fabs(sw - sum) < fabs(sum) * 1.0e-14) { break; }
            sw = sum;
        }
        r1 = reg;
        for (i = 2; i <= m + k; i++) { r1 *= i; }
        ck[k] = fac * sum / r1;
    }
}


void specfun_sdmn(int m, int n, double c, double cv, int kd, double *df) {

    // =====================================================
    // Purpose: Compute the expansion coefficients of the
    //          prolate and oblate spheroidal functions, dk
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          cv --- Characteristic value
    //          KD --- Function code
    //                 KD=1 for prolate; KD=-1 for oblate
    // Output:  DF(k) --- Expansion coefficients dk;
    //                    DF(1), DF(2), ... correspond to
    //                    d0, d2, ... for even n-m and d1,
    //                    d3, ... for odd n-m
    // =====================================================

    int nm, ip, k, kb;
    double cs, dk0, dk1, dk2, d2k, f, fs, f1, f0, fl, f2, su1,\
           su2, sw, r1, r3, r4, s0;

    nm = 25 + (int)(0.5 * (n - m) + c);

    if (c < 1e-10) {
        for (int i = 1; i <= nm; ++i) {
            df[i-1] = 0.0;
        }
        df[(n - m) / 2] = 1.0;
        return;
    }

    double *a = calloc(nm + 2, sizeof(double));
    double *d = calloc(nm + 2, sizeof(double));
    double *g = calloc(nm + 2, sizeof(double));
    cs = c*c*kd;
    ip = (n - m) % 2;

    for (int i = 1; i <= nm + 2; ++i) {
        k = (ip == 0 ? 2 * (i - 1) : 2 * i - 1);

        dk0 = m + k;
        dk1 = m + k + 1;
        dk2 = 2 * (m + k);
        d2k = 2 * m + k;

        a[i - 1] = (d2k + 2.0) * (d2k + 1.0) / ((dk2 + 3.0) * (dk2 + 5.0)) * cs;
        d[i - 1] = dk0 * dk1 + (2.0 * dk0 * dk1 - 2.0 * m * m - 1.0) / ((dk2 - 1.0) * (dk2 + 3.0)) * cs;
        g[i - 1] = k * (k - 1.0) / ((dk2 - 3.0) * (dk2 - 1.0)) * cs;
    }

    fs = 1.0;
    f1 = 0.0;
    f0 = 1.0e-100;
    kb = 0;
    df[nm] = 0.0;
    fl = 0.0;

    for (int k = nm; k >= 1; k--) {
        f = -((d[k] - cv) * f0 + a[k] * f1) / g[k];

        if (fabs(f) > fabs(df[k])) {
            df[k-1] = f;
            f1 = f0;
            f0 = f;

            if (fabs(f) > 1.0e+100) {
                for (int k1 = k; k1 <= nm; k1++)
                    df[k1 - 1] *= 1.0e-100;
                f1 *= 1.0e-100;
                f0 *= 1.0e-100;
            }
        } else {
            kb = k;
            fl = df[k];
            f1 = 1.0e-100;
            f2 = -((d[0] - cv) / a[0]) * f1;
            df[0] = f1;

            if (kb == 1) {
                fs = f2;
            } else if (kb == 2) {
                df[1] = f2;
                fs = -((d[1] - cv) * f2 + g[1] * f1) / a[1];
            } else {
                df[1] = f2;

                for (int j = 3; j <= kb + 1; j++) {
                    f = -((d[j - 2] - cv) * f2 + g[j - 2] * f1) / a[j - 2];
                    if (j <= kb) {
                        df[j-1] = f;
                    }
                    if (fabs(f) > 1.0e+100) {
                        for (int k1 = 1; k1 <= j; k1++) {
                            df[k1 - 1] *= 1.0e-100;
                        }
                        f *= 1.0e-100;
                        f2 *= 1.0e-100;
                    }
                    f1 = f2;
                    f2 = f;
                }
                fs = f;
            }
            break;
        }
    }

    su1 = 0.0;
    r1 = 1.0;

    for (int j = m + ip + 1; j <= 2 * (m + ip); j++) {
        r1 *= j;
    }
    su1 = df[0] * r1;

    for (int k = 2; k <= kb; k++) {
        r1 = -r1 * (k + m + ip - 1.5) / (k - 1.0);
        su1 += r1 * df[k - 1];
    }

    su2 = 0.0;
    sw = 0.0;

    for (int k = kb + 1; k <= nm; k++) {
        if (k != 1) {
            r1 = -r1 * (k + m + ip - 1.5) / (k - 1.0);
        }
        su2 += r1 * df[k - 1];

        if (fabs(sw - su2) < fabs(su2) * 1.0e-14) { break; }
        sw = su2;
    }
    r3 = 1.0;

    for (int j = 1; j <= (m + n + ip) / 2; j++) {
        r3 *= (j + 0.5 * (n + m + ip));
    }
    r4 = 1.0;

    for (int j = 1; j <= (n - m - ip) / 2; j++) {
        r4 *= -4.0 * j;
    }
    s0 = r3 / (fl * (su1 / fs) + su2) / r4;

    for (int k = 1; k <= kb; ++k) {
        df[k - 1] *= fl / fs * s0;
    }
    for (int k = kb + 1; k <= nm; ++k) {
        df[k - 1] *= s0;
    }
    free(a);free(d);free(g);
    return;
}


void specfun_segv(int m, int n, double c, int kd, double *cv, double *eg) {

    // =========================================================
    // Purpose: Compute the characteristic values of spheroidal
    //          wave functions
    // Input :  m  --- Mode parameter
    //          n  --- Mode parameter
    //          c  --- Spheroidal parameter
    //          KD --- Function code
    //                 KD=1 for Prolate; KD=-1 for Oblate
    // Output:  CV --- Characteristic value for given m, n and c
    //          EG(L) --- Characteristic value for mode m and n'
    //                    ( L = n' - m + 1 )
    // =========================================================


    int i, icm, j, k, k1, l, nm, nm1;
    double cs, dk0, dk1, dk2, d2k, s, t, t1, x1, xa, xb;
    // eg[<=200] is supplied by the caller

    if (c < 1e-10) {
        for (i = 1; i <= (n-m+1); i++) {
            eg[i-1] = (i+m) * (i + m -1);
        }
        *cv = eg[n-m];
        return;
    }

    // TODO: Following array sizes should be decided dynamically
    double *a = calloc(300, sizeof(double));
    double *b = calloc(100, sizeof(double));
    double *cv0 = calloc(100, sizeof(double));
    double *d = calloc(300, sizeof(double));
    double *e = calloc(300, sizeof(double));
    double *f = calloc(300, sizeof(double));
    double *g = calloc(300, sizeof(double));
    double *h = calloc(100, sizeof(double));
    icm = (n-m+2)/2;
    nm = 10 + (int)(0.5*(n-m)+c);
    cs = c*c*kd;
    k = 0;
    for (l = 0; l <= 1; l++) {
        for (i = 1; i <= nm; i++) {
            k = (l == 0 ? 2*(i - 1) : 2*i - 1);
            dk0 = m + k;
            dk1 = m + k + 1;
            dk2 = 2*(m + k);
            d2k = 2*m + k;
            a[i-1] = (d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs;
            d[i-1] = dk0*dk1+(2.0*dk0*dk1-2.0*m*m-1.0)/((dk2-1.0)*(dk2+3.0))*cs;
            g[i-1] = k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs;
        }
        for (k = 2; k <= nm; k++) {
            e[k-1] = sqrt(a[k-2]*g[k-1]);
            f[k-1] = e[k-1]*e[k-1];
        }
        f[0] = 0.0;
        e[0] = 0.0;
        xa = d[nm-1] + fabs(e[nm-1]);
        xb = d[nm-1] - fabs(e[nm-1]);
        nm1 = nm-1;
        for (i = 1; i <= nm1; i++) {
            t = fabs(e[i-1])+fabs(e[i]);
            t1 = d[i-1] + t;
            if (xa < t1) { xa = t1; }
            t1 = d[i-1] - t;
            if (t1 < xb) { xb = t1; }
        }
        for (i = 1; i <= icm; i++) {
            b[i-1] = xa;
            h[i-1] = xb;
        }
        for (k = 1; k <= icm; k++) {
            for (k1 = k; k1 <= icm; k1++) {
                if (b[k1-1] < b[k-1]) {
                    b[k-1] = b[k1-1];
                    break;
                }
            }
            if (k != 1) {
                if(h[k-1] < h[k-2]) { h[k-1] = h[k-2]; }
            }
            while (1) {
                x1 = (b[k-1]+h[k-1])/2.0;
                cv0[k-1] = x1;
                if (fabs((b[k-1] - h[k-1])/x1) < 1e-14) { break; }
                j = 0;
                s = 1.0;
                for (i = 1; i <= nm; i++) {
                    if (s == 0.0) { s += 1e-30; }
                    t = f[i-1]/s;
                    s = d[i-1] - t - x1;
                    if (s < 0.0) { j += 1; }
                }
                if (j < k) {
                    h[k-1] = x1;
                } else {
                    b[k-1] = x1;
                    if (j >= icm) {
                        b[icm - 1] = x1;
                    } else {
                        if (h[j] < x1) { h[j] = x1; }
                        if (x1 < b[j-1]) { b[j-1] = x1; }
                    }
                }
            }
            cv0[k-1] = x1;
            if (l == 0) eg[2*k-2] = cv0[k-1];
            if (l == 1) eg[2*k-1] = cv0[k-1];
        }
    }
    *cv = eg[n-m];
    free(a);free(b);free(cv0);free(d);free(e);free(f);free(g);free(h);
    return;
}


void specfun_sphj(double x, int n, int *nm, double *sj, double *dj) {

    //  MODIFIED to ALLOW N=0 CASE (ALSO IN SPHY)
    //
    //  =======================================================
    //  Purpose: Compute spherical Bessel functions jn(x) and
    //           their derivatives
    //  Input :  x --- Argument of jn(x)
    //           n --- Order of jn(x)  ( n = 0,1,… )
    //  Output:  SJ(n) --- jn(x)
    //           DJ(n) --- jn'(x)
    //           NM --- Highest order computed
    //  Routines called:
    //           MSTA1 and MSTA2 for computing the starting
    //           point for backward recurrence
    //  =======================================================

    int k, m;
    double cs, f, f0, f1, sa, sb;

    *nm = n;
    if (fabs(x) < 1e-100) {
        for (k = 0; k <= n; k++) {
            sj[k] = 0.0;
            dj[k] = 0.0;
        }
        sj[0] = 1.0;
        if (n > 0) {
            dj[0] = 1.0 / 3.0;
        }
        return;
    }
    sj[0] = sin(x)/x;
    dj[0] = (cos(x) - sin(x)/x)/x;
    if (n < 1) {
        return;
    }
    sj[1] = (sj[0] - cos(x))/x;
    if (n >= 2) {
        sa = sj[0];
        sb = sj[1];
        m = specfun_msta1(x, 200);
        if (m < n) {
            *nm = m;
        } else {
            m = specfun_msta2(x, n, 15);
        }
        f = 0.0;
        f0 = 0.0;
        f1 = 1e-100;
        for (k = m; k >= 0; k--) {
            f = (2.0*k + 3.0)*f1/x - f0;
            if (k <= *nm) { sj[k - 1] = f; }
            f0 = f1;
            f1 = f;
        }
        cs = (fabs(sa) > fabs(sb) ? sa/f : sb/f0);
        for (k = 0; k <= *nm; k++) {
            sj[k] *= cs;
        }
    }
    for (k = 1; k <= *nm; k++) {
        dj[k] = sj[k - 1] - (k + 1.0)*sj[k - 1]/x;
    }
    return;
}


void specfun_sphy(double x, int n, int *nm, double *sy, double *dy) {

    // ======================================================
    // Purpose: Compute spherical Bessel functions yn(x) and
    //          their derivatives
    // Input :  x --- Argument of yn(x) ( x ≥ 0 )
    //          n --- Order of yn(x) ( n = 0,1,… )
    // Output:  SY(n) --- yn(x)
    //          DY(n) --- yn'(x)
    //          NM --- Highest order computed
    // ======================================================

    double f, f0, f1;

    if (x < 1.0e-60) {
        for (int k = 0; k <= n; ++k) {
            sy[k] = -1.0e300;
            dy[k] = 1.0e300;
        }
        *nm = n;
        return;
    }
    sy[0] = -cos(x) / x;
    f0 = sy[0];
    dy[0] = (sin(x) + cos(x) / x) / x;

    if (n < 1) {
        *nm = n;
        return;
    }

    sy[1] = (sy[0] - sin(x)) / x;
    f1 = sy[1];

    for (int k = 2; k <= n; k++) {
        f = ((2.0 * k - 1.0) * f1 / x) - f0;
        sy[k] = f;
        if (fabs(f) >= 1.0e300) {
            *nm = k - 1;
            return;
        }
        f0 = f1;
        f1 = f;
    }
    *nm = n - 1;
    for (int k = 1; k <= *nm; k++) {
        dy[k] = sy[k - 1] - (k + 1.0) * sy[k] / x;
    }
    return;
}


double specfun_vvla(double x, double va) {

    // ===================================================
    // Purpose: Compute parabolic cylinder function Vv(x)
    //          for large argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PV --- Vv(x)
    // Routines called:
    //       (1) DVLA for computing Dv(x) for large |x|
    //       (2) GAMMA2 for computing Г(x)
    // ===================================================

    int k;
    double pv, qe, a0, r, x1, gl, dsl, pdl;

    const double pi = 3.141592653589793;
    const double eps = 1e-12;

    qe = exp(0.25*x*x);
    a0 = pow(fabs(x), -va-1)*sqrt(2.0/pi)*qe;
    r = 1.0;
    pv = 1.0;
    for (k = 1; k <= 18; k++) {
        r = 0.5 * r * (2.0*k + va - 1.0)*(2.0*k + va)/(k*x*x);
        pv += r;
        if (fabs(r/pv) < eps) { break; }
    }
    pv *= a0;
    if (x < 0.0) {
        x1 = -x;
        pdl = specfun_dvla(x1, va);
        gl = specfun_gamma2(-va);
        dsl = sin(pi*va)*sin(pi*va);
        pv = dsl*gl/pi*pdl - cos(pi*va)*pv;
    }
    return pv;
}


double specfun_vvsa(double x, double va) {

    // ===================================================
    // Purpose: Compute parabolic cylinder function Vv(x)
    //          for small argument
    // Input:   x  --- Argument
    //          va --- Order
    // Output:  PV --- Vv(x)
    // Routine called : GAMMA2 for computing Г(x)
    // ===================================================

    double a0, fac, g1, ga0, gm, gw, r, r1, sq2, sv, sv0, v1, vb0, vm, pv;

    const double eps = 1.0e-15;
    const double pi = 3.141592653589793;
    const double ep = exp(-0.25 * x * x);
    double va0 = 1.0 + 0.5 * va;

    if (x == 0.0) {
        if (((va0 <= 0.0) && (va0 == (int)va0)) || va == 0.0) {
            pv = 0.0;
        } else {
            vb0 = -0.5 * va;
            sv0 = sin(va0 * pi);
            ga0 = specfun_gamma2(va0);
            pv = pow(2.0, vb0) * sv0 / ga0;
        }
    } else {
        sq2 = sqrt(2.0);
        a0 = pow(2.0, -0.5 * va) * ep / (2.0 * pi);
        sv = sin(-(va + 0.5) * pi);
        v1 = -0.5 * va;
        g1 = specfun_gamma2(v1);
        pv = (sv + 1.0) * g1;
        r = 1.0;
        fac = 1.0;

        for (int m = 1; m <= 250; m++) {
            vm = 0.5 * (m - va);
            gm = specfun_gamma2(vm);
            r = r * sq2 * x / m;
            fac = -fac;
            gw = fac * sv + 1.0;
            r1 = gw * r * gm;
            pv += r1;
            if ((fabs(r1 / pv) < eps) && (gw != 0.0)) { break; }
        }
        pv *= a0;
    }
    return pv;
}
