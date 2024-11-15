/*
 *
 * This file accompanied with the header file cdflib.h is a C rewrite of
 * the Fortran code with the following as its original description:
 *
 * Cumulative distribution functions, inverses and parameters for
 * Beta, Binomial, Chi-square, noncentral Chi-square, F, noncentral F, Gamma,
 * negative Binomial, Normal, Poisson, Student's t distributions.
 * It uses various TOMS algorithms and Abramowitz & Stegun, also Bus Dekker
 * zero-finding algorithm.
 *
 * The original Fortran code can be found at Netlib
 * https://www.netlib.org/random/
 *
 *
 * References
 * ----------
 *
 *  J. C. P. Bus, T. J. Dekker, Two Efficient Algorithms with Guaranteed
 *  Convergence for Finding a Zero of a Function, ACM Trans. Math. Software 1:4
 *  (1975) 330-345, DOI:10.1145/355656.355659
 *
 *  M. Abramowitz and I. A. Stegun (Eds.) (1964) Handbook of Mathematical
 *  Functions with Formulas, Graphs, and Mathematical Tables. National Bureau
 *  of Standards Applied Mathematics Series, U.S. Government Printing Office,
 *  Washington, D.C..
 *
 */

/*
 *
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
 *
 */

#include "cdflib.h"

static const double PI = 3.1415926535897932384626433832795028841971693993751;

struct DinvrState
{
    double absstp;
    double abstol;
    double big;
    double fbig;
    double fx;
    double fsmall;
    double relstp;
    double reltol;
    double small;
    int status;
    double step;
    double stpmul;
    double x;
    double xhi;
    double xlb;
    double xlo;
    double xsave;
    double xub;
    double yy;
    double zx;
    double zy;
    double zz;
    int next_state;
    int qbdd;
    int qcond;
    int qdum1;
    int qdum2;
    int qhi;
    int qleft;
    int qincr;
    int qlim;
    int qok;
    int qup;
};

struct DzrorState
{
    double a;
    double atol;
    double b;
    double c;
    double d;
    double fa;
    double fb;
    double fc;
    double fd;
    double fda;
    double fdb;
    double fx;
    double m;
    double mb;
    double p;
    double q;
    double tol;
    double rtol;
    double w;
    double xhi;
    double xlo;
    double x;
    int ext;
    int status;
    int next_state;
    int first;
    int qrzero;
    int qleft;
    int qhi;
};

static const double spmpar[3] = {
    2.220446049250313e-16,                /* np.finfo(np.float64).eps    */
    2.2250738585072014e-308,              /* np.finfo(np.float64).tiny   */
    1.7976931348623157e+308,              /* np.finfo(np.float64).max    */
};

static double algdiv(double, double);
static double alngam(double);
static double alnrel(double);
static double apser(double, double, double, double);
static double basym(double, double, double, double);
static double bcorr(double, double);
static double betaln(double, double);
static double bfrac(double, double, double, double, double, double);
static struct TupleDI bgrat(double, double, double, double, double, double);
static double bpser(double, double, double, double);
static struct TupleDDI bratio(double, double, double, double);
static double brcmp1(int, double, double, double, double);
static double brcomp(double, double, double, double);
static double bup(double, double, double, double, int, double);
static struct TupleDD cumbet(double, double, double, double);
static struct TupleDD cumbin(double, double, double, double);
static struct TupleDD cumchi(double, double);
static struct TupleDD cumchn(double, double, double);
static struct TupleDD cumf(double, double, double);
static struct TupleDDI cumfnc(double, double, double, double);
static struct TupleDD cumgam(double, double);
static struct TupleDD cumnbn(double, double, double, double);
static struct TupleDD cumnor(double);
static struct TupleDD cumpoi(double, double);
static struct TupleDD cumt(double, double);
static struct TupleDD cumtnc(double, double, double);
static double devlpl(double *, int, double);
static double dinvnr(double, double);
static void dinvr(struct DinvrState *, struct DzrorState *);
static double dt1(double, double, double);
static void dzror(struct DzrorState *);
static double cdflib_erf(double);
static double erfc1(int, double);
static double esum(int, double);
static double fpser(double, double, double, double);
static double gam1(double);
static struct TupleDI gaminv(double, double, double, double);
static double gaminv_helper_30(double, double, double, double);
static double gamln(double);
static double gamln1(double);
static double cdflib_gamma(double);
static struct TupleDD grat1(double, double, double, double);
static struct TupleDD gratio(double, double, int);
static double gsumln(double, double);
static double psi(double);
static double rcomp(double, double);
static double rexp(double);
static double rlog(double);
static double rlog1(double);
static double stvaln(double);


double algdiv(double a, double b)
{
    //         Computation of ln(gamma(b)/gamma(a+b)) when b >= 8
    //
    //                             --------
    //
    //         In this algorithm, del(x) is the function defined by
    //         Ln(gamma(x)) = (x - 0.5)*ln(x) - x + 0.5*ln(2*pi) + del(x).
    //

    double c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2;
    double carr[6] = {0.833333333333333e-01, -0.277777777760991e-02,
                      0.793650666825390e-03, -0.595202931351870e-03,
                      0.837308034031215e-03, -0.165322962780713e-02};

    if (a > b) {
        h = b / a;
        c = 1./(1. + h);
        x = h/(1. + h);
        d = a + (b - 0.5);
    } else {
        h = a / b;
        c = h/(1. + h);
        x = 1./(1. + h);
        d = b + (a - 0.5);
    }
    // Set sn = (1 - x**n)/(1 - x)
    x2 = x*x;
    s3 = 1. + (x + x2);
    s5 = 1. + (x + x2*s3);
    s7 = 1. + (x + x2*s5);
    s9 = 1. + (x + x2*s7);
    s11 = 1. + (x + x2*s9);

    // Set w = del(b) - del(a + b)
    t = pow((1. / b), 2);
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0];
    w *= c / b;
    // Combine the results
    u = d * alnrel(a / b);
    v = a * (log(b) - 1.);
    return (u > v ? (w - v) - u : (w - u) - v);
}


double alngam(double x)
{
    //                   Double precision log of the gamma function
    //
    //
    //                                  Function
    //
    //
    //         Returns the natural logarithm of gamma(x).
    //
    //
    //                                  Arguments
    //
    //
    //         X --> value at which scaled log gamma is to be returned
    //                        x is double precision
    //
    //
    //                                  Method
    //
    //
    //         If x <= 6.0, then use recursion to get x below 3
    //         then apply rational approximation number 5236 of
    //         hart et al, computer approximations, john wiley and
    //         sons, ny, 1968.
    //
    //         If x > 6.0, then use recursion to get x to at least 12 and
    //         then use formula 5423 of the same source.

    double prod, xx, result, offset;
    int i, n;
    double scoefn[9] = {0.62003838007127258804e2, 0.36036772530024836321e2,
                        0.20782472531792126786e2, 0.6338067999387272343e1,
                        0.215994312846059073e1, 0.3980671310203570498e0,
                        0.1093115956710439502e0, 0.92381945590275995e-2,
                        0.29737866448101651e-2};
    double scoefd[4] = {0.62003838007126989331e2, 0.9822521104713994894e1,
                        -0.8906016659497461257e1, 0.1000000000000000000e1};
    double coef[5] = {0.83333333333333023564e-1, -0.27777777768818808e-2,
                      0.79365006754279e-3, -0.594997310889e-3, 0.8065880899e-3};

    if (x <= 6.0) {
        prod = 1.0;
        xx = x;

        if (x > 3.0) {
            while (xx > 3.0) {
                xx -= 1.0;
                prod *= xx;
            }
        }
        if (x < 2.0) {
            while (xx < 2.0) {
                prod /= xx;
                xx += 1.0;
            }
        }
        result = devlpl(scoefn, 9, xx - 2.) / devlpl(scoefd, 4, xx - 2.);
        // Compute rational approximation to gamma(x)
        return log(result * prod);
    }
    offset = 0.5*log(2.*PI);
    // If necessary make x at least 12 and carry correction in offset
    if (x <= 12.0) {
        n = (int)(12. - x);
        if (n > 0) {
            prod = 1.0;
            for (i = 0; i < n; i++) {
                prod *= x + i;
            }
            offset -= log(prod);
            xx = x + n;
        } else {
            xx = x;
        }
    } else {
        xx = x;
    }
    // Compute power series
    result = devlpl(coef, 5, pow((1./xx), 2)) / xx;
    result += offset + (xx - 0.5)*log(xx) - xx;
    return result;
}


double alnrel(double a)
{
    //    Evaluation of the function ln(1 + a)

    double p[3] = {-0.129418923021993e+01, 0.405303492862024e+00, -0.178874546012214e-01};
    double q[3] = {-0.162752256355323e+01, 0.747811014037616e+00, -0.845104217945565e-01};
    double t, t2, w;

    if (fabs(a) > 0.375) {
        return log(1. + a);
    } else {
        t = a / (a + 2.);
        t2 = t*t;
        w = ((p[2]*t2 + p[1])*t2 + p[0])*t2 + 1.;
        w /= ((q[2]*t2 + q[1])*t2 + q[0])*t2 + 1.;
        return 2*t*w;
    }
}


double apser(double a, double b, double x, double eps)
{
    //    Apser yields the incomplete beta ratio I_(1-x))(b,a) for
    //    a <= Min(eps,eps*b), b*x <= 1, And x <= 0.5. Used when
    //    a is very small. Use only if above inequalities are satisfied.

    double aj, bx, c, j, s, t, tol;
    double g = 0.577215664901532860606512090082;

    bx = b*x;
    t = x - bx;
    if (b*eps > 0.02) {
        c = log(bx) + g + t;
    } else {
        c = log(x) + psi(b) + g + t;
    }

    tol = 5.*eps*fabs(c);
    j = 1.0;
    s = 0.0;

    while (1) {
        j += 1;
        t *= (x - bx/j);
        aj = t / j;
        s += aj;
        if (fabs(aj) <= tol) { break; }
    }
    return -a * (c + s);
}


double basym(double a, double b, double lmbda, double eps)
{
    //    Asymptotic expansion for ix(a,b) for large a and b.
    //    Lambda = (a + b)*y - b  and eps is the tolerance used.
    //    It is assumed that lambda is nonnegative and that
    //    a and b are greater than or equal to 15.

    double a0[21] = { 0.0 };
    double b0[21] = { 0.0 };
    double c[21] = { 0.0 };
    double d[21] = { 0.0 };
    double bsum, dsum, f, h, h2, hn, j0, j1, r, r0, r1, s, ssum;
    double t, t0, t1, u, w, w0, z, z0, z2, zn, znm1;
    double e0 = 2. / sqrt(PI);
    double e1 = pow(2.0, (-3./2.));
    int i, imj, j, m, mmj, n, num;

    // ****** Num is the maximum value that n can take in the do loop
    //        ending at statement 50. It is required that num be even.
    //        The arrays a0, b0, c, d have dimension num + 1.
    num = 20;

    if (a < b) {
        h = a / b;
        r0 = 1./(1.+h);
        r1 = (b - a) / b;
        w0 = 1. / sqrt(a * (1. + h));
    } else {
        h = b / a;
        r0 = 1./(1.+h);
        r1 = (b - a) / a;
        w0 = 1. / sqrt(b * (1. + h));
    }
    f = a*rlog1(-lmbda/a) + b*rlog1(lmbda/b);
    t = exp(-f);
    if (t == 0.0) { return 0.0; }
    z0 = sqrt(f);
    z = 0.5*(z0/e1);
    z2 = f + f;

    a0[0] = (2./3.)*r1;
    c[0] = -0.5*a0[0];
    d[0] = -c[0];
    j0 = (0.5/e0)*erfc1(1, z0);
    j1 = e1;
    ssum = j0 + d[0]*w0*j1;

    s = 1.0;
    h2 = h*h;
    hn = 1.0;
    w = w0;
    znm1 = z;
    zn = z2;

    for (n = 2; n <= num; n += 2) {
        hn *= h2;
        a0[n-1] = 2.*r0*(1.+ h*hn)/(n + 2.);
        s += hn;
        a0[n] = 2.*r1*s/(n + 3.);

        for (i = n; i <= n + 1; i++) {
            r = -0.5*(i + 1.);
            b0[0] = r*a0[0];

            for (m = 2; m <= i; m++) {
                bsum = 0.0;
                for (j = 1; j < m; j++) {
                    mmj = m - j;
                    bsum += (j*r - mmj)*a0[j-1]*b0[mmj-1];
                }
                b0[m-1] = r*a0[m-1] + bsum/m;
            }
            c[i-1] = b0[i-1] / (i + 1.);
            dsum = 0.0;

            for (j = 1; j < i; j++) {
                imj = i - j;
                dsum += d[imj-1]*c[j-1];
            }
            d[i-1] = - (dsum+c[i-1]);
        }
        j0 = e1*znm1 + (n-1.)*j0;
        j1 = e1*zn + n*j1;
        znm1 *= z2;
        zn *= z2;
        w *= w0;
        t0 = d[n-1]*w*j0;
        w *= w0;
        t1 = d[n]*w*j1;
        ssum += t0 + t1;
        if ((fabs(t0) + fabs(t1)) <= eps*ssum) { break; }
    }
    u = exp(-bcorr(a, b));
    return e0*t*u*ssum;
}


double bcorr(double a0, double b0)
{
    //    Evaluation of  del(a0) + del(b0) - del(a0 + b0)  where
    //    ln(gamma(a)) = (a - 0.5)*ln(a) - a + 0.5*ln(2*pi) + del(a).
    //    It is assumed that a0 >= 8 And b0 >= 8.

    double a,b,c,h,s11,s3,s5,s7,s9,t,w,x,x2;
    double carr[6] = {0.833333333333333e-01, -0.277777777760991e-02,
                      0.793650666825390e-03, -0.595202931351870e-03,
                      0.837308034031215e-03, -0.165322962780713e-02};

    a = fmin(a0, b0);
    b = fmax(a0, b0);
    h = a / b;
    c = h/(1. + h);
    x = 1./(1. + h);
    x2 = x*x;
    //  Set sn = (1 - x**n)/(1 - x)
    s3 = 1. + (x + x2);
    s5 = 1. + (x + x2*s3);
    s7 = 1. + (x + x2*s5);
    s9 = 1. + (x + x2*s7);
    s11 = 1. + (x + x2*s9);
    // Set w = del(b) - del(a + b)
    t = pow((1. / b), 2);
    w = (((((carr[5]*s11
            )*t + carr[4]*s9
           )*t + carr[3]*s7
          )*t + carr[2]*s5
         )*t + carr[1]*s3
        )*t + carr[0];
    w *= c / b;
    // Compute  del(a) + w
    t = pow((1. / a), 2);
    return ((((((carr[5])*t + carr[4]
               )*t + carr[3]
              )*t + carr[2]
             )*t + carr[1]
            )*t + carr[0]
           )/a + w;
}


double betaln(double a0, double b0)
{
    //    Evaluation of the logarithm of the beta function

    double a, b, c, h, u, v, w, z;
    double e = 0.918938533204673;
    int i, n;

    a = fmin(a0, b0);
    b = fmax(a0, b0);

    if (a >= 8.0) {
        w = bcorr(a, b);
        h = a / b;
        c = h/(1. + h);
        u = -(a - 0.5)*log(c);
        v = b*alnrel(h);
        if (u > v) {
            return (((-0.5*log(b)+e)+w)-v) - u;
        } else {
            return (((-0.5*log(b)+e)+w)-u) - v;
        }
    }
    if (a < 1) {
        if (b > 8) {
            return gamln(a) + algdiv(a,b);
        } else {
            return gamln(a) + (gamln(b) - gamln(a+b));
        }
    }

    if (a <= 2) {
        if (b <= 2) {
            return gamln(a) + gamln(b) - gsumln(a, b);
        }
        if (b >= 8) {
            return gamln(a) + algdiv(a, b);
        }
        w = 0.;
    }

    if (a > 2) {
        if (b <= 1000) {
            n = (int)(a - 1.);
            w = 1.;
            for (i = 0; i < n; i++) {
                a -= 1.0;
                h = a / b;
                w *= h/(1.+h);
            }
            w = log(w);
            if (b >= 8.0) {
                return w + gamln(a) + algdiv(a, b);
            }
        } else {
            n = (int)(a - 1.);
            w = 1.0;
            for (i = 0; i < n; i++) {
                a -= 1.0;
                w *= a/(1. + (a/b));
            }
            return (log(w) - n*log(b)) + (gamln(a) + algdiv(a, b));
        }
    }
    n = (int)(b - 1.);
    z = 1.0;
    for (i = 0; i < n; i++) {
        b -= 1.0;
        z *= b / (a + b);
    }
    return w + log(z) + (gamln(a) + gamln(b) - gsumln(a, b));
}


double bfrac(double a, double b, double x, double y, double lmbda, double eps)
{
    //    Continued fraction expansion for ix(a,b) when a,b > 1.
    //    It is assumed that  lambda = (a + b)*y - b.

    double alpha, beta, e, r0, t, w, result;
    double c = 1. + lmbda;
    double c0 = b / a;
    double c1 = 1. + (1. / a);
    double yp1 = y + 1.;
    double n = 0.;
    double p = 1.;
    double s = a + 1.;
    double an = 0.;
    double bn = 1.;
    double anp1 = 1.;
    double bnp1 = c / c1;
    double r = c1 / c;

    result = brcomp(a, b, x, y);

    if (result == 0.0) { return 0; }

    // Continued fraction calculation
    while (1) {
        n += 1.0;
        t = n / a;
        w = n * (b - n)*x;
        e = a / s;
        alpha = (p*(p + c0)*e*e) * (w*x);
        e = (1. + t) / (c1 + t + t);
        beta = n + (w / s) + e*(c + n*yp1);
        p = 1. + t;
        s += 2.;
        // Update an, bn, anp1, and bnp1
        t = alpha*an + beta*anp1;
        an = anp1;
        anp1 = t;
        t = alpha*bn + beta*bnp1;
        bn = bnp1;
        bnp1 = t;
        r0 = r;
        r = anp1 / bnp1;
        if (!(fabs(r - r0) > eps*r)) { break; }
        // Rescale an, bn, anp1, and bnp1
        an /= bnp1;
        bn /= bnp1;
        anp1 = r;
        bnp1 = 1.0;
    }
    return result*r;
}


struct TupleDI bgrat(double a, double b, double x , double y, double w, double eps)
{
//    Asymptotic expansion for ix(a,b) when a is larger than b.
//    The result of the expansion is added to w. It is assumed
//    that a >= 15 And b <= 1.  Eps is the tolerance used.
//    Ierr is a variable that reports the status of the results.

    double bp2n, cn, coef, dj, j, l, n2, q, r, s, ssum, t, t2, u, v;
    double c[30] = { 0.0 };
    double d[30] = { 0.0 };
    double bm1 = (b - 0.5) - 0.5;
    double nu = a + bm1*0.5;
    double lnx = (y > 0.375 ? log(x) : alnrel(-y));
    double z = -nu*lnx;
    int i, n;

    if ((b*z) == 0.0) { return (struct TupleDI){.d1 = w, .i1 = 1}; }

    // COMPUTATION OF THE EXPANSION
    // SET R = EXP(-Z)*Z**B/GAMMA(B)
    r = b * (1. + gam1(b)) * exp(b*log(z));
    r *= exp(a*lnx) * exp(0.5*bm1*lnx);
    u = algdiv(b, a) + b*log(nu);
    u = r*exp(-u);
    if (u == 0.0) { return (struct TupleDI){.d1 = w, .i1 = 1}; }

    struct TupleDD ret = grat1(b, z, r, eps);
    q = ret.d2;
    v = 0.25 * pow((1 / nu), 2);
    t2 = 0.25*lnx*lnx;
    l = w / u;
    j = q / r;
    ssum = j;
    t = 1.;
    cn = 1.;
    n2 = 0.;

    for (n = 1; n <= 30; n++) {
        bp2n = b + n2;
        j = (bp2n*(bp2n + 1.)*j + (z + bp2n + 1.)*t)*v;
        n2 += 2.;
        t *= t2;
        cn *= 1/(n2*(n2 + 1.));
        c[n-1] = cn;
        s = 0.0;
        if (n > 1) {
            coef = b - n;
            for (i = 1; i < n; i++) {
                s += coef*c[i-1]*d[n-i-1];
                coef += b;
            }
        }
        d[n-1] = bm1*cn + s/n;
        dj = d[n-1]*j;
        ssum += dj;
        if (ssum <= 0.) { return (struct TupleDI){.d1 = w, .i1 = 1}; }
        if (!(fabs(dj) > eps*(ssum+l))) { break; }
    }
    return (struct TupleDI){.d1 = w + u*ssum, .i1 = 0};
}


double bpser(double a, double b, double x, double eps)
{
    //    Power series expansion for evaluating ix(a,b) when b <= 1
    //    Or b*x <= 0.7.  Eps is the tolerance used.

    double a0, apb, b0, c, n, ssum, t, tol, u, w, z;
    int i, m;
    double result = 0.0;

    if (x == 0.0) { return 0.0; }

    // Compute the factor x**a/(a*beta(a,b))
    a0 = fmin(a, b);

    if (a0 < 1.0) {
        b0 = fmax(a, b);
        if (b0 <= 1.0) {
            result = pow(x, a);
            if (result == 0.0) { return 0.0; }

            apb = a + b;
            if (apb > 1.0) {
                u = a + b - 1.0;
                z = (1. + gam1(u)) / apb;
            } else {
                z = 1. + gam1(apb);
            }
            c = (1. + gam1(a)) * (1. + gam1(b)) / z;
            result *= c * (b / apb);

        } else if (b0 < 8) {
            u = gamln1(a0);
            m = (int)(b0 - 1.0);
            if (m > 0) {
                c = 1.0;
                for (i = 0; i < m; i++) {
                    b0 -= 1.0;
                    c *= (b0 / (a0 + b0));
                }
                u += log(c);
            }
            z = a*log(x) - u;
            b0 -= 1.0;
            apb = a0 + b0;
            if (apb > 1.0) {
                u = a0 + b0 - 1.0;
                t = (1. + gam1(u)) / apb;
            } else {
                t = 1. + gam1(apb);
            }
            result = exp(z) * (a0 / a) * (1. + gam1(b0)) / t;

        } else if (b0 >= 8) {
            u = gamln1(a0) + algdiv(a0, b0);
            z = a*log(x) - u;
            result = (a0 / a) * exp(z);
        }
    } else {
        z = a*log(x) - betaln(a, b);
        result = exp(z) / a;
    }
    if ((result == 0.) || (a <= 0.1*eps)) { return result; }
    // Compute the series
    ssum = 0.0;
    n = 0;
    c = 1.0;
    tol = eps / a;
    while (1) {
        n += 1.0;
        c *= (0.5 + (0.5 - b/n))*x;
        w = c / (a + n);
        ssum += w;
        if (!(fabs(w) > tol)) { break; }
    }
    return result * (1. + a*ssum);
}


struct TupleDDI bratio(double a, double b, double x, double y)
{
    //           Evaluation of the incomplete beta function Ix(a,b)
    //
    //                    --------------------
    //
    //    It is assumed that a and b are nonnegative, and that x <= 1
    //    And y = 1 - x.  Bratio assigns w and w1 the values
    //
    //                     w  = Ix(a,b)
    //                     w1 = 1 - Ix(a,b)
    //
    //    Ierr is a variable that reports the status of the results.
    //    If no input errors are detected then ierr is set to 0 and
    //    w and w1 are computed. Otherwise, if an error is detected,
    //    then w and w1 are assigned the value 0 and ierr is set to
    //    one of the following values ...
    //
    //       Ierr = 1  if a or b is negative
    //       ierr = 2  if a = b = 0
    //       ierr = 3  if x < 0 Or x > 1
    //       Ierr = 4  if y < 0 Or y > 1
    //       Ierr = 5  if x + y .ne. 1
    //       Ierr = 6  if x = a = 0
    //       ierr = 7  if y = b = 0
    //

    double a0, b0, lmbda, x0, y0;
    int ind, n;
    double w = 0.;
    double w1 = 0.;
    double eps = fmax(spmpar[0], 1e-15);
    struct TupleDI bgratret;

    if ((a < 0.) || (b < 0.)) { return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 1}; }
    else if ((a == 0.) && (b == 0.)) { return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 2}; }
    else if ((x < 0.) || (x > 1.)) { return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 3}; }
    else if ((y < 0.) || (y > 1.)) { return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 4}; }
    else if (fabs(((x+y)-0.5) - 0.5) > 3.*eps) { return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 5}; }
    else if (x == 0.) {
        if (a == 0) {
            return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 6};
        } else {
            return (struct TupleDDI){.d1 = 0.0, .d2 = 1.0, .i1 = 0};
        }
    } else if (y == 0.) {
        if (b == 0) {
            return (struct TupleDDI){.d1 = w, .d2 = w1, .i1 = 7};
        } else {
            return (struct TupleDDI){.d1 = 1.0, .d2 = 0.0, .i1 = 0};
        }
    } else if (a == 0.) { return (struct TupleDDI){.d1 = 1.0, .d2 = 0.0, .i1 = 0}; }
    else if (b == 0.) { return (struct TupleDDI){.d1 = 0.0, .d2 = 1.0, .i1 = 0}; }
    else if (fmax(a, b) < 1e-3*eps) { return (struct TupleDDI){.d1 = b/(a+b), .d2 = a/(a+b), .i1 = 0}; }

    ind = 0;
    a0 = a;
    b0 = b;
    x0 = x;
    y0 = y;

    if (fmin(a0, b0) <= 1.0) {
        if (!(x <= 0.5)){
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
        }
        if (b0 < fmin(eps, eps*a0)) {
            w = fpser(a0, b0, x0, eps);
            w1 = 0.5 + (0.5 - w);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
        }
        if ((a0 < fmin(eps, eps*b0)) && (b0*x0 <= 1.)) {
            w1 = apser(a0, b0, x0, eps);
            w = 0.5 + (0.5 - w1);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
        }
        if (fmax(a0, b0) <= 1.0) {
            if ((pow(x0, a0) <= 0.9) || (a0 >= fmin(0.2, b0))) {
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            if (x0 >= 0.3) {
                w1 = bpser(b0, a0, y0, eps);
                w = 0.5 + (0.5 - w1);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            // 140
            n = 20;
            w1 = bup(b0, a0, y0, x0, n, eps);

            b0 += n;
            bgratret = bgrat(b0, a0, y0, x0, w1, 15.*eps);
            w1 = bgratret.d1;
            w = 0.5 + (0.5 - w1);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};

        } else {
            // 20
            if (b0 <= 1.0) {
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            if (x0 >= 0.3) {
                w1 = bpser(b0, a0, y0, eps);
                w = 0.5 + (0.5 - w1);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            if (x0 >= 0.1) {
                // 30
                if (b0 > 15.) {
                    bgratret = bgrat(b0, a0, y0, x0, w1, 15.*eps);
                    w1 = bgratret.d1;
                    w = 0.5 + (0.5 - w1);
                    return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
                } else {
                    // 140
                    n = 20;
                    w1 = bup(b0, a0, y0, x0, n, eps);
                    b0 += n;
                    bgratret = bgrat(b0, a0, y0, x0, w1, 15.*eps);
                    w1 = bgratret.d1;
                    w = 0.5 + (0.5 - w1);
                    return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
                }
            }
            if (pow((x0*b0),a0) <= 0.7) {
                w = bpser(a0, b0, x0, eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            if (b0 > 15.) {
                bgratret = bgrat(b0, a0, y0, x0, w1, 15.*eps);
                w1 = bgratret.d1;
                w = 0.5 + (0.5 - w1);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
            // 140
            n = 20;
            w1 = bup(b0, a0, y0, x0, n, eps);
            b0 += n;
            bgratret = bgrat(b0, a0, y0, x0, w1, 15.*eps);
            w1 = bgratret.d1;
            w = 0.5 + (0.5 - w1);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
        }
    } else {
        lmbda = (a > b ? ((a + b)*y - b) : (a - (a + b)*x));
        if (lmbda < 0.0) {
            ind = 1;
            a0 = b;
            b0 = a;
            x0 = y;
            y0 = x;
            lmbda = fabs(lmbda);
        }
        if ((b0 < 40.) && (b0*x0 <= 0.7)) {
            w = bpser(a0, b0, x0, eps);
            w1 = 0.5 + (0.5 - w);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
        }
        if (b0 < 40.) {
            n = (int)b0;
            b0 -= n;
            if (b0 == 0.) {
                n -= 1;
                b0 = 1.0;
            }
            w = bup(b0, a0, y0, x0, n, eps);
            if (x0 <= 0.7) {
                w += bpser(a0,b0,x0,eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            } else {
                if (a0 <= 15.) {
                    n = 20;
                    w += bup(a0, b0, x0, y0, n, eps);
                    a0 += n;
                }
                bgratret = bgrat(a0, b0, x0, y0, w, 15.*eps);
                w = bgratret.d1;
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
        }
        if (a0 > b0) {
            if ((b0 <= 100) || (lmbda > 0.03*b0)) {
                w = bfrac(a0, b0, x0, y0, lmbda, 15.0*eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            } else {
                // 200
                w = basym(a0, b0, lmbda, 100.*eps);
                w1 = 0.5 + (0.5 - w);
                return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
            }
        }
        if ((a0 <= 100.) || (lmbda > 0.03*a0)) {
            w = bfrac(a0, b0, x0, y0, lmbda, 15.0*eps);
            w1 = 0.5 + (0.5 - w);
            return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
        }
        // 200
        w = basym(a0, b0, lmbda, 100.*eps);
        w1 = 0.5 + (0.5 - w);
        return (struct TupleDDI){.d1 = (ind == 0? w : w1), .d2 = (ind == 0? w1 : w), .i1 = 0};
    }
}


double brcmp1(int mu, double a, double b, double x, double y)
{
    //    Evaluation of  exp(mu) * (x**a*y**b/beta(a,b))

    double a0, apb, b0, c, e, h, lmbda, t, u, v, x0, y0, z;
    double lnx = 1., lny = 1.;
    int i, n;
    double const r2pi = 1./sqrt(2.*PI);

    a0 = fmin(a, b);
    if (a0 >= 8.0) {
        if (a > b) {
            h = b / a;
            x0 = 1. / (1. + h);
            y0 = h / (1. + h);
            lmbda = (a + b)*y - b;
        } else {
            h = a / b;
            x0 = h / (1. + h);
            y0 = 1. / (1. + h);
            lmbda = a - (a + b)*x;
        }
        e = -lmbda / a;
        if (fabs(e) > 0.6) {
            u = e - log(x / x0);
        } else {
            u = rlog1(e);
        }
        e = lmbda / b;
        if (fabs(e) > 0.6) {
            v = e - log(y / y0);
        } else {
            v = rlog1(e);
        }
        z = esum(mu, -(a*u + b*v));
        return r2pi*sqrt(b*x0)*z*exp(-bcorr(a, b));
    }
    if ((x > 0.375) && (y > 0.375)) {
        lnx = log(x);
        lny = log(y);
    } else if ((x > 0.375) && (y <= 0.375)) {
        lnx = alnrel(-y);
        lny = log(y);
    } else if (x <= 0.375) {
        lnx = log(x);
        lny = alnrel(-x);
    }
    z = a*lnx + b*lny;

    if (a0 >= 1.0) {
        z -= betaln(a, b);
        return esum(mu, z);
    }

    b0 = fmax(a, b);
    if (b0 >= 8.0) {
        u = gamln1(a0) + algdiv(a0, b0);
        return a0*esum(mu, z - u);
    }

    if (b0 > 1.0) {
        u = gamln1(a0);
        n = (int)(b0 - 1.);
        if (n >= 1) {
            c = 1.;
            for (i = 0; i < n; i++) {
                b0 -= 1.;
                c *= b0 / (a0 + b0);
            }
            u += log(c);
        }
        z -= u;
        b0 -= 1.;
        apb = a0 + b0;

        if (apb > 1.0) {
            u = a0 + b0 - 1.;
            t = (1. + gam1(u)) / apb;
        } else {
            t = 1. + gam1(apb);
        }
        return a0*esum(mu, z)*(1. + gam1(b0)) / t;
    }
    if (esum(mu, z) == 0.0) {
        return 0.0;
    }
    apb = a + b;
    t = exp(z);
    if (apb > 1.0) {
        u = a + b - 1.0;
        z = (1. + gam1(u)) / apb;
    } else {
        z = 1. + gam1(apb);
    }

    c = (1. + gam1(a)) * (1. + gam1(b)) / z;
    return t*(a0*c) / (1. + a0 / b0);
}


double brcomp(double a, double b, double x, double y)
{
    //    Evaluation of x**a*y**b/beta(a,b)

    double a0, apb, b0, c, e, h, lmbda, lnx, lny, t, u, v, x0, y0, z;
    double const r2pi = 1. / sqrt(2 * PI);
    int i, n;

    if ((x == 0.) || (y == 0.)) {
        return 0.;
    }
    a0 = fmin(a, b);

    if (a0 >= 8.) {
        if (a > b) {
            h = b / a;
            x0 = 1. / (1. + h);
            y0 = h / (1. + h);
            lmbda = (a + b)*y - b;
        } else {
            h = a / b;
            x0 = h / (1. + h);
            y0 = 1. / (1. + h);
            lmbda = a - (a + b)*x;
        }
        e = -lmbda / a;
        if (fabs(e) > 0.6) {
            u = e - log(x / x0);
        } else {
            u = rlog1(e);
        }
        e = lmbda / b;
        if (fabs(e) > 0.6) {
            v = e - log(y / y0);
        } else {
            v = rlog1(e);
        }

        z = exp(-(a*u + b*v));
        return r2pi*sqrt(b*x0)*z*exp(-bcorr(a, b));
    }
    if (x <= 0.375) {
        lnx = log(x);
        lny = alnrel(-x);
    } else {
        lnx = (y > 0.375 ? log(x) : alnrel(-y));
        lny = log(y);
    }
    z = a*lnx + b*lny;
    if (a0 >= 1.0) {
        z -= betaln(a, b);
        return exp(z);
    }

    b0 = fmax(a, b);
    if (b0 >= 8.0) {
        u = gamln1(a0) + algdiv(a0, b0);
        return a0*exp(z - u);
    }

    if (b0 > 1.0) {
        u = gamln1(a0);
        n = (int)(b0 - 1.);
        if (n >= 1) {
            c = 1.0;
            for (i = 0; i < n; i++) {
                b0 -= 1.0;
                c *= b0 / (a0 + b0);
            }
            u += log(c);
        }
        z -= u;
        b0 -= 1.0;
        apb = a0 + b0;

        if (apb > 1.0) {
            u = a0 + b0 - 1.0;
            t = (1. + gam1(u)) / apb;
        } else {
            t = 1. + gam1(apb);
        }
        return a0*exp(z)*(1. + gam1(b0)) / t;
    }
    if (exp(z) == 0.) {
        return 0.0;
    }

    apb = a + b;
    t = exp(z);
    if (apb > 1.0) {
        u = a + b - 1.0;
        z = (1. + gam1(u)) / apb;
    } else {
        z = 1. + gam1(apb);
    }
    c = (1. + gam1(a)) * (1. + gam1(b)) / z;
    return t * (a0*c) / (1. + a0 / b0);
}


double bup(double a, double b, double x, double y, int n, double eps)
{
    //     Evaluation of Ix(a,b) - Ix(a+n,b) where n is a positive integer and
    //     eps is the tolerance used.
    double apb = a + b;
    double ap1 = a + 1.;
    double d = 1.;
    double r, t, w, result;
    int i, nm1;
    int k = 0;
    int mu = 0;

    //          Obtain the scaling factor exp(-mu) and
    //             exp(mu)*(x**a*y**b/beta(a,b))/a

    if (!((n == 1) || (a < 1) || (apb < 1.1*ap1))) {
        mu = 708;
        t = mu;
        d = exp(-708);
    }

    result = brcmp1(mu, a, b, x, y) / a;
    if ((n == 1) || (result == 0.)) {
        return result;
    }
    nm1 = n - 1;
    w = d;

    // Let k be the index of the maximum term
    k = 0;
    if (b <= 1.) {
        // 50
        for (i = 0; i < n-1; i++) {
            d *= ((apb + i) / (ap1 + i))*x;
            w += d;
            if (d <= eps*w) { break; }
        }
        return result*w;
    }
    if (y > 1e-4) {
        // 20
        r = (b - 1.)*x/y - a;
        if (r < 1.) {
            // 50
            for (i = 0; i < n-1; i++) {
                d *= ((apb + i) / (ap1 + i))*x;
                w += d;
                if (d <= eps*w) { break; }
            }
            return result*w;
        }
        k = nm1;
        t = nm1;
        if (r < t) {
            k = (int)r;
        }
    } else {
        k = nm1;
    }
    // 30
    // Add the increasing terms of the series
    for (i = 0; i < k; i++) {
        d *= ((apb + i) / (ap1 + i))*x;
        w += d;
    }

    if (k == nm1) {
        return result*w;
    }

    // 50
    // Add the remaining terms of the series
    for (i = k; i < n-1; i++) {
        d *= ((apb + i) / (ap1 + i))*x;
        w += d;
        if (d <= eps*w) { break; }
    }
    return result*w;
}

    //               Cumulative Distribution Function
    //                         BETa Distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the beta distribution given
    //     values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next four argument
    //               values is to be calculated from the others.
    //               Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from X,Y,A and B
    //               iwhich = 2 : Calculate X and Y from P,Q,A and B
    //               iwhich = 3 : Calculate A from P,Q,X,Y and B
    //               iwhich = 4 : Calculate B from P,Q,X,Y and A
    //
    //                    INTEGER WHICH
    //
    //     P <--> The integral from 0 to X of the chi-square
    //            distribution.
    //            Input range: [0, 1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: [0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //     X <--> Upper limit of integration of beta density.
    //            Input range: [0,1].
    //            Search range: [0,1]
    //                    DOUBLE PRECISION X
    //
    //     Y <--> 1-X.
    //            Input range: [0,1].
    //            Search range: [0,1]
    //            X + Y = 1.0.
    //                    DOUBLE PRECISION Y
    //
    //     A <--> The first parameter of the beta density.
    //            Input range: (0, +infinity).
    //            Search range: [1D-100,1D100]
    //                    DOUBLE PRECISION A
    //
    //     B <--> The second parameter of the beta density.
    //            Input range: (0, +infinity).
    //            Search range: [1D-100,1D100]
    //                    DOUBLE PRECISION B
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                4 if X + Y .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Cumulative distribution function  (P)  is calculated directly by
    //     code associated with the following reference.
    //
    //     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
    //     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
    //     Trans. Math.  Softw. 18 (1993), 360-373.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //                              Note
    //
    //
    //     The beta density is proportional to
    //               t^(A-1) * (1-t)^(B-1)
    //
    //**********************************************************************

struct TupleDDID cdfbet_which1(double x, double y, double a, double b)
{
    struct TupleDD betret;
    struct TupleDDID ret = {0};
    if (!((0 <= x) && (x <= 1.0))) {
        ret.i1 = -1;
        ret.d3 = (!(x > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= y) && (y <= 1.0))) {
        ret.i1 = -2;
        ret.d3 = (!(y > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < a)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    if (!(0 < b)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4, .d3 = 0.0};
    }
    if (((fabs(x+y)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 4,
                                  .d3 = (x+y < 0 ? 0.0 : 1.0)};
    }

    betret = cumbet(x, y, a, b);
    return (struct TupleDDID){.d1 = betret.d1, .d2 = betret.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDDID cdfbet_which2(double p, double q, double a, double b)
{
    double xx, yy;
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DzrorState DZ = {0};

    DZ.xlo = 0.0;
    DZ.xhi = 1.0;
    DZ.atol = atol;
    DZ.rtol = tol;
    DZ.x = 0.0;
    DZ.b = 0.0;
    struct TupleDD betret;
    struct TupleDDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d3 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d3 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < a)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    if (!(0 < b)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4, .d3 = 0.0};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 3,
                                  .d3 = (p+q < 0 ? 0.0 : 1.0)};
    }

    if (qporq) {
        dzror(&DZ);
        yy = 1.0 - DZ.x;
        while (DZ.status == 1) {
            betret = cumbet(DZ.x, yy, a, b);
            DZ.fx = betret.d1 - p;
            dzror(&DZ);
            yy = 1.0 - DZ.x;
        }
        xx = DZ.x;
    } else {
        dzror(&DZ);
        xx = 1.0 - DZ.x;
        while (DZ.status == 1) {
            betret = cumbet(xx, DZ.x, a, b);
            DZ.fx = betret.d2 - q;
            dzror(&DZ);
            xx = 1.0 - DZ.x;
        }
        yy = DZ.x;
    }

    if (DZ.status == -1) {
        ret.d1 = xx;
        ret.d2 = yy;
        ret.i1 = (DZ.qleft ? 1 : 2);
        ret.d3 = (DZ.qleft ? 0 : 1);
        return ret;
    } else {
        ret.d1 = xx;
        ret.d2 = yy;
        return ret;
    }
}


struct TupleDID cdfbet_which3(double p, double q, double x, double y, double b)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD betret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= x) && (x <= 1.0))) {
        ret.i1 = -3;
        ret.d2 = (!(x > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= y) && (y <= 1.0))) {
        ret.i1 = -4;
        ret.d2 = (!(y > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < b)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3, .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(x+y)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4, .d2 = (x+y < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        betret = cumbet(x, y, DS.x, b);
        DS.fx = (qporq ? betret.d1 - p : betret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfbet_which4(double p, double q, double x, double y, double a)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD betret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= x) && (x <= 1.0))) {
        ret.i1 = -3;
        ret.d2 = (!(x > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= y) && (y <= 1.0))) {
        ret.i1 = -4;
        ret.d2 = (!(y > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < a)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3, .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(x+y)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4, .d2 = (x+y < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        betret = cumbet(x, y, a, DS.x);
        DS.fx = (qporq ? betret.d1 - p : betret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //                         BINomial distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the binomial
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next four argument
    //               values is to be calculated from the others.
    //               Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
    //               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
    //               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
    //               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
    //                    INTEGER WHICH
    //
    //     P <--> The cumulation from 0 to S of the binomial distribution.
    //            (Probablility of S or fewer successes in XN trials each
    //            with probability of success PR.)
    //            Input range: [0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: [0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //     S <--> The number of successes observed.
    //            Input range: [0, XN]
    //            Search range: [0, XN]
    //                    DOUBLE PRECISION S
    //
    //     XN  <--> The number of binomial trials.
    //              Input range: (0, +infinity).
    //              Search range: [1E-100, 1E100]
    //                    DOUBLE PRECISION XN
    //
    //     PR  <--> The probability of success in each binomial trial.
    //              Input range: [0,1].
    //              Search range: [0,1]
    //                    DOUBLE PRECISION PR
    //
    //     OMPR  <--> 1-PR
    //              Input range: [0,1].
    //              Search range: [0,1]
    //              PR + OMPR = 1.0
    //                    DOUBLE PRECISION OMPR
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                4 if PR + OMPR .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
    //     Mathematical   Functions (1966) is   used  to reduce the  binomial
    //     distribution  to  the  cumulative incomplete    beta distribution.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //**********************************************************************

struct TupleDDID cdfbin_which1(double s, double xn, double pr, double ompr)
{
    if (!((0 <= s) && (s <= xn))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 < xn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!((0 <= pr) && (pr < 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3,
                                  .d3 = (pr > 0.0 ? 0.0 : 1.0)};
    }
    if (!((0 <= ompr) && (ompr < 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4,
                                  .d3 = (ompr > 0.0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 4,
                                  .d3 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }
    struct TupleDD res = cumbin(s, xn, pr, ompr);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfbin_which2(double p, double q, double xn, double pr, double ompr)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.0;
    DS.big = xn;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = xn / 2.0;
    struct TupleDD binret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < xn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!((0 <= pr) && (pr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -4, .d2 = (pr > 0.0 ? 0.0 : 1.0)};
    }
    if (!((0 <= ompr) && (ompr < 1))) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -5, .d2 = (ompr > 0.0 ? 0.0 : 1.0)};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3, .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4, .d2 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }
    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        binret = cumbin(DS.x, xn, pr, ompr);
        DS.fx = (qporq ? binret.d1 - p : binret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : xn);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfbin_which3(double p, double q, double s, double pr, double ompr)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD binret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 <= s)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!((0 <= pr) && (pr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -4, .d2 = (pr > 0.0 ? 0.0 : 1.0)};
    }
    if (!((0 <= ompr) && (ompr < 1))) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -5, .d2 = (ompr > 0.0 ? 0.0 : 1.0)};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3, .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4, .d2 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        binret = cumbin(s, DS.x, pr, ompr);
        DS.fx = (qporq ? binret.d1 - p : binret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDDID cdfbin_which4(double p, double q, double s, double xn)
{
    double pr, ompr;
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DzrorState DZ = {0};

    DZ.xlo = 0.0;
    DZ.xhi = 1.0;
    DZ.atol = atol;
    DZ.rtol = tol;
    DZ.x = 0.0;
    DZ.b = 0.0;
    struct TupleDD binret;
    struct TupleDDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d3 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d3 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= s) && (s <= xn))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    if (!(0 < xn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4, .d3 = 0.0};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 3,
                                  .d3 = (p+q < 0 ? 0.0 : 1.0)};
    }

    if (qporq) {
        dzror(&DZ);
        ompr = 1.0 - DZ.x;
        while (DZ.status == 1) {
            binret = cumbin(s, xn, DZ.x, ompr);
            DZ.fx = binret.d1 - p;
            dzror(&DZ);
            ompr = 1.0 - DZ.x;
        }
        pr = DZ.x;
    } else {
        dzror(&DZ);
        pr = 1.0 - DZ.x;
        while (DZ.status == 1) {
            binret = cumbin(s, xn, pr, DZ.x);
            DZ.fx = binret.d2 - q;
            dzror(&DZ);
            pr = 1.0 - DZ.x;
        }
        ompr = DZ.x;
    }

    if (DZ.status == -1) {
        ret.d1 = pr;
        ret.d2 = ompr;
        ret.i1 = (DZ.qleft ? 1 : 2);
        ret.d3 = (DZ.qleft ? 0 : 1);
        return ret;
    } else {
        ret.d1 = pr;
        ret.d2 = ompr;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //               CHI-Square distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the chi-square
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next three argument
    //               values is to be calculated from the others.
    //               Legal range: 1..3
    //               iwhich = 1 : Calculate P and Q from X and DF
    //               iwhich = 2 : Calculate X from P,Q and DF
    //               iwhich = 3 : Calculate DF from P,Q and X
    //                    INTEGER WHICH
    //
    //     P <--> The integral from 0 to X of the chi-square
    //            distribution.
    //            Input range: [0, 1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //     X <--> Upper limit of integration of the non-central
    //            chi-square distribution.
    //            Input range: [0, +infinity).
    //            Search range: [0,1E100]
    //                    DOUBLE PRECISION X
    //
    //     DF <--> Degrees of freedom of the
    //             chi-square distribution.
    //             Input range: (0, +infinity).
    //             Search range: [ 1E-100, 1E100]
    //                    DOUBLE PRECISION DF
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //               10 indicates error returned from cumgam.  See
    //                  references in cdfgam
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
    //     Mathematical Functions   (1966) is used   to reduce the chisqure
    //     distribution to the incomplete distribution.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //**********************************************************************

struct TupleDDID cdfchi_which1(double x, double df)
{
    struct TupleDD chiret;
    if (!(0 <= x)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 <= df)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    chiret = cumchi(x, df);
    return (struct TupleDDID){.d1 = chiret.d1, .d2 = chiret.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfchi_which2(double p, double q, double df)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    int porq = (qporq ? p : q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.0;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD chiret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 <= df)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        chiret = cumchi(DS.x, df);
        DS.fx = (qporq ? chiret.d1 - p : chiret.d2 - q);
        if (DS.fx + porq <= 1.5) {
            return (struct TupleDID){.d1 = DS.x, .i1 = 10, .d2 = 0.0};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfchi_which3(double p, double q, double x)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    int porq = (qporq ? p : q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD chiret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 <= x)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        chiret = cumchi(x, DS.x);
        DS.fx = (qporq ? chiret.d1 - p : chiret.d2 - q);
        if (DS.fx + porq > 1.5) {
            return (struct TupleDID){.d1 = DS.x, .i1 = 10, .d2 = 0.0};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //               Non-central Chi-Square
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the non-central chi-square
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next three argument
    //               values is to be calculated from the others.
    //               Input range: 1..4
    //               iwhich = 1 : Calculate P and Q from X and DF
    //               iwhich = 2 : Calculate X from P,DF and PNONC
    //               iwhich = 3 : Calculate DF from P,X and PNONC
    //               iwhich = 3 : Calculate PNONC from P,X and DF
    //                    INTEGER WHICH
    //
    //     P <--> The integral from 0 to X of the non-central chi-square
    //            distribution.
    //            Input range: [0, 1-1E-16).
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Q is not used by this subroutine and is only included
    //            for similarity with other cdf* routines.
    //                    DOUBLE PRECISION Q
    //
    //     X <--> Upper limit of integration of the non-central
    //            chi-square distribution.
    //            Input range: [0, +infinity).
    //            Search range: [0,1E300]
    //                    DOUBLE PRECISION X
    //
    //     DF <--> Degrees of freedom of the non-central
    //             chi-square distribution.
    //             Input range: (0, +infinity).
    //             Search range: [ 1E-300, 1E300]
    //                    DOUBLE PRECISION DF
    //
    //     PNONC <--> Non-centrality parameter of the non-central
    //                chi-square distribution.
    //                Input range: [0, +infinity).
    //                Search range: [0,1E4]
    //                    DOUBLE PRECISION PNONC
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of
    //     Mathematical  Functions (1966) is used to compute the cumulative
    //     distribution function.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //                            WARNING
    //
    //     The computation time  required for this  routine is proportional
    //     to the noncentrality  parameter  (PNONC).  Very large  values of
    //     this parameter can consume immense  computer resources.  This is
    //     why the search range is bounded by 1e9.
    //
    //**********************************************************************

struct TupleDDID cdfchn_which1(double x, double df, double pnonc)
{
    x = fmin(x, spmpar[2]);
    df = fmin(df, spmpar[2]);
    pnonc = fmin(pnonc, 1e9);
    struct TupleDD chnret;
    if (!(0 <= x)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 <= df)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!(0 <= pnonc)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }

    chnret = cumchn(x, df, pnonc);
    return (struct TupleDDID){.d1 = chnret.d1, .d2 = chnret.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfchn_which2(double p, double df, double pnonc)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e300;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD chnret;
    struct TupleDID ret = {0};

    df = fmin(df, spmpar[2]);
    pnonc = fmin(pnonc, 1.e9);

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= df)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= pnonc)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        chnret = cumchn(DS.x, df, pnonc);
        DS.fx = chnret.d1 - p;
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e300);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfchn_which3(double p, double x, double pnonc)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-300;
    DS.big = 1e300;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD chnret;
    struct TupleDID ret = {0};

    x = fmin(x, spmpar[2]);
    pnonc = fmin(pnonc, 1.e9);

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= x)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= pnonc)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        chnret = cumchn(x, DS.x, pnonc);
        DS.fx = chnret.d1 - p;
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-300 : 1e300);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfchn_which4(double p, double x, double df)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.0;
    DS.big = 1e9;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD chnret;
    struct TupleDID ret = {0};

    x = fmin(x, spmpar[2]);
    df = fmin(df, spmpar[2]);

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= x)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= df)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        chnret = cumchn(x, df, DS.x);
        DS.fx = chnret.d1 - p;
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e9);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //               F distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the F distribution
    //     given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next four argument
    //               values is to be calculated from the others.
    //               Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from F,DFN and DFD
    //               iwhich = 2 : Calculate F from P,Q,DFN and DFD
    //               iwhich = 3 : Calculate DFN from P,Q,F and DFD
    //               iwhich = 4 : Calculate DFD from P,Q,F and DFN
    //                    INTEGER WHICH
    //
    //       P <--> The integral from 0 to F of the f-density.
    //              Input range: [0,1].
    //                    DOUBLE PRECISION P
    //
    //       Q <--> 1-P.
    //              Input range: (0, 1].
    //              P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //       F <--> Upper limit of integration of the f-density.
    //              Input range: [0, +infinity).
    //              Search range: [0,1E100]
    //                    DOUBLE PRECISION F
    //
    //     DFN < --> Degrees of freedom of the numerator sum of squares.
    //               Input range: (0, +infinity).
    //               Search range: [ 1E-100, 1E100]
    //                    DOUBLE PRECISION DFN
    //
    //     DFD < --> Degrees of freedom of the denominator sum of squares.
    //               Input range: (0, +infinity).
    //               Search range: [ 1E-100, 1E100]
    //                    DOUBLE PRECISION DFD
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of
    //     Mathematical  Functions (1966) is used to reduce the computation
    //     of the  cumulative  distribution function for the  F  variate to
    //     that of an incomplete beta.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //                              WARNING
    //
    //     The value of the  cumulative  F distribution is  not necessarily
    //     monotone in  either degrees of freedom.  There  thus may  be two
    //     values  that  provide a given CDF  value.   This routine assumes
    //     monotonicity and will find an arbitrary one of the two values.
    //
    //**********************************************************************

struct TupleDDID cdff_which1(double f, double dfn, double dfd)
{
    struct TupleDD fret;
    if (!(0 <= f)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 < dfn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!(0 < dfd)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }

    fret = cumf(f, dfn, dfd);
    return (struct TupleDDID){.d1 = fret.d1, .d2 = fret.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdff_which2(double p, double q, double dfn, double dfd)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD fret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 < dfn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfd)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fret = cumf(DS.x, dfn, dfd);
        DS.fx = (qporq ? fret.d1 - p : fret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdff_which3(double p, double q, double f, double dfd)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD fret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 <= f)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfd)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fret = cumf(f, DS.x, dfd);
        DS.fx = (qporq ? fret.d1 - p : fret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdff_which4(double p, double q, double f, double dfn)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD fret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1.0))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1.0))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(0 <= f)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fret = cumf(f, dfn, DS.x);
        DS.fx = (qporq ? fret.d1 - p : fret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //               Non-central F distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the Non-central F
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next five argument
    //               values is to be calculated from the others.
    //               Legal range: 1..5
    //               iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
    //               iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
    //               iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
    //               iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
    //               iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD
    //                    INTEGER WHICH
    //
    //       P <--> The integral from 0 to F of the non-central f-density.
    //              Input range: [0,1-1E-16).
    //                    DOUBLE PRECISION P
    //
    //       Q <--> 1-P.
    //            Q is not used by this subroutine and is only included
    //            for similarity with other cdf* routines.
    //                    DOUBLE PRECISION Q
    //
    //       F <--> Upper limit of integration of the non-central f-density.
    //              Input range: [0, +infinity).
    //              Search range: [0,1E100]
    //                    DOUBLE PRECISION F
    //
    //     DFN < --> Degrees of freedom of the numerator sum of squares.
    //               Input range: (0, +infinity).
    //               Search range: [ 1E-100, 1E100]
    //                    DOUBLE PRECISION DFN
    //
    //     DFD < --> Degrees of freedom of the denominator sum of squares.
    //               Must be in range: (0, +infinity).
    //               Input range: (0, +infinity).
    //               Search range: [ 1E-100, 1E100]
    //                    DOUBLE PRECISION DFD
    //
    //     PNONC <-> The non-centrality parameter
    //               Input range: [0,infinity)
    //               Search range: [0,1E4]
    //                    DOUBLE PRECISION PHONC
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula  26.6.20   of   Abramowitz   and   Stegun,  Handbook  of
    //     Mathematical  Functions (1966) is used to compute the cumulative
    //     distribution function.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //                            WARNING
    //
    //     The computation time  required for this  routine is proportional
    //     to the noncentrality  parameter  (PNONC).  Very large  values of
    //     this parameter can consume immense  computer resources.  This is
    //     why the search range is bounded by 10,000.
    //
    //                              WARNING
    //
    //     The  value  of the  cumulative  noncentral F distribution is not
    //     necessarily monotone in either degrees  of freedom.  There  thus
    //     may be two values that provide a given  CDF value.  This routine
    //     assumes monotonicity  and will find  an arbitrary one of the two
    //     values.
    //
    //**********************************************************************

struct TupleDDID cdffnc_which1(double f, double dfn, double dfd, double phonc)
{
    struct TupleDDI fncret;
    if (!(0 <= f)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 < dfn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!(0 < dfd)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    if (!(0 <= phonc)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4, .d3 = 0.0};
    }
    fncret = cumfnc(f, dfn, dfd, phonc);
    if (fncret.i1 != 0) {
        return (struct TupleDDID){.d1 = fncret.d1, .d2 = fncret.d2, .i1 = 10, .d3 = 0.0};
    } else {
        return (struct TupleDDID){.d1 = fncret.d1, .d2 = fncret.d2, .i1 = 0, .d3 = 0.0};
    }
}


struct TupleDID cdffnc_which2(double p, double q, double dfn, double dfd, double phonc)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDDI fncret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 < dfn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfd)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!(0 <= phonc)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fncret = cumfnc(DS.x, dfn, dfd, phonc);
        DS.fx = fncret.d1 - p;
        if (fncret.i1 != 0) {
            return (struct TupleDID){.d1 = DS.x, .d2 = 0.0, .i1 = 10};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdffnc_which3(double p, double q, double f, double dfd, double phonc)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDDI fncret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= f)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfd)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!(0 <= phonc)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fncret = cumfnc(f, DS.x, dfd, phonc);
        DS.fx = fncret.d1 - p;
        if (fncret.i1 != 0) {
            return (struct TupleDID){.d1 = DS.x, .d2 = 0.0, .i1 = 10};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdffnc_which4(double p, double q, double f, double dfn, double phonc)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDDI fncret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= f)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!(0 <= phonc)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fncret = cumfnc(f, dfn, DS.x, phonc);
        DS.fx = fncret.d1 - p;
        if (fncret.i1 != 0) {
            return (struct TupleDID){.d1 = DS.x, .d2 = 0.0, .i1 = 10};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 1e-100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdffnc_which5(double p, double q, double f, double dfn, double dfd)
{
    double tol = 1e-10;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e4;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDDI fncret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(0 <= f)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < dfn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!(0 < dfd)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        fncret = cumfnc(f, dfn, dfd, DS.x);
        DS.fx = fncret.d1 - p;
        if (fncret.i1 != 0) {
            return (struct TupleDID){.d1 = DS.x, .d2 = 0.0, .i1 = 10};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0 : 1e4);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //                         GAMma Distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the gamma
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next four argument
    //               values is to be calculated from the others.
    //               Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE
    //               iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE
    //               iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE
    //               iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE
    //                    INTEGER WHICH
    //
    //     P <--> The integral from 0 to X of the gamma density.
    //            Input range: [0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //
    //     X <--> The upper limit of integration of the gamma density.
    //            Input range: [0, +infinity).
    //            Search range: [0,1E100]
    //                    DOUBLE PRECISION X
    //
    //     SHAPE <--> The shape parameter of the gamma density.
    //                Input range: (0, +infinity).
    //                Search range: [1E-100,1E100]
    //                  DOUBLE PRECISION SHAPE
    //
    //
    //     SCALE <--> The scale parameter of the gamma density.
    //                Input range: (0, +infinity).
    //                Search range: (1E-100,1E100]
    //                   DOUBLE PRECISION SCALE
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                10 if the gamma or inverse gamma routine cannot
    //                   compute the answer.  Usually happens only for
    //                   X and SHAPE very large (gt 1E10 or more)
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Cumulative distribution function (P) is calculated directly by
    //     the code associated with:
    //
    //     DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
    //     gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
    //     Softw. 12 (1986), 377-393.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //                              Note
    //
    //
    //
    //     The gamma density is proportional to
    //       T**(SHAPE - 1) * EXP(- SCALE * T)
    //
    //
    //**********************************************************************

struct TupleDDID cdfgam_which1(double x, double shape, double scale)
{
    struct TupleDD gamret;
    if (!(0 <= x)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 < shape)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!(0 < scale)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    gamret = cumgam(x*scale, shape);
    if (gamret.d1 >= 1.5) {
        return (struct TupleDDID){.d1 = gamret.d1, .d2 = gamret.d2, .i1 = 10, .d3 = 0.0};
    } else {
        return (struct TupleDDID){.d1 = gamret.d1, .d2 = gamret.d2, .i1 = 0, .d3 = 0.0};
    }
}


struct TupleDID cdfgam_which2(double p, double q, double shape, double scale)
{
    struct TupleDI invret;
    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 < shape)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < scale)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    invret = gaminv(shape, p, q, -1);
    if (invret.i1 < 0) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 10, .d2 = 0.0};
    } else {
        return (struct TupleDID){.d1 = invret.d1/scale, .i1 = 0, .d2 = 0.0};
    }
}


struct TupleDID cdfgam_which3(double p, double q, double x, double scale)
{
    double xscale = x*scale;
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = p <= q;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 1e-100;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.0;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD gamret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= x)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < scale)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        gamret = cumgam(xscale, DS.x);
        DS.fx = (qporq ? gamret.d1 - p : gamret.d2 - q);
        if (((qporq) && (gamret.d1 > 1.5)) || (!(qporq) && (gamret.d2 > 1.5))) {
            return (struct TupleDID){.d1 = DS.x, .d2 = 0.0, .i1 = 10};
        }
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }

}


struct TupleDID cdfgam_which4(double p, double q, double x, double shape)
{
    struct TupleDI invret;

    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= x)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!(0 < shape)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }

    invret = gaminv(shape, p, q, -1);
    if (invret.i1 < 0) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 10, .d2 = 0.0};
    } else {
        return (struct TupleDID){.d1 = invret.d1/x, .i1 = 0, .d2 = 0.0};
    }
}


    //               Cumulative Distribution Function
    //               Negative BiNomial distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the negative binomial
    //     distribution given values for the others.
    //
    //     The  cumulative  negative   binomial  distribution  returns  the
    //     probability that there  will be  F or fewer failures before  the
    //     XNth success in binomial trials each of which has probability of
    //     success PR.
    //
    //     The individual term of the negative binomial is the probability of
    //     S failures before XN successes and is
    //          Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which of the next four argument
    //               values is to be calculated from the others.
    //               Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
    //               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
    //               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
    //               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
    //                    INTEGER WHICH
    //
    //     P <--> The cumulation from 0 to S of the  negative
    //            binomial distribution.
    //            Input range: [0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //     S <--> The upper limit of cumulation of the binomial distribution.
    //            There are F or fewer failures before the XNth success.
    //            Input range: [0, +infinity).
    //            Search range: [0, 1E100]
    //                    DOUBLE PRECISION S
    //
    //     XN  <--> The number of successes.
    //              Input range: [0, +infinity).
    //              Search range: [0, 1E100]
    //                    DOUBLE PRECISION XN
    //
    //     PR  <--> The probability of success in each binomial trial.
    //              Input range: [0,1].
    //              Search range: [0,1].
    //                    DOUBLE PRECISION PR
    //
    //     OMPR  <--> 1-PR
    //              Input range: [0,1].
    //              Search range: [0,1]
    //              PR + OMPR = 1.0
    //                    DOUBLE PRECISION OMPR
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                4 if PR + OMPR .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula   26.5.26   of   Abramowitz  and  Stegun,  Handbook   of
    //     Mathematical Functions (1966) is used  to  reduce calculation of
    //     the cumulative distribution  function to that of  an  incomplete
    //     beta.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //**********************************************************************

struct TupleDDID cdfnbn_which1(double s, double xn, double pr, double ompr)
{
    if (!((0 <= s) && (s <= xn))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(0 <= xn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    if (!((0 <= pr) && (pr < 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3,
                                  .d3 = (pr > 0.0 ? 0.0 : 1.0)};
    }
    if (!((0 <= ompr) && (ompr < 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4,
                                  .d3 = (ompr > 0.0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 4,
                                  .d3 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }
    struct TupleDD res = cumnbn(s, xn, pr, ompr);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfnbn_which2(double p, double q, double xn, double pr, double ompr)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD nbnret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= xn)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!((0 <= pr) && (pr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!((0 <= ompr) && (ompr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4,
                                 .d2 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }
    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        nbnret = cumnbn(DS.x, xn, pr, ompr);
        DS.fx = (qporq ? nbnret.d1 - p : nbnret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfnbn_which3(double p, double q, double s, double pr, double ompr)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD nbnret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2};
    }
    if (!(0 <= s)) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3};
    }
    if (!((0 <= pr) && (pr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4};
    }
    if (!((0 <= ompr) && (ompr <= 1))) {
        return (struct TupleDID){.d1 = 0.0, .d2 = 0.0, .i1 = -5};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 3,
                                 .d2 = (p+q < 0 ? 0.0 : 1.0)};
    }
    if (((fabs(pr+ompr)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDID){.d1 = 0.0, .i1 = 4,
                                 .d2 = (pr+ompr < 0 ? 0.0 : 1.0)};
    }
    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        nbnret = cumnbn(s, DS.x, pr, ompr);
        DS.fx = (qporq ? nbnret.d1 - p : nbnret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDDID cdfnbn_which4(double p, double q, double s, double xn)
{
    double pr, ompr;
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DzrorState DZ = {0};

    DZ.xlo = 0.;
    DZ.xhi = 1.;
    DZ.atol = atol;
    DZ.rtol = tol;
    DZ.x = 0.;
    DZ.b = 0.;

    struct TupleDD nbnret;
    struct TupleDDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1,
                                  .d3 = (p > 0 ? 0.0 : 1.0)};
    }
    if (!((0 < q) && (q <= 1))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2,
                                  .d3 = (1 > 0 ? 0.0 : 1.0)};
    }
    if (!(0 <= s)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    if (!(0 <= xn)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -4, .d3 = 0.0};
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = 4,
                                 .d3 = (p+q < 0 ? 0.0 : 1.0)};
    }

    if (qporq) {
        dzror(&DZ);
        ompr = 1.0 - DZ.x;
        while (DZ.status == 1) {
            nbnret = cumnbn(s, xn, DZ.x, ompr);
            DZ.fx = nbnret.d1 - p;
            dzror(&DZ);
            ompr = 1.0 - DZ.x;
        }
        pr = DZ.x;
    } else {
        dzror(&DZ);
        pr = 1.0 - DZ.x;
        while (DZ.status == 1) {
            nbnret = cumnbn(s, xn, pr, DZ.x);
            DZ.fx = nbnret.d2 - q;
            dzror(&DZ);
            pr = 1.0 - DZ.x;
        }
        ompr = DZ.x;
    }
    if (DZ.status == -1) {
        ret.d1 = pr;
        ret.d2 = ompr;
        ret.i1 = (DZ.qleft ? 1 : 2);
        ret.d3 = (DZ.qleft ? 0 : 1);
        return ret;
    } else {
        ret.d1 = pr;
        ret.d2 = ompr;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //               NORmal distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the normal
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH  --> Integer indicating  which of the  next  parameter
    //     values is to be calculated using values  of the others.
    //     Legal range: 1..4
    //               iwhich = 1 : Calculate P and Q from X,MEAN and SD
    //               iwhich = 2 : Calculate X from P,Q,MEAN and SD
    //               iwhich = 3 : Calculate MEAN from P,Q,X and SD
    //               iwhich = 4 : Calculate SD from P,Q,X and MEAN
    //                    INTEGER WHICH
    //
    //     P <--> The integral from -infinity to X of the normal density.
    //            Input range: (0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //     X < --> Upper limit of integration of the normal-density.
    //             Input range: ( -infinity, +infinity)
    //                    DOUBLE PRECISION X
    //
    //     MEAN <--> The mean of the normal density.
    //               Input range: (-infinity, +infinity)
    //                    DOUBLE PRECISION MEAN
    //
    //     SD <--> Standard Deviation of the normal density.
    //             Input range: (0, +infinity).
    //                    DOUBLE PRECISION SD
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //
    //
    //     A slightly modified version of ANORM from
    //
    //     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
    //     Package of Special Function Routines and Test Drivers"
    //     acm Transactions on Mathematical Software. 19, 22-32.
    //
    //     is used to calculate the cumulative standard normal distribution.
    //
    //     The rational functions from pages  90-95  of Kennedy and Gentle,
    //     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
    //     starting values to Newton's Iterations which compute the inverse
    //     standard normal.  Therefore no  searches  are necessary for  any
    //     parameter.
    //
    //     For X < -15, the asymptotic expansion for the normal is used  as
    //     the starting value in finding the inverse standard normal.
    //     This is formula 26.2.12 of Abramowitz and Stegun.
    //
    //
    //                              Note
    //
    //
    //      The normal density is proportional to
    //      exp( - 0.5 * (( X - MEAN)/SD)**2)
    //
    //
    //**********************************************************************

struct TupleDDID cdfnor_which1(double x, double mean, double sd)
{
    if (!(sd > 0.0)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3, .d3 = 0.0};
    }
    double z = (x - mean)/sd;
    struct TupleDD res = cumnor(z);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfnor_which2(double p, double q, double mean, double sd)
{
    if (!(sd > 0.0)) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -4, .d2 = 0.0};
    }
    double z = dinvnr(p, q);
    return (struct TupleDID){.d1 = sd*z + mean, .i1 = 0, .d2 = 0.0};
}


struct TupleDID cdfnor_which3(double p, double q, double x, double sd)
{
    if (!(sd > 0.0)) {
        return (struct TupleDID){.d1 = 0.0, .i1 = -4, .d2 = 0.0};
    }
    double z = dinvnr(p, q);
    return (struct TupleDID){.d1 = x - sd*z, .i1 = 0, .d2 = 0.0};
}


struct TupleDID cdfnor_which4(double p, double q, double x, double mean)
{
    double z = dinvnr(p, q);
    return (struct TupleDID){.d1 = (x - mean)/z, .i1 = 0, .d2 = 0.0};
}


    //               Cumulative Distribution Function
    //               POIsson distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the Poisson
    //     distribution given values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which  argument
    //               value is to be calculated from the others.
    //               Legal range: 1..3
    //               iwhich = 1 : Calculate P and Q from S and XLAM
    //               iwhich = 2 : Calculate S from P,Q and XLAM
    //               iwhich = 3 : Calculate XLAM from P,Q and S
    //                    INTEGER WHICH
    //
    //        P <--> The cumulation from 0 to S of the poisson density.
    //               Input range: [0,1].
    //                    DOUBLE PRECISION P
    //
    //        Q <--> 1-P.
    //               Input range: (0, 1].
    //               P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //        S <--> Upper limit of cumulation of the Poisson.
    //               Input range: [0, +infinity).
    //               Search range: [0,1E100]
    //                    DOUBLE PRECISION S
    //
    //     XLAM <--> Mean of the Poisson distribution.
    //               Input range: [0, +infinity).
    //               Search range: [0,1E100]
    //                    DOUBLE PRECISION XLAM
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula   26.4.21  of   Abramowitz  and   Stegun,   Handbook  of
    //     Mathematical Functions (1966) is used  to reduce the computation
    //     of  the cumulative distribution function to that  of computing a
    //     chi-square, hence an incomplete gamma function.
    //
    //     Cumulative  distribution function  (P) is  calculated  directly.
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired value of  P.   The  search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //
    //**********************************************************************

struct TupleDDID cdfpoi_which1(double s, double xlam)
{
    if (!(s >= 0.0)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(xlam >= 0.0)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    struct TupleDD res = cumpoi(s, xlam);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdfpoi_which2(double p, double q, double xlam)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.;
    DS.big = 1e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.0;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD poiret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(xlam >= 0.0)) {
        ret.i1 = -3;
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    if ((xlam < 0.01) && (p < 0.975)) {
        // For sufficiently small xlam and p, the result is 0.
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        poiret = cumpoi(DS.x, xlam);
        DS.fx = (qporq ? poiret.d1 - p : poiret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdfpoi_which3(double p, double q, double s)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = 0.0;
    DS.big = 1e300;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.0;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;
    struct TupleDD poiret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 < q) && (q <= 1))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(s >= 0.0)) {
        ret.i1 = -3;
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        poiret = cumpoi(s, DS.x);
        DS.fx = (qporq ? poiret.d1 - p : poiret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0.0 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //                         T distribution
    //
    //
    //                              Function
    //
    //
    //     Calculates any one parameter of the t distribution given
    //     values for the others.
    //
    //
    //                              Arguments
    //
    //
    //     WHICH --> Integer indicating which  argument
    //               values is to be calculated from the others.
    //               Legal range: 1..3
    //               iwhich = 1 : Calculate P and Q from T and DF
    //               iwhich = 2 : Calculate T from P,Q and DF
    //               iwhich = 3 : Calculate DF from P,Q and T
    //                    INTEGER WHICH
    //
    //        P <--> The integral from -infinity to t of the t-density.
    //              Input range: (0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //        T <--> Upper limit of integration of the t-density.
    //               Input range: ( -infinity, +infinity).
    //               Search range: [ -1E100, 1E100 ]
    //                    DOUBLE PRECISION T
    //
    //        DF <--> Degrees of freedom of the t-distribution.
    //                Input range: (0 , +infinity).
    //                Search range: [1e-100, 1E10]
    //                    DOUBLE PRECISION DF
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //
    //                              Method
    //
    //
    //     Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
    //     Mathematical Functions  (1966) is used to reduce the computation
    //     of the cumulative distribution function to that of an incomplete
    //     beta.
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //**********************************************************************

struct TupleDDID cdft_which1(double t, double df)
{
    if (!(df > 0.0)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }
    struct TupleDD res = cumt(t, df);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdft_which2(double p, double q, double df)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = -1.e100;
    DS.big = 1.e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = dt1(p, q, df);
    struct TupleDD tret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!(df > 0.0)) {
        ret.i1 = -3;
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        tret = cumt(DS.x, df);
        DS.fx = (qporq ? tret.d1 - p : tret.d2 - q);
        dinvr(&DS, &DZ);
    }

    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? -1e100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdft_which3(double p, double q, double t)
{
    double tol = 1e-10;
    double atol = 1e-50;
    int qporq = (p <= q);
    DinvrState DS = {0};
    DzrorState DZ = {0};
    struct TupleDD tret;
    struct TupleDID ret = {0};

    DS.small = 1e-100;
    DS.big = 1e10;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.0;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.0;

    if (!((0 <= p) && (p <= 1))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (!((0 <= q) && (q <= 1))) {
        ret.i1 = -2;
        ret.d2 = (!(q > 0.0) ? 0.0 : 1.0);
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        tret = cumt(t, DS.x);
        DS.fx = (qporq ? tret.d1 - p : tret.d2 - q);
        dinvr(&DS, &DZ);
    }
    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? -1e100 : 1e10);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


    //               Cumulative Distribution Function
    //                  Non-Central T distribution
    //
    //                               Function
    //
    //     Calculates any one parameter of the noncentral t distribution give
    //     values for the others.
    //
    //                               Arguments
    //
    //     WHICH --> Integer indicating which  argument
    //               values is to be calculated from the others.
    //               Legal range: 1..3
    //               iwhich = 1 : Calculate P and Q from T,DF,PNONC
    //               iwhich = 2 : Calculate T from P,Q,DF,PNONC
    //               iwhich = 3 : Calculate DF from P,Q,T
    //               iwhich = 4 : Calculate PNONC from P,Q,DF,T
    //                    INTEGER WHICH
    //
    //        P <--> The integral from -infinity to t of the noncentral t-den
    //              Input range: (0,1].
    //                    DOUBLE PRECISION P
    //
    //     Q <--> 1-P.
    //            Input range: (0, 1].
    //            P + Q = 1.0.
    //                    DOUBLE PRECISION Q
    //
    //        T <--> Upper limit of integration of the noncentral t-density.
    //               Input range: ( -infinity, +infinity).
    //               Search range: [ -1E100, 1E100 ]
    //                    DOUBLE PRECISION T
    //
    //        DF <--> Degrees of freedom of the noncentral t-distribution.
    //                Input range: (0 , +infinity).
    //                Search range: [1e-100, 1E10]
    //                    DOUBLE PRECISION DF
    //
    //     PNONC <--> Noncentrality parameter of the noncentral t-distributio
    //                Input range: [-1e6, 1E6].
    //
    //     STATUS <-- 0 if calculation completed correctly
    //               -I if input parameter number I is out of range
    //                1 if answer appears to be lower than lowest
    //                  search bound
    //                2 if answer appears to be higher than greatest
    //                  search bound
    //                3 if P + Q .ne. 1
    //                    INTEGER STATUS
    //
    //     BOUND <-- Undefined if STATUS is 0
    //
    //               Bound exceeded by parameter number I if STATUS
    //               is negative.
    //
    //               Lower search bound if STATUS is 1.
    //
    //               Upper search bound if STATUS is 2.
    //
    //                                Method
    //
    //     Upper tail    of  the  cumulative  noncentral t is calculated usin
    //     formulae  from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou
    //     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
    //
    //     Computation of other parameters involve a search for a value that
    //     produces  the desired  value  of P.   The search relies  on  the
    //     monotinicity of P with the other parameter.
    //
    //***********************************************************************

struct TupleDDID cdftnc_which1(double t, double df, double pnonc)
{

    if (!(t == t)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -1, .d3 = 0.0};
    }
    if (!(df > 0.0)) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -2, .d3 = 0.0};
    }

    df = fmin(df, 1.e10);
    t = fmax(fmin(t, spmpar[2]), -spmpar[2]);

    if (!((-1.e6 <= pnonc) && (pnonc <= 1.e6))) {
        return (struct TupleDDID){.d1 = 0.0, .d2 = 0.0, .i1 = -3,
         .d3 = (pnonc > -1e6 ? 1.0e6 : -1.0e6)};
    }

    struct TupleDD res = cumtnc(t, df, pnonc);
    return (struct TupleDDID){.d1 = res.d1, .d2 = res.d2, .i1 = 0, .d3 = 0.0};
}


struct TupleDID cdftnc_which2(double p, double q, double df, double pnonc)
{

    double tol = 1e-8;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};

    DS.small = -1.e100;
    DS.big = 1.e100;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;
    struct TupleDD tncret;
    struct TupleDID ret = {0};

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }

    if (!(df > 0.)) {
        ret.i1 = -3;
        return ret;
    }
    df = fmin(df, 1.e10);

    if (!((-1e6 <= pnonc) && (pnonc <= 1e6))) {
        ret.i1 = -4;
        ret.d2 = (pnonc > -1e6 ? 1e6 : -1e6);
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        tncret = cumtnc(DS.x, df, pnonc);
        DS.fx = tncret.d1 - p;
        dinvr(&DS, &DZ);
    }
    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? -1e100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdftnc_which3(double p, double q, double t, double pnonc)
{

    double tol = 1e-8;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};
    struct TupleDD tncret;
    struct TupleDID ret = {0};

    DS.small = 1e-100;
    DS.big = 1.e10;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    t = fmax(fmin(t, spmpar[2]), -spmpar[2]);
    if (!(t == t)) {
        ret.i1 = -3;
        return ret;
    }

    if (!((-1e6 <= pnonc) && (pnonc <= 1e6))) {
        ret.i1 = -4;
        ret.d2 = (pnonc > -1e6 ? 1e6 : -1e6);
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        tncret = cumtnc(t, DS.x, pnonc);
        DS.fx = tncret.d1 - p;
        dinvr(&DS, &DZ);
    }
    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? -1e100 : 1e100);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDID cdftnc_which4(double p, double q, double t, double df)
{

    double tol = 1e-8;
    double atol = 1e-50;
    DinvrState DS = {0};
    DzrorState DZ = {0};
    struct TupleDD tncret;
    struct TupleDID ret = {0};

    DS.small = -1.e6;
    DS.big = 1.e6;
    DS.absstp = 0.5;
    DS.relstp = 0.5;
    DS.stpmul = 5.;
    DS.abstol = atol;
    DS.reltol = tol;
    DS.x = 5.;

    if (!((0 <= p) && (p <= (1. - 1e-16)))) {
        ret.i1 = -1;
        ret.d2 = (!(p > 0.0) ? 0.0 : (1. - 1e-16));
        return ret;
    }
    if (!(t == t)) {
        ret.i1 = -3;
        return ret;
    }
    if (!(df > 0.)) {
        ret.i1 = -4;
        return ret;
    }
    if (((fabs(p+q)-0.5)-0.5) > 3*spmpar[0]) {
        ret.i1 = 3;
        ret.d2 = (p+q < 0 ? 0.0 : 1.0);
        return ret;
    }

    t = fmax(fmin(t, spmpar[2]), -spmpar[2]);
    df = fmin(df, 1.e10);

    dinvr(&DS, &DZ);
    while (DS.status == 1) {
        tncret = cumtnc(t, df, DS.x);
        DS.fx = tncret.d1 - p;
        dinvr(&DS, &DZ);
    }
    if (DS.status == -1) {
        ret.d1 = DS.x;
        ret.i1 = (DS.qleft ? 1 : 2);
        ret.d2 = (DS.qleft ? 0 : 1e6);
        return ret;
    } else {
        ret.d1 = DS.x;
        return ret;
    }
}


struct TupleDD cumbet(double x, double y, double a, double b)
{
    //              Double precision cUMulative incomplete BETa distribution
    //
    //
    //                                  Function
    //
    //
    //         Calculates the cdf to X of the incomplete beta distribution
    //         with parameters a and b.  This is the integral from 0 to x
    //         of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
    //
    //
    //                                  Arguments
    //
    //
    //         X --> Upper limit of integration.
    //                                            X is DOUBLE PRECISION
    //
    //         Y --> 1 - X.
    //                                            Y is DOUBLE PRECISION
    //
    //         A --> First parameter of the beta distribution.
    //                                            A is DOUBLE PRECISION
    //
    //         B --> Second parameter of the beta distribution.
    //                                            B is DOUBLE PRECISION
    //
    //         CUM <-- Cumulative incomplete beta distribution.
    //                                            CUM is DOUBLE PRECISION
    //
    //         CCUM <-- Compliment of Cumulative incomplete beta distribution.
    //                                            CCUM is DOUBLE PRECISION
    //
    //
    //                                  Method
    //
    //
    //         Calls the routine BRATIO.
    //
    //                                       References
    //
    //         Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim
    //         708 Significant Digit Computation of the Incomplete Beta Function
    //         Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.

    if (x <= 0.0) {
        return (struct TupleDD){.d1 = 0.0, .d2 = 1.0};
    }
    if (y <= 0.0) {
        return (struct TupleDD){.d1 = 1.0, .d2 = 0.0};
    }
    struct TupleDDI res = bratio(a, b, x, y);
    return (struct TupleDD){.d1 = res.d1, .d2 = res.d2};
}


struct TupleDD cumbin(double s, double xn, double pr, double ompr)
{
    //                     CUmulative BINomial distribution
    //
    //
    //                               Function
    //
    //
    //      Returns the probability   of 0  to  S  successes in  XN   binomial
    //      trials, each of which has a probability of success, PBIN.
    //
    //
    //                               Arguments
    //
    //
    //      S --> The upper limit of cumulation of the binomial distribution.
    //                                                  S is DOUBLE PRECISION
    //
    //      XN --> The number of binomial trials.
    //                                                 XN is DOUBLE PRECISION
    //
    //      PBIN --> The probability of success in each binomial trial.
    //                                               PBIN is DOUBLE PRECISION
    //
    //      OMPR --> 1 - PBIN
    //                                               OMPR is DOUBLE PRECISION
    //
    //      CUM <-- Cumulative binomial distribution.
    //                                                CUM is DOUBLE PRECISION
    //
    //      CCUM <-- Compliment of Cumulative binomial distribution.
    //                                               CCUM is DOUBLE PRECISION
    //
    //
    //
    //                               Method
    //
    //
    //       Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
    //       Mathematical   Functions (1966) is   used  to reduce the  binomial
    //       distribution  to  the  cumulative    beta distribution.

    double cum, ccum;
    if (!(s < xn)) {
        return (struct TupleDD){.d1 = 1.0, .d2 = 0.0};
        cum = 1.0;
        ccum = 0.0;
    } else {
        struct TupleDD res = cumbet(pr, ompr, s + 1.0, xn - s);
        // Swap result
        return (struct TupleDD){.d1 = res.d2, .d2 = res.d1};
    }
}


struct TupleDD cumchi(double x, double df)
{
    //                 CUMulative of the CHi-square distribution
    //
    //
    //                                  Function
    //
    //
    //         Calculates the cumulative chi-square distribution.
    //
    //
    //                                  Arguments
    //
    //
    //         X       --> Upper limit of integration of the
    //                     chi-square distribution.
    //                                                    X is DOUBLE PRECISION
    //
    //         DF      --> Degrees of freedom of the
    //                     chi-square distribution.
    //                                                   DF is DOUBLE PRECISION
    //
    //         CUM <-- Cumulative chi-square distribution.
    //                                                  CUM is DOUBLE PRECISION
    //
    //         CCUM <-- Compliment of Cumulative chi-square distribution.
    //                                                 CCUM is DOUBLE PRECISION
    //
    //                                  Method
    //
    //         Calls incomplete gamma function (CUMGAM)
    //

    return cumgam(0.5*x, 0.5*df);
}


struct TupleDD cumchn(double x, double df, double pnonc)
{
    //             CUMulative of the Non-central CHi-square distribution
    //
    //                               Function
    //
    //     Calculates     the       cumulative      non-central    chi-square
    //     distribution, i.e.,  the probability   that  a   random   variable
    //     which    follows  the  non-central chi-square  distribution,  with
    //     non-centrality  parameter    PNONC  and   continuous  degrees   of
    //     freedom DF, is less than or equal to X.
    //
    //                              Arguments
    //
    //     X       --> Upper limit of integration of the non-central
    //                 chi-square distribution.
    //                                                 X is DOUBLE PRECISION
    //
    //     DF      --> Degrees of freedom of the non-central
    //                 chi-square distribution.
    //                                                 DF is DOUBLE PRECISION
    //
    //     PNONC   --> Non-centrality parameter of the non-central
    //                 chi-square distribution.
    //                                                 PNONC is DOUBLE PRECISION
    //
    //     CUM <-- Cumulative non-central chi-square distribution.
    //                                                 CUM is DOUBLE PRECISION
    //
    //     CCUM <-- Compliment of Cumulative non-central chi-square distribut
    //                                                 CCUM is DOUBLE PRECISION
    //
    //
    //                                Method
    //
    //     Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
    //     Mathematical    Functions,  US   NBS   (1966)    to calculate  the
    //     non-central chi-square.
    //
    //                                Variables
    //
    //     EPS     --- Convergence criterion.  The sum stops when a
    //                 term is less than EPS*SUM.
    //                                                 EPS is DOUBLE PRECISION
    //
    //     CCUM <-- Compliment of Cumulative non-central
    //              chi-square distribution.
    //                                                 CCUM is DOUBLE PRECISION

    double adj, centaj, centwt, chid2, dfd2, lcntaj, lcntwt;
    double lfact, pcent, pterm, ssum, sumadj, term, wt, xnonc;
    double eps = 1.e-15;
    double abstol = 1.e-300;
    int i, icent;

    if (!(x > 0.)) {
        return (struct TupleDD){.d1 = 0.0, .d2 = 1.0};
    }
    if (!(pnonc > 1e-10)) {
        return cumchi(x, df);
    }
    xnonc = pnonc/2.0;
    icent = (int)xnonc;
    if (icent == 0) {
        icent = 1;
    }
    chid2 = x / 2.0;

    lfact = alngam(icent + 1);
    lcntwt = -xnonc + icent*log(xnonc) - lfact;
    centwt = exp(lcntwt);

    struct TupleDD res = cumchi(x, df + 2.*icent);
    pcent = res.d1;
    dfd2 = (df + 2.*icent)/2.;
    lfact = alngam(1. + dfd2);
    lcntaj = dfd2*log(chid2) - chid2 - lfact;
    centaj = exp(lcntaj);
    ssum = centwt*pcent;

    sumadj = 0.;
    adj = centaj;
    wt = centwt;
    i = icent;

    while (1) {
        dfd2 = (df + 2.*i)/2.;
        adj *= dfd2 / chid2;
        sumadj += adj;
        pterm = pcent + sumadj;
        wt *= (i / xnonc);
        term = wt*pterm;
        ssum += term;
        i -= 1;
        if ((!((ssum >= abstol) && (term >= eps*ssum))) || (i == 0)) { break; }
    }
    sumadj = centaj;
    adj = centaj;
    wt = centwt;
    i = icent;

    while (1) {
        wt *= (xnonc / (i+1.));
        pterm = pcent - sumadj;
        term = wt*pterm;
        ssum += term;
        i += 1;
        dfd2 = (df + 2.*i)/2.0;
        adj *= chid2/dfd2;
        sumadj += adj;
        if (!((ssum >= abstol) && (term >= eps*ssum))) { break; }
    }

    return (struct TupleDD){.d1  = ssum, .d2 = 0.5 + (0.5 - ssum)};
}


struct TupleDD cumf(double f, double dfn, double dfd)
{
    //                    CUMulative F distribution
    //
    //
    //                              Function
    //
    //
    //     Computes  the  integral from  0  to  F of  the f-density  with DFN
    //     and DFD degrees of freedom.
    //
    //
    //                              Arguments
    //
    //
    //     F --> Upper limit of integration of the f-density.
    //                                                  F is DOUBLE PRECISION
    //
    //     DFN --> Degrees of freedom of the numerator sum of squares.
    //                                                  DFN is DOUBLE PRECISION
    //
    //     DFD --> Degrees of freedom of the denominator sum of squares.
    //                                                  DFD is DOUBLE PRECISION
    //
    //     CUM <-- Cumulative f distribution.
    //                                                  CUM is DOUBLE PRECISION
    //
    //     CCUM <-- Compliment of Cumulative f distribution.
    //                                                  CCUM is DOUBLE PRECISION
    //
    //
    //                              Method
    //
    //
    //     Formula  26.5.28 of  Abramowitz and   Stegun   is  used to  reduce
    //     the cumulative F to a cumulative beta distribution.
    //
    //
    //                              Note
    //
    //
    //     If F is less than or equal to 0, 0 is returned.

    double dsum, prod, xx, yy;

    if (f <= 0.0) {
        return (struct TupleDD){.d1 = 0.0, .d2 = 1.0};
    }
    prod = dfn*f;
    dsum = dfd + prod;
    xx = dfd / dsum;
    if (xx > 0.5) {
        yy = prod / dsum;
        xx = 1. - yy;
    } else {
        yy = 1. - xx;
    }
    struct TupleDDI res = bratio(dfd*0.5, dfn*0.5, xx, yy);
    return (struct TupleDD){.d1 = res.d2, .d2 = res.d1};
}


struct TupleDDI cumfnc(double f, double dfn, double dfd, double pnonc)
{
    //            F -NON- -C-ENTRAL F DISTRIBUTION
    //
    //
    //
    //                            Function
    //
    //
    //    COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD
    //    DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC
    //
    //
    //                            Arguments
    //
    //
    //    X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION
    //
    //    DFN --> DEGREES OF FREEDOM OF NUMERATOR
    //
    //    DFD -->  DEGREES OF FREEDOM OF DENOMINATOR
    //
    //    PNONC --> NONCENTRALITY PARAMETER.
    //
    //    CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION
    //
    //    CCUM <-- COMPLIMENT OF CUMMULATIVE
    //
    //
    //                            Method
    //
    //
    //    USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES.
    //    SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2
    //    (THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL
    //    THE CONVERGENCE CRITERION IS MET.
    //
    //    FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED
    //    BY FORMULA 26.5.16.
    //
    //
    //            REFERENCE
    //
    //
    //    HANDBOOD OF MATHEMATICAL FUNCTIONS
    //    EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN
    //    NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
    //    MARCH 1965
    //    P 947, EQUATIONS 26.6.17, 26.6.18
    //
    //
    //                            Note
    //
    //
    //    THE SUM CONTINUES UNTIL A SUCCEEDING TERM IS LESS THAN EPS
    //    TIMES THE SUM (OR THE SUM IS LESS THAN 1.0E-20).  EPS IS
    //    SET TO 1.0E-4 IN A DATA STATEMENT WHICH CAN BE CHANGED.

    double dsum, prod, xx, yy, adn, aup, b;
    double betdn, betup, centwt, dnterm, ssum;
    double upterm, xmult, xnonc;
    double eps = 1e-4;
    double abstol = 1e-300;
    int status = 0;
    int i, icent;

    if (!(f > 0.0)) {
        return (struct TupleDDI){.d1 = 0.0, .d2 = 1.0, .i1 = status};
    }

    if (!(pnonc > 1e-10)) {
        struct TupleDD res = cumf(f, dfn, dfd);
        return (struct TupleDDI){.d1 = res.d1, .d2 = res.d2, .i1=status};
    }

    xnonc = pnonc / 2.0;
    icent = (int)xnonc;
    if (fabs(xnonc - icent) >= 1.0) {
        return (struct TupleDDI){.d1 = 0.0, .d2 = 0.0, .i1 = 1};
    }
    if (icent == 0) {
        icent = 1;
    }
    centwt = exp(-xnonc + icent*log(xnonc) - alngam(icent + 1));
    prod = dfn * f;
    dsum = dfd + prod;
    yy = dfd / dsum;
    if (yy > 0.5) {
        xx = prod / dsum;
        yy = 1. - xx;
    } else {
        xx = 1. - yy;
    }
    struct TupleDDI resi = bratio(dfn*0.5 + icent, dfd*0.5, xx, yy);
    betdn = resi.d1;
    adn = dfn/2. + icent;
    aup = adn;
    b = dfd / 2.0;
    betup = betdn;
    ssum = centwt*betdn;

    xmult = centwt;
    i = icent;

    if (adn < 2.0) {
        dnterm = exp(alngam(adn+b) - alngam(adn+1.) - alngam(b) +
                     adn*log(xx) + b*log(yy));
    } else {
        dnterm = exp(-betaln(adn, b) - log(adn) + adn*log(xx) +
                     b*log(yy));
    }

    while (((ssum >= abstol) && (xmult*betdn >= eps*ssum)) && (i > 0)) {
        xmult *= (i/xnonc);
        i -= 1;
        adn -= 1;
        dnterm *= (adn + 1) / ((adn + b)*xx);
        betdn += dnterm;
        ssum += xmult*betdn;
    }

    i = icent + 1;
    xmult = centwt;

    if ((aup - 1 + b) == 0) {
          upterm = exp(-alngam(aup) - alngam(b) + (aup-1)*log(xx) + b*log(yy));
    } else {
        if (aup < 2) {
            upterm = exp(alngam(aup-1+b) - alngam(aup)
                         - alngam(b) + (aup-1)*log(xx) + b*log(yy));
        } else {
            // Same expression, but avoid problems for large aup
            upterm = exp(-betaln(aup-1, b) - log(aup-1) +
                         (aup-1)*log(xx) + b*log(yy));
        }
    }

    while (1) {
        xmult *= xnonc / i;
        i += 1;
        aup += 1;
        upterm *= (aup + b - 2.)*xx/(aup - 1.);
        betup -= upterm;
        ssum += xmult*betup;
        if (!((ssum >= abstol) && (xmult*betup >= eps*ssum))) { break; }
    }

    return (struct TupleDDI){.d1 = ssum, .d2 = 0.5 + (0.5 - ssum), .i1 = status};
}


struct TupleDD cumgam(double x, double a)
{
    //        Double precision cUMulative incomplete GAMma distribution
    //
    //
    //                            Function
    //
    //
    //    Computes   the  cumulative        of    the     incomplete   gamma
    //    distribution, i.e., the integral from 0 to X of
    //        (1/GAM(A))*EXP(-T)*T**(A-1) DT
    //    where GAM(A) is the complete gamma function of A, i.e.,
    //        GAM(A) = integral from 0 to infinity of
    //                EXP(-T)*T**(A-1) DT
    //
    //
    //                            Arguments
    //
    //
    //    X --> The upper limit of integration of the incomplete gamma.
    //                                            X is DOUBLE PRECISION
    //
    //    A --> The shape parameter of the incomplete gamma.
    //                                            A is DOUBLE PRECISION
    //
    //    CUM <-- Cumulative incomplete gamma distribution.
    //                                    CUM is DOUBLE PRECISION
    //
    //    CCUM <-- Compliment of Cumulative incomplete gamma distribution.
    //                                            CCUM is DOUBLE PRECISIO
    //
    //
    //                            Method
    //
    //
    //    Calls the routine GRATIO.
    if (x > 0.0) {
        return gratio(a, x, 0);
    } else {
        return (struct TupleDD){.d1 = 0.0, .d2 = 1.0};
    }
}


struct TupleDD cumnbn(double s, double xn, double pr, double ompr)
{
    //                CUmulative Negative BINomial distribution
    //
    //
    //                            Function
    //
    //
    //    Returns the probability that it there will be S or fewer failures
    //    before there are XN successes, with each binomial trial having
    //    a probability of success PR.
    //
    //    Prob(# failures = S | XN successes, PR)  =
    //                    ( XN + S - 1 )
    //                    (            ) * PR^XN * (1-PR)^S
    //                    (      S     )
    //
    //
    //                            Arguments
    //
    //
    //    S --> The number of failures
    //                                                S is DOUBLE PRECISION
    //
    //    XN --> The number of successes
    //                                                XN is DOUBLE PRECISION
    //
    //    PR --> The probability of success in each binomial trial.
    //                                                PR is DOUBLE PRECISION
    //
    //    OMPR --> 1 - PR
    //                                                OMPR is DOUBLE PRECISION
    //
    //    CUM <-- Cumulative negative binomial distribution.
    //                                                CUM is DOUBLE PRECISION
    //
    //    CCUM <-- Compliment of Cumulative negative binomial distribution.
    //                                                CCUM is DOUBLE PRECISION
    //
    //
    //                            Method
    //
    //    Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of
    //    Mathematical   Functions (1966) is   used  to reduce the  negative
    //    binomial distribution to the cumulative beta distribution.
    return cumbet(pr, ompr, xn, s+1.);
}


struct TupleDD cumnor(double x)
{
    //                                Function
    //        Computes the cumulative  of    the  normal   distribution,   i.e.,
    //        the integral from -infinity to x of
    //            (1/sqrt(2*pi)) exp(-u*u/2) du
    //        X --> Upper limit of integration.
    //                                            X is DOUBLE PRECISION
    //        RESULT <-- Cumulative normal distribution.
    //                                            RESULT is DOUBLE PRECISION
    //        CCUM <-- Compliment of Cumulative normal distribution.
    //                                            CCUM is DOUBLE PRECISION
    //        Renaming of function ANORM from:
    //        Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
    //        Package of Special Function Routines and Test Drivers"
    //        acm Transactions on Mathematical Software. 19, 22-32.
    //        with slight modifications to return ccum and to deal with
    //        machine constants.
    //    **********************************************************************
    //    Original Comments:
    //    ------------------------------------------------------------------
    //    This function evaluates the normal distribution function:
    //                                / x
    //                        1       |       -t*t/2
    //            P(x) = -----------  |      e       dt
    //                    sqrt(2 pi)  |
    //                                /-oo
    //    The main computation evaluates near-minimax approximations
    //    derived from those in "Rational Chebyshev approximations for
    //    the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
    //    This transportable program uses rational functions that
    //    theoretically approximate the normal distribution function to
    //    at least 18 significant decimal digits.  The accuracy achieved
    //    depends on the arithmetic system, the compiler, the intrinsic
    //    functions, and proper selection of the machine-dependent
    //    constants.

    double a[5] = {2.2352520354606839287e00, 1.6102823106855587881e02,
                   1.0676894854603709582e03, 1.8154981253343561249e04,
                   6.5682337918207449113e-2};
    double b[4] = {4.7202581904688241870e01, 9.7609855173777669322e02,
                   1.0260932208618978205e04, 4.5507789335026729956e04};
    double c[9] = {3.9894151208813466764e-1, 8.8831497943883759412e00,
                   9.3506656132177855979e01, 5.9727027639480026226e02,
                   2.4945375852903726711e03, 6.8481904505362823326e03,
                   1.1602651437647350124e04, 9.8427148383839780218e03,
                   1.0765576773720192317e-8};
    double d[8] = {2.2266688044328115691e01, 2.3538790178262499861e02,
                   1.5193775994075548050e03, 6.4855582982667607550e03,
                   1.8615571640885098091e04, 3.4900952721145977266e04,
                   3.8912003286093271411e04, 1.9685429676859990727e04};
    double p[6] = {2.1589853405795699e-01, 1.274011611602473639e-01,
                   2.2235277870649807e-02, 1.421619193227893466e-03,
                   2.9112874951168792e-05, 2.307344176494017303e-02};
    double q[5] = {1.28426009614491121e00, 4.68238212480865118e-01,
                   6.59881378689285515e-02, 3.78239633202758244e-03,
                   7.29751555083966205e-05};
    double eps = spmpar[0] * 0.5;
    double tiny = spmpar[1];
    double y = fabs(x);
    double threshold = 0.66291;
    double dl, result, xden, xnum, xsq, ccum;
    int i;

    if (y <= threshold) {
        // Evaluate  anorm  for  |X| <= 0.66291
        xsq = y > eps ? x*x : 0;
        xnum = a[4] * xsq;
        xden = xsq;
        for (int i=0; i<3; i++){
            xnum += a[i];
            xnum *= xsq;
            xden += b[i];
            xden *= xsq;
        }
        result = x*(xnum + a[3]) / (xden + b[3]);
        ccum = 0.5 - result;
        result += 0.5;
    } else if (y < sqrt(32)) {
        // Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
        xnum = c[8]*y;
        xden = y;
        for (i = 0; i < 7; i++){
            xnum += c[i];
            xnum *= y;
            xden += d[i];
            xden *= y;
        }
        result = (xnum + c[7]) / (xden + d[7]);
        xsq = (int)(y*1.6) / 1.6;
        dl = (y - xsq) * (y + xsq);
        result *= exp(-xsq*xsq*0.5)*exp(-0.5*dl);
        ccum = 1. - result;
    } else {
        // Evaluate  anorm  for |X| > sqrt(32)
        result = 0.;
        xsq = (1.0 / x / x);
        xnum = p[5] * xsq;
        xden = xsq;
        for(i = 0; i < 4; i++) {
            xnum += p[i];
            xnum *= xsq;
            xden += q[i];
            xden *= xsq;
        }
        result = xsq*(xnum + p[4]) / (xden + q[4]);
        result = (sqrt(1./(2.0*PI)) - result) / y;
        xsq = (int)(x*1.6) / 1.6;
        dl = (x - xsq) * (x + xsq);
        result *= exp(-xsq*xsq*0.5)*exp(-0.5*dl);
        ccum = 1. - result;
    }

    if (x > 0) {
        double tmp = result;
        result = ccum;
        ccum = tmp;
    }

    if (result < tiny){
        result = 0.0;
    }
    if (ccum < tiny){
        ccum = 0.0;
    }

    return (struct TupleDD){.d1 = result, .d2 = ccum};
}


struct TupleDD cumpoi(double s, double xlam)
{
    //                CUMulative POIsson distribution
    //
    //
    //                            Function
    //
    //
    //    Returns the  probability  of  S   or  fewer events in  a   Poisson
    //    distribution with mean XLAM.
    //
    //
    //                            Arguments
    //
    //
    //    S --> Upper limit of cumulation of the Poisson.
    //                                                S is DOUBLE PRECISION
    //
    //    XLAM --> Mean of the Poisson distribution.
    //                                                XLAM is DOUBLE PRECIS
    //
    //    CUM <-- Cumulative poisson distribution.
    //                                    CUM is DOUBLE PRECISION
    //
    //    CCUM <-- Compliment of Cumulative poisson distribution.
    //                                                CCUM is DOUBLE PRECIS
    //
    //
    //                            Method
    //
    //
    //    Uses formula  26.4.21   of   Abramowitz and  Stegun,  Handbook  of
    //    Mathematical   Functions  to reduce   the   cumulative Poisson  to
    //    the cumulative chi-square distribution.

    struct TupleDD res = cumchi(2*xlam, 2.*(s + 1.));
    return (struct TupleDD){.d1 = res.d2, .d2 = res.d1};
}


struct TupleDD cumt(double t, double df)
{
    //                CUMulative T-distribution
    //
    //
    //                            Function
    //
    //
    //    Computes the integral from -infinity to T of the t-density.
    //
    //
    //                            Arguments
    //
    //
    //    T --> Upper limit of integration of the t-density.
    //                                                T is DOUBLE PRECISION
    //
    //    DF --> Degrees of freedom of the t-distribution.
    //                                                DF is DOUBLE PRECISIO
    //
    //    CUM <-- Cumulative t-distribution.
    //                                                CCUM is DOUBLE PRECIS
    //
    //    CCUM <-- Compliment of Cumulative t-distribution.
    //                                                CCUM is DOUBLE PRECIS
    //
    //
    //                            Method
    //
    //
    //    Formula 26.5.27   of     Abramowitz  and   Stegun,    Handbook  of
    //    Mathematical Functions  is   used   to  reduce the  t-distribution
    //    to an incomplete beta.

    double a, oma, tt, dfptt, xx, yy, cum, ccum;

    tt = t*t;
    dfptt = df + tt;
    xx = df / dfptt;
    yy = tt / dfptt;
    struct TupleDD res = cumbet(xx, yy, 0.5*df, 0.5);
    a = res.d1;
    oma = res.d2;
    if (t > 0.0) {
        ccum = 0.5 * a;
        cum = oma + ccum;
    } else {
        cum = 0.5 * a;
        ccum = oma + cum;
    }
    return (struct TupleDD){.d1 = cum, .d2 = ccum};
}


struct TupleDD cumtnc(double t, double df, double pnonc)
{
    //                CUMulative Non-Central T-distribution
    //
    //
    //                            Function
    //
    //
    //    Computes the integral from -infinity to T of the non-central
    //    t-density.
    //
    //
    //                            Arguments
    //
    //
    //    T --> Upper limit of integration of the non-central t-density.
    //                                                T is DOUBLE PRECISION
    //
    //    DF --> Degrees of freedom of the non-central t-distribution.
    //                                                DF is DOUBLE PRECISION
    //
    //    PNONC --> Non-centrality parameter of the non-central t distibutio
    //                                                PNONC is DOUBLE PRECISION
    //
    //    CUM <-- Cumulative t-distribution.
    //                                                CCUM is DOUBLE PRECISION
    //
    //    CCUM <-- Compliment of Cumulative t-distribution.
    //                                                CCUM is DOUBLE PRECISION
    //
    //
    //                            Method
    //
    //    Upper tail    of  the  cumulative  noncentral t   using
    //    formulae from page 532  of Johnson, Kotz,  Balakrishnan, Coninuous
    //    Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
    //
    //    This implementation starts the calculation at i = lambda,
    //    which is near the largest Di.  It then sums forward and backward.
    double alghdf, b, bb, bbcent, bcent, cent, cum, ccum, d, dcent;
    double dpnonc, dum1, dum2, e, ecent, lmbda, lnomx, lnx, omx;
    double pnonc2, s, scent, ss, sscent, t2, term, tt, twoi, x;
    double xi, xlnd, xlne;
    double conv = 1.e-7;
    double tiny = 1.e-10;
    int qrevs;

    if (fabs(pnonc) <= tiny) { return cumt(t, df); }
    qrevs = t < 0.0;
    tt = (qrevs ? -t : t);
    dpnonc = (qrevs ? -pnonc : pnonc);
    pnonc2 = pow(dpnonc, 2);
    t2 = pow(tt, 2);
    if (fabs(tt) <= tiny) { return cumnor(-pnonc); }

    lmbda = 0.5*pnonc2;
    x = df / (df + t2);
    omx = 1. - x;
    lnx = log(x);
    lnomx = log(omx);
    alghdf = gamln(0.5*df);
    // Case : i = lambda
    cent = fmax(floor(lmbda), 1.);
    // Compute d=T(2i) in log space and offset by exp(-lambda)
    xlnd = cent*log(lmbda) - gamln(cent + 1.) - lmbda;
    dcent = exp(xlnd);
    // Compute e=t(2i+1) in log space offset by exp(-lambda)
    xlne = (cent + 0.5)*log(lmbda) - gamln(cent + 1.5) - lmbda;
    ecent = exp(xlne);
    if (dpnonc < 0.0) {
        ecent = -ecent;
    }
    // Compute bcent=B(2*cent)
    struct TupleDDI res1 = bratio(0.5*df, cent + 0.5, x, omx);
    bcent = res1.d1;
    dum1 = res1.d2;
    // Compute bbcent=B(2*cent+1)
    struct TupleDDI res2 = bratio(0.5*df, cent + 1., x, omx);
    bbcent = res2.d1;
    dum2 = res2.d2;

    // Case bcent and bbcent are essentially zero
    // Thus t is effectively infinite
    if ((bbcent + bcent) < tiny) {
        if (qrevs) {
            return (struct TupleDD){.d1 = 0.0, .d2 = 1.0};
        } else {
            return (struct TupleDD){.d1 = 1.0, .d2 = 0.0};
        }
    }
    // Case bcent and bbcent are essentially one
    // Thus t is effectively zero
    if ((dum1 + dum2) < tiny) { return cumnor(-pnonc); }

    // First term in ccum is D*B + E*BB
    ccum = dcent*bcent + ecent*bbcent;
    // Compute s(cent) = B(2*(cent+1)) - B(2*cent))
    scent = exp(gamln(0.5*df + cent + 0.5) - gamln(cent + 1.5) - alghdf
                + 0.5*df*lnx + (cent + 0.5)*lnomx);
    // Compute ss(cent) = B(2*cent+3) - B(2*cent+1)
    sscent = exp(gamln(0.5*df + cent + 1.) - gamln(cent + 2.) - alghdf
                 + 0.5*df*lnx + (cent + 1.)*lnomx);
    // Sum forward
    xi = cent + 1.0;
    twoi = 2.*xi;
    d = dcent;
    e = ecent;
    b = bcent;
    s = scent;
    bb = bbcent;
    ss = sscent;

    while (1) {
        b += s;
        bb += ss;
        d = (lmbda/xi)*d;
        e = (lmbda/(xi + 0.5))*e;
        term = d*b + e*bb;
        ccum += term;
        s *= omx*(df + twoi - 1.)/(twoi + 1.);
        ss *= omx*(df + twoi)/(twoi + 2.);
        xi += 1.0;
        twoi = 2.0*xi;
        if (fabs(term) <= conv*ccum) { break; }
    }
    // Sum Backward
    xi = cent;
    twoi = 2.0*xi;
    d = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s = scent*(1.0 + twoi)/((df + twoi - 1.)*omx);
    ss = sscent*(2.0 + twoi)/((df + twoi)*omx);

    while (1) {
        b -= s;
        bb -= ss;
        d *= (xi / lmbda);
        e *= (xi + 0.5)/lmbda;
        term = d*b + e*bb;
        ccum += term;
        xi -= 1.0;
        if (xi < 0.5) { break; }
        twoi = 2.0*xi;
        s *= (1. + twoi) / ((df + twoi - 1.)*omx);
        ss *= (2. + twoi) / ((df + twoi)*omx);
        if (fabs(term) <= conv*ccum) { break; }
    }
    // Due to roundoff error the answer may not lie between zero and one
    // Force it to do so
    if (qrevs) {
        cum = fmax(fmin(0.5*ccum, 1.), 0.);
        ccum = fmax(fmin(1.-cum, 1.), 0.);
    } else {
        ccum = fmax(fmin(0.5*ccum, 1.), 0.);
        cum = fmax(fmin(1.-ccum, 1.), 0.);
    }
    return (struct TupleDD){.d1 = cum, .d2 = ccum};
}


double devlpl(double *a, int n, double x)
{
    //            Double precision EVALuate a PoLynomial at X
    //
    //
    //                            Function
    //
    //
    //    returns
    //        A(1) + A(2)*X + ... + A(N)*X**(N-1)
    //
    //
    //                            Arguments
    //
    //
    //    A --> Array of coefficients of the polynomial.
    //                                    A is DOUBLE PRECISION(N)
    //
    //    N --> Length of A, also degree of polynomial - 1.
    //                                    N is INTEGER
    //
    //    X --> Point at which the polynomial is to be evaluated.
    //                                    X is DOUBLE PRECISION
    double temp = a[n-1];
    int i;

    for (i = n - 2; i >= 0; i--) {
        temp = a[i] + temp*x;
    }
    return temp;
}


double dinvnr(double p, double q)
{
    //    Double precision NoRmal distribution INVerse
    //
    //
    //                            Function
    //
    //
    //    Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
    //    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
    //
    //
    //                            Arguments
    //
    //
    //    P --> The probability whose normal deviate is sought.
    //                P is DOUBLE PRECISION
    //
    //    Q --> 1-P
    //                P is DOUBLE PRECISION
    //
    //
    //                            Method
    //
    //
    //    The  rational   function   on  page 95    of Kennedy  and  Gentle,
    //    Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
    //    value for the Newton method of finding roots.
    //
    //
    //                            Note
    //
    //
    //    If P or Q < machine EPS returns +/- DINVNR(EPS)

    int i, maxit = 100;
    double eps = 1e-13;
    double r2pi = sqrt(1. / (2.*PI));
    double strtx, xcur, cum, pp, dx;

    pp = (p > q ? q : p);
    strtx = stvaln(pp);
    xcur = strtx;

    for (i = 0; i < maxit; i++) {
        struct TupleDD res = cumnor(xcur);
        cum = res.d1;
        dx = (cum - pp) / (r2pi * exp(-0.5*xcur*xcur));
        xcur -= dx;
        if (fabs(dx / xcur) < eps) {
            return (p > q ? -xcur : xcur);
        }
    }
    return (p > q ? -strtx : strtx);
}


void dinvr(DinvrState *S, DzrorState *DZ)
{
    //        Double precision
    //        bounds the zero of the function and invokes zror
    //                Reverse Communication
    //
    //
    //                            Function
    //
    //
    //    Bounds the    function  and  invokes  ZROR   to perform the   zero
    //    finding.  STINVR  must  have   been  called  before this   routine
    //    in order to set its parameters.
    //
    //
    //                            Arguments
    //
    //
    //    STATUS <--> At the beginning of a zero finding problem, STATUS
    //                should be set to 0 and INVR invoked.  (The value
    //                of parameters other than X will be ignored on this cal
    //
    //                When INVR needs the function evaluated, it will set
    //                STATUS to 1 and return.  The value of the function
    //                should be set in FX and INVR again called without
    //                changing any of its other parameters.
    //
    //                When INVR has finished without error, it will return
    //                with STATUS 0.  In that case X is approximately a root
    //                of F(X).
    //
    //                If INVR cannot bound the function, it returns status
    //                -1 and sets QLEFT and QHI.
    //                        INTEGER STATUS
    //
    //    X <-- The value of X at which F(X) is to be evaluated.
    //                        DOUBLE PRECISION X
    //
    //    FX --> The value of F(X) calculated when INVR returns with
    //        STATUS = 1.
    //                        DOUBLE PRECISION FX
    //
    //    QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
    //        case it is .TRUE. If the stepping search terminated
    //        unsuccessfully at SMALL.  If it is .FALSE. the search
    //        terminated unsuccessfully at BIG.
    //                QLEFT is LOGICAL
    //
    //    QHI <-- Defined only if QMFINV returns .FALSE.  In that
    //        case it is .TRUE. if F(X) > Y at the termination
    //        of the search and .FALSE. if F(X) < Y at the
    //        termination of the search.
    //                QHI is LOGICAL
    while (1) {
        if (S->next_state == 0) {
            // See that SMALL and BIG bound the zero and set QINCR
            S->qcond = ((S->small <= S->x) && (S->x <= S->big));
            if (!(S->qcond)) {
                S->status = -2;
                return;
            }
            S->xsave = S->x;
            S->x = S->small;
            S->next_state = 10;
            S->status = 1;
            return;
        } else if (S->next_state == 10) {
            S->fsmall = S->fx;
            S->x = S->big;
            S->next_state = 20;
            S->status = 1;
            return;
        } else if (S->next_state == 20) {
            S->fbig = S->fx;
            S->qincr = (S->fbig > S->fsmall);
            S->status = -1;
            if (!(S->qincr)) {
                // 50
                if (S->fsmall >= 0.0) {
                    // 60
                    if (!(S->fbig <= 0.0)) {
                        S->qleft = 0;
                        S->qhi = 1;
                        return;
                    }
                } else {
                    S->qleft = 1;
                    S->qhi = 0;
                    return;
                }
            } else {
                if (S->fsmall <= 0.0) {
                    // 30
                    if (!(S->fbig >= 0.0)) {
                        S->qleft = 0;
                        S->qhi = 0;
                        return;
                    }
                } else {
                    S->qleft = 1;
                    S->qhi = 1;
                    return;
                }
            }
            S->status = 1;
            S->x = S->xsave;
            S->step = fmax(S->absstp, S->relstp*fabs(S->x));
            S->next_state = 90;
            return;
        } else if (S->next_state == 90) {
            S->yy = S->fx;
            if (S->yy == 0.0) {
                S->status = 0;
                S->qok = 1;
                return;
            }
            S->next_state = 100;
    //  Handle case in which we must step higher
        } else if (S->next_state == 100) {
            S->qup = (S->qincr & (S->yy < 0.)) || ((!(S->qincr)) && (S->yy > 0.));
            if (S->qup) {
                S->xlb = S->xsave;
                S->xub = fmin(S->xlb + S->step, S->big);
                S->next_state = 120;
            } else {
                // 170
                S->xub = S->xsave;
                S->xlb = fmax(S->xub - S->step, S->small);
                S->next_state = 190;
            }
        } else if (S->next_state == 120) {
            S->x = S->xub;
            S->status = 1;
            S->next_state = 130;
            return;
        } else if (S->next_state == 130) {
            S->yy = S->fx;
            S->qbdd = (S->qincr & (S->yy >= 0.)) || ((!(S->qincr)) && (S->yy <= 0.));
            S->qlim = (S->xub >= S->big);
            S->qcond = ((S->qbdd) || (S->qlim));
            if (S->qcond) {
                S->next_state = 150;
            } else {
                S->step *= S->stpmul;
                S->xlb = S->xub;
                S->xub = fmin(S->xlb + S->step, S->big);
                S->next_state = 120;
            }
        } else if (S->next_state == 150) {
            if (S->qlim & (!(S->qbdd))) {
                S->status = -1;
                S->qleft = 0;
                S->qhi = (!(S->qincr));
                S->x = S->big;
                return;
            } else {
                S->next_state = 240;
            }
        } else if (S->next_state == 190) {
            S->x = S->xlb;
            S->status = 1;
            S->next_state = 200;
            return;
    // Handle case in which we must step lower
        } else if (S->next_state == 200) {
            S->yy = S->fx;
            S->qbdd = (((S->qincr) && (S->yy <= 0.0)) || ((!(S->qincr)) && (S->yy >= 0.0)));
            S->qlim = ((S->xlb) <= (S->small));
            S->qcond = ((S->qbdd) || (S->qlim));
            if (S->qcond) {
                S->next_state = 220;
            } else {
                S->step *= S->stpmul;
                S->xub = S->xlb;
                S->xlb = fmax(S->xub - S->step, S->small);
                S->next_state = 190;
            }
        } else if (S->next_state == 220) {
            if ((S->qlim) && (!(S->qbdd))) {
                S->status = -1;
                S->qleft = 1;
                S->qhi = S->qincr;
                S->x = S->small;
                return;
            } else {
                S->next_state = 240;
            }
    // If we reach here, xlb and xub bound the zero of f.
        } else if (S->next_state == 240) {
            // Overwrite supplied DZ with the problem
            DZ->xhi = S->xub;
            DZ->xlo = S->xlb;
            DZ->atol = S->abstol;
            DZ->rtol = S->reltol;
            DZ->x = S->xlb;
            DZ->b = S->xlb;
            S->next_state = 250;
        } else if (S->next_state == 250) {
            dzror(DZ);
            if (DZ->status == 1) {
                S->next_state = 260;
                S->status = 1;
                S->x = DZ->x;
                return;
            } else {
                S->x = DZ->xlo;
                S->status = 0;
                return;
            }
        } else if (S->next_state == 260) {
            DZ->fx = S->fx;
            S->next_state = 250;
        } else {
    // Bad state, should not be possible to get here
            S->status = -9999;
            return;
        }
    }
}


double dt1(double p, double q, double df)
{
    //    Double precision Initialize Approximation to
    //        INVerse of the cumulative T distribution
    //
    //
    //                            Function
    //
    //
    //    Returns  the  inverse   of  the T   distribution   function, i.e.,
    //    the integral from 0 to INVT of the T density is P. This is an
    //    initial approximation
    //
    //
    //                            Arguments
    //
    //
    //    P --> The p-value whose inverse from the T distribution is
    //        desired.
    //                P is DOUBLE PRECISION
    //
    //    Q --> 1-P.
    //                Q is DOUBLE PRECISION
    //
    //    DF --> Degrees of freedom of the T distribution.
    //                DF is DOUBLE PRECISION

    double ssum, term, x, xx;
    double denpow = 1.0;
    double coef[4][5] = {{1., 1., 0., 0., 0.},
                              {3., 16., 5., 0., 0.},
                              {-15., 17., 19., 3., 0.},
                              {-945., -1920., 1482., 776., 79.}
                             };
    double denom[4] = {4., 96., 384.0, 92160.0};
    int i, ideg[4] = {2, 3, 4, 5};

    x = fabs(dinvnr(p, q));
    xx = x*x;
    ssum = x;
    for (i = 0; i < 4; i++){
        term = (devlpl(coef[i], ideg[i], xx))*x ;
        denpow *= df;
        ssum += term / (denpow*denom[i]);
    }

    return (p >= 0.5 ? ssum : -ssum);
}


void dzror(DzrorState *S)
{
    //    Double precision ZeRo of a function -- Reverse Communication
    //
    //
    //                            Function
    //
    //
    //    Performs the zero finding.  STZROR must have been called before
    //    this routine in order to set its parameters.
    //
    //
    //                            Arguments
    //
    //
    //    STATUS <--> At the beginning of a zero finding problem, STATUS
    //                should be set to 0 and ZROR invoked.  (The value
    //                of other parameters will be ignored on this call.)
    //
    //                When ZROR needs the function evaluated, it will set
    //                STATUS to 1 and return.  The value of the function
    //                should be set in FX and ZROR again called without
    //                changing any of its other parameters.
    //
    //                When ZROR has finished without error, it will return
    //                with STATUS 0.  In that case (XLO,XHI) bound the answe
    //
    //                If ZROR finds an error (which implies that F(XLO)-Y an
    //                F(XHI)-Y have the same sign, it returns STATUS -1.  In
    //                this case, XLO and XHI are undefined.
    //                        INTEGER STATUS
    //
    //    X <-- The value of X at which F(X) is to be evaluated.
    //                        DOUBLE PRECISION X
    //
    //    FX --> The value of F(X) calculated when ZROR returns with
    //        STATUS = 1.
    //                        DOUBLE PRECISION FX
    //
    //    XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
    //            inverval in X containing the solution below.
    //                        DOUBLE PRECISION XLO
    //
    //    XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
    //            inverval in X containing the solution above.
    //                        DOUBLE PRECISION XHI
    //
    //    QLEFT <-- .TRUE. if the stepping search terminated unsuccessfully
    //            at XLO.  If it is .FALSE. the search terminated
    //            unsuccessfully at XHI.
    //                QLEFT is LOGICAL
    //
    //    QHI <-- .TRUE. if F(X) > Y at the termination of the
    //            search and .FALSE. if F(X) < Y at the
    //            termination of the search.
    //                QHI is LOGICAL

    while (1) {
        if (S->next_state == 0) {
            S->next_state = 10;
            S->status = 1;
            return;
        } else if (S->next_state == 10) {
            S->fb = S->fx;
            S->xlo = S->xhi;
            S->a = S->xlo;
            S->x = S->xlo;
            S->next_state = 20;
            S->status = 1;
            return;

        } else if (S->next_state == 20) {
            if (!(S->fb > 0.) || (!(S->fx > 0.))) {
                if ((!(S->fb < 0.)) || (!(S->fx < 0.))) {
        // 60
                    S->fa = S->fx;
                    S->first = 1;
        // 70
                    S->c = S->a;
                    S->fc = S->fa;
                    S->ext = 0;
                    S->status = 1;
                    S->next_state = 80;
                } else {
        // 40
                    S->status = -1;
                    S->qleft = (S->fx > S->fb);
                    S->qhi = 1;
                    return;
                }
            } else {
                S->status = -1;
                S->qleft = (S->fx < S->fb);
                S->qhi = 0;
                return;
            }
        } else if (S->next_state == 80) {
            if (fabs(S->fc) < fabs(S->fb)) {
                if (!(S->c == S->a)) {
                    S->d = S->a;
                    S->fd = S->fa;
                }
            // 90
                S->a = S->b;
                S->fa = S->fb;
                S->xlo = S->c;
                S->b = S->xlo;
                S->fb = S->fc;
                S->c = S->a;
                S->fc = S->fa;
            }
            // 100
            S->tol = 0.5 * fmax(S->atol, S->rtol * fabs(S->xlo));
            S->m = (S->c + S->b) * 0.5;
            S->mb = S->m - S->b;

            if (!(fabs(S->mb) > S->tol)) {
                S->next_state = 240;
            } else {
                if (!(S->ext > 3)) {
                    S->tol = (S->mb < 0 ? -fabs(S->tol) : fabs(S->tol));
                    S->p = (S->b - S->a)*S->fb;
                    if (S->first) {
                        S->q = S->fa - S->fb;
                        S->first = 0;
                    } else {
            // 120
                        S->fdb = (S->fd - S->fb)/(S->d - S->b);
                        S->fda = (S->fd - S->fa)/(S->d - S->a);
                        S->p = S->fda*S->p;
                        S->q = S->fdb*S->fa - S->fda*S->fb;
                    }
            // 130
                    if (S->p < 0.) {
                        S->p = -S->p;
                        S->q = -S->q;
                    }
            // 140
                    if (S->ext == 3) {
                        S->p *= 2.;
                    }
                    if  (!((S->p*1. == 0.) || (S->p <= S->q*S->tol))) {
            // 150
                        if (!(S->p < (S->mb*S->q))) {
            // 160
                            S->w = S->mb;
                        } else {
                            S->w = S->p / S->q;
                        }
                    } else {
                        S->w = S->tol;
                    }
                } else {
                    S->w = S->mb;
                }
            // 190
                S->d = S->a;
                S->fd = S->fa;
                S->a = S->b;
                S->fa = S->fb;
                S->b = S->b + S->w;
                S->xlo = S->b;
                S->x = S->xlo;
                S->next_state = 200;
                S->status = 1;
                return;
            }
        } else if (S->next_state == 200) {
            S->fb = S->fx;
            if (!(S->fc*S->fb >= 0.)) {
            // 210
                if (!(S->w == S->mb)) {
            // 220
                    S->ext += 1;
                } else {
                    S->ext = 0;
                }
            } else {
                // 70
                S->c = S->a;
                S->fc = S->fa;
                S->ext = 0;
            }
            // 230
            S->next_state = 80;

        } else if (S->next_state == 240) {
            S->xhi = S->c;
            S->qrzero = (((S->fc >= 0.) && (S->fb <= 0.)) ||
                        ((S->fc < 0.) && (S->fb >= 0.)));

            if (S->qrzero) {
                S->status = 0;
            } else {
            // 250
                S->status = -1;
            }
            return;

        } else {
            // Bad state, should not be possible to get here
            S->status = -9999;
            return;
        }
    }
}


double cdflib_erf(double x)
{
    //    Evaluation of the real error function

    double ax, bot, t, top;
    double c = .564189583547756;
    double a[5] = {.771058495001320e-04, -.133733772997339e-02,
                   .323076579225834e-01, .479137145607681e-01,
                   .128379167095513e+00};
    double b[3] = {.301048631703895e-02, .538971687740286e-01,
                   .375795757275549e+00};
    double p[8] = {-1.36864857382717e-07, 5.64195517478974e-01,
                   7.21175825088309e+00, 4.31622272220567e+01,
                   1.52989285046940e+02, 3.39320816734344e+02,
                   4.51918953711873e+02, 3.00459261020162e+02};
    double q[8] = {1.00000000000000e+00, 1.27827273196294e+01,
                   7.70001529352295e+01, 2.77585444743988e+02,
                   6.38980264465631e+02, 9.31354094850610e+02,
                   7.90950925327898e+02, 3.00459260956983e+02};
    double r[5] = {2.10144126479064e+00, 2.62370141675169e+01,
                   2.13688200555087e+01, 4.65807828718470e+00,
                   2.82094791773523e-01};
    double s[4] = {9.41537750555460e+01, 1.87114811799590e+02,
                   9.90191814623914e+01, 1.80124575948747e+01};

    ax = fabs(x);
    if (ax <= 0.5) {
        t = x*x;
        top = ((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]) + 1.0;
        bot = ((b[0]*t+b[1])*t+b[2])*t + 1.0;
        return x*(top/bot);
    }
    if (ax <= 4.0) {
        top = (((((((p[0]
                    )*ax+p[1]
                   )*ax+p[2]
                  )*ax+p[3]
                 )*ax+p[4]
                )*ax+p[5]
               )*ax+p[6]
              )*ax + p[7];
        bot = (((((((q[0]
                    )*ax+q[1]
                   )*ax+q[2]
                  )*ax+q[3]
                 )*ax+q[4]
                )*ax+q[5]
               )*ax+q[6])*ax + q[7];
        t = 0.5 + (0.5 - exp(-x*x)*(top/bot));
        return (x < 0 ? -t : t);
    }
    if (ax < 5.8) {
        t = pow(1/x, 2);
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t + r[4];
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t + 1.0;
        t = 0.5 + (0.5 - exp(-x*x) * (c - top/(x*x*bot))/ax);
        return (x < 0 ? -t : t);
    }
    return (x < 0 ? -1 : 1);
}


double erfc1(int ind, double x)
{
//    Evaluation of the complementary error function
//
//    Erfc1(ind,x) = erfc(x)            if ind = 0
//    Erfc1(ind,x) = exp(x*x)*erfc(x)   otherwise

    double ax, bot, t, top, result;
    double c = 0.564189583547756;
    double a[5] = {.771058495001320e-04, -.133733772997339e-02,
                   .323076579225834e-01, .479137145607681e-01,
                   .128379167095513e+00};
    double b[3] = {.301048631703895e-02, .538971687740286e-01,
                   .375795757275549e+00};
    double p[8] = {-1.36864857382717e-07, 5.64195517478974e-01,
                   7.21175825088309e+00, 4.31622272220567e+01,
                   1.52989285046940e+02, 3.39320816734344e+02,
                   4.51918953711873e+02, 3.00459261020162e+02};
    double q[8] = {1.00000000000000e+00, 1.27827273196294e+01,
                   7.70001529352295e+01, 2.77585444743988e+02,
                   6.38980264465631e+02, 9.31354094850610e+02,
                   7.90950925327898e+02, 3.00459260956983e+02};
    double r[5] = {2.10144126479064e+00, 2.62370141675169e+01,
                   2.13688200555087e+01, 4.65807828718470e+00,
                   2.82094791773523e-01};
    double s[4] = {9.41537750555460e+01, 1.87114811799590e+02,
                   9.90191814623914e+01, 1.80124575948747e+01};

    if (x <= -5.6) { return (ind == 0 ? 2.0 : (2*exp(x*x))); }

    // sqrt(log(np.finfo(np.float64).max)) ~= 26.64
    if ((ind == 0) && (x > 26.64))  { return 0.0; }

    ax = fabs(x);

    if (ax <= 0.5) {
        t = x*x;
        top = (((((a[0])*t+a[1])*t+a[2])*t+a[3])*t+a[4]) + 1.0;
        bot = (((b[0])*t+b[1])*t+b[2])*t + 1.0;
        result = 0.5 + (0.5 - x*(top/bot));
        return (ind == 0 ? result : result*exp(t));
    } else if ((0.5 < ax) && (ax <= 4.0)) {
        top = (((((((p[0]
                    )*ax+p[1]
                   )*ax+p[2]
                  )*ax+p[3]
                 )*ax+p[4]
                )*ax+p[5]
               )*ax+p[6]
              )*ax + p[7];
        bot = (((((((q[0]
                  )*ax+q[1]
                 )*ax+q[2]
                )*ax+q[3]
               )*ax+q[4]
              )*ax+q[5]
             )*ax+q[6])*ax + q[7];
        result = top / bot;
    } else {
        t = pow(1 / x, 2);
        top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t + r[4];
        bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t + 1.0;
        result = (c - t*(top/bot)) / ax;
    }
    if (ind == 0) {
        result *= exp(-(x*x));
        return (x < 0 ? (2.0 - result) : result);
    } else {
        return (x < 0 ? (2.0*exp(x*x) - result) : result);
    }
}


double esum(int mu, double x)
{
    //    Evaluation of exp(mu + x)

    if (x > 0.0) {
        if ((mu > 0.) || (mu + x < 0)) {
            return exp(mu)*exp(x);
        } else {
            return exp(mu + x);
        }
    } else {
        if ((mu < 0.) || (mu + x > 0.)) {
            return exp(mu)*exp(x);
        } else {
            return exp(mu + x);
        }
    }
}


double fpser(double a, double b, double x, double eps)
{
    //           Evaluation of i_x(a,b)
    //
    //    for b < Min(eps,eps*a) and x <= 0.5.

    double an, c, s, t, tol, result = 1.0;

    if (!(a <= 1e-3*eps)){
        result = 0.0;
        t = a*log(x);
        if (t < -708.) {return result;}
        result = exp(t);
    }
    //  Note that 1/Beta(a,b) = b
    result *= (b / a);
    tol = eps /a;
    an = a + 1.0;
    t = x;
    s = t / an;
    while (1) {
        an += 1;
        t *= x;
        c = t / an;
        s += c;
        if (!(fabs(c) > tol)) { break; }
    }
    return result*(1. + a*s);
}


double gam1(double a)
{
    //    Computation of 1/gamma(a+1) - 1  for -0.5 <= A <= 1.5

    double bot, d, t, top, w;
    double p[7] = {.577215664901533e+00, -.409078193005776e+00,
                   -.230975380857675e+00, .597275330452234e-01,
                   .766968181649490e-02, -.514889771323592e-02,
                   .589597428611429e-03};
    double q[5] = {.100000000000000e+01, .427569613095214e+00,
                   .158451672430138e+00, .261132021441447e-01,
                   .423244297896961e-02};
    double r[9] = {-.422784335098468e+00, -.771330383816272e+00,
                   -.244757765222226e+00, .118378989872749e+00,
                   .930357293360349e-03, -.118290993445146e-01,
                   .223047661158249e-02, .266505979058923e-03,
                   -.132674909766242e-03};
    double s[2] = {.273076135303957e+00, .559398236957378e-01};

    d = a - 0.5;
    t = (d > 0 ? d - 0.5 : a);

    if (t == 0.0) { return 0.0; }

    if (t < 0) {
        top = ((((((((r[8]
                     )*t+r[7]
                    )*t+r[6]
                   )*t+r[5]
                  )*t+r[4]
                 )*t+r[3]
                )*t+r[2]
               )*t+r[1]
              )*t + r[0];
        bot = (s[1]*t + s[0])*t + 1.0;
        w = top / bot;
        if (d > 0.0) {
            return t*w/a;
        } else {
            return a * ((w + 0.5) + 0.5);
        }
    }
    top = ((((((p[6]
               )*t+p[5]
              )*t+p[4]
             )*t+p[3]
            )*t+p[2]
           )*t+p[1]
          )*t + p[0];
    bot = ((((q[4]
             )*t+q[3]
            )*t+q[2]
           )*t+q[1]
          )*t + 1.0;
    w = top / bot;
    if (d > 0.0) {
        return (t/a) * ((w - 0.5) - 0.5);
    } else {
        return a * w;
    }
}


struct TupleDI gaminv(double a, double p, double q, double x0)
{
    //        INVERSE INCOMPLETE GAMMA RATIO FUNCTION
    //
    //    GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
    //    THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
    //    ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
    //    TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
    //    PARTICULAR COMPUTER ARITHMETIC BEING USED.
    //
    //                    ------------
    //
    //    X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
    //    AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
    //    NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
    //    A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
    //    IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.
    //
    //    X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER
    //    DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET
    //    X0 <= 0.
    //
    //    IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
    //    WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
    //    VALUES ...
    //
    //    IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
    //                NOT USED.
    //    IERR>0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
    //                WERE PERFORMED.
    //    IERR = -2    (INPUT ERROR) A <= 0
    //    IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
    //                IS TOO LARGE.
    //    IERR = -4    (INPUT ERROR) P + Q .NE. 1
    //    IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
    //                RECENT VALUE OBTAINED FOR X IS GIVEN.
    //                THIS CANNOT OCCUR IF X0 <= 0.
    //    IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
    //                THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
    //    IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
    //                ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
    //                ITERATION CANNOT BE PERFORMED IN THIS
    //                CASE. IF X0 <= 0, THIS CAN OCCUR ONLY
    //                WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
    //                POSITIVE THEN THIS CAN OCCUR WHEN A IS
    //                EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
    //                LARGE (SAY A >= 1.E20).

    double act_val, ap1, ap2, ap3, apn, b, bot, d, g, h;
    double pn, qg, qn, r, rta, s, s2, ssum, t, top, tol_val, u, w, y, z;
    double am1 = 0.;
    double ln10 = log(10);
    double c = 0.57721566490153286060651209008;
    double tol = 1e-5;
    double e = spmpar[0];
    double e2 = 2*e;
    double amax = 0.4e-10 / (e*e);
    double xmin = spmpar[1];
    double xmax = spmpar[2];
    double xn = x0;
    double x = 0.;
    int ierr = 0;
    int iop = (e > 1e-10 ? 1 : 0);
    int use_p = (p > 0.5 ? 0 : 1);
    int skip140 = 0;
    double arr[4] = {3.31125922108741, 11.6616720288968,
                     4.28342155967104, 0.213623493715853};
    double barr[4] = {6.61053765625462, 6.40691597760039,
                      1.27364489782223, 0.036117081018842};
    double eps0[2] = {1e-10, 1e-8};
    double amin[2] = {500., 100.};
    double bmin[2] = {1e-28, 1e-13};
    double dmin[2] = {1e-6, 1e-4};
    double emin[2] = {2e-3, 6e-3};
    double eps = eps0[iop];

    if (!(a > 0.)){return (struct TupleDI){.d1=x, .i1=-2};}
    t = p + q - 1.0;
    if (!(fabs(t) <= e)) {return (struct TupleDI){.d1=x, .i1=-4};}
    ierr = 0;
    if (p == 0.0) {return (struct TupleDI){.d1=x, .i1=0};}
    if (q == 0.0) {return (struct TupleDI){.d1=xmax, .i1=0};}
    if (a == 1.) {
        return (struct TupleDI){.d1=(!(q >= 0.9)? -log(q) : -alnrel(-p)), .i1=0};
    }
    if (x0 > 0.0) {
        use_p = (p > 0.5 ? 0 : 1);
        am1 = (a - 0.5) - 0.5;
        if ((use_p ? p : q) <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
    } else if (a <= 1.0) {
        g = cdflib_gamma(a + 1.);
        qg = q*g;
        if (qg == 0.0) {return (struct TupleDI){.d1=xmax, .i1=-8};}
        b = qg / a;

        if ((qg > 0.6*a) || (((a >= 0.3) || (b < 0.35)) && (b >= 0.45))) {
            if (b*q > 1.e-8) {
    // 50
                if (p <= 0.9) {
                    xn = exp(log(p*g)/a);
                } else {
    // 60
                    xn = exp((alnrel(-q) + gamln1(a))/a);
                }
            } else {
    // 40
                xn = exp(-(q/a + c));
            }
    // 70
            if (xn == 0.0) {return (struct TupleDI){.d1=x, .i1=-3};}
            t = 0.5 + (0.5 - xn/(a + 1.));
            xn /= t;
            // 160
            use_p = 1;
            am1 = (a - 0.5) - 0.5;
            if ((use_p ? p : q) <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
        } else if ((a >= 0.3) || (b < 0.35)) {
    // 10
            if (b == 0.0) {return (struct TupleDI){.d1=xmax, .i1=-8};}
            y = -log(b);
            s = 0.5 + (0.5 - a);
            z = log(y);
            t = y - s*z;
            if (b < 0.15) {
    // 20
                if (b <= 0.01) {
                    xn = gaminv_helper_30(a, s, y, z);
                    if ((a <= 1.) || (b > bmin[iop])) {
                        return (struct TupleDI){.d1=xn, .i1=0};
                    }
                }
                u = ((t+2.*(3.-a))*t + (2.-a)* (3.-a))/((t+ (5.-a))*t+2.);
                xn = y - s*log(t) - log(u);
            } else {
                xn = y - s*log(t)-log(1.+s/(t+1.));
            }
            // 220
            use_p = 0;
            am1 = (a - 0.5) - 0.5;
            if (q <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
        } else {
            t = exp(-(b+c));
            u = t*exp(t);
            xn = t*exp(u);

            // 160
            use_p = 1;
            am1 = (a - 0.5) - 0.5;
            if (p <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
        }

    } else if (a > 1) {
        w = (q > 0.5 ? log(p) : log(q));
        t = sqrt(-2.*w);
        top = ((arr[3]*t+arr[2])*t+arr[1])*t+arr[0];
        bot = (((barr[3]*t+barr[2])*t+barr[1])*t+barr[0])*t+1.0;
        s = t - top/bot;
        if (q > 0.5) { s = -s; }
        rta = sqrt(a);
        s2 = s*s;
        xn = a + s*rta;
        xn += (s2 - 1.)/3. + s*(s2-7.)/(36.*rta);
        xn -= ((3.*s2+7.)*s2-16.)/(810.*a);
        xn += s*((9.*s2+256.)*s2-433.)/(38880.*a*rta);
        xn = (xn > 0.0 ? xn : 0);

        if (a >= amin[iop]) {
            x = xn;
            d = 0.5 + (0.5 - xn/a);
            if (fabs(d) <= dmin[iop]) {return (struct TupleDI){.d1=xn, .i1=0};}
        }
        if (p > 0.5) {
            // All branches that go to 220
            if (xn < 3.*a) {
                ; // Do nothing
            } else {
                y = - (w + gamln(a));
                d = fmax(2., a*(a - 1.));
                if (y < ln10*d) {
                    t = a - 1.;
                    // Note: recursion not duplicate
                    xn = y + t*log(xn) - alnrel(-t/(xn+1.));
                    xn = y + t*log(xn) - alnrel(-t/(xn+1.));
                } else {
                    s = 1. - a;
                    z = log(y);
                    xn = gaminv_helper_30(a, s, y, z);
                    if (a <= 1.) {return (struct TupleDI){.d1=xn, .i1=0};}
                }
            }
            // Go to 220
            use_p = 0;
            am1 = (a - 0.5) - 0.5;
            if (q <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
        } else {
            // all branches that go to 170
            ap1 = a + 1.0;
            if (xn > 0.7*ap1) {
                ; // Do nothing
            } else {
                // branches that go to 170 via 140
                w += gamln(ap1);
                if (xn <= 0.15*ap1) {
                    ap2 = a + 2.;
                    ap3 = a + 3.;
                    // Note: recursion not duplicate
                    x = exp((w+x)/a);
                    x = exp((w+x-log(1.+(x/ap1)*(1.+x/ap2)))/a);
                    x = exp((w+x-log(1.+(x/ap1)*(1.+x/ap2)))/a);
                    x = exp((w+x-log(1.+(x/ap1)*(1.+(x/ap2)*(1.+x/ap3))))/a);
                    xn = x;
                    if (xn <= 0.01*ap1){
                        if (x <= emin[iop]*ap1) {
                            return (struct TupleDI){.d1=x, .i1=0};
                        } else {
                            use_p = 1;
                            am1 = (a - 0.5) - 0.5;
                            if (p <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
                            skip140 = 1;
                        }
                    }
                }
                if (!(skip140)) {
                    // Go to 140
                    apn = ap1;
                    t = xn/apn;
                    ssum = 1. + t;
                    while (1) {
                        t *= xn/apn;
                        ssum += t;
                        if (!(t > 1e-4)) { break; }
                    }
                    t = w - log(ssum);
                    xn = exp((xn + t)/a);
                    xn *= (1. - (a*log(xn) - xn -t) / (a - xn));
                }
            }
            // Go to 170
            use_p = 1;
            am1 = (a - 0.5) - 0.5;
            if (p <= 1.e10*xmin) {return (struct TupleDI){.d1=xn, .i1=-8};}
        }
    }
    // Schroder iteration using P or Q
    for (ierr = 0; ierr < 20; ierr++)
    {
        if (a > amax) {
            d = 0.5 + (0.5 - xn / a);
            if (fabs(d) <= e2) {return (struct TupleDI){.d1=xn, .i1=-8};}
        }
        struct TupleDD pnqn = gratio(a, xn, 0);
        pn = pnqn.d1;
        qn = pnqn.d2;
        if ((pn == 0.) || (qn == 0.)) {return (struct TupleDI){.d1=xn, .i1=-8};}
        r = rcomp(a, xn);
        if (r == 0.) {return (struct TupleDI){.d1=xn, .i1=-8};}
        t =  (use_p ? (pn-p)/r : (q-qn)/r);
        w = 0.5 * (am1 - xn);

        if ((fabs(t) <= 0.1) && (fabs(w*t) <= 0.1)) {
            h = t * (1. + w*t);
            x = xn * (1. - h);
            if (x <= 0.) {return (struct TupleDI){.d1=x, .i1=-7};}
            if ((fabs(w) >= 1.) && (fabs(w)*t*t <= eps)) {
                return (struct TupleDI){.d1=x, .i1=ierr};
            }
            d = fabs(h);
        } else {
            x = xn * (1. - t);
            if (x <= 0.) {return (struct TupleDI){.d1=x, .i1=-7};}
            d = fabs(t);
        }

        xn = x;

        if (d <= tol) {
            act_val = (use_p ? fabs(pn-p) : fabs(qn-q));
            tol_val = tol*(use_p ? p : q);
            if ((d <= eps) || (act_val <= tol_val)){
                return (struct TupleDI){.d1=x, .i1=ierr};
            }
        }
    }
    return (struct TupleDI){.d1=x, .i1=-6};
}


double gaminv_helper_30(double a, double s, double y, double z)
{
    double c1, c2, c3, c4, c5;
    c1 = -s*z;
    c2 = -s*(1. + c1);
    c3 = s*((0.5*c1+ (2.-a))*c1 + (2.5-1.5*a));
    c4 = -s*(((c1/3. + (2.5-1.5*a))*c1 + ((a-6.)*a+7.))*c1 + ((11.*a-46.)*a+47.)/6.);
    c5 = -s*((((-c1/4.+ (11.*a-17.)/6.
               )*c1+ ((-3.*a+13.)*a-13.)
              )*c1 + 0.5*(((2.*a-25.)*a+72.)*a-61.)
             )*c1+ (((25.*a-195.)*a+477.)*a-379.)/12.);
    return ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y;
}


double gamln(double a)
{
    //    Evaluation of ln(gamma(a)) for positive a

    double t, w, d = .418938533204673;
    int i,n;
    const double c[6] = {.833333333333333e-01, -.277777777760991e-02,
                         .793650666825390e-03, -.595202931351870e-03,
                         .837308034031215e-03, -.165322962780713e-02};

    if (a <= 0.8) { return gamln1(a) - log(a); }

    if (a <= 2.25) {
        t = (a-0.5) - 0.5;
        return gamln1(t);
    }

    if (a < 10) {
        n = (int)(a - 1.25);
        t = a;
        w = 1.0;
        for (i = 0; i < n; i++)
        {
            t -= 1.0;
            w *= t;
        }
        return gamln1(t-1.) + log(w);
    }
    t = pow(1/a, 2);
    w = (((((c[5]*t+c[4])*t+c[3])*t+c[2])*t+c[1])*t+c[0])/a;
    return (d + w) + (a-0.5)*(log(a) - 1.);
}


double gamln1(double a)
{
    //    Evaluation of ln(gamma(1 + a)) for -0.2 <= A <= 1.25

    double bot, top, w, x;

    const double p[7] = { .577215664901533e+00,  .844203922187225e+00,
                         -.168860593646662e+00, -.780427615533591e+00,
                         -.402055799310489e+00, -.673562214325671e-01,
                         -.271935708322958e-02};
    const double q[6] = {.288743195473681e+01, .312755088914843e+01,
                         .156875193295039e+01, .361951990101499e+00,
                         .325038868253937e-01, .667465618796164e-03};
    const double r[6] = {.422784335098467e+00, .848044614534529e+00,
                         .565221050691933e+00, .156513060486551e+00,
                         .170502484022650e-01, .497958207639485e-03};
    const double s[5] = {.124313399877507e+01, .548042109832463e+00,
                         .101552187439830e+00, .713309612391000e-02,
                         .116165475989616e-03};

    if (a < 0.6) {
        top = ((((((p[6]
                   )*a+p[5]
                  )*a+p[4]
                 )*a+p[3]
                )*a+p[2]
               )*a+p[1]
              )*a+p[0];
        bot = ((((((q[5]
                   )*a+q[4]
                  )*a+q[3]
                 )*a+q[2]
                )*a+q[1]
               )*a+q[0]
              )*a+1.;
        w = top/bot;
        return -a*w;
    } else {
        x = (a - 0.5) - 0.5;
        top = (((((r[5]
                  )*x+r[4]
                 )*x+r[3]
                )*x+r[2]
               )*x+r[1]
              )*x+r[0];
        bot = (((((s[4]
                  )*x+s[3]
                 )*x+s[2]
                )*x+s[1]
               )*x+s[0]
              )*x+1.;
        w = top/bot;
        return x*w;
    }
}


double cdflib_gamma(double a)
{
    //        Evaluation of the gamma function for real arguments
    //
    //    Gamma(a) is assigned the value 0 when the gamma function cannot
    //    be computed.

    double bot, g, lnx, t, top, w, z, result;
    int i, j, m, n;
    double s = 0.0;
    double d = 0.5*(log(2.*PI) - 1);
    double x = a;
    double p[7] = {.539637273585445e-03, .261939260042690e-02,
                   .204493667594920e-01, .730981088720487e-01,
                   .279648642639792e+00, .553413866010467e+00,
                   1.0};
    double q[7] = {-.832979206704073e-03, .470059485860584e-02,
                   .225211131035340e-01, -.170458969313360e+00,
                   -.567902761974940e-01, .113062953091122e+01,
                   1.0};
    double r[5] = {.820756370353826e-03, -.595156336428591e-03,
                   .793650663183693e-03, -.277777777770481e-02,
                   .833333333333333e-01};

    result = 0.0;
    if (fabs(a) < 15) {
        t = 1.0;
        m = (int)(a) - 1;
        if (m > 0) {
            for (j = 0; j < m; j++)
            {
                x -= 1.0;
                t *= x;
            }
            x -= 1.0;
        } else if (m == 0) {
            x -= 1.0;
        } else {
            t = a;
            if (a <= 0.) {
                m = -m - 1;
                if (m != 0.) {
                    for (j = 0; j < m; j++)
                    {
                        x += 1.0;
                        t *= x;
                    }
                }
                x += 0.5;
                x += 0.5;
                t *= x;

                if (t == 0.) { return result; }
            }
            if (fabs(t) < 1e-30) {
                if (fabs(t)*spmpar[2] <= 1.0001) { return result; }
                return 1./t;
            }
        }
        top = p[0];
        bot = q[0];
        for (i = 1; i < 7; i++)
        {
            top *= x;
            top += p[i];
            bot *= x;
            bot += q[i];
        }
        result = top / bot;
        return (a < 1.0 ? result/t : result*t);
    }

    if (fabs(a) >= 1.e3) { return result; }

    if (a <= 0.0) {
        x = -a;
        n = (int)x;
        t = x - n;
        if (t > 0.9) {
            t = 1. - t;
        }
        s = sin(PI*t) / PI;
        if (n % 2 == 0) {
            s = -s;
        }
        if (s == 0.0) { return result; }
    }
    t = pow(1 / x, 2);
    g = ((((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4]) / x;
    lnx = log(x);
    z = x;
    g = (d + g) + (z -0.5)*(lnx - 1.);
    w = g;
    t = g - w;
    if (w > 0.99999*709) { return result; }
    result = exp(w)*(1. + t);
    return (a < 0.0 ? (1. / (result * s)) / x : result);
}


struct TupleDD grat1(double a, double x, double r, double eps)
{
    //        Evaluation of the incomplete gamma ratio functions
    //                    p(a,x) and q(a,x)
    //
    //    it is assumed that a <= 1.  Eps is the tolerance to be used.
    //    the input argument r has the value e**(-x)*x**a/gamma(a).

    double a2n, a2nm1, am0, an, an0, b2n, b2nm1, c, cma, g, h, j, l;
    double p, q, ssum, t, tol, w, z;

    if (a*x == 0.) {
        if (x > a) {
            return (struct TupleDD){.d1=1.0, .d2=0.0};
        } else {
            return (struct TupleDD){.d1=0.0, .d2=1.0};
        }
    }

    if (a == 0.5) {
        if (x < 0.25) {
            p = cdflib_erf(sqrt(x));
            return (struct TupleDD){.d1=p, .d2=0.5 + (0.5 - p)};
        } else {
            q = erfc1(0, sqrt(x));
            return (struct TupleDD){.d1=0.5 + (0.5 - q), .d2=q};
        }
    }

    if (x < 1.1) {
        //
        // Taylor series for p(a,x)/x**a
        //
        an = 3.0;
        c = x;
        ssum = x / (a + 3.);
        tol = 0.1*eps / (a + 1.);
        while (1) {
            an += 1;
            c *= -(x / an);
            t = c / (a + an);
            ssum += t;
            if (fabs(t) <= tol) { break; }
        }
        j = a*x*((ssum/6. - 0.5/(a+2.))*x + 1./(a+1.));

        z = a * log(x);
        h = gam1(a);
        g = 1. + h;

        if (((x >= 0.25) && (a >= x /2.59)) || ((x < 0.25) && (z <= -0.13394))) {
            w = exp(z);
            p = w*g*(0.5 + (0.5 - j));
            q = 0.5 + (0.5 - p);
            return (struct TupleDD){.d1=p, .d2=q};
        } else {
            l = rexp(z);
            w = 0.5 + (0.5 + l);
            q = (w*j - l)*g - h;
            if (q < 0.0) {
                return (struct TupleDD){.d1=1.0, .d2=0.0};
            }
            p = 0.5 + (0.5 - q);
            return (struct TupleDD){.d1=p, .d2=q};
        }
    }
    //  Continued fraction expansion
    a2nm1 = 1.0;
    a2n = 1.0;
    b2nm1 = x;
    b2n = x + (1. - a);
    c = 1.0;
    while (1) {
        a2nm1 = x*a2n + c*a2nm1;
        b2nm1 = x*b2n + c*b2nm1;
        am0 = a2nm1/b2nm1;
        c = c + 1.;
        cma = c - a;
        a2n = a2nm1 + cma*a2n;
        b2n = b2nm1 + cma*b2n;
        an0 = a2n/b2n;
        if (!(fabs(an0-am0) >= eps*an0)) { break; }
    }
    q = r*an0;
    p = 0.5 + (0.5 - q);
    return (struct TupleDD){.d1=p, .d2=q};
}


struct TupleDD gratio(double a, double x, int ind)
{
    //    Evaluation of the incomplete gamma ratio functions
    //                    P(a,x) and Q(a,x)
    //
    //                    ----------
    //
    //    It is assumed that a and x are nonnegative, where a and x
    //    Are not both 0.
    //
    //    Ans and qans are variables. Gratio assigns ans the value
    //    P(a,x) and qans the value q(a,x). Ind may be any integer.
    //    If ind = 0 then the user is requesting as much accuracy as
    //    Possible (up to 14 significant digits). Otherwise, if
    //    Ind = 1 then accuracy is requested to within 1 unit of the
    //    6-Th significant digit, and if ind .ne. 0,1 Then accuracy
    //    Is requested to within 1 unit of the 3rd significant digit.
    //
    //    Error return ...
    //    Ans is assigned the value 2 when a or x is negative,
    //    When a*x = 0, or when p(a,x) and q(a,x) are indeterminant.
    //    P(a,x) and q(a,x) are computationally indeterminant when
    //    X is exceedingly close to a and a is extremely large.

    double d10 = -.185185185185185e-02;
    double d20 = .413359788359788e-02;
    double d30 = .649434156378601e-03;
    double d40 = -.861888290916712e-03;
    double d50 = -.336798553366358e-03;
    double d60 = .531307936463992e-03;
    double d70 = .344367606892378e-03;
    double alog10 = log(10);
    double rt2pi = sqrt(1. / (2.*PI));
    double rtpi = sqrt(PI);
    double eps = spmpar[0];
    double acc, a2n, a2nm1, am0, amn, an, an0, ans, apn, b2n, b2nm1;
    double c, c0, c1, c2, c3, c4, c5, c6, cma, e0, g, h, j, l, qans, r;
    double rta, rtx, s, ssum, t, t1, tol, twoa, u, w, x0, y, z;
    int i, m, n, last_entry;

    double wk[20] = {0.0};
    double acc0[3] = {5.e-15, 5.e-7, 5.e-4};
    double big[3] = {20., 14., 10.};
    double e00[3] = {.00025, .025, .14};
    double x00[3] = {31., 17., 9.7};
    double d0[13] = {.833333333333333e-01, -.148148148148148e-01,
                     .115740740740741e-02, .352733686067019e-03,
                     -.178755144032922e-03, .391926317852244e-04,
                     -.218544851067999e-05, -.185406221071516e-05,
                     .829671134095309e-06, -.176659527368261e-06,
                     .670785354340150e-08, .102618097842403e-07,
                     -.438203601845335e-08};
    double d1[12] = {-.347222222222222e-02, .264550264550265e-02,
                     -.990226337448560e-03, .205761316872428e-03,
                     -.401877572016461e-06, -.180985503344900e-04,
                     .764916091608111e-05, -.161209008945634e-05,
                     .464712780280743e-08, .137863344691572e-06,
                     -.575254560351770e-07, .119516285997781e-07};
    double d2[10] = {-.268132716049383e-02, .771604938271605e-03,
                     .200938786008230e-05, -.107366532263652e-03,
                     .529234488291201e-04, -.127606351886187e-04,
                     .342357873409614e-07, .137219573090629e-05,
                     -.629899213838006e-06, .142806142060642e-06};
    double d3[8] = {.229472093621399e-03, -.469189494395256e-03,
                    .267720632062839e-03, -.756180167188398e-04,
                    -.239650511386730e-06, .110826541153473e-04,
                    -.567495282699160e-05, .142309007324359e-05};
    double d4[6] = {.784039221720067e-03, -.299072480303190e-03,
                    -.146384525788434e-05, .664149821546512e-04,
                    -.396836504717943e-04, .113757269706784e-04};
    double d5[4] = {-.697281375836586e-04, .277275324495939e-03,
                    -.199325705161888e-03, .679778047793721e-04};
    double d6[2] = {-.592166437353694e-03, .270878209671804e-03};

    if ((!(a >= 0.0)) || (!(x >= 0.0))) {
        return (struct TupleDD){.d1=2.0, .d2=0.0};
    }
    if ((a == 0.0) && (x == 0.0)) {
        return (struct TupleDD){.d1=2.0, .d2=0.0};
    }

    if (a*x == 0.0) {
        if (x > a) {
            return (struct TupleDD){.d1=1.0, .d2=0.0};
        } else {
            return (struct TupleDD){.d1=0.0, .d2=1.0};
        }
    }

    if ((!(ind == 0)) && (!(ind == 1))) {
        ind = 2;
    }

    acc = fmax(acc0[ind], eps);
    e0 = e00[ind];
    x0 = x00[ind];

    if (a >= 1.0) {
        if (a >= big[ind]) {
            // 30
            l = x / a;
            if (l == 0.0) {
                // 370
                return (struct TupleDD){.d1=0.0, .d2=1.0};
            }
            s = 0.5 + (0.5 - l);
            z = rlog(l);
            if (z >= (700. / a)) {
                // 410
                if (fabs(s) <= 2.*eps) {
                    return (struct TupleDD){.d1=2.0, .d2=0.0};
                }
                if (x > a) {
                    return (struct TupleDD){.d1=1.0, .d2=0.0};
                } else {
                    return (struct TupleDD){.d1=0.0, .d2=1.0};
                }
            }
            y = a*z;
            rta = sqrt(a);
            if (fabs(s) <= (e0 / rta)) {
                // 330
                // TEMME EXPANSION FOR L = 1
                //
                if (a*eps*eps > 3.28e-3) {
                    return (struct TupleDD){.d1=2.0, .d2=0.0};
                }
                c = 0.5 + (0.5 - y);
                w = (0.5 - sqrt(y) * (0.5 + (0.5 - y/3.))/rtpi)/c;
                u = 1. / a;
                z = sqrt(z+z);
                if (l < 1.0) {
                    z = -z;
                }
                if (ind == 0) {
                    c0 = (((((((d0[6]
                               )*z+d0[5]
                              )*z+d0[4]
                             )*z+d0[3]
                            )*z+d0[2]
                           )*z+d0[1]
                          )*z+d0[0]
                         )*z - (1./3.);
                    c1 = ((((((d1[5]
                              )*z+d1[4]
                             )*z+d1[3]
                            )*z+d1[2]
                           )*z+d1[1]
                          )*z+d1[0]
                         )*z + d10;
                    c2 = (((((d2[4])*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z + d20;
                    c3 = ((((d3[3])*z+d3[2])*z+d3[1])*z+d3[0])*z + d30;
                    c4 = (d4[1]*z+d4[0])*z + d40;
                    c5 = (d5[1]*z+d5[0])*z + d50;
                    c6 = d6[0]*z + d60;
                    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
                } else if (ind == 1) {
                    c0 = (d0[1]*z+d0[0])*z - (1. / 3.);
                    c1 = d1[0]*z + d10;
                    t = (d20*u+c1)*u + c0;

                } else {
                    t = d0[0]*z - (1. / 3.);
                }
                if (l < 1.) {
                    ans = c * (w - rt2pi*t/rta);
                    qans = 0.5 + (0.5 - ans);
                } else {
                    qans = c * (w + rt2pi*t/rta);
                    ans = 0.5 + (0.5 - qans);
                }
                return (struct TupleDD){.d1=ans, .d2=qans};
            }
            if (fabs(s) <= 0.4) {
                // 270
                // GENERAL TEMME EXPANSION
                //
                if ((fabs(s) <= 2.*eps) && (a*eps*eps > 3.28e-3)) {
                    return (struct TupleDD){.d1=2.0, .d2=0.0};
                }
                c = exp(-y);
                w = 0.5*erfc1(1, sqrt(y));
                u = 1./a;
                z = sqrt(z+z);
                if (l < 1.) {
                    z = -z;
                }
                if (ind == 0) {
                    if (fabs(s) <= 1e-3) {
                        c0 = (((((((d0[6]
                                   )*z+d0[5]
                                  )*z+d0[4]
                                 )*z+d0[3]
                                )*z+d0[2]
                               )*z+d0[1]
                              )*z+d0[0]
                             )*z - (1./3.);
                        c1 = ((((((d1[5]
                                  )*z+d1[4]
                                 )*z+d1[3]
                                )*z+d1[2]
                               )*z+d1[1]
                              )*z+d1[0]
                             )*z + d10;
                        c2 = (((((d2[4])*z+d2[3])*z+d2[2])*z+d2[1])*z+d2[0])*z + d20;
                        c3 = ((((d3[3])*z+d3[2])*z+d3[1])*z+d3[0])*z + d30;
                        c4 = (d4[1]*z+d4[0])*z + d40;
                        c5 = (d5[1]*z+d5[0])*z + d50;
                        c6 = d6[0]*z + d60;
                        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
                    } else {
                        c0 = (((((((((((((d0[12]
                                         )*z+d0[11]
                                        )*z+d0[10]
                                       )*z+d0[9]
                                      )*z+d0[8]
                                     )*z+d0[7]
                                    )*z+d0[6]
                                   )*z+d0[5]
                                  )*z+d0[4]
                                 )*z+d0[3]
                                )*z+d0[2]
                               )*z+d0[1]
                              )*z+d0[0]
                             )*z - (1./3.);
                        c1 = ((((((((((((d1[11]
                                        )*z+d1[10]
                                       )*z+d1[9]
                                      )*z+d1[8]
                                     )*z+d1[7]
                                    )*z+d1[6]
                                   )*z+d1[5]
                                  )*z+d1[4]
                                 )*z+d1[3]
                                )*z+d1[2]
                               )*z+d1[1]
                              )*z+d1[0]
                             )*z + d10;
                        c2 = ((((((((((d2[9]
                                      )*z+d2[8]
                                     )*z+d2[7]
                                    )*z+d2[6]
                                   )*z+d2[5]
                                  )*z+d2[4]
                                 )*z+d2[3]
                                )*z+d2[2]
                               )*z+d2[1]
                              )*z+d2[0]
                             )*z + d20;
                        c3 = ((((((((d3[7]
                                    )*z+d3[6]
                                   )*z+d3[5]
                                  )*z+d3[4]
                                 )*z+d3[3]
                                )*z+d3[2]
                               )*z+d3[1]
                              )*z+d3[0]
                             )*z + d30;
                        c4 = ((((((d4[5]
                                  )*z+d4[4]
                                 )*z+d4[3]
                                )*z+d4[2]
                               )*z+d4[1]
                              )*z+d4[0]
                             )*z + d40;
                        c5 = (((d5[3]*z+d5[2])*z+d5[1])*z+d5[0])*z + d50;
                        c6 = (d6[1]*z+d6[0])*z + d60;
                        t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0;
                    }
                } else if (ind == 1) {
                    c0 = ((((((d0[5]
                               )*z+d0[4]
                              )*z+d0[3]
                             )*z+d0[2]
                            )*z+d0[1]
                           )*z+d0[0]
                          )*z - (1./3.);
                    c1 = (((d1[3]*z+d1[2])*z+d1[1])*z+d1[0])*z + d10;
                    c2 = d2[0]*z + d20;
                    t = (c2*u+c1)*u + c0;
                } else {
                    t = ((d0[2]*z+d0[1])*z+d0[0])*z - (1./3.);
                }
                if (l < 1.) {
                    ans = c * (w - rt2pi*t/rta);
                    qans = 0.5 + (0.5 - ans);
                } else {
                    qans = c * (w + rt2pi*t/rta);
                    ans = 0.5 + (0.5 - qans);
                }
                return (struct TupleDD){.d1=ans, .d2=qans};
            }

            t = pow(1 / a, 2);
            t1 = (((0.75*t-1.)*t+3.5)*t-105.0)/ (a*1260.0);
            t1 -= y;
            r = rt2pi*rta*exp(t1);

            //
            // 40
            //
            if (r == 0.0) {
                // 420
                if (x > a) {
                    return (struct TupleDD){.d1=1.0, .d2=0.0};
                } else {
                    return (struct TupleDD){.d1=0.0, .d2=1.0};
                }
            }

            if (x <= fmax(a, alog10)) {
                // 50
                // TAYLOR SERIES FOR P/R
                //
                apn = a + 1.;
                t = x / apn;
                wk[0] = t;
                last_entry = 0;
                for (n = 1; n < 20; n++) {
                    apn += 1.0;
                    t *= (x / apn);
                    if (t <= 1e-3) { break; }
                    wk[n] = t;
                    last_entry = n;
                }
                ssum = t;
                tol = 0.5 * acc;
                while (1) {
                    apn += 1.0;
                    t *= x / apn;
                    ssum += t;
                    if (!(t > tol)) { break; }
                }

                for (m = last_entry; m >= 0; m--)
                {
                    ssum += wk[m];
                }
                ans = (r/a) * (1.0 + ssum);
                qans = 0.5 + (0.5 - ans);
                return (struct TupleDD){.d1=ans, .d2=qans};
            }

            if (x < x0) {
                // 250
                // CONTINUED FRACTION EXPANSION
                //
                tol = fmax(5.0*eps, acc);
                a2nm1 = 1.0;
                a2n = 1.0;
                b2nm1 = x;
                b2n = x + (1.0 - a);
                c = 1.0;

                while (1) {
                    a2nm1 = x*a2n + c*a2nm1;
                    b2nm1 = x*b2n + c*b2nm1;
                    am0 = a2nm1/b2nm1;
                    c += 1.0;
                    cma = c - a;
                    a2n = a2nm1 + cma*a2n;
                    b2n = b2nm1 + cma*b2n;
                    an0 = a2n/b2n;
                    if (!(fabs(an0-am0) >= tol*an0)) { break; }
                }

                qans = r*an0;
                ans = 0.5 + (0.5 - qans);
                return (struct TupleDD){.d1=ans, .d2=qans};
            }
            // 100
            // ASYMPTOTIC EXPANSION
            //
            amn = a - 1.;
            t = amn / x;
            wk[0] = t;
            last_entry = 0;
            for (n = 1; n < 20; n++)
            {
                amn -= 1.0;
                t *= amn / x;
                if (fabs(t) <= 1e-3) { break; }
                wk[n] = t;
                // F77 code was using "n" later in the code
                last_entry = n;
            }
            ssum = t;

            while (!(fabs(t) <= acc)) {
                amn -= 1.0;
                t *= amn / x;
                ssum += t;
            }
            for (m = last_entry; m >= 0; m--)
            {
                ssum += wk[m];
            }

            qans = (r/x) * (1. + ssum);
            ans = 0.5 + (0.5 - qans);
            return (struct TupleDD){.d1=ans, .d2=qans};
        }

        twoa = a + a;
        m = (int)twoa;

        if (((a > x) || (x >= x0)) || (twoa != m)) {
            t1 = a*log(x) - x;
            r = exp(t1)/cdflib_gamma(a);

            //
            // 40 Again - This time coming from 20
            //
            if (r == 0.0) {
                if (x > a) {
                    return (struct TupleDD){.d1=1.0, .d2=0.0};
                } else {
                    return (struct TupleDD){.d1=0.0, .d2=1.0};
                }
            }

            if (x <= (a > alog10 ? a : alog10)) {
                // 50
                // TAYLOR SERIES FOR P/R
                //
                apn = a + 1.;
                t = x / apn;
                wk[0] = t;
                last_entry = 0;
                for (n =1; n < 20; n++) {
                    apn += 1.0;
                    t *= (x / apn);
                    if (t <= 1e-3) { break; }
                    wk[n] = t;
                    last_entry = n;
                }
                ssum = t;
                tol = 0.5 * acc;
                while (!(t <= tol)) {
                    apn += 1.0;
                    t *= x / apn;
                    ssum += t;
                }
                for (m = last_entry; m >= 0; m--)
                {
                    ssum += wk[m];
                }
                ans = (r/a) * (1. + ssum);
                qans = 0.5 + (0.5 - ans);
                return (struct TupleDD){.d1=ans, .d2=qans};
            }

            if (x < x0) {
                // 250
                // CONTINUED FRACTION EXPANSION
                //
                tol = (5.0*eps > acc ? 5.0*eps : acc);
                a2nm1 = 1.0;
                a2n = 1.0;
                b2nm1 = x;
                b2n = x + (1.0 - a);
                c = 1.0;

                while (1) {
                    a2nm1 = x*a2n + c*a2nm1;
                    b2nm1 = x*b2n + c*b2nm1;
                    am0 = a2nm1/b2nm1;
                    c += 1.0;
                    cma = c - a;
                    a2n = a2nm1 + cma*a2n;
                    b2n = b2nm1 + cma*b2n;
                    an0 = a2n/b2n;
                    if (!(fabs(an0-am0) >= tol*an0)) { break; }
                }

                qans = r*an0;
                ans = 0.5 + (0.5 - qans);
                return (struct TupleDD){.d1=ans, .d2=qans};
            }
            // 100
            // ASYMPTOTIC EXPANSION
            //
            amn = a - 1.;
            t = amn / x;
            wk[0] = t;
            last_entry = 0;
            for (n = 1; n < 20; n++)
            {
                amn -= 1.0;
                t *= amn / x;
                if (fabs(t) <= 1e-3) { break; }
                wk[n] = t;
                last_entry = n;
            }
            ssum = t;

            while (!(fabs(t) <= acc)) {
                amn -= 1.0;
                t *= amn / x;
                ssum += t;
            }
            for (m = last_entry; m >= 0; m--)
            {
                ssum += wk[m];
            }

            qans = (r/x) * (1. + ssum);
            ans = 0.5 + (0.5 - qans);
            return (struct TupleDD){.d1=ans, .d2=qans};
        }

        i = m / 2;
        if (a == i) {
            // 210
            ssum = exp(-x);
            t = ssum;
            n = 1;
            c = 0.;
        } else {
            // 220
            rtx = sqrt(x);
            ssum = erfc1(0, rtx);
            t = exp(-x) / (rtpi*rtx);
            n = 0;
            c = -0.5;
        }

        while (n < i) {
            // 230
            n += 1;
            c += 1.0;
            t *= x/c;
            ssum += t;
        }

        qans = ssum;
        ans = 0.5 + (0.5 - qans);
        return (struct TupleDD){.d1=ans, .d2=qans};
    }

    if (a == 0.5) {
        // 390
        if (x >= 0.5) {
            qans = erfc1(0, sqrt(x));
            ans = 0.5 + (0.5 - qans);
        } else {
            ans = cdflib_erf(sqrt(x));
            qans = 0.5 + (0.5 - ans);
        }
        return (struct TupleDD){.d1=ans, .d2=qans};
    }

    if (x < 1.1) {
        // 160
        // TAYLOR SERIES FOR P(A,X)/X**A
        //
        an = 3.0;
        c = x;
        ssum = x / (a + 3.);
        tol = 3.*acc / (a + 1.);
        while (1) {
            an += 1.0;
            c *= -(x / an);
            t = c / (a + an);
            ssum += t;
            if (!(fabs(t) > tol)) { break; }
        }
        j = a*x*((ssum / 6. - 0.5 / (a + 2.))*x + 1./(a + 1.));
        z = a*log(x);
        h = gam1(a);
        g = 1. + h;

        if (((x < 0.25) && (z > -0.13394)) || (a < x/2.59)) {
            l = rexp(z);
            w = 0.5 + (0.5 + l);
            qans = (w*j - l)*g - h;
            if (qans < 0.) {
                return (struct TupleDD){.d1=1.0, .d2=0.0};
            }
            ans = 0.5 + (0.5 - qans);
        } else {
            w = exp(z);
            ans = w*g*(0.5 + (0.5 - j));
            qans = 0.5 + (0.5 - ans);
        }
        return (struct TupleDD){.d1=ans, .d2=qans};
    }

    t1 = a*log(x) - x;
    u = a*exp(t1);
    if (u == 0.0) {
        return (struct TupleDD){.d1=1.0, .d2=0.0};
    }

    r = u * (1. + gam1(a));
    // 250
    // CONTINUED FRACTION EXPANSION
    //
    tol = (5.0*eps > acc ? 5.0*eps : acc);
    a2nm1 = 1.0;
    a2n = 1.0;
    b2nm1 = x;
    b2n = x + (1.0 - a);
    c = 1.0;

    while (1) {
        a2nm1 = x*a2n + c*a2nm1;
        b2nm1 = x*b2n + c*b2nm1;
        am0 = a2nm1/b2nm1;
        c += 1.0;
        cma = c - a;
        a2n = a2nm1 + cma*a2n;
        b2n = b2nm1 + cma*b2n;
        an0 = a2n/b2n;
        if (!(fabs(an0-am0) >= tol*an0)) { break; }
    }
    qans = r*an0;
    ans = 0.5 + (0.5 - qans);
    return (struct TupleDD){.d1=ans, .d2=qans};
}


double gsumln(double a, double b)
{
    //     Evaluation of the function ln(gamma(a + b))
    //     for 1 <= A <= 2  And  1 <= B <= 2

    double x;

    x = a + b - 2;
    if (x <= 0.25) {
        return gamln1(1. + x);
    }

    if (x <= 1.25) {
        return gamln1(x) + alnrel(x);
    }

    return gamln1(x - 1.) + log(x*(1. + x));
}


double psi(double xx)
{
    //                    Evaluation of the digamma function
    //
    //                          -----------
    //
    //    Psi(xx) is assigned the value 0 when the digamma function cannot
    //    be computed.
    //
    //    The main computation involves evaluation of rational chebyshev
    //    approximations published in math. Comp. 27, 123-127(1973) By
    //    cody, strecok and thacher.
    //
    //    ----------------------------------------------------------------
    //    Psi was written at Argonne National Laboratory for the FUNPACK
    //    package of special function subroutines. Psi was modified by
    //    A.H. Morris (nswc).

    double aug, den, dx0, sgn, upper, w, x, xmax1, xmx0, xsmall, z;
    double p1[7] = {0.895385022981970e-02, 0.477762828042627e+01,
                    0.142441585084029e+03, 0.118645200713425e+04,
                    0.363351846806499e+04, 0.413810161269013e+04,
                    0.130560269827897e+04};
    double q1[6] = {0.448452573429826e+02, 0.520752771467162e+03,
                    0.221000799247830e+04, 0.364127349079381e+04,
                    0.190831076596300e+04, 0.691091682714533e-05};
    double p2[4] = {-0.212940445131011e+01, -0.701677227766759e+01,
                    -0.448616543918019e+01, -0.648157123766197e+00};
    double q2[4] = {0.322703493791143e+02, 0.892920700481861e+02,
                    0.546117738103215e+02, 0.777788548522962e+01};
    int nq, i;
    dx0 = 1.461632144968362341262659542325721325;
    xmax1 = 4503599627370496.0;
    xsmall = 1e-9;
    x = xx;
    aug = 0.0;

    if (x < 0.5) {
        if (fabs(x) <= xsmall) {
            if (x == 0.) {
                return 0.0;
            }
            aug = -1./x;
        } else {
            // 10
            w = -x;
            sgn = PI / 4;
            if (w <= 0.) {
                w = -w;
                sgn = -sgn;
            }
            // 20
            if (w >= xmax1) {
                return 0.0;
            }
            w -= (int)w;
            nq = (int)(w*4.0);
            w = 4.*(w - 0.25*nq);

            if (nq % 2 == 1) {
                w = 1. - w;
            }
            z = (PI / 4.)*w;

            if ((nq / 2) % 2 == 1) {
                sgn = -sgn;
            }
            if ((((nq + 1) / 2) % 2) == 1) {
                aug = sgn * (tan(z)*4.);
            } else {
                if (z == 0.) {
                    return 0.0;
                }
                aug = sgn * (4./tan(z));
            }
        }
        x = 1 - x;
    }

    if (x <= 3.0) {
        // 50
        den = x;
        upper = p1[0]*x;
        for (i = 0; i < 5; i++)
        {
            den = (den + q1[i])*x;
            upper = (upper + p1[i+1])*x;
        }
        den = (upper + p1[6]) / (den + q1[5]);
        xmx0 = x - dx0;
        return (den * xmx0) + aug;
    } else {
        // 70
        if (x < xmax1) {
            w = 1. / (x*x);
            den = w;
            upper = p2[0]*w;

            for (i = 0; i < 3; i++) {
                den = (den + q2[i])*w;
                upper = (upper + p2[i+1])*w;
            }
            aug += upper / (den + q2[3]) - 0.5/x;
        }
        return aug + log(x);
    }
}


double rcomp(double a, double x)
{
    // Evaluation of exp(-x)*x**a/gamma(a)

    double t, t1, u, r2pi = sqrt(1. / (2.0*PI));
    if (a < 20) {
        t = a*log(x) - x;
        return (a < 1 ? a*exp(t)*(1. + gam1(a)) : exp(t) / cdflib_gamma(a));
    } else {
        u = x / a;
        if (u == 0.) {
            return 0.;
        }
        t = pow((1 / a), 2);
        t1 = (((0.75*t - 1.)*t + 3.5)*t - 105.) / (a*1260);
        t1 -= a*rlog(u);
        return r2pi*sqrt(a)*exp(t1);
    }
}


double rexp(double x)
{
    // Evaluation of exp(-x)*x**a/gamma(a)

    double p[2] = {.914041914819518e-09, .238082361044469e-01};
    double q[4] = {-.499999999085958e+00, .107141568980644e+00,
                   -.119041179760821e-01, .595130811860248e-03};
    double w;

    if (fabs(x) <= 0.15) {
        return x*(((p[1]*x+p[0])*x+1.)/((((q[3]*x+q[2])*x+q[1])*x+q[0])*x + 1.));
    } else {
        w = exp(x);
        if (x > 0.) {
            return w*(0.5+ (0.5 - 1./w));
        }  else {
            return (w - 0.5) - 0.5;
        }
    }
}


double rlog(double x)
{
    // Computation of  x - 1 - ln(x)

    double r, t, u, w, w1;
    double a = .566749439387324e-01;
    double b = .456512608815524e-01;
    double p[3] = {.333333333333333, -.224696413112536, .620886815375787e-02};
    double q[2] = {-.127408923933623e+01, .354508718369557};

    if ((x < 0.61) || (x > 1.57))
    {
        return ((x - 0.5) - 0.5) - log(x);
    }

    if (x < 0.82)
    {
        u = (x - 0.7) / 0.7;
        w1 = a - u*0.3;
    }
    else if (x > 1.18)
    {
        u = 0.75*x - 1.;
        w1 = b + u/3.;
    }
    else
    {
        u = (x - 0.5) - 0.5;
        w1 = 0.0;
    }

    r = u / (u + 2.);
    t = r*r;
    w = ((p[2]*t+p[1])*t+p[0])/ ((q[1]*t+q[0])*t+1.);
    return 2.*t*(1. / (1. - r) - r*w) + w1;
}


double rlog1(double x)
{
    // Evaluation of the function x - ln(1 + x)

    double a = .566749439387324e-01;
    double b = .456512608815524e-01;
    double p0 = .333333333333333e+00;
    double p1 = -.224696413112536e+00;
    double p2 = .620886815375787e-02;
    double q1 = -.127408923933623e+01;
    double q2 = .354508718369557e+00;
    double h, r, t, w, w1;

    if ((-0.39 <= x) && (x <= 0.57)) {
        if ((-0.18 <= x) && (x <= 0.18))
        {
            h = x;
            w1 = 0.0;
        } else if (x < -0.18)
        {
            h = (x + 0.3)/0.7;
            w1 = a - h*0.3;
        } else
        {
            // 0.57 >= x > 0.18
            h = 0.75*x - 0.25;
            w1 = b + h/3.0;
        }

        r = h / (h + 2);
        t = r*r;
        w = ((p2*t + p1)*t + p0) / ((q2*t + q1)*t + 1.);
        return 2*t*(1./(1.-r) - r*w) + w1;

    } else {

       return x - log((x + 0.5) + 0.5);

    }
}


double stvaln(double p)
{
    //                STarting VALue for Newton-Raphon
    //            calculation of Normal distribution Inverse
    //
    //                          Function
    //
    // Returns X such that CUMNOR(X) = P,  i.e., the  integral from -
    // infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
    //
    //                          Arguments
    //
    // P --> The probability whose normal deviate is sought.
    //                P is DOUBLE PRECISION
    //
    //                          Method
    //
    // The  rational   function   on  page 95    of Kennedy  and  Gentle,
    // Statistical Computing, Marcel Dekker, NY , 1980.

    double y, z;
    double xnum[5] = {-0.322232431088, -1.000000000000, -0.342242088547,
                      -0.204231210245e-1, -0.453642210148e-4};
    double xden[5]  = {0.993484626060e-1, 0.588581570495, 0.531103462366,
                       0.103537752850, 0.38560700634e-2};
    z = (p > 0.5 ? 1.0 - p : p);
    y = sqrt(-2.0 * log(z));
    z = y + (devlpl(xnum, 5, y) / devlpl(xden, 5, y));
    return (p > 0.5 ? z : -z);
}
