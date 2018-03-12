/*                                                     ellpj.c
 *
 *     Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * double u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
 * and dn(u|m) of parameter m between 0 and 1, and real
 * argument u.
 *
 * These functions are periodic, with quarter-period on the
 * real axis equal to the complete elliptic integral
 * ellpk(m).
 *
 * Relation to incomplete elliptic integral:
 * If u = ellik(phi,m), then sn(u|m) = sin(phi),
 * and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
 *
 * Computation is by means of the arithmetic-geometric mean
 * algorithm, except when m is within 1e-9 of 0 or 1.
 *
 * In the latter case, when m is near 1, the functions `sn`, `cn`, and `dn`
 * are found using an approximation for `phi` given by 16.15.4 in [2].
 * As this approximation is only valid when `phi < pi/2` (i.e.,
 * when `u<K`), the computation only uses 16.15.4 to find the "principal"
 * part of `phi` or the remainder `dphi=phi%(pi/2)`. The returned value
 * value of `phi`, valid for all `u`, is `pi*n/2+dphi` where `n=floor(u/K)`
 * and `K` is the quarter-period.
 *
 *
 * ACCURACY:
 *
 * Tested at random points with u between 0 and 10, m between
 * 0 and 1.
 *
 *            Absolute error (* = relative error):
 * arithmetic   function   # trials      peak         rms
 *    IEEE      phi         10000       9.2e-16*    1.4e-16*
 *    IEEE      sn          50000       4.1e-15     4.6e-16
 *    IEEE      cn          40000       3.6e-15     4.4e-16
 *    IEEE      dn          10000       1.3e-12     1.8e-14
 *
 *  Peak error observed in consistency check using addition
 * theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
 * the above relation to the incomplete elliptic integral.
 * Accuracy deteriorates when u is large.
 *
 */

/*                                                     ellpj.c         */


/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Scipy changes:
 * - 07-18-2016: improve evaluation of dn near quarter periods
 */

#include "mconf.h"
extern double MACHEP;

int ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph)
{
    double ai, b, phi, t, twon, dnfac, K, du, uabs;
    double a[9], c[9];
    int i,j,r;

    /* Check for special cases */
    if (m < 0.0 || m > 1.0 || cephes_isnan(m)) {
        mtherr("ellpj", DOMAIN);
        *sn = NPY_NAN;
        *cn = NPY_NAN;
        *ph = NPY_NAN;
        *dn = NPY_NAN;
        return (-1);
    }
    if (m < 1.0e-9) {
        t = sin(u);
        b = cos(u);
        ai = 0.25 * m * (u - t * b);
        *sn = t - ai * b;
        *cn = b + ai * t;
        *ph = u - ai;
        *dn = 1.0 - 0.5 * m * t * t;
        return (0);
    }

    /* if m is near 1, use approximation for amplitude
       given in Abramowitz and Stegun, 16.15.4 */
    if (m >= 0.9999999999) {
        /* Evaluate the complete elliptic integral using the
           1-m<<1 approximation in ellpk.c */
        K = ellpk(1.0-m);
        uabs = fabs(u);
        j = (int)(uabs/K);
        du = uabs-((double)j)*K;
        r = j%4;

        /* When phi is in the I and III quadrants */
        if(r == 0 || r == 2) {
            t = cosh(du);
            phi = 0.25*(1.0-m)*(sinh(du)*t-du)/t;
            /* add Gudermannian term */
            phi += 2.0*atan(exp(du)) - NPY_PI_2;
        }
        /* When phi is in the II and IV quadrants */
        else{
            du = K-du;
            t = cosh(du);
            phi = 0.25*(1.0-m)*(sinh(du)*t-du)/t;
            /* add Gudermannian term */
            phi += 2.0*atan(exp(du)) - NPY_PI_2;
            phi = NPY_PI_2-phi;
        }
        phi += ((double)j)*NPY_PI_2;
        if (u < 0) {
            phi *= -1.0;
        }

        *sn = sin(phi);
        *cn = cos(phi);
        *dn = sqrt(1.0-m*(*sn)*(*sn));
        *ph = phi;

        return (0);
    }

    /* A. G. M. scale. See DLMF 22.20(ii) */
    a[0] = 1.0;
    b = sqrt(1.0 - m);
    c[0] = sqrt(m);
    twon = 1.0;
    i = 0;

    while (fabs(c[i] / a[i]) > MACHEP) {
        if (i > 7) {
            mtherr("ellpj", OVERFLOW);
            goto done;
        }
        ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        t = sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }

 done:

    /* backward recurrence */
    phi = twon * a[i] * u;
    do {
        t = c[i] * sin(phi) / a[i];
        b = phi;
        phi = (asin(t) + phi) / 2.0;
    }
    while (--i);

    *sn = sin(phi);
    t = cos(phi);
    *cn = t;
    dnfac = cos(phi - b);
    /* See discussion after DLMF 22.20.5 */
    if (fabs(dnfac) < 0.1) {
        *dn = sqrt(1 - m*(*sn)*(*sn));
    }
    else {
        *dn = t / dnfac;
    }
    *ph = phi;
    return (0);
}
