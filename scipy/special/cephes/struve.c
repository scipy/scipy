/*
 * Compute the Struve function.
 *
 * Notes
 * -----
 *
 * We use three expansions for the Struve function discussed in [1]:
 *
 * - power series
 * - expansion in Bessel functions
 * - asymptotic large-z expansion
 *
 * Rounding errors are estimated based on the largest terms in the sums.
 *
 * ``struve_convergence.py`` plots the convergence regions of the different
 * expansions.
 *
 * (i)
 *
 * Looking at the error in the asymptotic expansion, one finds that
 * it's not worth trying if z ~> 0.7 * v + 12 for v > 0.
 *
 * (ii)
 *
 * The Bessel function expansion tends to fail for |z| >~ |v| and is not tried
 * there.
 *
 * For Struve H it covers the quadrant v > z where the power series may fail to
 * produce reasonable results.
 *
 * (iii)
 *
 * The three expansions together cover for Struve H the region z > 0, v real.
 *
 * They also cover Struve L, except that some loss of precision may occur around
 * the transition region z ~ 0.7 |v|, v < 0, |v| >> 1 where the function changes
 * rapidly.
 *
 * (iv)
 *
 * The power series is evaluated in double-double precision. This fixes accuracy
 * issues in Struve H for |v| << |z| before the asymptotic expansion kicks in.
 * Moreover, it improves the Struve L behavior for negative v.
 *
 *
 * References
 * ----------
 * [1] NIST Digital Library of Mathematical Functions
 *     https://dlmf.nist.gov/11
 */

/*
 * Copyright (C) 2013  Pauli Virtanen
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Neither the name of Enthought nor the names of the SciPy Developers
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
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

#include "mconf.h"
#include "dd_real.h"

#include "special_wrappers.h"

#define STRUVE_MAXITER 10000
#define SUM_EPS 1e-16   /* be sure we are in the tail of the sum */
#define SUM_TINY 1e-100
#define GOOD_EPS 1e-12
#define ACCEPTABLE_EPS 1e-7
#define ACCEPTABLE_ATOL 1e-300

#define MIN(a, b) ((a) < (b) ? (a) : (b))

double struve_power_series(double v, double x, int is_h, double *err);
double struve_asymp_large_z(double v, double z, int is_h, double *err);
double struve_bessel_series(double v, double z, int is_h, double *err);

static double bessel_y(double v, double x);
static double bessel_j(double v, double x);
static double struve_hl(double v, double x, int is_h);

double struve_h(double v, double z)
{
    return struve_hl(v, z, 1);
}

double struve_l(double v, double z)
{
    return struve_hl(v, z, 0);
}

static double struve_hl(double v, double z, int is_h)
{
    double value[4], err[4], tmp;
    int n;

    if (z < 0) {
        n = v;
        if (v == n) {
            tmp = (n % 2 == 0) ? -1 : 1;
            return tmp * struve_hl(v, -z, is_h);
        }
        else {
            return NAN;
        }
    }
    else if (z == 0) {
        if (v < -1) {
            return gammasgn(v + 1.5) * INFINITY;
        }
        else if (v == -1) {
            return 2 / sqrt(M_PI) / Gamma(0.5);
        }
        else {
            return 0;
        }
    }

    n = -v - 0.5;
    if (n == -v - 0.5 && n > 0) {
        if (is_h) {
            return (n % 2 == 0 ? 1 : -1) * bessel_j(n + 0.5, z);
        }
        else {
            return iv(n + 0.5, z);
        }
    }

    /* Try the asymptotic expansion */
    if (z >= 0.7*v + 12) {
        value[0] = struve_asymp_large_z(v, z, is_h, &err[0]);
        if (err[0] < GOOD_EPS * fabs(value[0])) {
            return value[0];
        }
    }
    else {
        err[0] = INFINITY;
    }

    /* Try power series */
    value[1] = struve_power_series(v, z, is_h, &err[1]);
    if (err[1] < GOOD_EPS * fabs(value[1])) {
        return value[1];
    }

    /* Try bessel series */
    if (fabs(z) < fabs(v) + 20) {
        value[2] = struve_bessel_series(v, z, is_h, &err[2]);
        if (err[2] < GOOD_EPS * fabs(value[2])) {
            return value[2];
        }
    }
    else {
        err[2] = INFINITY;
    }

    /* Return the best of the three, if it is acceptable */
    n = 0;
    if (err[1] < err[n]) n = 1;
    if (err[2] < err[n]) n = 2;
    if (err[n] < ACCEPTABLE_EPS * fabs(value[n]) || err[n] < ACCEPTABLE_ATOL) {
        return value[n];
    }

    /* Maybe it really is an overflow? */
    tmp = -lgam(v + 1.5) + (v + 1)*log(z/2);
    if (!is_h) {
        tmp = fabs(tmp);
    }
    if (tmp > 700) {
        sf_error("struve", SF_ERROR_OVERFLOW, NULL);
        return INFINITY * gammasgn(v + 1.5);
    }

    /* Failure */
    sf_error("struve", SF_ERROR_NO_RESULT, NULL);
    return NAN;
}


/*
 * Power series for Struve H and L
 * https://dlmf.nist.gov/11.2.1
 *
 * Starts to converge roughly at |n| > |z|
 */
double struve_power_series(double v, double z, int is_h, double *err)
{
    int n, sgn;
    double term, sum, maxterm, scaleexp, tmp;
    double2 cterm, csum, cdiv, z2, c2v, ctmp;

    if (is_h) {
        sgn = -1;
    }
    else {
        sgn = 1;
    }

    tmp = -lgam(v + 1.5) + (v + 1)*log(z/2);
    if (tmp < -600 || tmp > 600) {
        /* Scale exponent to postpone underflow/overflow */
        scaleexp = tmp/2;
        tmp -= scaleexp;
    }
    else {
        scaleexp = 0;
    }

    term = 2 / sqrt(M_PI) * exp(tmp) * gammasgn(v + 1.5);
    sum = term;
    maxterm = 0;

    cterm = dd_create_d(term);
    csum = dd_create_d(sum);
    z2 = dd_create_d(sgn*z*z);
    c2v = dd_create_d(2*v);

    for (n = 0; n < STRUVE_MAXITER; ++n) {
        /* cdiv = (3 + 2*n) * (3 + 2*n + 2*v)) */
        cdiv = dd_create_d(3 + 2*n);
        ctmp = dd_create_d(3 + 2*n);
        ctmp = dd_add(ctmp, c2v);
        cdiv = dd_mul(cdiv, ctmp);

        /* cterm *= z2 / cdiv */
        cterm = dd_mul(cterm, z2);
        cterm = dd_div(cterm, cdiv);

        csum = dd_add(csum, cterm);

        term = dd_to_double(cterm);
        sum = dd_to_double(csum);

        if (fabs(term) > maxterm) {
            maxterm = fabs(term);
        }
        if (fabs(term) < SUM_TINY * fabs(sum) || term == 0 || !isfinite(sum)) {
            break;
        }
    }

    *err = fabs(term) + fabs(maxterm) * 1e-22;

    if (scaleexp != 0) {
        sum *= exp(scaleexp);
        *err *= exp(scaleexp);
    }

    if (sum == 0 && term == 0 && v < 0 && !is_h) {
        /* Spurious underflow */
        *err = INFINITY;
        return NAN;
    }

    return sum;
}


/*
 * Bessel series
 * https://dlmf.nist.gov/11.4.19
 */
double struve_bessel_series(double v, double z, int is_h, double *err)
{
    int n;
    double term, cterm, sum, maxterm;

    if (is_h && v < 0) {
        /* Works less reliably in this region */
        *err = INFINITY;
        return NAN;
    }

    sum = 0;
    maxterm = 0;

    cterm = sqrt(z / (2*M_PI));

    for (n = 0; n < STRUVE_MAXITER; ++n) {
        if (is_h) {
            term = cterm * bessel_j(n + v + 0.5, z) / (n + 0.5);
            cterm *= z/2 / (n + 1);
        }
        else {
            term = cterm * iv(n + v + 0.5, z) / (n + 0.5);
            cterm *= -z/2 / (n + 1);
        }
        sum += term;
        if (fabs(term) > maxterm) {
            maxterm = fabs(term);
        }
        if (fabs(term) < SUM_EPS * fabs(sum) || term == 0 || !isfinite(sum)) {
            break;
        }
    }

    *err = fabs(term) + fabs(maxterm) * 1e-16;

    /* Account for potential underflow of the Bessel functions */
    *err += 1e-300 * fabs(cterm);

    return sum;
}


/*
 * Large-z expansion for Struve H and L
 * https://dlmf.nist.gov/11.6.1
 */
double struve_asymp_large_z(double v, double z, int is_h, double *err)
{
    int n, sgn, maxiter;
    double term, sum, maxterm;
    double m;

    if (is_h) {
        sgn = -1;
    }
    else {
        sgn = 1;
    }

    /* Asymptotic expansion divergenge point */
    m = z/2;
    if (m <= 0) {
        maxiter = 0;
    }
    else if (m > STRUVE_MAXITER) {
        maxiter = STRUVE_MAXITER;
    }
    else {
        maxiter = (int)m;
    }
    if (maxiter == 0) {
        *err = INFINITY;
        return NAN;
    }

    if (z < v) {
        /* Exclude regions where our error estimation fails */
        *err = INFINITY;
        return NAN;
    }

    /* Evaluate sum */
    term = -sgn / sqrt(M_PI) * exp(-lgam(v + 0.5) + (v - 1) * log(z/2)) * gammasgn(v + 0.5);
    sum = term;
    maxterm = 0;

    for (n = 0; n < maxiter; ++n) {
        term *= sgn * (1 + 2*n) * (1 + 2*n - 2*v) / (z*z);
        sum += term;
        if (fabs(term) > maxterm) {
            maxterm = fabs(term);
        }
        if (fabs(term) < SUM_EPS * fabs(sum) || term == 0 || !isfinite(sum)) {
            break;
        }
    }

    if (is_h) {
        sum += bessel_y(v, z);
    }
    else {
        sum += iv(v, z);
    }

    /*
     * This error estimate is strictly speaking valid only for
     * n > v - 0.5, but numerical results indicate that it works
     * reasonably.
     */
    *err = fabs(term) + fabs(maxterm) * 1e-16;

    return sum;
}


static double bessel_y(double v, double x)
{
    return cbesy_wrap_real(v, x);
}

static double bessel_j(double v, double x)
{
    return cbesj_wrap_real(v, x);
}
