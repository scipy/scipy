/*                                                     expn.c
 *
 *             Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, expn();
 *
 * y = expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/*                                                     expn.c  */

/* Cephes Math Library Release 1.1:  March, 1985
 * Copyright 1985 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

/* Sources
 * [1] NIST, "The Digital Library of Mathematical Functions", dlmf.nist.gov
 */

/* Scipy changes:
 * - 09-10-2016: improved asymptotic expansion for large n
 */

#include "mconf.h"
#include "polevl.h"
#include "expn.h"

#define EUL 0.57721566490153286060
#define BIG  1.44115188075855872E+17
extern double MACHEP, MAXLOG;

static double expn_large_n(int, double);


double expn(int n, double x)
{
    double ans, r, t, yk, xk;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;
    double psi, z;
    int i, k;
    static double big = BIG;

    if (isnan(x)) {
	return NAN;
    }
    else if (n < 0 || x < 0) {
	sf_error("expn", SF_ERROR_DOMAIN, NULL);
	return NAN;
    }

    if (x > MAXLOG) {
	return (0.0);
    }

    if (x == 0.0) {
	if (n < 2) {
	    sf_error("expn", SF_ERROR_SINGULAR, NULL);
	    return (INFINITY);
	}
	else {
	    return (1.0 / (n - 1.0));
	}
    }

    if (n == 0) {
	return (exp(-x) / x);
    }

    /* Asymptotic expansion for large n, DLMF 8.20(ii) */
    if (n > 50) {
	ans = expn_large_n(n, x);
	goto done;
    }

    if (x > 1.0) {
	goto cfrac;
    }

    /* Power series expansion, DLMF 8.19.8 */
    psi = -EUL - log(x);
    for (i = 1; i < n; i++) {
	psi = psi + 1.0 / i;
    }

    z = -x;
    xk = 0.0;
    yk = 1.0;
    pk = 1.0 - n;
    if (n == 1) {
	ans = 0.0;
    } else {
	ans = 1.0 / pk;
    }
    do {
	xk += 1.0;
	yk *= z / xk;
	pk += 1.0;
	if (pk != 0.0) {
	    ans += yk / pk;
	}
	if (ans != 0.0)
	    t = fabs(yk / ans);
	else
	    t = 1.0;
    } while (t > MACHEP);
    k = xk;
    t = n;
    r = n - 1;
    ans = (pow(z, r) * psi / Gamma(t)) - ans;
    goto done;

    /* Continued fraction, DLMF 8.19.17 */
  cfrac:
    k = 1;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = 1.0;
    qkm1 = x + n;
    ans = pkm1 / qkm1;

    do {
	k += 1;
	if (k & 1) {
	    yk = 1.0;
	    xk = n + (k - 1) / 2;
	} else {
	    yk = x;
	    xk = k / 2;
	}
	pk = pkm1 * yk + pkm2 * xk;
	qk = qkm1 * yk + qkm2 * xk;
	if (qk != 0) {
	    r = pk / qk;
	    t = fabs((ans - r) / r);
	    ans = r;
	} else {
	    t = 1.0;
	}
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if (fabs(pk) > big) {
	    pkm2 /= big;
	    pkm1 /= big;
	    qkm2 /= big;
	    qkm1 /= big;
	}
    } while (t > MACHEP);

    ans *= exp(-x);

  done:
    return (ans);
}


/* Asymptotic expansion for large n, DLMF 8.20(ii) */
static double expn_large_n(int n, double x)
{
    int k;
    double p = n;
    double lambda = x/p;
    double multiplier = 1/p/(lambda + 1)/(lambda + 1);
    double fac = 1;
    double res = 1; /* A[0] = 1 */
    double expfac, term;

    expfac = exp(-lambda*p)/(lambda + 1)/p;
    if (expfac == 0) {
	sf_error("expn", SF_ERROR_UNDERFLOW, NULL);
	return 0;
    }

    /* Do the k = 1 term outside the loop since A[1] = 1 */
    fac *= multiplier;
    res += fac;

    for (k = 2; k < nA; k++) {
	fac *= multiplier;
	term = fac*polevl(lambda, A[k], Adegs[k]);
	res += term;
	if (fabs(term) < MACHEP*fabs(res)) {
	    break;
	}
    }

    return expfac*res;
}
