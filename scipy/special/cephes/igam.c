/*                                                     igam.c
 *
 *     Incomplete Gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igam();
 *
 * y = igam( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30       200000       3.6e-14     2.9e-15
 *    IEEE      0,100      300000       9.9e-14     1.5e-14
 */
/*							igamc()
 *
 *	Complemented incomplete Gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igamc();
 *
 * y = igamc( a, x );
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 * ACCURACY:
 *
 * Tested at random a, x.
 *                a         x                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
 *    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/* Scipy changes:
 * - 05-01-2016: added asymptotic expansion for igam to improve the
 *   a ~ x regime.
 */

#include "mconf.h"
#include "igam.h"

#ifdef MAXITER
#undef MAXITER
#endif
#define MAXITER 1000

#define SMALL 25
#define TOL 2.2204460492503131e-16

extern double MACHEP, MAXLOG;
static double big = 4.503599627370496e15;
static double biginv = 2.22044604925031308085e-16;

double igamc(double, double);
double igam(double, double);
double igam_pow(double, double);
double igam_asy(double, double);


double igamc(a, x)
double a, x;
{
    int i;
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;

    if ((x < 0) || (a <= 0)) {
	mtherr("gammaincc", DOMAIN);
	return (NPY_NAN);
    }

    if ((x < 1.0) || (x < a))
	return (1.0 - igam(a, x));
    
    if (cephes_isinf(x))
    return 0.0;

    ax = a * log(x) - x - lgam(a);
    if (ax < -MAXLOG) {
	mtherr("igamc", UNDERFLOW);
	return (0.0);
    }
    ax = exp(ax);

    /* continued fraction */
    y = 1.0 - a;
    z = x + y + 1.0;
    c = 0.0;
    pkm2 = 1.0;
    qkm2 = x;
    pkm1 = x + 1.0;
    qkm1 = z * x;
    ans = pkm1 / qkm1;

    for (i = 0; i < MAXITER; i++) {
	c += 1.0;
	y += 1.0;
	z += 2.0;
	yc = y * c;
	pk = pkm1 * z - pkm2 * yc;
	qk = qkm1 * z - qkm2 * yc;
	if (qk != 0) {
	    r = pk / qk;
	    t = fabs((ans - r) / r);
	    ans = r;
	}
	else
	    t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if (fabs(pk) > big) {
	    pkm2 *= biginv;
	    pkm1 *= biginv;
	    qkm2 *= biginv;
	    qkm1 *= biginv;
	}
	if (t <= MACHEP) {
	    break;
	}
    }

    return (ans * ax);
}


double igam(a, x)
double a, x;
{
    double lambda;

    /* Check zero integration limit first */
    if (x == 0)
	return (0.0);

    if ((x < 0) || (a <= 0)) {
	mtherr("gammainc", DOMAIN);
	return (NPY_NAN);
    }

    lambda = x / a;
    if (x > SMALL && a > SMALL && lambda > 0.7 && lambda < 1.3) {
	return igam_asy(a, x);
    }

    if ((x > 1.0) && (x > a))
	return (1.0 - igamc(a, x));

    return igam_pow(a, x);
}


/* Compute igam using DLMF 8.11.4. */
double igam_pow(a, x)
double a, x;
{
    int i;
    double ans, ax, c, r;

    /* Compute  x**a * exp(-x) / Gamma(a)  */
    ax = a * log(x) - x - lgam(a);
    if (ax < -MAXLOG) {
	mtherr("igam", UNDERFLOW);
	return (0.0);
    }
    ax = exp(ax);
    
    /* power series */
    r = a;
    c = 1.0;
    ans = 1.0;
    
    for (i = 0; i < MAXITER; i++) {
	r += 1.0;
	c *= x / r;
	ans += c;
	if (c <= MACHEP * ans) {
	    break;
	}
    }
    
    return (ans * ax / a);
}


/* Compute igam using DLMF 8.12.3. */
double igam_asy(a, x)
double a, x;
{
    int k, n;
    int maxpow = 0;
    double lambda = x/a;
    double eta, res, ck, ckterm, term, absterm;
    double absoldterm = NPY_INFINITY;
    double etapow[N] = {1};
    double sum = 0;
    double afac = 1;

    if (lambda > 1) {
	eta = sqrt(2*(lambda - 1 - log(lambda)));
    } else if (lambda < 1) {
	eta = -sqrt(2*(lambda - 1 - log(lambda)));
    } else {
	eta = 0;
    }
    res = 0.5*erfc(-eta*sqrt(a/2));

    for (k = 0; k < K; k++) {
	ck = d[k][0];
	for (n = 1; n < N; n++) {
	    if (n > maxpow) {
		etapow[n] = eta*etapow[n-1];
		maxpow += 1;
	    }
	    ckterm = d[k][n]*etapow[n];
	    ck += ckterm;
	    if (fabs(ckterm) < TOL*fabs(ck)) {
		break;
	    }
	}
	term = ck*afac;
	absterm = fabs(term);
	if (absterm > absoldterm) {
	    break;
	}
	sum += term;
	if (absterm < TOL*fabs(sum)) {
	    break;
	}
	absoldterm = absterm;
	afac /= a;
    }
    res -= exp(-0.5*a*eta*eta)*sum/sqrt(2*NPY_PI*a);
	    
    return res;
}
