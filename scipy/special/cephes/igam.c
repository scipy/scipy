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

/* Sources
 * [1] "The Digital Library of Mathematical Functions", dlmf.nist.gov
 * [2] Maddock et. al., "Incomplete Gamma Functions",
 *     https://www.boost.org/doc/libs/1_61_0/libs/math/doc/html/math_toolkit/sf_gamma/igamma.html
 */

/* Scipy changes:
 * - 05-01-2016: added asymptotic expansion for igam to improve the
 *   a ~ x regime.
 * - 06-19-2016: additional series expansion added for igamc to
 *   improve accuracy at small arguments.
 * - 06-24-2016: better choice of domain for the asymptotic series;
 *   improvements in accuracy for the asymptotic series when a and x
 *   are very close.
 */

#include "mconf.h"
#include "lanczos.h"
#include "igam.h"

#ifdef MAXITER
#undef MAXITER
#endif

#define MAXITER 2000
#define IGAM 1
#define IGAMC 0
#define SMALL 20
#define LARGE 200
#define SMALLRATIO 0.3
#define LARGERATIO 4.5

extern double MACHEP, MAXLOG;
static double big = 4.503599627370496e15;
static double biginv = 2.22044604925031308085e-16;

static double igamc_continued_fraction(double, double);
static double igam_series(double, double);
static double igamc_series(double, double);
static double asymptotic_series(double, double, int);


double igam(a, x)
double a, x;
{
    double absxma_a;

    /* Check zero integration limit first */
    if (x == 0)
	return (0.0);

    if ((x < 0) || (a <= 0)) {
	mtherr("gammainc", DOMAIN);
	return (NPY_NAN);
    }

    /* Asymptotic regime where a ~ x; see [2]. */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
	return asymptotic_series(a, x, IGAM);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
	return asymptotic_series(a, x, IGAM);
    }

    if ((x > 1.0) && (x > a)) {
	return (1.0 - igamc(a, x));
    }

    return igam_series(a, x);
}


double igamc(double a, double x)
{
    double absxma_a;

    if ((x < 0) || (a <= 0)) {
	mtherr("gammaincc", DOMAIN);
	return (NPY_NAN);
    } else if (x == 0) {
	return 1;
    } else if (cephes_isinf(x)) {
	return 0.0;
    }

    /* Asymptotic regime where a ~ x; see [2]. */
    absxma_a = fabs(x - a) / a;
    if ((a > SMALL) && (a < LARGE) && (absxma_a < SMALLRATIO)) {
	return asymptotic_series(a, x, IGAMC);
    } else if ((a > LARGE) && (absxma_a < LARGERATIO / sqrt(a))) {
	return asymptotic_series(a, x, IGAMC);
    }
    
    /* Everywhere else; see [2]. */
    if (x > 1.1) {
	if (x < a) {
	    return 1.0 - igam_series(a, x);
	} else {
	    return igamc_continued_fraction(a, x);
	}
    } else if (x <= 0.5) {
	if (-0.4 / log(x) < a) {
	    return 1.0 - igam_series(a, x);
	} else {
	    return igamc_series(a, x);
	}
    } else {
	if (x * 1.1 < a) {
	    return 1.0 - igam_series(a, x);
	} else {
	    return igamc_series(a, x);
	}
    }
}


/* Compute
 *
 * x^a * exp(-x) / gamma(a)
 *
 * corrected from (15) and (16) in [2] by replacing exp(x - a) with
 * exp(a - x).
 */
double igam_fac(double a, double x)
{
    double ax, fac, res, num;

    if (fabs(a - x) > 0.4 * fabs(a)) {
	ax = a * log(x) - x - lgam(a);
	if (ax < -MAXLOG) {
	    mtherr("igam", UNDERFLOW);
	    return 0.0;
	}
	return exp(ax);
    }
    
    fac = a + lanczos_g - 0.5;
    res = sqrt(fac / exp(1)) / lanczos_sum_expg_scaled(a);

    if ((a < 200) && (x < 200)) {
	res *= exp(a - x) * pow(x / fac, a);
    } else {
	num = x - a - lanczos_g + 0.5;
	res *= exp(a * log1pmx(num / fac) + x * (0.5 - lanczos_g) / fac);
    }
    
    return res;
}


/* Compute igamc using DLMF 8.9.2. */
static double igamc_continued_fraction(double a, double x)
{
    int i;
    double ans, ax, c, yc, r, t, y, z;
    double pk, pkm1, pkm2, qk, qkm1, qkm2;
    
    ax = igam_fac(a, x);
    if (ax == 0.0) {
	return 0.0;
    }

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


/* Compute igam using DLMF 8.11.4. */
static double igam_series(double a, double x)
{
    int i;
    double ans, ax, c, r;

    ax = igam_fac(a, x);
    if (ax == 0.0) {
	return 0.0;
    }
    
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


/* Compute igamc using DLMF 8.7.3. This is related to the series in
 * igam_series but extra care is taken to avoid cancellation.
 */
static double igamc_series(double a, double x)
{
    int n;
    double fac = 1;
    double sum = 0;
    double term, logx;

    for (n = 1; n < MAXITER; n++) {
	fac *= -x / n;
	term = fac / (a + n);
	sum += term;
	if (fabs(term) <= MACHEP * fabs(sum)) {
	    break;
	}
    }

    logx = log(x);
    term = -expm1(a * logx - lgam1p(a));
    return term - exp(a * logx - lgam(a)) * sum;
}


/* Compute igam/igamc using DLMF 8.12.3/8.12.4. */
static double asymptotic_series(double a, double x, int func)
{
    int k, n, sgn;
    int maxpow = 0;
    double lambda = x / a;
    double sigma = (x - a) / a;
    double eta, res, ck, ckterm, term, absterm;
    double absoldterm = NPY_INFINITY;
    double etapow[N] = {1};
    double sum = 0;
    double afac = 1;

    if (func == IGAM) {
	sgn = -1;
    } else {
	sgn = 1;
    }

    if (lambda > 1) {
 	eta = sqrt(-2 * log1pmx(sigma));
    } else if (lambda < 1) {
	eta = -sqrt(-2 * log1pmx(sigma));
    } else {
	eta = 0;
    }
    res = 0.5 * erfc(sgn * eta * sqrt(a / 2));

    for (k = 0; k < K; k++) {
	ck = d[k][0];
	for (n = 1; n < N; n++) {
	    if (n > maxpow) {
		etapow[n] = eta * etapow[n-1];
		maxpow += 1;
	    }
	    ckterm = d[k][n]*etapow[n];
	    ck += ckterm;
	    if (fabs(ckterm) < MACHEP * fabs(ck)) {
		break;
	    }
	}
	term = ck * afac;
	absterm = fabs(term);
	if (absterm > absoldterm) {
	    break;
	}
	sum += term;
	if (absterm < MACHEP * fabs(sum)) {
	    break;
	}
	absoldterm = absterm;
	afac /= a;
    }
    res += sgn * exp(-0.5 * a * eta * eta) * sum / sqrt(2 * NPY_PI * a);
	    
    return res;
}
