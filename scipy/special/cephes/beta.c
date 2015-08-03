/*                                                     beta.c
 *
 *     Beta function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, y, beta();
 *
 * y = beta( a, b );
 *
 *
 *
 * DESCRIPTION:
 *
 *                   -     -
 *                  | (a) | (b)
 * beta( a, b )  =  -----------.
 *                     -
 *                    | (a+b)
 *
 * For large arguments the logarithm of the function is
 * evaluated using lgam(), then exponentiated.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,30       30000       8.1e-14     1.1e-14
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * beta overflow    log(beta) > MAXLOG       0.0
 *                  a or b <0 integer        0.0
 *
 */

/*                                                     beta.c  */


/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include "mconf.h"

#define MAXGAM 171.624376956302725

extern double MAXLOG;
extern int sgngam;

#define ASYMP_FACTOR 1e6

static double lbeta_asymp(double a, double b, int *sgn);
static double lbeta_negint(int a, double b);
static double beta_negint(int a, double b);
double gammasgn(double x);

double beta(a, b)
double a, b;
{
    double y;
    int sign;

    sign = 1;

    if (a <= 0.0) {
	if (a == floor(a)) {
            if (a == (int)a) {
                return beta_negint((int)a, b);
            }
            else {
                goto over;
            }
        }
    }

    if (b <= 0.0) {
	if (b == floor(b)) {
            if (b == (int)b) {
                return beta_negint((int)b, a);
            }
            else {
                goto over;
            }
        }
    }

    if (fabs(a) < fabs(b)) {
        y = a; a = b; b = y;
    }

    if (fabs(a) > ASYMP_FACTOR * fabs(b) && a > ASYMP_FACTOR) {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        y = lbeta_asymp(a, b, &sign);
        return sign * exp(y);
    }

    y = a + b;
    if (fabs(y) > MAXGAM || fabs(a) > MAXGAM || fabs(b) > MAXGAM) {
	y = lgam(y);
	sign *= sgngam;		/* keep track of the sign */
	y = lgam(b) - y;
	sign *= sgngam;
	y = lgam(a) + y;
	sign *= sgngam;
	if (y > MAXLOG) {
	  over:
	    mtherr("beta", OVERFLOW);
	    return (sign * NPY_INFINITY);
	}
	return (sign * exp(y));
    }

    y = Gamma(y);
    a = Gamma(a);
    b = Gamma(b);
    if (y == 0.0)
	goto over;

    if (fabs(fabs(a) - fabs(y)) > fabs(fabs(b) - fabs(y))) {
        y = b / y;
        y *= a;
    }
    else {
        y = a / y;
        y *= b;
    }

    return (y);
}


/* Natural log of |beta|.  Return the sign of beta in sgngam.  */

double lbeta(a, b)
double a, b;
{
    double y;
    int sign;

    sign = 1;

    if (a <= 0.0) {
	if (a == floor(a)) {
            if (a == (int)a) {
                return lbeta_negint((int)a, b);
            }
            else {
                goto over;
            }
        }
    }

    if (b <= 0.0) {
	if (b == floor(b)) {
            if (b == (int)b) {
                return lbeta_negint((int)b, a);
            }
            else {
                goto over;
            }
        }
    }

    if (fabs(a) < fabs(b)) {
        y = a; a = b; b = y;
    }

    if (fabs(a) > ASYMP_FACTOR * fabs(b) && a > ASYMP_FACTOR) {
        /* Avoid loss of precision in lgam(a + b) - lgam(a) */
        y = lbeta_asymp(a, b, &sign);
        sgngam = sign;
        return y;
    }

    y = a + b;
    if (fabs(y) > MAXGAM || fabs(a) > MAXGAM || fabs(b) > MAXGAM) {
	y = lgam(y);
	sign *= sgngam;		/* keep track of the sign */
	y = lgam(b) - y;
	sign *= sgngam;
	y = lgam(a) + y;
	sign *= sgngam;
	sgngam = sign;
	return (y);
    }

    y = Gamma(y);
    a = Gamma(a);
    b = Gamma(b);
    if (y == 0.0) {
      over:
	mtherr("lbeta", OVERFLOW);
	return (sign * NPY_INFINITY);
    }

    if (fabs(fabs(a) - fabs(y)) > fabs(fabs(b) - fabs(y))) {
        y = b / y;
        y *= a;
    }
    else {
        y = a / y;
        y *= b;
    }

    if (y < 0) {
	sgngam = -1;
	y = -y;
    }
    else
	sgngam = 1;

    return (log(y));
}

/*
 * Asymptotic expansion for  ln(|B(a, b)|) for a > ASYMP_FACTOR*max(|b|, 1).
 */
static double lbeta_asymp(double a, double b, int *sgn)
{
    double r = lgam(b);
    *sgn = sgngam;
    r -= b * log(a);

    r += b*(1-b)/(2*a);
    r += b*(1-b)*(1-2*b)/(12*a*a);
    r += - b*b*(1-b)*(1-b)/(12*a*a*a);

    return r;
}


/*
 * Special case for a negative integer argument
 */

static double beta_negint(int a, double b)
{
    int sgn;
    if (b == (int)b && 1 - a - b > 0) {
        sgn = ((int)b % 2 == 0) ? 1 : -1;
        return sgn * beta(1 - a - b, b);
    }
    else {
	mtherr("lbeta", OVERFLOW);
        return NPY_INFINITY;
    }
}

static double lbeta_negint(int a, double b)
{
    double r;
    int sgn;
    if (b == (int)b && 1 - a - b > 0) {
        sgn = ((int)b % 2 == 0) ? 1 : -1;
        r = lbeta(1 - a - b, b);
        sgngam *= sgn;
        return r;
    }
    else {
	mtherr("lbeta", OVERFLOW);
        return NPY_INFINITY;
    }
}
