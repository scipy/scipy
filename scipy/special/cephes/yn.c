/*                                                     yn.c
 *
 *     Bessel function of second kind of integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, yn();
 * int n;
 *
 * y = yn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order n, where n is a
 * (possibly negative) integer.
 *
 * The function is evaluated by forward recurrence on
 * n, starting with values computed by the routines
 * y0() and y1().
 *
 * If n = 0 or 1 the routine for y0 or y1 is called
 * directly.
 *
 *
 *
 * ACCURACY:
 *
 *
 *                      Absolute error, except relative
 *                      when y > 1:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       3.4e-15     4.3e-16
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * yn singularity   x = 0              NPY_INFINITY
 * yn overflow                         NPY_INFINITY
 *
 * Spot checked against tables for x, n between 0 and 100.
 *
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */

#include "mconf.h"
extern double MAXLOG;

double yn(n, x)
int n;
double x;
{
    double an, anm1, anm2, r;
    int k, sign;

    if (n < 0) {
	n = -n;
	if ((n & 1) == 0)	/* -1**n */
	    sign = 1;
	else
	    sign = -1;
    }
    else
	sign = 1;


    if (n == 0)
	return (sign * y0(x));
    if (n == 1)
	return (sign * y1(x));

    /* test for overflow */
    if (x == 0.0) {
	sf_error("yn", SF_ERROR_SINGULAR, NULL);
	return -NPY_INFINITY * sign;
    }
    else if (x < 0.0) {
	sf_error("yn", SF_ERROR_DOMAIN, NULL);
	return NPY_NAN;
    }

    /* forward recurrence on n */

    anm2 = y0(x);
    anm1 = y1(x);
    k = 1;
    r = 2 * k;
    do {
	an = r * anm1 / x - anm2;
	anm2 = anm1;
	anm1 = an;
	r += 2.0;
	++k;
    }
    while (k < n);


    return (sign * an);
}
