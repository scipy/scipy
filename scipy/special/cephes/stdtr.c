/*                                                     stdtr.c
 *
 *     Student's t distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double t, stdtr();
 * short k;
 *
 * y = stdtr( k, t );
 *
 *
 * DESCRIPTION:
 *
 * Computes the integral from minus infinity to t of the Student
 * t distribution with integer k > 0 degrees of freedom:
 *
 *                                      t
 *                                      -
 *                                     | |
 *              -                      |         2   -(k+1)/2
 *             | ( (k+1)/2 )           |  (     x   )
 *       ----------------------        |  ( 1 + --- )        dx
 *                     -               |  (      k  )
 *       sqrt( k pi ) | ( k/2 )        |
 *                                   | |
 *                                    -
 *                                   -inf.
 *
 * Relation to incomplete beta integral:
 *
 *        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
 * where
 *        z = k/(k + t**2).
 *
 * For t < -2, this is the method of computation.  For higher t,
 * a direct method is derived from integration by parts.
 * Since the function is symmetric about t=0, the area under the
 * right tail of the density is found by calling the function
 * with -t instead of t.
 *
 * ACCURACY:
 *
 * Tested at random 1 <= k <= 25.  The "domain" refers to t.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -100,-2      50000       5.9e-15     1.4e-15
 *    IEEE     -2,100      500000       2.7e-15     4.9e-17
 */

/*                                                     stdtri.c
 *
 *     Functional inverse of Student's t distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * double p, t, stdtri();
 * int k;
 *
 * t = stdtri( k, p );
 *
 *
 * DESCRIPTION:
 *
 * Given probability p, finds the argument t such that stdtr(k,t)
 * is equal to p.
 *
 * ACCURACY:
 *
 * Tested at random 1 <= k <= 100.  The "domain" refers to p:
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    .001,.999     25000       5.7e-15     8.0e-16
 *    IEEE    10^-6,.001    25000       2.0e-12     2.9e-14
 */


/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1987, 1995 by Stephen L. Moshier
 */

#include "mconf.h"
#include <float.h>

extern double MACHEP;

double stdtr(k, t)
int k;
double t;
{
    double x, rk, z, f, tz, p, xsqk;
    int j;

    if (k <= 0) {
	sf_error("stdtr", SF_ERROR_DOMAIN, NULL);
	return (NPY_NAN);
    }

    if (t == 0)
	return (0.5);

    if (t < -2.0) {
	rk = k;
	z = rk / (rk + t * t);
	p = 0.5 * incbet(0.5 * rk, 0.5, z);
	return (p);
    }

    /*     compute integral from -t to + t */

    if (t < 0)
	x = -t;
    else
	x = t;

    rk = k;			/* degrees of freedom */
    z = 1.0 + (x * x) / rk;

    /* test if k is odd or even */
    if ((k & 1) != 0) {

	/*      computation for odd k   */

	xsqk = x / sqrt(rk);
	p = atan(xsqk);
	if (k > 1) {
	    f = 1.0;
	    tz = 1.0;
	    j = 3;
	    while ((j <= (k - 2)) && ((tz / f) > MACHEP)) {
		tz *= (j - 1) / (z * j);
		f += tz;
		j += 2;
	    }
	    p += f * xsqk / z;
	}
	p *= 2.0 / NPY_PI;
    }


    else {

	/*      computation for even k  */

	f = 1.0;
	tz = 1.0;
	j = 2;

	while ((j <= (k - 2)) && ((tz / f) > MACHEP)) {
	    tz *= (j - 1) / (z * j);
	    f += tz;
	    j += 2;
	}
	p = f * x / sqrt(z * rk);
    }

    /*     common exit     */


    if (t < 0)
	p = -p;			/* note destruction of relative accuracy */

    p = 0.5 + 0.5 * p;
    return (p);
}

double stdtri(k, p)
int k;
double p;
{
    double t, rk, z;
    int rflg;

    if (k <= 0 || p <= 0.0 || p >= 1.0) {
	sf_error("stdtri", SF_ERROR_DOMAIN, NULL);
	return (NPY_NAN);
    }

    rk = k;

    if (p > 0.25 && p < 0.75) {
	if (p == 0.5)
	    return (0.0);
	z = 1.0 - 2.0 * p;
	z = incbi(0.5, 0.5 * rk, fabs(z));
	t = sqrt(rk * z / (1.0 - z));
	if (p < 0.5)
	    t = -t;
	return (t);
    }
    rflg = -1;
    if (p >= 0.5) {
	p = 1.0 - p;
	rflg = 1;
    }
    z = incbi(0.5 * rk, 0.5, 2.0 * p);

    if (DBL_MAX * z < rk)
	return (rflg * NPY_INFINITY);
    t = sqrt(rk / z - rk);
    return (rflg * t);
}
