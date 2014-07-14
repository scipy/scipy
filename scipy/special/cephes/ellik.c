/*                                                     ellik.c
 *
 *     Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellik();
 *
 * y = ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi | m) =   |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */


/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

/*     Incomplete elliptic integral of first kind      */

#include "mconf.h"
extern double MACHEP;

double ellik(phi, m)
double phi, m;
{
    double a, b, c, e, temp, t, K, denom;
    int d, mod, sign, npio2;

    if (m == 0.0)
	return (phi);
    a = 1.0 - m;
    if (a == 0.0) {
	if (fabs(phi) >= NPY_PI_2) {
	    mtherr("ellik", SING);
	    return (NPY_INFINITY);
	}
	return (log(tan((NPY_PI_2 + phi) / 2.0)));
    }
    npio2 = floor(phi / NPY_PI_2);
    if (npio2 & 1)
	npio2 += 1;
    if (npio2) {
	K = ellpk(a);
	phi = phi - npio2 * NPY_PI_2;
    }
    else
	K = 0.0;
    if (phi < 0.0) {
	phi = -phi;
	sign = -1;
    }
    else
	sign = 0;
    b = sqrt(a);
    t = tan(phi);
    if (fabs(t) > 10.0) {
	/* Transform the amplitude */
	e = 1.0 / (b * t);
	/* ... but avoid multiple recursions.  */
	if (fabs(e) < 10.0) {
	    e = atan(e);
	    if (npio2 == 0)
		K = ellpk(a);
	    temp = K - ellik(e, m);
	    goto done;
	}
    }
    a = 1.0;
    c = sqrt(m);
    d = 1;
    mod = 0;

    while (fabs(c / a) > MACHEP) {
	temp = b / a;
	phi = phi + atan(t * temp) + mod * NPY_PI;
        denom = 1.0 - temp * t * t;
        if (fabs(denom) > 10*MACHEP) {
	    t = t * (1.0 + temp) / denom;
            mod = (phi + NPY_PI_2) / NPY_PI;
        }
        else {
            t = tan(phi);
            mod = (int)floor((phi - atan(t))/NPY_PI);
        }
	c = (a - b) / 2.0;
	temp = sqrt(a * b);
	a = (a + b) / 2.0;
	b = temp;
	d += d;
    }

    temp = (atan(t) + mod * NPY_PI) / (d * a);

  done:
    if (sign < 0)
	temp = -temp;
    temp += npio2 * K;
    return (temp);
}
