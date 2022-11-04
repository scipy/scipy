/*                                                     tandg.c
 *
 *     Circular tangent of argument in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, tandg();
 *
 * y = tandg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the argument x in degrees.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     0,10         30000      3.2e-16      8.4e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * tandg total loss   x > 1.0e14 (IEEE)     0.0
 * tandg singularity  x = 180 k  +  90     INFINITY
 */
/*							cotdg.c
 *
 *	Circular cotangent of argument in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cotdg();
 *
 * y = cotdg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular cotangent of the argument x in degrees.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * cotdg total loss   x > 1.0e14 (IEEE)     0.0
 * cotdg singularity  x = 180 k            INFINITY
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1984, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include "mconf.h"

static double PI180 = 1.74532925199432957692E-2;
static double lossth = 1.0e14;

static double tancot(double, int);

double tandg(double x)
{
    return (tancot(x, 0));
}


double cotdg(double x)
{
    return (tancot(x, 1));
}


static double tancot(double xx, int cotflg)
{
    double x;
    int sign;

    /* make argument positive but save the sign */
    if (xx < 0) {
	x = -xx;
	sign = -1;
    }
    else {
	x = xx;
	sign = 1;
    }

    if (x > lossth) {
	sf_error("tandg", SF_ERROR_NO_RESULT, NULL);
	return 0.0;
    }

    /* modulo 180 */
    x = x - 180.0 * floor(x / 180.0);
    if (cotflg) {
	if (x <= 90.0) {
	    x = 90.0 - x;
	}
	else {
	    x = x - 90.0;
	    sign *= -1;
	}
    }
    else {
	if (x > 90.0) {
	    x = 180.0 - x;
	    sign *= -1;
	}
    }
    if (x == 0.0) {
	return 0.0;
    }
    else if (x == 45.0) {
	return sign * 1.0;
    }
    else if (x == 90.0) {
	sf_error((cotflg ? "cotdg" : "tandg"), SF_ERROR_SINGULAR, NULL);
	return INFINITY;
    }
    /* x is now transformed into [0, 90) */
    return sign * tan(x * PI180);
}
