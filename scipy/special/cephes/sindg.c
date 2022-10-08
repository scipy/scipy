/*                                                     sindg.c
 *
 *     Circular sine of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, sindg();
 *
 * y = sindg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of 45 degrees.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *      x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *      1  -  x**2 P(x**2).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE      +-1000       30000      2.3e-16      5.6e-17
 *
 * ERROR MESSAGES:
 *
 *   message           condition        value returned
 * sindg total loss   x > 1.0e14 (IEEE)     0.0
 *
 */
/*							cosdg.c
 *
 *	Circular cosine of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cosdg();
 *
 * y = cosdg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of 45 degrees.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the cosine is approximated by
 *      1  -  x**2 P(x**2).
 * Between pi/4 and pi/2 the sine is represented as
 *      x  +  x**3 P(x**2).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE     +-1000        30000       2.1e-16     5.7e-17
 *  See also sin().
 *
 */

/* Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

#include "mconf.h"

static double sincof[] = {
    1.58962301572218447952E-10,
    -2.50507477628503540135E-8,
    2.75573136213856773549E-6,
    -1.98412698295895384658E-4,
    8.33333333332211858862E-3,
    -1.66666666666666307295E-1
};

static double coscof[] = {
    1.13678171382044553091E-11,
    -2.08758833757683644217E-9,
    2.75573155429816611547E-7,
    -2.48015872936186303776E-5,
    1.38888888888806666760E-3,
    -4.16666666666666348141E-2,
    4.99999999999999999798E-1
};

static double PI180 = 1.74532925199432957692E-2;	/* pi/180 */
static double lossth = 1.0e14;

double sindg(x)
double x;
{
    double y, z, zz;
    int j, sign;

    /* make argument positive but save the sign */
    sign = 1;
    if (x < 0) {
	x = -x;
	sign = -1;
    }

    if (x > lossth) {
	sf_error("sindg", SF_ERROR_NO_RESULT, NULL);
	return (0.0);
    }

    y = floor(x / 45.0);	/* integer part of x/NPY_PI_4 */

    /* strip high bits of integer part to prevent integer overflow */
    z = ldexp(y, -4);
    z = floor(z);		/* integer part of y/8 */
    z = y - ldexp(z, 4);	/* y - 16 * (y/16) */

    j = z;			/* convert to integer for tests on the phase angle */
    /* map zeros to origin */
    if (j & 1) {
	j += 1;
	y += 1.0;
    }
    j = j & 07;			/* octant modulo 360 degrees */
    /* reflect in x axis */
    if (j > 3) {
	sign = -sign;
	j -= 4;
    }

    z = x - y * 45.0;		/* x mod 45 degrees */
    z *= PI180;			/* multiply by pi/180 to convert to radians */
    zz = z * z;

    if ((j == 1) || (j == 2)) {
	y = 1.0 - zz * polevl(zz, coscof, 6);
    }
    else {
	y = z + z * (zz * polevl(zz, sincof, 5));
    }

    if (sign < 0)
	y = -y;

    return (y);
}


double cosdg(x)
double x;
{
    double y, z, zz;
    int j, sign;

    /* make argument positive */
    sign = 1;
    if (x < 0)
	x = -x;

    if (x > lossth) {
	sf_error("cosdg", SF_ERROR_NO_RESULT, NULL);
	return (0.0);
    }

    y = floor(x / 45.0);
    z = ldexp(y, -4);
    z = floor(z);		/* integer part of y/8 */
    z = y - ldexp(z, 4);	/* y - 16 * (y/16) */

    /* integer and fractional part modulo one octant */
    j = z;
    if (j & 1) {		/* map zeros to origin */
	j += 1;
	y += 1.0;
    }
    j = j & 07;
    if (j > 3) {
	j -= 4;
	sign = -sign;
    }

    if (j > 1)
	sign = -sign;

    z = x - y * 45.0;		/* x mod 45 degrees */
    z *= PI180;			/* multiply by pi/180 to convert to radians */

    zz = z * z;

    if ((j == 1) || (j == 2)) {
	y = z + z * (zz * polevl(zz, sincof, 5));
    }
    else {
	y = 1.0 - zz * polevl(zz, coscof, 6);
    }

    if (sign < 0)
	y = -y;

    return (y);
}


/* Degrees, minutes, seconds to radians: */

/* 1 arc second, in radians = 4.848136811095359935899141023579479759563533023727e-6 */
static double P64800 =
    4.848136811095359935899141023579479759563533023727e-6;

double radian(d, m, s)
double d, m, s;
{
    return (((d * 60.0 + m) * 60.0 + s) * P64800);
}
