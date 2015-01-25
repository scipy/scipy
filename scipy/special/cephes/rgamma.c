/*                                             rgamma.c
 *
 *     Reciprocal Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, rgamma();
 *
 * y = rgamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns one divided by the Gamma function of the argument.
 *
 * The function is approximated by a Chebyshev expansion in
 * the interval [0,1].  Range reduction is by recurrence
 * for arguments between -34.034 and +34.84425627277176174.
 * 0 is returned for positive arguments outside this
 * range.  For arguments less than -34.034 the cosecant
 * reflection formula is applied; lograrithms are employed
 * to avoid unnecessary overflow.
 *
 * The reciprocal Gamma function has no singularities,
 * but overflow and underflow may occur for large arguments.
 * These conditions return either NPY_INFINITY or 0 with
 * appropriate sign.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC      -30,+30       4000       1.2e-16     1.8e-17
 *    IEEE     -30,+30      30000       1.1e-15     2.0e-16
 * For arguments less than -34.034 the peak error is on the
 * order of 5e-15 (DEC), excepting overflow or underflow.
 */

/*
 * Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140
 */

#include "mconf.h"

/* Chebyshev coefficients for reciprocal Gamma function
 * in interval 0 to 1.  Function is 1/(x Gamma(x)) - 1
 */

#ifdef UNK
static double R[] = {
    3.13173458231230000000E-17,
    -6.70718606477908000000E-16,
    2.20039078172259550000E-15,
    2.47691630348254132600E-13,
    -6.60074100411295197440E-12,
    5.13850186324226978840E-11,
    1.08965386454418662084E-9,
    -3.33964630686836942556E-8,
    2.68975996440595483619E-7,
    2.96001177518801696639E-6,
    -8.04814124978471142852E-5,
    4.16609138709688864714E-4,
    5.06579864028608725080E-3,
    -6.41925436109158228810E-2,
    -4.98558728684003594785E-3,
    1.27546015610523951063E-1
};
#endif

#ifdef DEC
static unsigned short R[] = {
    0022420, 0066376, 0176751, 0071636,
    0123501, 0051114, 0042104, 0131153,
    0024036, 0107013, 0126504, 0033361,
    0025613, 0070040, 0035174, 0162316,
    0126750, 0037060, 0077775, 0122202,
    0027541, 0177143, 0037675, 0105150,
    0030625, 0141311, 0075005, 0115436,
    0132017, 0067714, 0125033, 0014721,
    0032620, 0063707, 0105256, 0152643,
    0033506, 0122235, 0072757, 0170053,
    0134650, 0144041, 0015617, 0016143,
    0035332, 0066125, 0000776, 0006215,
    0036245, 0177377, 0137173, 0131432,
    0137203, 0073541, 0055645, 0141150,
    0136243, 0057043, 0026226, 0017362,
    0037402, 0115554, 0033441, 0012310
};
#endif

#ifdef IBMPC
static unsigned short R[] = {
    0x2e74, 0xdfbd, 0x0d9f, 0x3c82,
    0x964d, 0x8888, 0x2a49, 0xbcc8,
    0x86de, 0x75a8, 0xd1c1, 0x3ce3,
    0x9c9a, 0x074f, 0x6e04, 0x3d51,
    0xb490, 0x0fff, 0x07c6, 0xbd9d,
    0xb14d, 0x67f7, 0x3fcc, 0x3dcc,
    0xb364, 0x2f40, 0xb859, 0x3e12,
    0x633a, 0x9543, 0xedf9, 0xbe61,
    0xdab4, 0xf155, 0x0cf8, 0x3e92,
    0xfe05, 0xaebd, 0xd493, 0x3ec8,
    0xe38c, 0x2371, 0x1904, 0xbf15,
    0xc192, 0xa03f, 0x4d8a, 0x3f3b,
    0x7663, 0xf7cf, 0xbfdf, 0x3f74,
    0xb84d, 0x2b74, 0x6eec, 0xbfb0,
    0xc3de, 0x6592, 0x6bc4, 0xbf74,
    0x2299, 0x86e4, 0x536d, 0x3fc0
};
#endif

#ifdef MIEEE
static unsigned short R[] = {
    0x3c82, 0x0d9f, 0xdfbd, 0x2e74,
    0xbcc8, 0x2a49, 0x8888, 0x964d,
    0x3ce3, 0xd1c1, 0x75a8, 0x86de,
    0x3d51, 0x6e04, 0x074f, 0x9c9a,
    0xbd9d, 0x07c6, 0x0fff, 0xb490,
    0x3dcc, 0x3fcc, 0x67f7, 0xb14d,
    0x3e12, 0xb859, 0x2f40, 0xb364,
    0xbe61, 0xedf9, 0x9543, 0x633a,
    0x3e92, 0x0cf8, 0xf155, 0xdab4,
    0x3ec8, 0xd493, 0xaebd, 0xfe05,
    0xbf15, 0x1904, 0x2371, 0xe38c,
    0x3f3b, 0x4d8a, 0xa03f, 0xc192,
    0x3f74, 0xbfdf, 0xf7cf, 0x7663,
    0xbfb0, 0x6eec, 0x2b74, 0xb84d,
    0xbf74, 0x6bc4, 0x6592, 0xc3de,
    0x3fc0, 0x536d, 0x86e4, 0x2299
};
#endif

static char name[] = "rgamma";

extern double MAXLOG;


double rgamma(x)
double x;
{
    double w, y, z;
    int sign;

    if (x > 34.84425627277176174) {
        return exp(-lgam(x));
    }
    if (x < -34.034) {
	w = -x;
	z = sin(NPY_PI * w);
	if (z == 0.0)
	    return (0.0);
	if (z < 0.0) {
	    sign = 1;
	    z = -z;
	}
	else
	    sign = -1;

	y = log(w * z) - log(NPY_PI) + lgam(w);
	if (y < -MAXLOG) {
	    mtherr(name, UNDERFLOW);
	    return (sign * 0.0);
	}
	if (y > MAXLOG) {
	    mtherr(name, OVERFLOW);
	    return (sign * NPY_INFINITY);
	}
	return (sign * exp(y));
    }
    z = 1.0;
    w = x;

    while (w > 1.0) {		/* Downward recurrence */
	w -= 1.0;
	z *= w;
    }
    while (w < 0.0) {		/* Upward recurrence */
	z /= w;
	w += 1.0;
    }
    if (w == 0.0)		/* Nonpositive integer */
	return (0.0);
    if (w == 1.0)		/* Other integer */
	return (1.0 / z);

    y = w * (1.0 + chbevl(4.0 * w - 2.0, R, 16)) / z;
    return (y);
}
