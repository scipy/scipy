/*                                                     exp2.c
 *
 *     Base 2 exponential function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, exp2();
 *
 * y = exp2( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 2 raised to the x power.
 *
 * Range reduction is accomplished by separating the argument
 * into an integer k and fraction f such that
 *     x    k  f
 *    2  = 2  2.
 *
 * A Pade' form
 *
 *   1 + 2x P(x**2) / (Q(x**2) - x P(x**2) )
 *
 * approximates 2**x in the basic range [-0.5, 0.5].
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    -1022,+1024   30000       1.8e-16     5.4e-17
 *
 *
 * See exp.c for comments on error amplification.
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * exp underflow    x < -MAXL2        0.0
 * exp overflow     x > MAXL2         NPY_INFINITY
 *
 * For DEC arithmetic, MAXL2 = 127.
 * For IEEE arithmetic, MAXL2 = 1024.
 */


/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 */



#include "mconf.h"

#ifdef UNK
static double P[] = {
    2.30933477057345225087E-2,
    2.02020656693165307700E1,
    1.51390680115615096133E3,
};

static double Q[] = {
    /* 1.00000000000000000000E0, */
    2.33184211722314911771E2,
    4.36821166879210612817E3,
};

#define MAXL2 1024.0
#define MINL2 -1024.0
#endif

#ifdef DEC
static unsigned short P[] = {
    0036675, 0027102, 0122327, 0053227,
    0041241, 0116724, 0115412, 0157355,
    0042675, 0036404, 0101733, 0132226,
};

static unsigned short Q[] = {
    /*0040200,0000000,0000000,0000000, */
    0042151, 0027450, 0077732, 0160744,
    0043210, 0100661, 0077550, 0056560,
};

#define MAXL2 127.0
#define MINL2 -127.0
#endif

#ifdef IBMPC
static unsigned short P[] = {
    0xead3, 0x549a, 0xa5c8, 0x3f97,
    0x5bde, 0x9361, 0x33ba, 0x4034,
    0x7693, 0x907b, 0xa7a0, 0x4097,
};

static unsigned short Q[] = {
    /*0x0000,0x0000,0x0000,0x3ff0, */
    0x5c3c, 0x0ffb, 0x25e5, 0x406d,
    0x0bae, 0x2fed, 0x1036, 0x40b1,
};

#define MAXL2 1024.0
#define MINL2 -1022.0
#endif

#ifdef MIEEE
static unsigned short P[] = {
    0x3f97, 0xa5c8, 0x549a, 0xead3,
    0x4034, 0x33ba, 0x9361, 0x5bde,
    0x4097, 0xa7a0, 0x907b, 0x7693,
};

static unsigned short Q[] = {
    /*0x3ff0,0x0000,0x0000,0x0000, */
    0x406d, 0x25e5, 0x0ffb, 0x5c3c,
    0x40b1, 0x1036, 0x2fed, 0x0bae,
};

#define MAXL2 1024.0
#define MINL2 -1022.0
#endif

double exp2(double x)
{
    double px, xx;
    short n;

    if (cephes_isnan(x))
	return (x);
    if (x > MAXL2) {
	return (NPY_INFINITY);
    }

    if (x < MINL2) {
	return (0.0);
    }

    xx = x;			/* save x */
    /* separate into integer and fractional parts */
    px = floor(x + 0.5);
    n = px;
    x = x - px;

    /* rational approximation
     * exp2(x) = 1 +  2xP(xx)/(Q(xx) - P(xx))
     * where xx = x**2
     */
    xx = x * x;
    px = x * polevl(xx, P, 2);
    x = px / (p1evl(xx, Q, 2) - px);
    x = 1.0 + ldexp(x, 1);

    /* scale by power of 2 */
    x = ldexp(x, n);
    return (x);
}
