/*                                                     j0.c
 *
 *     Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / NPY_PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 */

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

#include "mconf.h"

static double PP[7] = {
    7.96936729297347051624E-4,
    8.28352392107440799803E-2,
    1.23953371646414299388E0,
    5.44725003058768775090E0,
    8.74716500199817011941E0,
    5.30324038235394892183E0,
    9.99999999999999997821E-1,
};

static double PQ[7] = {
    9.24408810558863637013E-4,
    8.56288474354474431428E-2,
    1.25352743901058953537E0,
    5.47097740330417105182E0,
    8.76190883237069594232E0,
    5.30605288235394617618E0,
    1.00000000000000000218E0,
};

static double QP[8] = {
    -1.13663838898469149931E-2,
    -1.28252718670509318512E0,
    -1.95539544257735972385E1,
    -9.32060152123768231369E1,
    -1.77681167980488050595E2,
    -1.47077505154951170175E2,
    -5.14105326766599330220E1,
    -6.05014350600728481186E0,
};

static double QQ[7] = {
    /*  1.00000000000000000000E0, */
    6.43178256118178023184E1,
    8.56430025976980587198E2,
    3.88240183605401609683E3,
    7.24046774195652478189E3,
    5.93072701187316984827E3,
    2.06209331660327847417E3,
    2.42005740240291393179E2,
};

static double YP[8] = {
    1.55924367855235737965E4,
    -1.46639295903971606143E7,
    5.43526477051876500413E9,
    -9.82136065717911466409E11,
    8.75906394395366999549E13,
    -3.46628303384729719441E15,
    4.42733268572569800351E16,
    -1.84950800436986690637E16,
};

static double YQ[7] = {
    /* 1.00000000000000000000E0, */
    1.04128353664259848412E3,
    6.26107330137134956842E5,
    2.68919633393814121987E8,
    8.64002487103935000337E10,
    2.02979612750105546709E13,
    3.17157752842975028269E15,
    2.50596256172653059228E17,
};

/*  5.783185962946784521175995758455807035071 */
static double DR1 = 5.78318596294678452118E0;

/* 30.47126234366208639907816317502275584842 */
static double DR2 = 3.04712623436620863991E1;

static double RP[4] = {
    -4.79443220978201773821E9,
    1.95617491946556577543E12,
    -2.49248344360967716204E14,
    9.70862251047306323952E15,
};

static double RQ[8] = {
    /* 1.00000000000000000000E0, */
    4.99563147152651017219E2,
    1.73785401676374683123E5,
    4.84409658339962045305E7,
    1.11855537045356834862E10,
    2.11277520115489217587E12,
    3.10518229857422583814E14,
    3.18121955943204943306E16,
    1.71086294081043136091E18,
};

extern double SQ2OPI;

double j0(x)
double x;
{
    double w, z, p, q, xn;

    if (x < 0)
	x = -x;

    if (x <= 5.0) {
	z = x * x;
	if (x < 1.0e-5)
	    return (1.0 - z / 4.0);

	p = (z - DR1) * (z - DR2);
	p = p * polevl(z, RP, 3) / p1evl(z, RQ, 8);
	return (p);
    }

    w = 5.0 / x;
    q = 25.0 / (x * x);
    p = polevl(q, PP, 6) / polevl(q, PQ, 6);
    q = polevl(q, QP, 7) / p1evl(q, QQ, 7);
    xn = x - NPY_PI_4;
    p = p * cos(xn) - w * q * sin(xn);
    return (p * SQ2OPI / sqrt(x));
}

/*                                                     y0() 2  */
/* Bessel function of second kind, order zero  */

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / NPY_PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / NPY_PI
 * = 0.073804295108687225.
 */

/*
 * #define NPY_PI_4 .78539816339744830962
 * #define SQ2OPI .79788456080286535588
 */

double y0(x)
double x;
{
    double w, z, p, q, xn;

    if (x <= 5.0) {
	if (x == 0.0) {
	    sf_error("y0", SF_ERROR_SINGULAR, NULL);
	    return -INFINITY;
	}
	else if (x < 0.0) {
	    sf_error("y0", SF_ERROR_DOMAIN, NULL);
	    return NAN;
	}
	z = x * x;
	w = polevl(z, YP, 7) / p1evl(z, YQ, 7);
	w += NPY_2_PI * log(x) * j0(x);
	return (w);
    }

    w = 5.0 / x;
    z = 25.0 / (x * x);
    p = polevl(z, PP, 6) / polevl(z, PQ, 6);
    q = polevl(z, QP, 7) / p1evl(z, QQ, 7);
    xn = x - NPY_PI_4;
    p = p * sin(xn) + w * q * cos(xn);
    return (p * SQ2OPI / sqrt(x));
}
