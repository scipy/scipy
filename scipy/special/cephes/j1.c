/*                                                     j1.c
 *
 *     Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 */

/*
 * #define PIO4 .78539816339744830962
 * #define THPIO4 2.35619449019234492885
 * #define SQ2OPI .79788456080286535588
 */

#include "mconf.h"

static double RP[4] = {
    -8.99971225705559398224E8,
    4.52228297998194034323E11,
    -7.27494245221818276015E13,
    3.68295732863852883286E15,
};

static double RQ[8] = {
    /* 1.00000000000000000000E0, */
    6.20836478118054335476E2,
    2.56987256757748830383E5,
    8.35146791431949253037E7,
    2.21511595479792499675E10,
    4.74914122079991414898E12,
    7.84369607876235854894E14,
    8.95222336184627338078E16,
    5.32278620332680085395E18,
};

static double PP[7] = {
    7.62125616208173112003E-4,
    7.31397056940917570436E-2,
    1.12719608129684925192E0,
    5.11207951146807644818E0,
    8.42404590141772420927E0,
    5.21451598682361504063E0,
    1.00000000000000000254E0,
};

static double PQ[7] = {
    5.71323128072548699714E-4,
    6.88455908754495404082E-2,
    1.10514232634061696926E0,
    5.07386386128601488557E0,
    8.39985554327604159757E0,
    5.20982848682361821619E0,
    9.99999999999999997461E-1,
};

static double QP[8] = {
    5.10862594750176621635E-2,
    4.98213872951233449420E0,
    7.58238284132545283818E1,
    3.66779609360150777800E2,
    7.10856304998926107277E2,
    5.97489612400613639965E2,
    2.11688757100572135698E2,
    2.52070205858023719784E1,
};

static double QQ[7] = {
    /* 1.00000000000000000000E0, */
    7.42373277035675149943E1,
    1.05644886038262816351E3,
    4.98641058337653607651E3,
    9.56231892404756170795E3,
    7.99704160447350683650E3,
    2.82619278517639096600E3,
    3.36093607810698293419E2,
};

static double YP[6] = {
    1.26320474790178026440E9,
    -6.47355876379160291031E11,
    1.14509511541823727583E14,
    -8.12770255501325109621E15,
    2.02439475713594898196E17,
    -7.78877196265950026825E17,
};

static double YQ[8] = {
    /* 1.00000000000000000000E0, */
    5.94301592346128195359E2,
    2.35564092943068577943E5,
    7.34811944459721705660E7,
    1.87601316108706159478E10,
    3.88231277496238566008E12,
    6.20557727146953693363E14,
    6.87141087355300489866E16,
    3.97270608116560655612E18,
};


static double Z1 = 1.46819706421238932572E1;
static double Z2 = 4.92184563216946036703E1;

extern double THPIO4, SQ2OPI;

double j1(x)
double x;
{
    double w, z, p, q, xn;

    w = x;
    if (x < 0)
	return -j1(-x);

    if (w <= 5.0) {
	z = x * x;
	w = polevl(z, RP, 3) / p1evl(z, RQ, 8);
	w = w * x * (z - Z1) * (z - Z2);
	return (w);
    }

    w = 5.0 / x;
    z = w * w;
    p = polevl(z, PP, 6) / polevl(z, PQ, 6);
    q = polevl(z, QP, 7) / p1evl(z, QQ, 7);
    xn = x - THPIO4;
    p = p * cos(xn) - w * q * sin(xn);
    return (p * SQ2OPI / sqrt(x));
}


double y1(x)
double x;
{
    double w, z, p, q, xn;

    if (x <= 5.0) {
	if (x == 0.0) {
	    sf_error("y1", SF_ERROR_SINGULAR, NULL);
	    return -NPY_INFINITY;
	}
	else if (x <= 0.0) {
	    sf_error("y1", SF_ERROR_DOMAIN, NULL);
	    return NPY_NAN;
	}
	z = x * x;
	w = x * (polevl(z, YP, 5) / p1evl(z, YQ, 8));
	w += NPY_2_PI * (j1(x) * log(x) - 1.0 / x);
	return (w);
    }

    w = 5.0 / x;
    z = w * w;
    p = polevl(z, PP, 6) / polevl(z, PQ, 6);
    q = polevl(z, QP, 7) / p1evl(z, QQ, 7);
    xn = x - THPIO4;
    p = p * sin(xn) + w * q * cos(xn);
    return (p * SQ2OPI / sqrt(x));
}
