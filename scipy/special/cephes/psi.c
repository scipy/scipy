/*                                                     psi.c
 *
 *     Psi (digamma) function
 *
 *
 * SYNOPSIS:
 *
 * double x, y, psi();
 *
 * y = psi( x );
 *
 *
 * DESCRIPTION:
 *
 *              d      -
 *   psi(x)  =  -- ln | (x)
 *              dx
 *
 * is the logarithmic derivative of the gamma function.
 * For integer x,
 *                   n-1
 *                    -
 * psi(n) = -EUL  +   >  1/k.
 *                    -
 *                   k=1
 *
 * This formula is used for 0 < n <= 10.  If x is negative, it
 * is transformed to a positive argument by the reflection
 * formula  psi(1-x) = psi(x) + pi cot(pi x).
 * For general positive x, the argument is made greater than 10
 * using the recurrence  psi(x+1) = psi(x) + 1/x.
 * Then the following asymptotic expansion is applied:
 *
 *                           inf.   B
 *                            -      2k
 * psi(x) = log(x) - 1/2x -   >   -------
 *                            -        2k
 *                           k=1   2k x
 *
 * where the B2k are Bernoulli numbers.
 *
 * ACCURACY:
 *    Relative error (except absolute when |psi| < 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       1.3e-15     1.4e-16
 *    IEEE      -30,0       40000       1.5e-15     2.2e-16
 *
 * ERROR MESSAGES:
 *     message         condition      value returned
 * psi singularity    x integer <=0      NPY_INFINITY
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1992, 2000 by Stephen L. Moshier
 */

/*
 * Code for the rational approximation on [1, 2] is:
 *
 * (C) Copyright John Maddock 2006.
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "mconf.h"

static double A[] = {
    8.33333333333333333333E-2,
    -2.10927960927960927961E-2,
    7.57575757575757575758E-3,
    -4.16666666666666666667E-3,
    3.96825396825396825397E-3,
    -8.33333333333333333333E-3,
    8.33333333333333333333E-2
};


double digamma_imp_1_2(double x)
{
    /*
     * Rational approximation on [1, 2] taken from Boost.
     *
     * Now for the approximation, we use the form:
     *
     * digamma(x) = (x - root) * (Y + R(x-1))
     *
     * Where root is the location of the positive root of digamma,
     * Y is a constant, and R is optimised for low absolute error
     * compared to Y.
     *
     * Maximum Deviation Found:               1.466e-18
     * At double precision, max error found:  2.452e-17
     */
    double r, g;

    static const float Y = 0.99558162689208984f;

    static const double root1 = 1569415565.0 / 1073741824.0;
    static const double root2 = (381566830.0 / 1073741824.0) / 1073741824.0;
    static const double root3 = 0.9016312093258695918615325266959189453125e-19;

   static double P[] = {
       -0.0020713321167745952,
       -0.045251321448739056,
       -0.28919126444774784,
       -0.65031853770896507,
       -0.32555031186804491,
       0.25479851061131551
   };
   static double Q[] = {
       -0.55789841321675513e-6,
       0.0021284987017821144,
       0.054151797245674225,
       0.43593529692665969,
       1.4606242909763515,
       2.0767117023730469,
       1.0
   };
   g = x - root1;
   g -= root2;
   g -= root3;
   r = polevl(x - 1.0, P, 5) / polevl(x - 1.0, Q, 6);

   return g * Y + g * r;
}


double psi(x)
double x;
{
    double p, q, nz, s, w, y, z;
    int i, n, negative;

    negative = 0;
    nz = 0.0;

    if (x <= 0.0) {
	negative = 1;
	q = x;
	p = floor(q);
	if (p == q) {
	    mtherr("psi", SING);
	    return (NPY_INFINITY);
	}
	/* Remove the zeros of tan(NPY_PI x)
	 * by subtracting the nearest integer from x
	 */
	nz = q - p;
	if (nz != 0.5) {
	    if (nz > 0.5) {
		p += 1.0;
		nz = q - p;
	    }
	    nz = NPY_PI / tan(NPY_PI * nz);
	}
	else {
	    nz = 0.0;
	}
	x = 1.0 - x;
    }

    /* check for positive integer up to 10 */
    if ((x <= 10.0) && (x == floor(x))) {
	y = 0.0;
	n = x;
	for (i = 1; i < n; i++) {
	    w = i;
	    y += 1.0 / w;
	}
	y -= NPY_EULER;
	goto done;
    }

    if ((1.0 <= x) && (x <= 2.0)) {
	y = digamma_imp_1_2(x);
	goto done;
    }

    s = x;
    w = 0.0;
    while (s < 10.0) {
	w += 1.0 / s;
	s += 1.0;
    }

    if (s < 1.0e17) {
	z = 1.0 / (s * s);
	y = z * polevl(z, A, 6);
    }
    else
	y = 0.0;

    y = log(s) - (0.5 / s) - y - w;

  done:

    if (negative) {
	y -= nz;
    }

    return (y);
}
