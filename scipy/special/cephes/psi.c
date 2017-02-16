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
