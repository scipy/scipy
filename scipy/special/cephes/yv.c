/*
 * Cephes Math Library Release 2.8: June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */

#include "mconf.h"

#define DEBUG 0

static double stop = 1.37e-17;
extern double MACHEP;


double onef2(double a, double b, double c, double x, double *err)
{
    double n, a0, sum, t;
    double an, bn, cn, max, z;

    an = a;
    bn = b;
    cn = c;
    a0 = 1.0;
    sum = 1.0;
    n = 1.0;
    t = 1.0;
    max = 0.0;

    do {
	if (an == 0)
	    goto done;
	if (bn == 0)
	    goto error;
	if (cn == 0)
	    goto error;
	if ((a0 > 1.0e34) || (n > 200))
	    goto error;
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = fabs(a0);
	if (z > max)
	    max = z;
	if (sum != 0)
	    t = fabs(a0 / sum);
	else
	    t = z;
    }
    while (t > stop);

  done:

    *err = fabs(MACHEP * max / sum);

#if DEBUG
    printf(" onef2 cancellation error %.5E\n", *err);
#endif

    goto xit;

  error:
#if DEBUG
    printf("onef2 does not converge\n");
#endif
    *err = 1.0e38;

  xit:

#if DEBUG
    printf("onef2( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n,
	   sum);
#endif
    return (sum);
}


double threef0(double a, double b, double c, double x, double *err)
{
    double n, a0, sum, t, conv, conv1;
    double an, bn, cn, max, z;

    an = a;
    bn = b;
    cn = c;
    a0 = 1.0;
    sum = 1.0;
    n = 1.0;
    t = 1.0;
    max = 0.0;
    conv = 1.0e38;
    conv1 = conv;

    do {
	if (an == 0.0)
	    goto done;
	if (bn == 0.0)
	    goto done;
	if (cn == 0.0)
	    goto done;
	if ((a0 > 1.0e34) || (n > 200))
	    goto error;
	a0 *= (an * bn * cn * x) / n;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = fabs(a0);
	if (z > max)
	    max = z;
	if (z >= conv) {
	    if ((z < max) && (z > conv1))
		goto done;
	}
	conv1 = conv;
	conv = z;
	sum += a0;
	if (sum != 0)
	    t = fabs(a0 / sum);
	else
	    t = z;
    }
    while (t > stop);

  done:

    t = fabs(MACHEP * max / sum);
#if DEBUG
    printf(" threef0 cancellation error %.5E\n", t);
#endif

    max = fabs(conv / sum);
    if (max > t)
	t = max;
#if DEBUG
    printf(" threef0 convergence %.5E\n", max);
#endif

    goto xit;

  error:
#if DEBUG
    printf("threef0 does not converge\n");
#endif
    t = 1.0e38;

  xit:

#if DEBUG
    printf("threef0( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n,
	   sum);
#endif

    *err = t;
    return (sum);
}


/*
 * Bessel function of noninteger order
 */
double yv(double v, double x)
{
    double y, t;
    int n;

    n = v;
    if (n == v) {
	y = yn(n, x);
	return (y);
    }
    else if (v == floor(v)) {
        /* Zero in denominator. */
	mtherr("yv", DOMAIN);
        return NPY_NAN;
    }

    t = NPY_PI * v;
    y = (cos(t) * jv(v, x) - jv(-v, x)) / sin(t);

    if (cephes_isinf(y)) {
        if (v > 0) {
            mtherr("yv", OVERFLOW);
            return -NPY_INFINITY;
        }
        else if (v < -1e10) {
            /* Whether it's +inf or -inf is numerically ill-defined. */
            mtherr("yv", DOMAIN);
            return NPY_NAN;
        }
    }

    return (y);
}
