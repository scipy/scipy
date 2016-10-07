/*
 * Functions for adding two double precision numbers with rounding to
 * infinity or rounding to negative infinity without using <fenv.h>.
 *
 * References
 * ----------
 * [1] Hida, Li, Bailey, "Library for Double-Double and Quad-Double
 *     Arithmetic."
 */
#ifndef ROUND_H
#define ROUND_H

#include <numpy/npy_math.h>
#include "_c99compat.h"


double add_round_up(double a, double b)
{
    double s, e, v;

    /*
     * Use Algorithm 4 in [1] to compute s = fl(a + b) and 
     * e = err(a + b).
     */
    s = a + b;
    if (sc_isinf(s) && s < 0) {
	/* The sum is -oo; round up */
	return npy_nextafter(s, NPY_INFINITY);
    }
    v = s - a;
    e = (a - (s - v)) + (b - v);

    if (e > 0) {
	/* fl(a + b) rounded down */
	return npy_nextafter(s, NPY_INFINITY);
    }
    else {
	/* fl(a + b) rounded up or didn't round */
	return s;
    }
}


double add_round_down(double a, double b)
{
    double s, e, v;

    s = a + b;
    if (sc_isinf(s) && s > 0) {
	return npy_nextafter(s, -NPY_INFINITY);
    }
    v = s - a;
    e = (a - (s - v)) + (b - v);

    if (e < 0) {
	return npy_nextafter(s, -NPY_INFINITY);
    }
    else {
	return s;
    }
}


#endif /* _round.h */
