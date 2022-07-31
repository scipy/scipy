/*
 * Functions for adding two double precision numbers with rounding to
 * infinity or rounding to negative infinity without using <fenv.h>.
 */
#ifndef ROUND_H
#define ROUND_H

#include <numpy/npy_math.h>
#include "_c99compat.h"
#include "cephes/dd_idefs.h"


double add_round_up(double a, double b)
{
    double s, err;

    if (sc_isnan(a) || sc_isnan(b)) {
	return NPY_NAN;
    }

    s = two_sum(a, b, &err);
    if (err > 0) {
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
    double s, err;

    if (sc_isnan(a) || sc_isnan(b)) {
	return NPY_NAN;
    }

    s = two_sum(a, b, &err);
    if (err < 0) {
	return npy_nextafter(s, -NPY_INFINITY);
    }
    else {
	return s;
    }
}


/* Helper code for testing _round.h. */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) || defined(__cplusplus)
/* We have C99, or C++11 or higher; both have fenv.h */
#include <fenv.h>
#else

int fesetround(int round)
{
    return -1;
}

int fegetround()
{
    return -1;
}

#define FE_UPWARD -1
#define FE_DOWNWARD -1

#endif


#endif /* _round.h */
