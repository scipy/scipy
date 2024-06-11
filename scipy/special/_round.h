/*
 * Functions for adding two double precision numbers with rounding to
 * infinity or rounding to negative infinity without using <fenv.h>.
 */
#ifndef ROUND_H
#define ROUND_H

#include <math.h>


/* Computes fl(a+b) and err(a+b).  */
static inline double two_sum(double a, double b, double *err)
{
    volatile double s = a + b;
    volatile double c = s - a;
    volatile double d = b - c;
    volatile double e = s - c;
    *err = (a - e) + d;
    return s;
}


double add_round_up(double a, double b)
{
    double s, err;

    if (isnan(a) || isnan(b)) {
	return NAN;
    }

    s = two_sum(a, b, &err);
    if (err > 0) {
	/* fl(a + b) rounded down */
	return nextafter(s, INFINITY);
    }
    else {
	/* fl(a + b) rounded up or didn't round */
	return s;
    }
}


double add_round_down(double a, double b)
{
    double s, err;

    if (isnan(a) || isnan(b)) {
	return NAN;
    }

    s = two_sum(a, b, &err);
    if (err < 0) {
	return nextafter(s, -INFINITY);
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
