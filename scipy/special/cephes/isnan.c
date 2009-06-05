/*							isnan()
 *							signbit()
 *							isfinite()
 *
 *	Floating point numeric utilities
 *
 *
 *
 * SYNOPSIS:
 *
 * double ceil(), floor(), frexp(), ldexp(); --- gone
 * int signbit(), isnan(), isfinite();
 * double x, y;
 * int expnt, n;
 *
 * y = floor(x);             -gone 
 * y = ceil(x);              -gone
 * y = frexp( x, &expnt );   -gone
 * y = ldexp( x, n );        -gone
 * n = signbit(x);
 * n = isnan(x);
 * n = isfinite(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * All four routines return a double precision floating point
 * result.
 *
 * floor() returns the largest integer less than or equal to x.
 * It truncates toward minus infinity.
 *
 * ceil() returns the smallest integer greater than or equal
 * to x.  It truncates toward plus infinity.
 *
 * frexp() extracts the exponent from x.  It returns an integer
 * power of two to expnt and the significand between 0.5 and 1
 * to y.  Thus  x = y * 2**expn.
 *
 * ldexp() multiplies x by 2**n.
 *
 * signbit(x) returns 1 if the sign bit of x is 1, else 0.
 *
 * These functions are part of the standard C run time library
 * for many but not all C compilers.  The ones supplied are
 * written in C for either DEC or IEEE arithmetic.  They should
 * be used only if your compiler library does not already have
 * them.
 *
 * The IEEE versions assume that denormal numbers are implemented
 * in the arithmetic.  Some modifications will be required if
 * the arithmetic has abrupt rather than gradual underflow.
 */


/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include "mconf.h"

/* XXX: horrible hacks, but those cephes macros are buggy and just plain ugly anywa.
 * We should use npy_* macros instead once npy_math can be used reliably by
 * packages outside numpy
 */
#undef isnan
#undef signbit
#undef isfinite

#define isnan(x) ((x) != (x))

int cephes_isnan(double x)
{
	return isnan(x);
}

int isfinite(double x)
{
	return !isnan((x) + (-x));
}

static int isbigendian(void)
{
    const union {
        npy_uint32 i;
        char c[4];
    } bint = {0x01020304};

    if (bint.c[0] == 1) {
        return 1;
    }
    return 0;
}

int signbit(double x)
{
    union
    {
        double d;
        short s[4];
        int i[2];
    } u;

    u.d = x;

    /*
     * Tuis is stupid, we test for endianness every time, but that the easiest
     * way I can see without using platform checks - for scipy 0.8.0, we should
     * use npy_math
     */
#if SIZEOF_INT == 4
    if (isbigendian()) {
	return u.i[1] < 0;
    } else {
	return u.i[0] < 0;
    }

#else  /* SIZEOF_INT != 4 */

    if (isbigendian()) {
	return u.s[3] < 0;
    } else {
	return u.s[0] < 0;
    }
#endif  /* SIZEOF_INT */
}
