/*							log10.c
 *
 *	Common logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, log10();
 *
 * y = log10( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns logarithm to the base 10 of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  The logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)/Q(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0.5, 2.0     30000      1.5e-16     5.0e-17
 *    IEEE      0, MAXNUM    30000      1.4e-16     4.8e-17
 *    DEC       1, MAXNUM    50000      2.5e-17     6.0e-18
 *
 * In the tests over the interval [1, MAXNUM], the logarithms
 * of the random arguments were uniformly distributed over
 * [0, MAXLOG].
 *
 * ERROR MESSAGES:
 *
 * log10 singularity:  x = 0; returns -INFINITY
 * log10 domain:       x < 0; returns NAN
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include "mconf.h"
static char fname[] = {"log10"};

/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 */
#ifdef UNK
static double P[] = {
  4.58482948458143443514E-5,
  4.98531067254050724270E-1,
  6.56312093769992875930E0,
  2.97877425097986925891E1,
  6.06127134467767258030E1,
  5.67349287391754285487E1,
  1.98892446572874072159E1
};
static double Q[] = {
/* 1.00000000000000000000E0, */
  1.50314182634250003249E1,
  8.27410449222435217021E1,
  2.20664384982121929218E2,
  3.07254189979530058263E2,
  2.14955586696422947765E2,
  5.96677339718622216300E1
};
#endif

#ifdef DEC
static unsigned short P[] = {
0034500,0046473,0051374,0135174,
0037777,0037566,0145712,0150321,
0040722,0002426,0031543,0123107,
0041356,0046513,0170752,0004346,
0041562,0071553,0023536,0163343,
0041542,0170221,0024316,0114216,
0041237,0016454,0046611,0104602
};
static unsigned short Q[] = {
/*0040200,0000000,0000000,0000000,*/
0041160,0100260,0067736,0102424,
0041645,0075552,0036563,0147072,
0042134,0125025,0021132,0025320,
0042231,0120211,0046030,0103271,
0042126,0172241,0052151,0120426,
0041556,0125702,0072116,0047103
};
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x974f,0x6a5f,0x09a7,0x3f08,
0x5a1a,0xd979,0xe7ee,0x3fdf,
0x74c9,0xc66c,0x40a2,0x401a,
0x411d,0x7e3d,0xc9a9,0x403d,
0xdcdc,0x64eb,0x4e6d,0x404e,
0xd312,0x2519,0x5e12,0x404c,
0x3130,0x89b1,0xe3a5,0x4033
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xd0a2,0x0dfb,0x1016,0x402e,
0x79c7,0x47ae,0xaf6d,0x4054,
0x455a,0xa44b,0x9542,0x406b,
0x10d7,0x2983,0x3411,0x4073,
0x3423,0x2a8d,0xde94,0x406a,
0xc9c8,0x4e89,0xd578,0x404d
};
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3f08,0x09a7,0x6a5f,0x974f,
0x3fdf,0xe7ee,0xd979,0x5a1a,
0x401a,0x40a2,0xc66c,0x74c9,
0x403d,0xc9a9,0x7e3d,0x411d,
0x404e,0x4e6d,0x64eb,0xdcdc,
0x404c,0x5e12,0x2519,0xd312,
0x4033,0xe3a5,0x89b1,0x3130
};
static unsigned short Q[] = {
0x402e,0x1016,0x0dfb,0xd0a2,
0x4054,0xaf6d,0x47ae,0x79c7,
0x406b,0x9542,0xa44b,0x455a,
0x4073,0x3411,0x2983,0x10d7,
0x406a,0xde94,0x2a8d,0x3423,
0x404d,0xd578,0x4e89,0xc9c8
};
#endif

#define SQRTH 0.70710678118654752440
#define L102A 3.0078125E-1
#define L102B 2.48745663981195213739E-4
#define L10EA 4.3359375E-1
#define L10EB 7.00731903251827651129E-4

#ifndef ANSIPROT
double frexp(), ldexp(), polevl(), p1evl();
int isnan(), isfinite();
#endif
extern double LOGE2, SQRT2, INFINITY, NAN;

double log10(x)
double x;
{
VOLATILE double z;
double y;
#ifdef DEC
short *q;
#endif
int e;

#ifdef NANS
if( isnan(x) )
	return(x);
#endif
#ifdef INFINITIES
if( x == INFINITY )
	return(x);
#endif
/* Test for domain */
if( x <= 0.0 )
	{
	if( x == 0.0 )
	        {
		mtherr( fname, SING );
		return( -INFINITY );
	        }
	else
	        {
		mtherr( fname, DOMAIN );
		return( NAN );
	        }
	}

/* separate mantissa from exponent */

#ifdef DEC
q = (short *)&x;
e = *q;			/* short containing exponent */
e = ((e >> 7) & 0377) - 0200;	/* the exponent */
*q &= 0177;	/* strip exponent from x */
*q |= 040000;	/* x now between 0.5 and 1 */
#endif

#ifdef IBMPC
x = frexp( x, &e );
/*
q = (short *)&x;
q += 3;
e = *q;
e = ((e >> 4) & 0x0fff) - 0x3fe;
*q &= 0x0f;
*q |= 0x3fe0;
*/
#endif

/* Equivalent C language standard library function: */
#ifdef UNK
x = frexp( x, &e );
#endif

#ifdef MIEEE
x = frexp( x, &e );
#endif

/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */

if( x < SQRTH )
	{
	e -= 1;
	x = ldexp( x, 1 ) - 1.0; /*  2x - 1  */
	}	
else
	{
	x = x - 1.0;
	}


/* rational form */
z = x*x;
y = x * ( z * polevl( x, P, 6 ) / p1evl( x, Q, 6 ) );
y = y - ldexp( z, -1 );   /*  y - 0.5 * x**2  */

/* multiply log of fraction by log10(e)
 * and base 2 exponent by log10(2)
 */
z = (x + y) * L10EB;  /* accumulate terms in order of size */
z += y * L10EA;
z += x * L10EA;
z += e * L102B;
z += e * L102A;


return( z );
}
