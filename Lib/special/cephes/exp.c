/*							exp.c
 *
 *	Exponential function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, exp();
 *
 * y = exp( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns e (2.71828...) raised to the x power.
 *
 * Range reduction is accomplished by separating the argument
 * into an integer k and fraction f such that
 *
 *     x    k  f
 *    e  = 2  e.
 *
 * A Pade' form  1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
 * of degree 2/3 is used to approximate exp(f) in the basic
 * interval [-0.5, 0.5].
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       +- 88       50000       2.8e-17     7.0e-18
 *    IEEE      +- 708      40000       2.0e-16     5.6e-17
 *
 *
 * Error amplification in the exponential function can be
 * a serious matter.  The error propagation involves
 * exp( X(1+delta) ) = exp(X) ( 1 + X*delta + ... ),
 * which shows that a 1 lsb error in representing X produces
 * a relative error of X times 1 lsb in the function.
 * While the routine gives an accurate result for arguments
 * that are exactly represented by a double precision
 * computer number, the result contains amplified roundoff
 * error for large arguments not exactly represented.
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * exp underflow    x < MINLOG         0.0
 * exp overflow     x > MAXLOG         INFINITY
 *
 */

/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/


/*	Exponential function	*/

#include "mconf.h"

#ifdef UNK

static double P[] = {
 1.26177193074810590878E-4,
 3.02994407707441961300E-2,
 9.99999999999999999910E-1,
};
static double Q[] = {
 3.00198505138664455042E-6,
 2.52448340349684104192E-3,
 2.27265548208155028766E-1,
 2.00000000000000000009E0,
};
static double C1 = 6.93145751953125E-1;
static double C2 = 1.42860682030941723212E-6;
#endif

#ifdef DEC
static unsigned short P[] = {
0035004,0047156,0127442,0057502,
0036770,0033210,0063121,0061764,
0040200,0000000,0000000,0000000,
};
static unsigned short Q[] = {
0033511,0072665,0160662,0176377,
0036045,0070715,0124105,0132777,
0037550,0134114,0142077,0001637,
0040400,0000000,0000000,0000000,
};
static unsigned short sc1[] = {0040061,0071000,0000000,0000000};
#define C1 (*(double *)sc1)
static unsigned short sc2[] = {0033277,0137216,0075715,0057117};
#define C2 (*(double *)sc2)
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x4be8,0xd5e4,0x89cd,0x3f20,
0x2c7e,0x0cca,0x06d1,0x3f9f,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short Q[] = {
0x5fa0,0xbc36,0x2eb6,0x3ec9,
0xb6c0,0xb508,0xae39,0x3f64,
0xe074,0x9887,0x1709,0x3fcd,
0x0000,0x0000,0x0000,0x4000,
};
static unsigned short sc1[] = {0x0000,0x0000,0x2e40,0x3fe6};
#define C1 (*(double *)sc1)
static unsigned short sc2[] = {0xabca,0xcf79,0xf7d1,0x3eb7};
#define C2 (*(double *)sc2)
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3f20,0x89cd,0xd5e4,0x4be8,
0x3f9f,0x06d1,0x0cca,0x2c7e,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short Q[] = {
0x3ec9,0x2eb6,0xbc36,0x5fa0,
0x3f64,0xae39,0xb508,0xb6c0,
0x3fcd,0x1709,0x9887,0xe074,
0x4000,0x0000,0x0000,0x0000,
};
static unsigned short sc1[] = {0x3fe6,0x2e40,0x0000,0x0000};
#define C1 (*(double *)sc1)
static unsigned short sc2[] = {0x3eb7,0xf7d1,0xcf79,0xabca};
#define C2 (*(double *)sc2)
#endif

#ifndef ANSIPROT
double polevl(), p1evl(), floor(), ldexp();
int isnan(), isfinite();
#endif
extern double LOGE2, LOG2E, MAXLOG, MINLOG, MAXNUM;
#ifdef INFINITIES
extern double INFINITY;
#endif

double exp(x)
double x;
{
double px, xx;
int n;

#ifdef NANS
if( isnan(x) )
	return(x);
#endif
if( x > MAXLOG)
	{
#ifdef INFINITIES
	return( INFINITY );
#else
	mtherr( "exp", OVERFLOW );
	return( MAXNUM );
#endif
	}

if( x < MINLOG )
	{
#ifndef INFINITIES
	mtherr( "exp", UNDERFLOW );
#endif
	return(0.0);
	}

/* Express e**x = e**g 2**n
 *   = e**g e**( n loge(2) )
 *   = e**( g + n loge(2) )
 */
px = floor( LOG2E * x + 0.5 ); /* floor() truncates toward -infinity. */
n = px;
x -= px * C1;
x -= px * C2;

/* rational approximation for exponential
 * of the fractional part:
 * e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
 */
xx = x * x;
px = x * polevl( xx, P, 2 );
x =  px/( polevl( xx, Q, 3 ) - px );
x = 1.0 + 2.0 * x;

/* multiply by power of 2 */
x = ldexp( x, n );
return(x);
}
