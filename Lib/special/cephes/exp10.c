/*							exp10.c
 *
 *	Base 10 exponential function
 *      (Common antilogarithm)
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, exp10();
 *
 * y = exp10( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns 10 raised to the x power.
 *
 * Range reduction is accomplished by expressing the argument
 * as 10**x = 2**n 10**f, with |f| < 0.5 log10(2).
 * The Pade' form
 *
 *    1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
 *
 * is used to approximate 10**f.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -307,+307    30000       2.2e-16     5.5e-17
 * Test result from an earlier version (2.1):
 *    DEC       -38,+38     70000       3.1e-17     7.0e-18
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * exp10 underflow    x < -MAXL10        0.0
 * exp10 overflow     x > MAXL10       MAXNUM
 *
 * DEC arithmetic: MAXL10 = 38.230809449325611792.
 * IEEE arithmetic: MAXL10 = 308.2547155599167.
 *
 */

/*
Cephes Math Library Release 2.2:  January, 1991
Copyright 1984, 1991 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


#include "mconf.h"

#ifdef UNK
static double P[] = {
 4.09962519798587023075E-2,
 1.17452732554344059015E1,
 4.06717289936872725516E2,
 2.39423741207388267439E3,
};
static double Q[] = {
/* 1.00000000000000000000E0,*/
 8.50936160849306532625E1,
 1.27209271178345121210E3,
 2.07960819286001865907E3,
};
/* static double LOG102 = 3.01029995663981195214e-1; */
static double LOG210 = 3.32192809488736234787e0;
static double LG102A = 3.01025390625000000000E-1;
static double LG102B = 4.60503898119521373889E-6;
/* static double MAXL10 = 38.230809449325611792; */
static double MAXL10 = 308.2547155599167;
#endif

#ifdef DEC
static unsigned short P[] = {
0037047,0165657,0114061,0067234,
0041073,0166243,0123052,0144643,
0042313,0055720,0024032,0047443,
0043025,0121714,0070232,0050007,
};
static unsigned short Q[] = {
/*0040200,0000000,0000000,0000000,*/
0041652,0027756,0071216,0050075,
0042637,0001367,0077263,0136017,
0043001,0174673,0024157,0133416,
};
/*
static unsigned short L102[] = {0037632,0020232,0102373,0147770};
#define LOG102 *(double *)L102
*/
static unsigned short L210[] = {0040524,0115170,0045715,0015613};
#define LOG210 *(double *)L210
static unsigned short L102A[] = {0037632,0020000,0000000,0000000,};
#define LG102A *(double *)L102A
static unsigned short L102B[] = {0033632,0102373,0147767,0114220,};
#define LG102B *(double *)L102B
static unsigned short MXL[] = {0041430,0166131,0047761,0154130,};
#define MAXL10 ( *(double *)MXL )
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x2dd4,0xf306,0xfd75,0x3fa4,
0x5934,0x74c5,0x7d94,0x4027,
0x49e4,0x0503,0x6b7a,0x4079,
0x4a01,0x8e13,0xb479,0x40a2,
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xca08,0xce51,0x45fd,0x4055,
0x7782,0xefd6,0xe05e,0x4093,
0xf6e2,0x650d,0x3f37,0x40a0,
};
/*
static unsigned short L102[] = {0x79ff,0x509f,0x4413,0x3fd3};
#define LOG102 *(double *)L102
*/
static unsigned short L210[] = {0xa371,0x0979,0x934f,0x400a};
#define LOG210 *(double *)L210
static unsigned short L102A[] = {0x0000,0x0000,0x4400,0x3fd3,};
#define LG102A *(double *)L102A
static unsigned short L102B[] = {0xf312,0x79fe,0x509f,0x3ed3,};
#define LG102B *(double *)L102B
static double MAXL10 = 308.2547155599167;
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3fa4,0xfd75,0xf306,0x2dd4,
0x4027,0x7d94,0x74c5,0x5934,
0x4079,0x6b7a,0x0503,0x49e4,
0x40a2,0xb479,0x8e13,0x4a01,
};
static unsigned short Q[] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4055,0x45fd,0xce51,0xca08,
0x4093,0xe05e,0xefd6,0x7782,
0x40a0,0x3f37,0x650d,0xf6e2,
};
/*
static unsigned short L102[] = {0x3fd3,0x4413,0x509f,0x79ff};
#define LOG102 *(double *)L102
*/
static unsigned short L210[] = {0x400a,0x934f,0x0979,0xa371};
#define LOG210 *(double *)L210
static unsigned short L102A[] = {0x3fd3,0x4400,0x0000,0x0000,};
#define LG102A *(double *)L102A
static unsigned short L102B[] = {0x3ed3,0x509f,0x79fe,0xf312,};
#define LG102B *(double *)L102B
static double MAXL10 = 308.2547155599167;
#endif

#ifndef ANSIPROT
double floor(), ldexp(), polevl(), p1evl();
int isnan(), isfinite();
#endif
extern double MAXNUM;
#ifdef INFINITIES
extern double INFINITY;
#endif

double exp10(double x)
{
double px, xx;
short n;

#ifdef NANS
if( isnan(x) )
	return(x);
#endif
if( x > MAXL10 )
	{
#ifdef INFINITIES
	return( INFINITY );
#else
	mtherr( "exp10", OVERFLOW );
	return( MAXNUM );
#endif
	}

if( x < -MAXL10 )	/* Would like to use MINLOG but can't */
	{
	mtherr( "exp10", UNDERFLOW );
	return(0.0);
	}

/* Express 10**x = 10**g 2**n
 *   = 10**g 10**( n log10(2) )
 *   = 10**( g + n log10(2) )
 */
px = floor( LOG210 * x + 0.5 );
n = px;
x -= px * LG102A;
x -= px * LG102B;

/* rational approximation for exponential
 * of the fractional part:
 * 10**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
 */
xx = x * x;
px = x * polevl( xx, P, 3 );
x =  px/( p1evl( xx, Q, 3 ) - px );
x = 1.0 + ldexp( x, 1 );

/* multiply by power of 2 */
x = ldexp( x, n );

return(x);
}
