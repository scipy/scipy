/*							tandg.c
 *
 *	Circular tangent of argument in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, tandg();
 *
 * y = tandg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the argument x in degrees.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC      0,10          8000      3.4e-17      1.2e-17
 *    IEEE     0,10         30000      3.2e-16      8.4e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * tandg total loss   x > 8.0e14 (DEC)      0.0
 *                    x > 1.0e14 (IEEE)
 * tandg singularity  x = 180 k  +  90     MAXNUM
 */
/*							cotdg.c
 *
 *	Circular cotangent of argument in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cotdg();
 *
 * y = cotdg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular cotangent of the argument x in degrees.
 *
 * Range reduction is modulo pi/4.  A rational function
 *       x + x**3 P(x**2)/Q(x**2)
 * is employed in the basic interval [0, pi/4].
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * cotdg total loss   x > 8.0e14 (DEC)      0.0
 *                    x > 1.0e14 (IEEE)
 * cotdg singularity  x = 180 k            MAXNUM
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"

#ifdef UNK
static double P[] = {
-1.30936939181383777646E4,
 1.15351664838587416140E6,
-1.79565251976484877988E7
};
static double Q[] = {
/* 1.00000000000000000000E0,*/
 1.36812963470692954678E4,
-1.32089234440210967447E6,
 2.50083801823357915839E7,
-5.38695755929454629881E7
};
static double PI180 = 1.74532925199432957692E-2;
static double lossth = 1.0e14;
#endif

#ifdef DEC
static unsigned short P[] = {
0143514,0113306,0111171,0174674,
0045214,0147545,0027744,0167346,
0146210,0177526,0114514,0105660
};
static unsigned short Q[] = {
/*0040200,0000000,0000000,0000000,*/
0043525,0142457,0072633,0025617,
0145241,0036742,0140525,0162256,
0046276,0146176,0013526,0143573,
0146515,0077401,0162762,0150607
};
static unsigned short P1[] = {0036616,0175065,0011224,0164711};
#define PI180 *(double *)P1
static double lossth = 8.0e14;
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x3f38,0xd24f,0x92d8,0xc0c9,
0x9ddd,0xa5fc,0x99ec,0x4131,
0x9176,0xd329,0x1fea,0xc171
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x6572,0xeeb3,0xb8a5,0x40ca,
0xbc96,0x582a,0x27bc,0xc134,
0xd8ef,0xc2ea,0xd98f,0x4177,
0x5a31,0x3cbe,0xafe0,0xc189
};
static unsigned short P1[] = {0x9d39,0xa252,0xdf46,0x3f91};
#define PI180 *(double *)P1
static double lossth = 1.0e14;
#endif

#ifdef MIEEE
static unsigned short P[] = {
0xc0c9,0x92d8,0xd24f,0x3f38,
0x4131,0x99ec,0xa5fc,0x9ddd,
0xc171,0x1fea,0xd329,0x9176
};
static unsigned short Q[] = {
0x40ca,0xb8a5,0xeeb3,0x6572,
0xc134,0x27bc,0x582a,0xbc96,
0x4177,0xd98f,0xc2ea,0xd8ef,
0xc189,0xafe0,0x3cbe,0x5a31
};
static unsigned short P1[] = {
0x3f91,0xdf46,0xa252,0x9d39
};
#define PI180 *(double *)P1
static double lossth = 1.0e14;
#endif

#ifndef ANSIPROT
double polevl(), p1evl(), floor(), ldexp();
static double tancot();
#else
static double tancot(double, int);
#endif
extern double MAXNUM;
extern double PIO4;


double tandg(x)
double x;
{

return( tancot(x,0) );
}


double cotdg(x)
double x;
{

return( tancot(x,1) );
}


static double tancot( xx, cotflg )
double xx;
int cotflg;
{
double x, y, z, zz;
int j, sign;

/* make argument positive but save the sign */
if( xx < 0 )
	{
	x = -xx;
	sign = -1;
	}
else
	{
	x = xx;
	sign = 1;
	}

if( x > lossth )
	{
	mtherr( "tandg", TLOSS );
	return(0.0);
	}

/* compute x mod PIO4 */
y = floor( x/45.0 );

/* strip high bits of integer part */
z = ldexp( y, -3 );
z = floor(z);		/* integer part of y/8 */
z = y - ldexp( z, 3 );  /* y - 16 * (y/16) */

/* integer and fractional part modulo one octant */
j = z;

/* map zeros and singularities to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}

z = x - y * 45.0;
z *= PI180;

zz = z * z;

if( zz > 1.0e-14 )
	y = z  +  z * (zz * polevl( zz, P, 2 )/p1evl(zz, Q, 4));
else
	y = z;
	
if( j & 2 )
	{
	if( cotflg )
		y = -y;
	else
		{
		if( y != 0.0 )
			{
			y = -1.0/y;
			}
		else
			{
			mtherr( "tandg", SING );
			y = MAXNUM;
			}
		}
	}
else
	{
	if( cotflg )
		{
		if( y != 0.0 )
			y = 1.0/y;
		else
			{
			mtherr( "cotdg", SING );
			y = MAXNUM;
			}
		}
	}

if( sign < 0 )
	y = -y;

return( y );
}
