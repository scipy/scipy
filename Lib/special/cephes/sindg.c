/*							sindg.c
 *
 *	Circular sine of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, sindg();
 *
 * y = sindg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of 45 degrees.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *      x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *      1  -  x**2 P(x**2).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       +-1000        3100      3.3e-17      9.0e-18
 *    IEEE      +-1000       30000      2.3e-16      5.6e-17
 * 
 * ERROR MESSAGES:
 *
 *   message           condition        value returned
 * sindg total loss   x > 8.0e14 (DEC)      0.0
 *                    x > 1.0e14 (IEEE)
 *
 */
/*							cosdg.c
 *
 *	Circular cosine of angle in degrees
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, cosdg();
 *
 * y = cosdg( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of 45 degrees.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the cosine is approximated by
 *      1  -  x**2 P(x**2).
 * Between pi/4 and pi/2 the sine is represented as
 *      x  +  x**3 P(x**2).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC      +-1000         3400       3.5e-17     9.1e-18
 *    IEEE     +-1000        30000       2.1e-16     5.7e-17
 *  See also sin().
 *
 */

/* Cephes Math Library Release 2.0:  April, 1987
 * Copyright 1985, 1987 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

#include "mconf.h"

#ifdef UNK
static double sincof[] = {
 1.58962301572218447952E-10,
-2.50507477628503540135E-8,
 2.75573136213856773549E-6,
-1.98412698295895384658E-4,
 8.33333333332211858862E-3,
-1.66666666666666307295E-1
};
static double coscof[] = {
 1.13678171382044553091E-11,
-2.08758833757683644217E-9,
 2.75573155429816611547E-7,
-2.48015872936186303776E-5,
 1.38888888888806666760E-3,
-4.16666666666666348141E-2,
 4.99999999999999999798E-1
};
static double PI180 = 1.74532925199432957692E-2; /* pi/180 */
static double lossth = 1.0e14;
#endif

#ifdef DEC
static unsigned short sincof[] = {
0030056,0143750,0177170,0073013,
0131727,0027455,0044510,0132205,
0033470,0167432,0131752,0042263,
0135120,0006400,0146776,0174027,
0036410,0104210,0104207,0137202,
0137452,0125252,0125252,0125103
};
static unsigned short coscof[] = {
0027107,0176030,0153315,0110312,
0131017,0072476,0007450,0123243,
0032623,0171174,0070066,0146445,
0134320,0006400,0147355,0163313,
0035666,0005540,0133012,0165067,
0137052,0125252,0125252,0125206,
0040000,0000000,0000000,0000000
};
static unsigned short P1[] = {0036616,0175065,0011224,0164711};
#define PI180 *(double *)P1
static double lossth = 8.0e14;
#endif

#ifdef IBMPC
static unsigned short sincof[] = {
0x0ec1,0x1fcf,0xd8fd,0x3de5,
0x1691,0xa929,0xe5e5,0xbe5a,
0x4896,0x567d,0x1de3,0x3ec7,
0xdf03,0x19bf,0x01a0,0xbf2a,
0xf7d0,0x1110,0x1111,0x3f81,
0x5548,0x5555,0x5555,0xbfc5
};
static unsigned short coscof[] = {
0xb219,0x1ad9,0xff83,0x3da8,
0x14d4,0xc1e5,0xeea7,0xbe21,
0xd9a5,0x8e06,0x7e4f,0x3e92,
0xbcd9,0x19dd,0x01a0,0xbefa,
0x5d47,0x16c1,0xc16c,0x3f56,
0x5551,0x5555,0x5555,0xbfa5,
0x0000,0x0000,0x0000,0x3fe0
};

static unsigned short P1[] = {0x9d39,0xa252,0xdf46,0x3f91};
#define PI180 *(double *)P1
static double lossth = 1.0e14;
#endif

#ifdef MIEEE
static unsigned short sincof[] = {
0x3de5,0xd8fd,0x1fcf,0x0ec1,
0xbe5a,0xe5e5,0xa929,0x1691,
0x3ec7,0x1de3,0x567d,0x4896,
0xbf2a,0x01a0,0x19bf,0xdf03,
0x3f81,0x1111,0x1110,0xf7d0,
0xbfc5,0x5555,0x5555,0x5548
};
static unsigned short coscof[] = {
0x3da8,0xff83,0x1ad9,0xb219,
0xbe21,0xeea7,0xc1e5,0x14d4,
0x3e92,0x7e4f,0x8e06,0xd9a5,
0xbefa,0x01a0,0x19dd,0xbcd9,
0x3f56,0xc16c,0x16c1,0x5d47,
0xbfa5,0x5555,0x5555,0x5551,
0x3fe0,0x0000,0x0000,0x0000
};

static unsigned short P1[] = {
0x3f91,0xdf46,0xa252,0x9d39
};
#define PI180 *(double *)P1
static double lossth = 1.0e14;
#endif

#ifndef ANSIPROT
double polevl(), floor(), ldexp();
#else
extern double polevl (double, void *, int);
extern double floor(double);
extern double ldexp(double,int);
#endif
extern double PIO4;

double sindg(x)
double x;
{
double y, z, zz;
int j, sign;

/* make argument positive but save the sign */
sign = 1;
if( x < 0 )
	{
	x = -x;
	sign = -1;
	}

if( x > lossth )
	{
	mtherr( "sindg", TLOSS );
	return(0.0);
	}

y = floor( x/45.0 ); /* integer part of x/PIO4 */

/* strip high bits of integer part to prevent integer overflow */
z = ldexp( y, -4 );
z = floor(z);           /* integer part of y/8 */
z = y - ldexp( z, 4 );  /* y - 16 * (y/16) */

j = z; /* convert to integer for tests on the phase angle */
/* map zeros to origin */
if( j & 1 )
	{
	j += 1;
	y += 1.0;
	}
j = j & 07; /* octant modulo 360 degrees */
/* reflect in x axis */
if( j > 3)
	{
	sign = -sign;
	j -= 4;
	}

z = x - y * 45.0; /* x mod 45 degrees */
z *= PI180;	/* multiply by pi/180 to convert to radians */
zz = z * z;

if( (j==1) || (j==2) )
	{
	y = 1.0 - zz * polevl( zz, coscof, 6 );
	}
else
	{
	y = z  +  z * (zz * polevl( zz, sincof, 5 ));
	}

if(sign < 0)
	y = -y;

return(y);
}


double cosdg(x)
double x;
{
double y, z, zz;
int j, sign;

/* make argument positive */
sign = 1;
if( x < 0 )
	x = -x;

if( x > lossth )
	{
	mtherr( "cosdg", TLOSS );
	return(0.0);
	}

y = floor( x/45.0 );
z = ldexp( y, -4 );
z = floor(z);		/* integer part of y/8 */
z = y - ldexp( z, 4 );  /* y - 16 * (y/16) */

/* integer and fractional part modulo one octant */
j = z;
if( j & 1 )	/* map zeros to origin */
	{
	j += 1;
	y += 1.0;
	}
j = j & 07;
if( j > 3)
	{
	j -=4;
	sign = -sign;
	}

if( j > 1 )
	sign = -sign;

z = x - y * 45.0; /* x mod 45 degrees */
z *= PI180;	/* multiply by pi/180 to convert to radians */

zz = z * z;

if( (j==1) || (j==2) )
	{
	y = z  +  z * (zz * polevl( zz, sincof, 5 ));
	}
else
	{
	y = 1.0 - zz * polevl( zz, coscof, 6 );
	}

if(sign < 0)
	y = -y;

return(y);
}


/* Degrees, minutes, seconds to radians: */

/* 1 arc second, in radians = 4.848136811095359935899141023579479759563533023727e-6 */
static double P64800 = 4.848136811095359935899141023579479759563533023727e-6;

double radian(d,m,s)
double d,m,s;
{
return( ((d*60.0 + m)*60.0 + s)*P64800 );
}



