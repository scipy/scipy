/*							log2.c
 *
 *	Base 2 logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, log2();
 *
 * y = log2( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base 2 logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the base e
 * logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)/Q(x).
 *
 * Otherwise, setting  z = 2(x-1)/x+1),
 * 
 *     log(x) = z + z**3 P(z)/Q(z).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0.5, 2.0    30000       2.0e-16     5.5e-17
 *    IEEE      exp(+-700)  40000       1.3e-16     4.6e-17
 *
 * In the tests over the interval [exp(+-700)], the logarithms
 * of the random arguments were uniformly distributed.
 *
 * ERROR MESSAGES:
 *
 * log2 singularity:  x = 0; returns -INFINITY
 * log2 domain:       x < 0; returns NAN
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include "mconf.h"
static char fname[] = {"log2"};

/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 */
#ifdef UNK
static double P[] = {
 1.01875663804580931796E-4,
 4.97494994976747001425E-1,
 4.70579119878881725854E0,
 1.44989225341610930846E1,
 1.79368678507819816313E1,
 7.70838733755885391666E0,
};
static double Q[] = {
/* 1.00000000000000000000E0, */
 1.12873587189167450590E1,
 4.52279145837532221105E1,
 8.29875266912776603211E1,
 7.11544750618563894466E1,
 2.31251620126765340583E1,
};
#define LOG2EA 0.44269504088896340735992
#endif

#ifdef DEC
static unsigned short P[] = {
0037777,0127270,0162547,0057274,
0041001,0054665,0164317,0005341,
0041451,0034104,0031640,0105773,
0041677,0011276,0123617,0160135,
0041701,0126603,0053215,0117250,
0041420,0115777,0135206,0030232,
};
static unsigned short Q[] = {
/*0040200,0000000,0000000,0000000,*/
0041220,0144332,0045272,0174241,
0041742,0164566,0035720,0130431,
0042246,0126327,0166065,0116357,
0042372,0033420,0157525,0124560,
0042271,0167002,0066537,0172303,
0041730,0164777,0113711,0044407,
};
static unsigned short L[5] = {0037742,0124354,0122560,0057703};
#define LOG2EA (*(double *)(&L[0]))
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x1bb0,0x93c3,0xb4c2,0x3f1a,
0x52f2,0x3f56,0xd6f5,0x3fdf,
0x6911,0xed92,0xd2ba,0x4012,
0xeb2e,0xc63e,0xff72,0x402c,
0xc84d,0x924b,0xefd6,0x4031,
0xdcf8,0x7d7e,0xd563,0x401e,
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xef8e,0xae97,0x9320,0x4026,
0xc033,0x4e19,0x9d2c,0x4046,
0xbdbd,0xa326,0xbf33,0x4054,
0xae21,0xeb5e,0xc9e2,0x4051,
0x25b2,0x9e1f,0x200a,0x4037,
};
static unsigned short L[5] = {0x0bf8,0x94ae,0x551d,0x3fdc};
#define LOG2EA (*(double *)(&L[0]))
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3f1a,0xb4c2,0x93c3,0x1bb0,
0x3fdf,0xd6f5,0x3f56,0x52f2,
0x4012,0xd2ba,0xed92,0x6911,
0x402c,0xff72,0xc63e,0xeb2e,
0x4031,0xefd6,0x924b,0xc84d,
0x401e,0xd563,0x7d7e,0xdcf8,
};
static unsigned short Q[] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4026,0x9320,0xae97,0xef8e,
0x4046,0x9d2c,0x4e19,0xc033,
0x4054,0xbf33,0xa326,0xbdbd,
0x4051,0xc9e2,0xeb5e,0xae21,
0x4037,0x200a,0x9e1f,0x25b2,
};
static unsigned short L[5] = {0x3fdc,0x551d,0x94ae,0x0bf8};
#define LOG2EA (*(double *)(&L[0]))
#endif

/* Coefficients for log(x) = z + z**3 P(z)/Q(z),
 * where z = 2(x-1)/(x+1)
 * 1/sqrt(2) <= x < sqrt(2)
 */

#ifdef UNK
static double R[3] = {
-7.89580278884799154124E-1,
 1.63866645699558079767E1,
-6.41409952958715622951E1,
};
static double S[3] = {
/* 1.00000000000000000000E0,*/
-3.56722798256324312549E1,
 3.12093766372244180303E2,
-7.69691943550460008604E2,
};
/* log2(e) - 1 */
#define LOG2EA 0.44269504088896340735992
#endif
#ifdef DEC
static unsigned short R[12] = {
0140112,0020756,0161540,0072035,
0041203,0013743,0114023,0155527,
0141600,0044060,0104421,0050400,
};
static unsigned short S[12] = {
/*0040200,0000000,0000000,0000000,*/
0141416,0130152,0017543,0064122,
0042234,0006000,0104527,0020155,
0142500,0066110,0146631,0174731,
};
/* log2(e) - 1 */
#define LOG2EA 0.44269504088896340735992L
#endif
#ifdef IBMPC
static unsigned short R[12] = {
0x0e84,0xdc6c,0x443d,0xbfe9,
0x7b6b,0x7302,0x62fc,0x4030,
0x2a20,0x1122,0x0906,0xc050,
};
static unsigned short S[12] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x6d0a,0x43ec,0xd60d,0xc041,
0xe40e,0x112a,0x8180,0x4073,
0x3f3b,0x19b3,0x0d89,0xc088,
};
#endif
#ifdef MIEEE
static unsigned short R[12] = {
0xbfe9,0x443d,0xdc6c,0x0e84,
0x4030,0x62fc,0x7302,0x7b6b,
0xc050,0x0906,0x1122,0x2a20,
};
static unsigned short S[12] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0xc041,0xd60d,0x43ec,0x6d0a,
0x4073,0x8180,0x112a,0xe40e,
0xc088,0x0d89,0x19b3,0x3f3b,
};
#endif

#ifndef ANSIPROT
double frexp(), ldexp(), polevl(), p1evl();
int isnan(), isfinite();
#endif
#define SQRTH 0.70710678118654752440
extern double LOGE2, INFINITY, NAN;

double log2(x)
double x;
{
int e;
double y;
VOLATILE double z;
#ifdef DEC
short *q;
#endif

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

/* Note, frexp is used so that denormal numbers
 * will be handled properly.
 */
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


/* logarithm using log(x) = z + z**3 P(z)/Q(z),
 * where z = 2(x-1)/x+1)
 */

if( (e > 2) || (e < -2) )
{
if( x < SQRTH )
	{ /* 2( 2x-1 )/( 2x+1 ) */
	e -= 1;
	z = x - 0.5;
	y = 0.5 * z + 0.5;
	}	
else
	{ /*  2 (x-1)/(x+1)   */
	z = x - 0.5;
	z -= 0.5;
	y = 0.5 * x  + 0.5;
	}

x = z / y;
z = x*x;
y = x * ( z * polevl( z, R, 2 ) / p1evl( z, S, 3 ) );
goto ldone;
}



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

z = x*x;
#if DEC
y = x * ( z * polevl( x, P, 5 ) / p1evl( x, Q, 6 ) ) - ldexp( z, -1 );
#else
y = x * ( z * polevl( x, P, 5 ) / p1evl( x, Q, 5 ) ) - ldexp( z, -1 );
#endif

ldone:

/* Multiply log of fraction by log2(e)
 * and base 2 exponent by 1
 *
 * ***CAUTION***
 *
 * This sequence of operations is critical and it may
 * be horribly defeated by some compiler optimizers.
 */
z = y * LOG2EA;
z += x * LOG2EA;
z += y;
z += x;
z += e;
return( z );
}
