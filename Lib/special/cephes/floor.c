/*							ceil()
 *							floor()
 *							frexp()
 *							ldexp()
 *							signbit()
 *							isnan()
 *							isfinite()
 *
 *	Floating point numeric utilities
 *
 *
 *
 * SYNOPSIS:
 *
 * double ceil(), floor(), frexp(), ldexp();
 * int signbit(), isnan(), isfinite();
 * double x, y;
 * int expnt, n;
 *
 * y = floor(x);
 * y = ceil(x);
 * y = frexp( x, &expnt );
 * y = ldexp( x, n );
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


#include "mconf.h"

#ifdef UNK
/* ceil(), floor(), frexp(), ldexp() may need to be rewritten. */
#undef UNK
#if BIGENDIAN
#define MIEEE 1
#else
#define IBMPC 1
#endif
#endif

#ifdef DEC
#define EXPMSK 0x807f
#define MEXP 255
#define NBITS 56
#endif

#ifdef IBMPC
#define EXPMSK 0x800f
#define MEXP 0x7ff
#define NBITS 53
#endif

#ifdef MIEEE
#define EXPMSK 0x800f
#define MEXP 0x7ff
#define NBITS 53
#endif

extern double MAXNUM, NEGZERO;
#ifndef ANSIPROT
double floor();
#endif

double ceil(x)
double x;
{
double y;

#ifdef UNK
mtherr( "ceil", DOMAIN );
return(0.0);
#endif
#ifdef NANS
if( isnan(x) )
	return( x );
#endif
#ifdef INFINITIES
if(!isfinite(x))
	return(x);
#endif

y = floor(x);
if( y < x )
	y += 1.0;
#ifdef MINUSZERO
if( y == 0.0 && x < 0.0 )
	return( NEGZERO );
#endif
return(y);
}




/* Bit clearing masks: */

static unsigned short bmask[] = {
0xffff,
0xfffe,
0xfffc,
0xfff8,
0xfff0,
0xffe0,
0xffc0,
0xff80,
0xff00,
0xfe00,
0xfc00,
0xf800,
0xf000,
0xe000,
0xc000,
0x8000,
0x0000,
};





double floor(x)
double x;
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
unsigned short *p;
int e;

#ifdef UNK
mtherr( "floor", DOMAIN );
return(0.0);
#endif
#ifdef NANS
if( isnan(x) )
	return( x );
#endif
#ifdef INFINITIES
if(!isfinite(x))
	return(x);
#endif
#ifdef MINUSZERO
if(x == 0.0L)
	return(x);
#endif
u.y = x;
/* find the exponent (power of 2) */
#ifdef DEC
p = (unsigned short *)&u.sh[0];
e = (( *p  >> 7) & 0377) - 0201;
p += 3;
#endif

#ifdef IBMPC
p = (unsigned short *)&u.sh[3];
e = (( *p >> 4) & 0x7ff) - 0x3ff;
p -= 3;
#endif

#ifdef MIEEE
p = (unsigned short *)&u.sh[0];
e = (( *p >> 4) & 0x7ff) - 0x3ff;
p += 3;
#endif

if( e < 0 )
	{
	if( u.y < 0.0 )
		return( -1.0 );
	else
		return( 0.0 );
	}

e = (NBITS -1) - e;
/* clean out 16 bits at a time */
while( e >= 16 )
	{
#ifdef IBMPC
	*p++ = 0;
#endif

#ifdef DEC
	*p-- = 0;
#endif

#ifdef MIEEE
	*p-- = 0;
#endif
	e -= 16;
	}

/* clear the remaining bits */
if( e > 0 )
	*p &= bmask[e];

if( (x < 0) && (u.y != x) )
	u.y -= 1.0;

return(u.y);
}




double frexp( x, pw2 )
double x;
int *pw2;
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
int i;
#ifdef DENORMAL
int k;
#endif
short *q;

u.y = x;

#ifdef UNK
mtherr( "frexp", DOMAIN );
return(0.0);
#endif

#ifdef IBMPC
q = (short *)&u.sh[3];
#endif

#ifdef DEC
q = (short *)&u.sh[0];
#endif

#ifdef MIEEE
q = (short *)&u.sh[0];
#endif

/* find the exponent (power of 2) */
#ifdef DEC
i  = ( *q >> 7) & 0377;
if( i == 0 )
	{
	*pw2 = 0;
	return(0.0);
	}
i -= 0200;
*pw2 = i;
*q &= 0x807f;	/* strip all exponent bits */
*q |= 040000;	/* mantissa between 0.5 and 1 */
return(u.y);
#endif

#ifdef IBMPC
i  = ( *q >> 4) & 0x7ff;
if( i != 0 )
	goto ieeedon;
#endif

#ifdef MIEEE
i  =  *q >> 4;
i &= 0x7ff;
if( i != 0 )
	goto ieeedon;
#ifdef DENORMAL

#else
*pw2 = 0;
return(0.0);
#endif

#endif


#ifndef DEC
/* Number is denormal or zero */
#ifdef DENORMAL
if( u.y == 0.0 )
	{
	*pw2 = 0;
	return( 0.0 );
	}


/* Handle denormal number. */
do
	{
	u.y *= 2.0;
	i -= 1;
	k  = ( *q >> 4) & 0x7ff;
	}
while( k == 0 );
i = i + k;
#endif /* DENORMAL */

ieeedon:

i -= 0x3fe;
*pw2 = i;
*q &= 0x800f;
*q |= 0x3fe0;
return( u.y );
#endif
}







double ldexp( x, pw2 )
double x;
int pw2;
{
union
	{
	double y;
	unsigned short sh[4];
	} u;
short *q;
int e;

#ifdef UNK
mtherr( "ldexp", DOMAIN );
return(0.0);
#endif

u.y = x;
#ifdef DEC
q = (short *)&u.sh[0];
e  = ( *q >> 7) & 0377;
if( e == 0 )
	return(0.0);
#else

#ifdef IBMPC
q = (short *)&u.sh[3];
#endif
#ifdef MIEEE
q = (short *)&u.sh[0];
#endif
while( (e = (*q & 0x7ff0) >> 4) == 0 )
	{
	if( u.y == 0.0 )
		{
		return( 0.0 );
		}
/* Input is denormal. */
	if( pw2 > 0 )
		{
		u.y *= 2.0;
		pw2 -= 1;
		}
	if( pw2 < 0 )
		{
		if( pw2 < -53 )
			return(0.0);
		u.y /= 2.0;
		pw2 += 1;
		}
	if( pw2 == 0 )
		return(u.y);
	}
#endif /* not DEC */

e += pw2;

/* Handle overflow */
#ifdef DEC
if( e > MEXP )
	return( MAXNUM );
#else
if( e >= MEXP )
	return( 2.0*MAXNUM );
#endif

/* Handle denormalized results */
if( e < 1 )
	{
#ifdef DENORMAL
	if( e < -53 )
		return(0.0);
	*q &= 0x800f;
	*q |= 0x10;
	/* For denormals, significant bits may be lost even
	   when dividing by 2.  Construct 2^-(1-e) so the result
	   is obtained with only one multiplication.  */
	u.y *= ldexp(1.0, e-1);
	return(u.y);
#else
	return(0.0);
#endif
	}
else
	{
#ifdef DEC
	*q &= 0x807f;	/* strip all exponent bits */
	*q |= (e & 0xff) << 7;
#else
	*q &= 0x800f;
	*q |= (e & 0x7ff) << 4;
#endif
	return(u.y);
	}
}
