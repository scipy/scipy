/*							sqrt.c
 *
 *	Square root
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, sqrt();
 *
 * y = sqrt( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the square root of x.
 *
 * Range reduction involves isolating the power of two of the
 * argument and using a polynomial approximation to obtain
 * a rough value for the square root.  Then Heron's iteration
 * is used three times to converge to an accurate value.
 *
 *
 *
 * ACCURACY:
 *
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 10       60000       2.1e-17     7.9e-18
 *    IEEE      0,1.7e308   30000       1.7e-16     6.3e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * sqrt domain        x < 0            0.0
 *
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


#include "mconf.h"
#ifndef ANSIPROT
double frexp(), ldexp();
#endif
extern double SQRT2;  /*  SQRT2 = 1.41421356237309504880 */

double sqrt(x)
double x;
{
int e;
#ifndef UNK
short *q;
#endif
double z, w;

if( x <= 0.0 )
	{
	if( x < 0.0 )
		mtherr( "sqrt", DOMAIN );
	return( 0.0 );
	}
w = x;
/* separate exponent and significand */
#ifdef UNK
z = frexp( x, &e );
#endif
#ifdef DEC
q = (short *)&x;
e = ((*q >> 7) & 0377) - 0200;
*q &= 0177;
*q |= 040000;
z = x;
#endif

/* Note, frexp and ldexp are used in order to
 * handle denormal numbers properly.
 */
#ifdef IBMPC
z = frexp( x, &e );
q = (short *)&x;
q += 3;
/*
e = ((*q >> 4) & 0x0fff) - 0x3fe;
*q &= 0x000f;
*q |= 0x3fe0;
z = x;
*/
#endif
#ifdef MIEEE
z = frexp( x, &e );
q = (short *)&x;
/*
e = ((*q >> 4) & 0x0fff) - 0x3fe;
*q &= 0x000f;
*q |= 0x3fe0;
z = x;
*/
#endif

/* approximate square root of number between 0.5 and 1
 * relative error of approximation = 7.47e-3
 */
x = 4.173075996388649989089E-1 + 5.9016206709064458299663E-1 * z;

/* adjust for odd powers of 2 */
if( (e & 1) != 0 )
	x *= SQRT2;

/* re-insert exponent */
#ifdef UNK
x = ldexp( x, (e >> 1) );
#endif
#ifdef DEC
*q += ((e >> 1) & 0377) << 7;
*q &= 077777;
#endif
#ifdef IBMPC
x = ldexp( x, (e >> 1) );
/*
*q += ((e >>1) & 0x7ff) << 4;
*q &= 077777;
*/
#endif
#ifdef MIEEE
x = ldexp( x, (e >> 1) );
/*
*q += ((e >>1) & 0x7ff) << 4;
*q &= 077777;
*/
#endif

/* Newton iterations: */
#ifdef UNK
x = 0.5*(x + w/x);
x = 0.5*(x + w/x);
x = 0.5*(x + w/x);
#endif

/* Note, assume the square root cannot be denormal,
 * so it is safe to use integer exponent operations here.
 */
#ifdef DEC
x += w/x;
*q -= 0200;
x += w/x;
*q -= 0200;
x += w/x;
*q -= 0200;
#endif
#ifdef IBMPC
x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
#endif
#ifdef MIEEE
x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
x += w/x;
*q -= 0x10;
#endif

return(x);
}
