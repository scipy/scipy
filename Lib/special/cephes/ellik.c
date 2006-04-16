/*							ellik.c
 *
 *	Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double phi, m, y, ellik();
 *
 * y = ellik( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi | m) =   |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with m in [0, 1] and phi as indicated.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -10,10       200000      7.4e-16     1.0e-16
 *
 *
 */


/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/*	Incomplete elliptic integral of first kind	*/

#include "mconf.h"
#ifndef ANSIPROT
double sqrt(), fabs(), log(), tan(), atan(), floor(), ellpk();
#endif
extern double PI, PIO2, MACHEP, MAXNUM;

double ellik( phi, m )
double phi, m;
{
double a, b, c, e, temp, t, K;
int d, mod, sign, npio2;

if( m == 0.0 )
	return( phi );
a = 1.0 - m;
if( a == 0.0 )
	{
	if( fabs(phi) >= PIO2 )
		{
		mtherr( "ellik", SING );
		return( MAXNUM );
		}
	return(  log(  tan( (PIO2 + phi)/2.0 )  )   );
	}
npio2 = floor( phi/PIO2 );
if( npio2 & 1 )
	npio2 += 1;
if( npio2 )
	{
	K = ellpk( m );  /* Changed */
	phi = phi - npio2 * PIO2;
	}
else
	K = 0.0;
if( phi < 0.0 )
	{
	phi = -phi;
	sign = -1;
	}
else
	sign = 0;
b = sqrt(a);
t = tan( phi );
if( fabs(t) > 10.0 )
	{
	/* Transform the amplitude */
	e = 1.0/(b*t);
	/* ... but avoid multiple recursions.  */
	if( fabs(e) < 10.0 )
		{
		e = atan(e);
		if( npio2 == 0 )
		    K = ellpk( m ); /* Changed */
		temp = K - ellik( e, m );
		goto done;
		}
	}
a = 1.0;
c = sqrt(m);
d = 1;
mod = 0;

while( fabs(c/a) > MACHEP )
	{
	temp = b/a;
	phi = phi + atan(t*temp) + mod * PI;
	mod = (phi + PIO2)/PI;
	t = t * ( 1.0 + temp )/( 1.0 - temp * t * t );
	c = ( a - b )/2.0;
	temp = sqrt( a * b );
	a = ( a + b )/2.0;
	b = temp;
	d += d;
	}

temp = (atan(t) + mod * PI)/(d * a);

done:
if( sign < 0 )
	temp = -temp;
temp += npio2 * K;
return( temp );
}
