/*							igami()
 *
 *      Inverse of complemented imcomplete Gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, p, igami();
 *
 * x = igami( a, p );
 *
 * DESCRIPTION:
 *
 * Given p, the function finds x such that
 *
 *  igamc( a, x ) = p.
 *
 * Starting with the approximate value
 *
 *         3
 *  x = a t
 *
 *  where
 *
 *  t = 1 - d - ndtri(p) sqrt(d)
 * 
 * and
 *
 *  d = 1/9a,
 *
 * the routine performs up to 10 Newton iterations to find the
 * root of igamc(a,x) - p = 0.
 *
 * ACCURACY:
 *
 * Tested at random a, p in the intervals indicated.
 *
 *                a        p                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,0.5       100000       1.0e-14     1.7e-15
 *    IEEE     0.01,0.5  0,0.5       100000       9.0e-14     3.4e-15
 *    IEEE    0.5,10000  0,0.5        20000       2.3e-13     3.8e-14
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

#include "mconf.h"
#include <stdio.h>

extern double MACHEP, MAXNUM, MAXLOG, MINLOG;

double igami( a, y0 )
double a, y0;
{
double x0, x1, x, yl, yh, y, d, lgm, dithresh;
int i, dir;

/* bound the solution */
x0 = MAXNUM;
yl = 0;
x1 = 0;
yh = 1.0;
dithresh = 5.0 * MACHEP;

if ((y0<0.0) || (y0>1.0) || (a<=0)) {
   mtherr("igami", DOMAIN);
   return(NPY_NAN);
}

if (y0==0.0) {
  return(MAXNUM);
}

if (y0==1.0){
   return 0.0;
}

/* approximation to inverse function */
d = 1.0/(9.0*a);
y = ( 1.0 - d - ndtri(y0) * sqrt(d) );
x = a * y * y * y;

lgm = lgam(a);

for( i=0; i<10; i++ )
	{
	if( x > x0 || x < x1 )
		goto ihalve;
	y = igamc(a,x);
	if( y < yl || y > yh )
		goto ihalve;
	if( y < y0 )
		{
		x0 = x;
		yl = y;
		}
	else
		{
		x1 = x;
		yh = y;
		}
/* compute the derivative of the function at this point */
	d = (a - 1.0) * log(x) - x - lgm;
	if( d < -MAXLOG )
		goto ihalve;
	d = -exp(d);
/* compute the step to the next approximation of x */
	d = (y - y0)/d;
	if( fabs(d/x) < MACHEP )
		goto done;
	x = x - d;
	}

/* Resort to interval halving if Newton iteration did not converge. */
ihalve:

d = 0.0625;
if( x0 == MAXNUM )
	{
	if( x <= 0.0 )
		x = 1.0;
	while( x0 == MAXNUM )
		{
		x = (1.0 + d) * x;
		y = igamc( a, x );
		if( y < y0 )
			{
			x0 = x;
			yl = y;
			break;
			}
		d = d + d;
		}
	}
d = 0.5;
dir = 0;

for( i=0; i<400; i++ )
	{
	x = x1  +  d * (x0 - x1);
	y = igamc( a, x );
	lgm = (x0 - x1)/(x1 + x0);
	if( fabs(lgm) < dithresh )
		break;
	lgm = (y - y0)/y0;
	if( fabs(lgm) < dithresh )
		break;
	if( x <= 0.0 )
		break;
	if( y >= y0 )
		{
		x1 = x;
		yh = y;
		if( dir < 0 )
			{
			dir = 0;
			d = 0.5;
			}
		else if( dir > 1 )
			d = 0.5 * d + 0.5; 
		else
			d = (y0 - yl)/(yh - yl);
		dir += 1;
		}
	else
		{
		x0 = x;
		yl = y;
		if( dir > 0 )
			{
			dir = 0;
			d = 0.5;
			}
		else if( dir < -1 )
			d = 0.5 * d;
		else
			d = (y0 - yl)/(yh - yl);
		dir -= 1;
		}
	}
if( x == 0.0 )
	mtherr( "igami", UNDERFLOW );

done:
return( x );
}
