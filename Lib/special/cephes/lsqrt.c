/*							lsqrt.c
 *
 *	Integer square root
 *
 *
 *
 * SYNOPSIS:
 *
 * long x, y;
 * long lsqrt();
 *
 * y = lsqrt( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns a long integer square root of the long integer
 * argument.  The computation is by binary long division.
 *
 * The largest possible result is lsqrt(2,147,483,647)
 * = 46341.
 *
 * If x < 0, the square root of |x| is returned, and an
 * error message is printed.
 *
 *
 * ACCURACY:
 *
 * An extra, roundoff, bit is computed; hence the result
 * is the nearest integer to the actual square root.
 * NOTE: only DEC arithmetic is currently supported.
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"

long lsqrt(x)
long x;
{
long pq, num, sq;
long temp;
int i, j, k, n;

if( x < 0 )
	{
	mtherr( "lsqrt", DOMAIN );
	x = -x;
	}

num = 0;
sq = 0;
k = 24;
n = 4;

for( j=0; j<4; j++ )
	{
	num |= (x >> k) & 0xff;	/* bring in next byte of arg */
	if( j == 3 )		/* do roundoff bit at end */
		n = 5;
	for( i=0; i<n; i++ )
		{
		num <<= 2;		/* next 2 bits of arg */
		sq <<= 1;		/* shift up answer */
		temp = (sq << 1) + 256;	/* trial divisor */
		temp = num - temp;
		if( temp >= 0 )
			{
			num = temp;	/* it went in */
			sq += 256;	/* answer bit = 1 */
			}
		}
	k -= 8;	/* shift count to get next byte of arg */
	}

sq += 256;	/* add roundoff bit */
sq >>= 9;	/* truncate */
return( sq );
}
