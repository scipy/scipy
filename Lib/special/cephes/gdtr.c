/*							gdtr.c
 *
 *	Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, gdtr();
 *
 * y = gdtr( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from zero to x of the Gamma probability
 * density function:
 *
 *
 *                x
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               0
 *
 *  The incomplete Gamma integral is used, according to the
 * relation
 *
 * y = igam( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See igam().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtr domain         x < 0            0.0
 *
 */
/*							gdtrc.c
 *
 *	Complemented Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, b, x, y, gdtrc();
 *
 * y = gdtrc( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from x to infinity of the Gamma
 * probability density function:
 *
 *
 *               inf.
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               x
 *
 *  The incomplete Gamma integral is used, according to the
 * relation
 *
 * y = igamc( b, ax ).
 *
 *
 * ACCURACY:
 *
 * See igamc().
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtrc domain         x < 0            0.0
 *
 */

/*							gdtr()  */


/*
Cephes Math Library Release 2.3:  March,1995
Copyright 1984, 1987, 1995 by Stephen L. Moshier
*/

#include "mconf.h"
#ifndef ANSIPROT
double igam(), igamc();
#else
double gdtri(double,double,double);
#endif

extern double NAN;

double gdtr( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtr", DOMAIN );
	return( NAN );
	}
return(  igam( b, a * x )  );
}


double gdtrc( a, b, x )
double a, b, x;
{

if( x < 0.0 )
	{
	mtherr( "gdtrc", DOMAIN );
	return( NAN );
	}
return(  igamc( b, a * x )  );
}


double gdtri( a, b, y)
double a, b, y;
{

if ((y < 0.0) || (y > 1.0) || (a <= 0.0) || (b < 0.0))
  {
    mtherr("gdtri", DOMAIN);
    return( NAN );
  }

return ( igami (b, 1.0-y) / a);
}
