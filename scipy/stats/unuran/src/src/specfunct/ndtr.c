/*							ndtr.c
 *
 *	Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, ndtr();
 *
 * y = ndtr( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x:
 *
 *                            x
 *                             -
 *                   1        | |          2
 *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
 *               sqrt(2pi)  | |
 *                           -
 *                          -inf.
 *
 *             =  ( 1 + erf(z) ) / 2
 *             =  erfc(z) / 2
 *
 * where z = x/sqrt(2). Computation is via the functions
 * erf and erfc.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC      -13,0         8000       2.1e-15     4.8e-16
 *    IEEE     -13,0        30000       3.4e-14     6.7e-15
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition         value returned
 * erfc underflow    x > 37.519379347       0.0
 *
 */
/*							erf.c
 *
 *	Error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erf();
 *
 * y = erf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The integral is
 *
 *                           x 
 *                            -
 *                 2         | |          2
 *   erf(x)  =  --------     |    exp( - t  ) dt.
 *              sqrt(pi)   | |
 *                          -
 *                           0
 *
 * The magnitude of x is limited to 9.231948545 for DEC
 * arithmetic; 1 or -1 is returned outside this range.
 *
 * For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
 * erf(x) = 1 - erfc(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,1         14000       4.7e-17     1.5e-17
 *    IEEE      0,1         30000       3.7e-16     1.0e-16
 *
 */
/*							erfc.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erfc();
 *
 * y = erfc( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf. 
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 *
 * For small x, erfc(x) = 1 - erf(x); otherwise rational
 * approximations are computed.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 9.2319   12000       5.1e-16     1.2e-16
 *    IEEE      0,26.6417   30000       5.7e-14     1.5e-14
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition              value returned
 * erfc underflow    x > 9.231948545 (DEC)       0.0
 *
 *
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
*/

/* 
   Wed Jan  3, Josef Leydold:
   made ANSI compliant declaration
*/

#include "cephes_source.h"

/*---------------------------------------------------------------------------*/

static double P[] = {
 2.46196981473530512524E-10,
 5.64189564831068821977E-1,
 7.46321056442269912687E0,
 4.86371970985681366614E1,
 1.96520832956077098242E2,
 5.26445194995477358631E2,
 9.34528527171957607540E2,
 1.02755188689515710272E3,
 5.57535335369399327526E2
};
static double Q[] = {
/* 1.00000000000000000000E0,*/
 1.32281951154744992508E1,
 8.67072140885989742329E1,
 3.54937778887819891062E2,
 9.75708501743205489753E2,
 1.82390916687909736289E3,
 2.24633760818710981792E3,
 1.65666309194161350182E3,
 5.57535340817727675546E2
};
static double R[] = {
 5.64189583547755073984E-1,
 1.27536670759978104416E0,
 5.01905042251180477414E0,
 6.16021097993053585195E0,
 7.40974269950448939160E0,
 2.97886665372100240670E0
};
static double S[] = {
/* 1.00000000000000000000E0,*/
 2.26052863220117276590E0,
 9.39603524938001434673E0,
 1.20489539808096656605E1,
 1.70814450747565897222E1,
 9.60896809063285878198E0,
 3.36907645100081516050E0
};
static double T[] = {
 9.60497373987051638749E0,
 9.00260197203842689217E1,
 2.23200534594684319226E3,
 7.00332514112805075473E3,
 5.55923013010394962768E4
};
static double U[] = {
/* 1.00000000000000000000E0,*/
 3.35617141647503099647E1,
 5.21357949780152679795E2,
 4.59432382970980127987E3,
 2.26290000613890934246E4,
 4.92673942608635921086E4
};

/*---------------------------------------------------------------------------*/

double _unur_cephes_ndtr( double a )
{
double x, y, z;

x = a * SQRTH;
z = fabs(x);

if( z < SQRTH )
	y = 0.5 + 0.5 * _unur_cephes_erf(x);

else
	{
	y = 0.5 * _unur_cephes_erfc(z);

	if( x > 0 )
		y = 1.0 - y;
	}

return(y);
}

/*---------------------------------------------------------------------------*/

double _unur_cephes_erfc( double a )
{
double p,q,x,y,z;


if( a < 0.0 )
	x = -a;
else
	x = a;

if( x < 1.0 )
	return( 1.0 - _unur_cephes_erf(a) );

z = -a * a;

if( z < -MAXLOG )
	{
under:
	if( a < 0 )
		return( 2.0 );
	else
		return( 0.0 );
	}

z = exp(z);

if( x < 8.0 )
	{
	p = _unur_cephes_polevl( x, P, 8 );
	q = _unur_cephes_p1evl( x, Q, 8 );
	}
else
	{
	p = _unur_cephes_polevl( x, R, 5 );
	q = _unur_cephes_p1evl( x, S, 6 );
	}
y = (z * p)/q;

if( a < 0 )
	y = 2.0 - y;

if( _unur_iszero(y) )
	goto under;

return(y);
}

/*---------------------------------------------------------------------------*/

double _unur_cephes_erf( double x)
{
double y, z;

if( fabs(x) > 1.0 )
	return( 1.0 - _unur_cephes_erfc(x) );
z = x * x;
y = x * _unur_cephes_polevl( z, T, 4 ) / _unur_cephes_p1evl( z, U, 5 );
return( y );

}

/*---------------------------------------------------------------------------*/








