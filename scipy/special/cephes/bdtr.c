/*                                                     bdtr.c
 *
 *     Binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, bdtr();
 *
 * y = bdtr( k, n, p );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms 0 through k of the Binomial
 * probability density:
 *
 *   k
 *   --  ( n )   j      n-j
 *   >   (   )  p  (1-p)
 *   --  ( j )
 *  j=0
 *
 * The terms are not summed directly; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p), with p between 0 and 1.
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *  For p between 0.001 and 1:
 *    IEEE     0,100       100000      4.3e-15     2.6e-16
 * See also incbet.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtr domain         k < 0            0.0
 *                     n < k
 *                     x < 0, x > 1
 */
/*							bdtrc()
 *
 *	Complemented binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, bdtrc();
 *
 * y = bdtrc( k, n, p );
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 through n of the Binomial
 * probability density:
 *
 *   n
 *   --  ( n )   j      n-j
 *   >   (   )  p  (1-p)
 *   --  ( j )
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p).
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *  For p between 0.001 and 1:
 *    IEEE     0,100       100000      6.7e-15     8.2e-16
 *  For p between 0 and .001:
 *    IEEE     0,100       100000      1.5e-13     2.7e-15
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtrc domain      x<0, x>1, n<k       0.0
 */
/*							bdtri()
 *
 *	Inverse binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * double p, y, bdtri();
 *
 * p = bdtr( k, n, y );
 *
 * DESCRIPTION:
 *
 * Finds the event probability p such that the sum of the
 * terms 0 through k of the Binomial probability density
 * is equal to the given cumulative probability y.
 *
 * This is accomplished using the inverse beta integral
 * function and the relation
 *
 * 1 - p = incbi( n-k, k+1, y ).
 *
 * ACCURACY:
 *
 * Tested at random points (a,b,p).
 *
 *               a,b                     Relative error:
 * arithmetic  domain     # trials      peak         rms
 *  For p between 0.001 and 1:
 *    IEEE     0,100       100000      2.3e-14     6.4e-16
 *    IEEE     0,10000     100000      6.6e-12     1.2e-13
 *  For p between 10^-6 and 0.001:
 *    IEEE     0,100       100000      2.0e-12     1.3e-14
 *    IEEE     0,10000     100000      1.5e-12     3.2e-14
 * See also incbi.c.
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtri domain     k < 0, n <= k         0.0
 *                  x < 0, x > 1
 */

/*                                                             bdtr() */


/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1987, 1995 by Stephen L. Moshier
 */

#include "mconf.h"

double bdtrc(k, n, p)
int k, n;
double p;
{
    double dk, dn;

    if (npy_isnan(p)) {
	return NPY_NAN;
    }
    if ((p < 0.0) || (p > 1.0)) {
	goto domerr;
    }
    if (k < 0) {
	return 1.0;
    }

    if (n < k) {
      domerr:
	mtherr("bdtrc", DOMAIN);
	return NPY_NAN;
    }

    if (k == n)
	return (0.0);
    dn = n - k;
    if (k == 0) {
	if (p < .01)
	    dk = -expm1(dn * log1p(-p));
	else
	    dk = 1.0 - pow(1.0 - p, dn);
    }
    else {
	dk = k + 1;
	dk = incbet(dk, dn, p);
    }
    return (dk);
}



double bdtr(k, n, p)
int k, n;
double p;
{
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
	goto domerr;
    if ((k < 0) || (n < k)) {
      domerr:
	mtherr("bdtr", DOMAIN);
	return (NPY_NAN);
    }

    if (k == n)
	return (1.0);

    dn = n - k;
    if (k == 0) {
	dk = pow(1.0 - p, dn);
    }
    else {
	dk = k + 1;
	dk = incbet(dn, dk, 1.0 - p);
    }
    return (dk);
}


double bdtri(k, n, y)
int k, n;
double y;
{
    double dk, dn, p;

    if ((y < 0.0) || (y > 1.0))
	goto domerr;
    if ((k < 0) || (n <= k)) {
      domerr:
	mtherr("bdtri", DOMAIN);
	return (NPY_NAN);
    }

    dn = n - k;
    if (k == 0) {
	if (y > 0.8)
	    p = -expm1(log1p(y - 1.0) / dn);
	else
	    p = 1.0 - pow(y, 1.0 / dn);
    }
    else {
	dk = k + 1;
	p = incbet(dn, dk, 0.5);
	if (p > 0.5)
	    p = incbi(dk, dn, 1.0 - y);
	else
	    p = 1.0 - incbi(dn, dk, y);
    }
    return (p);
}
