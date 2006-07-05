/*							kn.c
 *
 *	Modified Bessel function, third kind, integer order
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, kn();
 * int n;
 *
 * y = kn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order n of the argument.
 *
 * The range is partitioned into the two intervals [0,9.55] and
 * (9.55, infinity).  An ascending power series is used in the
 * low range, and an asymptotic expansion in the high range.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         3000       1.3e-9      5.8e-11
 *    IEEE      0,30        90000       1.8e-8      3.0e-10
 *
 *  Error is high only near the crossover point x = 9.55
 * between the two expansions used.
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1988, 2000 by Stephen L. Moshier
*/


/*
Algorithm for Kn.
                       n-1 
                   -n   -  (n-k-1)!    2   k
K (x)  =  0.5 (x/2)     >  -------- (-x /4)
 n                      -     k!
                       k=0

                    inf.                                   2   k
       n         n   -                                   (x /4)
 + (-1)  0.5(x/2)    >  {p(k+1) + p(n+k+1) - 2log(x/2)} ---------
                     -                                  k! (n+k)!
                    k=0

where  p(m) is the psi function: p(1) = -EUL and

                      m-1
                       -
      p(m)  =  -EUL +  >  1/k
                       -
                      k=1

For large x,
                                         2        2     2
                                      u-1     (u-1 )(u-3 )
K (z)  =  sqrt(pi/2z) exp(-z) { 1 + ------- + ------------ + ...}
 v                                        1            2
                                    1! (8z)     2! (8z)
asymptotically, where

           2
    u = 4 v .

*/

#include "mconf.h"

#define EUL 5.772156649015328606065e-1
#define MAXFAC 31
#ifdef ANSIPROT
extern double fabs ( double );
extern double exp ( double );
extern double log ( double );
extern double sqrt ( double );
#else
double fabs(), exp(), log(), sqrt();
#endif
extern double MACHEP, MAXNUM, MAXLOG, PI;

double kn( nn, x )
int nn;
double x;
{
double k, kf, nk1f, nkf, zn, t, s, z0, z;
double ans, fn, pn, pk, zmn, tlg, tox;
int i, n;

if( nn < 0 )
	n = -nn;
else
	n = nn;

if( n > MAXFAC )
	{
overf:
	mtherr( "kn", OVERFLOW );
	return( MAXNUM );
	}

if( x <= 0.0 )
	{
	if( x < 0.0 )
		mtherr( "kn", DOMAIN );
	else
		mtherr( "kn", SING );
	return( MAXNUM );
	}


if( x > 9.55 )
	goto asymp;

ans = 0.0;
z0 = 0.25 * x * x;
fn = 1.0;
pn = 0.0;
zmn = 1.0;
tox = 2.0/x;

if( n > 0 )
	{
	/* compute factorial of n and psi(n) */
	pn = -EUL;
	k = 1.0;
	for( i=1; i<n; i++ )
		{
		pn += 1.0/k;
		k += 1.0;
		fn *= k;
		}

	zmn = tox;

	if( n == 1 )
		{
		ans = 1.0/x;
		}
	else
		{
		nk1f = fn/n;
		kf = 1.0;
		s = nk1f;
		z = -z0;
		zn = 1.0;
		for( i=1; i<n; i++ )
			{
			nk1f = nk1f/(n-i);
			kf = kf * i;
			zn *= z;
			t = nk1f * zn / kf;
			s += t;   
			if( (MAXNUM - fabs(t)) < fabs(s) )
				goto overf;
			if( (tox > 1.0) && ((MAXNUM/tox) < zmn) )
				goto overf;
			zmn *= tox;
			}
		s *= 0.5;
		t = fabs(s);
		if( (zmn > 1.0) && ((MAXNUM/zmn) < t) )
			goto overf;
		if( (t > 1.0) && ((MAXNUM/t) < zmn) )
			goto overf;
		ans = s * zmn;
		}
	}


tlg = 2.0 * log( 0.5 * x );
pk = -EUL;
if( n == 0 )
	{
	pn = pk;
	t = 1.0;
	}
else
	{
	pn = pn + 1.0/n;
	t = 1.0/fn;
	}
s = (pk+pn-tlg)*t;
k = 1.0;
do
	{
	t *= z0 / (k * (k+n));
	pk += 1.0/k;
	pn += 1.0/(k+n);
	s += (pk+pn-tlg)*t;
	k += 1.0;
	}
while( fabs(t/s) > MACHEP );

s = 0.5 * s / zmn;
if( n & 1 )
	s = -s;
ans += s;

return(ans);



/* Asymptotic expansion for Kn(x) */
/* Converges to 1.4e-17 for x > 18.4 */

asymp:

if( x > MAXLOG )
	{
	mtherr( "kn", UNDERFLOW );
	return(0.0);
	}
k = n;
pn = 4.0 * k * k;
pk = 1.0;
z0 = 8.0 * x;
fn = 1.0;
t = 1.0;
s = t;
nkf = MAXNUM;
i = 0;
do
	{
	z = pn - pk * pk;
	t = t * z /(fn * z0);
	nk1f = fabs(t);
	if( (i >= n) && (nk1f > nkf) )
		{
		goto adone;
		}
	nkf = nk1f;
	s += t;
	fn += 1.0;
	pk += 2.0;
	i += 1;
	}
while( fabs(t/s) > MACHEP );

adone:
ans = exp(-x) * sqrt( PI/(2.0*x) ) * s;
return(ans);
}
