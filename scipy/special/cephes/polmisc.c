
/* Square root, sine, cosine, and arctangent of polynomial.
 * See polyn.c for data structures and discussion.
 */

#include <stdio.h>
#include <stdlib.h>
#include "mconf.h"

/* Highest degree of polynomial to be handled
   by the polyn.c subroutine package.  */
#define N 16
/* Highest degree actually initialized at runtime.  */
extern int MAXPOL;

/* Taylor series coefficients for various functions
 */
double patan[N+1] = {
  0.0,     1.0,      0.0, -1.0/3.0,     0.0,
  1.0/5.0, 0.0, -1.0/7.0,      0.0, 1.0/9.0, 0.0, -1.0/11.0,
  0.0, 1.0/13.0, 0.0, -1.0/15.0, 0.0 };

double psin[N+1] = {
  0.0, 1.0, 0.0,   -1.0/6.0,  0.0, 1.0/120.0,  0.0,
  -1.0/5040.0, 0.0, 1.0/362880.0, 0.0, -1.0/39916800.0,
  0.0, 1.0/6227020800.0, 0.0, -1.0/1.307674368e12, 0.0};

double pcos[N+1] = {
  1.0, 0.0,   -1.0/2.0,  0.0, 1.0/24.0,  0.0,
  -1.0/720.0, 0.0, 1.0/40320.0, 0.0, -1.0/3628800.0, 0.0,
  1.0/479001600.0, 0.0, -1.0/8.7179291e10, 0.0, 1.0/2.0922789888e13};

double pasin[N+1] = {
  0.0,     1.0,  0.0, 1.0/6.0,  0.0,
  3.0/40.0, 0.0, 15.0/336.0, 0.0, 105.0/3456.0, 0.0, 945.0/42240.0,
  0.0, 10395.0/599040.0 , 0.0, 135135.0/9676800.0 , 0.0
};

/* Square root of 1 + x.  */
double psqrt[N+1] = {
  1.0, 1./2., -1./8., 1./16., -5./128., 7./256., -21./1024., 33./2048.,
  -429./32768., 715./65536., -2431./262144., 4199./524288., -29393./4194304.,
  52003./8388608., -185725./33554432., 334305./67108864.,
  -9694845./2147483648.};

/* Arctangent of the ratio num/den of two polynomials.
 */
void
polatn( num, den, ans, nn )
     double num[], den[], ans[];
     int nn;
{
  double a, t;
  double *polq, *polu, *polt;
  int i;

  if (nn > N)
    {
      mtherr ("polatn", OVERFLOW);
      return;
    }
  /* arctan( a + b ) = arctan(a) + arctan( b/(1 + ab + a**2) ) */
  t = num[0];
  a = den[0];
  if( (t == 0.0) && (a == 0.0 ) )
    {
      t = num[1];
      a = den[1];
    }
  t = atan2( t, a );  /* arctan(num/den), the ANSI argument order */
  polq = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polu = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polt = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polclr( polq, MAXPOL );
  i = poldiv( den, nn, num, nn, polq );
  a = polq[0]; /* a */
  polq[0] = 0.0; /* b */
  polmov( polq, nn, polu ); /* b */
  /* Form the polynomial
     1 + ab + a**2
     where a is a scalar.  */
  for( i=0; i<=nn; i++ )
    polu[i] *= a;
  polu[0] += 1.0 + a * a;
  poldiv( polu, nn, polq, nn, polt ); /* divide into b */
  polsbt( polt, nn, patan, nn, polu ); /* arctan(b)  */
  polu[0] += t; /* plus arctan(a) */
  polmov( polu, nn, ans );
  free( polt );
  free( polu );
  free( polq );
}



/* Square root of a polynomial.
 * Assumes the lowest degree nonzero term is dominant
 * and of even degree.  An error message is given
 * if the Newton iteration does not converge.
 */
void
polsqt( pol, ans, nn )
     double pol[], ans[];
     int nn;
{
  double t;
  double *x, *y;
  int i, n;
#if 0
  double z[N+1];
  double u;
#endif

  if (nn > N)
    {
      mtherr ("polatn", OVERFLOW);
      return;
    }
  x = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  y = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polmov( pol, nn, x );
  polclr( y, MAXPOL );

  /* Find lowest degree nonzero term.  */
  t = 0.0;
  for( n=0; n<nn; n++ )
    {
      if( x[n] != 0.0 )
	goto nzero;
    }
  polmov( y, nn, ans );
  return;

nzero:

  if( n > 0 )
    {
      if (n & 1)
        {
	  printf("error, sqrt of odd polynomial\n");
	  return;
	}
      /* Divide by x^n.  */
      y[n] = x[n];
      poldiv (y, nn, pol, N, x);
    }

  t = x[0];
  for( i=1; i<=nn; i++ )
    x[i] /= t;
  x[0] = 0.0;
  /* series development sqrt(1+x) = 1  +  x / 2  -  x**2 / 8  +  x**3 / 16
     hopes that first (constant) term is greater than what follows   */
  polsbt( x, nn, psqrt, nn, y);
  t = sqrt( t );
  for( i=0; i<=nn; i++ )
    y[i] *= t;

  /* If first nonzero coefficient was at degree n > 0, multiply by
     x^(n/2).  */
  if (n > 0)
    {
      polclr (x, MAXPOL);
      x[n/2] = 1.0;
      polmul (x, nn, y, nn, y);
    }
#if 0
/* Newton iterations */
for( n=0; n<10; n++ )
	{
	poldiv( y, nn, pol, nn, z );
	poladd( y, nn, z, nn, y );
	for( i=0; i<=nn; i++ )
		y[i] *= 0.5;
	for( i=0; i<=nn; i++ )
		{
		u = fabs( y[i] - z[i] );
		if( u > 1.0e-15 )
			goto more;
		}
	goto done;
more:	;
	}
printf( "square root did not converge\n" );
done:
#endif /* 0 */

polmov( y, nn, ans );
free( y );
free( x );
}



/* Sine of a polynomial.
 * The computation uses
 *     sin(a+b) = sin(a) cos(b) + cos(a) sin(b)
 * where a is the constant term of the polynomial and
 * b is the sum of the rest of the terms.
 * Since sin(b) and cos(b) are computed by series expansions,
 * the value of b should be small.
 */
void
polsin( x, y, nn )
     double x[], y[];
     int nn;
{
  double a, sc;
  double *w, *c;
  int i;

  if (nn > N)
    {
      mtherr ("polatn", OVERFLOW);
      return;
    }
  w = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  c = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polmov( x, nn, w );
  polclr( c, MAXPOL );
  polclr( y, nn );
  /* a, in the description, is x[0].  b is the polynomial x - x[0].  */
  a = w[0];
  /* c = cos (b) */
  w[0] = 0.0;
  polsbt( w, nn, pcos, nn, c );
  sc = sin(a);
  /* sin(a) cos (b) */
  for( i=0; i<=nn; i++ )
    c[i] *= sc;
  /* y = sin (b)  */
  polsbt( w, nn, psin, nn, y );
  sc = cos(a);
  /* cos(a) sin(b) */
  for( i=0; i<=nn; i++ )
    y[i] *= sc;
  poladd( c, nn, y, nn, y );
  free( c );
  free( w );
}


/* Cosine of a polynomial.
 * The computation uses
 *     cos(a+b) = cos(a) cos(b) - sin(a) sin(b)
 * where a is the constant term of the polynomial and
 * b is the sum of the rest of the terms.
 * Since sin(b) and cos(b) are computed by series expansions,
 * the value of b should be small.
 */
void
polcos( x, y, nn )
     double x[], y[];
     int nn;
{
  double a, sc;
  double *w, *c;
  int i;

  if (nn > N)
    {
      mtherr ("polatn", OVERFLOW);
      return;
    }
  w = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  c = (double * )malloc( (MAXPOL+1) * sizeof (double) );
  polmov( x, nn, w );
  polclr( c, MAXPOL );
  polclr( y, nn );
  a = w[0];
  w[0] = 0.0;
  /* c = cos(b)  */
  polsbt( w, nn, pcos, nn, c );
  sc = cos(a);
  /* cos(a) cos(b)  */
  for( i=0; i<=nn; i++ )
    c[i] *= sc;
  /* y = sin(b) */
  polsbt( w, nn, psin, nn, y );
  sc = sin(a);
  /* sin(a) sin(b) */
  for( i=0; i<=nn; i++ )
    y[i] *= sc;
  polsub( y, nn, c, nn, y );
  free( c );
  free( w );
}
