/*							struve.c
 *
 *      Struve function
 *
 *
 *
 * SYNOPSIS:
 *
 * double v, x, y, struve();
 *
 * y = struve( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the Struve function Hv(x) of order v, argument x.
 * Negative x is rejected unless v is an integer.
 *
 * This module also contains the hypergeometric functions 1F2
 * and 3F0 and a routine for the Bessel function Yv(x) with
 * noninteger v.
 *
 *
 *
 * ACCURACY:
 *
 * Not accurately characterized, but spot checked against tables.
 *
 */


/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#define ANSIPROT
#define DEBUG 0
#ifndef ANSIPROT
double Gamma(), pow(), sqrt(), yn(), yv(), jv(), fabs(), floor();
double sin(), cos();
#else
double onef2( double,double,double,double,double*);
double threef0( double,double,double,double,double*);
extern double fabs(double);
extern double floor(double);
extern double sqrt(double);
extern double sin ( double x );
extern double cos ( double x );
extern double pow ( double x, double y );
extern double Gamma ( double x );
extern double jv ( double n, double x );
extern double yn ( int n, double x );
double struve(double, double);
double yv(double, double);
#endif
static double stop = 1.37e-17;
extern double MACHEP, INFINITY;

double onef2( a, b, c, x, err )
double a, b, c, x;
double *err;
{
double n, a0, sum, t;
double an, bn, cn, max, z;

an = a;
bn = b;
cn = c;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
max = 0.0;

do
	{
	if( an == 0 )
		goto done;
	if( bn == 0 )
		goto error;
	if( cn == 0 )
		goto error;
	if( (a0 > 1.0e34) || (n > 200) )
		goto error;
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = fabs( a0 );
	if( z > max )
		max = z;
	if( sum != 0 )
		t = fabs( a0 / sum );
	else
		t = z;
	}
while( t > stop );

done:

*err = fabs( MACHEP*max /sum );

#if DEBUG
	printf(" onef2 cancellation error %.5E\n", *err );
#endif

goto xit;

error:
#if DEBUG
printf("onef2 does not converge\n");
#endif
*err = 1.0e38;

xit:

#if DEBUG
printf("onef2( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n, sum);
#endif
return(sum);
}




double threef0( a, b, c, x, err )
double a, b, c, x;
double *err;
{
double n, a0, sum, t, conv, conv1;
double an, bn, cn, max, z;

an = a;
bn = b;
cn = c;
a0 = 1.0;
sum = 1.0;
n = 1.0;
t = 1.0;
max = 0.0;
conv = 1.0e38;
conv1 = conv;

do
	{
	if( an == 0.0 )
		goto done;
	if( bn == 0.0 )
		goto done;
	if( cn == 0.0 )
		goto done;
	if( (a0 > 1.0e34) || (n > 200) )
		goto error;
	a0 *= (an * bn * cn * x) / n;
	an += 1.0;
	bn += 1.0;
	cn += 1.0;
	n += 1.0;
	z = fabs( a0 );
	if( z > max )
		max = z;
	if( z >= conv )
		{
		if( (z < max) && (z > conv1) )
			goto done;
		}
	conv1 = conv;
	conv = z;
	sum += a0;
	if( sum != 0 )
		t = fabs( a0 / sum );
	else
		t = z;
	}
while( t > stop );

done:

t = fabs( MACHEP*max/sum );
#if DEBUG
	printf(" threef0 cancellation error %.5E\n", t );
#endif

max = fabs( conv/sum );
if( max > t )
	t = max;
#if DEBUG
	printf(" threef0 convergence %.5E\n", max );
#endif

goto xit;

error:
#if DEBUG
printf("threef0 does not converge\n");
#endif
t = 1.0e38;

xit:

#if DEBUG
printf("threef0( %.2E %.2E %.2E %.5E ) =  %.3E  %.6E\n", a, b, c, x, n, sum);
#endif

*err = t;
return(sum);
}




extern double PI;

double struve( v, x )
double v, x;
{
double y, ya, f, g, h, t;
double onef2err, threef0err;

if (x == 0.0) {
  if ((v>-1) || ((floor(v)-v)==0.5)) return 0.0;
  if (v<-1) {
    if ((int)(floor(0.5-v)-1) % 2) return -INFINITY;
    else return INFINITY;
  }
  return 2.0/PI;
}

f = floor(v);
if( (v < 0) && ( v-f == 0.5 ) )
	{
	y = jv( -v, x );
	f = 1.0 - f;
	g =  2.0 * floor(f/2.0);
	if( g != f )
		y = -y;
	return(y);
	}
t = 0.25*x*x;
f = fabs(x);
g = 1.5 * fabs(v);
if( (f > 30.0) && (f > g) )
	{
	onef2err = 1.0e38;
	y = 0.0;
	}
else
	{
	y = onef2( 1.0, 1.5, 1.5+v, -t, &onef2err );
	}

if( (f < 18.0) || (x < 0.0) )
	{
	threef0err = 1.0e38;
	ya = 0.0;
	}
else
	{
	ya = threef0( 1.0, 0.5, 0.5-v, -1.0/t, &threef0err );
	}

f = sqrt( PI );
h = pow( 0.5*x, v-1.0 );

if( onef2err <= threef0err )
	{
	g = Gamma( v + 1.5 );
	y = y * h * t / ( 0.5 * f * g );
	return(y);
	}
else
	{
	g = Gamma( v + 0.5 );
	ya = ya * h / ( f * g );
	ya = ya + yv( v, x );
	return(ya);
	}
}




/* Bessel function of noninteger order
 */

double yv( v, x )
double v, x;
{
double y, t;
int n;

y = floor( v );
if( y == v )
	{
	n = v;
	y = yn( n, x );
	return( y );
	}
t = PI * v;
y = (cos(t) * jv( v, x ) - jv( -v, x ))/sin(t);
return( y );
}

/* Crossover points between ascending series and asymptotic series
 * for Struve function
 *
 *	 v	 x
 * 
 *	 0	19.2
 *	 1	18.95
 *	 2	19.15
 *	 3	19.3
 *	 5	19.7
 *	10	21.35
 *	20	26.35
 *	30	32.31
 *	40	40.0
 */
