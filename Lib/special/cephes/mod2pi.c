/* Program to test range reduction of trigonometry functions
 *
 * -- Steve Moshier
 */

#define TPI 6.283185307179586476925

main()
{
char s[40];
double a, n, t, x, y, z;
int lflg;
double floor(), ldexp(), sin();

x = TPI/4.0;
t = 1.0;

loop:

t = 2.0 * t;

/* Stop testing at a point beyond which the integer part of
 * x/2pi cannot be represented exactly by a double precision number.
 * The library trigonometry functions will probably give up long before
 * this point is reached.
 */
if( t > 1.0e16 )
	exit(0);

/* Adjust the following to choose a nontrivial x
 * where test function(x) has a slope of about 1 or more.
 */
x = TPI * t  + 0.5;

z = x;
lflg = 0;

inlup:

/* floor() returns the largest integer less than its argument.
 * If you do not have this, or AINT(), then you may convert x/TPI
 * to a long integer and then back to double; but in that case
 * x will be limited to the largest value that will fit into a
 * long integer.
 */
n = floor( z/TPI );

/* Carefully subtract 2 pi n from x.
 * This is done by subtracting n * 2**k in such a way that there
 * is no arithmetic cancellation error at any step.  The k are the
 * bits in the number 2 pi.
 *
 * If you do not have ldexp(), then you may multiply or
 * divide n by an appropriate power of 2 after each step.
 * For example:
 *  a = z - 4*n;
 *  a -= 2*n;
 *  n /= 4;
 *  a -= n;   n/4
 *  n /= 8;
 *  a -= n;   n/32
 * etc.
 * This will only work if division by a power of 2 is exact.
 */

a = z - ldexp(n, 2);	/* 4n */
a -= ldexp( n, 1);	/* 2n */
a -= ldexp( n, -2 );	/* n/4 */
a -= ldexp( n, -5 );	/* n/32 */
a -= ldexp( n, -9 );	/* n/512 */
a += ldexp( n, -15 );	/* add n/32768 */
a -= ldexp( n, -17 );	/* n/131072 */
a -= ldexp( n, -18 );
a -= ldexp( n, -20 );
a -= ldexp( n, -22 );
a -= ldexp( n, -24 );
a -= ldexp( n, -28 );
a -= ldexp( n, -32 );
a -= ldexp( n, -37 );
a -= ldexp( n, -39 );
a -= ldexp( n, -40 );
a -= ldexp( n, -42 );
a -= ldexp( n, -46 );
a -= ldexp( n, -47 );

/* Subtract what is left of 2 pi n after all the above reductions.
 */
a -= 2.44929359829470635445e-16 * n;

/* If the test is extended too far, it is possible
 * to have chosen the wrong value of n.  The following
 * will fix that, but at some reduction in accuracy.
 */
if( (a > TPI) || (a < -1e-11) )
	{
	z = a;
	lflg += 1;
	printf( "Warning! Reduction failed on first try.\n" );
	goto inlup;
	}
if( a < 0.0 )
	{
	printf( "Warning! Reduced value < 0\n" );
	a += TPI;
	}

/* Compute the test function at x and at a = x mod 2 pi.
 */
y = sin(x);
z = sin(a);
printf( "sin(%.15e) error = %.3e\n", x, y-z );
goto loop;
}

