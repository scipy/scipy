
/* Arithmetic operations on polynomials with rational coefficients
 *
 * In the following descriptions a, b, c are polynomials of degree
 * na, nb, nc respectively.  The degree of a polynomial cannot
 * exceed a run-time value MAXPOL.  An operation that attempts
 * to use or generate a polynomial of higher degree may produce a
 * result that suffers truncation at degree MAXPOL.  The value of
 * MAXPOL is set by calling the function
 *
 *     polini( maxpol );
 *
 * where maxpol is the desired maximum degree.  This must be
 * done prior to calling any of the other functions in this module.
 * Memory for internal temporary polynomial storage is allocated
 * by polini().
 *
 * Each polynomial is represented by an array containing its
 * coefficients, together with a separately declared integer equal
 * to the degree of the polynomial.  The coefficients appear in
 * ascending order; that is,
 *
 *                                        2                      na
 * a(x)  =  a[0]  +  a[1] * x  +  a[2] * x   +  ...  +  a[na] * x  .
 *
 *
 *
 * `a', `b', `c' are arrays of fracts.
 * poleva( a, na, &x, &sum );	Evaluate polynomial a(t) at t = x.
 * polprt( a, na, D );		Print the coefficients of a to D digits.
 * polclr( a, na );		Set a identically equal to zero, up to a[na].
 * polmov( a, na, b );		Set b = a.
 * poladd( a, na, b, nb, c );	c = b + a, nc = max(na,nb)
 * polsub( a, na, b, nb, c );	c = b - a, nc = max(na,nb)
 * polmul( a, na, b, nb, c );	c = b * a, nc = na+nb
 *
 *
 * Division:
 *
 * i = poldiv( a, na, b, nb, c );	c = b / a, nc = MAXPOL
 *
 * returns i = the degree of the first nonzero coefficient of a.
 * The computed quotient c must be divided by x^i.  An error message
 * is printed if a is identically zero.
 *
 *
 * Change of variables:
 * If a and b are polynomials, and t = a(x), then
 *     c(t) = b(a(x))
 * is a polynomial found by substituting a(x) for t.  The
 * subroutine call for this is
 *
 * polsbt( a, na, b, nb, c );
 *
 *
 * Notes:
 * poldiv() is an integer routine; poleva() is double.
 * Any of the arguments a, b, c may refer to the same array.
 *
 */

#include "mconf.h"
#include <malloc.h>
#ifndef NULL
#define NULL 0
#endif
typedef struct{
	double n;
	double d;
	}fract;

void radd(), rsub(), rmul(), rdiv();

/* near pointer version of malloc() */
/*
#define malloc _nmalloc
#define free _nfree
*/
/* Pointers to internal arrays.  Note poldiv() allocates
 * and deallocates some temporary arrays every time it is called.
 */
static fract *pt1 = 0;
static fract *pt2 = 0;
static fract *pt3 = 0;

/* Maximum degree of polynomial. */
int MAXPOL = 0;
extern int MAXPOL;

/* Number of bytes (chars) in maximum size polynomial. */
static int psize = 0;


/* Initialize max degree of polynomials
 * and allocate temporary storage.
 */
void polini( maxdeg )
int maxdeg;
{

MAXPOL = maxdeg;
psize = (maxdeg + 1) * sizeof(fract);

/* Release previously allocated memory, if any. */
if( pt3 )
	free(pt3);
if( pt2 )
	free(pt2);
if( pt1 )
	free(pt1);

/* Allocate new arrays */
pt1 = (fract * )malloc(psize); /* used by polsbt */
pt2 = (fract * )malloc(psize); /* used by polsbt */
pt3 = (fract * )malloc(psize); /* used by polmul */

/* Report if failure */
if( (pt1 == NULL) || (pt2 == NULL) || (pt3 == NULL) )
	{
	mtherr( "polini", ERANGE );
	exit(1);
	}
}



/* Print the coefficients of a, with d decimal precision.
 */
static char *form = "abcdefghijk";

void polprt( a, na, d )
fract a[];
int na, d;
{
int i, j, d1;
char *p;

/* Create format descriptor string for the printout.
 * Do this partly by hand, since sprintf() may be too
 * bug-ridden to accomplish this feat by itself.
 */
p = form;
*p++ = '%';
d1 = d + 8;
sprintf( p, "%d ", d1 );
p += 1;
if( d1 >= 10 )
	p += 1;
*p++ = '.';
sprintf( p, "%d ", d );
p += 1;
if( d >= 10 )
	p += 1;
*p++ = 'e';
*p++ = ' ';
*p++ = '\0';


/* Now do the printing.
 */
d1 += 1;
j = 0;
for( i=0; i<=na; i++ )
	{
/* Detect end of available line */
	j += d1;
	if( j >= 78 )
		{
		printf( "\n" );
		j = d1;
		}
	printf( form, a[i].n );
	j += d1;
	if( j >= 78 )
		{
		printf( "\n" );
		j = d1;
		}
	printf( form, a[i].d );
	}
printf( "\n" );
}



/* Set a = 0.
 */
void polclr( a, n )
fract a[];
int n;
{
int i;

if( n > MAXPOL )
	n = MAXPOL;
for( i=0; i<=n; i++ )
	{
	a[i].n = 0.0;
	a[i].d = 1.0;
	}
}



/* Set b = a.
 */
void polmov( a, na, b )
fract a[], b[];
int na;
{
int i;

if( na > MAXPOL )
	na = MAXPOL;

for( i=0; i<= na; i++ )
	{
	b[i].n = a[i].n;
	b[i].d = a[i].d;
	}
}


/* c = b * a.
 */
void polmul( a, na, b, nb, c )
fract a[], b[], c[];
int na, nb;
{
int i, j, k, nc;
fract temp;
fract *p;

nc = na + nb;
polclr( pt3, MAXPOL );

p = &a[0];
for( i=0; i<=na; i++ )
	{
	for( j=0; j<=nb; j++ )
		{
		k = i + j;
		if( k > MAXPOL )
			break;
		rmul( p, &b[j], &temp ); /*pt3[k] += a[i] * b[j];*/
		radd( &temp, &pt3[k], &pt3[k] );
		}
	++p;
	}

if( nc > MAXPOL )
	nc = MAXPOL;
for( i=0; i<=nc; i++ )
	{
	c[i].n = pt3[i].n;
	c[i].d = pt3[i].d;
	}
}



 
/* c = b + a.
 */
void poladd( a, na, b, nb, c )
fract a[], b[], c[];
int na, nb;
{
int i, n;


if( na > nb )
	n = na;
else
	n = nb;

if( n > MAXPOL )
	n = MAXPOL;

for( i=0; i<=n; i++ )
	{
	if( i > na )
		{
		c[i].n = b[i].n;
		c[i].d = b[i].d;
		}
	else if( i > nb )
		{
		c[i].n = a[i].n;
		c[i].d = a[i].d;
		}
	else
		{
		radd( &a[i], &b[i], &c[i] ); /*c[i] = b[i] + a[i];*/
		}
	}
}

/* c = b - a.
 */
void polsub( a, na, b, nb, c )
fract a[], b[], c[];
int na, nb;
{
int i, n;


if( na > nb )
	n = na;
else
	n = nb;

if( n > MAXPOL )
	n = MAXPOL;

for( i=0; i<=n; i++ )
	{
	if( i > na )
		{
		c[i].n = b[i].n;
		c[i].d = b[i].d;
		}
	else if( i > nb )
		{
		c[i].n = -a[i].n;
		c[i].d = a[i].d;
		}
	else
		{
		rsub( &a[i], &b[i], &c[i] ); /*c[i] = b[i] - a[i];*/
		}
	}
}



/* c = b/a
 */
int poldiv( a, na, b, nb, c )
fract a[], b[], c[];
int na, nb;
{
fract *ta, *tb, *tq;
fract quot;
fract temp;
int i, j, k, sing;

sing = 0;

/* Allocate temporary arrays.  This would be quicker
 * if done automatically on the stack, but stack space
 * may be hard to obtain on a small computer.
 */
ta = (fract * )malloc( psize );
polclr( ta, MAXPOL );
polmov( a, na, ta );

tb = (fract * )malloc( psize );
polclr( tb, MAXPOL );
polmov( b, nb, tb );

tq = (fract * )malloc( psize );
polclr( tq, MAXPOL );

/* What to do if leading (constant) coefficient
 * of denominator is zero.
 */
if( a[0].n == 0.0 )
	{
	for( i=0; i<=na; i++ )
		{
		if( ta[i].n != 0.0 )
			goto nzero;
		}
	mtherr( "poldiv", SING );
	goto done;

nzero:
/* Reduce the degree of the denominator. */
	for( i=0; i<na; i++ )
		{
		ta[i].n = ta[i+1].n;
		ta[i].d = ta[i+1].d;
		}
	ta[na].n = 0.0;
	ta[na].d = 1.0;

	if( b[0].n != 0.0 )
		{
/* Optional message:
		printf( "poldiv singularity, divide quotient by x\n" );
*/
		sing += 1;
		}
	else
		{
/* Reduce degree of numerator. */
		for( i=0; i<nb; i++ )
			{
			tb[i].n = tb[i+1].n;
			tb[i].d = tb[i+1].d;
			}
		tb[nb].n = 0.0;
		tb[nb].d = 1.0;
		}
/* Call self, using reduced polynomials. */
	sing += poldiv( ta, na, tb, nb, c );
	goto done;
	}

/* Long division algorithm.  ta[0] is nonzero.
 */
for( i=0; i<=MAXPOL; i++ )
	{
	rdiv( &ta[0], &tb[i], &quot ); /*quot = tb[i]/ta[0];*/
	for( j=0; j<=MAXPOL; j++ )
		{
		k = j + i;
		if( k > MAXPOL )
			break;

		rmul( &ta[j], &quot, &temp ); /*tb[k] -= quot * ta[j];*/
		rsub( &temp, &tb[k], &tb[k] );
		}
	tq[i].n = quot.n;
	tq[i].d = quot.d;
	}
/* Send quotient to output array. */
polmov( tq, MAXPOL, c );

done:

/* Restore allocated memory. */
free(tq);
free(tb);
free(ta);
return( sing );
}




/* Change of variables
 * Substitute a(y) for the variable x in b(x).
 * x = a(y)
 * c(x) = b(x) = b(a(y)).
 */

void polsbt( a, na, b, nb, c )
fract a[], b[], c[];
int na, nb;
{
int i, j, k, n2;
fract temp;
fract *p;

/* 0th degree term:
 */
polclr( pt1, MAXPOL );
pt1[0].n = b[0].n;
pt1[0].d = b[0].d;

polclr( pt2, MAXPOL );
pt2[0].n = 1.0;
pt2[0].d = 1.0;
n2 = 0;
p = &b[1];

for( i=1; i<=nb; i++ )
	{
/* Form ith power of a. */
	polmul( a, na, pt2, n2, pt2 );
	n2 += na;
/* Add the ith coefficient of b times the ith power of a. */
	for( j=0; j<=n2; j++ )
		{
		if( j > MAXPOL )
			break;
		rmul( &pt2[j], p, &temp ); /*pt1[j] += b[i] * pt2[j];*/
		radd( &temp, &pt1[j], &pt1[j] );
		}
	++p;
	}

k = n2 + nb;
if( k > MAXPOL )
	k = MAXPOL;
for( i=0; i<=k; i++ )
	{
	c[i].n = pt1[i].n;
	c[i].d = pt1[i].d;
	}
}




/* Evaluate polynomial a(t) at t = x.
 */
void poleva( a, na, x, s )
fract a[];
int na;
fract *x;
fract *s;
{
int i;
fract temp;

s->n = a[na].n;
s->d = a[na].d;
for( i=na-1; i>=0; i-- )
	{
	rmul( s, x, &temp ); /*s = s * x + a[i];*/
	radd( &a[i], &temp, s );
	}
}

