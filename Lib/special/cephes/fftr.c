/*							fftr.c
 *
 *	FFT of Real Valued Sequence
 *
 *
 *
 * SYNOPSIS:
 *
 * double x[], sine[];
 * int m;
 *
 * fftr( x, m, sine );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the (complex valued) discrete Fourier transform of
 * the real valued sequence x[].  The input sequence x[] contains
 * n = 2**m samples.  The program fills array sine[k] with
 * n/4 + 1 values of sin( 2 PI k / n ).
 *
 * Data format for complex valued output is real part followed
 * by imaginary part.  The output is developed in the input
 * array x[].
 *
 * The algorithm takes advantage of the fact that the FFT of an
 * n point real sequence can be obtained from an n/2 point
 * complex FFT.
 *
 * A radix 2 FFT algorithm is used.
 *
 * Execution time on an LSI-11/23 with floating point chip
 * is 1.0 sec for n = 256.
 *
 *
 *
 * REFERENCE:
 *
 * E. Oran Brigham, The Fast Fourier Transform;
 * Prentice-Hall, Inc., 1974
 *
 */
/*			fftr.c			*/


static short n0 = 0;
static short n4 = 0;
static short msav = 0;

extern double PI;

#ifndef ANSIPROT
double sin();
static int bitrv();
#else
static int bitrv(int, int);
#endif

fftr( x, m0, sine )
double x[];
int m0;
double sine[];
{
int th, nd, pth, nj, dth, m;
int n, n1, n2, j, k, l, r;
double xr, xi, tr, ti, co, si;
double a, b, c, d, bc, cs, bs, cc;
double *p, *q;

/* Array x assumed filled with real-valued data */
/* m0 = log2(n0)			*/
/* n0 is the number of real data samples */

if( m0 != msav )
	{
	msav = m0;

	/* Find n0 = 2**m0	*/
	n0 = 1;
	for( j=0; j<m0; j++ )
		n0 <<= 1;

	n4 = n0 >> 2;

	/* Calculate array of sines */
	xr = 2.0 * PI / n0;
	for( j=0; j<=n4; j++ )
		sine[j] = sin( j * xr );
	}

n = n0 >> 1;	/* doing half length transform */
m = m0 - 1;


/*							fftr.c	*/

/*  Complex Fourier Transform of n Complex Data Points */

/*	First, bit reverse the input data	*/

for( k=0; k<n; k++ )
	{
	j = bitrv( k, m );
	if( j > k )
		{ /* executed approx. n/2 times */
		p = &x[2*k];
		tr = *p++;
		ti = *p;
		q = &x[2*j+1];
		*p = *q;
		*(--p) = *(--q);
		*q++ = tr;
		*q = ti;
		}
	}
	
/*							fftr.c	*/
/*			Radix 2 Complex FFT			*/
n2 = n/2;
nj = 1;
pth = 1;
dth = 0;
th = 0;

for( l=0; l<m; l++ )
	{	/* executed log2(n) times, total */
	j = 0;
	do
		{	/* executed n-1 times, total */
		r = th << 1;
		si = sine[r];
		co = sine[ n4 - r ];
		if( j >= pth )
			{
			th -= dth;
			co = -co;
			}
		else
			th += dth;

		nd = j;

		do
			{ /* executed n/2 log2(n) times, total */
			r = (nd << 1) + (nj << 1);
			p = &x[ r ];
			xr = *p++;
			xi = *p;
			tr = xr * co + xi * si;
			ti = xi * co - xr * si;
			r = nd << 1;
			q = &x[ r ];
			xr = *q++;
			xi = *q;
			*p = xi - ti;
			*(--p) = xr - tr;
			*q = xi + ti;
			*(--q) = xr + tr;
			nd += nj << 1;
			}
		while( nd < n );
		}
	while( ++j < nj );

	n2 >>= 1;
	dth = n2;
	pth = nj;
	nj <<= 1;
	}

/*							fftr.c	*/

/*	Special trick algorithm			*/
/*	converts to spectrum of real series	*/

/* Highest frequency term; add space to input array if wanted */
/*
x[2*n] = x[0] - x[1];
x[2*n+1] = 0.0;
*/

/* Zero frequency term */
x[0] = x[0] + x[1];
x[1] = 0.0;
n2 = n/2;

for( j=1; j<=n2; j++ )
	{	/* executed n/2 times */
	si = sine[j];
	co = sine[ n4 - j ];
	p = &x[ 2*j ];
	xr = *p++;
	xi = *p;
	q = &x[ 2*(n-j) ];
	tr = *q++;
	ti = *q;
	a = xr + tr;
	b = xi + ti;
	c = xr - tr;
	d = xi - ti;
	bc = b * co;
	cs = c * si;
	bs = b * si;
	cc = c * co;
	*p = ( d - bs - cc )/2.0;
	*(--p) = ( a + bc - cs )/2.0;
	*q = -( d + bs + cc )/2.0;
	*(--q) = ( a - bc + cs )/2.0;
	}

return(0);
}

/*							fftr.c	*/

/*	Bit reverser	*/

int bitrv( j, m )
int j, m;
{
register int j1, ans;
short k;

ans = 0;
j1 = j;

for( k=0; k<m; k++ )
	{
	ans = (ans << 1) + (j1 & 1);
	j1 >>= 1;
	}

return( ans );
}
