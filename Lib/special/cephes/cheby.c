/*	cheby.c
 *
 * Program to calculate coefficients of the Chebyshev polynomial
 * expansion of a given input function.  The algorithm computes
 * the discrete Fourier cosine transform of the function evaluated
 * at unevenly spaced points.  Library routine chbevl.c uses the
 * coefficients to calculate an approximate value of the original
 * function.
 *    -- S. L. Moshier
 */

extern double PI;		/* 3.14159...	*/
extern double PIO2;
double cosi[33] = {0.0,};	/* cosine array for Fourier transform	*/
double func[65] = {0.0,};	/* values of the function		*/
double cos(), log(), exp(), sqrt();

main()
{
double c, r, s, t, x, y, z, temp;
double low, high, dtemp;
long n;
int i, ii, j, n2, k, rr, invflg;
short *p;
char st[40];

low = 0.0;		/* low end of approximation interval		*/
high = 1.0;		/* high end					*/
invflg = 0;		/* set to 1 if inverted interval, else zero	*/
/* Note: inverted interval goes from 1/high to 1/low	*/
z = 0.0;
n = 64;			/* will find 64 coefficients			*/
			/* but use only those greater than roundoff error */
n2 = n/2;
t = n;
t = PI/t;

/* calculate array of cosines */
puts("calculating cosines");
s = 1.0;
cosi[0] = 1.0;
i = 1;
while( i < 32 )
	{
	y = cos( s * t );
	cosi[i] = y;
	s += 1.0;
	++i;
	}
cosi[32] = 0.0;

/*							cheby.c 2 */

/* calculate function at special values of the argument */
puts("calculating function values");
x = low;
y = high;
if( invflg && (low != 0.0) )
	{	/* inverted interval */
	temp = 1.0/x;
	x = 1.0/y;
	y = temp;
	}
r = (x + y)/2.0;
printf( "center %.15E  ", r);
s = (y - x)/2.0;
printf( "width %.15E\n", s);
i = 0;
while( i < 65 )
	{
	if( i < n2 )
		c = cosi[i];
	else
		c = -cosi[64-i];
	temp = r + s * c;
/* if inverted interval, compute function(1/x)	*/
	if( invflg && (temp != 0.0) )
		temp = 1.0/temp;

	printf( "%.15E  ", temp );

/* insert call to function routine here:	*/
/**********************************/

	if( temp == 0.0 )
		y = 1.0;
	else
		y = exp( temp * log(2.0) );

/**********************************/
	func[i] = y;
	printf( "%.15E\n", y );
	++i;
	}

/*							cheby.c 3 */

puts( "calculating Chebyshev coefficients");
rr = 0;
while( rr < 65 )
	{
	z = func[0]/2.0;
	j = 1;
	while( j < 65 )
		{
		k = (rr * j)/n2;
		i = rr * j - n2 * k;
		k &= 3;
		if( k == 0 )
			c = cosi[i];
		if( k == 1 )
			{
			i = 32-i;
			c = -cosi[i];
			if( i == 32 )
				c = -c;
			}
		if( k == 2 )
			{
			c = -cosi[i];
			}
		if( k == 3 )
			{
			i = 32-i;
			c = cosi[i];
			}
		if( i != 32)
			{
			temp = func[j];
			temp = c * temp;
			z += temp;
			}
		++j;
		}

	if( i != 32 )
		{
		temp /= 2.0;
		z = z - temp;
		}
	z *= 2.0;
	temp = n;
	z /= temp;
	dtemp = z;
	++rr;
	sprintf( st, "/* %.16E */", dtemp );
	puts( st );
	}
}
