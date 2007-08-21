/*							cpmul.c
 *
 *	Multiply two polynomials with complex coefficients
 *
 *
 *
 * SYNOPSIS:
 *
 * typedef struct
 *		{
 *		double r;
 *		double i;
 *		}cmplx;
 *
 * cmplx a[], b[], c[];
 * int da, db, dc;
 *
 * cpmul( a, da, b, db, c, &dc );
 *
 *
 *
 * DESCRIPTION:
 *
 * The two argument polynomials are multiplied together, and
 * their product is placed in c.
 *
 * Each polynomial is represented by its coefficients stored
 * as an array of complex number structures (see the typedef).
 * The degree of a is da, which must be passed to the routine
 * as an argument; similarly the degree db of b is an argument.
 * Array a has da + 1 elements and array b has db + 1 elements.
 * Array c must have storage allocated for at least da + db + 1
 * elements.  The value da + db is returned in dc; this is
 * the degree of the product polynomial.
 *
 * Polynomial coefficients are stored in ascending order; i.e.,
 * a(x) = a[0]*x**0 + a[1]*x**1 + ... + a[da]*x**da.
 *
 *
 * If desired, c may be the same as either a or b, in which
 * case the input argument array is replaced by the product
 * array (but only up to terms of degree da + db).
 *
 */

/*							cpmul	*/

typedef struct
	{
	double r;
	double i;
	}cmplx;

void cpmul( cmplx*, int, cmplx*, int, cmplx*, int* );

void
cpmul( a, da, b, db, c, dc )
cmplx *a, *b, *c;
int da, db;
int *dc;
{
int i, j, k;
cmplx y;
register cmplx *pa, *pb, *pc;

if( da > db )	/* Know which polynomial has higher degree */
	{
	i = da;	/* Swapping is OK because args are on the stack */
	da = db;
	db = i;
	pa = a;
	a = b;
	b = pa;
	}
	
k = da + db;
*dc = k;		/* Output the degree of the product */
pc = &c[db+1];
for( i=db+1; i<=k; i++ )	/* Clear high order terms of output */
	{
	pc->r = 0;
	pc->i = 0;
	pc++;
	}
/* To permit replacement of input, work backward from highest degree */
pb = &b[db];
for( j=0; j<=db; j++ )
	{
	pa = &a[da];
	pc = &c[k-j];
	for( i=0; i<da; i++ )
		{
		y.r = pa->r * pb->r  -  pa->i * pb->i;	/* cmpx multiply */
		y.i = pa->r * pb->i  +  pa->i * pb->r;
		pc->r += y.r;	/* accumulate partial product */
		pc->i += y.i;
		pa--;
		pc--;
		}
	y.r = pa->r * pb->r  -  pa->i * pb->i;	/* replace last term,	*/
	y.i = pa->r * pb->i  +  pa->i * pb->r;	/* ...do not accumulate	*/
	pc->r = y.r;
	pc->i = y.i;
	pb--;
	}
}
