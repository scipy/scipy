/* fltest.c
 * Test program for floor(), frexp(), ldexp()
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/



#include "mconf.h"
extern double MACHEP;
#define UTH -1023

main()
{
double x, y, y0, z, f, x00, y00;
int i, j, k, e, e0;
int errfr, errld, errfl, underexp, err, errth, e00;
double frexp(), ldexp(), floor();


/*
if( 1 )
	goto flrtst;
*/

printf( "Testing frexp() and ldexp().\n" );
errfr = 0;
errld = 0;
underexp = 0;
f = 1.0;
x00 = 2.0;
y00 = 0.5;
e00 = 2;

for( j=0; j<20; j++ )
{
if( j == 10 )
	{
	f = 1.0;
	x00 = 2.0;
	e00 = 1;
/* Find 2**(2**10) / 2 */
#ifdef DEC
	for( i=0; i<5; i++ )
#else
	for( i=0; i<9; i++ )
#endif
		{
		x00 *= x00;
		e00 += e00;
		}
	y00 = x00/2.0;
	x00 = x00 * y00;
	e00 += e00;
	y00 = 0.5;
	}
x = x00 * f;
y0 = y00 * f;
e0 = e00;
for( i=0; i<2200; i++ )
	{
	x /= 2.0;
	e0 -= 1;
	if( x == 0.0 )
		{
		if( f == 1.0 )
			underexp = e0;
		y0 = 0.0;
		e0 = 0;
		}
	y = frexp( x, &e );
	if( (e0 < -1023) && (e != e0) )
		{
		if( e == (e0 - 1) )
			{
			e += 1;
			y /= 2.0;
			}
		if( e == (e0 + 1) )
			{
			e -= 1;
			y *= 2.0;
			}
		}
	err = y - y0;
	if( y0 != 0.0 )
		err /= y0;
	if( err < 0.0 )
		err = -err;
	if( e0 > -1023 )
		errth = 0.0;
	else
		{/* Denormal numbers may have rounding errors */
		if( e0 == -1023 )
			{
			errth = 2.0 * MACHEP;
			}
		else
			{
			errth *= 2.0;
			}
		}

	if( (x != 0.0) && ((err > errth) || (e != e0)) )
		{
		printf( "Test %d: ", j+1 );
		printf( " frexp( %.15e) =?= %.15e * 2**%d;", x, y, e );
		printf( " should be %.15e * 2**%d\n", y0, e0 );
		errfr += 1;
		}
	y = ldexp( x, 1-e0 );
	err = y - 1.0;
	if( err < 0.0 )
		err = -err;
	if( (err > errth) && ((x == 0.0) && (y != 0.0)) )
		{
		printf( "Test %d: ", j+1 );
		printf( "ldexp( %.15e, %d ) =?= %.15e;", x, 1-e0, y );
		if( x != 0.0 )
			printf( " should be %.15e\n", f );
		else
			printf( " should be %.15e\n", 0.0 );
		errld += 1;
		}
	if( x == 0.0 )
		{
		break;
		}
	}
f = f * 1.08005973889;
}


x = 2.22507385850720138309e-308;
for (i = 0; i < 52; i++)
  {
    y = ldexp (x, -i);
    z = ldexp (y, i);
    if (x != z)
      {
	printf ("x %.16e, i %d, y %.16e, z %.16e\n", x, i, y, z);
	errld += 1;
      }
  }


if( (errld == 0) && (errfr == 0) )
	{
	printf( "No errors found.\n" );
	}

flrtst:

printf( "Testing floor().\n" );
errfl = 0;

f = 1.0/MACHEP;
x00 = 1.0;
for( j=0; j<57; j++ )
{
x = x00 - 1.0;
for( i=0; i<128; i++ )
	{
	y = floor(x);
	if( y != x )
		{
		flierr( x, y, j );
		errfl += 1;
		}
/* Warning! the if() statement is compiler dependent,
 * since x-0.49 may be held in extra precision accumulator
 * so would never compare equal to x!  The subroutine call
 * y = floor() forces z to be stored as a double and reloaded
 * for the if() statement.
 */
	z = x - 0.49;
	y = floor(z);
	if( z == x )
		break;
	if( y != (x - 1.0) )
		{
		flierr( z, y, j );
		errfl += 1;
		}

	z = x + 0.49;
	y = floor(z);
	if( z != x )
		{
		if( y != x )
			{
			flierr( z, y, j );
			errfl += 1;
			}
		}
	x = -x;
	y = floor(x);
	if( z != x )
		{
		if( y != x )
			{
			flierr( x, y, j );
			errfl += 1;
			}
		}
	z = x + 0.49;
	y = floor(z);
	if( z != x )
		{
		if( y != x )
			{
			flierr( z, y, j );
			errfl += 1;
			}
		}
	z = x - 0.49;
	y = floor(z);
	if( z != x )
		{
		if( y != (x - 1.0) )
			{
			flierr( z, y, j );
			errfl += 1;
			}
		}
	x = -x;
	x += 1.0;
	}
x00 = x00 + x00;
}
y = floor(0.0);
if( y != 0.0 )
	{
	flierr( 0.0, y, 57 );
	errfl += 1;
	}
y = floor(-0.0);
if( y != 0.0 )
	{
	flierr( -0.0, y, 58 );
	errfl += 1;
	}
y = floor(-1.0);
if( y != -1.0 )
	{
	flierr( -1.0, y, 59 );
	errfl += 1;
	}
y = floor(-0.1);
if( y != -1.0 )
	{
	flierr( -0.1, y, 60 );
	errfl += 1;
	}

if( errfl == 0 )
	printf( "No errors found in floor().\n" );

}


flierr( x, y, k )
double x, y;
int k;
{
printf( "Test %d: ", k+1 );
printf( "floor(%.15e) =?= %.15e\n", x, y );
}
