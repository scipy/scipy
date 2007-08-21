/*							sincos.c
 *
 *	Circular sine and cosine of argument in degrees
 *	Table lookup and interpolation algorithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, sine, cosine, flg, sincos();
 *
 * sincos( x, &sine, &cosine, flg );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns both the sine and the cosine of the argument x.
 * Several different compile time options and minimax
 * approximations are supplied to permit tailoring the
 * tradeoff between computation speed and accuracy.
 * 
 * Since range reduction is time consuming, the reduction
 * of x modulo 360 degrees is also made optional.
 *
 * sin(i) is internally tabulated for 0 <= i <= 90 degrees.
 * Approximation polynomials, ranging from linear interpolation
 * to cubics in (x-i)**2, compute the sine and cosine
 * of the residual x-i which is between -0.5 and +0.5 degree.
 * In the case of the high accuracy options, the residual
 * and the tabulated values are combined using the trigonometry
 * formulas for sin(A+B) and cos(A+B).
 *
 * Compile time options are supplied for 5, 11, or 17 decimal
 * relative accuracy (ACC5, ACC11, ACC17 respectively).
 * A subroutine flag argument "flg" chooses betwen this
 * accuracy and table lookup only (peak absolute error
 * = 0.0087).
 *
 * If the argument flg = 1, then the tabulated value is
 * returned for the nearest whole number of degrees. The
 * approximation polynomials are not computed.  At
 * x = 0.5 deg, the absolute error is then sin(0.5) = 0.0087.
 *
 * An intermediate speed and precision can be obtained using
 * the compile time option LINTERP and flg = 1.  This yields
 * a linear interpolation using a slope estimated from the sine
 * or cosine at the nearest integer argument.  The peak absolute
 * error with this option is 3.8e-5.  Relative error at small
 * angles is about 1e-5.
 *
 * If flg = 0, then the approximation polynomials are computed
 * and applied.
 *
 *
 *
 * SPEED:
 *
 * Relative speed comparisons follow for 6MHz IBM AT clone
 * and Microsoft C version 4.0.  These figures include
 * software overhead of do loop and function calls.
 * Since system hardware and software vary widely, the
 * numbers should be taken as representative only.
 *
 *			flg=0	flg=0	flg=1	flg=1
 *			ACC11	ACC5	LINTERP	Lookup only
 * In-line 8087 (/FPi)
 * sin(), cos()		1.0	1.0	1.0	1.0
 *
 * In-line 8087 (/FPi)
 * sincos()		1.1	1.4	1.9	3.0
 *
 * Software (/FPa)
 * sin(), cos()		0.19	0.19	0.19	0.19
 *
 * Software (/FPa)
 * sincos()		0.39	0.50	0.73	1.7
 *
 *
 *
 * ACCURACY:
 *
 * The accurate approximations are designed with a relative error
 * criterion.  The absolute error is greatest at x = 0.5 degree.
 * It decreases from a local maximum at i+0.5 degrees to full
 * machine precision at each integer i degrees.  With the
 * ACC5 option, the relative error of 6.3e-6 is equivalent to
 * an absolute angular error of 0.01 arc second in the argument
 * at x = i+0.5 degrees.  For small angles < 0.5 deg, the ACC5
 * accuracy is 6.3e-6 (.00063%) of reading; i.e., the absolute
 * error decreases in proportion to the argument.  This is true
 * for both the sine and cosine approximations, since the latter
 * is for the function 1 - cos(x).
 *
 * If absolute error is of most concern, use the compile time
 * option ABSERR to obtain an absolute error of 2.7e-8 for ACC5
 * precision.  This is about half the absolute error of the
 * relative precision option.  In this case the relative error
 * for small angles will increase to 9.5e-6 -- a reasonable
 * tradeoff.
 */


#include "mconf.h"

/* Define one of the following to be 1:
 */
#define ACC5 1
#define ACC11 0
#define ACC17 0

/* Option for linear interpolation when flg = 1
 */
#define LINTERP 1

/* Option for absolute error criterion
 */
#define ABSERR 1

/* Option to include modulo 360 function:
 */
#define MOD360 0

/*
Cephes Math Library Release 2.1
Copyright 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


/* Table of sin(i degrees)
 * for 0 <= i <= 90
 */
static double sintbl[92] = {
  0.00000000000000000000E0,
  1.74524064372835128194E-2,
  3.48994967025009716460E-2,
  5.23359562429438327221E-2,
  6.97564737441253007760E-2,
  8.71557427476581735581E-2,
  1.04528463267653471400E-1,
  1.21869343405147481113E-1,
  1.39173100960065444112E-1,
  1.56434465040230869010E-1,
  1.73648177666930348852E-1,
  1.90808995376544812405E-1,
  2.07911690817759337102E-1,
  2.24951054343864998051E-1,
  2.41921895599667722560E-1,
  2.58819045102520762349E-1,
  2.75637355816999185650E-1,
  2.92371704722736728097E-1,
  3.09016994374947424102E-1,
  3.25568154457156668714E-1,
  3.42020143325668733044E-1,
  3.58367949545300273484E-1,
  3.74606593415912035415E-1,
  3.90731128489273755062E-1,
  4.06736643075800207754E-1,
  4.22618261740699436187E-1,
  4.38371146789077417453E-1,
  4.53990499739546791560E-1,
  4.69471562785890775959E-1,
  4.84809620246337029075E-1,
  5.00000000000000000000E-1,
  5.15038074910054210082E-1,
  5.29919264233204954047E-1,
  5.44639035015027082224E-1,
  5.59192903470746830160E-1,
  5.73576436351046096108E-1,
  5.87785252292473129169E-1,
  6.01815023152048279918E-1,
  6.15661475325658279669E-1,
  6.29320391049837452706E-1,
  6.42787609686539326323E-1,
  6.56059028990507284782E-1,
  6.69130606358858213826E-1,
  6.81998360062498500442E-1,
  6.94658370458997286656E-1,
  7.07106781186547524401E-1,
  7.19339800338651139356E-1,
  7.31353701619170483288E-1,
  7.43144825477394235015E-1,
  7.54709580222771997943E-1,
  7.66044443118978035202E-1,
  7.77145961456970879980E-1,
  7.88010753606721956694E-1,
  7.98635510047292846284E-1,
  8.09016994374947424102E-1,
  8.19152044288991789684E-1,
  8.29037572555041692006E-1,
  8.38670567945424029638E-1,
  8.48048096156425970386E-1,
  8.57167300702112287465E-1,
  8.66025403784438646764E-1,
  8.74619707139395800285E-1,
  8.82947592858926942032E-1,
  8.91006524188367862360E-1,
  8.98794046299166992782E-1,
  9.06307787036649963243E-1,
  9.13545457642600895502E-1,
  9.20504853452440327397E-1,
  9.27183854566787400806E-1,
  9.33580426497201748990E-1,
  9.39692620785908384054E-1,
  9.45518575599316810348E-1,
  9.51056516295153572116E-1,
  9.56304755963035481339E-1,
  9.61261695938318861916E-1,
  9.65925826289068286750E-1,
  9.70295726275996472306E-1,
  9.74370064785235228540E-1,
  9.78147600733805637929E-1,
  9.81627183447663953497E-1,
  9.84807753012208059367E-1,
  9.87688340595137726190E-1,
  9.90268068741570315084E-1,
  9.92546151641322034980E-1,
  9.94521895368273336923E-1,
  9.96194698091745532295E-1,
  9.97564050259824247613E-1,
  9.98629534754573873784E-1,
  9.99390827019095730006E-1,
  9.99847695156391239157E-1,
  1.00000000000000000000E0,
  9.99847695156391239157E-1,
};

#ifndef ANSIPROT
double floor();
#else
extern void sincos ( double x, double *s, double *c, int flg );
#endif

void
sincos(x, s, c, flg)
double x;
double *s, *c;
int flg;
{
int ix, ssign, csign, xsign;
double y, z, sx, sz, cx, cz;

/* Make argument nonnegative.
 */
xsign = 1;
if( x < 0.0 )
	{
	xsign = -1;
	x = -x;
	}


#if MOD360
x = x  -  360.0 * floor( x/360.0 );
#endif

/* Find nearest integer to x.
 * Note there should be a domain error test here,
 * but this is omitted to gain speed.
 */
ix = x + 0.5;
z = x - ix;		/* the residual */

/* Look up the sine and cosine of the integer.
 */
if( ix <= 180 )
	{
	ssign = 1;
	csign = 1;
	}
else
	{
	ssign = -1;
	csign = -1;
	ix -= 180;
	}

if( ix > 90 )
	{
	csign = -csign;
	ix = 180 - ix;
	}

sx = sintbl[ix];
if( ssign < 0 )
	sx = -sx;
cx = sintbl[ 90-ix ];
if( csign < 0 )
	cx = -cx;

/* If the flag argument is set, then just return
 * the tabulated values for arg to the nearest whole degree.
 */
if( flg )
	{
#if LINTERP
	y = sx + 1.74531263774940077459e-2 * z * cx;
	cx -= 1.74531263774940077459e-2 * z * sx;
	sx = y;
#endif
	if( xsign < 0 )
		sx = -sx;
	*s = sx;	/* sine */
	*c = cx;	/* cosine */
	return;
	}


if( ssign < 0 )
	sx = -sx;
if( csign < 0 )
	cx = -cx;

/* Find sine and cosine
 * of the residual angle between -0.5 and +0.5 degree.
 */
#if ACC5
#if ABSERR
/* absolute error = 2.769e-8: */
sz = 1.74531263774940077459e-2 * z;
/* absolute error = 4.146e-11: */
cz = 1.0 - 1.52307909153324666207e-4 * z * z;
#else
/* relative error = 6.346e-6: */
sz = 1.74531817576426662296e-2 * z;
/* relative error = 3.173e-6: */
cz = 1.0 - 1.52308226602566149927e-4 * z * z;
#endif
#else
y = z * z;
#endif


#if ACC11
sz = ( -8.86092781698004819918e-7 * y
      + 1.74532925198378577601e-2     ) * z;

cz = 1.0 - ( -3.86631403698859047896e-9 * y
            + 1.52308709893047593702e-4     ) * y;
#endif


#if ACC17
sz = ((  1.34959795251974073996e-11 * y
       - 8.86096155697856783296e-7     ) * y
       + 1.74532925199432957214e-2          ) * z;

cz = 1.0 - ((  3.92582397764340914444e-14 * y
             - 3.86632385155548605680e-9     ) * y
             + 1.52308709893354299569e-4          ) * y;
#endif


/* Combine the tabulated part and the calculated part
 * by trigonometry.
 */
y = sx * cz  +  cx * sz;
if( xsign < 0 )
	y = - y;
*s = y; /* sine */

*c = cx * cz  -  sx * sz; /* cosine */
}
