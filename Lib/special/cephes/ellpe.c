/*							ellpe.c
 *
 *	Complete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m, y, ellpe();
 *
 * y = ellpe( m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *            pi/2
 *             -
 *            | |                 2
 * E(m)  =    |    sqrt( 1 - m sin t ) dt
 *          | |    
 *           -
 *            0
 *
 * Where m = 1 - m1, using the approximation
 *
 *      P(x)  -  x log x Q(x).
 *
 * Though there are no singularities, the argument m1 is used
 * internally rather than m for compatibility with ellpk().
 *
 * E(1) = 1; E(0) = pi/2.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0, 1       13000       3.1e-17     9.4e-18
 *    IEEE       0, 1       10000       2.1e-16     7.3e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpe domain      x<0, x>1            0.0
 *
 */

/*							ellpe.c		*/

/* Elliptic integral of second kind */

/*
Cephes Math Library, Release 2.1:  February, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140

Feb, 2002:  altered by Travis Oliphant
            so that it is called with argument m 
            (which gets immediately converted to m1 = 1-m)
*/

#include "mconf.h"

#ifdef UNK
static double P[] = {
  1.53552577301013293365E-4,
  2.50888492163602060990E-3,
  8.68786816565889628429E-3,
  1.07350949056076193403E-2,
  7.77395492516787092951E-3,
  7.58395289413514708519E-3,
  1.15688436810574127319E-2,
  2.18317996015557253103E-2,
  5.68051945617860553470E-2,
  4.43147180560990850618E-1,
  1.00000000000000000299E0
};
static double Q[] = {
  3.27954898576485872656E-5,
  1.00962792679356715133E-3,
  6.50609489976927491433E-3,
  1.68862163993311317300E-2,
  2.61769742454493659583E-2,
  3.34833904888224918614E-2,
  4.27180926518931511717E-2,
  5.85936634471101055642E-2,
  9.37499997197644278445E-2,
  2.49999999999888314361E-1
};
#endif

#ifdef DEC
static unsigned short P[] = {
0035041,0001364,0141572,0117555,
0036044,0066032,0130027,0033404,
0036416,0053617,0064456,0102632,
0036457,0161100,0061177,0122612,
0036376,0136251,0012403,0124162,
0036370,0101316,0151715,0131613,
0036475,0105477,0050317,0133272,
0036662,0154232,0024645,0171552,
0037150,0126220,0047054,0030064,
0037742,0162057,0167645,0165612,
0040200,0000000,0000000,0000000
};
static unsigned short Q[] = {
0034411,0106743,0115771,0055462,
0035604,0052575,0155171,0045540,
0036325,0030424,0064332,0167756,
0036612,0052366,0063006,0115175,
0036726,0070430,0004533,0124654,
0037011,0022741,0030675,0030711,
0037056,0174452,0127062,0132122,
0037157,0177750,0142041,0072523,
0037277,0177777,0173137,0002627,
0037577,0177777,0177777,0101101
};
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x53ee,0x986f,0x205e,0x3f24,
0xe6e0,0x5602,0x8d83,0x3f64,
0xd0b3,0xed25,0xcaf1,0x3f81,
0xf4b1,0x0c4f,0xfc48,0x3f85,
0x750e,0x22a0,0xd795,0x3f7f,
0xb671,0xda79,0x1059,0x3f7f,
0xf6d7,0xea19,0xb167,0x3f87,
0xbe6d,0x4534,0x5b13,0x3f96,
0x8607,0x09c5,0x1592,0x3fad,
0xbd71,0xfdf4,0x5c85,0x3fdc,
0x0000,0x0000,0x0000,0x3ff0
};
static unsigned short Q[] = {
0x2b66,0x737f,0x31bc,0x3f01,
0x296c,0xbb4f,0x8aaf,0x3f50,
0x5dfe,0x8d1b,0xa622,0x3f7a,
0xd350,0xccc0,0x4a9e,0x3f91,
0x7535,0x012b,0xce23,0x3f9a,
0xa639,0x2637,0x24bc,0x3fa1,
0x568a,0x55c6,0xdf25,0x3fa5,
0x2eaa,0x1884,0xfffd,0x3fad,
0xe0b3,0xfecb,0xffff,0x3fb7,
0xf048,0xffff,0xffff,0x3fcf
};
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3f24,0x205e,0x986f,0x53ee,
0x3f64,0x8d83,0x5602,0xe6e0,
0x3f81,0xcaf1,0xed25,0xd0b3,
0x3f85,0xfc48,0x0c4f,0xf4b1,
0x3f7f,0xd795,0x22a0,0x750e,
0x3f7f,0x1059,0xda79,0xb671,
0x3f87,0xb167,0xea19,0xf6d7,
0x3f96,0x5b13,0x4534,0xbe6d,
0x3fad,0x1592,0x09c5,0x8607,
0x3fdc,0x5c85,0xfdf4,0xbd71,
0x3ff0,0x0000,0x0000,0x0000
};
static unsigned short Q[] = {
0x3f01,0x31bc,0x737f,0x2b66,
0x3f50,0x8aaf,0xbb4f,0x296c,
0x3f7a,0xa622,0x8d1b,0x5dfe,
0x3f91,0x4a9e,0xccc0,0xd350,
0x3f9a,0xce23,0x012b,0x7535,
0x3fa1,0x24bc,0x2637,0xa639,
0x3fa5,0xdf25,0x55c6,0x568a,
0x3fad,0xfffd,0x1884,0x2eaa,
0x3fb7,0xffff,0xfecb,0xe0b3,
0x3fcf,0xffff,0xffff,0xf048
};
#endif

#ifndef ANSIPROT
double polevl(), log();
#endif

extern double NAN;

double ellpe(x)
double x;
{
x = 1.0-x;
if( (x <= 0.0) || (x > 1.0) )
	{
	if( x == 0.0 )
		return( 1.0 );
	mtherr( "ellpe", DOMAIN );
	return( NAN );
	}
return( polevl(x,P,10) - log(x) * (x * polevl(x,Q,9)) );
}
