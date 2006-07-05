/*							k1.c
 *
 *	Modified Bessel function, third kind, order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1();
 *
 * y = k1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the modified Bessel function of the third kind
 * of order one of the argument.
 *
 * The range is partitioned into the two intervals [0,2] and
 * (2, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3300       8.9e-17     2.2e-17
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * k1 domain          x <= 0          MAXNUM
 *
 */
/*							k1e.c
 *
 *	Modified Bessel function, third kind, order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k1e();
 *
 * y = k1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order one of the argument:
 *
 *      k1e(x) = exp(x) * k1(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       7.8e-16     1.2e-16
 * See k1().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#include "mconf.h"

/* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
 * in the interval [0,2].
 * 
 * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
 */

#ifdef UNK
static double A[] =
{
-7.02386347938628759343E-18,
-2.42744985051936593393E-15,
-6.66690169419932900609E-13,
-1.41148839263352776110E-10,
-2.21338763073472585583E-8,
-2.43340614156596823496E-6,
-1.73028895751305206302E-4,
-6.97572385963986435018E-3,
-1.22611180822657148235E-1,
-3.53155960776544875667E-1,
 1.52530022733894777053E0
};
#endif

#ifdef DEC
static unsigned short A[] = {
0122001,0110501,0164746,0151255,
0124056,0165213,0150034,0147377,
0126073,0124026,0167207,0001044,
0130033,0030735,0141061,0033116,
0131676,0020350,0121341,0107175,
0133443,0046631,0062031,0070716,
0135065,0067427,0026435,0164022,
0136344,0112234,0165752,0006222,
0137373,0015622,0017016,0155636,
0137664,0150333,0125730,0067240,
0040303,0036411,0130200,0043120
};
#endif

#ifdef IBMPC
static unsigned short A[] = {
0xda56,0x3d3c,0x3228,0xbc60,
0x99e0,0x7a03,0xdd51,0xbce5,
0xe045,0xddd0,0x7502,0xbd67,
0x26ca,0xb846,0x663b,0xbde3,
0x31d0,0x145c,0xc41d,0xbe57,
0x2e3a,0x2c83,0x69b3,0xbec4,
0xbd02,0xe5a3,0xade2,0xbf26,
0x4192,0x9d7d,0x9293,0xbf7c,
0xdb74,0x43c1,0x6372,0xbfbf,
0x0dd4,0x757b,0x9a1b,0xbfd6,
0x08ca,0x3610,0x67a1,0x3ff8
};
#endif

#ifdef MIEEE
static unsigned short A[] = {
0xbc60,0x3228,0x3d3c,0xda56,
0xbce5,0xdd51,0x7a03,0x99e0,
0xbd67,0x7502,0xddd0,0xe045,
0xbde3,0x663b,0xb846,0x26ca,
0xbe57,0xc41d,0x145c,0x31d0,
0xbec4,0x69b3,0x2c83,0x2e3a,
0xbf26,0xade2,0xe5a3,0xbd02,
0xbf7c,0x9293,0x9d7d,0x4192,
0xbfbf,0x6372,0x43c1,0xdb74,
0xbfd6,0x9a1b,0x757b,0x0dd4,
0x3ff8,0x67a1,0x3610,0x08ca
};
#endif



/* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
 * in the interval [2,infinity].
 *
 * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
 */

#ifdef UNK
static double B[] =
{
-5.75674448366501715755E-18,
 1.79405087314755922667E-17,
-5.68946255844285935196E-17,
 1.83809354436663880070E-16,
-6.05704724837331885336E-16,
 2.03870316562433424052E-15,
-7.01983709041831346144E-15,
 2.47715442448130437068E-14,
-8.97670518232499435011E-14,
 3.34841966607842919884E-13,
-1.28917396095102890680E-12,
 5.13963967348173025100E-12,
-2.12996783842756842877E-11,
 9.21831518760500529508E-11,
-4.19035475934189648750E-10,
 2.01504975519703286596E-9,
-1.03457624656780970260E-8,
 5.74108412545004946722E-8,
-3.50196060308781257119E-7,
 2.40648494783721712015E-6,
-1.93619797416608296024E-5,
 1.95215518471351631108E-4,
-2.85781685962277938680E-3,
 1.03923736576817238437E-1,
 2.72062619048444266945E0
};
#endif

#ifdef DEC
static unsigned short B[] = {
0121724,0061352,0013041,0150076,
0022245,0074324,0016172,0173232,
0122603,0030250,0135670,0165221,
0023123,0165362,0023561,0060124,
0123456,0112436,0141654,0073623,
0024022,0163557,0077564,0006753,
0124374,0165221,0131014,0026524,
0024737,0017512,0144250,0175451,
0125312,0021456,0123136,0076633,
0025674,0077720,0020125,0102607,
0126265,0067543,0007744,0043701,
0026664,0152702,0033002,0074202,
0127273,0055234,0120016,0071733,
0027712,0133200,0042441,0075515,
0130346,0057000,0015456,0074470,
0031012,0074441,0051636,0111155,
0131461,0136444,0177417,0002101,
0032166,0111743,0032176,0021410,
0132674,0001224,0076555,0027060,
0033441,0077430,0135226,0106663,
0134242,0065610,0167155,0113447,
0035114,0131304,0043664,0102163,
0136073,0045065,0171465,0122123,
0037324,0152767,0147401,0017732,
0040456,0017275,0050061,0062120,
};
#endif

#ifdef IBMPC
static unsigned short B[] = {
0x3a08,0x42c4,0x8c5d,0xbc5a,
0x5ed3,0x838f,0xaf1a,0x3c74,
0x1d52,0x1777,0x6615,0xbc90,
0x2c0b,0x44ee,0x7d5e,0x3caa,
0x8ef2,0xd875,0xd2a3,0xbcc5,
0x81bd,0xefee,0x5ced,0x3ce2,
0x85ab,0x3641,0x9d52,0xbcff,
0x1f65,0x5915,0xe3e9,0x3d1b,
0xcfb3,0xd4cb,0x4465,0xbd39,
0xb0b1,0x040a,0x8ffa,0x3d57,
0x88f8,0x61fc,0xadec,0xbd76,
0x4f10,0x46c0,0x9ab8,0x3d96,
0xce7b,0x9401,0x6b53,0xbdb7,
0x2f6a,0x08a4,0x56d0,0x3dd9,
0xcf27,0x0365,0xcbc0,0xbdfc,
0xd24e,0x2a73,0x4f24,0x3e21,
0xe088,0x9fe1,0x37a4,0xbe46,
0xc461,0x668f,0xd27c,0x3e6e,
0xa5c6,0x8fad,0x8052,0xbe97,
0xd1b6,0x1752,0x2fe3,0x3ec4,
0xb2e5,0x1dcd,0x4d71,0xbef4,
0x908e,0x88f6,0x9658,0x3f29,
0xb48a,0xbe66,0x6946,0xbf67,
0x23fb,0xf9e0,0x9abe,0x3fba,
0x2c8a,0xaa06,0xc3d7,0x4005
};
#endif

#ifdef MIEEE
static unsigned short B[] = {
0xbc5a,0x8c5d,0x42c4,0x3a08,
0x3c74,0xaf1a,0x838f,0x5ed3,
0xbc90,0x6615,0x1777,0x1d52,
0x3caa,0x7d5e,0x44ee,0x2c0b,
0xbcc5,0xd2a3,0xd875,0x8ef2,
0x3ce2,0x5ced,0xefee,0x81bd,
0xbcff,0x9d52,0x3641,0x85ab,
0x3d1b,0xe3e9,0x5915,0x1f65,
0xbd39,0x4465,0xd4cb,0xcfb3,
0x3d57,0x8ffa,0x040a,0xb0b1,
0xbd76,0xadec,0x61fc,0x88f8,
0x3d96,0x9ab8,0x46c0,0x4f10,
0xbdb7,0x6b53,0x9401,0xce7b,
0x3dd9,0x56d0,0x08a4,0x2f6a,
0xbdfc,0xcbc0,0x0365,0xcf27,
0x3e21,0x4f24,0x2a73,0xd24e,
0xbe46,0x37a4,0x9fe1,0xe088,
0x3e6e,0xd27c,0x668f,0xc461,
0xbe97,0x8052,0x8fad,0xa5c6,
0x3ec4,0x2fe3,0x1752,0xd1b6,
0xbef4,0x4d71,0x1dcd,0xb2e5,
0x3f29,0x9658,0x88f6,0x908e,
0xbf67,0x6946,0xbe66,0xb48a,
0x3fba,0x9abe,0xf9e0,0x23fb,
0x4005,0xc3d7,0xaa06,0x2c8a
};
#endif

#ifdef ANSIPROT
extern double chbevl ( double, void *, int );
extern double exp ( double );
extern double i1 ( double );
extern double log ( double );
extern double sqrt ( double );
#else
double chbevl(), exp(), i1(), log(), sqrt();
#endif
extern double PI;
extern double MINLOG, MAXNUM;

double k1(x)
double x;
{
double y, z;

z = 0.5 * x;
if( z <= 0.0 )
	{
	mtherr( "k1", DOMAIN );
	return( MAXNUM );
	}

if( x <= 2.0 )
	{
	y = x * x - 2.0;
	y =  log(z) * i1(x)  +  chbevl( y, A, 11 ) / x;
	return( y );
	}

return(  exp(-x) * chbevl( 8.0/x - 2.0, B, 25 ) / sqrt(x) );
}




double k1e( x )
double x;
{
double y;

if( x <= 0.0 )
	{
	mtherr( "k1e", DOMAIN );
	return( MAXNUM );
	}

if( x <= 2.0 )
	{
	y = x * x - 2.0;
	y =  log( 0.5 * x ) * i1(x)  +  chbevl( y, A, 11 ) / x;
	return( y * exp(x) );
	}

return(  chbevl( 8.0/x - 2.0, B, 25 ) / sqrt(x) );
}
