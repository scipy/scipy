/*							k0.c
 *
 *	Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0();
 *
 * y = k0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at 2000 random points between 0 and 8.  Peak absolute
 * error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3100       1.3e-16     2.1e-17
 *    IEEE      0, 30       30000       1.2e-15     1.6e-16
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 *  K0 domain          x <= 0          MAXNUM
 *
 */
/*							k0e()
 *
 *	Modified Bessel function, third kind, order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, k0e();
 *
 * y = k0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order zero of the argument.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.4e-15     1.4e-16
 * See k0().
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

#include "mconf.h"

/* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
 * in the interval [0,2].  The odd order coefficients are all
 * zero; only the even order coefficients are listed.
 * 
 * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
 */

#ifdef UNK
static double A[] =
{
 1.37446543561352307156E-16,
 4.25981614279661018399E-14,
 1.03496952576338420167E-11,
 1.90451637722020886025E-9,
 2.53479107902614945675E-7,
 2.28621210311945178607E-5,
 1.26461541144692592338E-3,
 3.59799365153615016266E-2,
 3.44289899924628486886E-1,
-5.35327393233902768720E-1
};
#endif

#ifdef DEC
static unsigned short A[] = {
0023036,0073417,0032477,0165673,
0025077,0154126,0016046,0012517,
0027066,0011342,0035211,0005041,
0031002,0160233,0037454,0050224,
0032610,0012747,0037712,0173741,
0034277,0144007,0172147,0162375,
0035645,0140563,0125431,0165626,
0037023,0057662,0125124,0102051,
0037660,0043304,0004411,0166707,
0140011,0005467,0047227,0130370
};
#endif

#ifdef IBMPC
static unsigned short A[] = {
0xfd77,0xe6a7,0xcee1,0x3ca3,
0xc2aa,0xc384,0xfb0a,0x3d27,
0x2144,0x4751,0xc25c,0x3da6,
0x8a13,0x67e5,0x5c13,0x3e20,
0x5efc,0xe7f9,0x02bc,0x3e91,
0xfca0,0xfe8c,0xf900,0x3ef7,
0x3d73,0x7563,0xb82e,0x3f54,
0x9085,0x554a,0x6bf6,0x3fa2,
0x3db9,0x8121,0x08d8,0x3fd6,
0xf61f,0xe9d2,0x2166,0xbfe1
};
#endif

#ifdef MIEEE
static unsigned short A[] = {
0x3ca3,0xcee1,0xe6a7,0xfd77,
0x3d27,0xfb0a,0xc384,0xc2aa,
0x3da6,0xc25c,0x4751,0x2144,
0x3e20,0x5c13,0x67e5,0x8a13,
0x3e91,0x02bc,0xe7f9,0x5efc,
0x3ef7,0xf900,0xfe8c,0xfca0,
0x3f54,0xb82e,0x7563,0x3d73,
0x3fa2,0x6bf6,0x554a,0x9085,
0x3fd6,0x08d8,0x8121,0x3db9,
0xbfe1,0x2166,0xe9d2,0xf61f
};
#endif



/* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
 * in the inverted interval [2,infinity].
 * 
 * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
 */

#ifdef UNK
static double B[] = {
 5.30043377268626276149E-18,
-1.64758043015242134646E-17,
 5.21039150503902756861E-17,
-1.67823109680541210385E-16,
 5.51205597852431940784E-16,
-1.84859337734377901440E-15,
 6.34007647740507060557E-15,
-2.22751332699166985548E-14,
 8.03289077536357521100E-14,
-2.98009692317273043925E-13,
 1.14034058820847496303E-12,
-4.51459788337394416547E-12,
 1.85594911495471785253E-11,
-7.95748924447710747776E-11,
 3.57739728140030116597E-10,
-1.69753450938905987466E-9,
 8.57403401741422608519E-9,
-4.66048989768794782956E-8,
 2.76681363944501510342E-7,
-1.83175552271911948767E-6,
 1.39498137188764993662E-5,
-1.28495495816278026384E-4,
 1.56988388573005337491E-3,
-3.14481013119645005427E-2,
 2.44030308206595545468E0
};
#endif

#ifdef DEC
static unsigned short B[] = {
0021703,0106456,0076144,0173406,
0122227,0173144,0116011,0030033,
0022560,0044562,0006506,0067642,
0123101,0076243,0123273,0131013,
0023436,0157713,0056243,0141331,
0124005,0032207,0063726,0164664,
0024344,0066342,0051756,0162300,
0124710,0121365,0154053,0077022,
0025264,0161166,0066246,0077420,
0125647,0141671,0006443,0103212,
0026240,0076431,0077147,0160445,
0126636,0153741,0174002,0105031,
0027243,0040102,0035375,0163073,
0127656,0176256,0113476,0044653,
0030304,0125544,0006377,0130104,
0130751,0047257,0110537,0127324,
0031423,0046400,0014772,0012164,
0132110,0025240,0155247,0112570,
0032624,0105314,0007437,0021574,
0133365,0155243,0174306,0116506,
0034152,0004776,0061643,0102504,
0135006,0136277,0036104,0175023,
0035715,0142217,0162474,0115022,
0137000,0147671,0065177,0134356,
0040434,0026754,0175163,0044070
};
#endif

#ifdef IBMPC
static unsigned short B[] = {
0x9ee1,0xcf8c,0x71a5,0x3c58,
0x2603,0x9381,0xfecc,0xbc72,
0xcdf4,0x41a8,0x092e,0x3c8e,
0x7641,0x74d7,0x2f94,0xbca8,
0x785b,0x6b94,0xdbf9,0x3cc3,
0xdd36,0xecfa,0xa690,0xbce0,
0xdc98,0x4a7d,0x8d9c,0x3cfc,
0x6fc2,0xbb05,0x145e,0xbd19,
0xcfe2,0xcd94,0x9c4e,0x3d36,
0x70d1,0x21a4,0xf877,0xbd54,
0xfc25,0x2fcc,0x0fa3,0x3d74,
0x5143,0x3f00,0xdafc,0xbd93,
0xbcc7,0x475f,0x6808,0x3db4,
0xc935,0xd2e7,0xdf95,0xbdd5,
0xf608,0x819f,0x956c,0x3df8,
0xf5db,0xf22b,0x29d5,0xbe1d,
0x428e,0x033f,0x69a0,0x3e42,
0xf2af,0x1b54,0x0554,0xbe69,
0xe46f,0x81e3,0x9159,0x3e92,
0xd3a9,0x7f18,0xbb54,0xbebe,
0x70a9,0xcc74,0x413f,0x3eed,
0x9f42,0xe788,0xd797,0xbf20,
0x9342,0xfca7,0xb891,0x3f59,
0xf71e,0x2d4f,0x19f7,0xbfa0,
0x6907,0x9f4e,0x85bd,0x4003
};
#endif

#ifdef MIEEE
static unsigned short B[] = {
0x3c58,0x71a5,0xcf8c,0x9ee1,
0xbc72,0xfecc,0x9381,0x2603,
0x3c8e,0x092e,0x41a8,0xcdf4,
0xbca8,0x2f94,0x74d7,0x7641,
0x3cc3,0xdbf9,0x6b94,0x785b,
0xbce0,0xa690,0xecfa,0xdd36,
0x3cfc,0x8d9c,0x4a7d,0xdc98,
0xbd19,0x145e,0xbb05,0x6fc2,
0x3d36,0x9c4e,0xcd94,0xcfe2,
0xbd54,0xf877,0x21a4,0x70d1,
0x3d74,0x0fa3,0x2fcc,0xfc25,
0xbd93,0xdafc,0x3f00,0x5143,
0x3db4,0x6808,0x475f,0xbcc7,
0xbdd5,0xdf95,0xd2e7,0xc935,
0x3df8,0x956c,0x819f,0xf608,
0xbe1d,0x29d5,0xf22b,0xf5db,
0x3e42,0x69a0,0x033f,0x428e,
0xbe69,0x0554,0x1b54,0xf2af,
0x3e92,0x9159,0x81e3,0xe46f,
0xbebe,0xbb54,0x7f18,0xd3a9,
0x3eed,0x413f,0xcc74,0x70a9,
0xbf20,0xd797,0xe788,0x9f42,
0x3f59,0xb891,0xfca7,0x9342,
0xbfa0,0x19f7,0x2d4f,0xf71e,
0x4003,0x85bd,0x9f4e,0x6907
};
#endif

/*							k0.c	*/
#ifdef ANSIPROT 
extern double chbevl ( double, void *, int );
extern double exp ( double );
extern double i0 ( double );
extern double log ( double );
extern double sqrt ( double );
#else
double chbevl(), exp(), i0(), log(), sqrt();
#endif
extern double PI;
extern double MAXNUM;

double k0(x)
double x;
{
double y, z;

if( x <= 0.0 )
	{
	mtherr( "k0", DOMAIN );
	return( MAXNUM );
	}

if( x <= 2.0 )
	{
	y = x * x - 2.0;
	y = chbevl( y, A, 10 ) - log( 0.5 * x ) * i0(x);
	return( y );
	}
z = 8.0/x - 2.0;
y = exp(-x) * chbevl( z, B, 25 ) / sqrt(x);
return(y);
}




double k0e( x )
double x;
{
double y;

if( x <= 0.0 )
	{
	mtherr( "k0e", DOMAIN );
	return( MAXNUM );
	}

if( x <= 2.0 )
	{
	y = x * x - 2.0;
	y = chbevl( y, A, 10 ) - log( 0.5 * x ) * i0(x);
	return( y * exp(x) );
	}

y = chbevl( 8.0/x - 2.0, B, 25 ) / sqrt(x);
return(y);
}
