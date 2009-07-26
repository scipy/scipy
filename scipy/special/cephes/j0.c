/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        9400       7.0e-17     7.9e-18
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

#include "mconf.h"

#ifdef UNK
static double PP[7] = {
  7.96936729297347051624E-4,
  8.28352392107440799803E-2,
  1.23953371646414299388E0,
  5.44725003058768775090E0,
  8.74716500199817011941E0,
  5.30324038235394892183E0,
  9.99999999999999997821E-1,
};
static double PQ[7] = {
  9.24408810558863637013E-4,
  8.56288474354474431428E-2,
  1.25352743901058953537E0,
  5.47097740330417105182E0,
  8.76190883237069594232E0,
  5.30605288235394617618E0,
  1.00000000000000000218E0,
};
#endif
#ifdef DEC
static unsigned short PP[28] = {
0035520,0164604,0140733,0054470,
0037251,0122605,0115356,0107170,
0040236,0124412,0071500,0056303,
0040656,0047737,0045720,0045263,
0041013,0172143,0045004,0142103,
0040651,0132045,0026241,0026406,
0040200,0000000,0000000,0000000,
};
static unsigned short PQ[28] = {
0035562,0052006,0070034,0134666,
0037257,0057055,0055242,0123424,
0040240,0071626,0046630,0032371,
0040657,0011077,0032013,0012731,
0041014,0030307,0050331,0006414,
0040651,0145457,0065021,0150304,
0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short PP[28] = {
0x6b27,0x983b,0x1d30,0x3f4a,
0xd1cf,0xb35d,0x34b0,0x3fb5,
0x0b98,0x4e68,0xd521,0x3ff3,
0x0956,0xe97a,0xc9fb,0x4015,
0x9888,0x6940,0x7e8c,0x4021,
0x25a1,0xa594,0x3684,0x4015,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ[28] = {
0x9737,0xce03,0x4a80,0x3f4e,
0x54e3,0xab54,0xebc5,0x3fb5,
0x069f,0xc9b3,0x0e72,0x3ff4,
0x62bb,0xe681,0xe247,0x4015,
0x21a1,0xea1b,0x8618,0x4021,
0x3a19,0xed42,0x3965,0x4015,
0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short PP[28] = {
0x3f4a,0x1d30,0x983b,0x6b27,
0x3fb5,0x34b0,0xb35d,0xd1cf,
0x3ff3,0xd521,0x4e68,0x0b98,
0x4015,0xc9fb,0xe97a,0x0956,
0x4021,0x7e8c,0x6940,0x9888,
0x4015,0x3684,0xa594,0x25a1,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short PQ[28] = {
0x3f4e,0x4a80,0xce03,0x9737,
0x3fb5,0xebc5,0xab54,0x54e3,
0x3ff4,0x0e72,0xc9b3,0x069f,
0x4015,0xe247,0xe681,0x62bb,
0x4021,0x8618,0xea1b,0x21a1,
0x4015,0x3965,0xed42,0x3a19,
0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double QP[8] = {
-1.13663838898469149931E-2,
-1.28252718670509318512E0,
-1.95539544257735972385E1,
-9.32060152123768231369E1,
-1.77681167980488050595E2,
-1.47077505154951170175E2,
-5.14105326766599330220E1,
-6.05014350600728481186E0,
};
static double QQ[7] = {
/*  1.00000000000000000000E0,*/
  6.43178256118178023184E1,
  8.56430025976980587198E2,
  3.88240183605401609683E3,
  7.24046774195652478189E3,
  5.93072701187316984827E3,
  2.06209331660327847417E3,
  2.42005740240291393179E2,
};
#endif
#ifdef DEC
static unsigned short QP[32] = {
0136472,0035021,0142451,0141115,
0140244,0024731,0150620,0105642,
0141234,0067177,0124161,0060141,
0141672,0064572,0151557,0043036,
0142061,0127141,0003127,0043517,
0142023,0011727,0060271,0144544,
0141515,0122142,0126620,0143150,
0140701,0115306,0106715,0007344,
};
static unsigned short QQ[28] = {
/*0040200,0000000,0000000,0000000,*/
0041600,0121272,0004741,0026544,
0042526,0015605,0105654,0161771,
0043162,0123155,0165644,0062645,
0043342,0041675,0167576,0130756,
0043271,0052720,0165631,0154214,
0043000,0160576,0034614,0172024,
0042162,0000570,0030500,0051235,
};
#endif
#ifdef IBMPC
static unsigned short QP[32] = {
0x384a,0x38a5,0x4742,0xbf87,
0x1174,0x3a32,0x853b,0xbff4,
0x2c0c,0xf50e,0x8dcf,0xc033,
0xe8c4,0x5a6d,0x4d2f,0xc057,
0xe8ea,0x20ca,0x35cc,0xc066,
0x392d,0xec17,0x627a,0xc062,
0x18cd,0x55b2,0xb48c,0xc049,
0xa1dd,0xd1b9,0x3358,0xc018,
};
static unsigned short QQ[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x25ac,0x413c,0x1457,0x4050,
0x9c7f,0xb175,0xc370,0x408a,
0x8cb5,0xbd74,0x54cd,0x40ae,
0xd63e,0xbdef,0x4877,0x40bc,
0x3b11,0x1d73,0x2aba,0x40b7,
0x9e82,0xc731,0x1c2f,0x40a0,
0x0a54,0x0628,0x402f,0x406e,
};
#endif
#ifdef MIEEE
static unsigned short QP[32] = {
0xbf87,0x4742,0x38a5,0x384a,
0xbff4,0x853b,0x3a32,0x1174,
0xc033,0x8dcf,0xf50e,0x2c0c,
0xc057,0x4d2f,0x5a6d,0xe8c4,
0xc066,0x35cc,0x20ca,0xe8ea,
0xc062,0x627a,0xec17,0x392d,
0xc049,0xb48c,0x55b2,0x18cd,
0xc018,0x3358,0xd1b9,0xa1dd,
};
static unsigned short QQ[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4050,0x1457,0x413c,0x25ac,
0x408a,0xc370,0xb175,0x9c7f,
0x40ae,0x54cd,0xbd74,0x8cb5,
0x40bc,0x4877,0xbdef,0xd63e,
0x40b7,0x2aba,0x1d73,0x3b11,
0x40a0,0x1c2f,0xc731,0x9e82,
0x406e,0x402f,0x0628,0x0a54,
};
#endif


#ifdef UNK
static double YP[8] = {
 1.55924367855235737965E4,
-1.46639295903971606143E7,
 5.43526477051876500413E9,
-9.82136065717911466409E11,
 8.75906394395366999549E13,
-3.46628303384729719441E15,
 4.42733268572569800351E16,
-1.84950800436986690637E16,
};
static double YQ[7] = {
/* 1.00000000000000000000E0,*/
 1.04128353664259848412E3,
 6.26107330137134956842E5,
 2.68919633393814121987E8,
 8.64002487103935000337E10,
 2.02979612750105546709E13,
 3.17157752842975028269E15,
 2.50596256172653059228E17,
};
#endif
#ifdef DEC
static unsigned short YP[32] = {
0043563,0120677,0042264,0046166,
0146137,0140371,0113444,0042260,
0050241,0175707,0100502,0063344,
0152144,0125737,0007265,0164526,
0053637,0051621,0163035,0060546,
0155105,0004416,0107306,0060023,
0056035,0045133,0030132,0000024,
0155603,0065132,0144061,0131732,
};
static unsigned short YQ[28] = {
/*0040200,0000000,0000000,0000000,*/
0042602,0024422,0135557,0162663,
0045030,0155665,0044075,0160135,
0047200,0035432,0105446,0104005,
0051240,0167331,0056063,0022743,
0053223,0127746,0025764,0012160,
0055064,0044206,0177532,0145545,
0056536,0111375,0163715,0127201,
};
#endif
#ifdef IBMPC
static unsigned short YP[32] = {
0x898f,0xe896,0x7437,0x40ce,
0x8896,0x32e4,0xf81f,0xc16b,
0x4cdd,0xf028,0x3f78,0x41f4,
0xbd2b,0xe1d6,0x957b,0xc26c,
0xac2d,0x3cc3,0xea72,0x42d3,
0xcc02,0xd1d8,0xa121,0xc328,
0x4003,0x660b,0xa94b,0x4363,
0x367b,0x5906,0x6d4b,0xc350,
};
static unsigned short YQ[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xfcb6,0x576d,0x4522,0x4090,
0xbc0c,0xa907,0x1b76,0x4123,
0xd101,0x5164,0x0763,0x41b0,
0x64bc,0x2b86,0x1ddb,0x4234,
0x828e,0xc57e,0x75fc,0x42b2,
0x596d,0xdfeb,0x8910,0x4326,
0xb5d0,0xbcf9,0xd25f,0x438b,
};
#endif
#ifdef MIEEE
static unsigned short YP[32] = {
0x40ce,0x7437,0xe896,0x898f,
0xc16b,0xf81f,0x32e4,0x8896,
0x41f4,0x3f78,0xf028,0x4cdd,
0xc26c,0x957b,0xe1d6,0xbd2b,
0x42d3,0xea72,0x3cc3,0xac2d,
0xc328,0xa121,0xd1d8,0xcc02,
0x4363,0xa94b,0x660b,0x4003,
0xc350,0x6d4b,0x5906,0x367b,
};
static unsigned short YQ[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4090,0x4522,0x576d,0xfcb6,
0x4123,0x1b76,0xa907,0xbc0c,
0x41b0,0x0763,0x5164,0xd101,
0x4234,0x1ddb,0x2b86,0x64bc,
0x42b2,0x75fc,0xc57e,0x828e,
0x4326,0x8910,0xdfeb,0x596d,
0x438b,0xd25f,0xbcf9,0xb5d0,
};
#endif

#ifdef UNK
/*  5.783185962946784521175995758455807035071 */
static double DR1 = 5.78318596294678452118E0;
/* 30.47126234366208639907816317502275584842 */
static double DR2 = 3.04712623436620863991E1;
#endif

#ifdef DEC
static unsigned short R1[] = {0040671,0007734,0001061,0056734};
#define DR1 *(double *)R1
static unsigned short R2[] = {0041363,0142445,0030416,0165567};
#define DR2 *(double *)R2
#endif

#ifdef IBMPC
static unsigned short R1[] = {0x2bbb,0x8046,0x21fb,0x4017};
#define DR1 *(double *)R1
static unsigned short R2[] = {0xdd6f,0xa621,0x78a4,0x403e};
#define DR2 *(double *)R2
#endif

#ifdef MIEEE
static unsigned short R1[] = {0x4017,0x21fb,0x8046,0x2bbb};
#define DR1 *(double *)R1
static unsigned short R2[] = {0x403e,0x78a4,0xa621,0xdd6f};
#define DR2 *(double *)R2
#endif

#ifdef UNK
static double RP[4] = {
-4.79443220978201773821E9,
 1.95617491946556577543E12,
-2.49248344360967716204E14,
 9.70862251047306323952E15,
};
static double RQ[8] = {
/* 1.00000000000000000000E0,*/
 4.99563147152651017219E2,
 1.73785401676374683123E5,
 4.84409658339962045305E7,
 1.11855537045356834862E10,
 2.11277520115489217587E12,
 3.10518229857422583814E14,
 3.18121955943204943306E16,
 1.71086294081043136091E18,
};
#endif
#ifdef DEC
static unsigned short RP[16] = {
0150216,0161235,0064344,0014450,
0052343,0135216,0035624,0144153,
0154142,0130247,0003310,0003667,
0055411,0173703,0047772,0176635,
};
static unsigned short RQ[32] = {
/*0040200,0000000,0000000,0000000,*/
0042371,0144025,0032265,0136137,
0044451,0133131,0132420,0151466,
0046470,0144641,0072540,0030636,
0050446,0126600,0045042,0044243,
0052365,0172633,0110301,0071063,
0054215,0032424,0062272,0043513,
0055742,0005013,0171731,0072335,
0057275,0170646,0036663,0013134,
};
#endif
#ifdef IBMPC
static unsigned short RP[16] = {
0x8325,0xad1c,0xdc53,0xc1f1,
0x990d,0xc772,0x7751,0x427c,
0x00f7,0xe0d9,0x5614,0xc2ec,
0x5fb4,0x69ff,0x3ef8,0x4341,
};
static unsigned short RQ[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xb78c,0xa696,0x3902,0x407f,
0x1a67,0x36a2,0x36cb,0x4105,
0x0634,0x2eac,0x1934,0x4187,
0x4914,0x0944,0xd5b0,0x4204,
0x2e46,0x7218,0xbeb3,0x427e,
0x48e9,0x8c97,0xa6a2,0x42f1,
0x2e9c,0x7e7b,0x4141,0x435c,
0x62cc,0xc7b6,0xbe34,0x43b7,
};
#endif
#ifdef MIEEE
static unsigned short RP[16] = {
0xc1f1,0xdc53,0xad1c,0x8325,
0x427c,0x7751,0xc772,0x990d,
0xc2ec,0x5614,0xe0d9,0x00f7,
0x4341,0x3ef8,0x69ff,0x5fb4,
};
static unsigned short RQ[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x407f,0x3902,0xa696,0xb78c,
0x4105,0x36cb,0x36a2,0x1a67,
0x4187,0x1934,0x2eac,0x0634,
0x4204,0xd5b0,0x0944,0x4914,
0x427e,0xbeb3,0x7218,0x2e46,
0x42f1,0xa6a2,0x8c97,0x48e9,
0x435c,0x4141,0x7e7b,0x2e9c,
0x43b7,0xbe34,0xc7b6,0x62cc,
};
#endif

extern double TWOOPI, SQ2OPI;

double j0(x)
double x;
{
double w, z, p, q, xn;

if( x < 0 )
	x = -x;

if( x <= 5.0 )
	{
	z = x * x;
	if( x < 1.0e-5 )
		return( 1.0 - z/4.0 );

	p = (z - DR1) * (z - DR2);
	p = p * polevl( z, RP, 3)/p1evl( z, RQ, 8 );
	return( p );
	}

w = 5.0/x;
q = 25.0/(x*x);
p = polevl( q, PP, 6)/polevl( q, PQ, 6 );
q = polevl( q, QP, 7)/p1evl( q, QQ, 7 );
xn = x - NPY_PI_4;
p = p * cos(xn) - w * q * sin(xn);
return( p * SQ2OPI / sqrt(x) );
}

/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / PI
 * = 0.073804295108687225.
 */

/*
#define NPY_PI_4 .78539816339744830962
#define SQ2OPI .79788456080286535588
*/

double y0(x)
double x;
{
double w, z, p, q, xn;

if( x <= 5.0 )
	{
	if (x == 0.0) {
		mtherr("y0", SING);
		return -NPY_INFINITY;
	} else if (x < 0.0) {
		mtherr("y0", DOMAIN);
		return NPY_NAN;
	}
	z = x * x;
	w = polevl( z, YP, 7) / p1evl( z, YQ, 7 );
	w += TWOOPI * log(x) * j0(x);
	return( w );
	}

w = 5.0/x;
z = 25.0 / (x * x);
p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
xn = x - NPY_PI_4;
p = p * sin(xn) + w * q * cos(xn);
return( p * SQ2OPI / sqrt(x) );
}
