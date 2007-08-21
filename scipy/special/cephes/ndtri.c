/*							ndtri.c
 *
 *	Inverse of Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, ndtri();
 *
 * x = ndtri( y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the argument, x, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to x) is equal to y.
 *
 *
 * For small arguments 0 < y < exp(-2), the program computes
 * z = sqrt( -2.0 * log(y) );  then the approximation is
 * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
 * There are two rational functions P/Q, one for 0 < y < exp(-32)
 * and the other for y up to exp(-2).  For larger arguments,
 * w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain        # trials      peak         rms
 *    DEC      0.125, 1         5500       9.5e-17     2.1e-17
 *    DEC      6e-39, 0.135     3500       5.7e-17     1.3e-17
 *    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
 *    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition    value returned
 * ndtri domain       x <= 0        -MAXNUM
 * ndtri domain       x >= 1         MAXNUM
 *
 */


/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"
extern double MAXNUM;

#ifdef UNK
/* sqrt(2pi) */
static double s2pi = 2.50662827463100050242E0;
#endif

#ifdef DEC
static unsigned short s2p[] = {0040440,0066230,0177661,0034055};
#define s2pi *(double *)s2p
#endif

#ifdef IBMPC
static unsigned short s2p[] = {0x2706,0x1ff6,0x0d93,0x4004};
#define s2pi *(double *)s2p
#endif

#ifdef MIEEE
static unsigned short s2p[] = {
0x4004,0x0d93,0x1ff6,0x2706
};
#define s2pi *(double *)s2p
#endif

/* approximation for 0 <= |y - 0.5| <= 3/8 */
#ifdef UNK
static double P0[5] = {
-5.99633501014107895267E1,
 9.80010754185999661536E1,
-5.66762857469070293439E1,
 1.39312609387279679503E1,
-1.23916583867381258016E0,
};
static double Q0[8] = {
/* 1.00000000000000000000E0,*/
 1.95448858338141759834E0,
 4.67627912898881538453E0,
 8.63602421390890590575E1,
-2.25462687854119370527E2,
 2.00260212380060660359E2,
-8.20372256168333339912E1,
 1.59056225126211695515E1,
-1.18331621121330003142E0,
};
#endif
#ifdef DEC
static unsigned short P0[20] = {
0141557,0155170,0071360,0120550,
0041704,0000214,0172417,0067307,
0141542,0132204,0040066,0156723,
0041136,0163161,0157276,0007747,
0140236,0116374,0073666,0051764,
};
static unsigned short Q0[32] = {
/*0040200,0000000,0000000,0000000,*/
0040372,0026256,0110403,0123707,
0040625,0122024,0020277,0026661,
0041654,0134161,0124134,0007244,
0142141,0073162,0133021,0131371,
0042110,0041235,0043516,0057767,
0141644,0011417,0036155,0137305,
0041176,0076556,0004043,0125430,
0140227,0073347,0152776,0067251,
};
#endif
#ifdef IBMPC
static unsigned short P0[20] = {
0x142d,0x0e5e,0xfb4f,0xc04d,
0xedd9,0x9ea1,0x8011,0x4058,
0xdbba,0x8806,0x5690,0xc04c,
0xc1fd,0x3bd7,0xdcce,0x402b,
0xca7e,0x8ef6,0xd39f,0xbff3,
};
static unsigned short Q0[36] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x74f9,0xd220,0x4595,0x3fff,
0xe5b6,0x8417,0xb482,0x4012,
0x81d4,0x350b,0x970e,0x4055,
0x365f,0x56c2,0x2ece,0xc06c,
0xcbff,0xa8e9,0x0853,0x4069,
0xb7d9,0xe78d,0x8261,0xc054,
0x7563,0xc104,0xcfad,0x402f,
0xcdd5,0xfabf,0xeedc,0xbff2,
};
#endif
#ifdef MIEEE
static unsigned short P0[20] = {
0xc04d,0xfb4f,0x0e5e,0x142d,
0x4058,0x8011,0x9ea1,0xedd9,
0xc04c,0x5690,0x8806,0xdbba,
0x402b,0xdcce,0x3bd7,0xc1fd,
0xbff3,0xd39f,0x8ef6,0xca7e,
};
static unsigned short Q0[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x3fff,0x4595,0xd220,0x74f9,
0x4012,0xb482,0x8417,0xe5b6,
0x4055,0x970e,0x350b,0x81d4,
0xc06c,0x2ece,0x56c2,0x365f,
0x4069,0x0853,0xa8e9,0xcbff,
0xc054,0x8261,0xe78d,0xb7d9,
0x402f,0xcfad,0xc104,0x7563,
0xbff2,0xeedc,0xfabf,0xcdd5,
};
#endif


/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
 * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
 */
#ifdef UNK
static double P1[9] = {
 4.05544892305962419923E0,
 3.15251094599893866154E1,
 5.71628192246421288162E1,
 4.40805073893200834700E1,
 1.46849561928858024014E1,
 2.18663306850790267539E0,
-1.40256079171354495875E-1,
-3.50424626827848203418E-2,
-8.57456785154685413611E-4,
};
static double Q1[8] = {
/*  1.00000000000000000000E0,*/
 1.57799883256466749731E1,
 4.53907635128879210584E1,
 4.13172038254672030440E1,
 1.50425385692907503408E1,
 2.50464946208309415979E0,
-1.42182922854787788574E-1,
-3.80806407691578277194E-2,
-9.33259480895457427372E-4,
};
#endif
#ifdef DEC
static unsigned short P1[36] = {
0040601,0143074,0150744,0073326,
0041374,0031554,0113253,0146016,
0041544,0123272,0012463,0176771,
0041460,0051160,0103560,0156511,
0041152,0172624,0117772,0030755,
0040413,0170713,0151545,0176413,
0137417,0117512,0022154,0131671,
0137017,0104257,0071432,0007072,
0135540,0143363,0063137,0036166,
};
static unsigned short Q1[32] = {
/*0040200,0000000,0000000,0000000,*/
0041174,0075325,0004736,0120326,
0041465,0110044,0047561,0045567,
0041445,0042321,0012142,0030340,
0041160,0127074,0166076,0141051,
0040440,0046055,0040745,0150400,
0137421,0114146,0067330,0010621,
0137033,0175162,0025555,0114351,
0135564,0122773,0145750,0030357,
};
#endif
#ifdef IBMPC
static unsigned short P1[36] = {
0x8edb,0x9a3c,0x38c7,0x4010,
0x7982,0x92d5,0x866d,0x403f,
0x7fbf,0x42a6,0x94d7,0x404c,
0x1ba9,0x10ee,0x0a4e,0x4046,
0x463e,0x93ff,0x5eb2,0x402d,
0xbfa1,0x7a6c,0x7e39,0x4001,
0x9677,0x448d,0xf3e9,0xbfc1,
0x41c7,0xee63,0xf115,0xbfa1,
0xe78f,0x6ccb,0x18de,0xbf4c,
};
static unsigned short Q1[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xd41b,0xa13b,0x8f5a,0x402f,
0x296f,0x89ee,0xb204,0x4046,
0x461c,0x228c,0xa89a,0x4044,
0xd845,0x9d87,0x15c7,0x402e,
0xba20,0xa83c,0x0985,0x4004,
0x0232,0xcddb,0x330c,0xbfc2,
0xb31d,0x456d,0x7f4e,0xbfa3,
0x061e,0x797d,0x94bf,0xbf4e,
};
#endif
#ifdef MIEEE
static unsigned short P1[36] = {
0x4010,0x38c7,0x9a3c,0x8edb,
0x403f,0x866d,0x92d5,0x7982,
0x404c,0x94d7,0x42a6,0x7fbf,
0x4046,0x0a4e,0x10ee,0x1ba9,
0x402d,0x5eb2,0x93ff,0x463e,
0x4001,0x7e39,0x7a6c,0xbfa1,
0xbfc1,0xf3e9,0x448d,0x9677,
0xbfa1,0xf115,0xee63,0x41c7,
0xbf4c,0x18de,0x6ccb,0xe78f,
};
static unsigned short Q1[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x402f,0x8f5a,0xa13b,0xd41b,
0x4046,0xb204,0x89ee,0x296f,
0x4044,0xa89a,0x228c,0x461c,
0x402e,0x15c7,0x9d87,0xd845,
0x4004,0x0985,0xa83c,0xba20,
0xbfc2,0x330c,0xcddb,0x0232,
0xbfa3,0x7f4e,0x456d,0xb31d,
0xbf4e,0x94bf,0x797d,0x061e,
};
#endif

/* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
 * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
 */

#ifdef UNK
static double P2[9] = {
  3.23774891776946035970E0,
  6.91522889068984211695E0,
  3.93881025292474443415E0,
  1.33303460815807542389E0,
  2.01485389549179081538E-1,
  1.23716634817820021358E-2,
  3.01581553508235416007E-4,
  2.65806974686737550832E-6,
  6.23974539184983293730E-9,
};
static double Q2[8] = {
/*  1.00000000000000000000E0,*/
  6.02427039364742014255E0,
  3.67983563856160859403E0,
  1.37702099489081330271E0,
  2.16236993594496635890E-1,
  1.34204006088543189037E-2,
  3.28014464682127739104E-4,
  2.89247864745380683936E-6,
  6.79019408009981274425E-9,
};
#endif
#ifdef DEC
static unsigned short P2[36] = {
0040517,0033507,0036236,0125641,
0040735,0044616,0014473,0140133,
0040574,0012567,0114535,0102541,
0040252,0120340,0143474,0150135,
0037516,0051057,0115361,0031211,
0036512,0131204,0101511,0125144,
0035236,0016627,0043160,0140216,
0033462,0060512,0060141,0010641,
0031326,0062541,0101304,0077706,
};
static unsigned short Q2[32] = {
/*0040200,0000000,0000000,0000000,*/
0040700,0143322,0132137,0040501,
0040553,0101155,0053221,0140257,
0040260,0041071,0052573,0010004,
0037535,0066472,0177261,0162330,
0036533,0160475,0066666,0036132,
0035253,0174533,0027771,0044027,
0033502,0016147,0117666,0063671,
0031351,0047455,0141663,0054751,
};
#endif
#ifdef IBMPC
static unsigned short P2[36] = {
0xd574,0xe793,0xe6e8,0x4009,
0x780b,0xc327,0xa931,0x401b,
0xb0ac,0xf32b,0x82ae,0x400f,
0x9a0c,0x18e7,0x541c,0x3ff5,
0x2651,0xf35e,0xca45,0x3fc9,
0x354d,0x9069,0x5650,0x3f89,
0x1812,0xe8ce,0xc3b2,0x3f33,
0x2234,0x4c0c,0x4c29,0x3ec6,
0x8ff9,0x3058,0xccac,0x3e3a,
};
static unsigned short Q2[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xe828,0x568b,0x18da,0x4018,
0x3816,0xaad2,0x704d,0x400d,
0x6200,0x2aaf,0x0847,0x3ff6,
0x3c9b,0x5fd6,0xada7,0x3fcb,
0xc78b,0xadb6,0x7c27,0x3f8b,
0x2903,0x65ff,0x7f2b,0x3f35,
0xccf7,0xf3f6,0x438c,0x3ec8,
0x6b3d,0xb876,0x29e5,0x3e3d,
};
#endif
#ifdef MIEEE
static unsigned short P2[36] = {
0x4009,0xe6e8,0xe793,0xd574,
0x401b,0xa931,0xc327,0x780b,
0x400f,0x82ae,0xf32b,0xb0ac,
0x3ff5,0x541c,0x18e7,0x9a0c,
0x3fc9,0xca45,0xf35e,0x2651,
0x3f89,0x5650,0x9069,0x354d,
0x3f33,0xc3b2,0xe8ce,0x1812,
0x3ec6,0x4c29,0x4c0c,0x2234,
0x3e3a,0xccac,0x3058,0x8ff9,
};
static unsigned short Q2[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4018,0x18da,0x568b,0xe828,
0x400d,0x704d,0xaad2,0x3816,
0x3ff6,0x0847,0x2aaf,0x6200,
0x3fcb,0xada7,0x5fd6,0x3c9b,
0x3f8b,0x7c27,0xadb6,0xc78b,
0x3f35,0x7f2b,0x65ff,0x2903,
0x3ec8,0x438c,0xf3f6,0xccf7,
0x3e3d,0x29e5,0xb876,0x6b3d,
};
#endif

#ifndef ANSIPROT
double polevl(), p1evl(), log(), sqrt();
#endif

double ndtri(y0)
double y0;
{
double x, y, z, y2, x0, x1;
int code;

if( y0 <= 0.0 )
	{
	mtherr( "ndtri", DOMAIN );
	return( -MAXNUM );
	}
if( y0 >= 1.0 )
	{
	mtherr( "ndtri", DOMAIN );
	return( MAXNUM );
	}
code = 1;
y = y0;
if( y > (1.0 - 0.13533528323661269189) ) /* 0.135... = exp(-2) */
	{
	y = 1.0 - y;
	code = 0;
	}

if( y > 0.13533528323661269189 )
	{
	y = y - 0.5;
	y2 = y * y;
	x = y + y * (y2 * polevl( y2, P0, 4)/p1evl( y2, Q0, 8 ));
	x = x * s2pi; 
	return(x);
	}

x = sqrt( -2.0 * log(y) );
x0 = x - log(x)/x;

z = 1.0/x;
if( x < 8.0 ) /* y > exp(-32) = 1.2664165549e-14 */
	x1 = z * polevl( z, P1, 8 )/p1evl( z, Q1, 8 );
else
	x1 = z * polevl( z, P2, 8 )/p1evl( z, Q2, 8 );
x = x0 - x1;
if( code != 0 )
	x = -x;
return( x );
}
