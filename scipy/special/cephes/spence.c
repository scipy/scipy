/*							spence.c
 *
 *	Dilogarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, spence();
 *
 * y = spence( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes the integral
 *
 *                    x
 *                    -
 *                   | | log t
 * spence(x)  =  -   |   ----- dt
 *                 | |   t - 1
 *                  -
 *                  1
 *
 * for x >= 0.  A rational approximation gives the integral in
 * the interval (0.5, 1.5).  Transformation formulas for 1/x
 * and 1-x are employed outside the basic expansion range.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,4         30000       3.9e-15     5.4e-16
 *    DEC       0,4          3000       2.5e-16     4.5e-17
 *
 *
 */

/*							spence.c */


/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1985, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"

#ifdef UNK
static double A[8] = {
  4.65128586073990045278E-5,
  7.31589045238094711071E-3,
  1.33847639578309018650E-1,
  8.79691311754530315341E-1,
  2.71149851196553469920E0,
  4.25697156008121755724E0,
  3.29771340985225106936E0,
  1.00000000000000000126E0,
};
static double B[8] = {
  6.90990488912553276999E-4,
  2.54043763932544379113E-2,
  2.82974860602568089943E-1,
  1.41172597751831069617E0,
  3.63800533345137075418E0,
  5.03278880143316990390E0,
  3.54771340985225096217E0,
  9.99999999999999998740E-1,
};
#endif
#ifdef DEC
static unsigned short A[32] = {
0034503,0013315,0034120,0157771,
0036357,0135043,0016766,0150637,
0037411,0007533,0005212,0161475,
0040141,0031563,0023217,0120331,
0040455,0104461,0007002,0155522,
0040610,0034434,0065721,0120465,
0040523,0006674,0105671,0054427,
0040200,0000000,0000000,0000000,
};
static unsigned short B[32] = {
0035465,0021626,0032367,0144157,
0036720,0016326,0134431,0000406,
0037620,0161024,0133701,0120766,
0040264,0131557,0152055,0064512,
0040550,0152424,0051166,0034272,
0040641,0006233,0014672,0111572,
0040543,0006674,0105671,0054425,
0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short A[32] = {
0x1bff,0xa70a,0x62d9,0x3f08,
0xda34,0x63be,0xf744,0x3f7d,
0x5c68,0x6151,0x21eb,0x3fc1,
0xf41b,0x64d1,0x266e,0x3fec,
0x5b6a,0x21c0,0xb126,0x4005,
0x3427,0x8d7a,0x0723,0x4011,
0x2b23,0x9177,0x61b7,0x400a,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short B[32] = {
0xf90e,0xc69e,0xa472,0x3f46,
0x2021,0xd723,0x039a,0x3f9a,
0x343f,0x96f8,0x1c42,0x3fd2,
0xad29,0xfa85,0x966d,0x3ff6,
0xc717,0x8a4e,0x1aa2,0x400d,
0x526f,0x6337,0x2193,0x4014,
0x2b23,0x9177,0x61b7,0x400c,
0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short A[32] = {
0x3f08,0x62d9,0xa70a,0x1bff,
0x3f7d,0xf744,0x63be,0xda34,
0x3fc1,0x21eb,0x6151,0x5c68,
0x3fec,0x266e,0x64d1,0xf41b,
0x4005,0xb126,0x21c0,0x5b6a,
0x4011,0x0723,0x8d7a,0x3427,
0x400a,0x61b7,0x9177,0x2b23,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short B[32] = {
0x3f46,0xa472,0xc69e,0xf90e,
0x3f9a,0x039a,0xd723,0x2021,
0x3fd2,0x1c42,0x96f8,0x343f,
0x3ff6,0x966d,0xfa85,0xad29,
0x400d,0x1aa2,0x8a4e,0xc717,
0x4014,0x2193,0x6337,0x526f,
0x400c,0x61b7,0x9177,0x2b23,
0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifndef ANSIPROT
double fabs(), log(), polevl();
#endif
extern double PI, MACHEP, NAN;

double spence(x)
double x;
{
double w, y, z;
int flag;

if( x < 0.0 )
	{
	mtherr( "spence", DOMAIN );
	return(NAN);
	}

if( x == 1.0 )
	return( 0.0 );

if( x == 0.0 )
	return( PI*PI/6.0 );

flag = 0;

if( x > 2.0 )
	{
	x = 1.0/x;
	flag |= 2;
	}

if( x > 1.5 )
	{
	w = (1.0/x) - 1.0;
	flag |= 2;
	}

else if( x < 0.5 )
	{
	w = -x;
	flag |= 1;
	}

else
	w = x - 1.0;


y = -w * polevl( w, A, 7) / polevl( w, B, 7 );

if( flag & 1 )
	y = (PI * PI)/6.0  - log(x) * log(1.0-x) - y;

if( flag & 2 )
	{
	z = log(x);
	y = -0.5 * z * z  -  y;
	}

return( y );
}
