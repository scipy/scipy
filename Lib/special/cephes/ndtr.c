/*							ndtr.c
 *
 *	Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, ndtr();
 *
 * y = ndtr( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x:
 *
 *                            x
 *                             -
 *                   1        | |          2
 *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
 *               sqrt(2pi)  | |
 *                           -
 *                          -inf.
 *
 *             =  ( 1 + erf(z) ) / 2
 *             =  erfc(z) / 2
 *
 * where z = x/sqrt(2). Computation is via the functions
 * erf and erfc.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC      -13,0         8000       2.1e-15     4.8e-16
 *    IEEE     -13,0        30000       3.4e-14     6.7e-15
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition         value returned
 * erfc underflow    x > 37.519379347       0.0
 *
 */
/*							erf.c
 *
 *	Error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erf();
 *
 * y = erf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The integral is
 *
 *                           x 
 *                            -
 *                 2         | |          2
 *   erf(x)  =  --------     |    exp( - t  ) dt.
 *              sqrt(pi)   | |
 *                          -
 *                           0
 *
 * The magnitude of x is limited to 9.231948545 for DEC
 * arithmetic; 1 or -1 is returned outside this range.
 *
 * For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise
 * erf(x) = 1 - erfc(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,1         14000       4.7e-17     1.5e-17
 *    IEEE      0,1         30000       3.7e-16     1.0e-16
 *
 */
/*							erfc.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, erfc();
 *
 * y = erfc( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf. 
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 *
 * For small x, erfc(x) = 1 - erf(x); otherwise rational
 * approximations are computed.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 9.2319   12000       5.1e-16     1.2e-16
 *    IEEE      0,26.6417   30000       5.7e-14     1.5e-14
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition              value returned
 * erfc underflow    x > 9.231948545 (DEC)       0.0
 *
 *
 */


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


#include "mconf.h"

extern double SQRTH, NAN;
extern double MAXLOG;

#ifdef UNK
static double P[] = {
 2.46196981473530512524E-10,
 5.64189564831068821977E-1,
 7.46321056442269912687E0,
 4.86371970985681366614E1,
 1.96520832956077098242E2,
 5.26445194995477358631E2,
 9.34528527171957607540E2,
 1.02755188689515710272E3,
 5.57535335369399327526E2
};
static double Q[] = {
/* 1.00000000000000000000E0,*/
 1.32281951154744992508E1,
 8.67072140885989742329E1,
 3.54937778887819891062E2,
 9.75708501743205489753E2,
 1.82390916687909736289E3,
 2.24633760818710981792E3,
 1.65666309194161350182E3,
 5.57535340817727675546E2
};
static double R[] = {
 5.64189583547755073984E-1,
 1.27536670759978104416E0,
 5.01905042251180477414E0,
 6.16021097993053585195E0,
 7.40974269950448939160E0,
 2.97886665372100240670E0
};
static double S[] = {
/* 1.00000000000000000000E0,*/
 2.26052863220117276590E0,
 9.39603524938001434673E0,
 1.20489539808096656605E1,
 1.70814450747565897222E1,
 9.60896809063285878198E0,
 3.36907645100081516050E0
};
static double T[] = {
 9.60497373987051638749E0,
 9.00260197203842689217E1,
 2.23200534594684319226E3,
 7.00332514112805075473E3,
 5.55923013010394962768E4
};
static double U[] = {
/* 1.00000000000000000000E0,*/
 3.35617141647503099647E1,
 5.21357949780152679795E2,
 4.59432382970980127987E3,
 2.26290000613890934246E4,
 4.92673942608635921086E4
};

#define UTHRESH 37.519379347
#endif

#ifdef DEC
static unsigned short P[] = {
0030207,0054445,0011173,0021706,
0040020,0067272,0030661,0122075,
0040756,0151236,0173053,0067042,
0041502,0106175,0062555,0151457,
0042104,0102525,0047401,0003667,
0042403,0116176,0011446,0075303,
0042551,0120723,0061641,0123275,
0042600,0070651,0007264,0134516,
0042413,0061102,0167507,0176625
};
static unsigned short Q[] = {
/*0040200,0000000,0000000,0000000,*/
0041123,0123257,0165741,0017142,
0041655,0065027,0173413,0115450,
0042261,0074011,0021573,0004150,
0042563,0166530,0013662,0007200,
0042743,0176427,0162443,0105214,
0043014,0062546,0153727,0123772,
0042717,0012470,0006227,0067424,
0042413,0061103,0003042,0013254
};
static unsigned short R[] = {
0040020,0067272,0101024,0155421,
0040243,0037467,0056706,0026462,
0040640,0116017,0120665,0034315,
0040705,0020162,0143350,0060137,
0040755,0016234,0134304,0130157,
0040476,0122700,0051070,0015473
};
static unsigned short S[] = {
/*0040200,0000000,0000000,0000000,*/
0040420,0126200,0044276,0070413,
0041026,0053051,0007302,0063746,
0041100,0144203,0174051,0061151,
0041210,0123314,0126343,0177646,
0041031,0137125,0051431,0033011,
0040527,0117362,0152661,0066201
};
static unsigned short T[] = {
0041031,0126770,0170672,0166101,
0041664,0006522,0072360,0031770,
0043013,0100025,0162641,0126671,
0043332,0155231,0161627,0076200,
0044131,0024115,0021020,0117343
};
static unsigned short U[] = {
/*0040200,0000000,0000000,0000000,*/
0041406,0037461,0177575,0032714,
0042402,0053350,0123061,0153557,
0043217,0111227,0032007,0164217,
0043660,0145000,0004013,0160114,
0044100,0071544,0167107,0125471
};
#define UTHRESH 14.0
#endif

#ifdef IBMPC
static unsigned short P[] = {
0x6479,0xa24f,0xeb24,0x3df0,
0x3488,0x4636,0x0dd7,0x3fe2,
0x6dc4,0xdec5,0xda53,0x401d,
0xba66,0xacad,0x518f,0x4048,
0x20f7,0xa9e0,0x90aa,0x4068,
0xcf58,0xc264,0x738f,0x4080,
0x34d8,0x6c74,0x343a,0x408d,
0x972a,0x21d6,0x0e35,0x4090,
0xffb3,0x5de8,0x6c48,0x4081
};
static unsigned short Q[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x23cc,0xfd7c,0x74d5,0x402a,
0x7365,0xfee1,0xad42,0x4055,
0x610d,0x246f,0x2f01,0x4076,
0x41d0,0x02f6,0x7dab,0x408e,
0x7151,0xfca4,0x7fa2,0x409c,
0xf4ff,0xdafa,0x8cac,0x40a1,
0xede2,0x0192,0xe2a7,0x4099,
0x42d6,0x60c4,0x6c48,0x4081
};
static unsigned short R[] = {
0x9b62,0x5042,0x0dd7,0x3fe2,
0xc5a6,0xebb8,0x67e6,0x3ff4,
0xa71a,0xf436,0x1381,0x4014,
0x0c0c,0x58dd,0xa40e,0x4018,
0x960e,0x9718,0xa393,0x401d,
0x0367,0x0a47,0xd4b8,0x4007
};
static unsigned short S[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xce21,0x0917,0x1590,0x4002,
0x4cfd,0x21d8,0xcac5,0x4022,
0x2c4d,0x7f05,0x1910,0x4028,
0x7ff5,0x959c,0x14d9,0x4031,
0x26c1,0xaa63,0x37ca,0x4023,
0x2d90,0x5ab6,0xf3de,0x400a
};
static unsigned short T[] = {
0x5d88,0x1e37,0x35bf,0x4023,
0x067f,0x4e9e,0x81aa,0x4056,
0x35b7,0xbcb4,0x7002,0x40a1,
0xef90,0x3c72,0x5b53,0x40bb,
0x13dc,0xa442,0x2509,0x40eb
};
static unsigned short U[] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xa6ba,0x3fef,0xc7e6,0x4040,
0x3aee,0x14c6,0x4add,0x4080,
0xfd12,0xe680,0xf252,0x40b1,
0x7c0a,0x0101,0x1940,0x40d6,
0xf567,0x9dc8,0x0e6c,0x40e8
};
#define UTHRESH 37.519379347
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3df0,0xeb24,0xa24f,0x6479,
0x3fe2,0x0dd7,0x4636,0x3488,
0x401d,0xda53,0xdec5,0x6dc4,
0x4048,0x518f,0xacad,0xba66,
0x4068,0x90aa,0xa9e0,0x20f7,
0x4080,0x738f,0xc264,0xcf58,
0x408d,0x343a,0x6c74,0x34d8,
0x4090,0x0e35,0x21d6,0x972a,
0x4081,0x6c48,0x5de8,0xffb3
};
static unsigned short Q[] = {
0x402a,0x74d5,0xfd7c,0x23cc,
0x4055,0xad42,0xfee1,0x7365,
0x4076,0x2f01,0x246f,0x610d,
0x408e,0x7dab,0x02f6,0x41d0,
0x409c,0x7fa2,0xfca4,0x7151,
0x40a1,0x8cac,0xdafa,0xf4ff,
0x4099,0xe2a7,0x0192,0xede2,
0x4081,0x6c48,0x60c4,0x42d6
};
static unsigned short R[] = {
0x3fe2,0x0dd7,0x5042,0x9b62,
0x3ff4,0x67e6,0xebb8,0xc5a6,
0x4014,0x1381,0xf436,0xa71a,
0x4018,0xa40e,0x58dd,0x0c0c,
0x401d,0xa393,0x9718,0x960e,
0x4007,0xd4b8,0x0a47,0x0367
};
static unsigned short S[] = {
0x4002,0x1590,0x0917,0xce21,
0x4022,0xcac5,0x21d8,0x4cfd,
0x4028,0x1910,0x7f05,0x2c4d,
0x4031,0x14d9,0x959c,0x7ff5,
0x4023,0x37ca,0xaa63,0x26c1,
0x400a,0xf3de,0x5ab6,0x2d90
};
static unsigned short T[] = {
0x4023,0x35bf,0x1e37,0x5d88,
0x4056,0x81aa,0x4e9e,0x067f,
0x40a1,0x7002,0xbcb4,0x35b7,
0x40bb,0x5b53,0x3c72,0xef90,
0x40eb,0x2509,0xa442,0x13dc
};
static unsigned short U[] = {
0x4040,0xc7e6,0x3fef,0xa6ba,
0x4080,0x4add,0x14c6,0x3aee,
0x40b1,0xf252,0xe680,0xfd12,
0x40d6,0x1940,0x0101,0x7c0a,
0x40e8,0x0e6c,0x9dc8,0xf567
};
#define UTHRESH 37.519379347
#endif

#ifndef ANSIPROT
double polevl(), p1evl(), exp(), log(), fabs();
double erf(), erfc();
#endif

double ndtr(double a)
{
double x, y, z;

if (isnan(a)) {
  mtherr("ndtr", DOMAIN);
  return (NAN);
}

x = a * SQRTH;
z = fabs(x);

if( z < SQRTH )
	y = 0.5 + 0.5 * erf(x);

else
	{
	y = 0.5 * erfc(z);

	if( x > 0 )
		y = 1.0 - y;
	}

return(y);
}


double erfc(double a)
{
double p,q,x,y,z;

if (isnan(a)) {
  mtherr("erfc", DOMAIN);
  return (NAN);
}

if( a < 0.0 )
	x = -a;
else
	x = a;

if( x < 1.0 )
	return( 1.0 - erf(a) );

z = -a * a;

if( z < -MAXLOG )
	{
under:
	mtherr( "erfc", UNDERFLOW );
	if( a < 0 )
		return( 2.0 );
	else
		return( 0.0 );
	}

z = exp(z);

if( x < 8.0 )
	{
	p = polevl( x, P, 8 );
	q = p1evl( x, Q, 8 );
	}
else
	{
	p = polevl( x, R, 5 );
	q = p1evl( x, S, 6 );
	}
y = (z * p)/q;

if( a < 0 )
	y = 2.0 - y;

if( y == 0.0 )
	goto under;

return(y);
}



double erf(double x)
{
double y, z;

if (isnan(x)) {
  mtherr("erf", DOMAIN);
  return (NAN);
}

if( fabs(x) > 1.0 )
	return( 1.0 - erfc(x) );
z = x * x;

y = x * polevl( z, T, 4 ) / p1evl( z, U, 5 );
return( y );

}
