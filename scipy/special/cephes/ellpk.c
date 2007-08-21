/*							ellpk.c
 *
 *	Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * double m, y, ellpk();
 *
 * y = ellpk( m ); 
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used internally rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC        0,1        16000       3.5e-17     1.1e-17
 *    IEEE       0,1        30000       2.5e-16     6.8e-17
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpk domain       x<0, x>1           0.0
 *
 */

/*							ellpk.c */


/*
Cephes Math Library, Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140

Feb, 2002:  altered by Travis Oliphant 
            so that it is called with argument m 
            (which gets immediately converted to m1 = 1-m)
*/

#include "mconf.h"

#ifdef DEC
static unsigned short P[] =
{
0035020,0127576,0040430,0051544,
0036025,0070136,0042703,0153716,
0036402,0122614,0062555,0077777,
0036441,0102130,0072334,0025172,
0036341,0043320,0117242,0172076,
0036312,0146456,0077242,0154141,
0036420,0003467,0013727,0035407,
0036564,0137263,0110651,0020237,
0036775,0001330,0144056,0020305,
0037305,0144137,0157521,0141734,
0040261,0071027,0173721,0147572
};
static unsigned short Q[] =
{
0034366,0130371,0103453,0077633,
0035557,0122745,0173515,0113016,
0036302,0124470,0167304,0074473,
0036575,0132403,0117226,0117576,
0036703,0156271,0047124,0147733,
0036766,0137465,0002053,0157312,
0037031,0014423,0154274,0176515,
0037107,0177747,0143216,0016145,
0037217,0177777,0172621,0074000,
0037377,0177777,0177776,0156435,
0040000,0000000,0000000,0000000
};
static unsigned short ac1[] = {0040261,0071027,0173721,0147572};
#define C1 (*(double *)ac1)
#endif

#ifdef IBMPC
static unsigned short P[] =
{
0x0a6d,0xc823,0x15ef,0x3f22,
0x7afa,0xc8b8,0xae0b,0x3f62,
0xb000,0x8cad,0x54b1,0x3f80,
0x854f,0x0e9b,0x308b,0x3f84,
0x5e88,0x13d4,0x28da,0x3f7c,
0x5b0c,0xcfd4,0x59a5,0x3f79,
0xe761,0xe2fa,0x00e6,0x3f82,
0x2414,0x7235,0x97d6,0x3f8e,
0xc419,0x1905,0xa05b,0x3f9f,
0x387c,0xfbea,0xb90b,0x3fb8,
0x39ef,0xfefa,0x2e42,0x3ff6
};
static unsigned short Q[] =
{
0x6ff3,0x30e5,0xd61f,0x3efe,
0xb2c2,0xbee9,0xf4bc,0x3f4d,
0x8f27,0x1dd8,0x5527,0x3f78,
0xd3f0,0x73d2,0xb6a0,0x3f8f,
0x99fb,0x29ca,0x7b97,0x3f98,
0x7bd9,0xa085,0xd7e6,0x3f9e,
0x9faa,0x7b17,0x2322,0x3fa3,
0xc38d,0xf8d1,0xfffc,0x3fa8,
0x2f00,0xfeb2,0xffff,0x3fb1,
0xdba4,0xffff,0xffff,0x3fbf,
0x0000,0x0000,0x0000,0x3fe0
};
static unsigned short ac1[] = {0x39ef,0xfefa,0x2e42,0x3ff6};
#define C1 (*(double *)ac1)
#endif

#ifdef MIEEE
static unsigned short P[] =
{
0x3f22,0x15ef,0xc823,0x0a6d,
0x3f62,0xae0b,0xc8b8,0x7afa,
0x3f80,0x54b1,0x8cad,0xb000,
0x3f84,0x308b,0x0e9b,0x854f,
0x3f7c,0x28da,0x13d4,0x5e88,
0x3f79,0x59a5,0xcfd4,0x5b0c,
0x3f82,0x00e6,0xe2fa,0xe761,
0x3f8e,0x97d6,0x7235,0x2414,
0x3f9f,0xa05b,0x1905,0xc419,
0x3fb8,0xb90b,0xfbea,0x387c,
0x3ff6,0x2e42,0xfefa,0x39ef
};
static unsigned short Q[] =
{
0x3efe,0xd61f,0x30e5,0x6ff3,
0x3f4d,0xf4bc,0xbee9,0xb2c2,
0x3f78,0x5527,0x1dd8,0x8f27,
0x3f8f,0xb6a0,0x73d2,0xd3f0,
0x3f98,0x7b97,0x29ca,0x99fb,
0x3f9e,0xd7e6,0xa085,0x7bd9,
0x3fa3,0x2322,0x7b17,0x9faa,
0x3fa8,0xfffc,0xf8d1,0xc38d,
0x3fb1,0xffff,0xfeb2,0x2f00,
0x3fbf,0xffff,0xffff,0xdba4,
0x3fe0,0x0000,0x0000,0x0000
};
static unsigned short ac1[] = {
0x3ff6,0x2e42,0xfefa,0x39ef
};
#define C1 (*(double *)ac1)
#endif

#ifdef UNK
static double P[] =
{
 1.37982864606273237150E-4,
 2.28025724005875567385E-3,
 7.97404013220415179367E-3,
 9.85821379021226008714E-3,
 6.87489687449949877925E-3,
 6.18901033637687613229E-3,
 8.79078273952743772254E-3,
 1.49380448916805252718E-2,
 3.08851465246711995998E-2,
 9.65735902811690126535E-2,
 1.38629436111989062502E0
};

static double Q[] =
{
 2.94078955048598507511E-5,
 9.14184723865917226571E-4,
 5.94058303753167793257E-3,
 1.54850516649762399335E-2,
 2.39089602715924892727E-2,
 3.01204715227604046988E-2,
 3.73774314173823228969E-2,
 4.88280347570998239232E-2,
 7.03124996963957469739E-2,
 1.24999999999870820058E-1,
 4.99999999999999999821E-1
};
static double C1 = 1.3862943611198906188E0; /* log(4) */
#endif

#ifndef ANSIPROT
double polevl(), p1evl(), log();
#endif
extern double MACHEP, MAXNUM, NAN;

double ellpk(x)    /* Changed to use m argument rather than m1 = 1-m */
double x;
{

x = 1.0-x;  
if( (x < 0.0) || (x > 1.0) )
	{
	mtherr( "ellpk", DOMAIN );
	return( NAN );
	}

if( x > MACHEP )
	{
	return( polevl(x,P,10) - log(x) * polevl(x,Q,10) );
	}
else
	{
	if( x == 0.0 )
		{
		mtherr( "ellpk", SING );
		return( MAXNUM );
		}
	else
		{
		return( C1 - 0.5 * log(x) );
		}
	}
}
