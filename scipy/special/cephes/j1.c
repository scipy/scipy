/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */


/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

/*
#define PIO4 .78539816339744830962
#define THPIO4 2.35619449019234492885
#define SQ2OPI .79788456080286535588
*/

#include "mconf.h"

#ifdef UNK
static double RP[4] = {
-8.99971225705559398224E8,
 4.52228297998194034323E11,
-7.27494245221818276015E13,
 3.68295732863852883286E15,
};
static double RQ[8] = {
/* 1.00000000000000000000E0,*/
 6.20836478118054335476E2,
 2.56987256757748830383E5,
 8.35146791431949253037E7,
 2.21511595479792499675E10,
 4.74914122079991414898E12,
 7.84369607876235854894E14,
 8.95222336184627338078E16,
 5.32278620332680085395E18,
};
#endif
#ifdef DEC
static unsigned short RP[16] = {
0147526,0110742,0063322,0077052,
0051722,0112720,0065034,0061530,
0153604,0052227,0033147,0105650,
0055121,0055025,0032276,0022015,
};
static unsigned short RQ[32] = {
/*0040200,0000000,0000000,0000000,*/
0042433,0032610,0155604,0033473,
0044572,0173320,0067270,0006616,
0046637,0045246,0162225,0006606,
0050645,0004773,0157577,0053004,
0052612,0033734,0001667,0176501,
0054462,0054121,0173147,0121367,
0056237,0002777,0121451,0176007,
0057623,0136253,0131601,0044710,
};
#endif
#ifdef IBMPC
static unsigned short RP[16] = {
0x4fc5,0x4cda,0xd23c,0xc1ca,
0x8c6b,0x0d43,0x52ba,0x425a,
0xf175,0xe6cc,0x8a92,0xc2d0,
0xc482,0xa697,0x2b42,0x432a,
};
static unsigned short RQ[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x86e7,0x1b70,0x66b1,0x4083,
0x01b2,0x0dd7,0x5eda,0x410f,
0xa1b1,0xdc92,0xe954,0x4193,
0xeac1,0x7bef,0xa13f,0x4214,
0xffa8,0x8076,0x46fb,0x4291,
0xf45f,0x3ecc,0x4b0a,0x4306,
0x3f81,0xf465,0xe0bf,0x4373,
0x2939,0x7670,0x7795,0x43d2,
};
#endif
#ifdef MIEEE
static unsigned short RP[16] = {
0xc1ca,0xd23c,0x4cda,0x4fc5,
0x425a,0x52ba,0x0d43,0x8c6b,
0xc2d0,0x8a92,0xe6cc,0xf175,
0x432a,0x2b42,0xa697,0xc482,
};
static unsigned short RQ[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4083,0x66b1,0x1b70,0x86e7,
0x410f,0x5eda,0x0dd7,0x01b2,
0x4193,0xe954,0xdc92,0xa1b1,
0x4214,0xa13f,0x7bef,0xeac1,
0x4291,0x46fb,0x8076,0xffa8,
0x4306,0x4b0a,0x3ecc,0xf45f,
0x4373,0xe0bf,0xf465,0x3f81,
0x43d2,0x7795,0x7670,0x2939,
};
#endif

#ifdef UNK
static double PP[7] = {
 7.62125616208173112003E-4,
 7.31397056940917570436E-2,
 1.12719608129684925192E0,
 5.11207951146807644818E0,
 8.42404590141772420927E0,
 5.21451598682361504063E0,
 1.00000000000000000254E0,
};
static double PQ[7] = {
 5.71323128072548699714E-4,
 6.88455908754495404082E-2,
 1.10514232634061696926E0,
 5.07386386128601488557E0,
 8.39985554327604159757E0,
 5.20982848682361821619E0,
 9.99999999999999997461E-1,
};
#endif
#ifdef DEC
static unsigned short PP[28] = {
0035507,0144542,0061543,0024326,
0037225,0145105,0017766,0022661,
0040220,0043766,0010254,0133255,
0040643,0113047,0142611,0151521,
0041006,0144344,0055351,0074261,
0040646,0156520,0120574,0006416,
0040200,0000000,0000000,0000000,
};
static unsigned short PQ[28] = {
0035425,0142330,0115041,0165514,
0037214,0177352,0145105,0052026,
0040215,0072515,0141207,0073255,
0040642,0056427,0137222,0106405,
0041006,0062716,0166427,0165450,
0040646,0133352,0035425,0123304,
0040200,0000000,0000000,0000000,
};
#endif
#ifdef IBMPC
static unsigned short PP[28] = {
0x651b,0x4c6c,0xf92c,0x3f48,
0xc4b6,0xa3fe,0xb948,0x3fb2,
0x96d6,0xc215,0x08fe,0x3ff2,
0x3a6a,0xf8b1,0x72c4,0x4014,
0x2f16,0x8b5d,0xd91c,0x4020,
0x81a2,0x142f,0xdbaa,0x4014,
0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ[28] = {
0x3d69,0x1344,0xb89b,0x3f42,
0xaa83,0x5948,0x9fdd,0x3fb1,
0xeed6,0xb850,0xaea9,0x3ff1,
0x51a1,0xf7d2,0x4ba2,0x4014,
0xfd65,0xdda2,0xccb9,0x4020,
0xb4d9,0x4762,0xd6dd,0x4014,
0x0000,0x0000,0x0000,0x3ff0,
};
#endif
#ifdef MIEEE
static unsigned short PP[28] = {
0x3f48,0xf92c,0x4c6c,0x651b,
0x3fb2,0xb948,0xa3fe,0xc4b6,
0x3ff2,0x08fe,0xc215,0x96d6,
0x4014,0x72c4,0xf8b1,0x3a6a,
0x4020,0xd91c,0x8b5d,0x2f16,
0x4014,0xdbaa,0x142f,0x81a2,
0x3ff0,0x0000,0x0000,0x0000,
};
static unsigned short PQ[28] = {
0x3f42,0xb89b,0x1344,0x3d69,
0x3fb1,0x9fdd,0x5948,0xaa83,
0x3ff1,0xaea9,0xb850,0xeed6,
0x4014,0x4ba2,0xf7d2,0x51a1,
0x4020,0xccb9,0xdda2,0xfd65,
0x4014,0xd6dd,0x4762,0xb4d9,
0x3ff0,0x0000,0x0000,0x0000,
};
#endif

#ifdef UNK
static double QP[8] = {
 5.10862594750176621635E-2,
 4.98213872951233449420E0,
 7.58238284132545283818E1,
 3.66779609360150777800E2,
 7.10856304998926107277E2,
 5.97489612400613639965E2,
 2.11688757100572135698E2,
 2.52070205858023719784E1,
};
static double QQ[7] = {
/* 1.00000000000000000000E0,*/
 7.42373277035675149943E1,
 1.05644886038262816351E3,
 4.98641058337653607651E3,
 9.56231892404756170795E3,
 7.99704160447350683650E3,
 2.82619278517639096600E3,
 3.36093607810698293419E2,
};
#endif
#ifdef DEC
static unsigned short QP[32] = {
0037121,0037723,0055605,0151004,
0040637,0066656,0031554,0077264,
0041627,0122714,0153170,0161466,
0042267,0061712,0036520,0140145,
0042461,0133315,0131573,0071176,
0042425,0057525,0147500,0013201,
0042123,0130122,0061245,0154131,
0041311,0123772,0064254,0172650,
};
static unsigned short QQ[28] = {
/*0040200,0000000,0000000,0000000,*/
0041624,0074603,0002112,0101670,
0042604,0007135,0010162,0175565,
0043233,0151510,0157757,0172010,
0043425,0064506,0112006,0104276,
0043371,0164125,0032271,0164242,
0043060,0121425,0122750,0136013,
0042250,0005773,0053472,0146267,
};
#endif
#ifdef IBMPC
static unsigned short QP[32] = {
0xba40,0x6b70,0x27fa,0x3faa,
0x8fd6,0xc66d,0xedb5,0x4013,
0x1c67,0x9acf,0xf4b9,0x4052,
0x180d,0x47aa,0xec79,0x4076,
0x6e50,0xb66f,0x36d9,0x4086,
0x02d0,0xb9e8,0xabea,0x4082,
0xbb0b,0x4c54,0x760a,0x406a,
0x9eb5,0x4d15,0x34ff,0x4039,
};
static unsigned short QQ[28] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x5077,0x6089,0x8f30,0x4052,
0x5f6f,0xa20e,0x81cb,0x4090,
0xfe81,0x1bfd,0x7a69,0x40b3,
0xd118,0xd280,0xad28,0x40c2,
0x3d14,0xa697,0x3d0a,0x40bf,
0x1781,0xb4bd,0x1462,0x40a6,
0x5997,0x6ae7,0x017f,0x4075,
};
#endif
#ifdef MIEEE
static unsigned short QP[32] = {
0x3faa,0x27fa,0x6b70,0xba40,
0x4013,0xedb5,0xc66d,0x8fd6,
0x4052,0xf4b9,0x9acf,0x1c67,
0x4076,0xec79,0x47aa,0x180d,
0x4086,0x36d9,0xb66f,0x6e50,
0x4082,0xabea,0xb9e8,0x02d0,
0x406a,0x760a,0x4c54,0xbb0b,
0x4039,0x34ff,0x4d15,0x9eb5,
};
static unsigned short QQ[28] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4052,0x8f30,0x6089,0x5077,
0x4090,0x81cb,0xa20e,0x5f6f,
0x40b3,0x7a69,0x1bfd,0xfe81,
0x40c2,0xad28,0xd280,0xd118,
0x40bf,0x3d0a,0xa697,0x3d14,
0x40a6,0x1462,0xb4bd,0x1781,
0x4075,0x017f,0x6ae7,0x5997,
};
#endif

#ifdef UNK
static double YP[6] = {
 1.26320474790178026440E9,
-6.47355876379160291031E11,
 1.14509511541823727583E14,
-8.12770255501325109621E15,
 2.02439475713594898196E17,
-7.78877196265950026825E17,
};
static double YQ[8] = {
/* 1.00000000000000000000E0,*/
 5.94301592346128195359E2,
 2.35564092943068577943E5,
 7.34811944459721705660E7,
 1.87601316108706159478E10,
 3.88231277496238566008E12,
 6.20557727146953693363E14,
 6.87141087355300489866E16,
 3.97270608116560655612E18,
};
#endif
#ifdef DEC
static unsigned short YP[24] = {
0047626,0112763,0013715,0133045,
0152026,0134552,0142033,0024411,
0053720,0045245,0102210,0077565,
0155347,0000321,0136415,0102031,
0056463,0146550,0055633,0032605,
0157054,0171012,0167361,0054265,
};
static unsigned short YQ[32] = {
/*0040200,0000000,0000000,0000000,*/
0042424,0111515,0044773,0153014,
0044546,0005405,0171307,0075774,
0046614,0023575,0047105,0063556,
0050613,0143034,0101533,0156026,
0052541,0175367,0166514,0114257,
0054415,0014466,0134350,0171154,
0056164,0017436,0025075,0022101,
0057534,0103614,0103663,0121772,
};
#endif
#ifdef IBMPC
static unsigned short YP[24] = {
0xb6c5,0x62f9,0xd2be,0x41d2,
0x6521,0x5883,0xd72d,0xc262,
0x0fef,0xb091,0x0954,0x42da,
0xb083,0x37a1,0xe01a,0xc33c,
0x66b1,0x0b73,0x79ad,0x4386,
0x2b17,0x5dde,0x9e41,0xc3a5,
};
static unsigned short YQ[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x7ac2,0xa93f,0x9269,0x4082,
0xef7f,0xbe58,0xc160,0x410c,
0xacee,0xa9c8,0x84ef,0x4191,
0x7b83,0x906b,0x78c3,0x4211,
0x9316,0xfda9,0x3f5e,0x428c,
0x1e4e,0xd71d,0xa326,0x4301,
0xa488,0xc547,0x83e3,0x436e,
0x747f,0x90f6,0x90f1,0x43cb,
};
#endif
#ifdef MIEEE
static unsigned short YP[24] = {
0x41d2,0xd2be,0x62f9,0xb6c5,
0xc262,0xd72d,0x5883,0x6521,
0x42da,0x0954,0xb091,0x0fef,
0xc33c,0xe01a,0x37a1,0xb083,
0x4386,0x79ad,0x0b73,0x66b1,
0xc3a5,0x9e41,0x5dde,0x2b17,
};
static unsigned short YQ[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4082,0x9269,0xa93f,0x7ac2,
0x410c,0xc160,0xbe58,0xef7f,
0x4191,0x84ef,0xa9c8,0xacee,
0x4211,0x78c3,0x906b,0x7b83,
0x428c,0x3f5e,0xfda9,0x9316,
0x4301,0xa326,0xd71d,0x1e4e,
0x436e,0x83e3,0xc547,0xa488,
0x43cb,0x90f1,0x90f6,0x747f,
};
#endif


#ifdef UNK
static double Z1 = 1.46819706421238932572E1;
static double Z2 = 4.92184563216946036703E1;
#endif

#ifdef DEC
static unsigned short DZ1[] = {0041152,0164532,0006114,0010540};
static unsigned short DZ2[] = {0041504,0157663,0001625,0020621};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

#ifdef IBMPC
static unsigned short DZ1[] = {0x822c,0x4189,0x5d2b,0x402d};
static unsigned short DZ2[] = {0xa432,0x6072,0x9bf6,0x4048};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

#ifdef MIEEE
static unsigned short DZ1[] = {0x402d,0x5d2b,0x4189,0x822c};
static unsigned short DZ2[] = {0x4048,0x9bf6,0x6072,0xa432};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)
#endif

extern double TWOOPI, THPIO4, SQ2OPI;

double j1(x)
double x;
{
double w, z, p, q, xn;

w = x;
if( x < 0 )
        return -j1(-x);

if( w <= 5.0 )
	{
	z = x * x;	
	w = polevl( z, RP, 3 ) / p1evl( z, RQ, 8 );
	w = w * x * (z - Z1) * (z - Z2);
	return( w );
	}

w = 5.0/x;
z = w * w;
p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
xn = x - THPIO4;
p = p * cos(xn) - w * q * sin(xn);
return( p * SQ2OPI / sqrt(x) );
}


double y1(x)
double x;
{
double w, z, p, q, xn;

if( x <= 5.0 )
	{
	if (x == 0.0) {
		mtherr("y1", SING);
		return -NPY_INFINITY;
	} else if (x <= 0.0) {
		mtherr("y1", DOMAIN);
		return NPY_NAN;
	}
	z = x * x;
	w = x * (polevl( z, YP, 5 ) / p1evl( z, YQ, 8 ));
	w += TWOOPI * ( j1(x) * log(x)  -  1.0/x );
	return( w );
	}

w = 5.0/x;
z = w * w;
p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
xn = x - THPIO4;
p = p * sin(xn) + w * q * cos(xn);
return( p * SQ2OPI / sqrt(x) );
}
