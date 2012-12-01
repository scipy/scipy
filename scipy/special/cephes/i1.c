/*                                                     i1.c
 *
 *     Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1();
 *
 * y = i1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        3400       1.2e-16     2.3e-17
 *    IEEE      0, 30       30000       1.9e-15     2.1e-16
 *
 *
 */
/*							i1e.c
 *
 *	Modified Bessel function of order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i1e();
 *
 * y = i1e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order one of the argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       2.0e-15     2.0e-16
 * See i1().
 *
 */

/*                                                     i1.c 2          */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1985, 1987, 2000 by Stephen L. Moshier
 */

#include "mconf.h"

/* Chebyshev coefficients for exp(-x) I1(x) / x
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
 */

#ifdef UNK
static double A[] = {
    2.77791411276104639959E-18,
    -2.11142121435816608115E-17,
    1.55363195773620046921E-16,
    -1.10559694773538630805E-15,
    7.60068429473540693410E-15,
    -5.04218550472791168711E-14,
    3.22379336594557470981E-13,
    -1.98397439776494371520E-12,
    1.17361862988909016308E-11,
    -6.66348972350202774223E-11,
    3.62559028155211703701E-10,
    -1.88724975172282928790E-9,
    9.38153738649577178388E-9,
    -4.44505912879632808065E-8,
    2.00329475355213526229E-7,
    -8.56872026469545474066E-7,
    3.47025130813767847674E-6,
    -1.32731636560394358279E-5,
    4.78156510755005422638E-5,
    -1.61760815825896745588E-4,
    5.12285956168575772895E-4,
    -1.51357245063125314899E-3,
    4.15642294431288815669E-3,
    -1.05640848946261981558E-2,
    2.47264490306265168283E-2,
    -5.29459812080949914269E-2,
    1.02643658689847095384E-1,
    -1.76416518357834055153E-1,
    2.52587186443633654823E-1
};
#endif

#ifdef DEC
static unsigned short A[] = {
    0021514, 0174520, 0060742, 0000241,
    0122302, 0137206, 0016120, 0025663,
    0023063, 0017437, 0026235, 0176536,
    0123637, 0052523, 0170150, 0125632,
    0024410, 0165770, 0030251, 0044134,
    0125143, 0012160, 0162170, 0054727,
    0025665, 0075702, 0035716, 0145247,
    0126413, 0116032, 0176670, 0015462,
    0027116, 0073425, 0110351, 0105242,
    0127622, 0104034, 0137530, 0037364,
    0030307, 0050645, 0120776, 0175535,
    0131001, 0130331, 0043523, 0037455,
    0031441, 0026160, 0010712, 0100174,
    0132076, 0164761, 0022706, 0017500,
    0032527, 0015045, 0115076, 0104076,
    0133146, 0001714, 0015434, 0144520,
    0033550, 0161166, 0124215, 0077050,
    0134136, 0127715, 0143365, 0157170,
    0034510, 0106652, 0013070, 0064130,
    0135051, 0117126, 0117264, 0123761,
    0035406, 0045355, 0133066, 0175751,
    0135706, 0061420, 0054746, 0122440,
    0036210, 0031232, 0047235, 0006640,
    0136455, 0012373, 0144235, 0011523,
    0036712, 0107437, 0036731, 0015111,
    0137130, 0156742, 0115744, 0172743,
    0037322, 0033326, 0124667, 0124740,
    0137464, 0123210, 0021510, 0144556,
    0037601, 0051433, 0111123, 0177721
};
#endif

#ifdef IBMPC
static unsigned short A[] = {
    0x4014, 0x0c3c, 0x9f2a, 0x3c49,
    0x0576, 0xc38a, 0x57d0, 0xbc78,
    0xbfac, 0xe593, 0x63e3, 0x3ca6,
    0x1573, 0x7e0d, 0xeaaa, 0xbcd3,
    0x290c, 0x0615, 0x1d7f, 0x3d01,
    0x0b3b, 0x1c8f, 0x628e, 0xbd2c,
    0xd955, 0x4779, 0xaf78, 0x3d56,
    0x0366, 0x5fb7, 0x7383, 0xbd81,
    0x3154, 0xb21d, 0xcee2, 0x3da9,
    0x07de, 0x97eb, 0x5103, 0xbdd2,
    0xdf6c, 0xb43f, 0xea34, 0x3df8,
    0x67e6, 0x28ea, 0x361b, 0xbe20,
    0x5010, 0x0239, 0x258e, 0x3e44,
    0xc3e8, 0x24b8, 0xdd3e, 0xbe67,
    0xd108, 0xb347, 0xe344, 0x3e8a,
    0x992a, 0x8363, 0xc079, 0xbeac,
    0xafc5, 0xd511, 0x1c4e, 0x3ecd,
    0xbbcf, 0xb8de, 0xd5f9, 0xbeeb,
    0x0d0b, 0x42c7, 0x11b5, 0x3f09,
    0x94fe, 0xd3d6, 0x33ca, 0xbf25,
    0xdf7d, 0xb6c6, 0xc95d, 0x3f40,
    0xd4a4, 0x0b3c, 0xcc62, 0xbf58,
    0xa1b4, 0x49d3, 0x0653, 0x3f71,
    0xa26a, 0x7913, 0xa29f, 0xbf85,
    0x2349, 0xe7bb, 0x51e3, 0x3f99,
    0x9ebc, 0x537c, 0x1bbc, 0xbfab,
    0xf53c, 0xd536, 0x46da, 0x3fba,
    0x192e, 0x0469, 0x94d1, 0xbfc6,
    0x7ffa, 0x724a, 0x2a63, 0x3fd0
};
#endif

#ifdef MIEEE
static unsigned short A[] = {
    0x3c49, 0x9f2a, 0x0c3c, 0x4014,
    0xbc78, 0x57d0, 0xc38a, 0x0576,
    0x3ca6, 0x63e3, 0xe593, 0xbfac,
    0xbcd3, 0xeaaa, 0x7e0d, 0x1573,
    0x3d01, 0x1d7f, 0x0615, 0x290c,
    0xbd2c, 0x628e, 0x1c8f, 0x0b3b,
    0x3d56, 0xaf78, 0x4779, 0xd955,
    0xbd81, 0x7383, 0x5fb7, 0x0366,
    0x3da9, 0xcee2, 0xb21d, 0x3154,
    0xbdd2, 0x5103, 0x97eb, 0x07de,
    0x3df8, 0xea34, 0xb43f, 0xdf6c,
    0xbe20, 0x361b, 0x28ea, 0x67e6,
    0x3e44, 0x258e, 0x0239, 0x5010,
    0xbe67, 0xdd3e, 0x24b8, 0xc3e8,
    0x3e8a, 0xe344, 0xb347, 0xd108,
    0xbeac, 0xc079, 0x8363, 0x992a,
    0x3ecd, 0x1c4e, 0xd511, 0xafc5,
    0xbeeb, 0xd5f9, 0xb8de, 0xbbcf,
    0x3f09, 0x11b5, 0x42c7, 0x0d0b,
    0xbf25, 0x33ca, 0xd3d6, 0x94fe,
    0x3f40, 0xc95d, 0xb6c6, 0xdf7d,
    0xbf58, 0xcc62, 0x0b3c, 0xd4a4,
    0x3f71, 0x0653, 0x49d3, 0xa1b4,
    0xbf85, 0xa29f, 0x7913, 0xa26a,
    0x3f99, 0x51e3, 0xe7bb, 0x2349,
    0xbfab, 0x1bbc, 0x537c, 0x9ebc,
    0x3fba, 0x46da, 0xd536, 0xf53c,
    0xbfc6, 0x94d1, 0x0469, 0x192e,
    0x3fd0, 0x2a63, 0x724a, 0x7ffa
};
#endif

/*                                                     i1.c    */

/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
 */

#ifdef UNK
static double B[] = {
    7.51729631084210481353E-18,
    4.41434832307170791151E-18,
    -4.65030536848935832153E-17,
    -3.20952592199342395980E-17,
    2.96262899764595013876E-16,
    3.30820231092092828324E-16,
    -1.88035477551078244854E-15,
    -3.81440307243700780478E-15,
    1.04202769841288027642E-14,
    4.27244001671195135429E-14,
    -2.10154184277266431302E-14,
    -4.08355111109219731823E-13,
    -7.19855177624590851209E-13,
    2.03562854414708950722E-12,
    1.41258074366137813316E-11,
    3.25260358301548823856E-11,
    -1.89749581235054123450E-11,
    -5.58974346219658380687E-10,
    -3.83538038596423702205E-9,
    -2.63146884688951950684E-8,
    -2.51223623787020892529E-7,
    -3.88256480887769039346E-6,
    -1.10588938762623716291E-4,
    -9.76109749136146840777E-3,
    7.78576235018280120474E-1
};
#endif

#ifdef DEC
static unsigned short B[] = {
    0022012, 0125555, 0115227, 0043456,
    0021642, 0156127, 0052075, 0145203,
    0122526, 0072435, 0111231, 0011664,
    0122424, 0001544, 0161671, 0114403,
    0023252, 0144257, 0163532, 0142121,
    0023276, 0132162, 0174045, 0013204,
    0124007, 0077154, 0057046, 0110517,
    0124211, 0066650, 0116127, 0157073,
    0024473, 0133413, 0130551, 0107504,
    0025100, 0064741, 0032631, 0040364,
    0124675, 0045101, 0071551, 0012400,
    0125745, 0161054, 0071637, 0011247,
    0126112, 0117410, 0035525, 0122231,
    0026417, 0037237, 0131034, 0176427,
    0027170, 0100373, 0024742, 0025725,
    0027417, 0006417, 0105303, 0141446,
    0127246, 0163716, 0121202, 0060137,
    0130431, 0123122, 0120436, 0166000,
    0131203, 0144134, 0153251, 0124500,
    0131742, 0005234, 0122732, 0033006,
    0132606, 0157751, 0072362, 0121031,
    0133602, 0043372, 0047120, 0015626,
    0134747, 0165774, 0001125, 0046462,
    0136437, 0166402, 0117746, 0155137,
    0040107, 0050305, 0125330, 0124241
};
#endif

#ifdef IBMPC
static unsigned short B[] = {
    0xe8e6, 0xb352, 0x556d, 0x3c61,
    0xb950, 0xea87, 0x5b8a, 0x3c54,
    0x2277, 0xb253, 0xcea3, 0xbc8a,
    0x3320, 0x9c77, 0x806c, 0xbc82,
    0x588a, 0xfceb, 0x5915, 0x3cb5,
    0xa2d1, 0x5f04, 0xd68e, 0x3cb7,
    0xd22a, 0x8bc4, 0xefcd, 0xbce0,
    0xfbc7, 0x138a, 0x2db5, 0xbcf1,
    0x31e8, 0x762d, 0x76e1, 0x3d07,
    0x281e, 0x26b3, 0x0d3c, 0x3d28,
    0x22a0, 0x2e6d, 0xa948, 0xbd17,
    0xe255, 0x8e73, 0xbc45, 0xbd5c,
    0xb493, 0x076a, 0x53e1, 0xbd69,
    0x9fa3, 0xf643, 0xe7d3, 0x3d81,
    0x457b, 0x653c, 0x101f, 0x3daf,
    0x7865, 0xf158, 0xe1a1, 0x3dc1,
    0x4c0c, 0xd450, 0xdcf9, 0xbdb4,
    0xdd80, 0x5423, 0x34ca, 0xbe03,
    0x3528, 0x9ad5, 0x790b, 0xbe30,
    0x46c1, 0x94bb, 0x4153, 0xbe5c,
    0x5443, 0x2e9e, 0xdbfd, 0xbe90,
    0x0373, 0x49ca, 0x48df, 0xbed0,
    0xa9a6, 0x804a, 0xfd7f, 0xbf1c,
    0xdb4c, 0x53fc, 0xfda0, 0xbf83,
    0x1514, 0xb55b, 0xea18, 0x3fe8
};
#endif

#ifdef MIEEE
static unsigned short B[] = {
    0x3c61, 0x556d, 0xb352, 0xe8e6,
    0x3c54, 0x5b8a, 0xea87, 0xb950,
    0xbc8a, 0xcea3, 0xb253, 0x2277,
    0xbc82, 0x806c, 0x9c77, 0x3320,
    0x3cb5, 0x5915, 0xfceb, 0x588a,
    0x3cb7, 0xd68e, 0x5f04, 0xa2d1,
    0xbce0, 0xefcd, 0x8bc4, 0xd22a,
    0xbcf1, 0x2db5, 0x138a, 0xfbc7,
    0x3d07, 0x76e1, 0x762d, 0x31e8,
    0x3d28, 0x0d3c, 0x26b3, 0x281e,
    0xbd17, 0xa948, 0x2e6d, 0x22a0,
    0xbd5c, 0xbc45, 0x8e73, 0xe255,
    0xbd69, 0x53e1, 0x076a, 0xb493,
    0x3d81, 0xe7d3, 0xf643, 0x9fa3,
    0x3daf, 0x101f, 0x653c, 0x457b,
    0x3dc1, 0xe1a1, 0xf158, 0x7865,
    0xbdb4, 0xdcf9, 0xd450, 0x4c0c,
    0xbe03, 0x34ca, 0x5423, 0xdd80,
    0xbe30, 0x790b, 0x9ad5, 0x3528,
    0xbe5c, 0x4153, 0x94bb, 0x46c1,
    0xbe90, 0xdbfd, 0x2e9e, 0x5443,
    0xbed0, 0x48df, 0x49ca, 0x0373,
    0xbf1c, 0xfd7f, 0x804a, 0xa9a6,
    0xbf83, 0xfda0, 0x53fc, 0xdb4c,
    0x3fe8, 0xea18, 0xb55b, 0x1514
};
#endif

/*                                                     i1.c    */

double i1(x)
double x;
{
    double y, z;

    z = fabs(x);
    if (z <= 8.0) {
	y = (z / 2.0) - 2.0;
	z = chbevl(y, A, 29) * z * exp(z);
    }
    else {
	z = exp(z) * chbevl(32.0 / z - 2.0, B, 25) / sqrt(z);
    }
    if (x < 0.0)
	z = -z;
    return (z);
}

/*                                                     i1e()   */

double i1e(x)
double x;
{
    double y, z;

    z = fabs(x);
    if (z <= 8.0) {
	y = (z / 2.0) - 2.0;
	z = chbevl(y, A, 29) * z;
    }
    else {
	z = chbevl(32.0 / z - 2.0, B, 25) / sqrt(z);
    }
    if (x < 0.0)
	z = -z;
    return (z);
}
