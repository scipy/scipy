/*                                                     i0.c
 *
 *     Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0();
 *
 * y = i0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
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
 *    DEC       0,30         6000       8.2e-17     1.9e-17
 *    IEEE      0,30        30000       5.8e-16     1.4e-16
 *
 */
/*							i0e.c
 *
 *	Modified Bessel function of order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0e();
 *
 * y = i0e( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order zero of the argument.
 *
 * The function is defined as i0e(x) = exp(-|x|) j0( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       5.4e-16     1.2e-16
 * See i0().
 *
 */

/*                                                     i0.c            */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */

#include "mconf.h"

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */

#ifdef UNK
static double A[] = {
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1
};
#endif

#ifdef DEC
static unsigned short A[] = {
    0121642, 0162671, 0004646, 0103567,
    0022431, 0115424, 0135755, 0026104,
    0123214, 0023533, 0110365, 0156635,
    0023767, 0033304, 0117662, 0172716,
    0124522, 0100426, 0012277, 0157531,
    0025254, 0155062, 0054461, 0030465,
    0126010, 0131143, 0013560, 0153604,
    0026517, 0170577, 0006336, 0114437,
    0127227, 0162253, 0152243, 0052734,
    0027724, 0142766, 0061641, 0160200,
    0130416, 0123760, 0116564, 0125262,
    0031066, 0144035, 0021246, 0054641,
    0131537, 0053664, 0060131, 0102530,
    0032201, 0155664, 0165153, 0020652,
    0132617, 0061434, 0074423, 0176145,
    0033225, 0174444, 0136147, 0122542,
    0133624, 0031576, 0056453, 0020470,
    0034211, 0175305, 0172321, 0041314,
    0134561, 0054462, 0147040, 0165315,
    0035105, 0124333, 0120203, 0162532,
    0135427, 0013750, 0174257, 0055221,
    0035726, 0161654, 0050220, 0100162,
    0136215, 0131361, 0000325, 0041110,
    0036454, 0145417, 0117357, 0017352,
    0136702, 0072367, 0104415, 0133574,
    0037111, 0172126, 0072505, 0014544,
    0137302, 0055601, 0120550, 0033523,
    0037457, 0136543, 0136544, 0043002,
    0137633, 0177536, 0001276, 0066150,
    0040055, 0041164, 0100655, 0010521
};
#endif

#ifdef IBMPC
static unsigned short A[] = {
    0xd0ef, 0x2134, 0x5cb7, 0xbc54,
    0xa589, 0x977d, 0x3362, 0x3c83,
    0xbbb4, 0x721e, 0x84eb, 0xbcb1,
    0x5eba, 0x93f6, 0xe6d8, 0x3cde,
    0xfbeb, 0xc297, 0x5022, 0xbd0a,
    0x2627, 0x4b26, 0x9b46, 0x3d35,
    0x1af0, 0x62ee, 0x164c, 0xbd61,
    0xd324, 0xe19b, 0xfe2f, 0x3d89,
    0x6abc, 0x7a94, 0xfc95, 0xbdb2,
    0x3c10, 0xcc74, 0x98be, 0x3dda,
    0x9556, 0x13ae, 0xd4fe, 0xbe01,
    0xcb34, 0xa454, 0xd903, 0x3e26,
    0x30ab, 0x8c0b, 0xeaf6, 0xbe4b,
    0x6435, 0x9d4d, 0x3b76, 0x3e70,
    0x7f8d, 0x8f22, 0xec63, 0xbe91,
    0xf4ac, 0x978c, 0xbf24, 0x3eb2,
    0x6427, 0xcba5, 0x866f, 0xbed2,
    0x2859, 0xbe9a, 0x3f58, 0x3ef1,
    0x1d5a, 0x59c4, 0x2b26, 0xbf0e,
    0x7cab, 0x7410, 0xb51b, 0x3f28,
    0xeb52, 0x1f15, 0xe2fd, 0xbf42,
    0x100e, 0x8a12, 0xdc75, 0x3f5a,
    0xa849, 0x201a, 0xb65e, 0xbf71,
    0xe3dd, 0xf3dd, 0x9961, 0x3f85,
    0xb6f0, 0xf121, 0x4e9e, 0xbf98,
    0xa32d, 0xcea8, 0x3e8a, 0x3fa9,
    0x06ea, 0x342d, 0x4b70, 0xbfb8,
    0x88c0, 0x77ac, 0xf7ac, 0x3fc5,
    0xcd8d, 0xc057, 0x7feb, 0xbfd3,
    0xa22a, 0x9035, 0xa84e, 0x3fe5,
};
#endif

#ifdef MIEEE
static unsigned short A[] = {
    0xbc54, 0x5cb7, 0x2134, 0xd0ef,
    0x3c83, 0x3362, 0x977d, 0xa589,
    0xbcb1, 0x84eb, 0x721e, 0xbbb4,
    0x3cde, 0xe6d8, 0x93f6, 0x5eba,
    0xbd0a, 0x5022, 0xc297, 0xfbeb,
    0x3d35, 0x9b46, 0x4b26, 0x2627,
    0xbd61, 0x164c, 0x62ee, 0x1af0,
    0x3d89, 0xfe2f, 0xe19b, 0xd324,
    0xbdb2, 0xfc95, 0x7a94, 0x6abc,
    0x3dda, 0x98be, 0xcc74, 0x3c10,
    0xbe01, 0xd4fe, 0x13ae, 0x9556,
    0x3e26, 0xd903, 0xa454, 0xcb34,
    0xbe4b, 0xeaf6, 0x8c0b, 0x30ab,
    0x3e70, 0x3b76, 0x9d4d, 0x6435,
    0xbe91, 0xec63, 0x8f22, 0x7f8d,
    0x3eb2, 0xbf24, 0x978c, 0xf4ac,
    0xbed2, 0x866f, 0xcba5, 0x6427,
    0x3ef1, 0x3f58, 0xbe9a, 0x2859,
    0xbf0e, 0x2b26, 0x59c4, 0x1d5a,
    0x3f28, 0xb51b, 0x7410, 0x7cab,
    0xbf42, 0xe2fd, 0x1f15, 0xeb52,
    0x3f5a, 0xdc75, 0x8a12, 0x100e,
    0xbf71, 0xb65e, 0x201a, 0xa849,
    0x3f85, 0x9961, 0xf3dd, 0xe3dd,
    0xbf98, 0x4e9e, 0xf121, 0xb6f0,
    0x3fa9, 0x3e8a, 0xcea8, 0xa32d,
    0xbfb8, 0x4b70, 0x342d, 0x06ea,
    0x3fc5, 0xf7ac, 0x77ac, 0x88c0,
    0xbfd3, 0x7feb, 0xc057, 0xcd8d,
    0x3fe5, 0xa84e, 0x9035, 0xa22a
};
#endif


/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

#ifdef UNK
static double B[] = {
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1
};
#endif

#ifdef DEC
static unsigned short B[] = {
    0122005, 0066672, 0123124, 0054311,
    0121662, 0033323, 0030214, 0104602,
    0022515, 0170300, 0113314, 0020413,
    0022437, 0117350, 0035402, 0007146,
    0123243, 0000135, 0057220, 0177435,
    0123305, 0073476, 0144106, 0170702,
    0023777, 0071755, 0017527, 0154373,
    0024211, 0052214, 0102247, 0033270,
    0124454, 0017763, 0171453, 0012322,
    0125072, 0166316, 0075505, 0154616,
    0024612, 0133770, 0065376, 0025045,
    0025730, 0162143, 0056036, 0001632,
    0026112, 0015077, 0150464, 0063542,
    0126374, 0101030, 0014274, 0065457,
    0127150, 0077271, 0125763, 0157617,
    0127412, 0104350, 0040713, 0120445,
    0027121, 0023765, 0057500, 0001165,
    0030407, 0147146, 0003643, 0075644,
    0031151, 0061445, 0044422, 0156065,
    0031702, 0132224, 0003266, 0125551,
    0032534, 0000076, 0147153, 0005555,
    0033502, 0004536, 0004016, 0026055,
    0034620, 0076433, 0142314, 0171215,
    0036134, 0146145, 0013454, 0101104,
    0040115, 0171425, 0062500, 0047133
};
#endif

#ifdef IBMPC
static unsigned short B[] = {
    0x8b19, 0x54ca, 0xadb7, 0xbc60,
    0x9130, 0x6611, 0x46da, 0xbc56,
    0x8421, 0x12d9, 0xbe18, 0x3c89,
    0x41cd, 0x0760, 0xf3dd, 0x3c83,
    0x1fe4, 0xabd2, 0x600b, 0xbcb4,
    0xde38, 0xd908, 0xaee7, 0xbcb8,
    0xfb1f, 0xa3ea, 0xee7d, 0x3cdf,
    0xe6d7, 0x9094, 0x2a91, 0x3cf1,
    0x629a, 0x7e65, 0x83fe, 0xbd05,
    0xbb32, 0xcf68, 0x5d99, 0xbd27,
    0xc545, 0x0d5f, 0x56ff, 0x3d11,
    0xc073, 0x6b83, 0x1c8c, 0x3d5b,
    0x8cec, 0xfa26, 0x4347, 0x3d69,
    0x8d66, 0x0317, 0x9043, 0xbd7f,
    0x7bf2, 0x357e, 0x0fd7, 0xbdad,
    0x7425, 0x0839, 0x511d, 0xbdc1,
    0x004f, 0xabe8, 0x24fe, 0x3daa,
    0x6f75, 0xc0f4, 0xf9cc, 0x3e00,
    0x5b87, 0xa922, 0x2c64, 0x3e2d,
    0xd56d, 0x80d6, 0x5692, 0x3e58,
    0x616e, 0xd9cd, 0x8007, 0x3e8b,
    0xc586, 0xc101, 0x412b, 0x3ec8,
    0x9e52, 0x7899, 0x0fa3, 0x3f12,
    0x9049, 0xa2e5, 0x998c, 0x3f6b,
    0x09cb, 0xaca8, 0xbe62, 0x3fe9
};
#endif

#ifdef MIEEE
static unsigned short B[] = {
    0xbc60, 0xadb7, 0x54ca, 0x8b19,
    0xbc56, 0x46da, 0x6611, 0x9130,
    0x3c89, 0xbe18, 0x12d9, 0x8421,
    0x3c83, 0xf3dd, 0x0760, 0x41cd,
    0xbcb4, 0x600b, 0xabd2, 0x1fe4,
    0xbcb8, 0xaee7, 0xd908, 0xde38,
    0x3cdf, 0xee7d, 0xa3ea, 0xfb1f,
    0x3cf1, 0x2a91, 0x9094, 0xe6d7,
    0xbd05, 0x83fe, 0x7e65, 0x629a,
    0xbd27, 0x5d99, 0xcf68, 0xbb32,
    0x3d11, 0x56ff, 0x0d5f, 0xc545,
    0x3d5b, 0x1c8c, 0x6b83, 0xc073,
    0x3d69, 0x4347, 0xfa26, 0x8cec,
    0xbd7f, 0x9043, 0x0317, 0x8d66,
    0xbdad, 0x0fd7, 0x357e, 0x7bf2,
    0xbdc1, 0x511d, 0x0839, 0x7425,
    0x3daa, 0x24fe, 0xabe8, 0x004f,
    0x3e00, 0xf9cc, 0xc0f4, 0x6f75,
    0x3e2d, 0x2c64, 0xa922, 0x5b87,
    0x3e58, 0x5692, 0x80d6, 0xd56d,
    0x3e8b, 0x8007, 0xd9cd, 0x616e,
    0x3ec8, 0x412b, 0xc101, 0xc586,
    0x3f12, 0x0fa3, 0x7899, 0x9e52,
    0x3f6b, 0x998c, 0xa2e5, 0x9049,
    0x3fe9, 0xbe62, 0xaca8, 0x09cb
};
#endif

double i0(x)
double x;
{
    double y;

    if (x < 0)
	x = -x;
    if (x <= 8.0) {
	y = (x / 2.0) - 2.0;
	return (exp(x) * chbevl(y, A, 30));
    }

    return (exp(x) * chbevl(32.0 / x - 2.0, B, 25) / sqrt(x));

}




double i0e(x)
double x;
{
    double y;

    if (x < 0)
	x = -x;
    if (x <= 8.0) {
	y = (x / 2.0) - 2.0;
	return (chbevl(y, A, 30));
    }

    return (chbevl(32.0 / x - 2.0, B, 25) / sqrt(x));

}
