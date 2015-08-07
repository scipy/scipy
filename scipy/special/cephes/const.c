/*                                                     const.c
 *
 *     Globally declared constants
 *
 *
 *
 * SYNOPSIS:
 *
 * extern double nameofconstant;
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains a number of mathematical constants and
 * also some needed size parameters of the computer arithmetic.
 * The values are supplied as arrays of hexadecimal integers
 * for IEEE arithmetic, and in a normal decimal scientific notation for
 * other machines.  The particular notation used is determined
 * by a symbol (IBMPC, or UNK) defined in the include file
 * mconf.h.
 *
 * The default size parameters are as follows.
 *
 * For UNK mode:
 * MACHEP =  1.38777878078144567553E-17       2**-56
 * MAXLOG =  8.8029691931113054295988E1       log(2**127)
 * MINLOG = -8.872283911167299960540E1        log(2**-128)
 *
 * For IEEE arithmetic (IBMPC):
 * MACHEP =  1.11022302462515654042E-16       2**-53
 * MAXLOG =  7.09782712893383996843E2         log(2**1024)
 * MINLOG = -7.08396418532264106224E2         log(2**-1022)
 *
 * The global symbols for mathematical constants are
 * SQ2OPI =  7.9788456080286535587989E-1      sqrt( 2/pi )
 * LOGSQ2 =  3.46573590279972654709E-1        log(2)/2
 * THPIO4 =  2.35619449019234492885           3*pi/4
 *
 * These lists are subject to change.
 */

/*                                                     const.c */

/*
 * Cephes Math Library Release 2.3:  March, 1995
 * Copyright 1984, 1995 by Stephen L. Moshier
 */

#include "mconf.h"

#ifdef UNK
double MACHEP = 1.11022302462515654042E-16;	/* 2**-53 */

#ifdef DENORMAL
double MAXLOG = 7.09782712893383996732E2;	/* log(DBL_MAX) */

												      /* double MINLOG = -7.44440071921381262314E2; *//* log(2**-1074) */
double MINLOG = -7.451332191019412076235E2;	/* log(2**-1075) */
#else
double MAXLOG = 7.08396418532264106224E2;	/* log 2**1022 */
double MINLOG = -7.08396418532264106224E2;	/* log 2**-1022 */
#endif
double SQ2OPI = 7.9788456080286535587989E-1;	/* sqrt( 2/pi ) */
double LOGSQ2 = 3.46573590279972654709E-1;	/* log(2)/2 */
double THPIO4 = 2.35619449019234492885;	/* 3*pi/4 */

#endif

#ifdef IBMPC
			/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MACHEP[4] = { 0x0000, 0x0000, 0x0000, 0x3ca0 };

#ifdef DENORMAL
			/* log(DBL_MAX) =  7.09782712893383996732224E2 */
unsigned short MAXLOG[4] = { 0x39ef, 0xfefa, 0x2e42, 0x4086 };

			/* log(2**-1074) = - -7.44440071921381262314E2 */
/*unsigned short MINLOG[4] = {0x71c3,0x446d,0x4385,0xc087}; */
unsigned short MINLOG[4] = { 0x3052, 0xd52d, 0x4910, 0xc087 };
#else
			/* log(2**1022) =   7.08396418532264106224E2 */
unsigned short MAXLOG[4] = { 0xbcd2, 0xdd7a, 0x232b, 0x4086 };

			/* log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MINLOG[4] = { 0xbcd2, 0xdd7a, 0x232b, 0xc086 };
#endif
			/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short SQ2OPI[4] = { 0x3651, 0x33d4, 0x8845, 0x3fe9 };
unsigned short LOGSQ2[4] = { 0x39ef, 0xfefa, 0x2e42, 0x3fd6 };
unsigned short THPIO4[4] = { 0x21d2, 0x7f33, 0xd97c, 0x4002 };

#endif

#ifdef MIEEE
			/* 2**-53 =  1.11022302462515654042E-16 */
unsigned short MACHEP[4] = { 0x3ca0, 0x0000, 0x0000, 0x0000 };

#ifdef DENORMAL
			/* log(2**1024) =   7.09782712893383996843E2 */
unsigned short MAXLOG[4] = { 0x4086, 0x2e42, 0xfefa, 0x39ef };

			/* log(2**-1074) = - -7.44440071921381262314E2 */
/* unsigned short MINLOG[4] = {0xc087,0x4385,0x446d,0x71c3}; */
unsigned short MINLOG[4] = { 0xc087, 0x4910, 0xd52d, 0x3052 };
#else
			/* log(2**1022) =  7.08396418532264106224E2 */
unsigned short MAXLOG[4] = { 0x4086, 0x232b, 0xdd7a, 0xbcd2 };

			/* log(2**-1022) = - 7.08396418532264106224E2 */
unsigned short MINLOG[4] = { 0xc086, 0x232b, 0xdd7a, 0xbcd2 };
#endif
			/* 2**1024*(1-MACHEP) =  1.7976931348623158E308 */
unsigned short SQ2OPI[4] = { 0x3fe9, 0x8845, 0x33d4, 0x3651 };
unsigned short LOGSQ2[4] = { 0x3fd6, 0x2e42, 0xfefa, 0x39ef };
unsigned short THPIO4[4] = { 0x4002, 0xd97c, 0x7f33, 0x21d2 };

#endif

#ifndef UNK
extern unsigned short MACHEP[];
extern unsigned short MAXLOG[];
extern unsigned short UNDLOG[];
extern unsigned short MINLOG[];
extern unsigned short SQ2OPI[];
extern unsigned short LOGSQ2[];
extern unsigned short THPIO4[];
#endif
