/* 

   This file was considerably changed for UNU.RAN.
   As a consequence it has been renamed from "mconf.h" to
   "cephes_source.h" to avoid confusion.
   
   We were only interested in the files:
      gamma.c
      igamma.c
      incbet.c
      ndtr.c
      ndtri.c

   to use these files we needed the auxiliary files:
      isnan.c
      mtherr.c
      polevl.c

   The aim was to enhance portability for ANSI C without expert
   knowledge of the floating point unit.

   The main changes are
      numeric constants (in cephes - const.c) are replaced
      by macros defined in this file (source_mconf.h).

      We changed these constants using the constants provided by
      ANSI C in the math.h file.

      We changed also the definitions of MAXGAM (in incbet.c) and of
      MAXSTIR (in gamma.c).

      Everything concerning NANS and INFINITY was moved into isnan.c.
      (Also the defines to turn it off or on.)
      We have only tested this version with NAN and INFINITY turned
      off, but we have not changed anything concerning that
      question. So turning this support on should (could) work.


   March 26th, 2001, Josef Leydold and Wolfgang Hoermann.

*/

/*---------------------------------------------------------------------------*/

/*							mconf.h
 *
 *	Common include file for math routines
 *
 *
 *
 * SYNOPSIS:
 *
 * #include "mconf.h"
 *
 *
 *
 * DESCRIPTION:
 *
 * This file contains definitions for error codes that are
 * passed to the common error handling routine mtherr()
 * (which see).
 *
 * The file also includes a conditional assembly definition
 * for the type of computer arithmetic (IEEE, DEC, Motorola
 * IEEE, or UNKnown).
 * 
 * For Digital Equipment PDP-11 and VAX computers, certain
 * IBM systems, and others that use numbers with a 56-bit
 * significand, the symbol DEC should be defined.  In this
 * mode, most floating point constants are given as arrays
 * of octal integers to eliminate decimal to binary conversion
 * errors that might be introduced by the compiler.
 *
 * For little-endian computers, such as IBM PC, that follow the
 * IEEE Standard for Binary Floating Point Arithmetic (ANSI/IEEE
 * Std 754-1985), the symbol IBMPC should be defined.  These
 * numbers have 53-bit significands.  In this mode, constants
 * are provided as arrays of hexadecimal 16 bit integers.
 *
 * Big-endian IEEE format is denoted MIEEE.  On some RISC
 * systems such as Sun SPARC, double precision constants
 * must be stored on 8-byte address boundaries.  Since integer
 * arrays may be aligned differently, the MIEEE configuration
 * may fail on such machines.
 *
 * To accommodate other types of computer arithmetic, all
 * constants are also provided in a normal decimal radix
 * which one can hope are correctly converted to a suitable
 * format by the available C language compiler.  To invoke
 * this mode, define the symbol UNK.
 *
 * An important difference among these modes is a predefined
 * set of machine arithmetic constants for each.  The numbers
 * MACHEP (the machine roundoff error), MAXNUM (largest number
 * represented), and several other parameters are preset by
 * the configuration symbol.  Check the file const.c to
 * ensure that these values are correct for your computer.
 *
 * Configurations NANS, INFINITIES, MINUSZERO, and DENORMAL
 * may fail on many systems.  Verify that they are supposed
 * to work on your computer.
 */
/*
Cephes Math Library Release 2.3:  June, 1995
Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
*/

/*---------------------------------------------------------------------------*/
/* Include constants used in UNU.RAN.                                        */

#include <float.h>
#include <math.h>
#include <utils/unur_math_source.h>
#include <utils/unur_fp_source.h>
#include <utils/unur_fp_const_source.h>

/*---------------------------------------------------------------------------*/
/* Define mathematical constants (see x_math_source.h).                      */

#define PI      M_PI              /* Pi                                      */
#define SQRTH   M_SQRTH           /* sqrt(1/2)                               */

/*---------------------------------------------------------------------------*/
/* Define constant for floating point arithmetic.                            */

#define MAXNUM  DBL_MAX           /* largest number represented              */

/*---------------------------------------------------------------------------*/
/* Prototypes for functions from Cephes library                              */

/* gamma.c */
double _unur_cephes_gamma( double x );
double _unur_cephes_lgam( double x );

/* igam.c */
double _unur_cephes_igamc( double a, double x );
double _unur_cephes_igam( double a, double x );

/* incbet.c */
double _unur_cephes_incbet( double aa, double bb, double xx );

/* ndtr.c */
double _unur_cephes_ndtr( double a );
double _unur_cephes_erfc( double a );
double _unur_cephes_erf( double x );

/* ndtri.c */
double _unur_cephes_ndtri( double yval );

/* polevl.c */
double _unur_cephes_polevl( double x, double coef[], int N );
double _unur_cephes_p1evl( double x, double coef[], int N );

/*---------------------------------------------------------------------------*/
