/*                                                     mconf.h
 *
 *     Common include file for math routines
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
 * The file includes a conditional assembly definition for the type of
 * computer arithmetic (IEEE, Motorola IEEE, or UNKnown).
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
 * Cephes Math Library Release 2.3:  June, 1995
 * Copyright 1984, 1987, 1989, 1995 by Stephen L. Moshier
 */

#ifndef CEPHES_MCONF_H
#define CEPHES_MCONF_H

#include <Python.h>
#include <math.h>

#include "cephes_names.h"
#include "cephes.h"
#include "polevl.h"
#include "sf_error.h"

#define MAXITER        500
#define EDOM		33
#define ERANGE		34

/* Type of computer arithmetic */

/* UNKnown arithmetic, invokes coefficients given in
 * normal decimal format.  Beware of range boundary
 * problems (MACHEP, MAXLOG, etc. in const.c) and
 * roundoff problems in pow.c:
 * (Sun SPARCstation)
 */

/* SciPy note: by defining UNK, we prevent the compiler from
 * casting integers to floating point numbers.  If the Endianness
 * is detected incorrectly, this causes problems on some platforms.
 */
#define UNK 1

/* Define to support tiny denormal numbers, else undefine. */
#define DENORMAL 1

#define gamma Gamma

/*
 * Enable loop unrolling on GCC and use faster isnan et al.
 */
#if !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
#pragma GCC optimize("unroll-loops")
#define cephes_isnan(x) __builtin_isnan(x)
#define cephes_isinf(x) __builtin_isinf(x)
#define cephes_isfinite(x) __builtin_isfinite(x)
#endif
#endif
#ifndef cephes_isnan
#define cephes_isnan(x) isnan(x)
#define cephes_isinf(x) isinf(x)
#define cephes_isfinite(x) isfinite(x)
#endif

/* Constants needed that are not available in the C standard library */
#define SCIPY_EULER     0.577215664901532860606512090082402431   /* Euler constant */
#define SCIPY_El        2.718281828459045235360287471352662498L  /* e as long double */

#endif				/* CEPHES_MCONF_H */
