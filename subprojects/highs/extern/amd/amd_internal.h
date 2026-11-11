//------------------------------------------------------------------------------
// AMD/Include/amd_internal.h: internal definitions for AMD
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2023, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* This file is for internal use in AMD itself, and does not normally need to
 * be included in user code (it is included in UMFPACK, however).   All others
 * should use amd.h instead.
 */

#include "amd.h"

/* ------------------------------------------------------------------------- */
/* basic definitions */
/* ------------------------------------------------------------------------- */

#ifdef FLIP
#undef FLIP
#endif

#ifdef MAX
#undef MAX
#endif

#ifdef MIN
#undef MIN
#endif

#ifdef EMPTY
#undef EMPTY
#endif

#define PRIVATE static

/* FLIP is a "negation about -1", and is used to mark an integer i that is
 * normally non-negative.  FLIP (EMPTY) is EMPTY.  FLIP of a number > EMPTY
 * is negative, and FLIP of a number < EMTPY is positive.  FLIP (FLIP (i)) = i
 * for all integers i.  UNFLIP (i) is >= EMPTY. */
#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP (i) : (i))

/* for integer MAX/MIN, or for doubles when we don't care how NaN's behave: */
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* logical expression of p implies q: */
#define IMPLIES(p,q) (!(p) || (q))

/* Note that the IBM RS 6000 xlc predefines TRUE and FALSE in <types.h>. */
/* The Compaq Alpha also predefines TRUE and FALSE. */
#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

#define TRUE (1)
#define FALSE (0)
#define EMPTY (-1)

/* largest value of size_t */
#ifndef SIZE_T_MAX
#ifdef SIZE_MAX
/* C99 only */
#define SIZE_T_MAX SIZE_MAX
#else
#define SIZE_T_MAX ((size_t) (-1))
#endif
#endif

/* ------------------------------------------------------------------------- */
/* AMD routine definitions (not user-callable) */
/* ------------------------------------------------------------------------- */

size_t Highs_amd_aat
(
    amd_int n,
    const amd_int Ap [ ],
    const amd_int Ai [ ],
    amd_int Len [ ],
    amd_int Tp [ ],
    double Info [ ]
) ;

void Highs_amd_1
(
    amd_int n,
    const amd_int Ap [ ],
    const amd_int Ai [ ],
    amd_int P [ ],
    amd_int Pinv [ ],
    amd_int Len [ ],
    amd_int slen,
    amd_int S [ ],
    double Control [ ],
    double Info [ ]
) ;

void Highs_amd_postorder
(
    amd_int nn,
    amd_int Parent [ ],
    amd_int Npiv [ ],
    amd_int Fsize [ ],
    amd_int Order [ ],
    amd_int Child [ ],
    amd_int Sibling [ ],
    amd_int Stack [ ]
) ;

amd_int Highs_amd_post_tree
(
    amd_int root,
    amd_int k,
    amd_int Child [ ],
    const amd_int Sibling [ ],
    amd_int Order [ ],
    amd_int Stack [ ]
) ;

void Highs_amd_preprocess
(
    amd_int n,
    const amd_int Ap [ ],
    const amd_int Ai [ ],
    amd_int Rp [ ],
    amd_int Ri [ ],
    amd_int W [ ],
    amd_int Flag [ ]
) ;
