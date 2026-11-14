//------------------------------------------------------------------------------
// AMD/Include/amd.h:  approximate minimum degree ordering
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2024, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* AMD finds a symmetric ordering P of a matrix A so that the Cholesky
 * factorization of P*A*P' has fewer nonzeros and takes less work than the
 * Cholesky factorization of A.  If A is not symmetric, then it performs its
 * ordering on the matrix A+A'.  Two sets of user-callable routines are
 * provided, one for int32_t integers and the other for int64_t integers.
 *
 * The method is based on the approximate minimum degree algorithm, discussed
 * in Amestoy, Davis, and Duff, "An approximate degree ordering algorithm",
 * SIAM Journal of Matrix Analysis and Applications, vol. 17, no. 4, pp.
 * 886-905, 1996.  This package can perform both the AMD ordering (with
 * aggressive absorption), and the AMDBAR ordering (without aggressive
 * absorption) discussed in the above paper.  This package differs from the
 * Fortran codes discussed in the paper:
 *
 *       (1) it can ignore "dense" rows and columns, leading to faster run times
 *       (2) it computes the ordering of A+A' if A is not symmetric
 *       (3) it is followed by a depth-first post-ordering of the assembly tree
 *           (or supernodal elimination tree)
 *
 * For historical reasons, the Fortran versions, amd.f and amdbar.f, have
 * been left (nearly) unchanged.  They compute the identical ordering as
 * described in the above paper.
 */

#ifndef HIGHS_AMD_ORDERING_H
#define HIGHS_AMD_ORDERING_H

#include "SuiteSparse_config.h"
#include "HConfig.h"

#ifdef HIGHSINT64
typedef int64_t amd_int;
typedef uint64_t amd_uint;
#define amd_int_max INT64_MAX
#else
typedef int amd_int;
typedef unsigned int amd_uint;
#define amd_int_max INT_MAX
#endif

/* make it easy for C++ programs to include AMD */
#ifdef __cplusplus
extern "C" {
#endif

int Highs_amd_order  /* returns AMD_OK, AMD_OK_BUT_JUMBLED,
                                    * AMD_INVALID, or AMD_OUT_OF_MEMORY */
(
    amd_int n,                         /* A is n-by-n.  n must be >= 0. */
    const amd_int Ap [ ],              /* column pointers for A, of size n+1 */
    const amd_int Ai [ ],              /* row indices of A, of size nz = Ap [n] */
    amd_int P [ ],                     /* output permutation, of size n */
    double Control [ ],     /* input Control settings, of size AMD_CONTROL */
    double Info [ ]         /* output Info statistics, of size AMD_INFO */
) ;

/* Input arguments (not modified):
 *
 *      n: the matrix A is n-by-n.
 *      Ap: an int32_t/int64_t array of size n+1, containing column
 *          pointers of A.
 *      Ai: an int32_t/int64_t array of size nz, containing the row
 *          indices of A, where nz = Ap [n].
 *      Control:  a double array of size AMD_CONTROL, containing control
 *          parameters.  Defaults are used if Control is NULL.
 *
 * Output arguments (not defined on input):
 *
 *      P: an int32_t/int64_t array of size n, containing the output
 *          permutation. If row i is the kth pivot row, then P [k] = i.  In
 *          MATLAB notation, the reordered matrix is A (P,P).
 *      Info: a double array of size AMD_INFO, containing statistical
 *          information.  Ignored if Info is NULL.
 *
 * On input, the matrix A is stored in column-oriented form.  The row indices
 * of nonzero entries in column j are stored in Ai [Ap [j] ... Ap [j+1]-1].
 *
 * If the row indices appear in ascending order in each column, and there
 * are no duplicate entries, then amd_order is slightly more efficient in
 * terms of time and memory usage.  If this condition does not hold, a copy
 * of the matrix is created (where these conditions do hold), and the copy is
 * ordered.
 * 
 * Row indices must be in the range 0 to
 * n-1.  Ap [0] must be zero, and thus nz = Ap [n] is the number of nonzeros
 * in A.  The array Ap is of size n+1, and the array Ai is of size nz = Ap [n].
 * The matrix does not need to be symmetric, and the diagonal does not need to
 * be present (if diagonal entries are present, they are ignored except for
 * the output statistic Info [AMD_NZDIAG]).  The arrays Ai and Ap are not
 * modified.  This form of the Ap and Ai arrays to represent the nonzero
 * pattern of the matrix A is the same as that used internally by MATLAB.
 * If you wish to use a more flexible input structure, please see the
 * umfpack_*_triplet_to_col routines in the UMFPACK package, at
 * http://www.suitesparse.com.
 *
 * Restrictions:  n >= 0.  Ap [0] = 0.  Ap [j] <= Ap [j+1] for all j in the
 *      range 0 to n-1.  nz = Ap [n] >= 0.  Ai [0..nz-1] must be in the range 0
 *      to n-1.  Finally, Ai, Ap, and P must not be NULL.  If any of these
 *      restrictions are not met, AMD returns AMD_INVALID.
 *
 * AMD returns:
 *
 *      AMD_OK if the matrix is valid and sufficient memory can be allocated to
 *          perform the ordering.
 *
 *      AMD_OUT_OF_MEMORY if not enough memory can be allocated.
 *
 *      AMD_INVALID if the input arguments n, Ap, Ai are invalid, or if P is
 *          NULL.
 *
 *      AMD_OK_BUT_JUMBLED if the matrix had unsorted columns, and/or duplicate
 *          entries, but was otherwise valid.
 *
 * The AMD routine first forms the pattern of the matrix A+A', and then
 * computes a fill-reducing ordering, P.  If P [k] = i, then row/column i of
 * the original is the kth pivotal row.  In MATLAB notation, the permuted
 * matrix is A (P,P), except that 0-based indexing is used instead of the
 * 1-based indexing in MATLAB.
 *
 * The Control array is used to set various parameters for AMD.  If a NULL
 * pointer is passed, default values are used.  The Control array is not
 * modified.
 *
 *      Control [AMD_DENSE]:  controls the threshold for "dense" rows/columns.
 *          A dense row/column in A+A' can cause AMD to spend a lot of time in
 *          ordering the matrix.  If Control [AMD_DENSE] >= 0, rows/columns
 *          with more than Control [AMD_DENSE] * sqrt (n) entries are ignored
 *          during the ordering, and placed last in the output order.  The
 *          default value of Control [AMD_DENSE] is 10.  If negative, no
 *          rows/columns are treated as "dense".  Rows/columns with 16 or
 *          fewer off-diagonal entries are never considered "dense".
 *
 *      Control [AMD_AGGRESSIVE]: controls whether or not to use aggressive
 *          absorption, in which a prior element is absorbed into the current
 *          element if is a subset of the current element, even if it is not
 *          adjacent to the current pivot element (refer to Amestoy, Davis,
 *          & Duff, 1996, for more details).  The default value is nonzero,
 *          which means to perform aggressive absorption.  This nearly always
 *          leads to a better ordering (because the approximate degrees are
 *          more accurate) and a lower execution time.  There are cases where
 *          it can lead to a slightly worse ordering, however.  To turn it off,
 *          set Control [AMD_AGGRESSIVE] to 0.
 *
 *      Control [2..4] are not used in the current version, but may be used in
 *          future versions.
 *
 * The Info array provides statistics about the ordering on output.  If it is
 * not present, the statistics are not returned.  This is not an error
 * condition.
 * 
 *      Info [AMD_STATUS]:  the return value of AMD, either AMD_OK,
 *          AMD_OK_BUT_JUMBLED, AMD_OUT_OF_MEMORY, or AMD_INVALID.
 *
 *      Info [AMD_N]: n, the size of the input matrix
 *
 *      Info [AMD_NZ]: the number of nonzeros in A, nz = Ap [n]
 *
 *      Info [AMD_SYMMETRY]:  the symmetry of the matrix A.  It is the number
 *          of "matched" off-diagonal entries divided by the total number of
 *          off-diagonal entries.  An entry A(i,j) is matched if A(j,i) is also
 *          an entry, for any pair (i,j) for which i != j.  In MATLAB notation,
 *              S = spones (A) ;
 *              B = tril (S, -1) + triu (S, 1) ;
 *              symmetry = nnz (B & B') / nnz (B) ;
 *
 *      Info [AMD_NZDIAG]: the number of entries on the diagonal of A.
 *
 *      Info [AMD_NZ_A_PLUS_AT]:  the number of nonzeros in A+A', excluding the
 *          diagonal.  If A is perfectly symmetric (Info [AMD_SYMMETRY] = 1)
 *          with a fully nonzero diagonal, then Info [AMD_NZ_A_PLUS_AT] = nz-n
 *          (the smallest possible value).  If A is perfectly unsymmetric
 *          (Info [AMD_SYMMETRY] = 0, for an upper triangular matrix, for
 *          example) with no diagonal, then Info [AMD_NZ_A_PLUS_AT] = 2*nz
 *          (the largest possible value).
 *
 *      Info [AMD_NDENSE]: the number of "dense" rows/columns of A+A' that were
 *          removed from A prior to ordering.  These are placed last in the
 *          output order P.
 *
 *      Info [AMD_MEMORY]: the amount of memory used by AMD, in bytes.  In the
 *          current version, this is 1.2 * Info  [AMD_NZ_A_PLUS_AT] + 9*n
 *          times the size of an integer.  This is at most 2.4nz + 9n.  This
 *          excludes the size of the input arguments Ai, Ap, and P, which have
 *          a total size of nz + 2*n + 1 integers.
 *
 *      Info [AMD_NCMPA]: the number of garbage collections performed.
 *
 *      Info [AMD_LNZ]: the number of nonzeros in L (excluding the diagonal).
 *          This is a slight upper bound because mass elimination is combined
 *          with the approximate degree update.  It is a rough upper bound if
 *          there are many "dense" rows/columns.  The rest of the statistics,
 *          below, are also slight or rough upper bounds, for the same reasons.
 *          The post-ordering of the assembly tree might also not exactly
 *          correspond to a true elimination tree postordering.
 *
 *      Info [AMD_NDIV]: the number of divide operations for a subsequent LDL'
 *          or LU factorization of the permuted matrix A (P,P).
 *
 *      Info [AMD_NMULTSUBS_LDL]:  the number of multiply-subtract pairs for a
 *          subsequent LDL' factorization of A (P,P).
 *
 *      Info [AMD_NMULTSUBS_LU]:  the number of multiply-subtract pairs for a
 *          subsequent LU factorization of A (P,P), assuming that no numerical
 *          pivoting is required.
 *
 *      Info [AMD_DMAX]:  the maximum number of nonzeros in any column of L,
 *          including the diagonal.
 *
 *      Info [14..19] are not used in the current version, but may be used in
 *          future versions.
 */    

/* ------------------------------------------------------------------------- */
/* direct interface to AMD */
/* ------------------------------------------------------------------------- */

/* amd_2 is the primary AMD ordering routine.  It is not meant to be
 * user-callable because of its restrictive inputs and because it destroys
 * the user's input matrix.  It does not check its inputs for errors, either.
 * However, if you can work with these restrictions it can be faster than
 * amd_order and use less memory (assuming that you can create your own copy
 * of the matrix for AMD to destroy).  Refer to AMD/Source/amd_2.c for a
 * description of each parameter. */

void Highs_amd_2
(
    amd_int n,
    amd_int Pe [ ],
    amd_int Iw [ ],
    amd_int Len [ ],
    amd_int iwlen,
    amd_int pfree,
    amd_int Nv [ ],
    amd_int Next [ ], 
    amd_int Last [ ],
    amd_int Head [ ],
    amd_int Elen [ ],
    amd_int Degree [ ],
    amd_int W [ ],
    double Control [ ],
    double Info [ ]
) ;

/* ------------------------------------------------------------------------- */
/* amd_valid */
/* ------------------------------------------------------------------------- */

/* Returns AMD_OK or AMD_OK_BUT_JUMBLED if the matrix is valid as input to
 * amd_order; the latter is returned if the matrix has unsorted and/or
 * duplicate row indices in one or more columns.  Returns AMD_INVALID if the
 * matrix cannot be passed to amd_order.  For amd_order, the matrix must also
 * be square.  The first two arguments are the number of rows and the number
 * of columns of the matrix.  For its use in AMD, these must both equal n.
 */

int Highs_amd_valid
(
    amd_int n_row,                 /* # of rows */
    amd_int n_col,                 /* # of columns */
    const amd_int Ap [ ],          /* column pointers, of size n_col+1 */
    const amd_int Ai [ ]           /* row indices, of size Ap [n_col] */
) ;

/* ------------------------------------------------------------------------- */
/* AMD Control and Info arrays */
/* ------------------------------------------------------------------------- */

/* amd_defaults:  sets the default control settings */
void Highs_amd_defaults   (double Control [ ]) ;

/* amd_control: prints the control settings */
void Highs_amd_control    (double Control [ ]) ;

/* amd_info: prints the statistics */
void Highs_amd_info       (double Info [ ]) ;

#ifdef __cplusplus
}
#endif

#define AMD_CONTROL 5          /* size of Control array */
#define AMD_INFO 20            /* size of Info array */

/* contents of Control */
#define AMD_DENSE 0            /* "dense" if degree > Control [0] * sqrt (n) */
#define AMD_AGGRESSIVE 1    /* do aggressive absorption if Control [1] != 0 */

/* default Control settings */
#define AMD_DEFAULT_DENSE 10.0          /* default "dense" degree 10*sqrt(n) */
#define AMD_DEFAULT_AGGRESSIVE 1    /* do aggressive absorption by default */

/* contents of Info */
#define AMD_STATUS 0           /* return value of amd_order and amd_l_order */
#define AMD_N 1                /* A is n-by-n */
#define AMD_NZ 2      /* number of nonzeros in A */ 
#define AMD_SYMMETRY 3         /* symmetry of pattern (1 is sym., 0 is unsym.) */
#define AMD_NZDIAG 4           /* # of entries on diagonal */
#define AMD_NZ_A_PLUS_AT 5  /* nz in A+A' */
#define AMD_NDENSE 6           /* number of "dense" rows/columns in A */
#define AMD_MEMORY 7           /* amount of memory used by AMD */
#define AMD_NCMPA 8            /* number of garbage collections in AMD */
#define AMD_LNZ 9     /* approx. nz in L, excluding the diagonal */
#define AMD_NDIV 10            /* number of fl. point divides for LU and LDL' */
#define AMD_NMULTSUBS_LDL 11 /* number of fl. point (*,-) pairs for LDL' */
#define AMD_NMULTSUBS_LU 12  /* number of fl. point (*,-) pairs for LU */
#define AMD_DMAX 13             /* max nz. in any column of L, incl. diagonal */

/* ------------------------------------------------------------------------- */
/* return values of AMD */
/* ------------------------------------------------------------------------- */

#define AMD_OK 0           /* success */
#define AMD_OUT_OF_MEMORY -1        /* malloc failed, or problem too large */
#define AMD_INVALID -2              /* input arguments are not valid */
#define AMD_OK_BUT_JUMBLED 1        /* input matrix is OK for amd_order, but
    * columns were not sorted, and/or duplicate entries were present.  AMD had
    * to do extra work before ordering the matrix.  This is a warning, not an
    * error.  */

/* ========================================================================== */
/* === AMD version ========================================================== */
/* ========================================================================== */

/* AMD Version 1.2 and later include the following definitions.
 * As an example, to test if the version you are using is 1.2 or later:
 *
 * #ifdef AMD_VERSION
 *       if (AMD_VERSION >= AMD_VERSION_CODE (1,2)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *       #if defined(AMD_VERSION) && (AMD_VERSION >= AMD_VERSION_CODE (1,2))
 *           // This is version 1.2 or later
 *       #else
 *           // This is an early version
 *       #endif
 *
 * Versions 1.1 and earlier of AMD do not include a #define'd version number.
 */

#define AMD_DATE "July 25, 2025"
#define AMD_MAIN_VERSION   3
#define AMD_SUB_VERSION    3
#define AMD_SUBSUB_VERSION 4

#endif
