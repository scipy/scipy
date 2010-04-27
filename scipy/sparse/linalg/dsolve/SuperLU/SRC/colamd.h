/*! @file colamd.h
    \brief Colamd prototypes and definitions

	<pre> 
    ==========================================================================
    === colamd/symamd prototypes and definitions =============================
    ==========================================================================

    You must include this file (colamd.h) in any routine that uses colamd,
    symamd, or the related macros and definitions.

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (davis@cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Date:

	September 8, 2003.  Version 2.3.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Notice:

	Copyright (c) 1998-2003 by the University of Florida.
	All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

	Permission is hereby granted to use, copy, modify, and/or distribute
	this program, provided that the Copyright, this License, and the
	Availability of the original version is retained on all copies and made
	accessible to the end-user of any code or package that includes COLAMD
	or any modified version of COLAMD. 

    Availability:

	The colamd/symamd library is available at

	    http://www.cise.ufl.edu/research/sparse/colamd/

	This is the http://www.cise.ufl.edu/research/sparse/colamd/colamd.h
	file.  It is required by the colamd.c, colamdmex.c, and symamdmex.c
	files, and by any C code that calls the routines whose prototypes are
	listed below, or that uses the colamd/symamd definitions listed below.
 </pre>
*/

#ifndef COLAMD_H
#define COLAMD_H

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include <stdlib.h>

/* ========================================================================== */
/* === Knob and statistics definitions ====================================== */
/* ========================================================================== */

/* size of the knobs [ ] array.  Only knobs [0..1] are currently used. */
#define COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..6] are currently used. */
#define COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define COLAMD_DENSE_COL 1

/* stats [2]: memory defragmentation count output statistic */
#define COLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */ 
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

/* error codes returned in stats [3]: */
#define COLAMD_OK				(0)
#define COLAMD_OK_BUT_JUMBLED			(1)
#define COLAMD_ERROR_A_not_present		(-1)
#define COLAMD_ERROR_p_not_present		(-2)
#define COLAMD_ERROR_nrow_negative		(-3)
#define COLAMD_ERROR_ncol_negative		(-4)
#define COLAMD_ERROR_nnz_negative		(-5)
#define COLAMD_ERROR_p0_nonzero			(-6)
#define COLAMD_ERROR_A_too_small		(-7)
#define COLAMD_ERROR_col_length_negative	(-8)
#define COLAMD_ERROR_row_index_out_of_bounds	(-9)
#define COLAMD_ERROR_out_of_memory		(-10)
#define COLAMD_ERROR_internal_error		(-999)

/* ========================================================================== */
/* === Row and Column structures ============================================ */
/* ========================================================================== */

/* User code that makes use of the colamd/symamd routines need not directly */
/* reference these structures.  They are used only for the COLAMD_RECOMMENDED */
/* macro. */

typedef struct Colamd_Col_struct
{
    int start ;		/* index for A of first row in this column, or DEAD */
			/* if column is dead */
    int length ;	/* number of rows in this column */
    union
    {
	int thickness ;	/* number of original columns represented by this */
			/* col, if the column is alive */
	int parent ;	/* parent in parent tree super-column structure, if */
			/* the column is dead */
    } shared1 ;
    union
    {
	int score ;	/* the score used to maintain heap, if col is alive */
	int order ;	/* pivot ordering of this column, if col is dead */
    } shared2 ;
    union
    {
	int headhash ;	/* head of a hash bucket, if col is at the head of */
			/* a degree list */
	int hash ;	/* hash value, if col is not in a degree list */
	int prev ;	/* previous column in degree list, if col is in a */
			/* degree list (but not at the head of a degree list) */
    } shared3 ;
    union
    {
	int degree_next ;	/* next column, if col is in a degree list */
	int hash_next ;		/* next column, if col is in a hash list */
    } shared4 ;

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    int start ;		/* index for A of first col in this row */
    int length ;	/* number of principal columns in this row */
    union
    {
	int degree ;	/* number of principal & non-principal columns in row */
	int p ;		/* used as a row pointer in init_rows_cols () */
    } shared1 ;
    union
    {
	int mark ;	/* for computing set differences and marking dead rows*/
	int first_column ;/* first column in row (used in garbage collection) */
    } shared2 ;

} Colamd_Row ;

/* ========================================================================== */
/* === Colamd recommended memory size ======================================= */
/* ========================================================================== */

/*
    The recommended length Alen of the array A passed to colamd is given by
    the COLAMD_RECOMMENDED (nnz, n_row, n_col) macro.  It returns -1 if any
    argument is negative.  2*nnz space is required for the row and column
    indices of the matrix. COLAMD_C (n_col) + COLAMD_R (n_row) space is
    required for the Col and Row arrays, respectively, which are internal to
    colamd.  An additional n_col space is the minimal amount of "elbow room",
    and nnz/5 more space is recommended for run time efficiency.

    This macro is not needed when using symamd.

    Explicit typecast to int added Sept. 23, 2002, COLAMD version 2.2, to avoid
    gcc -pedantic warning messages.
*/

#define COLAMD_C(n_col) ((int) (((n_col) + 1) * sizeof (Colamd_Col) / sizeof (int)))
#define COLAMD_R(n_row) ((int) (((n_row) + 1) * sizeof (Colamd_Row) / sizeof (int)))

#define COLAMD_RECOMMENDED(nnz, n_row, n_col)                                 \
(                                                                             \
((nnz) < 0 || (n_row) < 0 || (n_col) < 0)                                     \
?                                                                             \
    (-1)                                                                      \
:                                                                             \
    (2 * (nnz) + COLAMD_C (n_col) + COLAMD_R (n_row) + (n_col) + ((nnz) / 5)) \
)

/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

int colamd_recommended		/* returns recommended value of Alen, */
				/* or (-1) if input arguments are erroneous */
(
    int nnz,			/* nonzeros in A */
    int n_row,			/* number of rows in A */
    int n_col			/* number of columns in A */
) ;

void colamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [COLAMD_KNOBS]	/* parameter settings for colamd */
) ;

int colamd			/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    int n_row,			/* number of rows in A */
    int n_col,			/* number of columns in A */
    int Alen,			/* size of the array A */
    int A [],			/* row indices of A, of size Alen */
    int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameter settings for colamd */
    int stats [COLAMD_STATS]	/* colamd output statistics and error codes */
) ;

int symamd				/* return (1) if OK, (0) otherwise */
(
    int n,				/* number of rows and columns of A */
    int A [],				/* row indices of A */
    int p [],				/* column pointers of A */
    int perm [],			/* output permutation, size n_col+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    int stats [COLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for MATLAB mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for MATLAB mexFunction) */
) ;

void colamd_report
(
    int stats [COLAMD_STATS]
) ;

void symamd_report
(
    int stats [COLAMD_STATS]
) ;

#endif /* COLAMD_H */
