/* ========================================================================== */
/* === umf_colamd.h ========================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*

Authors:

    The authors of the COLAMD code itself are Stefan I. Larimore and Timothy A.
    Davis, University of Florida.  The algorithm was developed in collaboration
    with John Gilbert, Xerox PARC, and Esmond Ng, Oak Ridge National Laboratory.

Date:

    UMFPACK Version: see above.
    COLAMD Version 2.0 was released on January 31, 2000.

Acknowledgements:

    This work was supported by the National Science Foundation, under
    grants DMS-9504974 and DMS-9803599.

UMFPACK:  Copyright (c) 2003 by Timothy A. Davis.  All Rights Reserved.

See the UMFPACK README file for the License for your use of this code.

Availability:

    Both UMFPACK and the original unmodified colamd/symamd library are
    available at http://www.cise.ufl.edu/research/sparse.

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

/* size of the knobs [ ] array.  Only knobs [0..2] are currently used. */
#define COLAMD_KNOBS 20

/* number of output statistics.  Only stats [0..8] are currently used. */
#define COLAMD_STATS 20

/* knobs [0] and stats [0]: dense row knob and output statistic. */
#define COLAMD_DENSE_ROW 0

/* knobs [1] and stats [1]: dense column knob and output statistic. */
#define COLAMD_DENSE_COL 1

/* knobs [2]: aggressive absorption option */
#define COLAMD_AGGRESSIVE 2

/* stats [2]: memory defragmentation count output statistic */
#define COLAMD_DEFRAG_COUNT 2

/* stats [3]: colamd status:  zero OK, > 0 warning or notice, < 0 error */
#define COLAMD_STATUS 3

/* stats [4..6]: error info, or info on jumbled columns */
#define COLAMD_INFO1 4
#define COLAMD_INFO2 5
#define COLAMD_INFO3 6

/* ------------------ */
/* added for UMFPACK: */
/* stats [7]: number of originally empty rows */
#define COLAMD_EMPTY_ROW 7
/* stats [8]: number of originally empty cols */
#define COLAMD_EMPTY_COL 8
/* stats [9]: number of rows with entries only in dense cols */
#define COLAMD_NEWLY_EMPTY_ROW 9
/* stats [10]: number of cols with entries only in dense rows */
#define COLAMD_NEWLY_EMPTY_COL 10
/* ------------------ */

/* error codes returned in stats [3]: */
#define COLAMD_OK				(0)
#define COLAMD_ERROR_jumbled_matrix		(-11)
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
    Int start ;		/* index for A of first row in this column, or DEAD */
			/* if column is dead */
    Int length ;	/* number of rows in this column */
    union
    {
	Int thickness ;	/* number of original columns represented by this */
			/* col, if the column is alive */
	Int parent ;	/* parent in parent tree super-column structure, if */
			/* the column is dead */
    } shared1 ;
    union
    {
	Int score ;	/* the score used to maintain heap, if col is alive */
	Int order ;	/* pivot ordering of this column, if col is dead */
    } shared2 ;
    union
    {
	Int headhash ;	/* head of a hash bucket, if col is at the head of */
			/* a degree list */
	Int hash ;	/* hash value, if col is not in a degree list */
	Int prev ;	/* previous column in degree list, if col is in a */
			/* degree list (but not at the head of a degree list) */
    } shared3 ;
    union
    {
	Int degree_next ;	/* next column, if col is in a degree list */
	Int hash_next ;		/* next column, if col is in a hash list */
    } shared4 ;

    /* ------------------ */
    /* added for UMFPACK: */
    Int nextcol ;	/* next column in this supercolumn */
    Int lastcol ;	/* last column in this supercolumn */
    /* ------------------ */

} Colamd_Col ;

typedef struct Colamd_Row_struct
{
    Int start ;		/* index for A of first col in this row */
    Int length ;	/* number of principal columns in this row */
    union
    {
	Int degree ;	/* number of principal & non-principal columns in row */
	Int p ;		/* used as a row pointer in init_rows_cols () */
    } shared1 ;
    union
    {
	Int mark ;	/* for computing set differences and marking dead rows*/
	Int first_column ;/* first column in row (used in garbage collection) */
    } shared2 ;

    /* ------------------ */
    /* added for UMFPACK: */
    Int thickness ;	/* number of original rows represented by this row */
			/* that are not yet pivotal */
    Int front ;		/* -1 if an original row */
			/* k if this row represents the kth frontal matrix */
			/* where k goes from 0 to at most n_col-1 */
    /* ------------------ */

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
*/

/* about 8*(n_col+1) integers: */
#define UMF_COLAMD_C(n_col) ((n_col + 1) * sizeof (Colamd_Col) / sizeof (Int))

/* about 6*(n_row+1) integers: */
#define UMF_COLAMD_R(n_row) ((n_row + 1) * sizeof (Colamd_Row) / sizeof (Int))

/* UMFPACK:  make sure Alen is >= 5*n_col + size of Col and Row structures.
 * Alen is typically about 2.2*nz + 9*n_col + 6*n_row, or 2.2nz+15n for
 * square matrices. */
#define UMF_COLAMD_RECOMMENDED(nnz, n_row, n_col)	\
(							\
((nnz) < 0 || (n_row) < 0 || (n_col) < 0)		\
?							\
    (-1)						\
:							\
    (MAX (2 * (nnz), 4 * (n_col)) +			\
    (Int) UMF_COLAMD_C (n_col) +			\
    (Int) UMF_COLAMD_R (n_row) + (n_col) + ((nnz) / 5))	\
)

/* ========================================================================== */
/* === Prototypes of user-callable routines ================================= */
/* ========================================================================== */

/* colamd_recommended removed for UMFPACK */

void UMF_colamd_set_defaults	/* sets default parameters */
(				/* knobs argument is modified on output */
    double knobs [COLAMD_KNOBS]	/* parameter settings for colamd */
) ;

Int UMF_colamd			/* returns (1) if successful, (0) otherwise*/
(				/* A and p arguments are modified on output */
    Int n_row,			/* number of rows in A */
    Int n_col,			/* number of columns in A */
    Int Alen,			/* size of the array A */
    Int A [],			/* row indices of A, of size Alen */
    Int p [],			/* column pointers of A, of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameter settings for colamd */
    Int stats [COLAMD_STATS]	/* colamd output statistics and error codes */
    /* ------------------ */
    /* added for UMFPACK: */
    , Int Front_npivcol [ ]
    , Int Front_nrows [ ]
    , Int Front_ncols [ ]
    , Int Front_parent [ ]
    , Int Front_cols [ ]
    , Int *p_nfr
    , Int InFront [ ]
    /* ------------------ */
) ;

/* symamd deleted for UMFPACK */

/* colamd_report deleted for UMFPACK */

/* symamd_report deleted for UMFPACK */

#endif /* COLAMD_H */
