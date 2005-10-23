/* ========================================================================== */
/* === umfpack_report_matrix ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_report_matrix
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ],
    int col_form,
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_dl_report_matrix
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ],
    long col_form,
    const double Control [UMFPACK_CONTROL]
) ;

int umfpack_zi_report_matrix
(
    int n_row,
    int n_col,
    const int Ap [ ],
    const int Ai [ ],
    const double Ax [ ], const double Az [ ],
    int col_form,
    const double Control [UMFPACK_CONTROL]
) ;

long umfpack_zl_report_matrix
(
    long n_row,
    long n_col,
    const long Ap [ ],
    const long Ai [ ],
    const double Ax [ ], const double Az [ ],
    long col_form,
    const double Control [UMFPACK_CONTROL]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    int n_row, n_col, *Ap, *Ai, status ;
    double *Ax, Control [UMFPACK_CONTROL] ;
    status = umfpack_di_report_matrix (n_row, n_col, Ap, Ai, Ax, 1, Control) ;
or:
    status = umfpack_di_report_matrix (n_row, n_col, Ap, Ai, Ax, 0, Control) ;

double long Syntax:

    #include "umfpack.h"
    long n_row, n_col, *Ap, *Ai, status ;
    double *Ax, Control [UMFPACK_CONTROL] ;
    status = umfpack_dl_report_matrix (n_row, n_col, Ap, Ai, Ax, 1, Control) ;
or:
    status = umfpack_dl_report_matrix (n_row, n_col, Ap, Ai, Ax, 0, Control) ;

complex int Syntax:

    #include "umfpack.h"
    int n_row, n_col, *Ap, *Ai, status ;
    double *Ax, *Az, Control [UMFPACK_CONTROL] ;
    status = umfpack_zi_report_matrix (n_row, n_col, Ap, Ai, Ax, Az, 1,
        Control) ;
or:
    status = umfpack_zi_report_matrix (n_row, n_col, Ap, Ai, Ax, Az, 0,
        Control) ;

complex long Syntax:

    #include "umfpack.h"
    long n_row, n_col, *Ap, *Ai, status ;
    double *Ax, Control [UMFPACK_CONTROL] ;
    status = umfpack_zl_report_matrix (n_row, n_col, Ap, Ai, Ax, Az, 1,
	Control) ;
or:
    status = umfpack_zl_report_matrix (n_row, n_col, Ap, Ai, Ax, Az, 0,
	Control) ;

Purpose:

    Verifies and prints a row or column-oriented sparse matrix.

Returns:

    UMFPACK_OK if Control [UMFPACK_PRL] <= 2 (the input is not checked).

    Otherwise (where n is n_col for the column form and n_row for row
    and let ni be n_row for the column form and n_col for row):

    UMFPACK_OK if the matrix is valid.

    UMFPACK_ERROR_n_nonpositive if n_row <= 0 or n_col <= 0.
    UMFPACK_ERROR_argument_missing if Ap and/or Ai are missing.
    UMFPACK_ERROR_invalid_matrix if Ap [n] < 0, if Ap [0] is not zero,
	if Ap [j+1] < Ap [j] for any j in the range 0 to n-1,
	if any row index in Ai is not in the range 0 to ni-1, or
	if the row indices in any column are not in
	ascending order, or contain duplicates.
    UMFPACK_ERROR_out_of_memory if out of memory.

Arguments:

    Int n_row ;		Input argument, not modified.
    Int n_col ;		Input argument, not modified.

	A is an n_row-by-n_row matrix.  Restriction: n_row > 0 and n_col > 0.

    Int Ap [n+1] ;	Input argument, not modified.

	n is n_row for a row-form matrix, and n_col for a column-form matrix.

	Ap is an integer array of size n+1.  If col_form is true (nonzero),
	then on input, it holds the "pointers" for the column form of the
	sparse matrix A.  The row indices of column j of the matrix A are held
	in Ai [(Ap [j]) ... (Ap [j+1]-1)].  Otherwise, Ap holds the
	row pointers, and the column indices of row j of the matrix are held
	in Ai [(Ap [j]) ... (Ap [j+1]-1)].

	The first entry, Ap [0], must be zero, and Ap [j] <= Ap [j+1] must hold
	for all j in the range 0 to n-1.  The value nz = Ap [n] is thus the
	total number of entries in the pattern of the matrix A.

    Int Ai [nz] ;	Input argument, not modified, of size nz = Ap [n].

	If col_form is true (nonzero), then the nonzero pattern (row indices)
	for column j is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)].  Row indices
	must be in the range 0 to n_row-1 (the matrix is 0-based).

	Otherwise, the nonzero pattern (column indices) for row j is stored in
	Ai [(Ap [j]) ... (Ap [j+1]-1)]. Column indices must be in the range 0
	to n_col-1 (the matrix is 0-based).

    double Ax [nz] ;	Input argument, not modified, of size nz = Ap [n].

	The numerical values of the sparse matrix A.

	If col_form is true (nonzero), then the nonzero pattern (row indices)
	for column j is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)], and the
	corresponding (real) numerical values are stored in
	Ax [(Ap [j]) ... (Ap [j+1]-1)].  The imaginary parts are stored in
	Az [(Ap [j]) ... (Ap [j+1]-1)], for the complex versions.

	Otherwise, the nonzero pattern (column indices) for row j
	is stored in Ai [(Ap [j]) ... (Ap [j+1]-1)], and the corresponding
	(real) numerical values are stored in Ax [(Ap [j]) ... (Ap [j+1]-1)].
	The imaginary parts are stored in Az [(Ap [j]) ... (Ap [j+1]-1)],
	for the complex versions.

	No numerical values are printed if Ax or Az are (double *) NULL.

    double Az [nz] ;	Input argument, not modified, for complex versions.

	The imaginary values of the sparse matrix A.   See the description
	of Ax, above.  No numerical values are printed if Az is NULL. 

	Future complex version:  if Ax is present and Az is NULL, then both real
	and imaginary parts will be contained in Ax[0..2*nz-1], with Ax[2*k]
	and Ax[2*k+1] being the real and imaginary part of the kth entry.

    Int col_form ;	Input argument, not modified.

	The matrix is in row-oriented form if form is col_form is false (0).
	Otherwise, the matrix is in column-oriented form.

    double Control [UMFPACK_CONTROL] ;	Input argument, not modified.

	If a (double *) NULL pointer is passed, then the default control
	settings are used.  Otherwise, the settings are determined from the
	Control array.  See umfpack_*_defaults on how to fill the Control
	array with the default settings.  If Control contains NaN's, the
	defaults are used.  The following Control parameters are used:

	Control [UMFPACK_PRL]:  printing level.

	    2 or less: no output.  returns silently without checking anything.
	    3: fully check input, and print a short summary of its status
	    4: as 3, but print first few entries of the input
	    5: as 3, but print all of the input
	    Default: 1
*/
