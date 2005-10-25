/* ========================================================================== */
/* === umfpack_get_numeric ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_get_numeric
(
    int Lp [ ],
    int Lj [ ],
    double Lx [ ],
    int Up [ ],
    int Ui [ ],
    double Ux [ ],
    int P [ ],
    int Q [ ],
    double Dx [ ],
    int *do_recip,
    double Rs [ ],
    void *Numeric
) ;

long umfpack_dl_get_numeric
(
    long Lp [ ],
    long Lj [ ],
    double Lx [ ],
    long Up [ ],
    long Ui [ ],
    double Ux [ ],
    long P [ ],
    long Q [ ],
    double Dx [ ],
    long *do_recip,
    double Rs [ ],
    void *Numeric
) ;

int umfpack_zi_get_numeric
(
    int Lp [ ],
    int Lj [ ],
    double Lx [ ], double Lz [ ],
    int Up [ ],
    int Ui [ ],
    double Ux [ ], double Uz [ ],
    int P [ ],
    int Q [ ],
    double Dx [ ], double Dz [ ],
    int *do_recip,
    double Rs [ ],
    void *Numeric
) ;

long umfpack_zl_get_numeric
(
    long Lp [ ],
    long Lj [ ],
    double Lx [ ], double Lz [ ],
    long Up [ ],
    long Ui [ ],
    double Ux [ ], double Uz [ ],
    long P [ ],
    long Q [ ],
    double Dx [ ], double Dz [ ],
    long *do_recip,
    double Rs [ ],
    void *Numeric
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int *Lp, *Lj, *Up, *Ui, *P, *Q, status, do_recip ;
    double *Lx, *Ux, *Dx, *Rs ;
    status = umfpack_di_get_numeric (Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx,
	&do_recip, Rs, Numeric) ;

double long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long *Lp, *Lj, *Up, *Ui, *P, *Q, status, do_recip ;
    double *Lx, *Ux, *Dx, *Rs ;
    status = umfpack_dl_get_numeric (Lp, Lj, Lx, Up, Ui, Ux, P, Q, Dx,
	&do_recip, Rs, Numeric) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int *Lp, *Lj, *Up, *Ui, *P, *Q, status, do_recip ;
    double *Lx, *Lz, *Ux, *Uz, *Dx, *Dz, *Rs ;
    status = umfpack_zi_get_numeric (Lp, Lj, Lx, Lz, Up, Ui, Ux, Uz, P, Q,
	Dx, Dz, &do_recip, Rs, Numeric) ;

complex long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    long *Lp, *Lj, *Up, *Ui, *P, *Q, status, do_recip ;
    double *Lx, *Lz, *Ux, *Uz, *Dx, *Dz, *Rs ;
    status = umfpack_zl_get_numeric (Lp, Lj, Lx, Lz, Up, Ui, Ux, Uz, P, Q,
	Dx, Dz, &do_recip, Rs, Numeric) ;

Purpose:

    This routine copies the LU factors and permutation vectors from the Numeric
    object into user-accessible arrays.  This routine is not needed to solve a
    linear system.  Note that the output arrays Lp, Lj, Lx, Up, Ui, Ux, P, Q,
    Dx, and Rs are not allocated by umfpack_*_get_numeric; they must exist on
    input.  Similarly, Lz, Uz and Dz must exist on input for the complex
    versions.

    All output arguments are optional.  If any of them are NULL
    on input, then that part of the LU factorization is not copied.  You can
    use this routine to extract just the parts of the LU factorization that
    you want.  For example, to retrieve just the column permutation Q, use:

    #define noD (double *) NULL
    #define noI (int *) NULL
    status = umfpack_di_get_numeric (noI, noI, noD, noI, noI, noD, noI,
	Q, noD, noI, noD, Numeric) ;

Returns:

    Returns UMFPACK_OK if successful.  Returns UMFPACK_ERROR_out_of_memory
    if insufficient memory is available for the 2*max(n_row,n_col) integer
    workspace that umfpack_*_get_numeric allocates to construct L and/or U.
    Returns UMFPACK_ERROR_invalid_Numeric_object if the Numeric object provided
    as input is invalid.

Arguments:

    Int Lp [n_row+1] ;	Output argument.
    Int Lj [lnz] ;	Output argument.
    double Lx [lnz] ;	Output argument.
    double Lz [lnz] ;	Output argument for complex versions.

	The n_row-by-min(n_row,n_col) matrix L is returned in compressed-row
	form.  The column indices of row i and corresponding numerical values
	are in:

	    Lj [Lp [i] ... Lp [i+1]-1]
	    Lx [Lp [i] ... Lp [i+1]-1]	real part
	    Lz [Lp [i] ... Lp [i+1]-1]	imaginary part (complex versions)

	respectively.  Each row is stored in sorted order, from low column
	indices to higher.  The last entry in each row is the diagonal, which
	is numerically equal to one.  The sizes of Lp, Lj, Lx, and Lz are
	returned by umfpack_*_get_lunz.    If Lp, Lj, or Ux (or Uz for the
	complex version) are not present, then the matrix L is not returned.
	This is not an error condition.  The L matrix can be printed if n_row,
	Lp, Lj, Lx (and Lz for the complex versions) are passed to
	umfpack_*_report_matrix (using the "row" form).

	Future complex version:  if Lx is present and Lz is NULL, then both real
	and imaginary parts will be returned in Lx[0..2*lnz-1], with Lx[2*k]
	and Lx[2*k+1] being the real and imaginary part of the kth entry.

    Int Up [n_col+1] ;	Output argument.
    Int Ui [unz] ;	Output argument.
    double Ux [unz] ;	Output argument.
    double Uz [unz] ;	Output argument for complex versions.

	The min(n_row,n_col)-by-n_col matrix U is returned in compressed-column
	form.  The row indices of column j and corresponding numerical values
	are in

	    Ui [Up [j] ... Up [j+1]-1]
	    Ux [Up [j] ... Up [j+1]-1]	real part
	    Uz [Up [j] ... Up [j+1]-1]	imaginary part (complex versions)

	respectively.  Each column is stored in sorted order, from low row
	indices to higher.  The last entry in each column is the diagonal
	(assuming that it is nonzero).  The sizes of Up, Ui, Ux, and Uz are
	returned by umfpack_*_get_lunz.  If Up, Ui, or Ux (or Uz for the complex
	version) are not present, then the matrix U is not returned.  This is
	not an error condition.  The U matrix can be printed if n_col, Up, Ui,
	Ux (and Uz for the complex versions) are passed to
	umfpack_*_report_matrix (using the "column" form).

	Future complex version:  if Ux is present and Uz is NULL, then both real
	and imaginary parts will be returned in Ux[0..2*unz-1], with Ux[2*k]
	and Ux[2*k+1] being the real and imaginary part of the kth entry.

    Int P [n_row] ;		Output argument.

	The permutation vector P is defined as P [k] = i, where the original
	row i of A is the kth pivot row in PAQ.  If you do not want the P vector
	to be returned, simply pass (Int *) NULL for P.  This is not an error
	condition.  You can print P and Q with umfpack_*_report_perm.

    Int Q [n_col] ;		Output argument.

	The permutation vector Q is defined as Q [k] = j, where the original
	column j of A is the kth pivot column in PAQ.  If you not want the Q
	vector to be returned, simply pass (Int *) NULL for Q.  This is not
	an error condition.  Note that Q is not necessarily identical to
	Qtree, the column pre-ordering held in the Symbolic object.  Refer to
	the description of Qtree and Front_npivcol in umfpack_*_get_symbolic for
	details.

    double Dx [min(n_row,n_col)] ;	Output argument.
    double Dz [min(n_row,n_col)] ;	Output argument for complex versions.

	The diagonal of U is also returned in Dx and Dz.  You can extract the
	diagonal of U without getting all of U by passing a non-NULL Dx (and
	Dz for the complex version) and passing Up, Ui, and Ux as NULL.  Dx is
	the real part of the diagonal, and Dz is the imaginary part.

	Future complex version:  if Dx is present and Dz is NULL, then both real
	and imaginary parts will be returned in Dx[0..2*min(n_row,n_col)-1],
	with Dx[2*k] and Dx[2*k+1] being the real and imaginary part of the kth
	entry.

    Int *do_recip ;		Output argument.

	If do_recip is returned as zero (false), then the scale factors Rs [i]
	are to be used by multiplying row i by Rs [i].  Otherwise, the entries
	in row i are to be divided by Rs [i].

	If UMFPACK has been compiled with gcc, or for MATLAB as either a
	built-in routine or as a mexFunction, then the NRECIPROCAL flag is
	set, and do_recip will always be zero (false).

	NOTE: this argument is new to version 4.1.

    double Rs [n_row] ;		Output argument.

	The row scale factors are returned in Rs [0..n_row-1].  Row i of A is
	scaled by dividing or multiplying its values by Rs [i].  If default
	scaling is in use, Rs [i] is the sum of the absolute values of row i
	(or its reciprocal).  If max row scaling is in use, then Rs [i] is the
	maximum absolute value in row i (or its reciprocal).
	Otherwise, Rs [i] = 1.  If row i is all zero, Rs [i] = 1 as well.  For
	the complex version, an approximate absolute value is used
	(|x_real|+|x_imag|).

	NOTE: this argument is new to version 4.1.

    void *Numeric ;	Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.
*/
