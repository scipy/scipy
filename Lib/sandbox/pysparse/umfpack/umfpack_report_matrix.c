/* ========================================================================== */
/* === UMFPACK_report_matrix ================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints a column or row-oriented matrix.  See
    umfpack_report_matrix.h for details.
*/

#include "umf_internal.h"
#include "umf_malloc.h"
#include "umf_free.h"

GLOBAL Int UMFPACK_report_matrix
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    Int col_form,		/* 1: column form, 0: row form */
    const double Control [UMFPACK_CONTROL]
)
{
    Int prl, i, k, length, ilast, p, nz, prl1, p1, p2, n, n_i, do_values ;
    char *vector, *index ;
    Entry a ;

    /* ---------------------------------------------------------------------- */
    /* determine the form, and check if inputs exist */
    /* ---------------------------------------------------------------------- */

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    if (col_form)
    {
	vector = "column" ;	/* column vectors */
	index = "row" ;		/* with row indices */
	n = n_col ;
	n_i = n_row ;
    }
    else
    {
	vector = "row" ;	/* row vectors */
	index = "column" ;	/* with column indices */
	n = n_row ;
	n_i = n_col ;
    }

    PRINTF (("%s-form matrix, n_row "ID" n_col "ID", ", vector, n_row, n_col)) ;

    if (n_row <= 0 || n_col <= 0)
    {
	PRINTF (("ERROR: n_row <= 0 or n_col <= 0\n\n")) ;
	return (UMFPACK_ERROR_n_nonpositive) ;
    }

    if (!Ap)
    {
	PRINTF (("ERROR: Ap missing\n\n")) ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    nz = Ap [n] ;
    PRINTF (("nz = "ID". ", nz)) ;
    if (nz < 0)
    {
	PRINTF (("ERROR: number of entries < 0\n\n")) ;
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    if (Ap [0] != 0)
    {
	PRINTF (("ERROR: Ap ["ID"] = "ID" must be "ID"\n\n",
	    (Int) INDEX (0), INDEX (Ap [0]), (Int) INDEX (0))) ;
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    if (!Ai)
    {
	PRINTF (("ERROR: Ai missing\n\n")) ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

#ifdef COMPLEX
    do_values = Ax && Az ;
#else
    do_values = Ax != (double *) NULL ;
#endif

    PRINTF4 (("\n")) ;

    /* ---------------------------------------------------------------------- */
    /* check the row/column pointers, Ap */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < n ; k++)
    {
	if (Ap [k] < 0)
	{
	    PRINTF (("ERROR: Ap ["ID"] < 0\n\n", INDEX (k))) ;
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
	if (Ap [k] > nz)
	{
	    PRINTF (("ERROR: Ap ["ID"] > size of Ai\n\n", INDEX (k))) ;
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
    }

    for (k = 0 ; k < n ; k++)
    {
	length = Ap [k+1] - Ap [k] ;
	if (length < 0)
	{
	    PRINTF (("ERROR: # entries in %s "ID" is < 0\n\n",
		vector, INDEX (k))) ;
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* print each vector */
    /* ---------------------------------------------------------------------- */

    prl1 = prl ;

    for (k = 0 ; k < n ; k++)
    {
	/* if prl is 4, print the first 10 entries of the first 10 vectors */
	if (k < 10)
	{
	    prl = prl1 ;
	}
	/* get the vector pointers */
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;
	length = p2 - p1 ;
	PRINTF4 (("\n    %s "ID": start: "ID" end: "ID" entries: "ID"\n",
	    vector, INDEX (k), p1, p2-1, length)) ;
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    PRINTF4 (("\t%s "ID" ", index, INDEX (i))) ;
	    if (do_values && prl >= 4)
	    {
		PRINTF ((":")) ;
		ASSIGN (a, Ax [p], Az [p]) ;
		PRINT_ENTRY (a) ;
	    }
	    if (i < 0 || i >= n_i)
	    {
		PRINTF ((" ERROR: %s index "ID" out of range in %s "ID"\n\n",
		    index, INDEX (i), vector, INDEX (k))) ;
		return (UMFPACK_ERROR_invalid_matrix) ;
	    }
	    if (i <= ilast)
	    {
		PRINTF ((" ERROR: %s index "ID" out of order (or duplicate) in "
		    "%s "ID"\n\n", index, INDEX (i), vector, INDEX (k))) ;
		return (UMFPACK_ERROR_invalid_matrix) ;
	    }
	    PRINTF4 (("\n")) ;
	    /* truncate printout, but continue to check matrix */
	    if (prl == 4 && (p - p1) == 9 && length > 10)
	    {
		PRINTF4 (("\t...\n")) ;
		prl-- ;
	    }
	    ilast = i ;
	}
	/* truncate printout, but continue to check matrix */
	if (prl == 4 && k == 9 && n > 10)
	{
	    PRINTF4 (("\n    ...\n")) ;
	    prl-- ;
	}
    }
    prl = prl1 ;

    /* ---------------------------------------------------------------------- */
    /* return the status of the matrix */
    /* ---------------------------------------------------------------------- */

    PRINTF4 (("    %s-form matrix ", vector)) ;
    PRINTF (("OK\n\n")) ;

    return (UMFPACK_OK) ;
}
