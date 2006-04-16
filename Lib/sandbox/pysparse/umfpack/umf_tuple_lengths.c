/* ========================================================================== */
/* === UMF_tuple_lengths ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Determine the tuple list lengths, and the amount of memory required for */
/* them.  Return the amount of memory needed to store all the tuples. */
/* This routine assumes that the tuple lists themselves are either already */
/* deallocated, or will be shortly (so Row[ ].tlen and Col[ ].tlen are */
/* overwritten) */

#include "umf_internal.h"

GLOBAL Int UMF_tuple_lengths	    /* return memory usage */
(
    NumericType *Numeric,
    WorkType *Work,
    double *p_dusage		    /* output argument */
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int e, nrows, ncols, nel, i, *Rows, *Cols, row, col, n_row, n_col, *E,
	*Row_degree, *Row_tlen, *Col_degree, *Col_tlen, usage, n1 ;
    double dusage ;
    Element *ep ;
    Unit *p ;

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    E = Work->E ;
    Row_degree = Numeric->Rperm ;   /* for NON_PIVOTAL_ROW macro only */
    Col_degree = Numeric->Cperm ;   /* for NON_PIVOTAL_COL macro only */
    Row_tlen   = Numeric->Uilen ;
    Col_tlen   = Numeric->Lilen ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;
    n1 = Work->n1 ;
    nel = Work->nel ;

    DEBUG3 (("TUPLE_LENGTHS: n_row "ID" n_col "ID" nel "ID"\n",
	n_row, n_col, nel)) ;
    ASSERT (nel < Work->elen) ;

    /* tuple list lengths already initialized to zero */

    /* ---------------------------------------------------------------------- */
    /* scan each element: count tuple list lengths (include element 0) */
    /* ---------------------------------------------------------------------- */

    for (e = 1 ; e <= nel ; e++)	/* for all elements, in any order */
    {
	if (E [e])
	{
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	    p = Numeric->Memory + E [e] ;
	    GET_ELEMENT_PATTERN (ep, p, Cols, Rows, ncols) ;
	    nrows = ep->nrows ;
	    for (i = 0 ; i < nrows ; i++)
	    {
		row = Rows [i] ;
		ASSERT (row == EMPTY || (row >= n1 && row < n_row)) ;
		if (row >= n1)
		{
		    ASSERT (NON_PIVOTAL_ROW (row)) ;
		    Row_tlen [row] ++ ;
		}
	    }
	    for (i = 0 ; i < ncols ; i++)
	    {
		col = Cols [i] ;
		ASSERT (col == EMPTY || (col >= n1 && col < n_col)) ;
		if (col >= n1)
		{
		    ASSERT (NON_PIVOTAL_COL (col)) ;
		    Col_tlen [col] ++ ;
		}
	    }
	}
    }

    /* note: tuple lengths are now modified, but the tuple lists are not */
    /* updated to reflect that fact. */

    /* ---------------------------------------------------------------------- */
    /* determine the required memory to hold all the tuple lists */
    /* ---------------------------------------------------------------------- */

    DEBUG0 (("UMF_build_tuples_usage\n")) ;

    usage = 0 ;
    dusage = 0 ;

    ASSERT (Col_tlen && Col_degree) ;

    for (col = n1 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col))
	{
	    usage  += 1 +  UNITS (Tuple, TUPLES (Col_tlen [col])) ;
	    dusage += 1 + DUNITS (Tuple, TUPLES (Col_tlen [col])) ;
	    DEBUG0 ((" col: "ID" tlen "ID" usage so far: "ID"\n",
		     col, Col_tlen [col], usage)) ;
	}
    }

    ASSERT (Row_tlen && Row_degree) ;

    for (row = n1 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row))
	{
	    usage  += 1 +  UNITS (Tuple, TUPLES (Row_tlen [row])) ;
	    dusage += 1 + DUNITS (Tuple, TUPLES (Row_tlen [row])) ;
	    DEBUG0 ((" row: "ID" tlen "ID" usage so far: "ID"\n",
		     row, Row_tlen [row], usage)) ;
	}
    }

    DEBUG0 (("UMF_build_tuples_usage "ID" %g\n", usage, dusage)) ;

    *p_dusage = dusage ;
    return (usage) ;
}
