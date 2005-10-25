/* ========================================================================== */
/* === UMF_scale_column ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Scale the current pivot column, move the pivot row and column into place,
    and log the permutation.
*/

#include "umf_internal.h"
#include "umf_mem_free_tail_block.h"
#include "umf_scale.h"

/* ========================================================================== */
/* === shift_pivot_row ====================================================== */
/* ========================================================================== */

/* Except for the BLAS, most of the time is typically spent in the following
 * shift_pivot_row routine.  It copies the pivot row into the U block, and
 * then fills in the whole in the C block by shifting the last row of C into
 * the row vacated by the pivot row.
 */

PRIVATE void shift_pivot_row (Entry *Fd, Entry *Fs, Entry *Fe, Int len, Int d)
{
    Int j ;
#pragma ivdep
    for (j = 0 ; j < len ; j++)
    {
	Fd [j]   = Fs [j*d] ;
	Fs [j*d] = Fe [j*d] ;
    }
}

/* ========================================================================== */
/* === UMF_scale_column ===================================================== */
/* ========================================================================== */

GLOBAL void UMF_scale_column
(
    NumericType *Numeric,
    WorkType *Work
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int k, k1, fnr_curr, fnrows, fncols, *Frpos, *Fcpos, pivrow, pivcol,
	*Frows, *Fcols, fnc_curr, fnpiv, *Row_tuples, nb,
	*Col_tuples, *Rperm, *Cperm, fspos, col2, row2 ;
    Entry pivot_value, *Fcol,
	*Flublock, *Flblock, *Fublock, *Fcblock ;
#ifndef NDEBUG
    Int *Col_degree, *Row_degree ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    fnrows = Work->fnrows ;
    fncols = Work->fncols ;
    fnpiv = Work->fnpiv ;

    /* ---------------------------------------------------------------------- */

    Rperm = Numeric->Rperm ;
    Cperm = Numeric->Cperm ;

    /* ---------------------------------------------------------------------- */

    Flublock = Work->Flublock ;
    Flblock  = Work->Flblock ;
    Fublock  = Work->Fublock ;
    Fcblock  = Work->Fcblock ;

    fnr_curr = Work->fnr_curr ;
    fnc_curr = Work->fnc_curr ;
    Frpos = Work->Frpos ;
    Fcpos = Work->Fcpos ;
    Frows = Work->Frows ;
    Fcols = Work->Fcols ;
    pivrow = Work->pivrow ;
    pivcol = Work->pivcol ;

    ASSERT (pivrow >= 0 && pivrow < Work->n_row) ;
    ASSERT (pivcol >= 0 && pivcol < Work->n_col) ;

#ifndef NDEBUG
    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
#endif

    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;
    nb = Work->nb ;

#ifndef NDEBUG
    ASSERT (fnrows == Work->fnrows_new + 1) ;
    ASSERT (fncols == Work->fncols_new + 1) ;
    DEBUG1 (("SCALE COL: fnrows "ID" fncols "ID"\n", fnrows, fncols)) ;
    DEBUG2 (("\nFrontal matrix, including all space:\n"
		"fnr_curr "ID" fnc_curr "ID" nb    "ID"\n"
		"fnrows   "ID" fncols   "ID" fnpiv "ID"\n",
		fnr_curr, fnc_curr, nb, fnrows, fncols, fnpiv)) ;
    DEBUG2 (("\nJust the active part:\n")) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Fcblock,  fnr_curr, fnrows, fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Flblock,  fnr_curr, fnrows, fnpiv);
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Fublock,  fnc_curr, fncols, fnpiv) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Flublock, nb, fnpiv, fnpiv) ;
#endif

    /* ====================================================================== */
    /* === Shift pivot row and column ======================================= */
    /* ====================================================================== */

    /* ---------------------------------------------------------------------- */
    /* move pivot column into place */
    /* ---------------------------------------------------------------------- */

    /* Note that the pivot column is already in place.  Just shift the last
     * column into the position vacated by the pivot column. */

    fspos = Fcpos [pivcol] ;

    /* one less column in the contribution block */
    fncols = --(Work->fncols) ;

    if (fspos != fncols * fnr_curr)
    {

	Int fs = fspos / fnr_curr ;

	DEBUG6 (("Shift pivot column in front\n")) ;
	DEBUG6 (("fspos: "ID" flpos: "ID"\n", fspos, fncols * fnr_curr)) ;

	/* ------------------------------------------------------------------ */
	/* move Fe => Fs */
	/* ------------------------------------------------------------------ */

	/* column of the contribution block: */
	{
	    /* Fs: current position of pivot column in contribution block */
	    /* Fe: position of last column in contribution block */
	    Int i ;
	    Entry *Fs, *Fe ;
	    Fs = Fcblock + fspos ;
	    Fe = Fcblock + fncols * fnr_curr ;
#pragma ivdep
	    for (i = 0 ; i < fnrows ; i++)
	    {
		Fs [i] = Fe [i] ;
	    }
	}

	/* column of the U2 block */
	{
	    /* Fs: current position of pivot column in U block */
	    /* Fe: last column in U block */
	    Int i ;
	    Entry *Fs, *Fe ;
	    Fs = Fublock + fs ;
	    Fe = Fublock + fncols ;
#pragma ivdep
	    for (i = 0 ; i < fnpiv ; i++)
	    {
		Fs [i * fnc_curr] = Fe [i * fnc_curr] ;
	    }
	}

	/* move column Fe to Fs in the Fcols pattern */
	col2 = Fcols [fncols] ;
	Fcols [fs] = col2 ;
	Fcpos [col2] = fspos ;
    }

    /* pivot column is no longer in the frontal matrix */
    Fcpos [pivcol] = EMPTY ;

#ifndef NDEBUG
    DEBUG2 (("\nFrontal matrix after col swap, including all space:\n"
		"fnr_curr "ID" fnc_curr "ID" nb    "ID"\n"
		"fnrows   "ID" fncols   "ID" fnpiv "ID"\n",
		fnr_curr, fnc_curr, nb,
		fnrows, fncols, fnpiv)) ;
    DEBUG2 (("\nJust the active part:\n")) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Fcblock,  fnr_curr, fnrows, fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Flblock,  fnr_curr, fnrows, fnpiv+1);
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Fublock,  fnc_curr, fncols, fnpiv) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Flublock, nb, fnpiv, fnpiv+1) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* move pivot row into place */
    /* ---------------------------------------------------------------------- */

    fspos = Frpos [pivrow] ;

    /* one less row in the contribution block */
    fnrows = --(Work->fnrows) ;

    DEBUG6 (("Swap/shift pivot row in front:\n")) ;
    DEBUG6 (("fspos: "ID" flpos: "ID"\n", fspos, fnrows)) ;

    if (fspos == fnrows)
    {

	/* ------------------------------------------------------------------ */
	/* move Fs => Fd */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("row case 1\n")) ;

	/* row of the contribution block: */
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Fublock + fnpiv * fnc_curr ;
	    Fs = Fcblock + fspos ;
#pragma ivdep
	    for (j = 0 ; j < fncols ; j++)
	    {
		Fd [j] = Fs [j * fnr_curr] ;
	    }
	}

	/* row of the L2 block: */
	if (Work->pivrow_in_front)
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
#pragma ivdep
	    for (j = 0 ; j <= fnpiv ; j++)
	    {
		Fd [j * nb] = Fs [j * fnr_curr] ;
	    }
	}
	else
	{
	    Int j ;
	    Entry *Fd, *Fs ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
#pragma ivdep
	    for (j = 0 ; j < fnpiv ; j++)
	    {
		ASSERT (IS_ZERO (Fs [j * fnr_curr])) ;
		CLEAR (Fd [j * nb]) ;
	    }
	    Fd [fnpiv * nb] = Fs [fnpiv * fnr_curr] ;
	}
    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* move Fs => Fd */
	/* move Fe => Fs */
	/* ------------------------------------------------------------------ */

	DEBUG6 (("row case 2\n")) ;
	/* this is the most common case, by far */

	/* row of the contribution block: */
	{
	    /* Fd: destination of pivot row on U block */
	    /* Fs: current position of pivot row in contribution block */
	    /* Fe: position of last row in contribution block */
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Fublock + fnpiv * fnc_curr ;
	    Fs = Fcblock + fspos ;
	    Fe = Fcblock + fnrows ;
	    shift_pivot_row (Fd, Fs, Fe, fncols, fnr_curr) ;
	}

	/* row of the L2 block: */
	if (Work->pivrow_in_front)
	{
	    /* Fd: destination of pivot row in LU block */
	    /* Fs: current position of pivot row in L block */
	    /* Fe: last row in L block */
	    Int j ;
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
	    Fe = Flblock  + fnrows ;
#pragma ivdep
	    for (j = 0 ; j <= fnpiv ; j++)
	    {
		Fd [j * nb]       = Fs [j * fnr_curr] ;
		Fs [j * fnr_curr] = Fe [j * fnr_curr] ;
	    }
	}
	else
	{
	    Int j ;
	    Entry *Fd, *Fs, *Fe ;
	    Fd = Flublock + fnpiv ;
	    Fs = Flblock  + fspos ;
	    Fe = Flblock  + fnrows ;
#pragma ivdep
	    for (j = 0 ; j < fnpiv ; j++)
	    {
		ASSERT (IS_ZERO (Fs [j * fnr_curr])) ;
		CLEAR (Fd [j * nb]) ;
		Fs [j * fnr_curr] = Fe [j * fnr_curr] ;
	    }
	    Fd [fnpiv * nb]       = Fs [fnpiv * fnr_curr] ;
	    Fs [fnpiv * fnr_curr] = Fe [fnpiv * fnr_curr] ;
	}

	/* move row Fe to Fs in the Frows pattern */
	row2 = Frows [fnrows] ;
	Frows [fspos] = row2 ;
	Frpos [row2] = fspos ;

    }
    /* pivot row is no longer in the frontal matrix */
    Frpos [pivrow] = EMPTY ;

#ifndef NDEBUG
    DEBUG2 (("\nFrontal matrix after row swap, including all space:\n"
		"fnr_curr "ID" fnc_curr "ID" nb    "ID"\n"
		"fnrows   "ID" fncols   "ID" fnpiv "ID"\n",
		Work->fnr_curr, Work->fnc_curr, Work->nb,
		Work->fnrows, Work->fncols, Work->fnpiv)) ;
    DEBUG2 (("\nJust the active part:\n")) ;
    DEBUG7 (("C  block: ")) ;
    UMF_dump_dense (Fcblock,  fnr_curr, fnrows, fncols) ;
    DEBUG7 (("L  block: ")) ;
    UMF_dump_dense (Flblock,  fnr_curr, fnrows, fnpiv+1);
    DEBUG7 (("U' block: ")) ;
    UMF_dump_dense (Fublock,  fnc_curr, fncols, fnpiv+1) ;
    DEBUG7 (("LU block: ")) ;
    UMF_dump_dense (Flublock, nb, fnpiv+1, fnpiv+1) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* Frpos [row] >= 0 for each row in pivot column pattern.   */
    /* offset into pattern is given by:				*/
    /* Frpos [row] == offset - 1				*/
    /* Frpos [pivrow] is EMPTY */

    /* Fcpos [col] >= 0 for each col in pivot row pattern.	*/
    /* Fcpos [col] == (offset - 1) * fnr_curr			*/
    /* Fcpos [pivcol] is EMPTY */

    /* Fcols [0..fncols-1] is the pivot row pattern (excl pivot cols) */
    /* Frows [0..fnrows-1] is the pivot col pattern (excl pivot rows) */

    /* ====================================================================== */
    /* === scale pivot column =============================================== */
    /* ====================================================================== */

    /* pivot column (except for pivot entry itself) */
    Fcol = Flblock + fnpiv * fnr_curr ;
    /* fnpiv-th pivot in frontal matrix located in Flublock (fnpiv, fnpiv) */
    pivot_value = Flublock [fnpiv + fnpiv * nb] ;

    /* this is the kth global pivot */
    k = Work->npiv + fnpiv ;

    DEBUG4 (("Pivot value: ")) ;
    EDEBUG4 (pivot_value) ;
    DEBUG4 (("\n")) ;

    UMF_scale (fnrows, pivot_value, Fcol) ;

    /* ---------------------------------------------------------------------- */
    /* deallocate the pivot row and pivot column tuples */
    /* ---------------------------------------------------------------------- */

    UMF_mem_free_tail_block (Numeric, Row_tuples [pivrow]) ;
    UMF_mem_free_tail_block (Numeric, Col_tuples [pivcol]) ;

    Row_tuples [pivrow] = 0 ;
    Col_tuples [pivcol] = 0 ;

    DEBUG5 (("number of pivots prior to this one: "ID"\n", k)) ;
    ASSERT (NON_PIVOTAL_ROW (pivrow)) ;
    ASSERT (NON_PIVOTAL_COL (pivcol)) ;

    /* save row and column inverse permutation */
    k1 = ONES_COMPLEMENT (k) ;
    Rperm [pivrow] = k1 ;			/* aliased with Row_degree */
    Cperm [pivcol] = k1 ;			/* aliased with Col_degree */

    ASSERT (!NON_PIVOTAL_ROW (pivrow)) ;
    ASSERT (!NON_PIVOTAL_COL (pivcol)) ;

    /* ---------------------------------------------------------------------- */
    /* Keep track of the pivot order.  This is the kth pivot row and column. */
    /* ---------------------------------------------------------------------- */

    /* keep track of pivot rows and columns in the LU, L, and U blocks */
    ASSERT (fnpiv < MAXNB) ;
    Work->Pivrow [fnpiv] = pivrow ;
    Work->Pivcol [fnpiv] = pivcol ;

    /* ====================================================================== */
    /* === one step in the factorization is done ============================ */
    /* ====================================================================== */

    /* One more step is done, except for pending updates to the U and C blocks
     * of this frontal matrix.  Those are saved up, and applied by
     * UMF_blas3_update when enough pivots have accumulated.   Also, the
     * LU factors for these pending pivots have not yet been stored. */

    Work->fnpiv++ ;

#ifndef NDEBUG
    DEBUG7 (("Current frontal matrix: (after pivcol scale)\n")) ;
    UMF_dump_current_front (Numeric, Work, TRUE) ;
#endif

}
