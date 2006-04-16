/* ========================================================================== */
/* === UMFPACK_report_symbolic ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Prints the Symbolic object. See umfpack_report_symbolic.h
    for details.  Does not print new Cdeg, Rdeg, Esize, and the Diagonal_map.

    Dynamic memory usage:  Allocates a size MAX (n_row,n_col)*sizeof(Int)
    workspace via a single call to UMF_malloc and then frees all of it via
    UMF_free on return.  The workspace is not allocated if an early error
    return occurs before the workspace is needed.
*/

#include "umf_internal.h"
#include "umf_valid_symbolic.h"
#include "umf_report_perm.h"
#include "umf_malloc.h"
#include "umf_free.h"

GLOBAL Int UMFPACK_report_symbolic
(
    void *SymbolicHandle,
    const double Control [UMFPACK_CONTROL]
)
{
    Int n_row, n_col, nz, nchains, nfr, maxnrows, maxncols, prl,
	k, chain, frontid, frontid1, frontid2, kk, *Chain_start, *W,
	*Chain_maxrows, *Chain_maxcols, *Front_npivcol, *Front_1strow,
	*Front_leftmostdesc, *Front_parent, done, status1, status2 ;
    SymbolicType *Symbolic ;

    prl = GET_CONTROL (UMFPACK_PRL, UMFPACK_DEFAULT_PRL) ;

    if (prl <= 2)
    {
	return (UMFPACK_OK) ;
    }

    PRINTF (("Symbolic object: ")) ;

    Symbolic = (SymbolicType *) SymbolicHandle ;
    if (!UMF_valid_symbolic (Symbolic))
    {
	PRINTF (("ERROR: invalid\n")) ;
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    n_row = Symbolic->n_row ;
    n_col = Symbolic->n_col ;

    nz = Symbolic->nz ;

    nchains = Symbolic->nchains ;
    nfr = Symbolic->nfr ;
    maxnrows = Symbolic->maxnrows ;
    maxncols = Symbolic->maxncols ;

    Chain_start = Symbolic->Chain_start ;
    Chain_maxrows = Symbolic->Chain_maxrows ;
    Chain_maxcols = Symbolic->Chain_maxcols ;
    Front_npivcol = Symbolic->Front_npivcol ;
    Front_1strow = Symbolic->Front_1strow ;
    Front_leftmostdesc = Symbolic->Front_leftmostdesc ;
    Front_parent = Symbolic->Front_parent ;

    if (prl >= 4)
    {

	PRINTF (("\n    matrix to be factorized:\n")) ;
	PRINTF (("\tn_row: "ID" n_col: "ID"\n", n_row, n_col)) ;
	PRINTF (("\tnumber of entries: "ID"\n", nz)) ;
	PRINTF (("    block size used for dense matrix kernels:   "ID"\n",
	Symbolic->nb)) ;

	PRINTF (("    strategy used:                              ")) ;
	/* strategy cannot be auto */
	if (Symbolic->strategy == UMFPACK_STRATEGY_SYMMETRIC)
	{
	    PRINTF (("symmetric")) ;
	}
	else if (Symbolic->strategy == UMFPACK_STRATEGY_UNSYMMETRIC)
	{
	    PRINTF (("unsymmetric")) ;
	}
	else if (Symbolic->strategy == UMFPACK_STRATEGY_2BY2)
	{
	    PRINTF (("symmetric 2-by-2")) ;
	}
	PRINTF (("\n")) ;

	PRINTF (("    ordering used:                              ")) ;
	if (Symbolic->ordering == UMFPACK_ORDERING_COLAMD)
	{
	    PRINTF (("colamd on A\n")) ;
	}
	else if (Symbolic->ordering == UMFPACK_ORDERING_AMD)
	{
	    PRINTF (("amd on A+A'\n")) ;
	}
	else if (Symbolic->ordering == UMFPACK_ORDERING_GIVEN)
	{
	    PRINTF (("provided by user")) ;
	}
	PRINTF (("\n")) ;

	PRINTF (("    performn column etree postorder:            ")) ;
	if (Symbolic->fixQ)
	{
	    PRINTF (("no\n")) ;
	}
	else
	{
	    PRINTF (("yes\n")) ;
	}

	PRINTF (("    prefer diagonal pivoting (attempt P=Q):     ")) ;
	if (Symbolic->prefer_diagonal)
	{
	    PRINTF (("yes\n")) ;
	}
	else
	{
	    PRINTF (("no\n")) ;
	}

	PRINTF (("    variable-size part of Numeric object:\n")) ;
	PRINTF (("\tminimum initial size (Units): %.20g  (MBytes): %.1f\n",
	    Symbolic->dnum_mem_init_usage,
	    MBYTES (Symbolic->dnum_mem_init_usage))) ;
	PRINTF (("\testimated peak size (Units):  %.20g  (MBytes): %.1f\n",
	    Symbolic->num_mem_usage_est,
	    MBYTES (Symbolic->num_mem_usage_est))) ;
	PRINTF (("\testimated final size (Units): %.20g  (MBytes): %.1f\n",
	    Symbolic->num_mem_size_est,
	    MBYTES (Symbolic->num_mem_size_est))) ;
	PRINTF (("    symbolic factorization memory usage (Units):"
	    " %.20g  (MBytes): %.1f\n",
	    Symbolic->peak_sym_usage,
	    MBYTES (Symbolic->peak_sym_usage))) ;
	PRINTF (("    frontal matrices / supercolumns:\n")) ;
	PRINTF (("\tnumber of frontal chains: "ID"\n", nchains)) ;
	PRINTF (("\tnumber of frontal matrices: "ID"\n", nfr)) ;
	PRINTF (("\tlargest frontal matrix row dimension: "ID"\n", maxnrows)) ;
	PRINTF (("\tlargest frontal matrix column dimension: "ID"\n",maxncols));
    }

    k = 0 ;
    done = FALSE ;

    for (chain = 0 ; chain < nchains ; chain++)
    {
	frontid1 = Chain_start [chain] ;
	frontid2 = Chain_start [chain+1] - 1 ;
	PRINTF4 (("\n    Frontal chain: "ID".  Frontal matrices "ID" to "ID"\n",
	    INDEX (chain), INDEX (frontid1), INDEX (frontid2))) ;
	PRINTF4 (("\tLargest frontal matrix in Frontal chain: "ID"-by-"ID"\n",
	    Chain_maxrows [chain], Chain_maxcols [chain])) ;
	for (frontid = frontid1 ; frontid <= frontid2 ; frontid++)
	{
	    kk = Front_npivcol [frontid] ;
	    PRINTF4 (("\tFront: "ID"  pivot cols: "ID" (pivot columns "ID" to "
		ID")\n", INDEX (frontid), kk, INDEX (k), INDEX (k+kk-1))) ;
	    PRINTF4 (("\t    pivot row candidates: "ID" to "ID"\n",
		INDEX (Front_1strow [Front_leftmostdesc [frontid]]),
		INDEX (Front_1strow [frontid+1]-1))) ;
	    PRINTF4 (("\t    leftmost descendant: "ID"\n",
		INDEX (Front_leftmostdesc [frontid]))) ;
	    PRINTF4 (("\t    1st new candidate row : "ID"\n",
		INDEX (Front_1strow [frontid]))) ;
	    PRINTF4 (("\t    parent:")) ;
	    if (Front_parent [frontid] == EMPTY)
	    {
		PRINTF4 ((" (none)\n")) ;
	    }
	    else
	    {
		PRINTF4 ((" "ID"\n", INDEX (Front_parent [frontid]))) ;
	    }
	    done = (frontid == 20 && frontid < nfr-1 && prl == 4) ;
	    if (done)
	    {
		PRINTF4 (("\t...\n")) ;
		break ;
	    }
	    k += kk ;
	}
	if (Front_npivcol [nfr] != 0)
	{
	    PRINTF4 (("\tFront: "ID" placeholder for "ID" empty columns\n",
		INDEX (nfr), Front_npivcol [nfr])) ;
	}
	if (done)
	{
	    break ;
	}
    }

    W = (Int *) UMF_malloc (MAX (n_row, n_col), sizeof (Int)) ;
    if (!W)
    {
	PRINTF (("ERROR: out of memory to check Symbolic object\n\n")) ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    PRINTF4 (("\nInitial column permutation, Q1: ")) ;
    status1 = UMF_report_perm (n_col, Symbolic->Cperm_init, W, prl, 0) ;

    PRINTF4 (("\nInitial row permutation, P1: ")) ;
    status2 = UMF_report_perm (n_row, Symbolic->Rperm_init, W, prl, 0) ;

    (void) UMF_free ((void *) W) ;

    if (status1 != UMFPACK_OK || status2 != UMFPACK_OK)
    {
	return (UMFPACK_ERROR_invalid_Symbolic_object) ;
    }

    PRINTF4 (("    Symbolic object:  ")) ;
    PRINTF (("OK\n\n")) ;
    return (UMFPACK_OK) ;
}
