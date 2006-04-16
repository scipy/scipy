/* ========================================================================== */
/* === AMD_2 ================================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* AMD_2:  performs the AMD ordering on a symmetric sparse matrix A, followed
 * by a postordering (via depth-first search) of the assembly tree using the
 * AMD_postorder routine.
 */

#include "amd_internal.h"

GLOBAL void AMD_2
(
    Int n,		/* A is n-by-n, where n > 0 */
    Int Pe [ ],		/* Pe [0..n-1]: index in Iw of row i on input */
    Int Iw [ ],		/* workspace of size iwlen. Iw [0..pfree-1]
			 * holds the matrix on input */
    Int Len [ ],	/* Len [0..n-1]: length for row/column i on input */
    Int iwlen,		/* length of Iw. iwlen >= pfree + n */
    Int pfree,		/* Iw [pfree ... iwlen-1] is empty on input */

    /* 7 size-n workspaces, not defined on input: */
    Int Nv [ ],		/* the size of each supernode on output */
    Int Next [ ],	/* the output inverse permutation */
    Int Last [ ],	/* the output permutation */
    Int Head [ ],
    Int Elen [ ],	/* the size columns of L for each supernode */
    Int Degree [ ],
    Int W [ ],

    /* control parameters and output statistics */
    double Control [ ],	/* array of size AMD_CONTROL */
    double Info [ ]	/* array of size AMD_INFO */
)
{

/*
 * Given a representation of the nonzero pattern of a symmetric matrix, A, 
 * (excluding the diagonal) perform an approximate minimum (UMFPACK/MA38-style)
 * degree ordering to compute a pivot order such that the introduction of
 * nonzeros (fill-in) in the Cholesky factors A = LL' is kept low.  At each
 * step, the pivot selected is the one with the minimum UMFAPACK/MA38-style
 * upper-bound on the external degree.  This routine can optionally perform
 * aggresive absorption (as done by MC47B in the Harwell Subroutine
 * Library).
 * 
 * The approximate degree algorithm implemented here is the symmetric analog of
 * the degree update algorithm in MA38 and UMFPACK (the Unsymmetric-pattern
 * MultiFrontal PACKage, both by Davis and Duff).  The routine is based on the
 * MA27 minimum degree ordering algorithm by Iain Duff and John Reid.
 * 
 * This routine is a translation of the original AMDBAR and MC47B routines,
 * in Fortran, with the following modifications:
 * 
 * (1) dense rows/columns are removed prior to ordering the matrix, and placed
 *	last in the output order.  The presence of a dense row/column can
 *	increase the ordering time by up to O(n^2), unless they are removed
 *	prior to ordering.
 *
 * (2) the minimum degree ordering is followed by a postordering (depth-first
 *	search) of the assembly tree.  Note that mass elimination (discussed
 *	below) combined with the approximate degree update can lead to the mass
 *	elimination of nodes with lower exact degree than the current pivot
 *	element.  No additional fill-in is caused in the representation of the
 *	Schur complement.  The mass-eliminated nodes merge with the current
 *	pivot element.  They are ordered prior to the current pivot element.
 *	Because they can have lower exact degree than the current element, the
 *	merger of two or more of these nodes in the current pivot element can
 *	lead to a single element that is not a "fundamental supernode".  The
 *	diagonal block can have zeros in it.  Thus, the assembly tree used here
 *	is not guaranteed to be the precise supernodal elemination tree (with
 *	"funadmental" supernodes), and the postordering performed by this
 *	routine is not guaranteed to be a precise postordering of the
 *	elimination tree.
 *
 * (3) input parameters are added, to control aggressive absorption and the
 *	detection of "dense" rows/columns of A.
 *
 * (4) additional statistical information is returned, such as the number of
 *	nonzeros in L, and the flop counts for subsequent LDL' and LU
 *	factorizations.  These are slight upper bounds, because of the mass
 *	elimination issue discussed above.
 *
 * (5) additional routines are added to interface this routine to MATLAB
 *	to provide a simple C-callable user-interface, to check inputs for
 *	errors, compute the symmetry of the pattern of A and the number of
 *	nonzeros in each row/column of A+A', to compute the pattern of A+A',
 *	to perform the assembly tree postordering, and to provide debugging
 *	ouput.  Many of these functions are also provided by the Fortran
 *	Harwell Subroutine Library routine MC47A.
 *
 * (6) both "int" and "long" versions are provided.  In the descriptions below
 *	and integer is and "int" or "long", depending on which version is
 *	being used.

 **********************************************************************
 ***** CAUTION:  ARGUMENTS ARE NOT CHECKED FOR ERRORS ON INPUT.  ******
 **********************************************************************
 ** If you want error checking, a more versatile input format, and a **
 ** simpler user interface, use amd_order or amd_l_order instead.    **
 ** This routine is not meant to be user-callable.                   **
 **********************************************************************

 * -----------------------------------------------------------------------------
 * References:
 * -----------------------------------------------------------------------------
 *
 *  [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern multifrontal
 *	method for sparse LU factorization", SIAM J. Matrix Analysis and
 *	Applications, vol. 18, no. 1, pp. 140-158.  Discusses UMFPACK / MA38,
 *	which first introduced the approximate minimum degree used by this
 *	routine.
 *
 *  [2] Patrick Amestoy, Timothy A. Davis, and Iain S. Duff, "An approximate
 *	minimum degree ordering algorithm," SIAM J. Matrix Analysis and
 *	Applications, vol. 17, no. 4, pp. 886-905, 1996.  Discusses AMDBAR and
 *	MC47B, which are the Fortran versions of this routine.
 *
 *  [3] Alan George and Joseph Liu, "The evolution of the minimum degree
 *	ordering algorithm," SIAM Review, vol. 31, no. 1, pp. 1-19, 1989.
 *	We list below the features mentioned in that paper that this code
 *	includes:
 *
 *	mass elimination:
 *	    Yes.  MA27 relied on supervariable detection for mass elimination.
 *
 *	indistinguishable nodes:
 *	    Yes (we call these "supervariables").  This was also in the MA27
 *	    code - although we modified the method of detecting them (the
 *	    previous hash was the true degree, which we no longer keep track
 *	    of).  A supervariable is a set of rows with identical nonzero
 *	    pattern.  All variables in a supervariable are eliminated together.
 *	    Each supervariable has as its numerical name that of one of its
 *	    variables (its principal variable).
 *
 *	quotient graph representation:
 *	    Yes.  We use the term "element" for the cliques formed during
 *	    elimination.  This was also in the MA27 code.  The algorithm can
 *	    operate in place, but it will work more efficiently if given some
 *	    "elbow room."
 *
 *	element absorption:
 *	    Yes.  This was also in the MA27 code.
 *
 *	external degree:
 *	    Yes.  The MA27 code was based on the true degree.
 *
 *	incomplete degree update and multiple elimination:
 *	    No.  This was not in MA27, either.  Our method of degree update
 *	    within MC47B is element-based, not variable-based.  It is thus
 *	    not well-suited for use with incomplete degree update or multiple
 *	    elimination.
 *
 * Authors, and Copyright (C) 2003 by:
 * Timothy A. Davis, Patrick Amestoy, Iain S. Duff, John K. Reid.
 *
 * Acknowledgements: This work (and the UMFPACK package) was supported by the
 * National Science Foundation (ASC-9111263, DMS-9223088, and DMS-0203270).
 * The UMFPACK/MA38 approximate degree update algorithm, the unsymmetric analog
 * which forms the basis of AMD, was developed while Tim Davis was supported by
 * CERFACS (Toulouse, France) in a post-doctoral position.  This C version, and
 * the etree postorder, were written while Tim Davis was on sabbatical at
 * Stanford University and Lawrence Berkeley National Laboratory.

 * -----------------------------------------------------------------------------
 * INPUT ARGUMENTS (unaltered):
 * -----------------------------------------------------------------------------
 
 * n:  The matrix order.  Restriction:  n >= 1.
 *
 * iwlen:  The size of the Iw array.  On input, the matrix is stored in
 *	Iw [0..pfree-1].  However, Iw [0..iwlen-1] should be slightly larger
 *	than what is required to hold the matrix, at least iwlen >= pfree + n.
 *	Otherwise, excessive compressions will take place.  The recommended
 *	value of iwlen is 1.2 * pfree + n, which is the value used in the
 *	user-callable interface to this routine (amd_order.c).  The algorithm
 *	will not run at all if iwlen < pfree.  Restriction: iwlen >= pfree + n.
 *	Note that this is slightly more restrictive than the actual minimum
 *	(iwlen >= pfree), but AMD_2 will be very slow with no elbow room.
 *	Thus, this routine enforces a bare minimum elbow room of size n.
 *
 * pfree: On input the tail end of the array, Iw [pfree..iwlen-1], is empty,
 *	and the matrix is stored in Iw [0..pfree-1].  During execution,
 *	additional data is placed in Iw, and pfree is modified so that
 *	Iw [pfree..iwlen-1] is always the unused part of Iw.
 *
 * Control:  A double array of size AMD_CONTROL containing input parameters that
 *	affect how the ordering is computed.  If NULL, then default settings
 *	are used.
 *
 *	Control [AMD_DENSE] is used to determine whether or not a given input
 *	row is "dense".  A row is "dense" if the number of entries in the row
 *	exceeds Control [AMD_DENSE] times sqrt (n), except that rows with 16 or
 *	fewer entries are never considered "dense".  To turn off the detection
 *	of dense rows, set Control [AMD_DENSE] to a negative number, or to a
 *	number larger than sqrt (n).  The default value of Control [AMD_DENSE]
 *	is AMD_DEFAULT_DENSE, which is defined in amd.h as 10.
 *
 *	Control [AMD_AGGRESSIVE] is used to determine whether or not aggressive
 *	absorption is to be performed.  If nonzero, then aggressive absorption
 *	is performed (this is the default).

 * -----------------------------------------------------------------------------
 * INPUT/OUPUT ARGUMENTS:
 * -----------------------------------------------------------------------------
 *
 * Pe:  An integer array of size n.  On input, Pe [i] is the index in Iw of
 *	the start of row i.  Pe [i] is ignored if row i has no off-diagonal
 *	entries.  Thus Pe [i] must be in the range 0 to pfree-1 for non-empty
 *	rows.
 *
 *	During execution, it is used for both supervariables and elements:
 *
 *	Principal supervariable i:  index into Iw of the description of
 *	    supervariable i.  A supervariable represents one or more rows of
 *	    the matrix with identical nonzero pattern.  In this case,
 *	    Pe [i] >= 0.
 *
 *	Non-principal supervariable i:  if i has been absorbed into another
 *	    supervariable j, then Pe [i] = FLIP (j), where FLIP (j) is defined
 *	    as (-(j)-2).  Row j has the same pattern as row i.  Note that j
 *	    might later be absorbed into another supervariable j2, in which
 *	    case Pe [i] is still FLIP (j), and Pe [j] = FLIP (j2) which is
 *	    < EMPTY, where EMPTY is defined as (-1) in amd_internal.h.
 *
 *	Unabsorbed element e:  the index into Iw of the description of element
 *	    e, if e has not yet been absorbed by a subsequent element.  Element
 *	    e is created when the supervariable of the same name is selected as
 *	    the pivot.  In this case, Pe [i] >= 0.
 *
 *	Absorbed element e:  if element e is absorbed into element e2, then
 *	    Pe [e] = FLIP (e2).  This occurs when the pattern of e (which we
 *	    refer to as Le) is found to be a subset of the pattern of e2 (that
 *	    is, Le2).  In this case, Pe [i] < EMPTY.  If element e is "null"
 *	    (it has no nonzeros outside its pivot block), then Pe [e] = EMPTY,
 *	    and e is the root of an assembly subtree (or the whole tree if
 *	    there is just one such root).
 *
 *	Dense variable i:  if i is "dense", then Pe [i] = EMPTY.
 *
 *	On output, Pe holds the assembly tree/forest, which implicitly
 *	represents a pivot order with identical fill-in as the actual order
 *	(via a depth-first search of the tree), as follows.  If Nv [i] > 0,
 *	then i represents a node in the assembly tree, and the parent of i is
 *	Pe [i], or EMPTY if i is a root.  If Nv [i] = 0, then (i, Pe [i])
 *	represents an edge in a subtree, the root of which is a node in the
 *	assembly tree.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 * Info:  A double array of size AMD_INFO.  If present, (that is, not NULL),
 *	then statistics about the ordering are returned in the Info array.
 *	See amd.h for a description.

 * -----------------------------------------------------------------------------
 * INPUT/MODIFIED (undefined on output):
 * -----------------------------------------------------------------------------
 *
 * Len:  An integer array of size n.  On input, Len [i] holds the number of
 *	entries in row i of the matrix, excluding the diagonal.  The contents
 *	of Len are undefined on output.
 *
 * Iw:  An integer array of size iwlen.  On input, Iw [0..pfree-1] holds the
 *	description of each row i in the matrix.  The matrix must be symmetric,
 *	and both upper and lower triangular parts must be present.  The
 *	diagonal must not be present.  Row i is held as follows:
 *
 *	    Len [i]:  the length of the row i data structure in the Iw array.
 *	    Iw [Pe [i] ... Pe [i] + Len [i] - 1]:
 *		the list of column indices for nonzeros in row i (simple
 *		supervariables), excluding the diagonal.  All supervariables
 *		start with one row/column each (supervariable i is just row i).
 *		If Len [i] is zero on input, then Pe [i] is ignored on input.
 *
 *	    Note that the rows need not be in any particular order, and there
 *	    may be empty space between the rows.
 *
 *	During execution, the supervariable i experiences fill-in.  This is
 *	represented by placing in i a list of the elements that cause fill-in
 *	in supervariable i:
 *
 *	    Len [i]:  the length of supervariable i in the Iw array.
 *	    Iw [Pe [i] ... Pe [i] + Elen [i] - 1]:
 *		the list of elements that contain i.  This list is kept short
 *		by removing absorbed elements.
 *	    Iw [Pe [i] + Elen [i] ... Pe [i] + Len [i] - 1]:
 *		the list of supervariables in i.  This list is kept short by
 *		removing nonprincipal variables, and any entry j that is also
 *		contained in at least one of the elements (j in Le) in the list
 *		for i (e in row i).
 *
 *	When supervariable i is selected as pivot, we create an element e of
 *	the same name (e=i):
 *
 *	    Len [e]:  the length of element e in the Iw array.
 *	    Iw [Pe [e] ... Pe [e] + Len [e] - 1]:
 *		the list of supervariables in element e.
 *
 *	An element represents the fill-in that occurs when supervariable i is
 *	selected as pivot (which represents the selection of row i and all
 *	non-principal variables whose principal variable is i).  We use the
 *	term Le to denote the set of all supervariables in element e.  Absorbed
 *	supervariables and elements are pruned from these lists when
 *	computationally convenient.
 *
 *  CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
 *  The contents of Iw are undefined on output.

 * -----------------------------------------------------------------------------
 * OUTPUT (need not be set on input):
 * -----------------------------------------------------------------------------
 *
 * Nv:  An integer array of size n.  During execution, ABS (Nv [i]) is equal to
 *	the number of rows that are represented by the principal supervariable
 *	i.  If i is a nonprincipal or dense variable, then Nv [i] = 0.
 *	Initially, Nv [i] = 1 for all i.  Nv [i] < 0 signifies that i is a
 *	principal variable in the pattern Lme of the current pivot element me.
 *	After element me is constructed, Nv [i] is set back to a positive
 *	value.
 *
 *	On output, Nv [i] holds the number of pivots represented by super
 *	row/column i of the original matrix, or Nv [i] = 0 for non-principal
 *	rows/columns.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 * Elen:  An integer array of size n.  See the description of Iw above.  At the
 *	start of execution, Elen [i] is set to zero for all rows i.  During
 *	execution, Elen [i] is the number of elements in the list for
 *	supervariable i.  When e becomes an element, Elen [e] = FLIP (esize) is
 *	set, where esize is the size of the element (the number of pivots, plus
 *	the number of nonpivotal entries).  Thus Elen [e] < EMPTY.
 *	Elen (i) = EMPTY set when variable i becomes nonprincipal.
 *
 *	For variables, Elen (i) >= EMPTY holds until just before the
 *	postordering and permutation vectors are computed.  For elements,
 *	Elen [e] < EMPTY holds.
 *
 *	On output, Elen [i] is the degree of the row/column in the Cholesky
 *	factorization of the permuted matrix, corresponding to the original row
 *	i, if i is a super row/column.  It is equal to EMPTY if i is
 *	non-principal.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 *	Note that the contents of Elen on output differ from the Fortran
 *	version (Elen holds the inverse permutation in the Fortran version,
 *	which is instead returned in the Next array in this C version,
 *	described below).
 *
 * Last: In a degree list, Last [i] is the supervariable preceding i, or EMPTY
 *	if i is the head of the list.  In a hash bucket, Last [i] is the hash
 *	key for i.
 *
 *	Last [Head [hash]] is also used as the head of a hash bucket if
 *	Head [hash] contains a degree list (see the description of Head,
 *	below).
 *
 *	On output, Last [0..n-1] holds the permutation.  That is, if
 *	i = Last [k], then row i is the kth pivot row (where k ranges from 0 to
 *	n-1).  Row Last [k] of A is the kth row in the permuted matrix, PAP'.
 *
 * Next: Next [i] is the supervariable following i in a link list, or EMPTY if
 *	i is the last in the list.  Used for two kinds of lists:  degree lists
 *	and hash buckets (a supervariable can be in only one kind of list at a
 *	time).
 *
 *	On output Next [0..n-1] holds the inverse permutation. 	That is, if
 *	k = Next [i], then row i is the kth pivot row. Row i of A appears as
 *	the (Next[i])-th row in the permuted matrix, PAP'.
 *
 *	Note that the contents of Next on output differ from the Fortran
 *	version (Next is undefined on output in the Fortran version).

 * -----------------------------------------------------------------------------
 * LOCAL WORKSPACE (not input or output - used only during execution):
 * -----------------------------------------------------------------------------
 *
 * Degree:  An integer array of size n.  If i is a supervariable, then
 *	Degree [i] holds the current approximation of the external degree of
 *	row i (an upper bound).  The external degree is the number of nonzeros
 *	in row i, minus ABS (Nv [i]), the diagonal part.  The bound is equal to
 *	the exact external degree if Elen [i] is less than or equal to two.
 *
 *	We also use the term "external degree" for elements e to refer to
 *	|Le \ Lme|.  If e is an element, then Degree [e] is |Le|, which is the
 *	degree of the off-diagonal part of the element e (not including the
 *	diagonal part).
 *
 * Head:   An integer array of size n.  Head is used for degree lists.
 *	Head [deg] is the first supervariable in a degree list.  All
 *	supervariables i in a degree list Head [deg] have the same approximate
 *	degree, namely, deg = Degree [i].  If the list Head [deg] is empty then
 *	Head [deg] = EMPTY.
 *
 *	During supervariable detection Head [hash] also serves as a pointer to
 *	a hash bucket.  If Head [hash] >= 0, there is a degree list of degree
 *	hash.  The hash bucket head pointer is Last [Head [hash]].  If
 *	Head [hash] = EMPTY, then the degree list and hash bucket are both
 *	empty.  If Head [hash] < EMPTY, then the degree list is empty, and
 *	FLIP (Head [hash]) is the head of the hash bucket.  After supervariable
 *	detection is complete, all hash buckets are empty, and the
 *	(Last [Head [hash]] = EMPTY) condition is restored for the non-empty
 *	degree lists.
 *
 * W:  An integer array of size n.  The flag array W determines the status of
 *	elements and variables, and the external degree of elements.
 *
 *	for elements:
 *	    if W [e] = 0, then the element e is absorbed.
 *	    if W [e] >= wflg, then W [e] - wflg is the size of the set
 *		|Le \ Lme|, in terms of nonzeros (the sum of ABS (Nv [i]) for
 *		each principal variable i that is both in the pattern of
 *		element e and NOT in the pattern of the current pivot element,
 *		me).
 *	    if wflg > W [e] > 0, then e is not absorbed and has not yet been
 *		seen in the scan of the element lists in the computation of
 *		|Le\Lme| in Scan 1 below.
 *
 *	for variables:
 *	    during supervariable detection, if W [j] != wflg then j is
 *	    not in the pattern of variable i.
 *
 *	The W array is initialized by setting W [i] = 1 for all i, and by
 *	setting wflg = 2.  It is reinitialized if wflg becomes too large (to
 *	ensure that wflg+n does not cause integer overflow).

 * -----------------------------------------------------------------------------
 * LOCAL INTEGERS:
 * -----------------------------------------------------------------------------
 */

    Int deg, degme, dext, lemax, e, elenme, eln, i, ilast, inext, j,
	jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
	nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, x, ok, ndense, ncmpa,
	dense, aggressive ;

    unsigned Int hash ;	    /* unsigned, so that hash % n is well defined.*/

/*
 * deg:		the degree of a variable or element
 * degme:	size, |Lme|, of the current element, me (= Degree [me])
 * dext:	external degree, |Le \ Lme|, of some element e
 * lemax:	largest |Le| seen so far (called dmax in Fortran version)
 * e:		an element
 * elenme:	the length, Elen [me], of element list of pivotal variable
 * eln:		the length, Elen [...], of an element list
 * hash:	the computed value of the hash function
 * i:		a supervariable
 * ilast:	the entry in a link list preceding i
 * inext:	the entry in a link list following i
 * j:		a supervariable
 * jlast:	the entry in a link list preceding j
 * jnext:	the entry in a link list, or path, following j
 * k:		the pivot order of an element or variable
 * knt1:	loop counter used during element construction
 * knt2:	loop counter used during element construction
 * knt3:	loop counter used during compression
 * lenj:	Len [j]
 * ln:		length of a supervariable list
 * me:		current supervariable being eliminated, and the current
 *		    element created by eliminating that supervariable
 * mindeg:	current minimum degree
 * nel:		number of pivots selected so far
 * nleft:	n - nel, the number of nonpivotal rows/columns remaining
 * nvi:		the number of variables in a supervariable i (= Nv [i])
 * nvj:		the number of variables in a supervariable j (= Nv [j])
 * nvpiv:	number of pivots in current element
 * slenme:	number of variables in variable list of pivotal variable
 * wbig:	= INT_MAX - n for the "int" version, LONG_MAX - n for the
 *		    "long" version.  wflg is not allowed to be >= wbig.
 * we:		W [e]
 * wflg:	used for flagging the W array.  See description of Iw.
 * wnvi:	wflg - Nv [i]
 * x:		either a supervariable or an element
 *
 * ok:		true if supervariable j can be absorbed into i
 * ndense:	number of "dense" rows/columns
 * dense:	rows/columns with initial degree > dense are considered "dense" 
 * aggressive:	true if aggressive absorption is being performed
 * ncmpa:	number of garbage collections

 * -----------------------------------------------------------------------------
 * LOCAL DOUBLES, used for statistical output only (except for alpha):
 * -----------------------------------------------------------------------------
 */

    double f, r, ndiv, s, nms_lu, nms_ldl, dmax, alpha, lnz, lnzme ;

/*
 * f:		nvpiv
 * r:		degme + nvpiv
 * ndiv:	number of divisions for LU or LDL' factorizations
 * s:		number of multiply-subtract pairs for LU factorization, for the
 *		    current element me
 * nms_lu	number of multiply-subtract pairs for LU factorization
 * nms_ldl	number of multiply-subtract pairs for LDL' factorization
 * dmax:	the largest number of entries in any column of L, including the
 *		    diagonal
 * alpha:	"dense" degree ratio
 * lnz:		the number of nonzeros in L (excluding the diagonal)
 * lnzme:	the number of nonzeros in L (excl. the diagonal) for the
 *		    current element me

 * -----------------------------------------------------------------------------
 * LOCAL "POINTERS" (indices into the Iw array)
 * -----------------------------------------------------------------------------
*/

    Int p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, psrc ;

/*
 * Any parameter (Pe [...] or pfree) or local variable starting with "p" (for
 * Pointer) is an index into Iw, and all indices into Iw use variables starting
 * with "p."  The only exception to this rule is the iwlen input argument.
 *
 * p:           pointer into lots of things
 * p1:          Pe [i] for some variable i (start of element list)
 * p2:          Pe [i] + Elen [i] -  1 for some variable i
 * p3:          index of first supervariable in clean list
 * p4:		
 * pdst:        destination pointer, for compression
 * pend:        end of memory to compress
 * pj:          pointer into an element or variable
 * pme:         pointer into the current element (pme1...pme2)
 * pme1:        the current element, me, is stored in Iw [pme1...pme2]
 * pme2:        the end of the current element
 * pn:          pointer into a "clean" variable, also used to compress
 * psrc:        source pointer, for compression
*/

/* ========================================================================== */
/*  INITIALIZATIONS */
/* ========================================================================== */

    /* Note that this restriction on iwlen is slightly more restrictive than
     * what is actually required in AMD_2.  AMD_2 can operate with no elbow
     * room at all, but it will be slow.  For better performance, at least
     * size-n elbow room is enforced. */
    ASSERT (iwlen >= pfree + n) ;
    ASSERT (n > 0) ;

    /* initialize output statistics */
    lnz = 0 ;
    ndiv = 0 ;
    nms_lu = 0 ;
    nms_ldl = 0 ;
    dmax = 1 ;
    me = EMPTY ;

    wflg = 2 ;
    mindeg = 0 ;
    ncmpa = 0 ;
    nel = 0 ;
    lemax = 0 ;		/* this is called dmax in the Fortran version */

#ifdef TEST_FOR_INTEGER_OVERFLOW
    /* for testing only */
    wbig = 3*n ;
#else
    /* normal operation */
    wbig = Int_MAX - n ;
#endif

    /* get control parameters */
    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = (Control [AMD_AGGRESSIVE] != 0) ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }
    if (alpha < 0)
    {
	/* no dense rows/columns */
	dense = n ;
    }
    else
    {
	dense = alpha * sqrt ((double) n) ;
    }
    dense = MAX (16, dense) ;
    dense = MIN (n,  dense) ;
    AMD_DEBUG1 (("AMD (debug), alpha %g, aggr. "ID"\n", alpha, aggressive)) ;

    for (i = 0 ; i < n ; i++)
    {
	Last [i] = EMPTY ;
	Head [i] = EMPTY ;
	Next [i] = EMPTY ;
	/* if seperate Hhead array is used for hash buckets: *
	Hhead [i] = EMPTY ;
	*/
	Nv [i] = 1 ;
	W [i] = 1 ;
	Elen [i] = 0 ;
	Degree [i] = Len [i] ;
    }

#ifndef NDEBUG
    AMD_DEBUG1 (("\n======Nel "ID"\n", nel)) ;
    AMD_dump (n, Pe, Iw, Len, iwlen, pfree, Nv, Next, Last,
		Head, Elen, Degree, W, -1) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* initialize degree lists and eliminate dense and empty rows */
    /* ---------------------------------------------------------------------- */

    ndense = 0 ;

    /* for (i = n-1 ; i >= 0 ; i--) */
    for (i = 0 ; i < n ; i++)
    {
	deg = Degree [i] ;
	ASSERT (deg >= 0 && deg < n) ;
	if (deg == 0)
	{

	    /* --------------------------------------------------------------
	     * we have a variable that can be eliminated at once because
	     * there is no off-diagonal non-zero in its row.  Note that
	     * Nv [i] = 1 for an empty variable i.  It is treated just
	     * the same as an eliminated element i.
	     * -------------------------------------------------------------- */

	    Elen [i] = FLIP (1) ;
	    nel++ ;
	    Pe [i] = EMPTY ;
	    W [i] = 0 ;

	}
	else if (deg > dense)
	{

	    /* --------------------------------------------------------------
	     * Dense variables are not treated as elements, but as unordered,
	     * non-principal variables that have no parent.  They do not take
	     * part in the postorder, since Nv [i] = 0.  Note that the Fortran
	     * version does not have this option.
	     * -------------------------------------------------------------- */

	    AMD_DEBUG1 (("Dense node "ID" degree "ID"\n", i, deg)) ;
	    ndense++ ;
	    Nv [i] = 0 ;		/* do not postorder this node */
	    Elen [i] = EMPTY ;
	    nel++ ;
	    Pe [i] = EMPTY ;

	}
	else
	{

	    /* --------------------------------------------------------------
	     * place i in the degree list corresponding to its degree
	     * -------------------------------------------------------------- */

	    inext = Head [deg] ;
	    ASSERT (inext >= EMPTY && inext < n) ;
	    if (inext != EMPTY) Last [inext] = i ;
	    Next [i] = inext ;
	    Head [deg] = i ;

	}
    }

/* ========================================================================== */
/* WHILE (selecting pivots) DO */
/* ========================================================================== */

    while (nel < n)
    {

#ifndef NDEBUG
	AMD_DEBUG1 (("\n======Nel "ID"\n", nel)) ;
	if (AMD_debug >= 2) AMD_dump (n, Pe, Iw, Len, iwlen, pfree, Nv, Next,
	    Last, Head, Elen, Degree, W, nel) ;
#endif

/* ========================================================================== */
/* GET PIVOT OF MINIMUM DEGREE */
/* ========================================================================== */

	/* ------------------------------------------------------------------ */
	/* find next supervariable for elimination */
	/* ------------------------------------------------------------------ */

	ASSERT (mindeg >= 0 && mindeg < n) ;
	for (deg = mindeg ; deg < n ; deg++)
	{
	    me = Head [deg] ;
	    if (me != EMPTY) break ;
	}
	mindeg = deg ;
	ASSERT (me >= 0 && me < n) ;
	AMD_DEBUG1 (("=================me: "ID"\n", me)) ;

	/* ------------------------------------------------------------------ */
	/* remove chosen variable from link list */
	/* ------------------------------------------------------------------ */

	inext = Next [me] ;
	ASSERT (inext >= EMPTY && inext < n) ;
	if (inext != EMPTY) Last [inext] = EMPTY ;
	Head [deg] = inext ;

	/* ------------------------------------------------------------------ */
	/* me represents the elimination of pivots nel to nel+Nv[me]-1. */
	/* place me itself as the first in this set. */
	/* ------------------------------------------------------------------ */

	elenme = Elen [me] ;
	nvpiv = Nv [me] ;
	ASSERT (nvpiv > 0) ;
	nel += nvpiv ;

/* ========================================================================== */
/* CONSTRUCT NEW ELEMENT */
/* ========================================================================== */

	/* ------------------------------------------------------------------
	 * At this point, me is the pivotal supervariable.  It will be
	 * converted into the current element.  Scan list of the pivotal
	 * supervariable, me, setting tree pointers and constructing new list
	 * of supervariables for the new element, me.  p is a pointer to the
	 * current position in the old list.
	 * ------------------------------------------------------------------ */

	/* flag the variable "me" as being in Lme by negating Nv [me] */
	Nv [me] = -nvpiv ;
	degme = 0 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	if (elenme == 0)
	{

	    /* -------------------------------------------------------------- */
	    /* construct the new element in place */
	    /* -------------------------------------------------------------- */

	    pme1 = Pe [me] ;
	    pme2 = pme1 - 1 ;

	    for (p = pme1 ; p <= pme1 + Len [me] - 1 ; p++)
	    {
		i = Iw [p] ;
		ASSERT (i >= 0 && i < n && Nv [i] >= 0) ;
		nvi = Nv [i] ;
		if (nvi > 0)
		{

		    /* ------------------------------------------------------ */
		    /* i is a principal variable not yet placed in Lme. */
		    /* store i in new list */
		    /* ------------------------------------------------------ */

                    /* flag i as being in Lme by negating Nv [i] */
		    degme += nvi ;
		    Nv [i] = -nvi ;
		    Iw [++pme2] = i ;

		    /* ------------------------------------------------------ */
		    /* remove variable i from degree list. */
		    /* ------------------------------------------------------ */

		    ilast = Last [i] ;
		    inext = Next [i] ;
		    ASSERT (ilast >= EMPTY && ilast < n) ;
		    ASSERT (inext >= EMPTY && inext < n) ;
		    if (inext != EMPTY) Last [inext] = ilast ;
		    if (ilast != EMPTY)
		    {
			Next [ilast] = inext ;
		    }
		    else
		    {
                        /* i is at the head of the degree list */
			ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			Head [Degree [i]] = inext ;
		    }
		}
	    }
	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* construct the new element in empty space, Iw [pfree ...] */
	    /* -------------------------------------------------------------- */

	    p = Pe [me] ;
	    pme1 = pfree ;
	    slenme = Len [me] - elenme ;

	    for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)
	    {

		if (knt1 > elenme)
		{
		    /* search the supervariables in me. */
		    e = me ;
		    pj = p ;
		    ln = slenme ;
		    AMD_DEBUG2 (("Search sv: "ID" "ID" "ID"\n", me,pj,ln)) ;
		}
		else
		{
                    /* search the elements in me. */
		    e = Iw [p++] ;
		    ASSERT (e >= 0 && e < n) ;
		    pj = Pe [e] ;
		    ln = Len [e] ;
		    AMD_DEBUG2 (("Search element e "ID" in me "ID"\n", e,me)) ;
		    ASSERT (Elen [e] < EMPTY && W [e] > 0 && pj >= 0) ;
		}
		ASSERT (ln >= 0 && (ln == 0 || (pj >= 0 && pj < iwlen))) ;

		/* ----------------------------------------------------------
		 * search for different supervariables and add them to the
		 * new list, compressing when necessary. this loop is
		 * executed once for each element in the list and once for
		 * all the supervariables in the list.
		 * ---------------------------------------------------------- */

		for (knt2 = 1 ; knt2 <= ln ; knt2++)
		{
		    i = Iw [pj++] ;
		    ASSERT (i >= 0 && i < n && (i == me || Elen [i] >= EMPTY)) ;
		    nvi = Nv [i] ;
		    AMD_DEBUG2 ((": "ID" "ID" "ID" "ID"\n",
				i, Elen [i], Nv [i], wflg)) ;

		    if (nvi > 0)
		    {

			/* -------------------------------------------------- */
			/* compress Iw, if necessary */
			/* -------------------------------------------------- */

			if (pfree >= iwlen)
			{

			    AMD_DEBUG1 (("GARBAGE COLLECTION\n")) ;

			    /* prepare for compressing Iw by adjusting pointers
			     * and lengths so that the lists being searched in
			     * the inner and outer loops contain only the
			     * remaining entries. */

			    Pe [me] = p ;
			    Len [me] -= knt1 ;
                            /* check if nothing left of supervariable me */
			    if (Len [me] == 0) Pe [me] = EMPTY ;
			    Pe [e] = pj ;
			    Len [e] = ln - knt2 ;
                            /* nothing left of element e */
			    if (Len [e] == 0) Pe [e] = EMPTY ;

			    ncmpa++ ;	/* one more garbage collection */

                            /* store first entry of each object in Pe */
                            /* FLIP the first entry in each object */
			    for (j = 0 ; j < n ; j++)
			    {
				pn = Pe [j] ;
				if (pn >= 0)
				{
				    ASSERT (pn >= 0 && pn < iwlen) ;
				    Pe [j] = Iw [pn] ;
				    Iw [pn] = FLIP (j) ;
				}
			    }

                            /* psrc/pdst point to source/destination */
			    psrc = 0 ;
			    pdst = 0 ;
			    pend = pme1 - 1 ;

			    while (psrc <= pend)
			    {
				/* search for next FLIP'd entry */
				j = FLIP (Iw [psrc++]) ;
				if (j >= 0)
				{
				    AMD_DEBUG2 (("Got object j: "ID"\n", j)) ;
				    Iw [pdst] = Pe [j] ;
				    Pe [j] = pdst++ ;
				    lenj = Len [j] ;
				    /* copy from source to destination */
				    for (knt3 = 0 ; knt3 <= lenj - 2 ; knt3++)
				    {
					Iw [pdst++] = Iw [psrc++] ;
				    }
				}
			    }

                            /* move the new partially-constructed element */
			    p1 = pdst ;
			    for (psrc = pme1 ; psrc <= pfree-1 ; psrc++)
			    {
				Iw [pdst++] = Iw [psrc] ;
			    }
			    pme1 = p1 ;
			    pfree = pdst ;
			    pj = Pe [e] ;
			    p = Pe [me] ;

			}

			/* -------------------------------------------------- */
			/* i is a principal variable not yet placed in Lme */
			/* store i in new list */
			/* -------------------------------------------------- */

                        /* flag i as being in Lme by negating Nv [i] */
			degme += nvi ;
			Nv [i] = -nvi ;
			Iw [pfree++] = i ;
			AMD_DEBUG2 (("     s: "ID"     nv "ID"\n", i, Nv [i])) ;

			/* -------------------------------------------------- */
			/* remove variable i from degree link list */
			/* -------------------------------------------------- */

			ilast = Last [i] ;
			inext = Next [i] ;
			ASSERT (ilast >= EMPTY && ilast < n) ;
			ASSERT (inext >= EMPTY && inext < n) ;
			if (inext != EMPTY) Last [inext] = ilast ;
			if (ilast != EMPTY)
			{
			    Next [ilast] = inext ;
			}
			else
			{
                            /* i is at the head of the degree list */
			    ASSERT (Degree [i] >= 0 && Degree [i] < n) ;
			    Head [Degree [i]] = inext ;
			}
		    }
		}

		if (e != me)
		{
                    /* set tree pointer and flag to indicate element e is
                     * absorbed into new element me (the parent of e is me) */
		    AMD_DEBUG1 ((" Element "ID" => "ID"\n", e, me)) ;
		    Pe [e] = FLIP (me) ;
		    W [e] = 0 ;
		}
	    }

	    pme2 = pfree - 1 ;
	}

	/* ------------------------------------------------------------------ */
	/* me has now been converted into an element in Iw [pme1..pme2] */
	/* ------------------------------------------------------------------ */

        /* degme holds the external degree of new element */
	Degree [me] = degme ;
	Pe [me] = pme1 ;
	Len [me] = pme2 - pme1 + 1 ;
	ASSERT (Pe [me] >= 0 && Pe [me] < iwlen) ;

	Elen [me] = FLIP (nvpiv + degme) ;
        /* FLIP (Elen (me)) is now the degree of pivot (including
	 * diagonal part). */

#ifndef NDEBUG
	AMD_DEBUG2 (("New element structure: length= "ID"\n", pme2-pme1+1)) ;
	for (pme = pme1 ; pme <= pme2 ; pme++) AMD_DEBUG3 ((" "ID"", Iw [pme]));
	AMD_DEBUG3 (("\n")) ;
#endif

	/* ------------------------------------------------------------------ */
	/* make sure that wflg is not too large. */
	/* ------------------------------------------------------------------ */

	/* With the current value of wflg, wflg+n must not cause integer
	 * overflow */

	if (wflg >= wbig)
	{
	    for (x = 0 ; x < n ; x++)
	    {
		if (W [x] != 0) W [x] = 1 ;
	    }
	    wflg = 2 ;
	}

/* ========================================================================== */
/* COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS */
/* ========================================================================== */

	/* ------------------------------------------------------------------
	 * Scan 1:  compute the external degrees of previous elements with
	 * respect to the current element.  That is:
	 *       (W [e] - wflg) = |Le \ Lme|
	 * for each element e that appears in any supervariable in Lme.  The
	 * notation Le refers to the pattern (list of supervariables) of a
	 * previous element e, where e is not yet absorbed, stored in
	 * Iw [Pe [e] + 1 ... Pe [e] + Iw [Pe [e]]].  The notation Lme
	 * refers to the pattern of the current element (stored in
	 * Iw [pme1..pme2]).   If aggressive absorption is enabled, and
	 * (W [e] - wflg) becomes zero, then the element e will be absorbed
	 * in Scan 2.
	 * ------------------------------------------------------------------ */

	AMD_DEBUG2 (("me: ")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    eln = Elen [i] ;
	    AMD_DEBUG3 ((""ID" Elen "ID": \n", i, eln)) ;
	    if (eln > 0)
	    {
                /* note that Nv [i] has been negated to denote i in Lme: */
		nvi = -Nv [i] ;
		ASSERT (nvi > 0 && Pe [i] >= 0 && Pe [i] < iwlen) ;
		wnvi = wflg - nvi ;
		for (p = Pe [i] ; p <= Pe [i] + eln - 1 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    AMD_DEBUG4 (("    e "ID" we "ID" ", e, we)) ;
		    if (we >= wflg)
		    {
                        /* unabsorbed element e has been seen in this loop */
			AMD_DEBUG4 (("    unabsorbed, first time seen")) ;
			we -= nvi ;
		    }
		    else if (we != 0)
		    {
                        /* e is an unabsorbed element */
                        /* this is the first we have seen e in all of Scan 1 */
			AMD_DEBUG4 (("    unabsorbed")) ;
			we = Degree [e] + wnvi ;
		    }
		    AMD_DEBUG4 (("\n")) ;
		    W [e] = we ;
		}
	    }
	}
	AMD_DEBUG2 (("\n")) ;

/* ========================================================================== */
/* DEGREE UPDATE AND ELEMENT ABSORPTION */
/* ========================================================================== */

	/* ------------------------------------------------------------------
	 * Scan 2:  for each i in Lme, sum up the degree of Lme (which is
	 * degme), plus the sum of the external degrees of each Le for the
	 * elements e appearing within i, plus the supervariables in i.
	 * Place i in hash list.
	 * ------------------------------------------------------------------ */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n && Nv [i] < 0 && Elen [i] >= 0) ;
	    AMD_DEBUG2 (("Updating: i "ID" "ID" "ID"\n", i, Elen[i], Len [i])) ;
	    p1 = Pe [i] ;
	    p2 = p1 + Elen [i] - 1 ;
	    pn = p1 ;
	    hash = 0 ;
	    deg = 0 ;
	    ASSERT (p1 >= 0 && p1 < iwlen && p2 >= -1 && p2 < iwlen) ;

	    /* -------------------------------------------------------------- */
	    /* scan the element list associated with supervariable i */
	    /* -------------------------------------------------------------- */

            /* UMFPACK/MA38-style approximate degree: */
	    if (aggressive)
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* dext = | Le \ Lme | */
			dext = we - wflg ;
			if (dext > 0)
			{
			    deg += dext ;
			    Iw [pn++] = e ;
			    hash += e ;
			    AMD_DEBUG4 ((" e: "ID" hash = "ID"\n",e,hash)) ;
			}
			else
			{
			    /* external degree of e is zero, absorb e into me */
			    AMD_DEBUG1 ((" Element "ID" => "ID" (aggressive)\n",
				e, me)) ;
			    ASSERT (dext == 0) ;
			    Pe [e] = FLIP (me) ;
			    W [e] = 0 ;
			}
		    }
		}
	    }
	    else
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    ASSERT (e >= 0 && e < n) ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			dext = we - wflg ;
			ASSERT (dext >= 0) ;
			deg += dext ;
			Iw [pn++] = e ;
			hash += e ;
			AMD_DEBUG4 (("	e: "ID" hash = "ID"\n",e,hash)) ;
		    }
		}
	    }

            /* count the number of elements in i (including me): */
	    Elen [i] = pn - p1 + 1 ;

	    /* -------------------------------------------------------------- */
	    /* scan the supervariables in the list associated with i */
	    /* -------------------------------------------------------------- */

	    /* The bulk of the AMD run time is typically spent in this loop,
	     * particularly if the matrix has many dense rows that are not
	     * removed prior to ordering. */
	    p3 = pn ;
	    p4 = p1 + Len [i] ;
	    for (p = p2 + 1 ; p < p4 ; p++)
	    {
		j = Iw [p] ;
		ASSERT (j >= 0 && j < n) ;
		nvj = Nv [j] ;
		if (nvj > 0)
		{
                    /* j is unabsorbed, and not in Lme. */
                    /* add to degree and add to new list */
		    deg += nvj ;
		    Iw [pn++] = j ;
		    hash += j ;
		    AMD_DEBUG4 (("  s: "ID" hash "ID" Nv[j]= "ID"\n",
				j, hash, nvj)) ;
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* update the degree and check for mass elimination */
	    /* -------------------------------------------------------------- */

	    /* with aggressive absorption, deg==0 is identical to the 
	     * Elen [i] == 1 && p3 == pn test, below. */
	    ASSERT (IMPLIES (aggressive, (deg==0) == (Elen[i]==1 && p3==pn))) ;

	    if (Elen [i] == 1 && p3 == pn)
	    {

		/* ---------------------------------------------------------- */
	        /* mass elimination */
		/* ---------------------------------------------------------- */

		/* There is nothing left of this node except for an edge to
		 * the current pivot element.  Elen [i] is 1, and there are
		 * no variables adjacent to node i.  Absorb i into the
		 * current pivot element, me.  Note that if there are two or
		 * more mass eliminations, fillin due to mass elimination is
		 * possible within the nvpiv-by-nvpiv pivot block.  It is this
		 * step that causes AMD's analysis to be an upper bound.
		 *
		 * The reason is that the selected pivot has a lower
		 * approximate degree than the true degree of the two mass
		 * eliminated nodes.  There is no edge between the two mass
		 * eliminated nodes.  They are merged with the current pivot
		 * anyway.
		 *
		 * No fillin occurs in the Schur complement, in any case,
		 * and this effect does not decrease the quality of the
		 * ordering itself, just the quality of the nonzero and
		 * flop count analysis.  It also means that the post-ordering
		 * is not an exact elimination tree post-ordering. */

		AMD_DEBUG1 (("  MASS i "ID" => parent e "ID"\n", i, me)) ;
		Pe [i] = FLIP (me) ;
		nvi = -Nv [i] ;
		degme -= nvi ;
		nvpiv += nvi ;
		nel += nvi ;
		Nv [i] = 0 ;
		Elen [i] = EMPTY ;

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* update the upper-bound degree of i */
		/* ---------------------------------------------------------- */

                /* the following degree does not yet include the size
                 * of the current element, which is added later: */

		Degree [i] = MIN (Degree [i], deg) ;

		/* ---------------------------------------------------------- */
		/* add me to the list for i */
		/* ---------------------------------------------------------- */

                /* move first supervariable to end of list */
		Iw [pn] = Iw [p3] ;
                /* move first element to end of element part of list */
		Iw [p3] = Iw [p1] ;
                /* add new element, me, to front of list. */
		Iw [p1] = me ;
                /* store the new length of the list in Len [i] */
		Len [i] = pn - p1 + 1 ;

		/* ---------------------------------------------------------- */
		/* place in hash bucket.  Save hash key of i in Last [i]. */
		/* ---------------------------------------------------------- */

		/* NOTE: this can fail if hash is negative, because the ANSI C
		 * standard does not define a % b when a and/or b are negative.
		 * That's why hash is defined as an unsigned Int, to avoid this
		 * problem. */
		hash = hash % n ;
		ASSERT (((Int) hash) >= 0 && ((Int) hash) < n) ;

		/* if the Hhead array is not used: */
		j = Head [hash] ;
		if (j <= EMPTY)
		{
		    /* degree list is empty, hash head is FLIP (j) */
		    Next [i] = FLIP (j) ;
		    Head [hash] = FLIP (i) ;
		}
		else
		{
		    /* degree list is not empty, use Last [Head [hash]] as
		     * hash head. */
		    Next [i] = Last [j] ;
		    Last [j] = i ;
		}

		/* if a seperate Hhead array is used: *
		Next [i] = Hhead [hash] ;
		Hhead [hash] = i ;
		*/

		Last [i] = hash ;
	    }
	}

	Degree [me] = degme ;

	/* ------------------------------------------------------------------ */
	/* Clear the counter array, W [...], by incrementing wflg. */
	/* ------------------------------------------------------------------ */

        /* make sure that wflg+n does not cause integer overflow */
	lemax =  MAX (lemax, degme) ;
	wflg += lemax ;
	if (wflg >= wbig)
	{
	    for (x = 0 ; x < n ; x++)
	    {
		if (W [x] != 0) W [x] = 1 ;
	    }
	    wflg = 2 ;
	}
        /*  at this point, W [0..n-1] < wflg holds */

/* ========================================================================== */
/* SUPERVARIABLE DETECTION */
/* ========================================================================== */

	AMD_DEBUG1 (("Detecting supervariables:\n")) ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    AMD_DEBUG2 (("Consider i "ID" nv "ID"\n", i, Nv [i])) ;
	    if (Nv [i] < 0)
	    {
                /* i is a principal variable in Lme */

		/* ----------------------------------------------------------
		 * examine all hash buckets with 2 or more variables.  We do
		 * this by examing all unique hash keys for supervariables in
		 * the pattern Lme of the current element, me
		 * ---------------------------------------------------------- */

                /* let i = head of hash bucket, and empty the hash bucket */
		ASSERT (Last [i] >= 0 && Last [i] < n) ;
		hash = Last [i] ;

		/* if Hhead array is not used: */
		j = Head [hash] ;
		if (j == EMPTY)
		{
		    /* hash bucket and degree list are both empty */
		    i = EMPTY ;
		}
		else if (j < EMPTY)
		{
		    /* degree list is empty */
		    i = FLIP (j) ;
		    Head [hash] = EMPTY ;
		}
		else
		{
		    /* degree list is not empty, restore Last [j] of head j */
		    i = Last [j] ;
		    Last [j] = EMPTY ;
		}

		/* if seperate Hhead array is used: *
		i = Hhead [hash] ;
		Hhead [hash] = EMPTY ;
		*/

		ASSERT (i >= EMPTY && i < n) ;
		AMD_DEBUG2 (("----i "ID" hash "ID"\n", i, hash)) ;

		while (i != EMPTY && Next [i] != EMPTY)
		{

		    /* ------------------------------------------------------
		     * this bucket has one or more variables following i.
		     * scan all of them to see if i can absorb any entries
		     * that follow i in hash bucket.  Scatter i into w.
		     * ------------------------------------------------------ */

		    ln = Len [i] ;
		    eln = Elen [i] ;
		    ASSERT (ln >= 0 && eln >= 0) ;
		    ASSERT (Pe [i] >= 0 && Pe [i] < iwlen) ;
                    /* do not flag the first element in the list (me) */
		    for (p = Pe [i] + 1 ; p <= Pe [i] + ln - 1 ; p++)
		    {
			ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			W [Iw [p]] = wflg ;
		    }

		    /* ------------------------------------------------------ */
		    /* scan every other entry j following i in bucket */
		    /* ------------------------------------------------------ */

		    jlast = i ;
		    j = Next [i] ;
		    ASSERT (j >= EMPTY && j < n) ;

		    while (j != EMPTY)
		    {
			/* -------------------------------------------------- */
			/* check if j and i have identical nonzero pattern */
			/* -------------------------------------------------- */

			AMD_DEBUG3 (("compare i "ID" and j "ID"\n", i,j)) ;

			/* check if i and j have the same Len and Elen */
			ASSERT (Len [j] >= 0 && Elen [j] >= 0) ;
			ASSERT (Pe [j] >= 0 && Pe [j] < iwlen) ;
			ok = (Len [j] == ln) && (Elen [j] == eln) ;
                        /* skop the first element in the list (me) */
			for (p = Pe [j] + 1 ; ok && p <= Pe [j] + ln - 1 ; p++)
			{
			    ASSERT (Iw [p] >= 0 && Iw [p] < n) ;
			    if (W [Iw [p]] != wflg) ok = 0 ;
			}
			if (ok)
			{
			    /* ---------------------------------------------- */
			    /* found it!  j can be absorbed into i */
			    /* ---------------------------------------------- */

			    AMD_DEBUG1 (("found it! j "ID" => i "ID"\n", j,i)) ;
			    Pe [j] = FLIP (i) ;
			    /* both Nv [i] and Nv [j] are negated since they */
			    /* are in Lme, and the absolute values of each */
			    /* are the number of variables in i and j: */
			    Nv [i] += Nv [j] ;
			    Nv [j] = 0 ;
			    Elen [j] = EMPTY ;
			    /* delete j from hash bucket */
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			    Next [jlast] = j ;

			}
			else
			{
			    /* j cannot be absorbed into i */
			    jlast = j ;
			    ASSERT (j != Next [j]) ;
			    j = Next [j] ;
			}
			ASSERT (j >= EMPTY && j < n) ;
		    }

		    /* ------------------------------------------------------
		     * no more variables can be absorbed into i
		     * go to next i in bucket and clear flag array
		     * ------------------------------------------------------ */

		    wflg++ ;
		    i = Next [i] ;
		    ASSERT (i >= EMPTY && i < n) ;

		}
	    }
	}
	AMD_DEBUG2 (("detect done\n")) ;

/* ========================================================================== */
/* RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT */
/* ========================================================================== */

	p = pme1 ;
	nleft = n - nel ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    ASSERT (i >= 0 && i < n) ;
	    nvi = -Nv [i] ;
	    AMD_DEBUG3 (("Restore i "ID" "ID"\n", i, nvi)) ;
	    if (nvi > 0)
	    {
                /* i is a principal variable in Lme */
                /* restore Nv [i] to signify that i is principal */
		Nv [i] = nvi ;

		/* ---------------------------------------------------------- */
		/* compute the external degree (add size of current element) */
		/* ---------------------------------------------------------- */

		deg = Degree [i] + degme - nvi ;
		deg = MIN (deg, nleft - nvi) ;
		ASSERT (IMPLIES (aggressive, deg > 0) && deg >= 0 && deg < n) ;

		/* ---------------------------------------------------------- */
		/* place the supervariable at the head of the degree list */
		/* ---------------------------------------------------------- */

		inext = Head [deg] ;
		ASSERT (inext >= EMPTY && inext < n) ;
		if (inext != EMPTY) Last [inext] = i ;
		Next [i] = inext ;
		Last [i] = EMPTY ;
		Head [deg] = i ;

		/* ---------------------------------------------------------- */
		/* save the new degree, and find the minimum degree */
		/* ---------------------------------------------------------- */

		mindeg = MIN (mindeg, deg) ;
		Degree [i] = deg ;

		/* ---------------------------------------------------------- */
		/* place the supervariable in the element pattern */
		/* ---------------------------------------------------------- */

		Iw [p++] = i ;

	    }
	}
	AMD_DEBUG2 (("restore done\n")) ;

/* ========================================================================== */
/* FINALIZE THE NEW ELEMENT */
/* ========================================================================== */

	AMD_DEBUG2 (("ME = "ID" DONE\n", me)) ;
	Nv [me] = nvpiv ;
        /* save the length of the list for the new element me */
	Len [me] = p - pme1 ;
	if (Len [me] == 0)
	{
            /* there is nothing left of the current pivot element */
	    /* it is a root of the assembly tree */
	    Pe [me] = EMPTY ;
	    W [me] = 0 ;
	}
	if (elenme != 0)
	{
	    /* element was not constructed in place: deallocate part of */
	    /* it since newly nonprincipal variables may have been removed */
	    pfree = p ;
	}

	/* The new element has nvpiv pivots and the size of the contribution
	 * block for a multifrontal method is degme-by-degme, not including
	 * the "dense" rows/columns.  If the "dense" rows/columns are included,
	 * the frontal matrix is no larger than
	 * (degme+ndense)-by-(degme+ndense).
	 */

	if (Info != (double *) NULL)
	{
	    f = nvpiv ;
	    r = degme + ndense ;
	    dmax = MAX (dmax, f + r) ;

	    /* number of nonzeros in L (excluding the diagonal) */
	    lnzme = f*r + (f-1)*f/2 ;
	    lnz += lnzme ;

	    /* number of divide operations for LDL' and for LU */
	    ndiv += lnzme ;

	    /* number of multiply-subtract pairs for LU */
	    s = f*r*r + r*(f-1)*f + (f-1)*f*(2*f-1)/6 ;
	    nms_lu += s ;

	    /* number of multiply-subtract pairs for LDL' */
	    nms_ldl += (s + lnzme)/2 ;
	}

#ifndef NDEBUG
	AMD_DEBUG2 (("finalize done nel "ID" n "ID"\n   ::::\n", nel, n)) ;
	for (pme = Pe [me] ; pme <= Pe [me] + Len [me] - 1 ; pme++)
	{
	      AMD_DEBUG3 ((" "ID"", Iw [pme])) ;
	}
	AMD_DEBUG3 (("\n")) ;
#endif

    }

/* ========================================================================== */
/* DONE SELECTING PIVOTS */
/* ========================================================================== */

    if (Info != (double *) NULL)
    {

	/* count the work to factorize the ndense-by-ndense submatrix */
	f = ndense ;
	dmax = MAX (dmax, (double) ndense) ;

	/* number of nonzeros in L (excluding the diagonal) */
	lnzme = (f-1)*f/2 ;
	lnz += lnzme ;

	/* number of divide operations for LDL' and for LU */
	ndiv += lnzme ;

	/* number of multiply-subtract pairs for LU */
	s = (f-1)*f*(2*f-1)/6 ;
	nms_lu += s ;

	/* number of multiply-subtract pairs for LDL' */
	nms_ldl += (s + lnzme)/2 ;

	/* number of nz's in L (excl. diagonal) */
	Info [AMD_LNZ] = lnz ;

	/* number of divide ops for LU and LDL' */
	Info [AMD_NDIV] = ndiv ;

	/* number of multiply-subtract pairs for LDL' */
	Info [AMD_NMULTSUBS_LDL] = nms_ldl ;

	/* number of multiply-subtract pairs for LU */
	Info [AMD_NMULTSUBS_LU] = nms_lu ;

	/* number of "dense" rows/columns */
	Info [AMD_NDENSE] = ndense ;

	/* largest front is dmax-by-dmax */
	Info [AMD_DMAX] = dmax ;

	/* number of garbage collections in AMD */
	Info [AMD_NCMPA] = ncmpa ;

	/* successful ordering */
	Info [AMD_STATUS] = AMD_OK ;
    }

/* --------------------------------------------------------------------------
 * Variables at this point:
 *
 * Pe: holds the elimination tree.  The parent of j is FLIP (Pe [j]),
 *	or EMPTY if j is a root.  The tree holds both elements and
 *	non-principal (unordered) variables absorbed into them.
 *	Dense variables are non-principal and unordered.
 *
 * Elen: holds the size of each element, including the diagonal part.
 *	FLIP (Elen [e]) > 0 if e is an element.  For unordered
 *	variables i, Elen [i] is EMPTY.
 *
 * Nv: Nv [e] > 0 is the number of pivots represented by the element e.
 *	For unordered variables i, Nv [i] is zero.
 *
 * Contents no longer needed:
 *	W, Iw, Len, Degree, Head, Next, Last.
 *
 * The matrix itself has been destroyed.
 *
 * n: the size of the matrix.
 * No other scalars needed (pfree, iwlen, etc.)
 * -------------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	Pe [i] = FLIP (Pe [i]) ;
	Elen [i] = FLIP (Elen [i]) ;
    }

/* Now the parent of j is Pe [j], or EMPTY if j is a root.  Elen [e] > 0
 * is the size of element e.  Elen [i] is EMPTY for unordered variable i. */

#ifndef NDEBUG
    AMD_DEBUG2 (("\nTree:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	AMD_DEBUG2 ((" "ID" parent: "ID"   ", i, Pe [i])) ;
	ASSERT (Pe [i] >= EMPTY && Pe [i] < n) ;
	if (Nv [i] > 0)
	{
	    /* this is an element */
	    e = i ;
	    AMD_DEBUG2 ((" element, size is "ID"\n", Elen [i])) ;
	    ASSERT (Elen [e] > 0) ;
	}
	AMD_DEBUG2 (("\n")) ;
    }
    AMD_DEBUG2 (("\nelements:\n")) ;
    for (e = 0 ; e < n ; e++)
    {
	if (Nv [e] > 0)
	{
	    AMD_DEBUG3 (("Element e= "ID" size "ID" nv "ID" \n", e,
		Elen [e], Nv [e])) ;
	}
    }
    AMD_DEBUG2 (("\nvariables:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	Int cnt ;
	if (Nv [i] == 0)
	{
	    AMD_DEBUG3 (("i unordered: "ID"\n", i)) ;
	    j = Pe [i] ;
	    cnt = 0 ;
	    AMD_DEBUG3 (("  j: "ID"\n", j)) ;
	    if (j == EMPTY)
	    {
		AMD_DEBUG3 (("	i is a dense variable\n")) ;
	    }
	    else
	    {
		ASSERT (j >= 0 && j < n) ;
		while (Nv [j] == 0)
		{
		    AMD_DEBUG3 (("	j : "ID"\n", j)) ;
		    j = Pe [j] ;
		    AMD_DEBUG3 (("	j:: "ID"\n", j)) ;
		    cnt++ ;
		    if (cnt > n) break ;
		}
		e = j ;
		AMD_DEBUG3 (("	got to e: "ID"\n", e)) ;
	    }
	}
    }
#endif

/* ========================================================================== */
/* compress the paths of the variables */
/* ========================================================================== */

    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{

	    /* --------------------------------------------------------------
	     * i is an un-ordered row.  Traverse the tree from i until
	     * reaching an element, e.  The element, e, was the principal
	     * supervariable of i and all nodes in the path from i to when e
	     * was selected as pivot.
	     * -------------------------------------------------------------- */

	    AMD_DEBUG1 (("Path compression, i unordered: "ID"\n", i)) ;
	    j = Pe [i] ;
	    ASSERT (j >= EMPTY && j < n) ;
	    AMD_DEBUG3 (("	j: "ID"\n", j)) ;
	    if (j == EMPTY)
	    {
		/* Skip a dense variable.  It has no parent. */
		AMD_DEBUG3 (("      i is a dense variable\n")) ;
		continue ;
	    }

	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		AMD_DEBUG3 (("		j : "ID"\n", j)) ;
		j = Pe [j] ;
		AMD_DEBUG3 (("		j:: "ID"\n", j)) ;
		ASSERT (j >= 0 && j < n) ;
	    }
	    /* got to an element e */
	    e = j ;
	    AMD_DEBUG3 (("got to e: "ID"\n", e)) ;

	    /* --------------------------------------------------------------
	     * traverse the path again from i to e, and compress the path
	     * (all nodes point to e).  Path compression allows this code to
	     * compute in O(n) time.
	     * -------------------------------------------------------------- */

	    j = i ;
	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		jnext = Pe [j] ;
		AMD_DEBUG3 (("j "ID" jnext "ID"\n", j, jnext)) ;
		Pe [j] = e ;
		j = jnext ;
		ASSERT (j >= 0 && j < n) ;
	    }
	}
    }

/* ========================================================================== */
/* postorder the assembly tree */
/* ========================================================================== */

    AMD_postorder (n, Pe, Nv, Elen,
	W,			/* output order */
	Head, Next, Last) ;	/* workspace */

/* ========================================================================== */
/* compute output permutation and inverse permutation */
/* ========================================================================== */

    /* W [e] = k means that element e is the kth element in the new
     * order.  e is in the range 0 to n-1, and k is in the range 0 to
     * the number of elements.  Use Head for inverse order. */

    for (k = 0 ; k < n ; k++)
    {
	Head [k] = EMPTY ;
	Next [k] = EMPTY ;
    }
    for (e = 0 ; e < n ; e++)
    {
	k = W [e] ;
	ASSERT ((k == EMPTY) == (Nv [e] == 0)) ;
	if (k != EMPTY)
	{
	    ASSERT (k >= 0 && k < n) ;
	    Head [k] = e ;
	}
    }

    /* construct output inverse permutation in Next,
     * and permutation in Last */
    nel = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	e = Head [k] ;
	if (e == EMPTY) break ;
	ASSERT (e >= 0 && e < n && Nv [e] > 0) ;
	Next [e] = nel ;
	nel += Nv [e] ;
    }
    ASSERT (nel == n - ndense) ;

    /* order non-principal variables (dense, and those merged into supervar's */
    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{
	    e = Pe [i] ;
	    ASSERT (e >= EMPTY && e < n) ;
	    if (e != EMPTY)
	    {
		/* This is an unordered variable that was merged
		 * into element e via supernode detection or mass
		 * elimination of i when e became the pivot element.
		 * Place i in order just before e. */
		ASSERT (Next [i] == EMPTY && Nv [e] > 0) ;
		Next [i] = Next [e] ;
		Next [e]++ ;
	    }
	    else
	    {
		/* This is a dense unordered variable, with no parent.
		 * Place it last in the output order. */
		Next [i] = nel++ ;
	    }
	}
    }
    ASSERT (nel == n) ; 

    AMD_DEBUG2 (("\n\nPerm:\n")) ;
    for (i = 0 ; i < n ; i++)
    {
	k = Next [i] ;
	ASSERT (k >= 0 && k < n) ;
	Last [k] = i ;
	AMD_DEBUG2 (("   perm ["ID"] = "ID"\n", k, i)) ;
    }
}
