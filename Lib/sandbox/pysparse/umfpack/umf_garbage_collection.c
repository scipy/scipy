/* ========================================================================== */
/* === UMF_garbage_collection =============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Compress the elements at the tail of Numeric->Memory, and delete the tuples.
    Elements are renumbered.  The new numbering space is compressed, and
    in the order of element creation (original elements of A first, followed
    by the new elements in the order that they were formed).

    Only called by UMF_get_memory.

    There are 5 ways in which garbage collection can be performed:

	Allocate a new working array for the current frontal matrix.  In this
	case, there are never any pivot rows/columns in the current frontal
	matrix (fnpiv = 0), and the old working array for the current frontal
	matrix can always be fully compacted, to fnrows-by-fncols.

	    UMF_kernel : UMF_extend : UMF_grow_front : UMF_get_memory
	    UMF_kernel : UMF_init_front : UMF_grow_front : UMF_get_memory
	    UMF_kernel : UMF_start_front : UMF_grow_front : UMF_get_memory

	Allocate a new element.  In this case, UMF_grow_front may or may not
	be subsequently called, depending on Work->do_grow.  There are never
	any pivot rows/columns in the current frontal matrix (fnpiv=0), but one
	may be added if UMF_init_front is to be called just after
	UMF_create_element.  If do_grow is true, then the current front can be
	fully compacted, to fnrows-by-fncols.  Otherwise, it can only be
	partially compacted, to MAX (fnrows, fnrows_new + 1) -by-
	MAX (fncols, fncols_new + 1).

	    UMF_kernel : UMF_create_element : UMF_get_memory

	Allocate rows of L and columns of U.  In this case, the current
	frontal matrix is only partially compacted, to (fnrows_new + 1)-by-
	(fncols_new + 1).  There are pivots in the frontal matrix (fnpiv > 0).

	    UMF_kernel : UMF_store_lu : UMF_get_memory
*/

#include "umf_internal.h"

GLOBAL void UMF_garbage_collection
(
    NumericType *Numeric,
    WorkType *Work,
    Int drnew,	    /* compact current front to drnew-by-dcnew */
    Int dcnew,
    Int do_Fcpos
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int size, e, n_row, n_col, nrows, ncols, nrowsleft, ncolsleft, prevsize,
	csize, size2, i2, j2, i, j, cdeg, rdeg, *E, row, col,
	*Rows, *Cols, *Rows2, *Cols2, nel, e2, *Row_tuples, *Col_tuples,
	*Row_degree, *Col_degree ;
    Entry *C, *C1, *C3, *C2 ;
    Unit *psrc, *pdest, *p, *pnext ;
    Element *epsrc, *epdest ;

#ifndef NDEBUG
    Int nmark ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get parameters */
    /* ---------------------------------------------------------------------- */

    Col_degree = Numeric->Cperm ;	/* for NON_PIVOTAL_COL macro */
    Row_degree = Numeric->Rperm ;	/* for NON_PIVOTAL_ROW macro */
    Row_tuples = Numeric->Uip ;
    Col_tuples = Numeric->Lip ;
    E = Work->E ;
    n_row = Work->n_row ;
    n_col = Work->n_col ;

    /* note that the tuple lengths (Col_tlen and Row_tlen) are updated, but */
    /* the tuple lists themselves are stale and are about to be destroyed */
    /* and recreated.  Do not attempt to scan them until they are recreated. */

#ifndef NDEBUG
    DEBUGm1 (("::::GARBAGE COLLECTION::::\n")) ;
    UMF_dump_memory (Numeric) ;
#endif

    Numeric->ngarbage++ ;

    /* ---------------------------------------------------------------------- */
    /* delete the tuple lists by marking the blocks as free */
    /* ---------------------------------------------------------------------- */

    /* do not modify Row_tlen and Col_tlen */
    /* those are needed for UMF_build_tuples */

    for (row = 0 ; row < n_row ; row++)
    {
	if (NON_PIVOTAL_ROW (row) && Row_tuples [row])
	{
	    DEBUG2 (("row "ID" tuples "ID"\n", row, Row_tuples [row])) ;
	    p = Numeric->Memory + Row_tuples [row] - 1 ;
	    DEBUG2 (("Freeing tuple list row "ID", p-S "ID", size "ID"\n",
		row, (Int) (p-Numeric->Memory), p->header.size)) ;
	    ASSERT (p->header.size > 0) ;
	    ASSERT (p >= Numeric->Memory + Numeric->itail) ;
	    ASSERT (p < Numeric->Memory + Numeric->size) ;
	    p->header.size = -p->header.size ;
	    Row_tuples [row] = 0 ;
	}
    }

    for (col = 0 ; col < n_col ; col++)
    {
	if (NON_PIVOTAL_COL (col) && Col_tuples [col])
	{
	    DEBUG2 (("col "ID" tuples "ID"\n", col, Col_tuples [col])) ;
	    p = Numeric->Memory + Col_tuples [col] - 1 ;
	    DEBUG2 (("Freeing tuple list col "ID", p-S "ID", size "ID"\n",
		col, (Int) (p-Numeric->Memory), p->header.size)) ;
	    ASSERT (p->header.size > 0) ;
	    ASSERT (p >= Numeric->Memory + Numeric->itail) ;
	    ASSERT (p < Numeric->Memory + Numeric->size) ;
	    p->header.size = -p->header.size ;
	    Col_tuples [col] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* mark the elements, and compress the name space */
    /* ---------------------------------------------------------------------- */

    nel = Work->nel ;
    ASSERT (nel < Work->elen) ;

#ifndef NDEBUG
    nmark = 0 ;
    UMF_dump_current_front (Numeric, Work, FALSE) ;
    DEBUGm1 (("E [0] "ID"  \n", E [0])) ;
    ASSERT (IMPLIES (E [0],
		Work->Flublock == (Entry *) (Numeric->Memory + E [0]))) ;
    ASSERT (IMPLIES (Work->Flublock,
		Work->Flublock == (Entry *) (Numeric->Memory + E [0]))) ;
    ASSERT ((E [0] != 0) == (Work->Flublock != (Entry *) NULL)) ;
#endif

    e2 = 0 ;

    for (e = 0 ; e <= nel ; e++) /* for all elements in order of creation */
    {
	if (E [e])
	{
	    psrc = Numeric->Memory + E [e] ;
	    psrc-- ;		/* get the header of this block */
	    if (e > 0)
	    {
		e2++ ;	/* do not renumber element zero */
	    }
	    ASSERT (psrc->header.size > 0) ;
	    psrc->header.size = e2  ;	/* store the new name in the header */
#ifndef NDEBUG
	    nmark++ ;
#endif
	    DEBUG7 ((ID":: Mark e "ID" at psrc-S "ID", new e "ID"\n",
		nmark, e, (Int) (psrc-Numeric->Memory), e2)) ;
	    E [e] = 0 ;
	    if (e == Work->prior_element)
	    {
		Work->prior_element = e2 ;
	    }
	}
    }

    /* all 1..e2 are now in use (element zero may or may not be in use) */
    Work->nel = e2 ;
    nel = Work->nel ;

#ifndef NDEBUG
    for (e = 0 ; e < Work->elen ; e++)
    {
	ASSERT (!E [e]) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* compress the elements */
    /* ---------------------------------------------------------------------- */

    /* point to tail marker block of size 1 + header */
    psrc = Numeric->Memory + Numeric->size - 2 ;
    pdest = psrc ;
    prevsize = psrc->header.prevsize ;
    DEBUG7 (("Starting the compression:\n")) ;

    while (prevsize > 0)
    {

	/* ------------------------------------------------------------------ */
	/* move up to the next element above the current header, and */
	/* get the element name and size */
	/* (if it is an element, the name will be positive) */
	/* ------------------------------------------------------------------ */

	size = prevsize ;
	psrc -= (size + 1) ;
	e = psrc->header.size ;
	prevsize = psrc->header.prevsize ;
	/* top block at tail has prevsize of 0 */

	/* a free block will have a negative size, so skip it */
	/* otherwise, if size >= 0, it holds the element name, not the size */

	DEBUG8 (("psrc-S: "ID" prevsize: "ID" size: "ID,
	    (Int) (psrc-Numeric->Memory), prevsize, size)) ;

	if (e == 0)
	{
	    /* -------------------------------------------------------------- */
	    /* this is the current frontal matrix */
	    /* -------------------------------------------------------------- */

	    Entry *F1, *F2, *Fsrc, *Fdst ;
	    Int c, r, k, dr, dc, gap, gap1, gap2, nb ;

	    /* shift the frontal matrix down */
	    F1 = (Entry *) (psrc + 1) ;

	    /* get the size of the current front.  r and c could be zero */
	    k = Work->fnpiv ;
	    dr = Work->fnr_curr ;
	    dc = Work->fnc_curr ;
	    r = Work->fnrows ;
	    c = Work->fncols ;
	    nb = Work->nb ;

	    ASSERT ((dr >= 0 && (dr % 2) == 1) || dr == 0) ;
	    ASSERT (drnew >= 0) ;
	    if (drnew % 2 == 0)
	    {
		/* make sure leading frontal matrix dimension is always odd */
		drnew++ ;
	    }
	    drnew = MIN (dr, drnew) ;
	    ASSERT ((drnew >= 0 && (drnew % 2) == 1) || drnew == 0) ;

	    pnext = pdest ;

#ifndef NDEBUG
	    DEBUGm2 (("move front: dr "ID" dc "ID" r "ID" drnew "ID" c "ID
		" dcnew " ID" k "ID"\n", dr, dc, r, drnew, c, dcnew, k)) ;
	    DEBUG7 (("\n")) ;
	    DEBUG7 ((ID":: Move current frontal matrix from: psrc-S: "ID" \n",
		nmark, (Int) (psrc-Numeric->Memory))) ;
	    nmark-- ;
	    ASSERT (E [e] == 0) ;
	    ASSERT (Work->Flublock == F1) ;
	    ASSERT (Work->Flblock  == Work->Flublock + nb*nb) ;
	    ASSERT (Work->Fublock  == Work->Flblock  + dr*nb) ;
	    ASSERT (Work->Fcblock  == Work->Fublock  + nb*dc) ;
	    DEBUG7 (("C  block: ")) ;
	    UMF_dump_dense (Work->Fcblock,  dr, r, c) ;
	    DEBUG7 (("L  block: ")) ;
	    UMF_dump_dense (Work->Flblock,  dr, r, k);
	    DEBUG7 (("U' block: ")) ;
	    UMF_dump_dense (Work->Fublock,  dc, c, k) ;
	    DEBUG7 (("LU block: ")) ;
	    UMF_dump_dense (Work->Flublock, nb, k, k) ;
	    ASSERT (r <= drnew && c <= dcnew && drnew <= dr && dcnew <= dc) ;
#endif

	    /* compact frontal matrix to drnew-by-dcnew before moving it */

	    /* do not compact the LU block (nb-by-nb) */

	    /* compact the columns of L (from dr-by-nb to drnew-by-nb) */
	    Fsrc = Work->Flblock ;
	    Fdst = Work->Flblock ;
	    ASSERT (Fdst == F1 + nb*nb) ;
	    gap1 = dr - r ;
	    gap2 = drnew - r ;
	    ASSERT (gap1 >= 0) ;
	    for (j = 0 ; j < k ; j++)
	    {
		for (i = 0 ; i < r ; i++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*k) ;
	    Fdst += drnew * (nb - k) ;

	    /* compact the rows of U (U' from dc-by-nb to dcnew-by-nb) */
	    Fsrc = Work->Fublock ;
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb) ;
	    gap1 = dc - c ;
	    gap2 = dcnew - c ;
	    for (i = 0 ; i < k ; i++)
	    {
		for (j = 0 ; j < c ; j++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + dcnew*k) ;
	    Fdst += dcnew * (nb - k) ;

	    /* compact the columns of C (from dr-by-dc to drnew-by-dcnew) */
	    Fsrc = Work->Fcblock ;
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + nb*dcnew) ;
	    gap1 = dr - r ;
	    gap2 = drnew - r ;
	    for (j = 0 ; j < c ; j++)
	    {
		for (i = 0 ; i < r ; i++)
		{
		    *Fdst++ = *Fsrc++ ;
		}
		Fsrc += gap1 ;
		Fdst += gap2 ;
	    }
	    ASSERT (Fdst == F1 + nb*nb + drnew*nb + nb*dcnew + drnew*c) ;

	    /* recompute Fcpos, if necessary */
	    if (do_Fcpos)
	    {
		Int *Fcols, *Fcpos ;
		Fcols = Work->Fcols ;
		Fcpos = Work->Fcpos ;
		for (j = 0 ; j < c ; j++)
		{
		    col = Fcols [j] ;
		    ASSERT (col >= 0 && col < Work->n_col) ;
		    ASSERT (Fcpos [col] == j * dr) ;
		    Fcpos [col] = j * drnew ;
		}
#ifndef NDEBUG
		{
		    Int cnt = 0 ;
		    for (j = 0 ; j < Work->n_col ; j++)
		    {
			if (Fcpos [j] != EMPTY) cnt++ ;
		    }
		    DEBUGm2 (("Recompute Fcpos cnt "ID" c "ID"\n", cnt, c)) ;
		    ASSERT (cnt == c) ;
		}
#endif
	    }

#ifndef NDEBUG
	    DEBUGm2 (("Compacted front, drnew "ID" dcnew "ID"\n", drnew, dcnew)) ;
	    DEBUG7 (("C  block: ")) ;
	    UMF_dump_dense (F1 + nb*nb + drnew*nb + nb*dcnew, drnew, r, c) ;
	    DEBUG7 (("L  block: ")) ;
	    UMF_dump_dense (F1 + nb*nb, drnew, r, k) ;
	    DEBUG7 (("U  block: ")) ;
	    UMF_dump_dense (F1 + nb*nb + drnew*nb, nb, k, c) ;
	    DEBUG7 (("LU block: ")) ;
	    UMF_dump_dense (F1, nb, k, k) ;
#endif

	    /* Compacted dimensions of the new frontal matrix. */
	    Work->fnr_curr = drnew ;
	    Work->fnc_curr = dcnew ;
	    Work->fcurr_size = (drnew + nb) * (dcnew + nb) ;
	    size = UNITS (Entry, Work->fcurr_size) ;

	    /* make sure the object doesn't evaporate.  The front can have
	     * zero size (Work->fcurr_size = 0), but the size of the memory
	     * block containing it cannot have zero size. */
	    size = MAX (1, size) ;

	    /* get the destination of frontal matrix */
	    pnext->header.prevsize = size ;
	    pdest -= (size + 1) ;
	    F2 = (Entry *) (pdest + 1) ;

	    ASSERT ((unsigned Int) psrc + 1 + size <= (unsigned Int) pnext) ;
	    ASSERT (psrc <= pdest) ;
	    ASSERT (F1 <= F2) ;

	    /* move the C block first */
	    Fsrc = F1 + nb*nb + drnew*nb + nb*dcnew + drnew*c ;
	    Fdst = F2 + nb*nb + drnew*nb + nb*dcnew + drnew*c ;
	    gap = drnew - r ;
	    for (j = c-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		/* move column j of C */
		for (i = r-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb + nb*dcnew) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb + nb*dcnew) ;

	    /* move the U block */
	    Fsrc -= dcnew * (nb - k) ;
	    Fdst -= dcnew * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb + dcnew*k) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb + dcnew*k) ;
	    gap = dcnew - c ;
	    for (i = k-1 ; i >= 0 ; i--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (j = c-1 ; j >= 0 ; j--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb + drnew*nb) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*nb) ;

	    /* move the L block */
	    Fsrc -= drnew * (nb - k) ;
	    Fdst -= drnew * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*nb + drnew*k) ;
	    ASSERT (Fdst == F2 + nb*nb + drnew*k) ;
	    gap = drnew - r ;
	    for (j = k-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (i = r-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1 + nb*nb) ;
	    ASSERT (Fdst == F2 + nb*nb) ;

	    /* move the LU block */
	    Fsrc -= nb * (nb - k) ;
	    Fdst -= nb * (nb - k) ;
	    ASSERT (Fsrc == F1 + nb*k) ;
	    ASSERT (Fdst == F2 + nb*k) ;
	    gap = nb - k ;
	    for (j = k-1 ; j >= 0 ; j--)
	    {
		Fsrc -= gap ;
		Fdst -= gap ;
		for (i = k-1 ; i >= 0 ; i--)
		{
		    *--Fdst = *--Fsrc ;
		}
	    }
	    ASSERT (Fsrc == F1) ;
	    ASSERT (Fdst == F2) ;

	    E [0] = (pdest + 1) - Numeric->Memory ;

	    Work->Flublock = (Entry *) (Numeric->Memory + E [0]) ;
	    ASSERT (Work->Flublock == F2) ;
	    Work->Flblock  = Work->Flublock + nb * nb ;
	    Work->Fublock  = Work->Flblock  + drnew * nb ;
	    Work->Fcblock  = Work->Fublock  + nb * dcnew ;

	    pdest->header.prevsize = 0 ;
	    pdest->header.size = size ;

#ifndef NDEBUG
	    DEBUG7 (("After moving compressed current frontal matrix:\n")) ;
	    DEBUG7 (("C  block: ")) ;
	    UMF_dump_dense (Work->Fcblock,  drnew, r, c) ;
	    DEBUG7 (("L  block: ")) ;
	    UMF_dump_dense (Work->Flblock,  drnew, r, k);
	    DEBUG7 (("U' block: ")) ;
	    UMF_dump_dense (Work->Fublock,  dcnew, c, k) ;
	    DEBUG7 (("LU block: ")) ;
	    UMF_dump_dense (Work->Flublock, nb, k, k) ;
#endif

	}
	else if (e > 0)
	{

	    /* -------------------------------------------------------------- */
	    /* this is an element, compress and move from psrc down to pdest */
	    /* -------------------------------------------------------------- */

#ifndef NDEBUG
	    DEBUG7 (("\n")) ;
	    DEBUG7 ((ID":: Move element "ID": from: "ID" \n",
		nmark, e, (Int) (psrc-Numeric->Memory))) ;
	    nmark-- ;
	    ASSERT (e <= nel) ;
	    ASSERT (E [e] == 0) ;
#endif

	    /* -------------------------------------------------------------- */
	    /* get the element scalars, and pointers to C, Rows, and Cols: */
	    /* -------------------------------------------------------------- */

	    p = psrc + 1 ;
	    GET_ELEMENT (epsrc, p, Cols, Rows, ncols, nrows, C) ;
	    nrowsleft = epsrc->nrowsleft ;
	    ncolsleft = epsrc->ncolsleft ;
	    cdeg = epsrc->cdeg ;
	    rdeg = epsrc->rdeg ;

#ifndef NDEBUG
	    DEBUG7 ((" nrows "ID" nrowsleft "ID"\n", nrows, nrowsleft)) ;
	    DEBUG7 ((" ncols "ID" ncolsleft "ID"\n", ncols, ncolsleft)) ;
	    DEBUG8 ((" Rows:")) ;
	    for (i = 0 ; i < nrows ; i++) DEBUG8 ((" "ID, Rows [i])) ;
	    DEBUG8 (("\n Cols:")) ;
	    for (j = 0 ; j < ncols ; j++) DEBUG8 ((" "ID, Cols [j])) ;
	    DEBUG8 (("\n")) ;
#endif

	    /* -------------------------------------------------------------- */
	    /* determine the layout of the new element */
	    /* -------------------------------------------------------------- */

	    csize = nrowsleft * ncolsleft ;
	    size2 = UNITS (Element, 1)
		  + UNITS (Int, nrowsleft + ncolsleft)
		  + UNITS (Entry, csize) ;

	    DEBUG7 (("Old size "ID" New size "ID"\n", size, size2)) ;

	    pnext = pdest ;
	    pnext->header.prevsize = size2 ;
	    pdest -= (size2 + 1) ;

	    ASSERT (size2 <= size) ;
	    ASSERT ((unsigned Int) psrc + 1 + size <= (unsigned Int) pnext) ;
	    ASSERT (psrc <= pdest) ;

	    p = pdest + 1 ;
	    epdest = (Element *) p ;
	    p += UNITS (Element, 1) ;
	    Cols2 = (Int *) p ;
	    Rows2 = Cols2 + ncolsleft ;
	    p += UNITS (Int, nrowsleft + ncolsleft) ;
	    C2 = (Entry *) p ;

	    ASSERT (epdest >= epsrc) ;
	    ASSERT (Rows2 >= Rows) ;
	    ASSERT (Cols2 >= Cols) ;
	    ASSERT (C2 >= C) ;
	    ASSERT (p + UNITS (Entry, csize) == pnext) ;

	    /* -------------------------------------------------------------- */
	    /* move the contribution block */
	    /* -------------------------------------------------------------- */

	    /* overlap = psrc + size + 1 > pdest ; */

	    if (nrowsleft < nrows || ncolsleft < ncols)
	    {

		/* ---------------------------------------------------------- */
		/* compress contribution block in place prior to moving it */
		/* ---------------------------------------------------------- */

		DEBUG7 (("Compress C in place prior to move:\n"));
#ifndef NDEBUG
		UMF_dump_dense (C, nrows, nrows, ncols) ;
#endif
		C1 = C ;
		C3 = C ;
		for (j = 0 ; j < ncols ; j++)
		{
		    if (Cols [j] >= 0)
		    {
			for (i = 0 ; i < nrows ; i++)
			{
			    if (Rows [i] >= 0)
			    {
				*C3++ = C1 [i] ;
			    }
			}
		    }
		    C1 += nrows ;
		}
		ASSERT (C3-C == csize) ;
		DEBUG8 (("Newly compressed contrib. block (all in use):\n")) ;
#ifndef NDEBUG
		UMF_dump_dense (C, nrowsleft, nrowsleft, ncolsleft) ;
#endif
	    }

	    /* shift the contribution block down */
	    C += csize ;
	    C2 += csize ;
	    for (i = 0 ; i < csize ; i++)
	    {
		*--C2 = *--C ;
	    }

	    /* -------------------------------------------------------------- */
	    /* move the row indices */
	    /* -------------------------------------------------------------- */

	    i2 = nrowsleft ;
	    for (i = nrows - 1 ; i >= 0 ; i--)
	    {
		ASSERT (Rows2+i2 >= Rows+i) ;
		if (Rows [i] >= 0)
		{
		    Rows2 [--i2] = Rows [i] ;
		}
	    }
	    ASSERT (i2 == 0) ;

	    j2 = ncolsleft ;
	    for (j = ncols - 1 ; j >= 0 ; j--)
	    {
		ASSERT (Cols2+j2 >= Cols+j) ;
		if (Cols [j] >= 0)
		{
		    Cols2 [--j2] = Cols [j] ;
		}
	    }
	    ASSERT (j2 == 0) ;

	    /* -------------------------------------------------------------- */
	    /* construct the new header */
	    /* -------------------------------------------------------------- */

	    /* E [0...e] is now valid */
	    E [e] = (pdest + 1) - Numeric->Memory ;
	    epdest = (Element *) (pdest + 1) ;

	    epdest->next = EMPTY ;	/* destroys the son list */
	    epdest->ncols = ncolsleft ;
	    epdest->nrows = nrowsleft ;
	    epdest->ncolsleft = ncolsleft ;
	    epdest->nrowsleft = nrowsleft ;
	    epdest->rdeg = rdeg ;
	    epdest->cdeg = cdeg ;

	    ASSERT (size2 <= size) ;
	    pdest->header.prevsize = 0 ;
	    pdest->header.size = size2 ;

	    DEBUG7 (("After moving it:\n")) ;
#ifndef NDEBUG
	    UMF_dump_element (Numeric, Work, e, FALSE) ;
#endif
	}

#ifndef NDEBUG
	else
	{
	    DEBUG8 ((" free\n")) ;
	}
#endif
	DEBUG7 (("psrc "ID"  tail "ID"\n",
	(Int) (psrc-Numeric->Memory), Numeric->itail)) ;
    }

    ASSERT (psrc == Numeric->Memory + Numeric->itail) ;
    ASSERT (nmark == 0) ;

    /* ---------------------------------------------------------------------- */
    /* final tail pointer */
    /* ---------------------------------------------------------------------- */

    ASSERT (pdest >= Numeric->Memory + Numeric->itail) ;
    Numeric->itail = pdest - Numeric->Memory ;
    pdest->header.prevsize = 0 ;
    Numeric->ibig = EMPTY ;
    Numeric->tail_usage = Numeric->size - Numeric->itail ;

    /* ---------------------------------------------------------------------- */
    /* clear the unused E [nel+1 .. Work->elen - 1] */
    /* ---------------------------------------------------------------------- */

    for (e = nel+1 ; e < Work->elen ; e++)
    {
	E [e] = 0 ;
    }

#ifndef NDEBUG
    UMF_dump_packed_memory (Numeric, Work) ;
#endif

    DEBUG8 (("::::GARBAGE COLLECTION DONE::::\n")) ;
}
