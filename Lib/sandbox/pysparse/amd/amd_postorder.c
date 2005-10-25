/* ========================================================================== */
/* === AMD_postorder ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* Perform a postordering (via depth-first search) of an assembly tree. */

#include "amd_internal.h"

GLOBAL void AMD_postorder
(
    /* inputs, not modified on output: */
    Int nn,		/* nodes are in the range 0..nn-1 */
    Int Parent [ ],	/* Parent [j] is the parent of j, or EMPTY if root */
    Int Nv [ ],		/* Nv [j] > 0 number of pivots represented by node j,
			 * or zero if j is not a node. */
    Int Fsize [ ],	/* Fsize [j]: size of node j */

    /* output, not defined on input: */
    Int Order [ ],	/* output post-order */

    /* workspaces of size nn: */
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
)
{
    Int i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev, bigf, fnext ;

    for (j = 0 ; j < nn ; j++)
    {
	Child [j] = EMPTY ;
	Sibling [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* place the children in link lists - bigger elements tend to be last */
    /* ---------------------------------------------------------------------- */

    for (j = nn-1 ; j >= 0 ; j--)
    {
	if (Nv [j] > 0)
	{
	    /* this is an element */
	    parent = Parent [j] ;
	    if (parent != EMPTY)
	    {
		/* place the element in link list of the children its parent */
		/* bigger elements will tend to be at the end of the list */
		Sibling [j] = Child [parent] ;
		Child [parent] = j ;
	    }
	}
    }

#ifndef NDEBUG
    {
	Int nels, ff, nchild ;
	AMD_DEBUG1 (("\n\n================================ AMD_postorder:\n")) ;
	nels = 0 ;
	for (j = 0 ; j < nn ; j++)
	{
	    if (Nv [j] > 0)
	    {
		AMD_DEBUG1 (( ""ID" :  nels "ID" npiv "ID" size "ID
		    " parent "ID" maxfr "ID"\n", j, nels,
		    Nv [j], Fsize [j], Parent [j], Fsize [j])) ;
		/* this is an element */
		/* dump the link list of children */
		nchild = 0 ;
		AMD_DEBUG1 (("    Children: ")) ;
		for (ff = Child [j] ; ff != EMPTY ; ff = Sibling [ff])
		{
		    AMD_DEBUG1 ((ID" ", ff)) ;
		    ASSERT (Parent [ff] == j) ;
		    nchild++ ;
		    ASSERT (nchild < nn) ;
		}
		AMD_DEBUG1 (("\n")) ;
		parent = Parent [j] ;
		if (parent != EMPTY)
		{
		    ASSERT (Nv [parent] > 0) ;
		}
		nels++ ;
	    }
	}
    }
    AMD_DEBUG1 (("\n\nGo through the children of each node, and put\n"
		 "the biggest child last in each list:\n")) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* place the largest child last in the list of children for each node */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	if (Nv [i] > 0 && Child [i] != EMPTY)
	{

#ifndef NDEBUG
	    Int nchild ;
	    AMD_DEBUG1 (("Before partial sort, element "ID"\n", i)) ;
	    nchild = 0 ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		AMD_DEBUG1 (("      f: "ID"  size: "ID"\n", f, Fsize [f])) ;
		nchild++ ;
		ASSERT (nchild <= nn) ;
	    }
#endif

	    /* find the biggest element in the child list */
	    fprev = EMPTY ;
	    maxfrsize = EMPTY ;
	    bigfprev = EMPTY ;
	    bigf = EMPTY ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		frsize = Fsize [f] ;
		if (frsize >= maxfrsize)
		{
		    /* this is the biggest seen so far */
		    maxfrsize = frsize ;
		    bigfprev = fprev ;
		    bigf = f ;
		}
		fprev = f ;
	    }
	    ASSERT (bigf != EMPTY) ;

	    fnext = Sibling [bigf] ;

	    AMD_DEBUG1 (("bigf "ID" maxfrsize "ID" bigfprev "ID" fnext "ID
		" fprev " ID"\n", bigf, maxfrsize, bigfprev, fnext, fprev)) ;

	    if (fnext != EMPTY)
	    {
		/* if fnext is EMPTY, then bigf is already at the end of list */

		if (bigfprev == EMPTY)
		{
		    /* delete bigf from the element of the list */
		    Child [i] = fnext ;
		}
		else
		{
		    /* delete bigf from the middle of the list */
		    Sibling [bigfprev] = fnext ;
		}

		/* put bigf at the end of the list */
		Sibling [bigf] = EMPTY ;
		ASSERT (Child [i] != EMPTY) ;
		ASSERT (fprev != bigf) ;
		ASSERT (fprev != EMPTY) ;
		Sibling [fprev] = bigf ;
	    }

#ifndef NDEBUG
	    AMD_DEBUG1 (("After partial sort, element "ID"\n", i)) ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (f >= 0 && f < nn) ;
		AMD_DEBUG1 (("        "ID"  "ID"\n", f, Fsize [f])) ;
		ASSERT (Nv [f] > 0) ;
		nchild-- ;
	    }
	    ASSERT (nchild == 0) ;
#endif

	}
    }

    /* ---------------------------------------------------------------------- */
    /* postorder the assembly tree */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	Order [i] = EMPTY ;
    }

    k = 0 ;

    for (i = 0 ; i < nn ; i++)
    {
	if (Parent [i] == EMPTY && Nv [i] > 0)
	{
	    AMD_DEBUG1 (("Root of assembly tree "ID"\n", i)) ;
	    k = AMD_post_tree (i, k, Child, Sibling, Order, Stack
#ifndef NDEBUG
		, nn
#endif
		) ;
	}
    }
}
