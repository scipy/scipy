/* ========================================================================== */
/* === AMD_post_tree ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* AMD Version 1.0 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A. Davis,   */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.          */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/amd                           */
/* -------------------------------------------------------------------------- */

/* Post-ordering of a supernodal elimination tree.  */

#include "amd_internal.h"

GLOBAL Int AMD_post_tree
(
    Int root,			/* root of the tree */
    Int k,			/* start numbering at k */
    Int Child [ ],		/* input argument of size nn, undefined on
				 * output.  Child [i] is the head of a link
				 * list of all nodes that are children of node
				 * i in the tree. */
    const Int Sibling [ ],	/* input argument of size nn, not modified.
				 * If f is a node in the link list of the
				 * children of node i, then Sibling [f] is the
				 * next child of node i.
				 */
    Int Order [ ],		/* output order, of size nn.  Order [i] = k
				 * if node i is the kth node of the reordered
				 * tree. */
    Int Stack [ ]		/* workspace of size nn */
#ifndef NDEBUG
    , Int nn			/* nodes are in the range 0..nn-1. */
#endif
)
{
    Int f, head, h, i ;

#if 0
    /* ---------------------------------------------------------------------- */
    /* recursive version (Stack [ ] is not used): */
    /* ---------------------------------------------------------------------- */

    /* this is simple, but can caouse stack overflow if nn is large */
    i = root ;
    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
    {
	k = AMD_post_tree (f, k, Child, Sibling, Order, Stack, nn) ;
    }
    Order [i] = k++ ;
    return (k) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* ---------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	ASSERT (head < nn) ;
	i = Stack [head] ;
	AMD_DEBUG1 (("head of stack "ID" \n", i)) ;
	ASSERT (i >= 0 && i < nn) ;

	if (Child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		head++ ;
		ASSERT (head < nn) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    h = head ;
	    ASSERT (head < nn) ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		ASSERT (h > 0) ;
		Stack [h--] = f ;
		AMD_DEBUG1 (("push "ID" on stack\n", f)) ;
		ASSERT (f >= 0 && f < nn) ;
	    }
	    ASSERT (Stack [h] == i) ;

	    /* delete child list so that i gets ordered next time we see it */
	    Child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    AMD_DEBUG1 (("pop "ID" order "ID"\n", i, k)) ;
	    Order [i] = k++ ;
	    ASSERT (k <= nn) ;
	}

#ifndef NDEBUG
	AMD_DEBUG1 (("\nStack:")) ;
	for (h = head ; h >= 0 ; h--)
	{
	    Int j = Stack [h] ;
	    AMD_DEBUG1 ((" "ID, j)) ;
	    ASSERT (j >= 0 && j < nn) ;
	}
	AMD_DEBUG1 (("\n\n")) ;
	ASSERT (head < nn) ;
#endif

    }
    return (k) ;
}
