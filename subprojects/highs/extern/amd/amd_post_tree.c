//------------------------------------------------------------------------------
// AMD/Source/amd_post_tree: post-ordering of a single etree
//------------------------------------------------------------------------------

// AMD, Copyright (c) 1996-2022, Timothy A. Davis, Patrick R. Amestoy, and
// Iain S. Duff.  All Rights Reserved.
// SPDX-License-Identifier: BSD-3-clause

//------------------------------------------------------------------------------

/* Post-ordering of a supernodal elimination tree.  */

#include "amd_internal.h"

amd_int Highs_amd_post_tree
(
    amd_int root,			/* root of the tree */
    amd_int k,			/* start numbering at k */
    amd_int Child [ ],		/* input argument of size nn, undefined on
				 * output.  Child [i] is the head of a link
				 * list of all nodes that are children of node
				 * i in the tree. */
    const amd_int Sibling [ ],	/* input argument of size nn, not modified.
				 * If f is a node in the link list of the
				 * children of node i, then Sibling [f] is the
				 * next child of node i.
				 */
    amd_int Order [ ],		/* output order, of size nn.  Order [i] = k
				 * if node i is the kth node of the reordered
				 * tree. */
    amd_int Stack [ ]		/* workspace of size nn */
)
{
    amd_int f, head, h, i ;

#if 0
    /* --------------------------------------------------------------------- */
    /* recursive version (Stack [ ] is not used): */
    /* --------------------------------------------------------------------- */

    /* this is simple, but can cause stack overflow if nn is large */
    i = root ;
    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
    {
	k = Highs_amd_post_tree (f, k, Child, Sibling, Order, Stack, nn) ;
    }
    Order [i] = k++ ;
    return (k) ;
#endif

    /* --------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* --------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	
	i = Stack [head] ;
	
	

	if (Child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		head++ ;
		
		
	    }
	    h = head ;
	    
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		
		Stack [h--] = f ;
		
		
	    }
	    

	    /* delete child list so that i gets ordered next time we see it */
	    Child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    
	    Order [i] = k++ ;
	    
	}

    }
    return (k) ;
}
