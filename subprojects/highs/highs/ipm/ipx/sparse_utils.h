#ifndef IPX_SPARSE_UTILS_H_
#define IPX_SPARSE_UTILS_H_

#include "ipm/ipx/ipx_internal.h"

namespace ipx {

// Depth-first search in the graph of matrix A.
//
// @istart starting node, must be unmarked on entry.
// @Ap, @Ai pattern of matrix A in CSC format.
// @colmap The neighbours of node i are the entries in col j = @colmap[i] of A.
//         If j is negative, node i has no neighbours. @colmap can be NULL, in
//         which case the identity mapping is assumed.
// @top, @istack On return @istack[newtop..@top-1] holds the reached nodes that
//               have previously been unmarked; they are marked now and newtop
//               is returned.
// @marked, @marker Node i is "marked" iff @marked[i] == @marker.
// @work workspace of size # rows of A.
//
// Nodes are written to @istack in the order in which their search finishes
// (postorder), from position @top-1 downwards. Hence @istack[newtop..@top-1]
// is a topological order of the reached nodes. @istack[0..@top-1] is used as
// scratch space during the search.
//
Int DepthFirstSearch(Int istart, const Int* Ap, const Int* Ai,
                     const Int* colmap, Int top, Int* istack, Int* marked,
                     Int marker, Int* work);

// Alternating augmenting path for extending a matching for an m-by-n matrix A.
//
// @jstart column of matrix A that should be matched.
// @Ap, @Ai pattern of matrix A in CSC format.
// @jmatch array of size m.
//         Row i is matched to column j if @jmatch[i] = j >= 0.
//         Row i is yet unmatched if @jmatch[i] == -1.
//         Row i is unmatched and excluded from being matched if jmatch[i] < -1.
// @cheap array of size n. On the first call, @cheap must hold @Ap[0..n-1].
//        It must be unchanged between subsequent calls.
// @marked array of size n. On entry @marked[j] != @jstart for all j.
//         On return some entries were set to @jstart.
// @work size m workspace.
// @work2 size m+1 workspace.
// @work3 size m+1 workspace.
//
// Returns true if the matching was extended.
//
// The search is a depth-first search for an alternating path from @jstart to
// a free row, with a "cheap assignment" check for a directly adjacent free
// row whenever a column is first visited, as in [1].
//
// [1] I. S. Duff, "On algorithms for obtaining a maximum transversal",
//     ACM Trans. Math. Software 7 (1981), pp. 315-330.
//
bool AugmentingPath(Int jstart, const Int* Ap, const Int* Ai, Int* jmatch,
                    Int* cheap, Int* marked, Int* work, Int* work2, Int* work3);

}  // namespace ipx

#endif // IPX_SPARSE_UTILS_H_
