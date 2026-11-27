#include "ipm/ipx/sparse_utils.h"
#include <cassert>

namespace ipx {

// Iterative depth-first search. The recursion is simulated with two parallel
// stacks indexed by depth d:
//   work[d]   node being expanded at depth d,
//   istack[d] position in Ai of the next neighbour of work[d] to examine.
// Finished nodes are written to istack from position top-1 downwards, so the
// output occupies istack[top..] while the cursor stack grows from istack[0].
// At any time (depth + #finished) <= #newly reached nodes <= top, so the two
// regions of istack never overlap.
Int DepthFirstSearch(Int istart, const Int* Ap, const Int* Ai,
                     const Int* colmap, Int top, Int* istack, Int* marked,
                     Int marker, Int* work) {
    assert(marked[istart] != marker);

    // Range [begin,end) of neighbours of node i in Ai. Empty if i has none.
    auto neighbours = [&](Int i, Int& begin, Int& end) {
        const Int j = colmap ? colmap[i] : i;
        if (j < 0) {
            begin = end = 0;
        } else {
            begin = Ap[j];
            end = Ap[j+1];
        }
    };

    Int depth = 0;
    Int begin, end;
    marked[istart] = marker;
    work[0] = istart;
    neighbours(istart, begin, end);
    istack[0] = begin;

    while (depth >= 0) {
        const Int node = work[depth];
        neighbours(node, begin, end);
        Int pos = istack[depth];

        // Advance to the first neighbour that has not been reached yet.
        while (pos < end && marked[Ai[pos]] == marker)
            pos++;

        if (pos < end) {
            // Descend into an unvisited neighbour. Remember where to resume.
            const Int child = Ai[pos];
            istack[depth] = pos + 1;
            depth++;
            assert(depth < top);
            marked[child] = marker;
            work[depth] = child;
            neighbours(child, begin, end);
            istack[depth] = begin;
        } else {
            // All neighbours done: node is finished.
            depth--;
            istack[--top] = node;
        }
    }
    return top;
}

// Searches for an alternating path that starts at column jstart and ends at a
// free row; on success the matching is flipped along the path.
//
// The search is a depth-first search over columns. Column c is entered from a
// row r that is currently matched to c; the path then continues from c via
// one of its other rows. Stacks indexed by depth d:
//   work2[d]  column at depth d (work2[0] == jstart),
//   work3[d]  next position in Ai to examine for column work2[d],
//   work[d]   row used to step from column work2[d] to column work2[d+1].
// Each column on the path other than jstart is matched to a distinct row, so
// the depth is at most m, which fits the sizes m+1 and m.
//
// When a column is entered for the first time it is checked for a directly
// adjacent free row. Since a matched row never becomes free again, that check
// resumes where it stopped last time (cheap[c]), giving amortised linear cost
// over all calls.
bool AugmentingPath(Int jstart, const Int* Ap, const Int* Ai, Int* jmatch,
                    Int* cheap, Int* marked, Int* work, Int* work2,
                    Int* work3) {
    Int* col_stack = work2;
    Int* pos_stack = work3;
    Int* row_stack = work;

    // Looks for a free row in column c starting from cheap[c]. Returns the row
    // or -1.
    auto take_free_row = [&](Int c) -> Int {
        Int p = cheap[c];
        const Int end = Ap[c+1];
        while (p < end && jmatch[Ai[p]] != -1)
            p++;
        if (p == end) {
            cheap[c] = end;
            return -1;
        }
        // The row found here becomes matched, so the next scan of column c
        // can start after it.
        cheap[c] = p + 1;
        return Ai[p];
    };

    Int depth = 0;
    Int free_row = -1;
    assert(marked[jstart] != jstart);
    marked[jstart] = jstart;
    col_stack[0] = jstart;
    pos_stack[0] = Ap[jstart];
    free_row = take_free_row(jstart);

    while (free_row < 0 && depth >= 0) {
        const Int c = col_stack[depth];
        const Int end = Ap[c+1];
        Int p = pos_stack[depth];
        Int next_col = -1;

        // Find a row of column c whose partner column is not yet on the
        // search tree.
        for (; p < end; p++) {
            const Int partner = jmatch[Ai[p]];
            if (partner >= 0 && marked[partner] != jstart) {
                next_col = partner;
                break;
            }
        }

        if (next_col < 0) {
            depth--;  // column c is exhausted; backtrack
            continue;
        }

        row_stack[depth] = Ai[p];
        pos_stack[depth] = p + 1;
        depth++;
        marked[next_col] = jstart;
        col_stack[depth] = next_col;
        pos_stack[depth] = Ap[next_col];
        free_row = take_free_row(next_col);
    }

    if (free_row < 0)
        return false;

    // Flip the matching along the path: the free row goes to the deepest
    // column, and each row used to step down goes to the column above it.
    jmatch[free_row] = col_stack[depth];
    for (Int d = depth - 1; d >= 0; d--)
        jmatch[row_stack[d]] = col_stack[d];
    return true;
}

}  // namespace ipx
