#include "ipm/ipx/sparse_utils.h"
#include <cassert>

namespace ipx {

Int DepthFirstSearch(Int istart, const Int* Ap, const Int* Ai,
                     const Int* colmap, Int top, Int* istack, Int* marked,
                     Int marker, Int* work) {
    assert(marked[istart] != marker);
    Int* pstack = work;
    Int head = 0;
    istack[head] = istart;
    while (head >= 0) {
        const Int i = istack[head];
        const Int j = colmap ? colmap[i] : i;
        if (marked[i] != marker) {
            marked[i] = marker;
            pstack[head] = j >= 0 ? Ap[j] : 0;
        }
        bool done = true;
        const Int pbeg = pstack[head];
        const Int pend = j >= 0 ? Ap[j+1] : 0;
        for (Int p = pbeg; p < pend; p++) {
            Int inext = Ai[p];
            if (marked[inext] == marker)
                continue;
            pstack[head] = p+1;
            istack[++head] = inext;
            done = false;
            break;
        }
        if (done) {
            head--;
            istack[--top] = i;
        }
    }
    return top;
}

bool AugmentingPath(Int jstart, const Int* Ap, const Int* Ai, Int* jmatch,
                    Int* cheap, Int* marked, Int* work, Int* work2, Int* work3){
    // istack holds row indices without duplicates.
    // jstack holds jstart and up to m indices from jmatch.
    // pstack holds the corresponding column pointers.
    Int* istack = work;
    Int* jstack = work2;
    Int* pstack = work3;
    bool found = false;
    Int head = 0;
    jstack[head] = jstart;
    while (head >= 0) {
        const Int j = jstack[head];
        if (marked[j] != jstart) { // first time j visited in this path
            marked[j] = jstart;
            Int i, p;
            for (p = cheap[j]; p < Ap[j+1] && !found; p++) {
                i = Ai[p];
                found = jmatch[i] == -1;
            }
            cheap[j] = p;
            if (found) {
                istack[head] = i;
                break;
            }
            pstack[head] = Ap[j]; // start depth-first search from j
        }
        const Int pbeg = pstack[head];
        const Int pend = Ap[j+1];
        Int p;
        for (p = pbeg; p < pend; p++) {
            Int i = Ai[p];
            if (jmatch[i] < -1) // row to be ignored
                continue;
            assert(jmatch[i] >= 0);
            if (marked[jmatch[i]] == jstart)
                continue;
            pstack[head] = p+1;
            istack[head] = i;
            jstack[++head] = jmatch[i]; // continue augmenting path here
            break;
        }
        if (p == Ap[j+1])
            head--;
    }
    if (found) {
        for (Int p = head; p >= 0; p--)
            jmatch[istack[p]] = jstack[p];
    }
    return found;
}

}  // namespace ipx
