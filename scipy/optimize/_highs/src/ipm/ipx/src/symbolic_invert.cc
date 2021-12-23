// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include "symbolic_invert.h"
#include <algorithm>
#include <cassert>
#include <random>
#include "sparse_matrix.h"
#include "sparse_utils.h"

namespace ipx {

static std::vector<Int> RandomPermute(const std::vector<Int>& basis) {
    const Int m = basis.size();
    std::default_random_engine re;
    std::uniform_int_distribution<Int> dist(0,m-1);
    std::vector<Int> permuted_basis(basis);
    for (Int k = 0; k < m; k++)
        std::swap(permuted_basis[k], permuted_basis[dist(re)]);
    return permuted_basis;
}

static std::vector<Int> Matching(const Model& model,
                                 const std::vector<Int>& basis) {
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const Int* Ap = AI.colptr();
    const Int* Ai = AI.rowidx();
    std::vector<Int> jmatch(m, -1);
    std::vector<Int> cheap(Ap, Ap+n+m);
    std::vector<Int> marked(n+m, -1);
    std::vector<Int> work(m), work2(m+1), work3(m+1);

    // Match singletons first, since then they are all cheap matches.
    for (Int j : basis) {
        if (Ap[j+1] == Ap[j]+1) {
            bool matched = AugmentingPath(j, Ap, Ai, jmatch.data(),
                                          cheap.data(), marked.data(),
                                          work.data(), work2.data(),
                                          work3.data());
            assert(matched);
            (void)(matched);
        }
    }
    for (Int j : basis) {
        if (Ap[j+1] != Ap[j]+1) {
            bool matched = AugmentingPath(j, Ap, Ai, jmatch.data(),
                                          cheap.data(), marked.data(),
                                          work.data(), work2.data(),
                                          work3.data());
            assert(matched);
            (void)(matched);
        }
    }
    return jmatch;
}

static std::vector<std::vector<Int> >
Blockperm(const SparseMatrix& AI, const std::vector<Int>& jmatch,
          const SparseMatrix& BT) {
    const Int m = AI.rows();
    const Int* Ap = AI.colptr();
    const Int* Ai = AI.rowidx();
    const Int* BTp = BT.colptr();
    const Int* BTi = BT.rowidx();

    // Put row indices into istack, working from the bottom of the stack (m-1)
    // to the top (0), so that indices whose depth-first search finished first
    // appear further down the stack.
    std::vector<Int> istack(m), marked(m), work(m);
    Int top = m;
    for (Int i = 0; i < m; i++) {
        if (marked[i] != 1)
            top = DepthFirstSearch(i, Ap, Ai, jmatch.data(), top, istack.data(),
                                   marked.data(), 1, work.data());
    }
    assert(top == 0);

    std::vector<std::vector<Int> > blocks;
    std::vector<Int> rowperm(m);
    top = m;
    for (Int i : istack) {
        if (marked[i] != 2) {
            // Node i is not part of any block yet.
            Int end = top;
            top = DepthFirstSearch(i, BTp, BTi, nullptr, top, rowperm.data(),
                                   marked.data(), 2, work.data());
            assert(top < end);
            blocks.push_back(std::vector<Int>(rowperm.begin() + top,
                                              rowperm.begin() + end));
        }
    }
    assert(top == 0);
    std::reverse(blocks.begin(), blocks.end());
    return blocks;
}

static SparseMatrix
CoarsenedGraph(const SparseMatrix& BT,
               const std::vector<std::vector<Int> >& blocks) {
    const Int m = BT.rows();
    const Int nb = blocks.size();

    std::vector<Int> map2block(m, -1);
    for (Int b = 0; b < nb; b++) {
        for (Int i : blocks[b]) {
            assert(map2block[i] == -1);
            map2block[i] = b;
        }
    }
    for (Int i = 0; i < m; i++)
        assert(map2block[i] >= 0);

    SparseMatrix C(nb, 0);
    std::vector<Int> marked(m, -1);
    for (Int k = 0; k < nb; k++) {
        for (Int i : blocks[k]) {
            for (Int p = BT.begin(i); p < BT.end(i); p++) {
                Int b = map2block[BT.index(p)];
                if (marked[b] != k) {
                    marked[b] = k;
                    C.push_back(b, 1.0);
                    assert(b >= k);
                }
            }
        }
        C.add_column();
    }
    return C;
}

void SymbolicInvert(const Model& model, const std::vector<Int>& basis,
                    Int* rowcounts, Int* colcounts) {
    const SparseMatrix& AI = model.AI();
    const Int m = AI.rows();
    assert((int)basis.size() == m);

    // jmatch is a permutation of basis such that B = AI[:,jmatch] has a
    // zero-free diagonal.
    std::vector<Int> jmatch = Matching(model, RandomPermute(basis));
    SparseMatrix BT = CopyColumns(AI, jmatch);
    BT = Transpose(BT);

    std::vector<std::vector<Int> > blocks = Blockperm(AI, jmatch, BT);
    SparseMatrix C = CoarsenedGraph(BT, blocks);

    const Int nb = blocks.size();
    std::vector<Int> stack(nb), marked(nb), work(nb);

    if (rowcounts) {
        std::vector<Int> jcount(AI.cols(), -1);
        std::fill(marked.begin(), marked.end(), -1);
        const Int* Cp = C.colptr();
        const Int* Ci = C.rowidx();
        for (Int b = 0; b < nb; b++) {
            Int top = DepthFirstSearch(b, Cp, Ci, nullptr, nb, stack.data(),
                                       marked.data(), b, work.data());
            Int nz = 0;
            for (Int t = top; t < nb; t++)
                nz += blocks[stack[t]].size();
            for (Int i : blocks[b])
                jcount[jmatch[i]] = nz;
        }
        for (Int p = 0; p < m; p++) {
            rowcounts[p] = jcount[basis[p]];
            assert(rowcounts[p] >= 0);
        }
    }

    if (colcounts) {
        C = Transpose(C);
        std::fill(marked.begin(), marked.end(), -1);
        const Int* Cp = C.colptr();
        const Int* Ci = C.rowidx();
        for (Int b = 0; b < nb; b++) {
            Int top = DepthFirstSearch(b, Cp, Ci, nullptr, nb, stack.data(),
                                       marked.data(), b, work.data());
            Int nz = 0;
            for (Int t = top; t < nb; t++)
                nz += blocks[stack[t]].size();
            for (Int i : blocks[b])
                colcounts[i] = nz;
        }
    }
}

}  // namespace ipx
