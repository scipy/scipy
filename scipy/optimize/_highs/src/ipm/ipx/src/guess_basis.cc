// Copyright (c) 2018 ERGO-Code. See license.txt for license.

#include <algorithm>
#include "guess_basis.h"
#include <cassert>
#include <cmath>
#include "sparse_utils.h"
#include "utils.h"

namespace ipx {

// Computes the pattern of L\AI[:,j] in topological order in pattern[top..m-1].
static Int ComputePattern(const SparseMatrix& L, const SparseMatrix& AI,
                          Int j, const Int* rownumber, Int* pattern,
                          Int* pstack, Int* marked, Int marker) {
    const Int m = L.rows();
    const Int* Lp = L.colptr();
    const Int* Li = L.rowidx();
    Int top = m;
    for (Int p = AI.begin(j); p < AI.end(j); p++) {
        Int i = AI.index(p);
        if (marked[i] != marker)
            top = DepthFirstSearch(i, Lp, Li, rownumber, top, pattern, marked,
                                   marker, pstack);
    }
    return top;
}

// Computes the values of L\AI[:,j] scattered into lhs. Positions that are not
// in the pattern are not accessed. Returns the position of the maximum entry
// in lhs that has not been pivotal yet (-1 if no such nonzero).
static Int ComputeValues(const SparseMatrix& L, const SparseMatrix& AI,
                         Int j, const Int* rownumber, const Int* pattern,
                         Int top, Vector& lhs){
    const Int m = L.rows();
    const Int* Lp = L.colptr();
    const Int* Li = L.rowidx();
    const double* Lx = L.values();
    for (Int t = top; t < m; t++)
        lhs[pattern[t]] = 0.0;
    for (Int p = AI.begin(j); p < AI.end(j); p++) // scatter RHS into lhs
        lhs[AI.index(p)] = AI.value(p);
    double lhsmax = 0.0;        // maximum entry that has not been pivotal
    Int imax = -1;              // corresponding position in lhs
    for (Int t = top; t < m; t++) {
        Int i = pattern[t];
        double temp = lhs[i];
        Int k = rownumber[i];
        if (temp != 0.0) {
            if (k >= 0) {
                // Row i has been pivotal; update lhs.
                for (Int p = Lp[k]; p < Lp[k+1]; p++)
                    lhs[Li[p]] -= Lx[p] * temp;
            } else if (std::abs(temp) > lhsmax) {
                // Row i has not been pivotal; update lhsmax.
                lhsmax = std::abs(temp);
                imax = i;
            }
        }
    }
    return imax;
}

// Performs incomplete left-looking LU factorization of columns with infinite
// weight. If a candidate column has a sufficiently large entry after being
// updated by the columns in L, it is added to the basis. The row with maximum
// entry is chosen as pivot row.
static void ProcessFreeColumns(const Control& control, const Model& model,
                               const double* weights, std::vector<Int>* basis,
                               Int* rownumber, int* active) {
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();

    // Workspace and output of left-looking incomplete LU factorization.
    // L is stored without the unit diagonal, row indices are (unpermuted) row
    // indices from AI. The upper triangular factor is not stored.
    std::vector<Int> pattern(m), work(m), marked(m, -1);
    Vector lhs(m);
    SparseMatrix L(m,0);

    Int num_free = 0;
    for (Int j = 0; j < n+m; j++) {
        if (weights[j] != INFINITY)
            continue;
        Int top = ComputePattern(L, AI, j, rownumber, pattern.data(),
                                 work.data(), marked.data(), j);
        Int imax = ComputeValues(L, AI, j, rownumber, pattern.data(), top, lhs);
        double pivot = imax >= 0 ? lhs[imax] : 0.0;
        if (std::abs(pivot) > kLuDependencyTol) {
            assert(rownumber[imax] == -1);
            rownumber[imax] = basis->size();
            basis->push_back(j);
            // Add new column to L. Discard off-diagonal entries that were not
            // in the pattern of AI[:,j].
            for (Int t = AI.begin(j); t < AI.end(j); t++) {
                Int i = AI.index(t);
                if (rownumber[i] < 0 && lhs[i] != 0.0)
                    L.push_back(i, lhs[i]/pivot);
            }
            L.add_column();
            num_free++;
        }
        active[j] = false;
    }
    control.Debug()
        << Textline("Number of free variables in starting basis:")
        << num_free << '\n';
}

// Adds singleton columns to the basis if their entry is maximum up to a factor
// 2.0 in its row.
static void ProcessSingletons(const Control& control, const Model& model,
                              const double* weights, std::vector<Int>* basis,
                              Int* rownumber, int* active) {
    const Int m = model.rows();
    const SparseMatrix& AI = model.AI();
    const SparseMatrix& AT = model.AIt();
    Int num_singletons = 0;
    for (Int i = 0; i < m; i++) {
        if (rownumber[i] >= 0)
            continue;
        double rowmax = 0.0;
        double max_singleton = 0.0;
        Int jsingleton = -1;
        for (Int p = AT.begin(i); p < AT.end(i); p++) {
            Int j = AT.index(p);
            if (!active[j])
                continue;
            double a = std::abs(AT.value(p)) * weights[j];
            rowmax = std::max(rowmax, a);
            if (a > max_singleton && AI.end(j) == AI.begin(j)+1) {
                max_singleton = a;
                jsingleton = j;
            }
        }
        if (max_singleton > 0.0 && max_singleton >= 0.5*rowmax) {
            rownumber[i] = basis->size();
            basis->push_back(jsingleton);
            active[jsingleton] = false;
            num_singletons++;
        }
    }
    control.Debug()
        << Textline("Number of singletons in starting basis:")
        << num_singletons << '\n';
}

// Computes a matching of active columns with rows that have not been pivotal.
// Each matched column is added to the basis. Columns are processed in
// decreasing order of their weight.
static void ProcessRemaining(const Control& control, const Model& model,
                             const double* weights, std::vector<Int>* basis,
                             Int* rownumber, int* active) {
    const Int m = model.rows();
    const Int n = model.cols();
    const SparseMatrix& AI = model.AI();
    const std::vector<Int> colperm = Sortperm(n+m, weights, true);

    // jmatch[i] == -1  if row i is unmatched but eligible for matching
    // jmatch[i] == -2  if row i is excluded from being matched
    // jmatch[i] == j >= 0 if row i is matched to column j
    std::vector<Int> jmatch(m, -1);
    for (Int i = 0; i < m; i++) {
        if (rownumber[i] >= 0)
            jmatch[i] = -2;
    }
    std::vector<Int> marked(n+m, -1);
    const Int* Ap = AI.colptr();
    const Int* Ai = AI.rowidx();
    std::vector<Int> cheap(Ap, Ap+n+m);
    std::vector<Int> work(m), work2(m+1), work3(m+1);
    Int num_matched = 0;
    Int num_failed = 0;

    for (Int j : colperm) {
        if (!active[j])
            continue;
        if (weights[j] == 0.0)
            break;
        bool matched = AugmentingPath(j, Ap, Ai, jmatch.data(), cheap.data(),
                                      marked.data(), work.data(), work2.data(),
                                      work3.data());
        if (matched) {
            basis->push_back(j);
            num_matched++;
        } else {
            num_failed++;
        }
        if (num_failed >= 10*(m - (Int) basis->size()))
            break;
    }
    for (Int i = 0; i < m; i++) {
        if (jmatch[i] >= 0) {
            // rownumber[i] = m indicates that row i was matched
            assert(rownumber[i] < 0);
            rownumber[i] = m;
        }
    }
    control.Debug()
        << Textline("Number of other columns matched:")
        << num_matched << '\n'
        << Textline("Number of other columns failed:")
        << num_failed << '\n';
}

std::vector<Int> GuessBasis(const Control& control, const Model& model,
                            const double* colweights) {
    const Int m = model.rows();
    const Int n = model.cols();

    // basis starts empty and is filled one index at a time. rownumber[i] >= 0
    // iff row i was pivot row when a column was added to the basis. (The
    // specific value has a different meaning in each method.) A column is
    // "active" if it is eligible for being added to the basis.
    std::vector<Int> basis, rownumber(m, -1);
    std::vector<int> active(n+m, 1);

    ProcessFreeColumns(control, model, colweights, &basis, rownumber.data(),
                       active.data());
    ProcessSingletons(control, model, colweights, &basis, rownumber.data(),
                      active.data());
    ProcessRemaining(control, model, colweights, &basis, rownumber.data(),
                     active.data());

    // Complete basis with unit columns.
    for (Int i = 0; i < m; i++) {
        if (rownumber[i] < 0)
            basis.push_back(n+i);
    }
    assert((int)basis.size() == m);
    return basis;
}

}  // namespace ipx
