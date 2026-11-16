/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsGFkLU.h
 * @brief linear system solve in GF(k) for mod-k cut separation
 */

#ifndef HIGHS_GFk_SOLVE_H_
#define HIGHS_GFk_SOLVE_H_

#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <tuple>
#include <vector>

#include "lp_data/HConst.h"

// helper struct to compute the multiplicative inverse by using fermats
// theorem and recursive repeated squaring.
// Under the assumption that k is a small prime and an 32bit HighsInt is enough
// to hold the number (k-1)^(k-2) good compilers should be able to optimize this
// code by inlining and unfolding the recursion to code that uses the fewest
// amount of integer multiplications and optimize the single integer division
// away as k is a compile time constant. Since the class was developed for the
// purpose of separating maximally violated mod-k cuts the assumption that k is
// a small constant prime number is not restrictive.

template <HighsInt k>
struct HighsGFk;

template <>
struct HighsGFk<2> {
  static constexpr unsigned int powk(unsigned int a) { return a * a; }
  static constexpr unsigned int inverse(unsigned int) { return 1; }
};

template <>
struct HighsGFk<3> {
  static constexpr unsigned int powk(unsigned int a) { return a * a * a; }
  static constexpr unsigned int inverse(unsigned int a) { return a; }
};

template <HighsInt k>
struct HighsGFk {
  static constexpr unsigned int powk(unsigned int a) {
    return (k & 1) == 0 ? HighsGFk<2>::powk(HighsGFk<k / 2>::powk(a))
                        : HighsGFk<k - 1>::powk(a) * a;
  }

  static unsigned int inverse(unsigned int a) {
    return HighsGFk<k - 2>::powk(a) % k;
  }
};

class HighsGFkSolve {
  HighsInt numCol;
  HighsInt numRow;

  // triplet format
  std::vector<HighsInt> Arow;
  std::vector<HighsInt> Acol;
  std::vector<unsigned int> Avalue;

  // sizes of rows and columns
  std::vector<HighsInt> rowsize;
  std::vector<HighsInt> colsize;

  // linked list links for column based links for each nonzero
  std::vector<HighsInt> colhead;
  std::vector<HighsInt> Anext;
  std::vector<HighsInt> Aprev;

  // splay tree links for row based iteration and lookup
  std::vector<HighsInt> rowroot;
  std::vector<HighsInt> ARleft;
  std::vector<HighsInt> ARright;

  // right hand side vector
  std::vector<unsigned int> rhs;

  // column permutation for the factorization required for backwards solve
  std::vector<HighsInt> factorColPerm;
  std::vector<HighsInt> factorRowPerm;
  std::vector<int8_t> colBasisStatus;
  std::vector<HighsBool> rowUsed;

  // working memory
  std::vector<HighsInt> iterstack;
  std::vector<HighsInt> rowpositions;
  std::vector<HighsInt> rowposColsizes;

  // priority queue to reuse free slots
  std::priority_queue<HighsInt, std::vector<HighsInt>, std::greater<HighsInt>>
      freeslots;

  void link(HighsInt pos);

  void unlink(HighsInt pos);

  void dropIfZero(HighsInt pos) {
    if (Avalue[pos] == 0) unlink(pos);
  }

  void storeRowPositions(HighsInt pos);

  void addNonzero(HighsInt row, HighsInt col, unsigned int val);

 public:
  struct SolutionEntry {
    HighsInt index;
    HighsUInt weight;

    bool operator<(const SolutionEntry& other) const {
      return index < other.index;
    }
  };

  // access to triplets and find nonzero function for unit test
  const std::vector<HighsInt>& getArow() const { return Arow; }
  const std::vector<HighsInt>& getAcol() const { return Acol; }
  const std::vector<unsigned>& getAvalue() const { return Avalue; }
  HighsInt numNonzeros() const { return int(Avalue.size() - freeslots.size()); }
  HighsInt findNonzero(HighsInt row, HighsInt col);

  template <unsigned int k, int kNumRhs = 1, typename T>
  void fromCSC(const std::vector<T>& Aval, const std::vector<HighsInt>& Aindex,
               const std::vector<HighsInt>& Astart, HighsInt numRow) {
    Avalue.clear();
    Acol.clear();
    Arow.clear();

    freeslots = decltype(freeslots)();

    numCol = Astart.size() - 1;
    this->numRow = numRow;

    colhead.assign(numCol, -1);
    colsize.assign(numCol, 0);

    rhs.assign(kNumRhs * numRow, 0);
    rowroot.assign(numRow, -1);
    rowsize.assign(numRow, 0);

    Avalue.reserve(Aval.size());
    Acol.reserve(Aval.size());
    Arow.reserve(Aval.size());

    for (HighsInt i = 0; i != numCol; ++i) {
      for (HighsInt j = Astart[i]; j != Astart[i + 1]; ++j) {
        assert(Aval[j] == (int64_t)Aval[j]);
        int64_t val = ((int64_t)Aval[j]) % k;
        if (val == 0) continue;

        if (val < 0) val += k;
        assert(val >= 0);

        Avalue.push_back(val);
        Acol.push_back(i);
        Arow.push_back(Aindex[j]);
      }
    }

    HighsInt nnz = Avalue.size();
    Anext.resize(nnz);
    Aprev.resize(nnz);
    ARleft.resize(nnz);
    ARright.resize(nnz);
    for (HighsInt pos = 0; pos != nnz; ++pos) link(pos);
  }

  template <unsigned int k, int kNumRhs = 1, typename T>
  void setRhs(HighsInt row, T val, int rhsIndex = 0) {
    rhs[kNumRhs * row + rhsIndex] = ((unsigned int)std::abs(val)) % k;
  }

  template <unsigned int k, int kNumRhs = 1, typename ReportSolution>
  void solve(ReportSolution&& reportSolution) {
    auto cmpPrio = [](const std::pair<HighsInt, HighsInt>& a,
                      const std::pair<HighsInt, HighsInt>& b) {
      return a.first > b.first;
    };
    std::priority_queue<std::pair<HighsInt, HighsInt>,
                        std::vector<std::pair<HighsInt, HighsInt>>,
                        decltype(cmpPrio)>
        pqueue(cmpPrio);

    for (HighsInt i = 0; i != numCol; ++i) pqueue.emplace(colsize[i], i);

    HighsInt maxPivot = std::min(numRow, numCol);
    factorColPerm.clear();
    factorRowPerm.clear();
    factorColPerm.reserve(maxPivot);
    factorRowPerm.reserve(maxPivot);
    colBasisStatus.assign(numCol, 0);
    rowUsed.assign(numRow, false);
    HighsInt numPivot = 0;

    while (!pqueue.empty()) {
      HighsInt pivotCol;
      HighsInt oldColSize;

      std::tie(oldColSize, pivotCol) = pqueue.top();
      pqueue.pop();

      if (colsize[pivotCol] == 0) continue;

      assert(colBasisStatus[pivotCol] == 0);

      if (colsize[pivotCol] != oldColSize) {
        pqueue.emplace(colsize[pivotCol], pivotCol);
        continue;
      }

      HighsInt pivot = -1;
      HighsInt pivotRow = -1;
      HighsInt pivotRowLen = kHighsIInf;
      for (HighsInt coliter = colhead[pivotCol]; coliter != -1;
           coliter = Anext[coliter]) {
        HighsInt row = Arow[coliter];
        if (rowUsed[row]) continue;
        if (rowsize[row] < pivotRowLen) {
          pivotRowLen = rowsize[row];
          pivotRow = row;
          pivot = coliter;
        }
      }

      assert(pivot != -1);
      assert(Acol[pivot] == pivotCol);
      assert(Arow[pivot] == pivotRow);
      assert(Avalue[pivot] > 0);
      assert(Avalue[pivot] < k);

      unsigned int pivotInverse = HighsGFk<k>::inverse(Avalue[pivot]);
      assert((Avalue[pivot] * pivotInverse) % k == 1);

      rowpositions.clear();
      rowposColsizes.clear();
      storeRowPositions(rowroot[pivotRow]);
      assert(pivotRowLen == (HighsInt)rowpositions.size());
      HighsInt next;
      for (HighsInt coliter = colhead[pivotCol]; coliter != -1;
           coliter = next) {
        next = Anext[coliter];
        if (coliter == pivot) continue;

        assert(Acol[coliter] == pivotCol);
        assert(Arow[coliter] != pivotRow);
        assert(Avalue[coliter] != 0);

        HighsInt row = Arow[coliter];
        if (rowUsed[row]) continue;

        unsigned int pivotRowScale = pivotInverse * (k - Avalue[coliter]);

        for (HighsInt i = 0; i < kNumRhs; ++i)
          rhs[kNumRhs * row + i] =
              (rhs[kNumRhs * row + i] +
               pivotRowScale * rhs[kNumRhs * pivotRow + i]) %
              k;

        for (HighsInt pivotRowPos : rowpositions) {
          HighsInt nonzeroPos = findNonzero(Arow[coliter], Acol[pivotRowPos]);

          if (nonzeroPos == -1) {
            assert(Acol[pivotRowPos] != pivotCol);
            unsigned int val = (pivotRowScale * Avalue[pivotRowPos]) % k;
            if (val != 0) addNonzero(row, Acol[pivotRowPos], val);
          } else {
            Avalue[nonzeroPos] =
                (Avalue[nonzeroPos] + pivotRowScale * Avalue[pivotRowPos]) % k;
            assert(Acol[pivotRowPos] != pivotCol || Avalue[nonzeroPos] == 0);
            dropIfZero(nonzeroPos);
          }
        }
      }

      ++numPivot;
      factorColPerm.push_back(pivotCol);
      factorRowPerm.push_back(pivotRow);
      colBasisStatus[pivotCol] = 1;
      rowUsed[pivotRow] = true;
      if (numPivot == maxPivot) break;

      for (HighsInt i = 0; i != pivotRowLen; ++i) {
        assert(Arow[rowpositions[i]] == pivotRow);
        HighsInt col = Acol[rowpositions[i]];
        HighsInt oldsize = rowposColsizes[i];

        // we only want to count rows that are not used so far, so we need to
        // decrease the counter for all columns in the pivot row by one
        --colsize[col];

        // the column size should never get negative
        assert(colsize[col] >= 0);

        if (colsize[col] == 0) continue;

        // the pivot column should occur in zero unused rows
        assert(col != pivotCol);

        // only reinsert the column if the size is smaller, as it otherwise
        // either is already at the correct position in the queue, or is at a
        // too good position and will be lazily reinserted when it is extracted
        // from the queue and the size does not match
        if (colsize[col] < oldsize) pqueue.emplace(colsize[col], col);
      }
    }

    // check if a solution exists by scanning the linearly dependent rows for
    // nonzero right hand sides
    std::array<bool, kNumRhs> hasSolution;
    HighsInt numRhsWithSolution = 0;
    for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex) {
      hasSolution[rhsIndex] = true;
      for (HighsInt i = 0; i != numRow; ++i) {
        // if the row was used it is linearly independent
        if (rowUsed[i]) continue;

        // if the row is linearly dependent, the right hand side must be zero,
        // otherwise no solution exists
        if (rhs[kNumRhs * i + rhsIndex] != 0) {
          hasSolution[rhsIndex] = false;
          break;
        }
      }

      numRhsWithSolution += hasSolution[rhsIndex];
    }

    if (numRhsWithSolution == 0) return;

    // now iterate a subset of the basic solutions.
    // When a column leaves the basis we do not allow it to enter again so that
    // we iterate at most one solution for each nonbasic column
    std::array<std::vector<SolutionEntry>, kNumRhs> solution;
    for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex)
      if (hasSolution[rhsIndex]) solution[rhsIndex].reserve(numCol);

    HighsInt numFactorRows = factorRowPerm.size();

    // create vector for swapping different columns into the basis
    // For each column we want to iterate one basic solution where the
    // column is basic
    std::vector<std::pair<HighsInt, HighsInt>> basisSwaps;
    assert(iterstack.empty());
    for (HighsInt i = numFactorRows - 1; i >= 0; --i) {
      HighsInt row = factorRowPerm[i];
      iterstack.push_back(rowroot[row]);

      while (!iterstack.empty()) {
        HighsInt rowpos = iterstack.back();
        iterstack.pop_back();
        assert(rowpos != -1);

        if (ARleft[rowpos] != -1) iterstack.push_back(ARleft[rowpos]);
        if (ARright[rowpos] != -1) iterstack.push_back(ARright[rowpos]);

        HighsInt col = Acol[rowpos];
        if (colBasisStatus[col] != 0) continue;

        colBasisStatus[col] = -1;
        basisSwaps.emplace_back(i, col);
      }
    }

    HighsInt basisSwapPos = 0;

    bool performedBasisSwap;
    do {
      performedBasisSwap = false;

      for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex)
        solution[rhsIndex].clear();

      for (HighsInt i = numFactorRows - 1; i >= 0; --i) {
        HighsInt row = factorRowPerm[i];

        std::array<unsigned int, kNumRhs> solval;

        for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex) {
          if (!hasSolution[rhsIndex]) continue;
          solval[rhsIndex] = 0;

          for (const SolutionEntry& solentry : solution[rhsIndex]) {
            HighsInt pos = findNonzero(row, solentry.index);
            if (pos != -1) solval[rhsIndex] += Avalue[pos] * solentry.weight;
          }

          solval[rhsIndex] =
              rhs[kNumRhs * row + rhsIndex] + k - (solval[rhsIndex] % k);
        }

        HighsInt col = factorColPerm[i];
        HighsInt pos = findNonzero(row, col);
        assert(pos != -1);
        unsigned int colValInverse = HighsGFk<k>::inverse(Avalue[pos]);

        for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex) {
          if (!hasSolution[rhsIndex]) continue;
          assert(solval[rhsIndex] >= 0);
          assert(colValInverse != 0);

          solval[rhsIndex] = (solval[rhsIndex] * colValInverse) % k;

          assert(solval[rhsIndex] >= 0 && solval[rhsIndex] < k);

          // only record nonzero solution values
          if (solval[rhsIndex] != 0)
            solution[rhsIndex].emplace_back(
                SolutionEntry{col, solval[rhsIndex]});
        }
      }

      for (int rhsIndex = 0; rhsIndex < kNumRhs; ++rhsIndex)
        if (hasSolution[rhsIndex]) reportSolution(solution[rhsIndex], rhsIndex);

      if (basisSwapPos < (HighsInt)basisSwaps.size()) {
        HighsInt basisIndex = basisSwaps[basisSwapPos].first;
        HighsInt enteringCol = basisSwaps[basisSwapPos].second;
        HighsInt leavingCol = factorColPerm[basisIndex];
        assert(colBasisStatus[leavingCol] == 1);
        factorColPerm[basisIndex] = enteringCol;
        colBasisStatus[enteringCol] = 1;
        colBasisStatus[leavingCol] = 0;
        performedBasisSwap = true;
        ++basisSwapPos;
      }
    } while (performedBasisSwap);
  }
};

#endif
