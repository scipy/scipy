/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsDynamicRowMatrix.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <numeric>

HighsDynamicRowMatrix::HighsDynamicRowMatrix(HighsInt ncols) {
  AheadPos_.resize(ncols, -1);
  AheadNeg_.resize(ncols, -1);
}
/// adds a row to the matrix with the given values and returns its index
HighsInt HighsDynamicRowMatrix::addRow(HighsInt* Rindex, double* Rvalue,
                                       HighsInt Rlen, bool linkCols) {
  HighsInt rowindex;
  HighsInt start;
  HighsInt end;

  // insert the row in an existing empty space or append the values to the end
  // if no space that is large enough exists
  std::set<std::pair<HighsInt, HighsInt>>::iterator it;
  if (freespaces_.empty() || (it = freespaces_.lower_bound(std::make_pair(
                                  Rlen, HighsInt{-1}))) == freespaces_.end()) {
    start = ARindex_.size();
    end = start + Rlen;

    ARindex_.resize(end);
    ARvalue_.resize(end);
    ARrowindex_.resize(end);
    AprevPos_.resize(end, -1);
    AnextPos_.resize(end, -1);
    AprevNeg_.resize(end, -1);
    AnextNeg_.resize(end, -1);
  } else {
    std::pair<HighsInt, HighsInt> freeslot = *it;
    freespaces_.erase(it);

    start = freeslot.second;
    end = start + Rlen;
    // if the space was not completely occupied, we register the remainder of
    // it again in the priority queue
    if (freeslot.first > Rlen) {
      freespaces_.emplace(freeslot.first - Rlen, end);
    }
  }

  // register the range of values for this row with a reused or a new index
  if (deletedrows_.empty()) {
    rowindex = ARrange_.size();
    ARrange_.emplace_back(start, end);
    colsLinked.push_back(linkCols);
  } else {
    rowindex = deletedrows_.back();
    deletedrows_.pop_back();
    ARrange_[rowindex].first = start;
    ARrange_[rowindex].second = end;
    colsLinked[rowindex] = linkCols;
  }

  // now add the nonzeros in the order sorted by the index value
  for (HighsInt i = start; i != end; ++i) {
    ARindex_[i] = Rindex[i - start];
    ARvalue_[i] = Rvalue[i - start];
    ARrowindex_[i] = rowindex;
  }

  // link the row values to the columns
  if (!linkCols) return rowindex;

  for (HighsInt i = start; i != end; ++i) {
    HighsInt col = ARindex_[i];

    if (ARvalue_[i] > 0) {
      AprevPos_[i] = -1;
      HighsInt head = AheadPos_[col];
      AheadPos_[col] = i;
      AnextPos_[i] = head;

      if (head != -1) AprevPos_[head] = i;
    } else {
      AprevNeg_[i] = -1;
      HighsInt head = AheadNeg_[col];
      AheadNeg_[col] = i;
      AnextNeg_[i] = head;

      if (head != -1) AprevNeg_[head] = i;
    }
  }

  return rowindex;
}

void HighsDynamicRowMatrix::unlinkColumns(HighsInt rowindex) {
  if (!colsLinked[rowindex]) return;

  colsLinked[rowindex] = false;
  HighsInt start = ARrange_[rowindex].first;
  HighsInt end = ARrange_[rowindex].second;
  for (HighsInt i = start; i != end; ++i) {
    HighsInt col = ARindex_[i];

    if (ARvalue_[i] > 0) {
      HighsInt prev = AprevPos_[i];
      HighsInt next = AnextPos_[i];

      if (next != -1) {
        assert(AprevPos_[next] == i);
        AprevPos_[next] = prev;
      }

      if (prev != -1) {
        assert(AnextPos_[prev] == i);
        AnextPos_[prev] = next;
      } else {
        assert(AheadPos_[col] == i);
        AheadPos_[col] = next;
      }
    } else {
      HighsInt prev = AprevNeg_[i];
      HighsInt next = AnextNeg_[i];

      if (next != -1) {
        assert(AprevNeg_[next] == i);
        AprevNeg_[next] = prev;
      }

      if (prev != -1) {
        assert(AnextNeg_[prev] == i);
        AnextNeg_[prev] = next;
      } else {
        assert(AheadNeg_[col] == i);
        AheadNeg_[col] = next;
      }
    }
  }
}

/// removes the row with the given index from the matrix, afterwards the index
/// can be reused for new rows
void HighsDynamicRowMatrix::removeRow(HighsInt rowindex) {
  HighsInt start = ARrange_[rowindex].first;
  HighsInt end = ARrange_[rowindex].second;

  if (colsLinked[rowindex]) {
    for (HighsInt i = start; i != end; ++i) {
      HighsInt col = ARindex_[i];

      if (ARvalue_[i] > 0) {
        HighsInt prev = AprevPos_[i];
        HighsInt next = AnextPos_[i];

        if (next != -1) {
          assert(AprevPos_[next] == i);
          AprevPos_[next] = prev;
        }

        if (prev != -1) {
          assert(AnextPos_[prev] == i);
          AnextPos_[prev] = next;
        } else {
          assert(AheadPos_[col] == i);
          AheadPos_[col] = next;
        }
      } else {
        HighsInt prev = AprevNeg_[i];
        HighsInt next = AnextNeg_[i];

        if (next != -1) {
          assert(AprevNeg_[next] == i);
          AprevNeg_[next] = prev;
        }

        if (prev != -1) {
          assert(AnextNeg_[prev] == i);
          AnextNeg_[prev] = next;
        } else {
          assert(AheadNeg_[col] == i);
          AheadNeg_[col] = next;
        }
      }
    }
  }

  // register the space of the deleted row and the index so that it can be
  // reused
  deletedrows_.push_back(rowindex);
  freespaces_.emplace(end - start, start);

  // set the range to -1,-1 to indicate a deleted row
  ARrange_[rowindex].first = -1;
  ARrange_[rowindex].second = -1;
}
