/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_DYNAMIC_ROW_MATRIX_H_
#define HIGHS_DYNAMIC_ROW_MATRIX_H_

#include <set>
#include <utility>
#include <vector>

#include "util/HighsType.h"

class HighsDynamicRowMatrix {
 private:
  /// vector of index ranges in the index and value arrays of AR for each row
  std::vector<std::pair<HighsInt, HighsInt>> ARrange_;

  /// column indices for each nonzero in AR
  std::vector<HighsInt> ARindex_;
  /// values for each nonzero in AR
  std::vector<double> ARvalue_;

  std::vector<HighsInt> ARrowindex_;
  std::vector<HighsInt> AnextPos_;
  std::vector<HighsInt> AprevPos_;
  std::vector<HighsInt> AnextNeg_;
  std::vector<HighsInt> AprevNeg_;

  /// vector of pointers to the head/tail of the nonzero block list for each
  /// column
  std::vector<HighsInt> AheadPos_;
  std::vector<HighsInt> AheadNeg_;

  std::vector<HighsBool> colsLinked;

  /// vector of column sizes

  /// keep an ordered set of free spaces in the row arrays so that they can be
  /// reused efficiently
  std::set<std::pair<HighsInt, HighsInt>> freespaces_;

  /// vector of deleted rows so that their indices can be reused
  std::vector<HighsInt> deletedrows_;

 public:
  HighsDynamicRowMatrix(HighsInt ncols);

  bool columnsLinked(HighsInt rowindex) const { return colsLinked[rowindex]; }

  void unlinkColumns(HighsInt rowindex);

  /// adds a row to the matrix with the given values and returns its index
  HighsInt addRow(HighsInt* Rindex, double* Rvalue, HighsInt Rlen,
                  bool linkCols = true);

  /// removes the row with the given index from the matrix, afterwards the index
  /// can be reused for new rows
  void removeRow(HighsInt rowindex);

  std::size_t nonzeroCapacity() const { return ARvalue_.size(); }

  /// calls the given function object for each entry in the given column.
  /// The function object should accept the row index as first argument and
  /// the nonzero value of the column in that row as the second argument.
  template <typename Func>
  void forEachPositiveColumnEntry(HighsInt col, Func&& f) const {
    HighsInt iter = AheadPos_[col];

    while (iter != -1) {
      if (!f(ARrowindex_[iter], ARvalue_[iter])) break;
      iter = AnextPos_[iter];
    }
  }

  template <typename Func>
  void forEachNegativeColumnEntry(HighsInt col, Func&& f) const {
    HighsInt iter = AheadNeg_[col];

    while (iter != -1) {
      if (!f(ARrowindex_[iter], ARvalue_[iter])) break;
      iter = AnextNeg_[iter];
    }
  }

  HighsInt getNumRows() const { return ARrange_.size(); }

  HighsInt getNumDelRows() const { return deletedrows_.size(); }

  HighsInt getRowStart(HighsInt row) const { return ARrange_[row].first; }

  HighsInt getRowEnd(HighsInt row) const { return ARrange_[row].second; }

  const HighsInt* getARindex() const { return ARindex_.data(); }

  const double* getARvalue() const { return ARvalue_.data(); }
};

#endif
