/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNlaProductForm.cpp
 *
 * @brief Product form update methods for simplex NLA
 */
#include "simplex/HSimplexNla.h"

// #include <stdio.h>
#include <algorithm>

using std::fabs;
const HighsInt kProductFormExtraEntries = 1000;
const HighsInt kProductFormMaxUpdates = 50;
const double kProductFormPivotTolerance = 1e-8;

void ProductFormUpdate::clear() {
  valid_ = false;
  num_row_ = 0;
  update_count_ = 0;
  pivot_index_.clear();
  pivot_value_.clear();
  start_.clear();
  index_.clear();
  value_.clear();
}

void ProductFormUpdate::setup(const HighsInt num_row,
                              const double expected_density) {
  valid_ = true;
  num_row_ = num_row;
  update_count_ = 0;
  start_.push_back(0);
  HighsInt reserve_entry_space =
      kProductFormExtraEntries +
      kProductFormMaxUpdates * num_row * expected_density;
  index_.reserve(reserve_entry_space);
  value_.reserve(reserve_entry_space);
}

HighsInt ProductFormUpdate::update(HVector* aq, HighsInt* pivot_row) {
  assert(0 <= *pivot_row && *pivot_row < num_row_);
  if (update_count_ >= kProductFormMaxUpdates)
    return kRebuildReasonUpdateLimitReached;
  double pivot = aq->array[*pivot_row];
  if (fabs(pivot) < kProductFormPivotTolerance)
    return kRebuildReasonPossiblySingularBasis;
  pivot_index_.push_back(*pivot_row);
  pivot_value_.push_back(pivot);
  for (HighsInt iX = 0; iX < aq->count; iX++) {
    HighsInt iRow = aq->index[iX];
    if (iRow == *pivot_row) continue;
    index_.push_back(iRow);
    value_.push_back(aq->array[iRow]);
  }
  start_.push_back(index_.size());
  update_count_++;
  return kRebuildReasonNo;
}

void ProductFormUpdate::btran(HVector& rhs) const {
  if (!valid_) return;
  assert(rhs.size == num_row_);
  assert((int)start_.size() == update_count_ + 1);
  for (HighsInt iX = update_count_ - 1; iX >= 0; iX--) {
    const HighsInt pivot_index = pivot_index_[iX];
    double pivot_value = rhs.array[pivot_index];
    for (HighsInt iEl = start_[iX]; iEl < start_[iX + 1]; iEl++)
      pivot_value -= value_[iEl] * rhs.array[index_[iEl]];
    pivot_value /= pivot_value_[iX];
    if (rhs.array[pivot_index] == 0) rhs.index[rhs.count++] = pivot_index;
    rhs.array[pivot_index] =
        (fabs(pivot_value) < kHighsTiny) ? 1e-100 : pivot_value;
  }
}

void ProductFormUpdate::ftran(HVector& rhs) const {
  if (!valid_) return;
  assert(rhs.size == num_row_);
  assert((int)start_.size() == update_count_ + 1);
  // Use the zeroed rhs.cwork to record whether a row is in the index
  // list. If RHS fill-in occurs in a row, then we have to add it to
  // the list. We're not tracking cancellation, so we don't need to
  // know where a row appears in the list
  vector<HighsBool>& in_index = rhs.cwork;
  for (HighsInt iX = 0; iX < rhs.count; iX++) in_index[rhs.index[iX]] = true;

  for (HighsInt iX = 0; iX < update_count_; iX++) {
    const HighsInt pivot_index = pivot_index_[iX];
    double pivot_value = rhs.array[pivot_index];
    if (fabs(pivot_value) > kHighsTiny) {
      assert(in_index[pivot_index]);
      pivot_value /= pivot_value_[iX];
      rhs.array[pivot_index] = pivot_value;
      for (HighsInt iEl = start_[iX]; iEl < start_[iX + 1]; iEl++) {
        HighsInt iRow = index_[iEl];
        rhs.array[iRow] -= pivot_value * value_[iEl];
        if (in_index[iRow]) continue;
        in_index[iRow] = true;
        rhs.index[rhs.count++] = iRow;
      }
    } else {
      rhs.array[pivot_index] = 0;
    }
  }
  // Reset the in_index entries used to point into the index list
  for (HighsInt iX = 0; iX < rhs.count; iX++) in_index[rhs.index[iX]] = false;
}
