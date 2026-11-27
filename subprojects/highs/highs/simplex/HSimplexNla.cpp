/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HSimplexNla.cpp
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#include "simplex/HSimplexNla.h"

#include <stdio.h>

#include "../extern/pdqsort/pdqsort.h"
#include "parallel/HighsParallel.h"
#include "simplex/HSimplex.h"

using std::vector;

void HSimplexNla::setup(const HighsLp* lp, HighsInt* basic_index,
                        const HighsOptions* options, HighsTimer* timer,
                        HighsSimplexAnalysis* analysis,
                        const HighsSparseMatrix* factor_a_matrix,
                        const double factor_pivot_threshold) {
  this->setLpAndScalePointers(lp);
  this->basic_index_ = basic_index;
  this->options_ = options;
  this->timer_ = timer;
  this->analysis_ = analysis;
  this->report_ = false;
  this->factor_.setupGeneral(
      this->lp_->num_col_, this->lp_->num_row_, this->lp_->num_row_,
      factor_a_matrix->start_.data(), factor_a_matrix->index_.data(),
      factor_a_matrix->value_.data(), this->basic_index_,
      factor_pivot_threshold, this->options_->factor_pivot_tolerance,
      this->options_->highs_debug_level, &(this->options_->log_options));
  assert(debugCheckData("After HSimplexNla::setup") == HighsDebugStatus::kOk);
}

void HSimplexNla::setLpAndScalePointers(const HighsLp* for_lp) {
  this->lp_ = for_lp;
  this->scale_ = NULL;
  if (for_lp->scale_.has_scaling && !for_lp->is_scaled_)
    this->scale_ = &(for_lp->scale_);
}

void HSimplexNla::setBasicIndexPointers(HighsInt* basic_index) {
  this->basic_index_ = basic_index;
  this->factor_.basic_index = basic_index;
}

void HSimplexNla::setPointers(const HighsLp* for_lp,
                              const HighsSparseMatrix* factor_a_matrix,
                              HighsInt* basic_index,
                              const HighsOptions* options, HighsTimer* timer,
                              HighsSimplexAnalysis* analysis) {
  this->setLpAndScalePointers(for_lp);
  if (factor_a_matrix) factor_.setupMatrix(factor_a_matrix);
  if (basic_index) basic_index_ = basic_index;
  if (options) options_ = options;
  if (timer) timer_ = timer;
  if (analysis) analysis_ = analysis;
}

void HSimplexNla::clear() {
  lp_ = NULL;
  scale_ = NULL;
  basic_index_ = NULL;
  options_ = NULL;
  timer_ = NULL;
  analysis_ = NULL;
  report_ = false;
  build_synthetic_tick_ = 0;
}

HighsInt HSimplexNla::invert() {
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  if (analysis_->analyse_factor_time) {
    HighsInt thread_id = highs::parallel::thread_num();
#if 0  // def OPENMP
    thread_id = omp_get_thread_num();
#endif
    factor_timer_clock_pointer =
        analysis_->getThreadFactorTimerClockPtr(thread_id);
  }
  HighsInt rank_deficiency = factor_.build(factor_timer_clock_pointer);
  // Must not have timed out
  assert(rank_deficiency >= 0);
  build_synthetic_tick_ = factor_.build_synthetic_tick;
  return rank_deficiency;
}

void HSimplexNla::btran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  applyBasisMatrixColScale(rhs);
  btranInScaledSpace(rhs, expected_density, factor_timer_clock_pointer);
  applyBasisMatrixRowScale(rhs);
}

void HSimplexNla::ftran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  applyBasisMatrixRowScale(rhs);
  ftranInScaledSpace(rhs, expected_density, factor_timer_clock_pointer);
  applyBasisMatrixColScale(rhs);
}

void HSimplexNla::btranInScaledSpace(
    HVector& rhs, const double expected_density,
    HighsTimerClock* factor_timer_clock_pointer) const {
  factor_.btranCall(rhs, expected_density, factor_timer_clock_pointer);
}

void HSimplexNla::ftranInScaledSpace(
    HVector& rhs, const double expected_density,
    HighsTimerClock* factor_timer_clock_pointer) const {
  factor_.ftranCall(rhs, expected_density, factor_timer_clock_pointer);
}

void HSimplexNla::update(HVector* aq, HVector* ep, HighsInt* iRow,
                         HighsInt* hint) {
  reportPackValue("  pack: aq Bf ", aq);
  reportPackValue("  pack: ep Bf ", ep);
  factor_.refactor_info_.clear();
  if (!update_.valid_) {
    factor_.update(aq, ep, iRow, hint);
  } else {
    *hint = update_.update(aq, iRow);
  }
}

double HSimplexNla::rowEp2NormInScaledSpace(const HighsInt iRow,
                                            const HVector& row_ep) const {
  if (scale_ == NULL) {
    return row_ep.norm2();
  }
  const vector<double>& row_scale = scale_->row;
  // Get the 2-norm of row_ep in the scaled space otherwise for
  // checking
  const bool DSE_check = false;
  double alt_row_ep_2norm = 0;
  if (DSE_check) {
    HVector alt_row_ep;
    alt_row_ep.setup(lp_->num_row_);
    alt_row_ep.clear();
    alt_row_ep.count = 1;
    alt_row_ep.index[0] = iRow;
    alt_row_ep.array[iRow] = 1;
    alt_row_ep.packFlag = false;
    factor_.btranCall(alt_row_ep, 0);
    alt_row_ep_2norm = alt_row_ep.norm2();
  }
  // Get the 2-norm of row_ep in the scaled space
  //
  // Determine the scaling that was applied to the unit RHS before
  // scaled BTRAN. This must be unapplied to all components of the
  // result.
  double col_scale_value = basicColScaleFactor(iRow);
  // Now compute the 2-norm of row_ep in the scaled space, unapplying
  // the scaling that was applied after BTRAN
  double row_ep_2norm = 0;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(row_ep.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_row_indices ? row_ep.index[iEntry] : iEntry;
    const double value_in_scaled_space =
        row_ep.array[iRow] / (row_scale[iRow] * col_scale_value);
    row_ep_2norm += value_in_scaled_space * value_in_scaled_space;
  }
  if (DSE_check) {
    const double error = std::fabs(row_ep_2norm - alt_row_ep_2norm) /
                         std::max(1.0, row_ep_2norm);
    if (error > 1e-4)
      printf(
          "rowEp2NormInScaledSpace: iRow = %2d has deduced norm = %10.4g and "
          "alt "
          "norm = %10.4g, giving error %10.4g\n",
          (int)iRow, row_ep_2norm, alt_row_ep_2norm, error);
  }
  return row_ep_2norm;
}

void HSimplexNla::transformForUpdate(HVector* aq, HVector* ep,
                                     const HighsInt variable_in,
                                     const HighsInt row_out) {
  if (scale_ == NULL) return;
  // For (\hat)aq, UPDATE needs packValue and array[row_out] to
  // correspond to \bar{B}^{-1}(R.aq.cq), but CB.\bar{B}^{-1}(R.aq)
  // has been computed.
  //
  // Hence packValue and array[row_out] need to be scaled by cq;
  //
  // array[row_out] has to be unscaled by the corresponding entry of
  // CB
  //
  reportPackValue("pack aq Bf ", aq);
  double cq_scale_factor = variableScaleFactor(variable_in);

  for (HighsInt ix = 0; ix < aq->packCount; ix++)
    aq->packValue[ix] *= cq_scale_factor;
  reportPackValue("pack aq Af ", aq);
  //
  // Now focus on the pivot value, aq->array[row_out]
  double pivot_in_scaled_space = pivotInScaledSpace(aq, variable_in, row_out);
  // First scale by cq
  aq->array[row_out] *= cq_scale_factor;
  //
  // Also have to unscale by cp
  double cp_scale_factor = basicColScaleFactor(row_out);
  aq->array[row_out] /= cp_scale_factor;
  assert(pivot_in_scaled_space == aq->array[row_out]);
  // For (\hat)ep, UPDATE needs packValue to correspond to
  // \bar{B}^{-T}ep, but R.\bar{B}^{-T}(CB.ep) has been computed.
  //
  // Hence packValue needs to be unscaled by cp
  for (HighsInt ix = 0; ix < ep->packCount; ix++)
    ep->packValue[ix] /= cp_scale_factor;
}

double HSimplexNla::variableScaleFactor(const HighsInt iVar) const {
  if (scale_ == NULL) return 1.0;
  return iVar < lp_->num_col_ ? scale_->col[iVar]
                              : 1.0 / scale_->row[iVar - lp_->num_col_];
}

double HSimplexNla::basicColScaleFactor(const HighsInt iCol) const {
  if (scale_ == NULL) return 1.0;
  return variableScaleFactor(basic_index_[iCol]);
}

double HSimplexNla::pivotInScaledSpace(const HVector* aq,
                                       const HighsInt variable_in,
                                       const HighsInt row_out) const {
  return aq->array[row_out] * variableScaleFactor(variable_in) /
         variableScaleFactor(basic_index_[row_out]);
}

void HSimplexNla::setPivotThreshold(const double new_pivot_threshold) {
  factor_.setPivotThreshold(new_pivot_threshold);
}

void HSimplexNla::applyBasisMatrixRowScale(HVector& rhs) const {
  if (scale_ == NULL) return;
  const vector<double>& row_scale = scale_->row;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(rhs.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_row_indices ? rhs.index[iEntry] : iEntry;
    rhs.array[iRow] *= row_scale[iRow];
  }
}

void HSimplexNla::applyBasisMatrixColScale(HVector& rhs) const {
  if (scale_ == NULL) return;
  const vector<double>& col_scale = scale_->col;
  const vector<double>& row_scale = scale_->row;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(rhs.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iCol = use_row_indices ? rhs.index[iEntry] : iEntry;
    HighsInt iVar = basic_index_[iCol];
    if (iVar < lp_->num_col_) {
      rhs.array[iCol] *= col_scale[iVar];
    } else {
      rhs.array[iCol] /= row_scale[iVar - lp_->num_col_];
    }
  }
}

void HSimplexNla::unapplyBasisMatrixRowScale(HVector& rhs) const {
  if (scale_ == NULL) return;
  const vector<double>& row_scale = scale_->row;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(rhs.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    const HighsInt iRow = use_row_indices ? rhs.index[iEntry] : iEntry;
    rhs.array[iRow] /= row_scale[iRow];
  }
}

void HSimplexNla::addCols(const HighsLp* updated_lp) {
  // Adding columns is easy, since they are nonbasic
  //
  // Set the pointers for the LP and scaling. The pointer to the
  // vector of basic variables isn't updated, since it hasn't been
  // resized. The HFactor matrix isn't needed until reinversion has to
  // be performed
  setLpAndScalePointers(updated_lp);
}

void HSimplexNla::addRows(const HighsLp* updated_lp, HighsInt* basic_index,
                          const HighsSparseMatrix* scaled_ar_matrix) {
  // Adding rows is not so easy, since their slacks are basic
  //
  // Set the pointers for the LP, scaling and basic variables. The
  // HFactor matrix isn't needed until reinversion has to be performed
  setLpAndScalePointers(updated_lp);
  basic_index_ = basic_index;
  factor_.basic_index = basic_index;
  factor_.addRows(scaled_ar_matrix);
}

bool HSimplexNla::sparseLoopStyle(const HighsInt count, const HighsInt dim,
                                  HighsInt& to_entry) const {
  // Parameter to decide whether to use just the values in a HVector, or
  // use the indices of their nonzeros
  const bool use_indices = count >= 0 && count < kDensityForIndexing * dim;
  if (use_indices) {
    to_entry = count;
  } else {
    to_entry = dim;
  }
  return use_indices;
}

void HSimplexNla::reportArray(const std::string& message, const HVector* vector,
                              const bool force) const {
  reportArray(message, 0, vector, force);
}

void HSimplexNla::reportArray(const std::string& message, const HighsInt offset,
                              const HVector* vector, const bool force) const {
  if (!report_ && !force) return;
  const HighsInt num_row = lp_->num_row_;
  if (num_row > kReportItemLimit) {
    reportArraySparse(message, offset, vector, force);
  } else {
    printf("%s", message.c_str());
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow > 0 && iRow % 10 == 0)
        printf("\n                                 ");
      printf("%11.4g ", vector->array[iRow]);
    }
    printf("\n");
  }
}

void HSimplexNla::reportVector(const std::string& message,
                               const HighsInt num_index,
                               const vector<double>& vector_value,
                               const vector<HighsInt>& vector_index,
                               const bool force) const {
  if (!report_ && !force) return;
  assert((int)vector_value.size() >= num_index);
  if (num_index <= 0) return;
  const HighsInt num_row = lp_->num_row_;
  if (num_index > kReportItemLimit) {
    analyseVectorValues(nullptr, message, num_row, vector_value, true);
  } else {
    printf("%s", message.c_str());
    for (HighsInt iX = 0; iX < num_index; iX++) {
      if (iX % 5 == 0) printf("\n");
      printf("[%4d %11.4g] ", (int)vector_index[iX], vector_value[iX]);
    }
    printf("\n");
  }
}

void HSimplexNla::reportArraySparse(const std::string& message,
                                    const HVector* vector,
                                    const bool force) const {
  reportArraySparse(message, 0, vector, force);
}

void HSimplexNla::reportArraySparse(const std::string& message,
                                    const HighsInt offset,
                                    const HVector* vector,
                                    const bool force) const {
  if (!report_ && !force) return;
  const HighsInt num_row = lp_->num_row_;
  if (vector->count > kReportItemLimit) {
    analyseVectorValues(nullptr, message, num_row, vector->array, true);
  } else if (vector->count < num_row) {
    std::vector<HighsInt> sorted_index = vector->index;
    pdqsort(sorted_index.begin(), sorted_index.begin() + vector->count);
    printf("%s", message.c_str());
    for (HighsInt en = 0; en < vector->count; en++) {
      HighsInt iRow = sorted_index[en];
      if (en % 5 == 0) printf("\n");
      printf("[%4d ", (int)(iRow));
      if (offset) printf("(%4d)", (int)(offset + iRow));
      printf("%11.4g] ", vector->array[iRow]);
    }
  } else {
    if (num_row > kReportItemLimit) {
      analyseVectorValues(nullptr, message, num_row, vector->array, true);
      return;
    }
    printf("%s", message.c_str());
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow % 5 == 0) printf("\n");
      printf("%11.4g ", vector->array[iRow]);
    }
  }
  printf("\n");
}

void HSimplexNla::reportPackValue(const std::string& message,
                                  const HVector* vector,
                                  const bool force) const {
  if (!report_ && !force) return;
  if (vector->packCount > kReportItemLimit) {
    analyseVectorValues(nullptr, message, vector->packCount, vector->packValue,
                        true);
    return;
  }
  printf("%s", message.c_str());
  std::vector<HighsInt> sorted_index = vector->packIndex;
  pdqsort(sorted_index.begin(), sorted_index.begin() + vector->packCount);
  for (HighsInt en = 0; en < vector->packCount; en++) {
    HighsInt iRow = sorted_index[en];
    if (en % 5 == 0) printf("\n");
    printf("[%4d %11.4g] ", (int)iRow, vector->packValue[en]);
  }
  printf("\n");
}

HighsDebugStatus HSimplexNla::debugCheckData(const std::string& message) const {
  if (options_->highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kOk;
  std::string scale_status;
  if (scale_ == NULL) {
    scale_status = "NULL";
  } else {
    scale_status = "non-NULL";
  }
  HighsLp check_lp = *lp_;
  bool error0_found = false;
  bool error1_found = false;
  bool error2_found = false;
  bool error_found = false;
  const HighsInt* factor_Astart = factor_.getAstart();
  const HighsInt* factor_Aindex = factor_.getAindex();
  const double* factor_Avalue = factor_.getAvalue();

  if (scale_ == NULL) {
    if (factor_Astart != lp_->a_matrix_.start_.data()) error0_found = true;
    if (factor_Aindex != lp_->a_matrix_.index_.data()) error1_found = true;
    if (factor_Avalue != lp_->a_matrix_.value_.data()) error2_found = true;
    error_found = error0_found || error1_found || error2_found;
    if (error_found) {
      highsLogUser(options_->log_options, HighsLogType::kError,
                   "CheckNlaData: (%s) scale_ is %s lp_ - factor_ matrix "
                   "pointer errors\n",
                   message.c_str(), scale_status.c_str());
      if (error0_found)
        printf("a_matrix_.start_ pointer error: %p vs %p\n",
               (void*)factor_Astart, (void*)lp_->a_matrix_.start_.data());
      if (error1_found) printf("a_matrix_.index pointer error\n");
      if (error2_found) printf("a_matrix_.value pointer error\n");
      assert(!error_found);
      return HighsDebugStatus::kLogicalError;
    }
  } else {
    check_lp.applyScale();
  }
  HighsInt error_col = -1;
  for (HighsInt iCol = 0; iCol < check_lp.num_col_ + 1; iCol++) {
    if (check_lp.a_matrix_.start_[iCol] != factor_Astart[iCol]) {
      error_col = iCol;
      break;
    }
  }
  error_found = error_col >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.start_ != "
                 "factor_Astart for col %d (%d != %d)\n",
                 message.c_str(), scale_status.c_str(), (int)error_col,
                 (int)check_lp.a_matrix_.start_[error_col],
                 (int)factor_Astart[error_col]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  HighsInt nnz = check_lp.numNz();
  HighsInt error_el = -1;
  for (HighsInt iEl = 0; iEl < nnz; iEl++) {
    if (check_lp.a_matrix_.index_[iEl] != factor_Aindex[iEl]) {
      error_el = iEl;
      break;
    }
  }
  error_found = error_el >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.index_ != "
                 "factor_Aindex for el %d (%d != %d)\n",
                 message.c_str(), scale_status.c_str(), (int)error_el,
                 (int)check_lp.a_matrix_.index_[error_el],
                 (int)factor_Aindex[error_el]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  for (HighsInt iEl = 0; iEl < nnz; iEl++) {
    if (check_lp.a_matrix_.value_[iEl] != factor_Avalue[iEl]) {
      error_el = iEl;
      break;
    }
  }
  error_found = error_el >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.value_ != "
                 "factor_Avalue for el %d (%g != %g)\n",
                 message.c_str(), scale_status.c_str(), (int)error_el,
                 check_lp.a_matrix_.value_[error_el], factor_Avalue[error_el]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  return HighsDebugStatus::kOk;
}
