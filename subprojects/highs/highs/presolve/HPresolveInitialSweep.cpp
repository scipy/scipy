/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolveInitialSweep.h"

#include "presolve/HPresolveUtils.h"
#include "presolve/HighsPostsolveStack.h"
#include "util/HighsCDouble.h"

#define CHECKED_CALL(call)                                \
  do {                                                    \
    HPresolveInitialSweep::Result __result = call;        \
    if (__result != HPresolveInitialSweep::Result::kOk) { \
      return __result;                                    \
    }                                                     \
  } while (0)

namespace presolve {

HPresolveInitialSweep::HPresolveInitialSweep(
    HighsLp& lp, const HighsOptions& options,
    const std::vector<HighsBool>& allow_rule, const double primal_feastol)
    : lp_(lp),
      options_(options),
      allow_rule_(allow_rule),
      primal_feastol_(primal_feastol),
      num_deleted_rows_(0),
      num_deleted_cols_(0) {}

HPresolveInitialSweep::Result HPresolveInitialSweep::checkColBounds(
    HighsInt col, bool& isFixed) {
  double boundDiff = lp_.col_upper_[col] - lp_.col_lower_[col];
  double max_abs_col_value = 0;
  for (HighsInt iEl = lp_.a_matrix_.start_[col];
       iEl < lp_.a_matrix_.start_[col + 1]; iEl++)
    max_abs_col_value =
        std::max(std::abs(lp_.a_matrix_.value_[iEl]), max_abs_col_value);
  isFixed = false;
  // Check for simple infeasibility in the original LP should
  // already have been carried out in
  // HPresolve::checkOriginalLpBounds()
  assert(boundDiff >= 0);
  if (boundDiff <= primal_feastol_ &&
      (boundDiff <= options_.small_matrix_value ||
       max_abs_col_value * boundDiff <= primal_feastol_)) {
    // check for unboundedness
    if (std::abs(lp_.col_lower_[col]) == kHighsInf)
      return Result::kDualInfeasible;
    // column is fixed
    isFixed = true;
  }
  return Result::kOk;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::emptyCol(
    HighsPostsolveStack& postsolve_stack, HighsInt col) {
  const HighsInt col_nnz =
      lp_.a_matrix_.start_[col + 1] - lp_.a_matrix_.start_[col];
  assert(col_nnz == 0);
  double cost = lp_.col_cost_[col];
  const double lower = lp_.col_lower_[col];
  const double upper = lp_.col_upper_[col];

  if ((cost > 0 && lower == -kHighsInf) || (cost < 0 && upper == kHighsInf)) {
    if (std::abs(cost) <= options_.dual_feasibility_tolerance)
      cost = 0;
    else
      return Result::kDualInfeasible;
  }
  double fixval = kHighsInf;
  if (cost > 0) {
    fixval = lower;
    if (fixval == -kHighsInf) return Result::kDualInfeasible;
  } else if (cost < 0 || std::abs(upper) < std::abs(lower)) {
    fixval = upper;
    if (fixval == kHighsInf) return Result::kDualInfeasible;
  } else if (lower != -kHighsInf) {
    fixval = lower;
    if (fixval == -kHighsInf) return Result::kDualInfeasible;
  } else {
    fixval = 0.0;
  }
  assert(fixval != kHighsInf);
  postsolve_stack.removedModelFixedCol(col, fixval, cost, col_nnz, nullptr,
                                       nullptr);
  lp_.offset_ += cost * fixval;
  if (col == lp_.fme_obj_col_) lp_.fme_obj_col_ = -1;
  return Result::kOk;
}

void HPresolveInitialSweep::removeFixedCol(HighsInt col) {
  double fixval = lp_.col_lower_[col];
  if (col == lp_.fme_obj_col_) lp_.fme_obj_col_ = -1;
  for (HighsInt iEl = lp_.a_matrix_.start_[col];
       iEl < lp_.a_matrix_.start_[col + 1]; iEl++) {
    HighsInt colrow = lp_.a_matrix_.index_[iEl];
    double colval = lp_.a_matrix_.value_[iEl];
    if (lp_.row_lower_[colrow] != -kHighsInf)
      lp_.row_lower_[colrow] -= colval * fixval;
    if (lp_.row_upper_[colrow] != kHighsInf)
      lp_.row_upper_[colrow] -= colval * fixval;
  }
  lp_.offset_ += lp_.col_cost_[col] * fixval;
  assert(std::isfinite(lp_.offset_));
  lp_.col_cost_[col] = 0;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::emptyRow(
    HighsPostsolveStack& postsolve_stack, HighsInt row) {
  if (lp_.row_upper_[row] < -primal_feastol_ ||
      lp_.row_lower_[row] > primal_feastol_)
    return Result::kPrimalInfeasible;
  postsolve_stack.redundantRow(row);
  return Result::kOk;
}

double HPresolveInitialSweep::getMaxAbsColVal(HighsInt col) const {
  double maxVal = 0.0;
  for (HighsInt iEl = lp_.a_matrix_.start_[col];
       iEl < lp_.a_matrix_.start_[col + 1]; iEl++)
    maxVal = std::max(std::abs(lp_.a_matrix_.value_[iEl]), maxVal);
  return maxVal;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::singletonRow(
    HighsPostsolveStack& postsolve_stack, HighsInt row, HighsInt col,
    double val) {
  assert(std::abs(val) > options_.small_matrix_value);

  double lb, ub;
  bool lowerTightened, upperTightened;
  SingletonRowResult sr = computeSingletonRowBounds(
      val, lp_.row_lower_[row], lp_.row_upper_[row], lp_.col_lower_[col],
      lp_.col_upper_[col], primal_feastol_, getMaxAbsColVal(col),
      lp_.integrality_[col] != HighsVarType::kContinuous, lb, ub,
      lowerTightened, upperTightened);
  if (sr == SingletonRowResult::kRedundant) {
    postsolve_stack.redundantRow(row);
    return Result::kOk;
  }
  if (sr == SingletonRowResult::kPrimalInfeasible)
    return Result::kPrimalInfeasible;

  postsolve_stack.singletonRow(row, col, val, lowerTightened, upperTightened);

  lp_.col_lower_[col] = lb;
  lp_.col_upper_[col] = ub;
  return Result::kOk;
}

HPresolveInitialSweep::Result HPresolveInitialSweep::run(
    HighsPostsolveStack& postsolve_stack) {
  const bool have_col_names = lp_.col_names_.size() > 0;
  const bool have_row_names = lp_.row_names_.size() > 0;
  const HighsInt original_num_col = lp_.num_col_;
  const HighsInt original_num_row = lp_.num_row_;

  HighsInt num_fixed_col = 0;
  HighsInt num_empty_col = 0;
  HighsInt num_empty_row = 0;
  HighsInt num_singleton_row = 0;
  HighsInt num_redundant_row = 0;
  HighsInt num_col = 0;
  HighsInt nnz = 0;
  bool isFixed;

  std::vector<HighsInt> newColIndex(lp_.num_col_);
  std::vector<HighsInt> row_count(lp_.num_row_, 0);
  // Col of row is used to identify the column containing each
  // singleton row, and val_of_row the matrix entry of the singleton
  std::vector<HighsInt> col_of_row(lp_.num_row_, -1);
  std::vector<double> val_of_row(lp_.num_row_, 0);
  const HighsBool compute_implied_bounds =
      allow_rule_[kPresolveRuleRedundantRow];
  // Compute the implied bounds on rows
  std::vector<HighsCDouble> implied_row_lower;
  std::vector<HighsCDouble> implied_row_upper;
  if (compute_implied_bounds) {
    implied_row_lower.assign(lp_.num_row_, 0);
    implied_row_upper.assign(lp_.num_row_, 0);
  }

  // Pass through the columns, identifying any that are empty or
  // fixed, so can be removed, updating the LP in place.
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    HighsInt col_nnz =
        lp_.a_matrix_.start_[iCol + 1] - lp_.a_matrix_.start_[iCol];
    CHECKED_CALL(checkColBounds(iCol, isFixed));
    if (allow_rule_[kPresolveRuleEmptyCol] && col_nnz == 0) {
      // Remove empty column
      CHECKED_CALL(emptyCol(postsolve_stack, iCol));
      newColIndex[iCol] = -1;
      num_empty_col++;
      num_deleted_cols_++;
    } else if (allow_rule_[kPresolveRuleFixedCol] && isFixed) {
      // Remove fixed column
      // use pointer arithmetic since iEl equals the number of nonzeros if
      // this column and all remaining columns are empty
      HighsInt iEl = lp_.a_matrix_.start_[iCol];
      postsolve_stack.removedModelFixedCol(
          iCol, lp_.col_lower_[iCol], lp_.col_cost_[iCol], col_nnz,
          lp_.a_matrix_.index_.data() + iEl, lp_.a_matrix_.value_.data() + iEl);
      removeFixedCol(iCol);
      newColIndex[iCol] = -1;
      num_fixed_col++;
      num_deleted_cols_++;
    } else {
      // Column retained: update the LP in place by shifting the cost,
      // bounds, any names, and the matrix data. Also compute this
      // column's contribution to the implied row bounds.
      newColIndex[iCol] = num_col;
      if (num_col < iCol) {
        lp_.col_cost_[num_col] = lp_.col_cost_[iCol];
        lp_.col_lower_[num_col] = lp_.col_lower_[iCol];
        lp_.col_upper_[num_col] = lp_.col_upper_[iCol];
        lp_.integrality_[num_col] = lp_.integrality_[iCol];
        if (have_col_names)
          lp_.col_names_[num_col] = std::move(lp_.col_names_[iCol]);
      }
      HighsInt from_os = lp_.a_matrix_.start_[iCol];
      HighsInt new_col_start = nnz;
      for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
        HighsInt iRow = lp_.a_matrix_.index_[from_os + iEl];
        double value = lp_.a_matrix_.value_[from_os + iEl];
        row_count[iRow]++;
        col_of_row[iRow] = num_col;
        val_of_row[iRow] = value;
        lp_.a_matrix_.index_[nnz] = iRow;
        lp_.a_matrix_.value_[nnz] = value;
        nnz++;
        if (compute_implied_bounds) {
          double row_lower_bnd =
              value > 0 ? lp_.col_lower_[num_col] : lp_.col_upper_[num_col];
          double row_upper_bnd =
              value > 0 ? lp_.col_upper_[num_col] : lp_.col_lower_[num_col];
          if (std::abs(row_lower_bnd) == kHighsInf)
            implied_row_lower[iRow] = static_cast<HighsCDouble>(-kHighsInf);
          else if (static_cast<double>(implied_row_lower[iRow]) > -kHighsInf)
            implied_row_lower[iRow] +=
                static_cast<HighsCDouble>(value) * row_lower_bnd;
          if (std::abs(row_upper_bnd) == kHighsInf)
            implied_row_upper[iRow] = static_cast<HighsCDouble>(kHighsInf);
          else if (static_cast<double>(implied_row_upper[iRow]) < kHighsInf)
            implied_row_upper[iRow] +=
                static_cast<HighsCDouble>(value) * row_upper_bnd;
        }
      }
      lp_.a_matrix_.start_[num_col] = new_col_start;
      num_col++;
    }
  }
  lp_.a_matrix_.start_[num_col] = nnz;
  HighsInt num_removed_cols = num_empty_col + num_fixed_col;
  assert(num_deleted_cols_ == num_removed_cols);
  assert(num_col + num_removed_cols == lp_.num_col_);
  lp_.col_cost_.resize(num_col);
  lp_.col_lower_.resize(num_col);
  lp_.col_upper_.resize(num_col);
  lp_.integrality_.resize(num_col);
  if (have_col_names) lp_.col_names_.resize(num_col);
  lp_.num_col_ = num_col;
  lp_.a_matrix_.num_col_ = num_col;
  lp_.a_matrix_.start_.resize(num_col + 1);
  lp_.a_matrix_.index_.resize(nnz);
  lp_.a_matrix_.value_.resize(nnz);
  if (lp_.fme_obj_col_ >= 0) lp_.fme_obj_col_ = newColIndex[lp_.fme_obj_col_];
  postsolve_stack.compressColIndexMap(newColIndex);

  auto isRedundant = [&](const HighsInt iRow) {
    return static_cast<double>(implied_row_lower[iRow]) >=
               lp_.row_lower_[iRow] - primal_feastol_ &&
           static_cast<double>(implied_row_upper[iRow]) <=
               lp_.row_upper_[iRow] + primal_feastol_;
  };

  // Row pass: count empty, singleton, and redundant rows
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    if (allow_rule_[kPresolveRuleEmptyRow] && row_count[iRow] == 0)
      num_empty_row++;
    else if (allow_rule_[kPresolveRuleSingletonRow] && row_count[iRow] == 1)
      num_singleton_row++;
    else if (allow_rule_[kPresolveRuleRedundantRow] && isRedundant(iRow))
      num_redundant_row++;
  }

  HighsInt num_removed_rows =
      num_empty_row + num_singleton_row + num_redundant_row;
  if (num_empty_row || num_singleton_row || num_redundant_row) {
    HighsInt num_row = 0;
    std::vector<HighsBool> has_singleton_row(lp_.num_col_, false);
    std::vector<HighsInt> newRowIndex(lp_.num_row_);
    for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
      if (allow_rule_[kPresolveRuleEmptyRow] && row_count[iRow] == 0) {
        // Remove empty row
        CHECKED_CALL(emptyRow(postsolve_stack, iRow));
        newRowIndex[iRow] = -1;
        num_deleted_rows_++;
      } else if (allow_rule_[kPresolveRuleSingletonRow] &&
                 row_count[iRow] == 1) {
        // Remove singleton row
        has_singleton_row[col_of_row[iRow]] = true;
        assert(val_of_row[iRow]);
        CHECKED_CALL(singletonRow(postsolve_stack, iRow, col_of_row[iRow],
                                  val_of_row[iRow]));
        newRowIndex[iRow] = -1;
        num_deleted_rows_++;
      } else if (allow_rule_[kPresolveRuleRedundantRow] && isRedundant(iRow)) {
        // Remove redundant row
        postsolve_stack.redundantRow(iRow);
        newRowIndex[iRow] = -1;
        num_deleted_rows_++;
      } else {
        newRowIndex[iRow] = num_row;
        if (num_row < iRow) {
          lp_.row_lower_[num_row] = lp_.row_lower_[iRow];
          lp_.row_upper_[num_row] = lp_.row_upper_[iRow];
          if (have_row_names)
            lp_.row_names_[num_row] = std::move(lp_.row_names_[iRow]);
        }
        num_row++;
      }
    }
    assert(num_deleted_rows_ == num_removed_rows);
    assert(num_row + num_removed_rows == lp_.num_row_);

    if (num_redundant_row == 0) {
      // Only removing entries corresponding to singleton rows so
      // there are few to remove and it can be done efficiently
      nnz = 0;
      HighsInt from_col = 0;
      // Lambda for shifting column data and updating row indices
      auto shiftCols = [&](const HighsInt to_col) {
        for (HighsInt iCol = from_col; iCol < to_col; iCol++) {
          HighsInt from_os = lp_.a_matrix_.start_[iCol];
          HighsInt col_nnz = lp_.a_matrix_.start_[iCol + 1] - from_os;
          HighsInt new_col_start = nnz;
          for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
            HighsInt iRow = lp_.a_matrix_.index_[from_os + iEl];
            HighsInt newRow = newRowIndex[iRow];
            assert(newRow >= 0);
            lp_.a_matrix_.index_[nnz] = newRow;
            lp_.a_matrix_.value_[nnz] = lp_.a_matrix_.value_[from_os + iEl];
            nnz++;
          }
          lp_.a_matrix_.start_[iCol] = new_col_start;
        }
      };
      for (HighsInt iCol0 = 0; iCol0 < lp_.num_col_; iCol0++) {
        if (!has_singleton_row[iCol0]) continue;
        // Column iCol0 contains a row singleton, so update the matrix
        // entries for the columns since the last with a row singleton
        shiftCols(iCol0);
        HighsInt from_os = lp_.a_matrix_.start_[iCol0];
        HighsInt col_nnz = lp_.a_matrix_.start_[iCol0 + 1] - from_os;
        HighsInt new_col_start = nnz;
        bool found_row_singleton = false;
        for (HighsInt iEl = 0; iEl < col_nnz; iEl++) {
          HighsInt iRow = lp_.a_matrix_.index_[from_os + iEl];
          HighsInt newRow = newRowIndex[iRow];
          if (newRow >= 0) {
            lp_.a_matrix_.index_[nnz] = newRow;
            lp_.a_matrix_.value_[nnz] = lp_.a_matrix_.value_[from_os + iEl];
            nnz++;
          } else {
            assert(row_count[iRow] == 1);
            assert(col_of_row[iRow] == iCol0);
            assert(val_of_row[iRow] == lp_.a_matrix_.value_[from_os + iEl]);
            found_row_singleton = true;
          }
        }
        assert(found_row_singleton);
        lp_.a_matrix_.start_[iCol0] = new_col_start;
        from_col = iCol0 + 1;
      }
      // Update the matrix entries for the columns since the last with a
      // row singleton
      shiftCols(lp_.num_col_);
      lp_.a_matrix_.start_[num_col] = nnz;
    } else {
      // Also removing redundant rows, so make the matrix rowwise and
      // remove rows simply above
      nnz = 0;
      HighsInt from_row = 0;
      num_row = 0;
      // Lambda for shifting row data and updating row indices
      auto shiftRows = [&](const HighsInt to_row) {
        for (HighsInt iRow = from_row; iRow < to_row; iRow++) {
          HighsInt new_row_start = nnz;
          for (HighsInt iEl = lp_.a_matrix_.start_[iRow];
               iEl < lp_.a_matrix_.start_[iRow + 1]; iEl++) {
            lp_.a_matrix_.index_[nnz] = lp_.a_matrix_.index_[iEl];
            lp_.a_matrix_.value_[nnz] = lp_.a_matrix_.value_[iEl];
            nnz++;
          }
          lp_.a_matrix_.start_[num_row] = new_row_start;
          num_row++;
        }
      };
      lp_.a_matrix_.ensureRowwise();
      for (HighsInt iRow0 = 0; iRow0 < lp_.num_row_; iRow0++) {
        if (newRowIndex[iRow0] >= 0) continue;
        // Row iRow0 is removed, so update the matrix entries for the
        // rows since the last removed
        shiftRows(iRow0);
        from_row = iRow0 + 1;
      }
      // Update the matrix entries for the rows since the last removed
      shiftRows(lp_.num_row_);
      assert(num_row + num_removed_rows == lp_.num_row_);
      lp_.a_matrix_.start_[num_row] = nnz;
      lp_.a_matrix_.num_row_ = num_row;
      lp_.a_matrix_.ensureColwise();
    }
    lp_.row_lower_.resize(num_row);
    lp_.row_upper_.resize(num_row);
    if (have_row_names) lp_.row_names_.resize(num_row);
    lp_.num_row_ = num_row;
    lp_.a_matrix_.num_row_ = num_row;
    lp_.a_matrix_.index_.resize(nnz);
    lp_.a_matrix_.value_.resize(nnz);
    postsolve_stack.compressRowIndexMap(newRowIndex);
  }
  // Add doubleton equations, column singletons, variable locks

  if (num_fixed_col || num_empty_col)
    highsLogUser(
        options_.log_options, HighsLogType::kInfo,
        "Initial sweep removes %d + %d = %d / %d empty + fixed columns\n",
        int(num_empty_col), int(num_fixed_col), int(num_removed_cols),
        int(original_num_col));
  if (num_empty_row || num_singleton_row || num_redundant_row)
    highsLogUser(options_.log_options, HighsLogType::kInfo,
                 "Initial sweep identifies %d + %d + %d = %d / %d empty + "
                 "singleton + redundant rows\n",
                 int(num_empty_row), int(num_singleton_row),
                 int(num_redundant_row), int(num_removed_rows),
                 int(original_num_row));

  return Result::kOk;
}

}  // namespace presolve
