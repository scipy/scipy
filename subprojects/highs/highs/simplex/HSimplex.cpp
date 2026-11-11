/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplex.cpp
 * @brief
 */

#include "simplex/HSimplex.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

// #include "lp_data/HighsLpUtils.h"
// #include "util/HighsSort.h"

using std::max;
using std::min;
// using std::runtime_error;

void SimplexBasis::clear() {
  this->hash = 0;
  this->basicIndex_.clear();
  this->nonbasicFlag_.clear();
  this->nonbasicMove_.clear();
  this->debug_id = -1;
  this->debug_update_count = -1;
  this->debug_origin_name = "None";
}

void SimplexBasis::setup(const HighsInt num_col, const HighsInt num_row) {
  this->hash = 0;
  this->basicIndex_.resize(num_row);
  this->nonbasicFlag_.resize(num_col + num_row);
  this->nonbasicMove_.resize(num_col + num_row);
  this->debug_id = -1;
  this->debug_update_count = -1;
  this->debug_origin_name = "None";
}

void appendNonbasicColsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                               HighsInt XnumNewCol) {
  assert(highs_basis.valid);
  if (!highs_basis.valid) {
    printf("\n!!Appending columns to invalid basis!!\n\n");
  }
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  HighsInt newNumCol = lp.num_col_ + XnumNewCol;
  highs_basis.col_status.resize(newNumCol);
  // Make any new columns nonbasic
  for (HighsInt iCol = lp.num_col_; iCol < newNumCol; iCol++) {
    if (!highs_isInfinity(-lp.col_lower_[iCol])) {
      highs_basis.col_status[iCol] = HighsBasisStatus::kLower;
    } else if (!highs_isInfinity(lp.col_upper_[iCol])) {
      highs_basis.col_status[iCol] = HighsBasisStatus::kUpper;
    } else {
      highs_basis.col_status[iCol] = HighsBasisStatus::kZero;
    }
  }
}

void appendNonbasicColsToBasis(HighsLp& lp, SimplexBasis& basis,
                               HighsInt XnumNewCol) {
  // Add nonbasic structurals
  if (XnumNewCol == 0) return;
  HighsInt newNumCol = lp.num_col_ + XnumNewCol;
  HighsInt newNumTot = newNumCol + lp.num_row_;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  // Shift the row data in basicIndex, nonbasicFlag and nonbasicMove if
  // necessary
  for (HighsInt iRow = lp.num_row_ - 1; iRow >= 0; iRow--) {
    HighsInt iCol = basis.basicIndex_[iRow];
    if (iCol >= lp.num_col_) {
      // This basic variable is a row, so shift its index
      basis.basicIndex_[iRow] += XnumNewCol;
    }
    basis.nonbasicFlag_[newNumCol + iRow] =
        basis.nonbasicFlag_[lp.num_col_ + iRow];
    basis.nonbasicMove_[newNumCol + iRow] =
        basis.nonbasicMove_[lp.num_col_ + iRow];
  }
  // Make any new columns nonbasic
  for (HighsInt iCol = lp.num_col_; iCol < newNumCol; iCol++) {
    basis.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    int8_t move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (fabs(lower) < fabs(upper)) {
          move = kNonbasicMoveUp;
        } else {
          move = kNonbasicMoveDn;
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis.nonbasicMove_[iCol] = move;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, HighsBasis& highs_basis,
                            HighsInt XnumNewRow) {
  assert(highs_basis.valid);
  if (!highs_basis.valid) {
    printf("\n!!Appending rows to invalid basis!!\n\n");
  }
  // Add basic logicals
  if (XnumNewRow == 0) return;
  HighsInt newNumRow = lp.num_row_ + XnumNewRow;
  highs_basis.row_status.resize(newNumRow);
  // Make the new rows basic
  for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++) {
    highs_basis.row_status[iRow] = HighsBasisStatus::kBasic;
  }
}

void appendBasicRowsToBasis(HighsLp& lp, SimplexBasis& basis,
                            HighsInt XnumNewRow) {
  // Add basic logicals
  if (XnumNewRow == 0) return;

  HighsInt newNumRow = lp.num_row_ + XnumNewRow;
  HighsInt newNumTot = lp.num_col_ + newNumRow;
  basis.nonbasicFlag_.resize(newNumTot);
  basis.nonbasicMove_.resize(newNumTot);
  basis.basicIndex_.resize(newNumRow);
  // Make the new rows basic
  for (HighsInt iRow = lp.num_row_; iRow < newNumRow; iRow++) {
    basis.nonbasicFlag_[lp.num_col_ + iRow] = kNonbasicFlagFalse;
    basis.nonbasicMove_[lp.num_col_ + iRow] = 0;
    basis.basicIndex_[iRow] = lp.num_col_ + iRow;
  }
}

void getUnscaledInfeasibilities(const HighsOptions& options,
                                const HighsScale& scale,
                                const SimplexBasis& basis,
                                const HighsSimplexInfo& info,
                                HighsInfo& highs_info) {
  const double primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  const double dual_feasibility_tolerance = options.dual_feasibility_tolerance;

  HighsInt& num_primal_infeasibilities = highs_info.num_primal_infeasibilities;
  double& max_primal_infeasibility = highs_info.max_primal_infeasibility;
  double& sum_primal_infeasibilities = highs_info.sum_primal_infeasibilities;
  HighsInt& num_dual_infeasibilities = highs_info.num_dual_infeasibilities;
  double& max_dual_infeasibility = highs_info.max_dual_infeasibility;
  double& sum_dual_infeasibilities = highs_info.sum_dual_infeasibilities;

  // Zero the counts of unscaled primal and dual infeasibilities
  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;
  num_dual_infeasibilities = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibilities = 0;

  double scale_mu = 1.0;
  assert(int(scale.col.size()) == scale.num_col);
  assert(int(scale.row.size()) == scale.num_row);
  for (HighsInt iVar = 0; iVar < scale.num_col + scale.num_row; iVar++) {
    // Look at the dual infeasibilities of nonbasic variables
    if (basis.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
    // No dual infeasibility for fixed rows and columns
    if (info.workLower_[iVar] == info.workUpper_[iVar]) continue;
    bool col = iVar < scale.num_col;
    HighsInt iCol = 0;
    HighsInt iRow = 0;
    if (col) {
      iCol = iVar;
      assert(int(scale.col.size()) > iCol);
      scale_mu = 1 / (scale.col[iCol] / scale.cost);
    } else {
      iRow = iVar - scale.num_col;
      assert(int(scale.row.size()) > iRow);
      scale_mu = scale.row[iRow] * scale.cost;
    }
    const double dual = info.workDual_[iVar];
    const double lower = info.workLower_[iVar];
    const double upper = info.workUpper_[iVar];
    const double unscaled_dual = dual * scale_mu;

    double dual_infeasibility;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(unscaled_dual);
    } else {
      // Not fixed: any dual infeasibility is given by value signed by
      // nonbasicMove. This assumes that nonbasicMove=0 for fixed
      // variables
      dual_infeasibility = -basis.nonbasicMove_[iVar] * unscaled_dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= dual_feasibility_tolerance)
        num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
  }
  // Look at the primal infeasibilities of basic variables
  for (HighsInt ix = 0; ix < scale.num_row; ix++) {
    HighsInt iVar = basis.basicIndex_[ix];
    bool col = iVar < scale.num_col;
    HighsInt iCol = 0;
    HighsInt iRow = 0;
    if (col) {
      iCol = iVar;
      scale_mu = scale.col[iCol];
    } else {
      iRow = iVar - scale.num_col;
      scale_mu = 1 / scale.row[iRow];
    }
    double unscaled_lower = info.baseLower_[ix] * scale_mu;
    double unscaled_value = info.baseValue_[ix] * scale_mu;
    double unscaled_upper = info.baseUpper_[ix] * scale_mu;
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (unscaled_value < unscaled_lower - primal_feasibility_tolerance) {
      primal_infeasibility = unscaled_lower - unscaled_value;
    } else if (unscaled_value > unscaled_upper + primal_feasibility_tolerance) {
      primal_infeasibility = unscaled_value - unscaled_upper;
    }
    if (primal_infeasibility > 0) {
      num_primal_infeasibilities++;
      max_primal_infeasibility =
          max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibilities += primal_infeasibility;
    }
  }
  setSolutionStatus(highs_info);
}

void setSolutionStatus(HighsInfo& highs_info) {
  if (highs_info.num_primal_infeasibilities < 0) {
    highs_info.primal_solution_status = kSolutionStatusNone;
  } else if (highs_info.num_primal_infeasibilities == 0) {
    highs_info.primal_solution_status = kSolutionStatusFeasible;
  } else {
    highs_info.primal_solution_status = kSolutionStatusInfeasible;
  }
  if (highs_info.num_dual_infeasibilities < 0) {
    highs_info.dual_solution_status = kSolutionStatusNone;
  } else if (highs_info.num_dual_infeasibilities == 0) {
    highs_info.dual_solution_status = kSolutionStatusFeasible;
  } else {
    highs_info.dual_solution_status = kSolutionStatusInfeasible;
  }
}
// SCALING

void scaleSimplexCost(const HighsOptions& options, HighsLp& lp,
                      double& cost_scale) {
  // Scale the costs by no less than minAlwCostScale
  double max_allowed_cost_scale = pow(2.0, options.allowed_cost_scale_factor);
  double max_nonzero_cost = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (lp.col_cost_[iCol]) {
      max_nonzero_cost = max(fabs(lp.col_cost_[iCol]), max_nonzero_cost);
    }
  }
  // Scaling the costs up effectively increases the dual tolerance to
  // which the problem is solved - so, if the max cost is small the
  // scaling factor pushes it up by a power of 2 so it's close to 1
  // Scaling the costs down effectively decreases the dual tolerance
  // to which the problem is solved - so this can't be done too much
  cost_scale = 1;
  const double ln2 = log(2.0);
  // Scale if the max cost is positive and outside the range [1/16, 16]
  if ((max_nonzero_cost > 0) &&
      ((max_nonzero_cost < (1.0 / 16)) || (max_nonzero_cost > 16))) {
    cost_scale = max_nonzero_cost;
    cost_scale = pow(2.0, floor(log(cost_scale) / ln2 + 0.5));
    cost_scale = min(cost_scale, max_allowed_cost_scale);
  }
  if (cost_scale == 1) {
    highsLogUser(options.log_options, HighsLogType::kInfo,
                 "LP cost vector not scaled down: max cost is %g\n",
                 max_nonzero_cost);
    return;
  }
  // Scale the costs (and record of max_nonzero_cost) by cost_scale, being at
  // most max_allowed_cost_scale
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    lp.col_cost_[iCol] /= cost_scale;
  }
  max_nonzero_cost /= cost_scale;
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "LP cost vector scaled down by %g: max cost is %g\n", cost_scale,
               max_nonzero_cost);
}

void unscaleSimplexCost(HighsLp& lp, double cost_scale) {
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.col_cost_[iCol] *= cost_scale;
}

bool isBasisRightSize(const HighsLp& lp, const SimplexBasis& basis) {
  bool right_size = true;
  right_size =
      (HighsInt)basis.nonbasicFlag_.size() == lp.num_col_ + lp.num_row_ &&
      right_size;
  right_size =
      (HighsInt)basis.nonbasicMove_.size() == lp.num_col_ + lp.num_row_ &&
      right_size;
  right_size = (HighsInt)basis.basicIndex_.size() == lp.num_row_ && right_size;
  return right_size;
}
