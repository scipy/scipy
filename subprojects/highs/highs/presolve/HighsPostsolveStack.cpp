/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HighsPostsolveStack.h"

#include <numeric>

#include "../extern/pdqsort/pdqsort.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsModelUtils.h"  // For debugging #2001
#include "lp_data/HighsOptions.h"
#include "util/HighsCDouble.h"
#include "util/HighsUtils.h"

namespace presolve {

void HighsPostsolveStack::initializeIndexMaps(HighsInt numRow,
                                              HighsInt numCol) {
  origNumRow = numRow;
  origNumCol = numCol;
  nextRowIndex = numRow;
  nextColIndex = numCol;

  origRowIndex.resize(numRow);
  std::iota(origRowIndex.begin(), origRowIndex.end(), 0);

  origRowType.resize(numRow, OrigRowType::kOriginal);

  origColIndex.resize(numCol);
  std::iota(origColIndex.begin(), origColIndex.end(), 0);

  linearlyTransformable.resize(numCol, true);
}

void HighsPostsolveStack::compressIndexMaps(
    const std::vector<HighsInt>& newRowIndex,
    const std::vector<HighsInt>& newColIndex) {
  compressColIndexMap(newColIndex);
  compressRowIndexMap(newRowIndex);
}

void HighsPostsolveStack::compressRowIndexMap(
    const std::vector<HighsInt>& newRowIndex) {
  compressIndexMap(newRowIndex, this->origRowIndex);
  compressIndexMap(newRowIndex, this->origRowType);
}

void HighsPostsolveStack::compressColIndexMap(
    const std::vector<HighsInt>& newColIndex) {
  compressIndexMap(newColIndex, this->origColIndex);
}

void HighsPostsolveStack::LinearTransform::undo(const HighsOptions& options,
                                                HighsSolution& solution) const {
  solution.col_value[col] *= scale;
  solution.col_value[col] += constant;

  if (solution.dual_valid) solution.col_dual[col] /= scale;
}

void HighsPostsolveStack::LinearTransform::transformToPresolvedSpace(
    std::vector<double>& primalSol) const {
  primalSol[col] -= constant;
  primalSol[col] /= scale;
}

static HighsBasisStatus computeRowStatus(double dual,
                                         HighsPostsolveStack::RowType rowType) {
  if (rowType == HighsPostsolveStack::RowType::kEq)
    return dual < 0 ? HighsBasisStatus::kUpper : HighsBasisStatus::kLower;
  else if (rowType == HighsPostsolveStack::RowType::kGeq)
    return HighsBasisStatus::kLower;
  else
    return HighsBasisStatus::kUpper;
}

void HighsPostsolveStack::FreeColSubstitution::undo(
    const HighsOptions& options, const std::vector<Nonzero>& rowValues,
    const std::vector<Nonzero>& colValues, HighsSolution& solution,
    HighsBasis& basis) {
  // compute primal values
  double colCoef = 0;
  HighsCDouble rowValue = 0;
  for (const auto& rowVal : rowValues) {
    if (rowVal.index == col)
      colCoef = rowVal.value;
    else
      rowValue += rowVal.value * solution.col_value[rowVal.index];
  }

  assert(colCoef != 0);
  // Row values aren't fully postsolved, so why do this?
  solution.row_value[row] =
      static_cast<double>(rowValue + colCoef * solution.col_value[col]);
  solution.col_value[col] = static_cast<double>((rhs - rowValue) / colCoef);

  // if no dual values requested, return here
  if (!solution.dual_valid) return;

  // compute the row dual value such that reduced cost of basic column is 0
  solution.row_dual[row] = 0;
  HighsCDouble dualval = colCost;
  for (const auto& colVal : colValues) {
    dualval -= colVal.value * solution.row_dual[colVal.index];
  }
  solution.row_dual[row] = static_cast<double>(dualval / colCoef);

  solution.col_dual[col] = 0;

  // set basis status if necessary
  if (!basis.valid) return;

  basis.col_status[col] = HighsBasisStatus::kBasic;
  basis.row_status[row] = computeRowStatus(solution.row_dual[row], rowType);
}

static HighsBasisStatus computeStatus(double dual, HighsBasisStatus& status,
                                      double dual_feasibility_tolerance) {
  if (dual > dual_feasibility_tolerance)
    status = HighsBasisStatus::kLower;
  else if (dual < -dual_feasibility_tolerance)
    status = HighsBasisStatus::kUpper;

  return status;
}

static HighsBasisStatus computeStatus(double dual,
                                      double dual_feasibility_tolerance) {
  if (dual > dual_feasibility_tolerance)
    return HighsBasisStatus::kLower;
  else if (dual < -dual_feasibility_tolerance)
    return HighsBasisStatus::kUpper;
  else
    return HighsBasisStatus::kBasic;
}

void HighsPostsolveStack::DoubletonEquation::undo(
    const HighsOptions& options, const std::vector<Nonzero>& colValues,
    HighsSolution& solution, HighsBasis& basis) const {
  // retrieve the row and column index, the row side and the two
  // coefficients then compute the primal values
  solution.col_value[colSubst] = static_cast<double>(
      (rhs - static_cast<HighsCDouble>(coef) * solution.col_value[col]) /
      coefSubst);

  // can only do primal postsolve, stop here
  if (row == -1 || !solution.dual_valid) return;

  const HighsBasisStatus colStatus =
      !basis.valid
          ? computeStatus(solution.col_dual[col],
                          options.dual_feasibility_tolerance)
          : computeStatus(solution.col_dual[col], basis.col_status[col],
                          options.dual_feasibility_tolerance);

  // compute the current dual values of the row and the substituted column
  // before deciding on which column becomes basic
  // for each entry in a row i of the substituted column we added the
  // doubleton equation row with scale -a_i/substCoef. Therefore the dual
  // multiplier of this row i implicitly increases the dual multiplier of this
  // doubleton equation row with that scale.
  HighsCDouble rowDual = 0.0;
  solution.row_dual[row] = 0;
  for (const auto& colVal : colValues) {
    rowDual -= colVal.value * solution.row_dual[colVal.index];
  }
  rowDual /= coefSubst;
  solution.row_dual[row] = static_cast<double>(rowDual);
  // the equation was also added to the objective, so the current values need
  // to be changed
  solution.col_dual[colSubst] = substCost;
  solution.col_dual[col] += substCost * coef / coefSubst;

  if ((upperTightened && colStatus == HighsBasisStatus::kUpper) ||
      (lowerTightened && colStatus == HighsBasisStatus::kLower)) {
    // column must get zero reduced cost as the current bound cannot be used
    // so alter the dual multiplier of the row to make the dual multiplier of
    // column zero
    double rowDualDelta = solution.col_dual[col] / coef;
    solution.row_dual[row] = static_cast<double>(rowDual + rowDualDelta);
    solution.col_dual[col] = 0.0;
    solution.col_dual[colSubst] = static_cast<double>(
        static_cast<HighsCDouble>(solution.col_dual[colSubst]) -
        rowDualDelta * coefSubst);

    if (basis.valid) {
      if ((std::signbit(coef) == std::signbit(coefSubst) &&
           basis.col_status[col] == HighsBasisStatus::kUpper) ||
          (std::signbit(coef) != std::signbit(coefSubst) &&
           basis.col_status[col] == HighsBasisStatus::kLower))
        basis.col_status[colSubst] = HighsBasisStatus::kLower;
      else
        basis.col_status[colSubst] = HighsBasisStatus::kUpper;
      basis.col_status[col] = HighsBasisStatus::kBasic;
    }
  } else {
    // otherwise make the reduced cost of the substituted column zero and make
    // that column basic
    double rowDualDelta = solution.col_dual[colSubst] / coefSubst;
    solution.row_dual[row] = static_cast<double>(rowDual + rowDualDelta);
    solution.col_dual[colSubst] = 0.0;
    solution.col_dual[col] =
        static_cast<double>(static_cast<HighsCDouble>(solution.col_dual[col]) -
                            rowDualDelta * coef);
    if (basis.valid) basis.col_status[colSubst] = HighsBasisStatus::kBasic;
  }

  if (!basis.valid) return;

  basis.row_status[row] = computeRowStatus(solution.row_dual[row], rowType);
}

void HighsPostsolveStack::EqualityRowAddition::undo(
    const HighsOptions& options, const std::vector<Nonzero>& eqRowValues,
    HighsSolution& solution, HighsBasis& basis) const {
  // nothing more to do if the row is zero in the dual solution or there is
  // no dual solution
  if (!solution.dual_valid || solution.row_dual[row] == 0.0) return;

  // the dual multiplier of the row implicitly increases the dual multiplier
  // of the equation with the scale the equation was added with
  solution.row_dual[addedEqRow] = static_cast<double>(
      static_cast<HighsCDouble>(eqRowScale) * solution.row_dual[row] +
      solution.row_dual[addedEqRow]);

  assert(!basis.valid);
}

void HighsPostsolveStack::EqualityRowAdditions::undo(
    const HighsOptions& options, const std::vector<Nonzero>& eqRowValues,
    const std::vector<Nonzero>& targetRows, HighsSolution& solution,
    HighsBasis& basis) const {
  // nothing more to do if the row is zero in the dual solution or there is
  // no dual solution
  if (!solution.dual_valid) return;

  // the dual multiplier of the rows where the eq row was added implicitly
  // increases the dual multiplier of the equation with the scale that was
  // used for adding the equation
  HighsCDouble eqRowDual = solution.row_dual[addedEqRow];
  for (const auto& targetRow : targetRows) {
    eqRowDual += static_cast<HighsCDouble>(targetRow.value) *
                 solution.row_dual[targetRow.index];
  }
  solution.row_dual[addedEqRow] = static_cast<double>(eqRowDual);

  assert(!basis.valid);
}

void HighsPostsolveStack::ForcingColumn::undo(
    const HighsOptions& options, const std::vector<Nonzero>& colValues,
    HighsSolution& solution, HighsBasis& basis) const {
  HighsInt nonbasicRow = -1;
  HighsBasisStatus nonbasicRowStatus = HighsBasisStatus::kNonbasic;
  double colValFromNonbasicRow = colBound;

  HighsInt debug_num_use_row_value = 0;
  const bool debug_report = false;

  auto computeColVal = [&](HighsInt direction) {
    // choose solution value that makes all rows feasible
    for (const auto& colVal : colValues) {
      // Row values aren't fully postsolved, so how can this work?
      debug_num_use_row_value++;
      double colValFromRow = solution.row_value[colVal.index] / colVal.value;
      if (direction * colValFromRow > direction * colValFromNonbasicRow) {
        nonbasicRow = colVal.index;
        colValFromNonbasicRow = colValFromRow;
        nonbasicRowStatus = direction * colVal.value > 0
                                ? HighsBasisStatus::kLower
                                : HighsBasisStatus::kUpper;
      }
    }
    // round solution value if column is integer-constrained
    if (nonbasicRow != -1 && colIntegral)
      colValFromNonbasicRow =
          direction * std::ceil(direction * colValFromNonbasicRow -
                                options.mip_feasibility_tolerance);
  };

  if (atInfiniteUpper) {
    // choose largest value as then all rows are feasible
    computeColVal(HighsInt{1});
  } else {
    // choose smallest value, as then all rows are feasible
    computeColVal(HighsInt{-1});
  }
  if (debug_num_use_row_value && debug_report) {
    printf(
        "HighsPostsolveStack::ForcingColumn::undo Using %d unknown row "
        "activit%s\n",
        int(debug_num_use_row_value),
        debug_num_use_row_value > 1 ? "ies" : "y");
  }

  solution.col_value[col] = colValFromNonbasicRow;

  if (!solution.dual_valid) return;

  solution.col_dual[col] = 0.0;

  if (!basis.valid) return;
  if (nonbasicRow == -1) {
    basis.col_status[col] =
        atInfiniteUpper ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
  } else {
    basis.col_status[col] = HighsBasisStatus::kBasic;
    basis.row_status[nonbasicRow] = nonbasicRowStatus;
  }
}

void HighsPostsolveStack::ForcingColumnRemovedRow::undo(
    const HighsOptions& options, const std::vector<Nonzero>& rowValues,
    HighsSolution& solution, HighsBasis& basis) const {
  // we use the row value as storage for the scaled value implied on the
  // column dual
  HighsCDouble val = rhs;
  for (const auto& rowVal : rowValues)
    val -= rowVal.value * solution.col_value[rowVal.index];

  // Row values aren't fully postsolved, so why do this?
  solution.row_value[row] = static_cast<double>(val);

  if (solution.dual_valid) solution.row_dual[row] = 0.0;
  if (basis.valid) basis.row_status[row] = HighsBasisStatus::kBasic;
}

void HighsPostsolveStack::SingletonRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) const {
  // nothing to do if the row's dual value is zero in the dual
  // solution or there is no dual solution
  if (!solution.dual_valid) return;

  const HighsBasisStatus colStatus =
      !basis.valid
          ? computeStatus(solution.col_dual[col],
                          options.dual_feasibility_tolerance)
          : computeStatus(solution.col_dual[col], basis.col_status[col],
                          options.dual_feasibility_tolerance);

  if ((!colLowerTightened || colStatus != HighsBasisStatus::kLower) &&
      (!colUpperTightened || colStatus != HighsBasisStatus::kUpper)) {
    // the tightened bound is not used in the basic solution
    // hence we simply make the row basic and give it a dual multiplier of 0
    if (basis.valid) basis.row_status[row] = HighsBasisStatus::kBasic;
    solution.row_dual[row] = 0;
    return;
  }

  // choose the row dual value such that the columns reduced cost becomes
  // zero
  solution.row_dual[row] = solution.col_dual[col] / coef;
  solution.col_dual[col] = 0;

  if (!basis.valid) return;

  switch (colStatus) {
    case HighsBasisStatus::kLower:
      assert(colLowerTightened);
      if (coef > 0)
        // tightened lower bound comes from row lower bound
        basis.row_status[row] = HighsBasisStatus::kLower;
      else
        // tightened lower bound comes from row upper bound
        basis.row_status[row] = HighsBasisStatus::kUpper;

      break;
    case HighsBasisStatus::kUpper:
      if (coef > 0)
        // tightened upper bound comes from row upper bound
        basis.row_status[row] = HighsBasisStatus::kUpper;
      else
        // tightened upper bound comes from row lower bound
        basis.row_status[row] = HighsBasisStatus::kLower;
      break;
    default:
      assert(false);
  }

  // column becomes basic
  basis.col_status[col] = HighsBasisStatus::kBasic;
}

// column fixed to lower or upper bound
void HighsPostsolveStack::FixedCol::undo(const HighsOptions& options,
                                         const std::vector<Nonzero>& colValues,
                                         HighsSolution& solution,
                                         HighsBasis& basis) const {
  // set solution value
  solution.col_value[col] = fixValue;

  if (!solution.dual_valid) return;

  // compute reduced cost

  HighsCDouble reducedCost = colCost;
  for (const auto& colVal : colValues) {
    reducedCost -= colVal.value * solution.row_dual[colVal.index];
  }

  solution.col_dual[col] = static_cast<double>(reducedCost);

  // set basis status
  if (basis.valid) {
    basis.col_status[col] = fixType;
    if (basis.col_status[col] == HighsBasisStatus::kNonbasic)
      basis.col_status[col] = solution.col_dual[col] >= 0
                                  ? HighsBasisStatus::kLower
                                  : HighsBasisStatus::kUpper;
  }
}

void HighsPostsolveStack::RedundantRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) const {
  // set row dual to zero if dual solution requested
  if (!solution.dual_valid) return;

  solution.row_dual[row] = 0.0;

  if (basis.valid) basis.row_status[row] = HighsBasisStatus::kBasic;
}

void HighsPostsolveStack::ImpliedEquation::undo(
    const std::vector<Nonzero>& rowValues, HighsSolution& solution) const {
  if (!solution.dual_valid) return;
  double oldDual = solution.row_dual[row];
  if (atLower ? (oldDual < 0) : (oldDual > 0)) {
    solution.row_dual[row] = 0;
    for (const auto& nz : rowValues)
      solution.col_dual[nz.index] += nz.value * oldDual;
  }
}

void HighsPostsolveStack::ForcingRow::undo(
    const HighsOptions& options, const std::vector<Nonzero>& rowValues,
    HighsSolution& solution, HighsBasis& basis) const {
  if (!solution.dual_valid) return;

  // compute the row dual multiplier and determine the new basic column
  HighsInt basicCol = -1;
  double dualDelta = 0;
  HighsInt direction = rowType == RowType::kLeq ? 1 : -1;
  for (const auto& rowVal : rowValues) {
    double colDual = solution.col_dual[rowVal.index] - rowVal.value * dualDelta;
    if (direction * colDual * rowVal.value < 0) {
      // column is dual infeasible, decrease the row dual such that its
      // reduced cost become zero and remember this column as the new basic
      // column for this row
      dualDelta = solution.col_dual[rowVal.index] / rowVal.value;
      basicCol = rowVal.index;
    }
  }

  if (basicCol != -1) {
    solution.row_dual[row] = solution.row_dual[row] + dualDelta;
    for (const auto& rowVal : rowValues) {
      solution.col_dual[rowVal.index] = static_cast<double>(
          solution.col_dual[rowVal.index] -
          static_cast<HighsCDouble>(dualDelta) * rowVal.value);
    }
    solution.col_dual[basicCol] = 0;

    if (basis.valid) {
      basis.row_status[row] =
          (rowType == RowType::kGeq ? HighsBasisStatus::kLower
                                    : HighsBasisStatus::kUpper);

      basis.col_status[basicCol] = HighsBasisStatus::kBasic;
    }
  }
}

void HighsPostsolveStack::DuplicateRow::undo(const HighsOptions& options,
                                             HighsSolution& solution,
                                             HighsBasis& basis) const {
  if (!solution.dual_valid) return;
  if (!rowUpperTightened && !rowLowerTightened) {
    // simple case of row2 being redundant, in which case it just gets a
    // dual multiplier of 0 and is made basic
    solution.row_dual[duplicateRow] = 0.0;
    if (basis.valid) basis.row_status[duplicateRow] = HighsBasisStatus::kBasic;
    return;
  }

  const HighsBasisStatus rowStatus =
      !basis.valid
          ? computeStatus(solution.row_dual[row],
                          options.dual_feasibility_tolerance)
          : computeStatus(solution.row_dual[row], basis.row_status[row],
                          options.dual_feasibility_tolerance);

  auto computeRowDualAndStatus = [&](bool tightened, const HighsInt dir) {
    if (tightened) {
      solution.row_dual[duplicateRow] =
          solution.row_dual[row] * duplicateRowScale;
      if (basis.valid) {
        if (dir * duplicateRowScale > 0)
          basis.row_status[duplicateRow] = HighsBasisStatus::kUpper;
        else
          basis.row_status[duplicateRow] = HighsBasisStatus::kLower;
      }
      solution.row_dual[row] = 0.0;
      if (basis.valid) basis.row_status[row] = HighsBasisStatus::kBasic;
    } else {
      solution.row_dual[duplicateRow] = 0.0;
      if (basis.valid)
        basis.row_status[duplicateRow] = HighsBasisStatus::kBasic;
    }
  };

  // at least one bound of the row was tightened by using the bound of the
  // scaled parallel row, hence we might need to make the parallel row
  // nonbasic and the row basic

  switch (rowStatus) {
    case HighsBasisStatus::kBasic:
      // if row is basic the parallel row is also basic
      solution.row_dual[duplicateRow] = 0.0;
      if (basis.valid)
        basis.row_status[duplicateRow] = HighsBasisStatus::kBasic;
      break;
    case HighsBasisStatus::kUpper:
      // if row sits on its upper bound, and the row upper bound was
      // tightened using the parallel row we make the row basic and
      // transfer its dual value to the parallel row with the proper scale
      computeRowDualAndStatus(rowUpperTightened, 1);
      break;
    case HighsBasisStatus::kLower:
      computeRowDualAndStatus(rowLowerTightened, -1);
      break;
    default:
      assert(false);
  }
}

void HighsPostsolveStack::DuplicateColumn::undo(const HighsOptions& options,
                                                HighsSolution& solution,
                                                HighsBasis& basis) const {
  const bool debug_report = false;
  const double mergeVal = solution.col_value[col];

  auto okResidual = [&](const double x, const double y) {
    const double check_mergeVal = x + colScale * y;
    const double residual = std::fabs(check_mergeVal - mergeVal);
    const bool ok_residual = residual <= options.primal_feasibility_tolerance;
    if (!ok_residual && debug_report) {
      printf(
          "HighsPostsolveStack::DuplicateColumn::undo %g + %g.%g = %g != %g: "
          "residual = %g\n",
          x, colScale, y, check_mergeVal, mergeVal, residual);
    }
    return ok_residual;
  };

  auto isAtBound = [&](const double value, const double bound) {
    if (value < bound - options.primal_feasibility_tolerance) return false;
    if (value <= bound + options.primal_feasibility_tolerance) return true;
    return false;
  };

  //  const bool ok_merge = okMerge(options.mip_feasibility_tolerance);
  //  assert(ok_merge);
  //
  // the column dual of the duplicate column is easily computed by scaling
  // since col * colScale yields the coefficient values and cost of the
  // duplicate column.
  if (solution.dual_valid)
    solution.col_dual[duplicateCol] = solution.col_dual[col] * colScale;

  if (basis.valid) {
    // do postsolve using basis status if a basis is available: if the
    // merged column is nonbasic, we can just set both columns to
    // appropriate nonbasic status and value
    //
    // Undoing z = x + a.y
    //
    // Since x became z, its basis status is unchanged
    //
    // For a > 0, z\in [x_l + a.y_l, x_u + a.y_u]
    //
    // If z is nonbasic at its lower (upper) bound, set y to be
    // nonbasic at its lower (upper) bound
    //
    // For a < 0, z\in [x_l + a.y_u, x_u + a.y_l]
    //
    // If z is nonbasic at lower (upper) bound, set y to be nonbasic
    // at its upper (lower) bounds
    //
    // Check for perturbations
    switch (basis.col_status[col]) {
      case HighsBasisStatus::kLower: {
        solution.col_value[col] = colLower;
        if (colScale > 0) {
          basis.col_status[duplicateCol] = HighsBasisStatus::kLower;
          solution.col_value[duplicateCol] = duplicateColLower;
        } else {
          basis.col_status[duplicateCol] = HighsBasisStatus::kUpper;
          solution.col_value[duplicateCol] = duplicateColUpper;
        }
        // nothing else to do
        assert(okResidual(solution.col_value[col],
                          solution.col_value[duplicateCol]));
        return;
      }
      case HighsBasisStatus::kUpper: {
        solution.col_value[col] = colUpper;
        if (colScale > 0) {
          basis.col_status[duplicateCol] = HighsBasisStatus::kUpper;
          solution.col_value[duplicateCol] = duplicateColUpper;
        } else {
          basis.col_status[duplicateCol] = HighsBasisStatus::kLower;
          solution.col_value[duplicateCol] = duplicateColLower;
        }
        // nothing else to do
        assert(okResidual(solution.col_value[col],
                          solution.col_value[duplicateCol]));
        return;
      }
      case HighsBasisStatus::kZero: {
        solution.col_value[col] = 0.0;
        basis.col_status[duplicateCol] = HighsBasisStatus::kZero;
        solution.col_value[duplicateCol] = 0.0;
        // nothing else to do
        assert(okResidual(solution.col_value[col],
                          solution.col_value[duplicateCol]));
        return;
      }
      case HighsBasisStatus::kBasic:
      case HighsBasisStatus::kNonbasic:;
    }
    // Nonbasic cases should have been considered; basic case
    // considered later
    assert(basis.col_status[col] == HighsBasisStatus::kBasic);
  }

  // There is either no basis for postsolve, or column status is
  // basic. One of the two columns must become nonbasic. In case of
  // integrality it is simpler to choose col, since it has a
  // coefficient of +1 in the equation z = col + colScale *
  // duplicateCol where the merged column is z and is currently using
  // the index of col. The duplicateCol can have a positive or
  // negative coefficient. So for postsolve, we first start out with
  // col sitting at the lower bound and compute the corresponding
  // value for the duplicate column as (z - colLower)/colScale. Then
  // the following things might happen:
  //
  // - case 1: the value computed for duplicateCol is within the bounds
  // - case 1.1: duplicateCol is continuous -> accept value, make col nonbasic
  // at lower and duplicateCol basic
  // - case 1.2: duplicateCol is integer -> accept value if integer feasible,
  // otherwise round down and compute value of col as
  // col = z - colScale * duplicateCol
  // - case 2: the value for duplicateCol violates the column bounds: make it
  // sit at the bound that is violated
  //           and compute the value of col as col = z - colScale *
  //           duplicateCol for basis postsolve col is basic and duplicateCol
  //           nonbasic at lower/upper depending on which bound is violated.

  if (colLower != -kHighsInf)
    solution.col_value[col] = colLower;
  else
    solution.col_value[col] = std::min(0.0, colUpper);
  // Determine a prospective value for col_value[duplicateCol]
  solution.col_value[duplicateCol] = static_cast<double>(
      (static_cast<HighsCDouble>(mergeVal) - solution.col_value[col]) /
      colScale);

  bool recomputeCol = false;
  bool duplicateColAtUpper = false;

  // Set any basis status for duplicateCol to kNonbasic to check that
  // it is set
  if (basis.valid) basis.col_status[duplicateCol] = HighsBasisStatus::kNonbasic;

  if (solution.col_value[duplicateCol] > duplicateColUpper) {
    // Prospective value exceeds the upper bound, so trim it to the
    // upper bound and force recalculation of col_value[col]
    solution.col_value[duplicateCol] = duplicateColUpper;
    recomputeCol = true;
    duplicateColAtUpper = true;
    if (basis.valid) basis.col_status[duplicateCol] = HighsBasisStatus::kUpper;
  } else if (solution.col_value[duplicateCol] < duplicateColLower) {
    // Prospective value exceeds the lower bound, so trim it to the
    // lower bound and force recalculation of col_value[col]
    solution.col_value[duplicateCol] = duplicateColLower;
    recomputeCol = true;
    if (basis.valid) basis.col_status[duplicateCol] = HighsBasisStatus::kLower;
  } else if (duplicateColIntegral) {
    assert(!basis.valid);
    if (fractionality(solution.col_value[duplicateCol]) >
        options.mip_feasibility_tolerance) {
      // Prospective value is not integer, so trim it to the integer
      // value and force recalculation of col_value[col]
      solution.col_value[duplicateCol] =
          std::floor(solution.col_value[duplicateCol]);
      recomputeCol = true;
    }
  }

  if (recomputeCol) {
    // Determine the value of col_value[col] corresponding to
    // col_value[duplicateCol]
    solution.col_value[col] =
        mergeVal - colScale * solution.col_value[duplicateCol];
    if (!duplicateColIntegral && colIntegral) {
      // If column is integral and duplicateCol is not we need to make sure
      // we split the values into an integral one for col
      assert(!basis.valid);
      solution.col_value[col] =
          duplicateColAtUpper == (colScale > 0)
              ? std::ceil(solution.col_value[col] -
                          options.mip_feasibility_tolerance)
              : std::floor(solution.col_value[col] +
                           options.mip_feasibility_tolerance);
      solution.col_value[duplicateCol] = static_cast<double>(
          (static_cast<HighsCDouble>(mergeVal) - solution.col_value[col]) /
          colScale);
    }
  } else {
    // Have set col_value[col] to its lower bound (if finite),
    // otherwise min(0, upper), and if there is a basis, then
    // col_status[col] is kBasic.
    if (basis.valid) {
      // Now set the basis status of col and duplicateCol.
      //
      // Set col_status[col] to be nonbasic corresponding to its
      // value. This is fine if it's at its lower or upper bound, or
      // at zero and free, but what if it's at zero with a finite
      // upper bound?
      //
      // The basis is maintained by setting duplicateCol to be basic
      assert(basis.col_status[col] == HighsBasisStatus::kBasic);
      basis.col_status[duplicateCol] = basis.col_status[col];
      if (colLower != -kHighsInf) {
        assert(solution.col_value[col] == colLower);
        basis.col_status[col] = HighsBasisStatus::kLower;
      } else if (colUpper <= 0.0) {
        assert(solution.col_value[col] == colUpper);
        basis.col_status[col] = HighsBasisStatus::kUpper;
      } else {
        assert(solution.col_value[col] == 0.0);
        basis.col_status[col] = HighsBasisStatus::kZero;
        if (colUpper < kHighsInf) {
          // Nonbasic at zero with bounds (-inf, colUpper)
          highsLogUser(
              options.log_options, HighsLogType::kWarning,
              "HighsPostsolveStack::DuplicateColumn::undo Col is nonbasic at "
              "zero with upper bound of %g\n",
              colUpper);
        }
      }
      assert(basis.col_status[duplicateCol] == HighsBasisStatus::kBasic);
    }
  }
  // Check that any basis status for duplicateCol has been set
  if (basis.valid)
    assert(basis.col_status[duplicateCol] != HighsBasisStatus::kNonbasic);

  auto hasError = [&]() {
    bool illegal_duplicateCol_lower =
        solution.col_value[duplicateCol] <
        duplicateColLower - options.mip_feasibility_tolerance;
    bool illegal_duplicateCol_upper =
        solution.col_value[duplicateCol] >
        duplicateColUpper + options.mip_feasibility_tolerance;
    bool illegal_col_lower =
        solution.col_value[col] < colLower - options.mip_feasibility_tolerance;
    bool illegal_col_upper =
        solution.col_value[col] > colUpper + options.mip_feasibility_tolerance;
    bool illegal_residual =
        !okResidual(solution.col_value[col], solution.col_value[duplicateCol]);
    return (illegal_duplicateCol_lower || illegal_duplicateCol_upper ||
            illegal_col_lower || illegal_col_upper || illegal_residual);
  };

  if (hasError()) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: col = %d(%g), duplicateCol = %d(%g)\n"
          "%g\n%g\n%g %g %d\n%g %g %d\n",
          int(col), solution.col_value[col], int(duplicateCol),
          solution.col_value[duplicateCol], mergeVal, colScale, colLower,
          colUpper, colIntegral, duplicateColLower, duplicateColUpper,
          duplicateColIntegral);
    // Fix error due to undo
    undoFix(options, solution, mergeVal);
  } else {
    return;
  }
  const bool allow_assert = false;
  if (allow_assert) assert(!hasError());
  // Following undoFix, set any basis status, ideally keeping col basic
  if (basis.valid) {
    bool duplicateCol_basic = false;
    if (duplicateColLower <= -kHighsInf && duplicateColUpper >= kHighsInf) {
      // duplicateCol is free, so may be zero
      if (solution.col_value[duplicateCol] == 0) {
        basis.col_status[col] = HighsBasisStatus::kBasic;
        basis.col_status[duplicateCol] = HighsBasisStatus::kZero;
      } else {
        duplicateCol_basic = true;
      }
    } else if (isAtBound(solution.col_value[duplicateCol], duplicateColLower)) {
      basis.col_status[col] = HighsBasisStatus::kBasic;
      basis.col_status[duplicateCol] = HighsBasisStatus::kLower;
    } else if (isAtBound(solution.col_value[duplicateCol], duplicateColUpper)) {
      basis.col_status[col] = HighsBasisStatus::kBasic;
      basis.col_status[duplicateCol] = HighsBasisStatus::kUpper;
    } else {
      // duplicateCol is not free or at a bound, so must be basic
      duplicateCol_basic = true;
    }
    if (duplicateCol_basic) {
      // duplicateCol must be basic
      basis.col_status[duplicateCol] = HighsBasisStatus::kBasic;
      // Hopefully col can be nonbasic
      if (isAtBound(solution.col_value[col], colLower)) {
        basis.col_status[col] = HighsBasisStatus::kLower;
      } else if (isAtBound(solution.col_value[col], colUpper)) {
        basis.col_status[col] = HighsBasisStatus::kUpper;
      } else {
        basis.col_status[col] = HighsBasisStatus::kNonbasic;
        if (debug_report)
          printf(
              "When demerging, neither col nor duplicateCol can be "
              "nonbasic\n");
        if (kAllowDeveloperAssert) assert(666 == 999);
      }
    }
  }
}

bool HighsPostsolveStack::DuplicateColumn::okMerge(
    const double tolerance) const {
  // When merging x and y to x+a.y, not all values of a are permitted,
  // since it must be possible to map back onto feasible values of x
  // and y.
  //
  // Assume WLOG that a > 0, x\in[x_l, x_u], y\in[y_l, y_u]
  //
  // Let z = x + a.y
  //
  // Range for z is [x_l+a.y_l, x_u+a.y_u]
  //
  // * If x and y are both integer:
  //
  // z will be integer and x+a.y must generate all integer values in
  // [x_l+a.y_l, x_u+a.y_u]. Hence a must be an integer. If a >=
  // (x_u-x_l)+2 then, since [a.y_l, a.y_u] contains integer multiples
  // of a, some of the intervening integers don't correspond to a
  // value of x. Hence a must be an integer and a <= (x_u-x_l)+1
  //
  // For example, if x and y are binary, then x+a.y is [0, 1, a,
  // 1+a]. For this to be a continuous sequence of integers, we must
  // have a <= 2.
  //
  // * If x is integer and y is continuous:
  //
  // z will be continuous and x+a.y must generate all values in
  // [x_l+a.y_l, x_u+a.y_u]. Since [x_l, x_u] are integers, [a.y_l,
  // a.y_u] = a[y_l, y_u] must be of length at least 1. Hence a must
  // be at least 1/(y_u-y_l) in magnitude.
  //
  // * If x is continuous and y is integer:
  //
  // z will be continuous and x+a.y must generate all values in
  // [x_l+a.y_l, x_u+a.y_u]. Since [a.y_l, a.y_u] contains integer
  // multiples of a, the gaps between them must not exceed the length
  // of [x_l, x_u]. Hence a must be at most x_u-x_l in
  // magnitude.
  //
  // Observe that this is equivalent to requiring 1/a to be at least
  // 1/(x_u-x_l) in magnitude, the symmetric result corresponding to
  // the merge (1/a)x+y.
  //
  //  * If x and y are both continuous
  //
  // z will be continuous and x+a.y naturally generates all values in
  // [x_l+a.y_l, x_u+a.y_u].

  const double scale = colScale;
  const bool x_int = colIntegral;
  const bool y_int = duplicateColIntegral;
  const double x_lo = x_int ? std::ceil(colLower - tolerance) : colLower;
  const double x_up = x_int ? std::floor(colUpper + tolerance) : colUpper;
  const double y_lo =
      y_int ? std::ceil(duplicateColLower - tolerance) : duplicateColLower;
  const double y_up =
      y_int ? std::floor(duplicateColUpper + tolerance) : duplicateColUpper;
  const double x_len = x_up - x_lo;
  const double y_len = y_up - y_lo;
  std::string newline = "\n";
  bool ok_merge = true;
  const bool debug_report = false;
  if (scale == 0) {
    if (debug_report)
      printf("%sDuplicateColumn::checkMerge: Scale cannot be zero\n",
             newline.c_str());
    newline = "";
    ok_merge = false;
  }
  const double abs_scale = std::fabs(scale);
  if (x_int) {
    if (y_int) {
      // Scale must be integer and not exceed (x_u-x_l)+1 in magnitude
      if (fractionality(scale) > tolerance) {
        if (debug_report)
          printf(
              "%sDuplicateColumn::checkMerge: scale must be integer, but is "
              "%g\n",
              newline.c_str(), scale);
        newline = "";
        ok_merge = false;
      }
      double scale_limit = x_len + 1 + tolerance;
      if (abs_scale > scale_limit) {
        if (debug_report)
          printf(
              "%sDuplicateColumn::checkMerge: scale = %g, but |scale| cannot "
              "exceed %g since x is [%g, %g]\n",
              newline.c_str(), scale, scale_limit, x_lo, x_up);
        newline = "";
        ok_merge = false;
      }
    } else {  // y is continuous
      if (debug_report)
        printf("DuplicateColumn::checkMerge: x-integer; y-continuous\n");
      // Scale must be at least 1/(y_u-y_l) in magnitude
      if (y_len == 0) {
        if (debug_report)
          printf(
              "%sDuplicateColumn::checkMerge: scale = %g is too small in "
              "magnitude, as y is [%g, %g]\n",
              newline.c_str(), scale, y_lo, y_up);
        newline = "";
        ok_merge = false;
      } else {
        double scale_limit = 1 / y_len;
        if (abs_scale < scale_limit) {
          if (debug_report)
            printf(
                "%sDuplicateColumn::checkMerge: scale = %g, but |scale| must "
                "be "
                "at least %g since y is [%g, %g]\n",
                newline.c_str(), scale, scale_limit, y_lo, y_up);
          newline = "";
          ok_merge = false;
        }
      }
    }
  } else {
    if (y_int) {
      if (debug_report)
        printf("DuplicateColumn::checkMerge: x-continuous; y-integer\n");
      // Scale must be at most (x_u-x_l) in magnitude
      double scale_limit = x_len;
      if (abs_scale > scale_limit) {
        if (debug_report)
          printf(
              "%sDuplicateColumn::checkMerge: scale = %g, but |scale| must "
              "be "
              "at "
              "most %g since x is [%g, %g]\n",
              newline.c_str(), scale, scale_limit, x_lo, x_up);
        newline = "";
        ok_merge = false;
      }
    } else {
      // x and y are continuous
      //	if (debug_report) printf("DuplicateColumn::checkMerge:
      // x-continuous ;
      // y-continuous\n");
    }
  }
  return ok_merge;
}

void HighsPostsolveStack::DuplicateColumn::undoFix(
    const HighsOptions& options, HighsSolution& solution,
    const double mergeValue) const {
  const double mip_feasibility_tolerance = options.mip_feasibility_tolerance;
  const double primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  std::vector<double>& col_value = solution.col_value;
  const bool allow_assert = false;
  const bool debug_report = false;
  //=============================================================================================

  auto isInteger = [&](const double v) {
    return (fractionality(v) <= mip_feasibility_tolerance);
  };

  auto isFeasible = [&](const double l, const double v, const double u) {
    return v >= l - primal_feasibility_tolerance &&
           v <= u + primal_feasibility_tolerance;
  };
  const double scale = colScale;
  const bool x_int = colIntegral;
  const bool y_int = duplicateColIntegral;
  const int x_ix = col;
  const int y_ix = duplicateCol;
  const double x_lo =
      x_int ? std::ceil(colLower - mip_feasibility_tolerance) : colLower;
  const double x_up =
      x_int ? std::floor(colUpper + mip_feasibility_tolerance) : colUpper;
  const double y_lo =
      y_int ? std::ceil(duplicateColLower - mip_feasibility_tolerance)
            : duplicateColLower;
  const double y_up =
      y_int ? std::floor(duplicateColUpper + mip_feasibility_tolerance)
            : duplicateColUpper;
  if (kAllowDeveloperAssert) assert(scale);
  double x_v = mergeValue;
  double y_v = -kHighsInf;

  auto checkIntVar = [](double z_lo, double z_up, double& z_0, double& z_d,
                        double& z_1, bool& z_free) {
    const double value_max = 1000;
    z_0 = 0;
    z_d = 0;
    z_1 = 0;
    z_free = false;
    if (z_lo <= -kHighsInf) {
      if (z_up >= kHighsInf) {
        // z is free
        z_free = true;
        z_0 = 0;
        z_d = 1.0;
        z_1 = value_max;
      } else {
        // z is (-int, u]
        z_0 = z_up;
        z_d = -1.0;
        z_1 = -value_max;
      }
    } else {
      // z is [l, u] or [l, inf) respectively
      z_0 = z_lo;
      z_d = 1.0;
      z_1 = z_up >= kHighsInf ? value_max : z_up;
    }
  };

  auto computeValue = [mergeValue, scale](double value) {
    return static_cast<double>((static_cast<HighsCDouble>(mergeValue) - value) /
                               scale);
  };

  auto computeInvValue = [mergeValue, scale](double value) {
    return static_cast<double>(static_cast<HighsCDouble>(mergeValue) -
                               static_cast<HighsCDouble>(value) * scale);
  };

  auto findValue = [&](double& z_value, double z_0, double z_1, double z_delta,
                       double& other_value, double other_lower,
                       double other_upper, bool other_int) {
    const double eps = 1e-8;
    for (z_value = z_0;; z_value += z_delta) {
      other_value = computeValue(z_value);
      if (isFeasible(other_lower, other_value, other_upper) &&
          (!other_int || isInteger(other_value)))
        return true;
      if (z_delta > 0 && z_value + z_delta >= z_1 + eps) return false;
      if (z_delta < 0 && z_value + z_delta <= z_1 - eps) return false;
    }
    return false;
  };

  auto setXValue = [&](double value) {
    x_v = value;
    y_v = computeValue(value);
  };

  auto setYValue = [&](double value) {
    y_v = value;
    x_v = computeInvValue(value);
  };

  if (x_int) {
    double x_0;
    double x_d;
    double x_1;
    bool x_free;
    checkIntVar(x_lo, x_up, x_0, x_d, x_1, x_free);
    // x is integer, so look through its possible values to find a
    // suitable y
    if (x_free && debug_report) printf("DuplicateColumn::undo x is free\n");
    if (debug_report)
      printf("DuplicateColumn::undo Using x (%g; %g; %g)\n", x_0, x_d, x_1);
    bool found_y = findValue(x_v, x_0, x_1, x_d, y_v, y_lo, y_up, y_int);
    if (allow_assert) assert(found_y);
  } else if (y_int) {
    double y_0;
    double y_d;
    double y_1;
    bool y_free;
    checkIntVar(y_lo, y_up, y_0, y_d, y_1, y_free);
    // y is integer, so look through its possible values to find a
    // suitable x
    if (y_free && debug_report) printf("DuplicateColumn::undo y is free\n");
    if (debug_report)
      printf("DuplicateColumn::undo Using y (%g; %g; %g)\n", y_0, y_d, y_1);
    bool found_x = findValue(y_v, y_0, y_1, y_d, x_v, x_lo, x_up, x_int);
    if (allow_assert) assert(found_x);
  } else {
    // x and y are both continuous
    double v_m_a_ylo = 0;
    double v_m_a_yup = 0;
    if (y_lo <= -kHighsInf) {
      v_m_a_ylo = scale > 0 ? kHighsInf : -kHighsInf;
    } else {
      v_m_a_ylo = computeInvValue(y_lo);
    }
    if (y_up >= kHighsInf) {
      v_m_a_yup = scale > 0 ? -kHighsInf : kHighsInf;
    } else {
      v_m_a_yup = computeInvValue(y_up);
    }
    // Need to ensure that y puts x in [x_l, x_u]
    if (scale > 0) {
      if (debug_report)
        printf("DuplicateColumn::undo [V-a(y_u), V-a(y_l)] == [%g, %g]\n",
               v_m_a_yup, v_m_a_ylo);
      // V-ay is in [V-a(y_u), V-a(y_l)] == [v_m_a_yup, v_m_a_ylo]
      if (y_up < kHighsInf) {
        // If v_m_a_yup is right of x_up+eps then [v_m_a_yup, v_m_a_ylo] is
        // right of [x_lo-eps, x_up+eps] so there's no solution. [Could
        // try v_m_a_ylo computed from y_lo-eps.]
        if (kAllowDeveloperAssert)
          assert(x_up + primal_feasibility_tolerance >= v_m_a_yup);
        // This assignment is OK unless x_v < x_lo-eps
        setYValue(y_up);
        if (x_v < x_lo - primal_feasibility_tolerance) {
          // Try y_v corresponding to x_lo
          setXValue(x_lo);
          if (y_v < y_lo - primal_feasibility_tolerance) {
            // Very tight: use x_v on its margin and hope!
            setXValue(x_lo - primal_feasibility_tolerance);
          }
        }
      } else if (y_lo > -kHighsInf) {
        // If v_m_a_ylo is left of x_lo-eps then [v_m_a_yup, v_m_a_ylo] is
        // left of [x_lo-eps, x_up+eps] so there's no solution. [Could
        // try v_m_a_yup computed from y_up+eps.]
        if (kAllowDeveloperAssert)
          assert(x_lo - primal_feasibility_tolerance <= v_m_a_ylo);
        // This assignment is OK unless x_v > x_up-eps
        setYValue(y_lo);
        if (x_v > x_up + primal_feasibility_tolerance) {
          // Try y_v corresponding to x_up
          setXValue(x_up);
          if (y_v > y_up + primal_feasibility_tolerance) {
            // Very tight: use x_v on its margin and hope!
            setXValue(x_up + primal_feasibility_tolerance);
          }
        }
      } else {
        // y is free, so use x_v = max(0, x_lo)
        setXValue(std::max(0.0, x_lo));
      }
    } else {  // scale < 0
      if (debug_report)
        printf("DuplicateColumn::undo [V-a(y_l), V-a(y_u)] == [%g, %g]\n",
               v_m_a_ylo, v_m_a_yup);
      // V-ay is in [V-a(y_l), V-a(y_u)] == [v_m_a_ylo, v_m_a_yup]
      if (y_lo > -kHighsInf) {
        // If v_m_a_ylo is right of x_up+eps then [v_m_a_ylo, v_m_a_yup] is
        // right of [x_lo-eps, x_up+eps] so there's no solution. [Could
        // try v_m_a_ylo computed from y_up+eps.]
        if (kAllowDeveloperAssert)
          assert(x_up + primal_feasibility_tolerance >= v_m_a_ylo);
        // This assignment is OK unless x_v < x_lo-eps
        setYValue(y_lo);
        if (x_v < x_lo - primal_feasibility_tolerance) {
          // Try y_v corresponding to x_lo
          setXValue(x_lo);
          if (y_v > y_up + primal_feasibility_tolerance) {
            // Very tight: use x_v on its margin and hope!
            setXValue(x_lo - primal_feasibility_tolerance);
          }
        }
      } else if (y_up < kHighsInf) {
        // If v_m_a_yup is left of x_lo-eps then [v_m_a_ylo, v_m_a_yup] is
        // left of [x_lo-eps, x_up+eps] so there's no solution. [Could
        // try v_m_a_yup computed from y_lo-eps.]
        if (kAllowDeveloperAssert)
          assert(x_lo - primal_feasibility_tolerance <= v_m_a_yup);
        // This assignment is OK unless x_v < x_lo-eps
        setYValue(y_up);
        if (x_v > x_up + primal_feasibility_tolerance) {
          // Try y_v corresponding to x_up
          setXValue(x_up);
          if (y_v < y_lo - primal_feasibility_tolerance) {
            // Very tight: use x_v on its margin and hope!
            setXValue(x_up + primal_feasibility_tolerance);
          }
        }
      } else {
        // y is free, so use x_v = max(0, x_lo)
        setXValue(std::max(0.0, x_lo));
      }
    }
  }
  const double residual_tolerance = 1e-12;
  double residual = std::fabs(static_cast<double>(
      static_cast<HighsCDouble>(x_v) - computeInvValue(y_v)));
  const bool x_y_ok =
      isFeasible(x_lo, x_v, x_up) && isFeasible(y_lo, y_v, y_up) &&
      (!x_int || isInteger(x_v)) && (!y_int || isInteger(y_v)) &&
      (std::fabs(x_v) < kHighsInf) && (std::fabs(y_v) < kHighsInf) &&
      (residual <= residual_tolerance);

  bool check;
  check = isFeasible(x_lo, x_v, x_up);
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: isFeasible(x_lo, x_v, x_up) is "
          "false\n");
    if (allow_assert) assert(check);
  }
  check = isFeasible(y_lo, y_v, y_up);
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: isFeasible(y_lo, y_v, y_up) is "
          "false\n");
    if (allow_assert) assert(check);
  }
  check = !x_int || isInteger(x_v);
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: !x_int || isInteger(x_v) is false\n");
    if (allow_assert) assert(check);
  }
  check = !y_int || isInteger(y_v);
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: !y_int || isInteger(y_v) is false\n");
    if (allow_assert) assert(check);
  }
  check = std::fabs(x_v) < kHighsInf;
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: std::fabs(x_v) < kHighsInf is "
          "false\n");
    if (allow_assert) assert(check);
  }
  check = std::fabs(y_v) < kHighsInf;
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: std::fabs(y_v) < kHighsInf is "
          "false\n");
    if (allow_assert) assert(check);
  }
  check = residual <= residual_tolerance;
  if (!check) {
    if (debug_report)
      printf(
          "DuplicateColumn::undo error: residual <= residual_tolerance is "
          "false\n");
    if (allow_assert) assert(check);
  }
  check = residual <= residual_tolerance;
  if (debug_report)
    printf("DuplicateColumn::undo%s x = %g; y = %g to give x + (%g)y = %g",
           x_y_ok ? "" : " ERROR", x_v, y_v, scale, mergeValue);
  if (x_y_ok) {
    if (debug_report) printf(": FIXED\n");
  } else if (check) {
    if (debug_report) printf("\n");
  } else {
    if (debug_report) printf(": residual = %g\n", residual);
  }
  //=============================================================================================
  if (x_y_ok) {
    col_value[x_ix] = x_v;
    col_value[y_ix] = y_v;
  }
}

void HighsPostsolveStack::DuplicateColumn::transformToPresolvedSpace(
    std::vector<double>& primalSol) const {
  primalSol[col] = primalSol[col] + colScale * primalSol[duplicateCol];
}

void HighsPostsolveStack::FourierMotzkinObjCol::transformToPresolvedSpace(
    const std::vector<Nonzero>& costEntries,
    std::vector<double>& primalSol) const {
  double val = offset;
  for (const Nonzero& entry : costEntries)
    val += entry.value * primalSol[entry.index];
  primalSol[col] = val;
}

void HighsPostsolveStack::FourierMotzkinObjCol::undo(
    const std::vector<Nonzero>& costEntries, HighsSolution& solution) const {
  if (!solution.dual_valid) return;
  double zDual = solution.col_dual[col];
  for (const Nonzero& entry : costEntries)
    solution.col_dual[entry.index] += entry.value * zDual;
  solution.col_dual[col] = 0.0;
}

std::vector<HighsPostsolveStack::FmeStepData>
HighsPostsolveStack::popFourierMotzkinBlock(HighsDataStack& stack) {
  HighsInt numSteps;
  stack.pop(numSteps);

  std::vector<FmeStepData> steps(numSteps);

  // step headers
  for (HighsInt s = numSteps - 1; s >= 0; --s) stack.pop(steps[s].header);

  // new row origins
  for (HighsInt s = numSteps - 1; s >= 0; --s) stack.pop(steps[s].newRows);

  // descendants
  for (HighsInt s = numSteps - 1; s >= 0; --s) {
    HighsInt numMinus = steps[s].header.numMinus;
    steps[s].minusRows.resize(numMinus);
    for (HighsInt m = numMinus - 1; m >= 0; --m)
      stack.pop(steps[s].minusRows[m].descendants);
    HighsInt numPlus = steps[s].header.numPlus;
    steps[s].plusRows.resize(numPlus);
    for (HighsInt p = numPlus - 1; p >= 0; --p)
      stack.pop(steps[s].plusRows[p].descendants);
  }

  // row data
  for (HighsInt s = numSteps - 1; s >= 0; --s) {
    // minus row data
    std::vector<FmeRowHeader> minusHeaders;
    std::vector<double> minusCoefs;
    stack.pop(minusHeaders);
    stack.pop(minusCoefs);
    HighsInt numMinus = static_cast<HighsInt>(minusCoefs.size());
    for (HighsInt r = 0; r < numMinus; ++r) {
      steps[s].minusRows[r].header = minusHeaders[r];
      steps[s].minusRows[r].coef = minusCoefs[r];
    }
    for (HighsInt r = numMinus - 1; r >= 0; --r)
      stack.pop(steps[s].minusRows[r].entries);

    // plus row data
    std::vector<FmeRowHeader> plusHeaders;
    std::vector<double> plusCoefs;
    stack.pop(plusHeaders);
    stack.pop(plusCoefs);
    HighsInt numPlus = static_cast<HighsInt>(plusCoefs.size());
    for (HighsInt r = 0; r < numPlus; ++r) {
      steps[s].plusRows[r].header = plusHeaders[r];
      steps[s].plusRows[r].coef = plusCoefs[r];
    }
    for (HighsInt r = numPlus - 1; r >= 0; --r)
      stack.pop(steps[s].plusRows[r].entries);
  }

  return steps;
}

void HighsPostsolveStack::undoFourierMotzkinBlock(
    const std::vector<FmeStepData>& steps, const HighsOptions& options,
    HighsSolution& solution, HighsBasis& basis) {
  const double tol = options.primal_feasibility_tolerance;
  const double dual_tol = options.dual_feasibility_tolerance;

  HighsInt numSteps = static_cast<HighsInt>(steps.size());

  // primal postsolve (Algorithm 3): process in reverse elimination order
  for (HighsInt s = numSteps - 1; s >= 0; --s) {
    const auto& step = steps[s];
    HighsInt col = step.header.col;
    double lower = step.header.colLower;
    double upper = step.header.colUpper;

    auto tightenBounds = [&](const std::vector<FmeParentRow>& rows,
                             double& lowerBound, double& upperBound) {
      for (const auto& row : rows) {
        double aij = row.coef;
        HighsCDouble sum = 0.0;
        for (const auto& nz : row.entries)
          sum += static_cast<HighsCDouble>(nz.value) *
                 solution.col_value[nz.index];
        HighsInt direction = aij > 0 ? HighsInt{1} : HighsInt{-1};
        double rhs_upper =
            direction > 0 ? row.header.rowUpper : row.header.rowLower;
        double rhs_lower =
            direction > 0 ? row.header.rowLower : row.header.rowUpper;
        if (direction * rhs_upper != kHighsInf) {
          double bound = static_cast<double>((rhs_upper - sum) / aij);
          upperBound = std::min(upperBound, bound);
        }
        if (direction * rhs_lower != -kHighsInf) {
          double bound = static_cast<double>((rhs_lower - sum) / aij);
          lowerBound = std::max(lowerBound, bound);
        }
      }
    };

    tightenBounds(step.plusRows, lower, upper);
    tightenBounds(step.minusRows, lower, upper);

    if (lower <= tol && upper >= -tol)
      solution.col_value[col] = 0.0;
    else if (lower > 0.0)
      solution.col_value[col] = lower;
    else
      solution.col_value[col] = upper;
  }

  if (!solution.dual_valid) return;

  // dual postsolve (Algorithm 4): process in reverse elimination order
  for (HighsInt s = numSteps - 1; s >= 0; --s) {
    const auto& step = steps[s];
    HighsInt col = step.header.col;
    HighsInt numPlus = step.header.numPlus;
    HighsInt numMinus = step.header.numMinus;

    // u_i = Σ_{k ∈ K^j_i} λ_k * scaleFactor
    auto recoverDual = [&](const std::vector<FmeParentRow>& rows) {
      for (const auto& row : rows) {
        HighsCDouble dual = 0.0;
        for (const auto& desc : row.descendants)
          dual += static_cast<HighsCDouble>(solution.row_dual[desc.row]) *
                  desc.scaleFactor;
        solution.row_dual[row.header.row] += static_cast<double>(dual);
      }
    };
    recoverDual(step.plusRows);
    recoverDual(step.minusRows);

    // col_dual = -Σ a_{ij} * row_dual[i] (cost is zero after reformulation)
    HighsCDouble colDual = 0.0;
    std::vector<HighsBool> visited(solution.row_dual.size(), false);
    for (const auto& row : step.plusRows) {
      colDual -= static_cast<HighsCDouble>(row.coef) *
                 solution.row_dual[row.header.row];
      visited[row.header.row] = true;
    }
    for (const auto& row : step.minusRows) {
      if (visited[row.header.row]) continue;
      colDual -= static_cast<HighsCDouble>(row.coef) *
                 solution.row_dual[row.header.row];
    }
    solution.col_dual[col] = static_cast<double>(colDual);
  }

  // basis postsolve: use dual solution to determine basis status
  if (!basis.valid) return;

  // pre-compute lower and upper slacks for each row
  auto computeSlacks = [&](HighsInt col, const std::vector<FmeParentRow>& rows,
                           std::vector<double>& lowerSlacks,
                           std::vector<double>& upperSlacks) {
    HighsInt n = static_cast<HighsInt>(rows.size());
    lowerSlacks.resize(n);
    upperSlacks.resize(n);
    for (HighsInt r = 0; r < n; ++r) {
      HighsCDouble activity =
          static_cast<HighsCDouble>(rows[r].coef) * solution.col_value[col];
      for (const auto& nz : rows[r].entries)
        activity +=
            static_cast<HighsCDouble>(nz.value) * solution.col_value[nz.index];
      double act = static_cast<double>(activity);
      lowerSlacks[r] = rows[r].header.rowLower != -kHighsInf
                           ? act - rows[r].header.rowLower
                           : kHighsInf;
      upperSlacks[r] = rows[r].header.rowUpper != kHighsInf
                           ? rows[r].header.rowUpper - act
                           : kHighsInf;
    }
  };

  // row must be basic if it has zero dual and activity strictly
  // between bounds (complementary slackness)
  auto rowMustBeBasic = [&](HighsInt row, double lowerSlack,
                            double upperSlack) {
    return std::abs(solution.row_dual[row]) <= dual_tol && lowerSlack > tol &&
           upperSlack > tol;
  };

  // assign row as basic
  auto assignBasicRowStatus = [&](HighsInt row, HighsInt& basicAssigned) {
    if (basis.row_status[row] == HighsBasisStatus::kBasic ||
        std::abs(solution.row_dual[row]) > dual_tol)
      return false;
    basis.row_status[row] = HighsBasisStatus::kBasic;
    basicAssigned++;
    return true;
  };

  // assign row as non-basic
  auto assignNonBasicRowStatus = [&](HighsInt row, double lowerSlack,
                                     double upperSlack) {
    if (solution.row_dual[row] > dual_tol)
      basis.row_status[row] = HighsBasisStatus::kLower;
    else if (solution.row_dual[row] < -dual_tol)
      basis.row_status[row] = HighsBasisStatus::kUpper;
    else
      basis.row_status[row] = upperSlack < lowerSlack
                                  ? HighsBasisStatus::kUpper
                                  : HighsBasisStatus::kLower;
  };

  // assign row status
  auto assignRowStatus = [&](HighsInt row, double lowerSlack, double upperSlack,
                             HighsInt& basicAssigned,
                             bool forceNonBasic = false) {
    if (forceNonBasic || !assignBasicRowStatus(row, basicAssigned))
      assignNonBasicRowStatus(row, lowerSlack, upperSlack);
  };

  // collect a single candidate for basic assignment
  auto collectCandidate =
      [&](HighsInt row, bool forcedNonBasic, double lowerSlack,
          double upperSlack, const std::vector<Nonzero>& entries,
          HighsInt parentIndex, bool isMinus,
          std::vector<std::tuple<HighsInt, HighsInt, bool>>& candidates) {
        if (basis.row_status[row] == HighsBasisStatus::kBasic) return;
        if (forcedNonBasic || std::abs(solution.row_dual[row]) > dual_tol) {
          assignNonBasicRowStatus(row, lowerSlack, upperSlack);
          return;
        }
        HighsInt nonBasicCount = 0;
        for (const auto& nz : entries)
          if (basis.col_status[nz.index] != HighsBasisStatus::kBasic)
            nonBasicCount++;
        candidates.emplace_back(nonBasicCount, parentIndex, isMinus);
      };

  for (HighsInt s = numSteps - 1; s >= 0; --s) {
    const auto& step = steps[s];
    HighsInt col = step.header.col;
    HighsInt numPlus = step.header.numPlus;
    HighsInt numMinus = step.header.numMinus;

    // compute slacks
    std::vector<double> plusLowerSlack;
    std::vector<double> plusUpperSlack;
    std::vector<double> minusLowerSlack;
    std::vector<double> minusUpperSlack;
    computeSlacks(col, step.plusRows, plusLowerSlack, plusUpperSlack);
    computeSlacks(col, step.minusRows, minusLowerSlack, minusUpperSlack);

    // non-basic propagation: if a generated row is non-basic (with nonzero
    // dual), both its parents are forced non-basic. mark them so the greedy
    // passes skip them. only force if the parent doesn't must-be-basic.
    std::vector<HighsBool> forcedNonBasicPlus(numPlus, false);
    std::vector<HighsBool> forcedNonBasicMinus(numMinus, false);
    for (const auto& nr : step.newRows) {
      // get indices of parent rows
      HighsInt p = nr.plusParentIdx;
      HighsInt m = nr.minusParentIdx;
      // skip basic rows (zero dual) and degenerate non-basic rows (with zero
      // dual)
      if (p < 0 || m < 0 || std::abs(solution.row_dual[nr.row]) <= dual_tol)
        continue;
      // mark rows that do not have to be basic
      if (!rowMustBeBasic(step.plusRows[p].header.row, plusLowerSlack[p],
                          plusUpperSlack[p]))
        forcedNonBasicPlus[p] = true;
      if (!rowMustBeBasic(step.minusRows[m].header.row, minusLowerSlack[m],
                          minusUpperSlack[m]))
        forcedNonBasicMinus[m] = true;
    }

    // mark ranged rows (appearing in both plus and minus sets)
    std::vector<HighsBool> isMinusRowRanged(numMinus, false);
    HighsInt numRanged = 0;
    for (HighsInt m = 0; m < numMinus; ++m)
      for (HighsInt p = 0; p < numPlus; ++p)
        if (step.minusRows[m].header.row == step.plusRows[p].header.row) {
          isMinusRowRanged[m] = true;
          numRanged++;
          break;
        }

    // count number of basic new rows
    HighsInt numNewRows = static_cast<HighsInt>(step.newRows.size());
    HighsInt numBasicNewRows = 0;
    for (const auto& nr : step.newRows)
      if (basis.row_status[nr.row] == HighsBasisStatus::kBasic)
        numBasicNewRows++;

    // how many basic variables are needed?
    HighsInt basicNeeded =
        (numPlus + numMinus - numRanged - numNewRows) + numBasicNewRows;
    HighsInt basicAssigned = 0;

    // determine col status
    bool colMustBeBasic =
        solution.col_value[col] > step.header.colLower + tol &&
        solution.col_value[col] < step.header.colUpper - tol;
    bool colCanBeBasic =
        colMustBeBasic || std::abs(solution.col_dual[col]) <= dual_tol;

    // pass 1: assign all must-be-basic (col and rows)
    if (colMustBeBasic) {
      basis.col_status[col] = HighsBasisStatus::kBasic;
      basicAssigned++;
    }
    for (HighsInt p = 0; p < numPlus; ++p) {
      if (forcedNonBasicPlus[p]) continue;
      if (rowMustBeBasic(step.plusRows[p].header.row, plusLowerSlack[p],
                         plusUpperSlack[p]))
        assignBasicRowStatus(step.plusRows[p].header.row, basicAssigned);
    }
    for (HighsInt m = 0; m < numMinus; ++m) {
      if (isMinusRowRanged[m] || forcedNonBasicMinus[m]) continue;
      if (rowMustBeBasic(step.minusRows[m].header.row, minusLowerSlack[m],
                         minusUpperSlack[m]))
        assignBasicRowStatus(step.minusRows[m].header.row, basicAssigned);
    }

    // pass 2: assign can-be-basic col (if not already assigned)
    if (!colMustBeBasic) {
      if (colCanBeBasic && basicAssigned < basicNeeded) {
        basis.col_status[col] = HighsBasisStatus::kBasic;
        basicAssigned++;
      } else if (solution.col_value[col] <= step.header.colLower + tol) {
        basis.col_status[col] = HighsBasisStatus::kLower;
      } else {
        basis.col_status[col] = HighsBasisStatus::kUpper;
      }
    }

    // pass 3: assign can-be-basic rows, sorted by non-basic support count
    // to reduce risk of rank deficiency in degenerate cases
    std::vector<std::tuple<HighsInt, HighsInt, bool>> candidates;
    for (HighsInt p = 0; p < numPlus; ++p)
      collectCandidate(step.plusRows[p].header.row, forcedNonBasicPlus[p],
                       plusLowerSlack[p], plusUpperSlack[p],
                       step.plusRows[p].entries, p, false, candidates);
    for (HighsInt m = 0; m < numMinus; ++m) {
      if (isMinusRowRanged[m]) continue;
      collectCandidate(step.minusRows[m].header.row, forcedNonBasicMinus[m],
                       minusLowerSlack[m], minusUpperSlack[m],
                       step.minusRows[m].entries, m, true, candidates);
    }
    // sort descending by non-basic support count
    pdqsort(candidates.begin(), candidates.end(),
            [](const std::tuple<HighsInt, HighsInt, bool>& a,
               const std::tuple<HighsInt, HighsInt, bool>& b) {
              return std::get<0>(a) > std::get<0>(b);
            });
    for (const auto& cand : candidates) {
      HighsInt parentIndex = std::get<1>(cand);
      if (std::get<2>(cand)) {
        assignRowStatus(step.minusRows[parentIndex].header.row,
                        minusLowerSlack[parentIndex],
                        minusUpperSlack[parentIndex], basicAssigned,
                        basicAssigned >= basicNeeded);
      } else {
        assignRowStatus(step.plusRows[parentIndex].header.row,
                        plusLowerSlack[parentIndex],
                        plusUpperSlack[parentIndex], basicAssigned,
                        basicAssigned >= basicNeeded);
      }
    }
  }
}

void HighsPostsolveStack::ZeroObjSingletonContinuousCol::undo(
    const HighsOptions& options, const std::vector<Nonzero>& rowValues,
    HighsSolution& solution, HighsBasis& basis) {
  assert(origRowLower != -kHighsInf && origRowUpper != kHighsInf);

  const double primal_tol = options.primal_feasibility_tolerance;
  const double dual_tol = options.dual_feasibility_tolerance;
  // Get activity of row without removed singleton
  HighsCDouble act = 0;
  solution.col_value[col] = 0.0;
  for (const auto& rowVal : rowValues)
    act += rowVal.value * solution.col_value[rowVal.index];

  // Find a suitable bound within the domain
  double col_value = kHighsInf;
  // Detemine whether the row was at its lower or upper bound by
  // virtue of basis status or value
  bool row_at_lower =
      (basis.valid && basis.row_status[row] == HighsBasisStatus::kLower) ||
      static_cast<double>(act - primal_tol) <= new_row_lb;
  bool row_at_upper =
      (basis.valid && basis.row_status[row] == HighsBasisStatus::kUpper) ||
      static_cast<double>(act + primal_tol) >= new_row_ub;
  bool col_at_lower = false;
  bool col_at_upper = false;
  if (row_at_lower) {
    if (coef > 0) {
      col_at_upper = true;
    } else {
      col_at_lower = true;
    }
  } else if (row_at_upper) {
    if (coef > 0) {
      col_at_lower = true;
    } else {
      col_at_upper = true;
    }
  } else {
    // Determine whether the column can be at its lower or upper bound
    // with the row between its bound basic
    double act_with_col_at_lower = static_cast<double>(act + coef * lb);
    double act_with_col_at_upper = static_cast<double>(act + coef * ub);
    if (act_with_col_at_lower >= origRowLower - primal_tol &&
        act_with_col_at_lower <= origRowUpper + primal_tol) {
      col_at_lower = true;
    } else if (act_with_col_at_upper >= origRowLower - primal_tol &&
               act_with_col_at_upper <= origRowUpper + primal_tol) {
      col_at_upper = true;
    } else {
      // Otherwise, there will be a column value for which the row is at a bound
      double col_value_for_lower =
          static_cast<double>((origRowLower - act) / coef);
      double col_value_for_upper =
          static_cast<double>((origRowUpper - act) / coef);
      if (col_value_for_lower >= lb - primal_tol &&
          col_value_for_lower <= ub + primal_tol) {
        // Can set column basic at this value
        col_value = col_value_for_lower;
        row_at_lower = true;
      } else if (col_value_for_upper >= lb - primal_tol &&
                 col_value_for_upper <= ub + primal_tol) {
        col_value = col_value_for_upper;
        row_at_upper = true;
      }
    }
  }
  if (col_at_lower || col_at_upper) {
    col_value = col_at_lower ? lb : ub;
    if (solution.dual_valid)
      solution.col_dual[col] = -coef * solution.row_dual[row];
    if (basis.valid)
      basis.col_status[col] =
          col_at_lower ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
  } else {
    // Column takes value between its bounds and inherits any basic status from
    // the row
    if (solution.dual_valid) {
      solution.col_dual[col] = 0;
      solution.row_dual[row] = 0;
    }
    if (basis.valid) {
      basis.row_status[row] =
          row_at_lower ? HighsBasisStatus::kLower : HighsBasisStatus::kUpper;
      basis.col_status[col] = HighsBasisStatus::kBasic;
    }
  }
  solution.col_value[col] = col_value;
}

}  // namespace presolve
