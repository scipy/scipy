/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::dualFixing(HighsPostsolveStack& postsolve_stack,
                                        HighsInt col) {
  // fix variables or tighten bounds using dual arguments
  // see section 4.4 "Dual fixing, substitution and bound strengthening",
  // Achterberg et al., Presolve Reductions in Mixed Integer Programming,
  // INFORMS Journal on Computing 32(2):473-506.
  assert(!colDeleted[col]);

  // return if variable is already fixed
  if (model->col_lower_[col] == model->col_upper_[col]) return Result::kOk;

  // struct for single equation handling
  struct equationNonZero {
    HighsInt col;
    double val;
    HighsInt mark;
  };

  // lambda for variable substitution
  auto substituteCol = [&](HighsInt col, HighsInt row, HighsInt direction,
                           double colBound, double otherColBound) {
    assert(model->col_lower_[col] != -kHighsInf &&
           model->col_upper_[col] != kHighsInf);
    constexpr HighsInt maxSubstRowSize = 1000;
    if (rowsizeInteger[row] == 0 || rowsize[row] > maxSubstRowSize)
      return Result::kOk;
    // check lhs and rhs for finiteness
    const bool lhsFinite = model->row_lower_[row] != -kHighsInf;
    const bool rhsFinite = model->row_upper_[row] != kHighsInf;
    if (lhsFinite && rhsFinite) {
      // Row cannot be redundant if the column to be substituted can
      // account for all row activity
      HighsInt pos = findNonzero(row, col);
      assert(pos != -1);
      const double colCoef = Avalue[pos];
      const HighsCDouble activityRange =
          std::abs(colCoef) *
          (static_cast<HighsCDouble>(model->col_upper_[col]) -
           model->col_lower_[col]);
      if (activityRange > static_cast<HighsCDouble>(model->row_upper_[row]) -
                              model->row_lower_[row] + 2 * primal_feastol) {
        return Result::kOk;
      }
    }

    // use storeRow and getStoredRow since getRowVector's rowroot[row] would be
    // overwritten by subsequent findNonZero calls, which would produce
    // undefined behavior
    storeRow(row);
    for (const auto& rowNz : getStoredRow()) {
      // skip column index that was passed to this lambda
      if (rowNz.index() == col) continue;

      // only consider non-fixed binary variables
      if (model->integrality_[rowNz.index()] != HighsVarType::kInteger ||
          model->col_lower_[rowNz.index()] != 0.0 ||
          model->col_upper_[rowNz.index()] != 1.0)
        continue;

      // skip binary variable if setting it to its lower bound does not make the
      // row redundant
      if ((rhsFinite && impliedRowBounds.getResidualSumUpperOrig(
                            row, rowNz.index(), rowNz.value()) >
                            model->row_upper_[row] + primal_feastol) ||
          (lhsFinite && impliedRowBounds.getResidualSumLowerOrig(
                            row, rowNz.index(), rowNz.value()) <
                            model->row_lower_[row] - primal_feastol))
        continue;

      // now compute the implied lower bound (direction = 1) or implied upper
      // bound (direction = -1) provided that the binary variable is set to its
      // upper bound.
      double bestBound =
          direction > 0
              ? computeImpliedLowerBound(col, rowNz.index(),
                                         model->col_upper_[rowNz.index()])
              : -computeImpliedUpperBound(col, rowNz.index(),
                                          model->col_upper_[rowNz.index()]);

      // check if lower / upper bound is implied
      if (bestBound >= direction * colBound - primal_feastol) {
        // substitute variable
        double offset = otherColBound;
        double scale = colBound - otherColBound;
        postsolve_stack.doubletonEquation(
            -1, col, rowNz.index(), 1.0, -scale, offset, model->col_lower_[col],
            model->col_upper_[col], 0.0, false, false,
            HighsPostsolveStack::RowType::kEq, HighsEmptySlice());
        HPRESOLVE_CHECKED_CALL(substitute(col, rowNz.index(), offset, scale));
        HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
        break;
      }
    }
    return Result::kOk;
  };

  // lambda that checks column for single equation handling
  auto checkColumn = [&](HighsInt col, HighsInt colDirection, HighsInt row) {
    if (colDirection * model->col_cost_[col] < 0) return false;
    for (const auto& colNz : getColumnVector(col)) {
      // skip equation we are inspecting; coefficient is positive (see sign
      // adjustment performed by caller)
      if (colNz.index() == row) continue;
      // skip redundant rows
      if (isRedundant(colNz.index())) continue;
      // another ranged row blocks variable in both directions -> cannot be a
      // candidate for fixing
      if (isRanged(colNz.index())) return false;
      // compute row direction
      HighsInt rowDirection =
          model->row_lower_[colNz.index()] == -kHighsInf &&
                  model->row_upper_[colNz.index()] != kHighsInf
              ? 1
              : -1;
      // coefficient must have the correct sign
      if (colDirection * rowDirection * colNz.value() < 0) return false;
    }
    return true;
  };

  // lambda that computes activities for single equation handling
  auto computeActivity = [&](HighsInt col, double val, HighsCDouble& activity,
                             bool& activityFinite, HighsInt direction) {
    // direction =  1 -> update upper bound on activity (supremum)
    // direction = -1 -> update lower bound on activity (infimum)
    double bound =
        direction * val > 0 ? model->col_upper_[col] : model->col_lower_[col];
    activityFinite = activityFinite && std::abs(bound) != kHighsInf;
    if (activityFinite) activity += static_cast<HighsCDouble>(val) * bound;
  };

  // lambda for fixing variables
  auto fixCols = [&](const std::vector<equationNonZero>& equationRowVector,
                     HighsInt direction) {
    for (const auto& rowNz : equationRowVector) {
      if (direction * rowNz.mark >= 0) continue;
      if (direction * rowNz.val > 0)
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, rowNz.col));
      else
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, rowNz.col));
    }
    return Result::kOk;
  };

  // lambda for handling single equations
  auto handleSingleEquation = [&](HighsInt row) {
    assert(isEquation(row));
    std::vector<equationNonZero> equationRowVector;
    equationRowVector.reserve(rowsize[row]);
    HighsInt numSPlus = 0;
    HighsInt numSMinus = 0;
    for (const auto& rowNz : getRowVector(row)) {
      // flip column direction if coefficient is negative
      // (equivalent to negating the column's coefficients)
      HighsInt colDirection = std::copysign(HighsInt{1}, rowNz.value());
      HighsInt mark = 0;
      if (checkColumn(rowNz.index(), colDirection, row)) {
        // store in S_+
        numSPlus++;
        mark = 1;
      } else if (checkColumn(rowNz.index(), -colDirection, row)) {
        // store in S_-
        numSMinus++;
        mark = -1;
      }
      // store non-zero
      equationRowVector.emplace_back(
          equationNonZero{rowNz.index(), rowNz.value(), mark});
    }
    // return if both sets are empty
    if (numSPlus == 0 && numSMinus == 0) return Result::kOk;
    // compute activities
    HighsCDouble activityTPlus = 0.0;
    HighsCDouble activityTMinus = 0.0;
    HighsCDouble activitySCPlus = 0.0;
    HighsCDouble activitySCMinus = 0.0;
    bool activityTPlusFinite = true;
    bool activityTMinusFinite = true;
    bool activitySCPlusFinite = true;
    bool activitySCMinusFinite = true;
    for (const auto& rowNz : equationRowVector) {
      if (rowNz.mark <= 0 ||
          model->integrality_[rowNz.col] != HighsVarType::kContinuous)
        // T_+
        computeActivity(rowNz.col, rowNz.val, activityTPlus,
                        activityTPlusFinite, HighsInt{1});
      else
        // S^C_+
        computeActivity(rowNz.col, rowNz.val, activitySCPlus,
                        activitySCPlusFinite, HighsInt{-1});
      if (rowNz.mark >= 0 ||
          model->integrality_[rowNz.col] != HighsVarType::kContinuous)
        // T_-
        computeActivity(rowNz.col, rowNz.val, activityTMinus,
                        activityTMinusFinite, HighsInt{-1});
      else
        // S^C_-
        computeActivity(rowNz.col, rowNz.val, activitySCMinus,
                        activitySCMinusFinite, HighsInt{1});
      // break if activities are not finite
      if ((numSMinus == 0 || !activityTPlusFinite || !activitySCPlusFinite) &&
          (numSPlus == 0 || !activityTMinusFinite || !activitySCMinusFinite))
        break;
    }
    // fix variables
    if (numSMinus > 0 && activityTPlusFinite && activitySCPlusFinite &&
        activityTPlus + activitySCPlus <=
            model->row_lower_[row] + primal_feastol)
      // fix all variables in S_-
      HPRESOLVE_CHECKED_CALL(fixCols(equationRowVector, HighsInt{1}));
    else if (numSPlus > 0 && activityTMinusFinite && activitySCMinusFinite &&
             activityTMinus + activitySCMinus >=
                 model->row_lower_[row] - primal_feastol)
      // fix all variables in S_+
      HPRESOLVE_CHECKED_CALL(fixCols(equationRowVector, HighsInt{-1}));
    return Result::kOk;
  };

  // lambda for computing tighter bounds
  auto hasTighterBound = [&](HighsInt col, HighsInt direction,
                             double currentBound, double& newBound) {
    // return if objective coefficient has wrong sign
    if (direction * model->col_cost_[col] < 0) return false;

    // do not accept huge bounds
    double hugeBound = primal_feastol / kHighsTiny;

    // compute worst-case bounds (direction = 1: lower bound,
    // direction = -1: upper bound)
    newBound = direction > 0 ? computeWorstCaseLowerBound(col)
                             : -computeWorstCaseUpperBound(col);

    // return if bound is too large
    if (newBound != -kHighsInf &&
        (newBound >= direction * currentBound - primal_feastol ||
         std::abs(newBound) > hugeBound))
      return false;

    // flip sign
    newBound *= direction;
    return true;
  };

  // lambda callback for lock computation
  HighsInt numDownLocks = 0;
  HighsInt numUpLocks = 0;
  HighsInt downLockRow = -1;
  HighsInt upLockRow = -1;
  auto lockCallback = [&](HighsInt row, bool hasDownLock, bool hasUpLock) {
    // count locks and remember row index
    if (hasUpLock) {
      numUpLocks++;
      upLockRow = row;
    }
    if (hasDownLock) {
      numDownLocks++;
      downLockRow = row;
    }
    // stop early if there are locks in both directions, since the variable
    // cannot be fixed in this case.
    return numDownLocks > 1 && numUpLocks > 1;
  };
  // compute locks
  computeLocks(col, true, lockCallback);

  // If there are no up (down) locks, the variable must/can be fixed
  // at its upper (lower) bound. Note that if there are no up (down)
  // locks, then cost <= 0 (>= 0).
  //
  // If |cost| > dual_feasibility_tolerance then fixing is forced, and
  // may identify unboundedness
  //
  // If |cost| <= dual_feasibility_tolerance, then the variable is
  // dual feasible at any value, so only fix if the corresponding
  // bound is finite
  if (numDownLocks == 0 || numUpLocks == 0) {
    if (numDownLocks == 0) {
      if (model->col_cost_[col] > options->dual_feasibility_tolerance ||
          model->col_lower_[col] > -kHighsInf)
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
    } else {
      assert(numUpLocks == 0);
      if (model->col_cost_[col] < -options->dual_feasibility_tolerance ||
          model->col_upper_[col] < kHighsInf)
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
    }
  } else {
    bool hasSingleDownLock = numDownLocks == 1 && downLockRow != -1;
    bool hasSingleUpLock = numUpLocks == 1 && upLockRow != -1;
    if (hasSingleDownLock || hasSingleUpLock) {
      HighsInt equationRow =
          hasSingleDownLock && isEquation(downLockRow)
              ? downLockRow
              : (hasSingleUpLock && isEquation(upLockRow) ? upLockRow : -1);
      if (equationRow != -1 && !singleEquationChecked[equationRow]) {
        // see section 6.1 "Extension of dual fixing for single equations",
        // Achterberg et al., Presolve Reductions in Mixed Integer
        // Programming, INFORMS Journal on Computing 32(2):473-506.
        HPRESOLVE_CHECKED_CALL(handleSingleEquation(equationRow));
        assert(!rowDeleted[equationRow]);
        singleEquationChecked[equationRow] = true;
        ++numSingleEquationChecked;
        if (colDeleted[col]) return Result::kOk;
      } else if (mipsolver != nullptr && model->col_lower_[col] != -kHighsInf &&
                 model->col_upper_[col] != kHighsInf) {
        // try substitution
        if (hasSingleDownLock) {
          HPRESOLVE_CHECKED_CALL(substituteCol(col, downLockRow, HighsInt{1},
                                               model->col_upper_[col],
                                               model->col_lower_[col]));
          if (colDeleted[col]) return Result::kOk;
        }
        if (hasSingleUpLock) {
          HPRESOLVE_CHECKED_CALL(substituteCol(col, upLockRow, HighsInt{-1},
                                               model->col_lower_[col],
                                               model->col_upper_[col]));
          if (colDeleted[col]) return Result::kOk;
        }
      }
    }
    // try to strengthen bounds
    double newBound = 0.0;
    if (hasTighterBound(col, HighsInt{1}, model->col_upper_[col], newBound)) {
      // do not make bounds inconsistent
      newBound = std::max(newBound, model->col_lower_[col]);
      // update upper bound
      if (newBound < model->col_upper_[col] - primal_feastol) {
        if (newBound == model->col_lower_[col])
          HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
        else if (model->integrality_[col] != HighsVarType::kContinuous)
          HPRESOLVE_CHECKED_CALL(changeColUpper(col, newBound));
      }
    } else if (hasTighterBound(col, HighsInt{-1}, model->col_lower_[col],
                               newBound)) {
      // do not make bounds inconsistent
      newBound = std::min(newBound, model->col_upper_[col]);
      if (newBound > model->col_lower_[col] + primal_feastol) {
        // update lower bound
        if (newBound == model->col_upper_[col])
          HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
        else if (model->integrality_[col] != HighsVarType::kContinuous)
          HPRESOLVE_CHECKED_CALL(changeColLower(col, newBound));
      }
    }
  }
  return Result::kOk;
}

}  // namespace presolve
