/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolverData.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::dominatedColumns(
    HighsPostsolveStack& postsolve_stack) {
  // See section 6.4 "Dominated columns", Achterberg et al., Presolve Reductions
  // in Mixed Integer Programming, INFORMS Journal on Computing 32(2):473-506.
  // See also Gamrath, G., Koch, T., Martin, A. et al., Progress in presolving
  // for mixed integer programming, Math. Prog. Comp. 7, 367–398 (2015).
  //
  // Given a domination relationship x_j ≻ x_k (Definition 1 in Gamrath
  // et al.), there exists an optimal solution where x_j = u_j or x_k = l_k
  // (Theorem 2). If this can be verified (e.g. via implied bounds or
  // cliques), the dominated variable is fixed. Otherwise, predictive bound
  // analysis (Theorem 3) derives tighter bounds on x_j and x_k by computing
  // implied bounds conditional on the other variable being at its extreme.
  // Worst-case bounds (MINU/MAXU) use the worst-case row activity (maximal
  // for ≤-constraints, minimal for ≥-constraints) to obtain bounds that are
  // feasible independent of other variables' values.

  // non-zero signatures for comparing columns
  std::vector<std::pair<uint32_t, uint32_t>> signatures(model->num_col_);

  // count overall number of domination checks and number of checks performed
  // for predictive bound analysis
  size_t numDomChecks = 0;
  size_t numDomChecksPredBndAnalysis = 0;

  auto addSignature = [&](HighsInt row, HighsInt col, uint32_t rowLowerFinite,
                          uint32_t rowUpperFinite) {
    HighsInt rowHashedPos = (HighsHashHelpers::hash(row) >> 59);
    assert(rowHashedPos < 32);
    signatures[col].first |= rowLowerFinite << rowHashedPos;
    signatures[col].second |= rowUpperFinite << rowHashedPos;
  };

  auto checkDominationNonZero = [&](HighsInt row, double aj, double ak) {
    if (isRanged(row)) {
      // the row is an equality or ranged row, therefore the coefficients must
      // be parallel, otherwise one of the inequalities given by the row rules
      // out domination
      if (std::abs(aj - ak) > options->small_matrix_value) return false;
      return true;
    }

    // normalize row to a <= constraint
    if (model->row_upper_[row] == kHighsInf) {
      aj = -aj;
      ak = -ak;
    }

    // the coefficient of the dominating column needs to be smaller than or
    // equal to the coefficient of the dominated column
    if (aj > ak + options->small_matrix_value) return false;
    return true;
  };

  auto checkDomination = [&](HighsInt scalj, HighsInt j, HighsInt scalk,
                             HighsInt k) {
    // increment counter
    numDomChecks++;

    // rule out domination from integers to continuous variables
    if (model->integrality_[j] == HighsVarType::kInteger &&
        model->integrality_[k] == HighsVarType::kContinuous)
      return false;

    // check the signatures
    uint32_t sjMinus = signatures[j].first;
    uint32_t sjPlus = signatures[j].second;
    if (scalj == -1) std::swap(sjPlus, sjMinus);

    uint32_t skMinus = signatures[k].first;
    uint32_t skPlus = signatures[k].second;
    if (scalk == -1) std::swap(skPlus, skMinus);

    // the set of rows with a negative coefficient must be a superset of the
    // dominated columns set of rows with a negative coefficient
    if ((~sjMinus & skMinus) != 0) return false;

    // the set of rows with a positive coefficient must be a subset of the
    // dominated columns set of rows with a positive coefficient
    if ((sjPlus & ~skPlus) != 0) return false;

    // the dominating columns cost must be smaller or equal to the dominated
    // columns cost
    if (scalj * model->col_cost_[j] >
        scalk * model->col_cost_[k] + options->small_matrix_value)
      return false;

    // finally check the column vectors
    for (const HighsSliceNonzero& nonz : getColumnVector(j)) {
      HighsInt row = nonz.index();
      HighsInt akPos = findNonzero(row, k);
      if (!checkDominationNonZero(row, scalj * nonz.value(),
                                  scalk * (akPos == -1 ? 0.0 : Avalue[akPos])))
        return false;
    }

    // check row only occurring in the column vector of k
    for (const HighsSliceNonzero& nonz : getColumnVector(k)) {
      HighsInt row = nonz.index();
      HighsInt ajPos = findNonzero(row, j);
      // only rows in which aj does not occur are left to check
      if (ajPos != -1) continue;
      if (!checkDominationNonZero(row, 0.0, scalk * nonz.value())) return false;
    }

    return true;
  };

  for (size_t i = 0; i != Avalue.size(); ++i) {
    if (Avalue[i] == 0) continue;

    HighsInt row = Arow[i];
    HighsInt col = Acol[i];
    bool rowLowerFinite = model->row_lower_[row] != -kHighsInf;
    bool rowUpperFinite = model->row_upper_[row] != kHighsInf;

    if (Avalue[i] > 0)
      addSignature(row, col, rowLowerFinite, rowUpperFinite);
    else
      addSignature(row, col, rowUpperFinite, rowLowerFinite);
  }

  // count number of fixed columns and modified bounds
  HighsInt numCols = 0;
  HighsInt numFixedCols = 0;
  HighsInt numFixedColsPredBndAnalysis = 0;
  HighsInt numModifiedBndsPredBndAnalysis = 0;

  // parameters for predictive bound analysis
  const size_t maxAverageNumDomChecksPredBndAnalysis = 10000;
  const double minAverageNumRedsPredBndAnalysis = 1e-2;

  // perform predictive bound analysis?
  bool allowPredBndAnalysis = true;
  std::vector<HighsInt> check_rowpositions;

  for (HighsInt j = 0; j < model->num_col_; ++j) {
    // skip deleted columns
    if (colDeleted[j]) continue;

    // increment counter for number of columns
    numCols++;

    // initialise
    HighsInt bestRowPlus = -1;
    HighsInt bestRowPlusLen = kHighsIInf;
    double ajBestRowPlus = 0.0;
    HighsInt bestRowMinus = -1;
    HighsInt bestRowMinusLen = kHighsIInf;
    double ajBestRowMinus = 0.0;

    for (const HighsSliceNonzero& nonz : getColumnVector(j)) {
      HighsInt row = nonz.index();
      HighsInt scale = model->row_upper_[row] != kHighsInf ? 1 : -1;

      double val = scale * nonz.value();
      if (val > 0.0 && rowsize[row] < bestRowPlusLen) {
        bestRowPlus = row;
        bestRowPlusLen = rowsize[row];
        ajBestRowPlus = nonz.value();
      }

      if (val < 0.0 && rowsize[row] < bestRowMinusLen) {
        bestRowMinus = row;
        bestRowMinusLen = rowsize[row];
        ajBestRowMinus = nonz.value();
      }
    }

    // lambda for fixing variables
    auto fixCol = [&](HighsInt col, HighsInt direction) {
      numFixedCols++;
      if (direction > 0) {
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
      } else {
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
      }
      // remove row singletons and doubleton equations
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
      HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
      return Result::kOk;
    };

    // Predictive bound analysis (Theorem 3 from Gamrath et al. 2015).
    // The paper defines x_j ≻ x_k (j dominates k). The code checks
    // (direction * x_col) ≻ (direction_k * x_otherCol):
    //   direction = +1, multiplier = +1: x_j =  x_col,      x_k =  x_otherCol
    //   direction = +1, multiplier = -1: x_j =  x_col,      x_k = -x_otherCol
    //   direction = -1, multiplier = +1: x_j =  x_otherCol, x_k =  x_col
    //   direction = -1, multiplier = -1: x_j = -x_col,      x_k =  x_otherCol
    // col is the column being tightened, otherCol is the conditioning
    // column (whose value otherColBound is substituted).
    // colIsAtUpper = true : col = x_j, apply (i)/(iii)/(v)
    // colIsAtUpper = false: col = x_k, apply (ii)/(iv)/(vi)
    // otherColCoeffPattern: +1 = same-sign, -1 = opposite-sign coefficient
    //                       filter (opposite signs from negated-column
    //                       domination)
    auto tightenBounds = [&](HighsInt col, double colBound, bool colIsAtUpper,
                             HighsInt otherCol, double otherColBound,
                             HighsInt otherColCoeffPattern) {
      // bound should be finite
      assert(std::abs(otherColBound) != kHighsInf);
      // return if variable is already fixed
      if (model->col_lower_[col] == model->col_upper_[col]) return Result::kOk;
      // initialise bounds
      double lowerBound = -kHighsInf;
      double upperBound = kHighsInf;
      if (colIsAtUpper) {
        // For negated-column domination the paper's x_k is negated,
        // so (i)/(iii)/(v) below correspond to (ii)/(iv)/(vi) in the
        // paper with negated bounds and cost.
        // (i) x_j <= MINL^k_j(otherColBound)
        upperBound = computeImpliedUpperBound(col, otherCol, otherColBound,
                                              otherColCoeffPattern);
        // (iii) x_j >= min{colBound, MAXL^k_j(otherColBound)}
        lowerBound = std::min(
            colBound, computeImpliedLowerBound(col, otherCol, otherColBound,
                                               otherColCoeffPattern));
        if (model->col_cost_[col] <= 0) {
          // (v) if c_j <= 0, then x_j >= min{colBound,
          //                                  MINU^k_j(otherColBound)}
          double worstCaseUpper = computeWorstCaseUpperBound(
              col, otherCol, otherColBound, otherColCoeffPattern);
          lowerBound = std::max(lowerBound, std::min(colBound, worstCaseUpper));
        }
      } else {
        // For negated-column domination the paper's x_j is negated,
        // so (ii)/(iv)/(vi) below correspond to (i)/(iii)/(v) in the
        // paper with negated bounds and cost.
        // (ii) x_k >= MAXL^j_k(otherColBound)
        lowerBound = computeImpliedLowerBound(col, otherCol, otherColBound,
                                              otherColCoeffPattern);
        // (iv) x_k <= max{colBound, MINL^j_k(otherColBound)}
        upperBound = std::max(
            colBound, computeImpliedUpperBound(col, otherCol, otherColBound,
                                               otherColCoeffPattern));
        if (model->col_cost_[col] >= 0) {
          // (vi) if c_k >= 0, then x_k <= max{colBound,
          //                                   MAXU^j_k(otherColBound)}
          double worstCaseLower = computeWorstCaseLowerBound(
              col, otherCol, otherColBound, otherColCoeffPattern);
          upperBound = std::min(upperBound, std::max(colBound, worstCaseLower));
        }
      }
      // update bounds
      if (lowerBound < kHighsInf &&
          lowerBound > model->col_lower_[col] + primal_feastol) {
        if (lowerBound == model->col_upper_[col]) {
          numFixedColsPredBndAnalysis++;
          return fixCol(col, HighsInt{1});
        } else if (model->integrality_[col] != HighsVarType::kContinuous) {
          numModifiedBndsPredBndAnalysis++;
          HPRESOLVE_CHECKED_CALL(changeColLower(col, lowerBound));
        }
      }
      if (upperBound > -kHighsInf &&
          upperBound < model->col_upper_[col] - primal_feastol) {
        if (upperBound == model->col_lower_[col]) {
          numFixedColsPredBndAnalysis++;
          return fixCol(col, HighsInt{-1});
        } else if (model->integrality_[col] != HighsVarType::kContinuous) {
          numModifiedBndsPredBndAnalysis++;
          HPRESOLVE_CHECKED_CALL(changeColUpper(col, upperBound));
        }
      }
      return Result::kOk;
    };

    // lambda for checking whether two variables have the same type
    auto varsHaveSameType = [&](HighsInt col1, HighsInt col2) {
      return ((model->integrality_[col1] != HighsVarType::kContinuous &&
               model->integrality_[col2] != HighsVarType::kContinuous) ||
              (model->integrality_[col1] != HighsVarType::kInteger &&
               model->integrality_[col2] != HighsVarType::kInteger));
    };

    // lambda for (1) checking whether one of the two columns is dominated by
    // the other one and (2) fixing variables or strengthening bounds
    auto checkCols = [&](HighsInt row, HighsInt col, HighsInt k,
                         HighsInt direction, HighsInt multiplier,
                         bool boundImplied, bool hasCliques, bool sameVarType) {
      // compute direction for dominated variable
      // (1: lower bound; -1: upper bound)
      HighsInt direction_k = multiplier * direction;
      // get bounds
      double dominatingBound =
          direction > 0 ? model->col_upper_[col] : model->col_lower_[col];
      double dominatedBound =
          direction_k > 0 ? model->col_lower_[k] : model->col_upper_[k];
      // check if bounds are finite
      bool isDominatingBoundFinite = direction * dominatingBound != kHighsInf;
      bool isDominatedBoundFinite = direction_k * dominatedBound != -kHighsInf;
      // check whether variable 'k' can potentially be fixed. check if there are
      // cliques in advance instead of directly searching for a common clique.
      bool tryToFix = isDominatedBoundFinite && (boundImplied || hasCliques);
      // check whether predictive bound analysis can be performed. both
      // variables need to have the same type.
      bool tryToStrengthenBounds =
          (isDominatingBoundFinite || isDominatedBoundFinite) && sameVarType;
      if (tryToFix || tryToStrengthenBounds) {
        // increment counter for number of domination checks due to predictive
        // bound analysis
        if (!tryToFix) numDomChecksPredBndAnalysis++;
        // check for domination
        if (checkDomination(direction, col, direction_k, k)) {
          if (tryToFix &&
              (boundImplied ||
               mipsolver->mipdata_->cliquetable.haveCommonClique(
                   HighsCliqueTable::CliqueVar(col, direction > 0 ? 1 : 0),
                   HighsCliqueTable::CliqueVar(k, direction_k > 0 ? 1 : 0)))) {
            // direction =  1, multiplier =  1:
            // case (i)   ub(x_j) =  inf,  x_j >  x_k: set x_k = lb(x_k)
            // direction =  1, multiplier = -1:
            // case (ii)  ub(x_j) =  inf,  x_j > -x_k: set x_k = ub(x_k)
            // direction = -1, multiplier =  1:
            // case (iii) lb(x_j) = -inf, -x_j > -x_k: set x_k = ub(x_k)
            // direction = -1, multiplier = -1:
            // case (iv)  lb(x_j) = -inf, -x_j >  x_k: set x_k = lb(x_k)
            HPRESOLVE_CHECKED_CALL(fixCol(k, -direction_k));
          } else if (tryToStrengthenBounds) {
            // tighten bounds via predictive bound analysis, see Theorem 3
            // from Gamrath et al.'s paper
            if (isDominatedBoundFinite)
              HPRESOLVE_CHECKED_CALL(
                  tightenBounds(col, dominatingBound, direction > 0, k,
                                dominatedBound, direction * direction_k));
            if (!colDeleted[col] && isDominatingBoundFinite)
              HPRESOLVE_CHECKED_CALL(
                  tightenBounds(k, dominatedBound, direction_k < 0, col,
                                dominatingBound, direction * direction_k));
          }
        }
      }
      return Result::kOk;
    };

    // lambda for finding a domination relationship in the given row
    auto checkRow = [&](HighsInt row, HighsInt col, HighsInt direction,
                        double bestVal, bool hasCliques) {
      if (rowDeleted[row]) return Result::kOk;
      getRowPositions(row, check_rowpositions);
      for (HighsInt pos : check_rowpositions) {
        if (Avalue[pos] == 0.0 || Arow[pos] != row) continue;

        // get column index
        HighsInt k = Acol[pos];

        // skip column 'col' and deleted columns
        if (k == col || colDeleted[k]) continue;

        // get coefficient
        double ak = Avalue[pos];

        // check if variables have the same type
        bool sameVarType = varsHaveSameType(col, k);

        // check if bound is implied (computed fresh due to earlier fixings)
        bool boundImplied =
            direction > 0 ? isUpperImplied(col) : isLowerImplied(col);

        // skip checks if nothing to do
        if (!boundImplied && !hasCliques && !sameVarType) continue;

        // try to fix variables or strengthen bounds
        // check already known non-zeros in respective columns in advance to
        // avoid (potentially slow) element-wise comparison if possible.
        if (checkDominationNonZero(row, direction * bestVal, direction * ak)) {
          HPRESOLVE_CHECKED_CALL(checkCols(row, col, k, direction, HighsInt{1},
                                           boundImplied, hasCliques,
                                           sameVarType));
          if (colDeleted[col] || rowDeleted[row]) break;
        }

        if (!colDeleted[k] && Avalue[pos] != 0.0 && Arow[pos] == row &&
            checkDominationNonZero(row, direction * bestVal, -direction * ak)) {
          HPRESOLVE_CHECKED_CALL(checkCols(row, col, k, direction, HighsInt{-1},
                                           boundImplied, hasCliques,
                                           sameVarType));
          if (colDeleted[col] || rowDeleted[row]) break;
        }
      }
      return Result::kOk;
    };

    // check if bounds are implied or there are cliques
    bool lowerImplied = isLowerImplied(j);
    bool upperImplied = isUpperImplied(j);
    bool hasNegCliques =
        isBinary(j) && mipsolver->mipdata_->cliquetable.numCliques(j, 0) > 0;
    bool hasPosCliques =
        isBinary(j) && mipsolver->mipdata_->cliquetable.numCliques(j, 1) > 0;

    // use row 'bestRowMinus'
    if (bestRowMinus != -1 &&
        (allowPredBndAnalysis || lowerImplied || hasNegCliques))
      HPRESOLVE_CHECKED_CALL(checkRow(bestRowMinus, j, HighsInt{-1},
                                      ajBestRowMinus, hasNegCliques));

    // use row 'bestRowPlus'
    if (!colDeleted[j] && bestRowPlus != -1 &&
        (allowPredBndAnalysis || upperImplied || hasPosCliques))
      HPRESOLVE_CHECKED_CALL(
          checkRow(bestRowPlus, j, HighsInt{1}, ajBestRowPlus, hasPosCliques));

    // do not use predictive bound analysis if it requires many domination
    // checks and only yields few fixings or improved bounds on average
    size_t averageNumDomChecksPredBndAnalysis =
        numDomChecksPredBndAnalysis / static_cast<size_t>(numCols);
    double averageNumRedsPredBndAnalysis =
        (numFixedColsPredBndAnalysis + numModifiedBndsPredBndAnalysis) /
        static_cast<double>(numCols);
    allowPredBndAnalysis =
        allowPredBndAnalysis &&
        (numDomChecksPredBndAnalysis <=
             30 * maxAverageNumDomChecksPredBndAnalysis ||
         (averageNumDomChecksPredBndAnalysis <=
              maxAverageNumDomChecksPredBndAnalysis &&
          averageNumRedsPredBndAnalysis >= minAverageNumRedsPredBndAnalysis));
  }

  if (numFixedCols > 0 || numModifiedBndsPredBndAnalysis > 0)
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "Fixed %d dominated columns and strengthened %d bounds\n",
                static_cast<int>(numFixedCols),
                static_cast<int>(numModifiedBndsPredBndAnalysis));

  return Result::kOk;
}

}  // namespace presolve
