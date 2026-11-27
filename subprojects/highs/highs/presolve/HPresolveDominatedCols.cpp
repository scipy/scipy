/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
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

  // the two shortest rows of a column that the columns it dominates must also
  // appear in
  struct DominationRows {
    HighsInt shortestRow = -1;
    HighsInt secondShortestRow = -1;
    double coefShortestRow = 0.0;
    double coefSecondShortestRow = 0.0;
    // whether the column has a positive normalised coefficient (multiplied by
    // the direction) in some row. only then the direction is searched.
    bool hasPositiveCoef = false;
  };

  // a column and a direction to be searched, together with the two shortest
  // required rows
  struct ColToCheck {
    HighsInt shortestRow;
    HighsInt secondShortestRow;
    HighsInt col;
    HighsInt direction;
  };
  std::vector<ColToCheck> colsToCheck;

  // non-zero signatures for comparing columns
  std::vector<std::pair<uint64_t, uint64_t>> signatures(model->num_col_);

  // positions of the candidates and of the row they are collected from
  std::vector<HighsInt> check_rowpositions;
  std::vector<HighsInt> collect_rowpositions;

  // the rows for which the current candidates were collected: the shortest
  // required row (-1 if the candidates have to be collected again, e.g. after
  // fixing variables, since this modifies the matrix) and the second shortest
  // required row they were filtered by (-1 if they were not filtered)
  HighsInt candidatesRow = -1;
  HighsInt candidatesFilterRow = -1;

  // count number of fixed columns and modified bounds
  HighsInt numFixedCols = 0;
  HighsInt numModifiedBnds = 0;

  // check for timeout according to this frequency (number of columns)
  const size_t checkTimeFrequency = 100;

  // work limit: the work (candidates scanned, positions collected and
  // nonzeros looked up) may not exceed a base, which is proportional to the
  // number of nonzeros, plus a bonus for each reduction
  const int64_t numNonzeros =
      static_cast<int64_t>(Avalue.size() - freeslots.size());
  const int64_t workLimitBase =
      std::max(int64_t{10} * numNonzeros, int64_t{10000000});
  const int64_t workBonusPerReduction = 10000;
  int64_t work = 0;
  auto currentWorkLimit = [&]() {
    return workLimitBase +
           workBonusPerReduction * (numFixedCols + numModifiedBnds);
  };
  auto workLimitReached = [&]() { return work > currentWorkLimit(); };

  auto addSignature = [&](HighsInt row, HighsInt col, uint64_t rowLowerFinite,
                          uint64_t rowUpperFinite) {
    HighsInt rowHashedPos = (HighsHashHelpers::hash(row) >> 58);
    assert(rowHashedPos < 64);
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

  auto checkSignatures = [&](HighsInt scalj, HighsInt j, HighsInt scalk,
                             HighsInt k) {
    uint64_t sjMinus = signatures[j].first;
    uint64_t sjPlus = signatures[j].second;
    if (scalj == -1) std::swap(sjPlus, sjMinus);

    uint64_t skMinus = signatures[k].first;
    uint64_t skPlus = signatures[k].second;
    if (scalk == -1) std::swap(skPlus, skMinus);

    // the set of rows with a negative coefficient must be a superset of the
    // dominated columns set of rows with a negative coefficient
    if ((~sjMinus & skMinus) != 0) return false;

    // the set of rows with a positive coefficient must be a subset of the
    // dominated columns set of rows with a positive coefficient
    if ((sjPlus & ~skPlus) != 0) return false;

    return true;
  };

  auto checkDomination = [&](HighsInt scalj, HighsInt j, HighsInt scalk,
                             HighsInt k) {
    // rule out domination from integers to continuous variables
    if (model->integrality_[j] == HighsVarType::kInteger &&
        model->integrality_[k] == HighsVarType::kContinuous)
      return false;

    // the signatures are checked before calling this function
    assert(checkSignatures(scalj, j, scalk, k));

    // the dominating columns cost must be smaller or equal to the dominated
    // columns cost
    if (scalj * model->col_cost_[j] >
        scalk * model->col_cost_[k] + options->small_matrix_value)
      return false;

    // finally check the column vectors
    for (const HighsSliceNonzero& nonz : getColumnVector(j)) {
      work++;
      HighsInt row = nonz.index();
      HighsInt akPos = findNonzero(row, k);
      if (!checkDominationNonZero(row, scalj * nonz.value(),
                                  scalk * (akPos == -1 ? 0.0 : Avalue[akPos])))
        return false;
    }

    // check row only occurring in the column vector of k
    for (const HighsSliceNonzero& nonz : getColumnVector(k)) {
      work++;
      HighsInt row = nonz.index();
      HighsInt ajPos = findNonzero(row, j);
      // only rows in which aj does not occur are left to check
      if (ajPos != -1) continue;
      if (!checkDominationNonZero(row, 0.0, scalk * nonz.value())) return false;
    }

    return true;
  };

  // lambda for fixing variables
  auto fixCol = [&](HighsInt col, HighsInt direction) {
    numFixedCols++;
    candidatesRow = -1;
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
        return fixCol(col, HighsInt{1});
      } else if (model->integrality_[col] != HighsVarType::kContinuous) {
        numModifiedBnds++;
        HPRESOLVE_CHECKED_CALL(changeColLower(col, lowerBound));
      }
    }
    if (upperBound > -kHighsInf &&
        upperBound < model->col_upper_[col] - primal_feastol) {
      if (upperBound == model->col_lower_[col]) {
        return fixCol(col, HighsInt{-1});
      } else if (model->integrality_[col] != HighsVarType::kContinuous) {
        numModifiedBnds++;
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
  auto checkCols = [&](HighsInt col, HighsInt k, HighsInt direction,
                       HighsInt multiplier, bool boundImplied, bool hasCliques,
                       bool sameVarType) {
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

  // lambda for finding a domination relationship with the candidates (positions
  // in the shortest required row)
  auto checkRow = [&](const DominationRows& rows, HighsInt col,
                      HighsInt direction, bool hasCliques) {
    HighsInt row = rows.shortestRow;
    if (rowDeleted[row]) return Result::kOk;
    for (HighsInt pos : check_rowpositions) {
      work++;
      if (Avalue[pos] == 0.0 || Arow[pos] != row) continue;

      // get column index
      HighsInt k = Acol[pos];

      // skip column 'col' and deleted columns
      if (k == col || colDeleted[k]) continue;

      // check whether 'col' dominates k (multiplier 1) or its negation
      // (multiplier -1)
      for (HighsInt multiplier : {HighsInt{1}, HighsInt{-1}}) {
        // the signatures rule out domination for most candidates, so check
        // them first
        if (!checkSignatures(direction, col, multiplier * direction, k))
          continue;

        // the non-zero might have been modified by a reduction
        if (Avalue[pos] == 0.0 || Arow[pos] != row) break;

        // check the coefficients in this row in advance to avoid the
        // (potentially slow) element-wise comparison of the columns
        if (!checkDominationNonZero(row, direction * rows.coefShortestRow,
                                    multiplier * direction * Avalue[pos]))
          continue;

        // the signatures may pass due to hash collisions, so also check the
        // coefficients in the second shortest required row
        HighsInt secondRow = rows.secondShortestRow;
        if (secondRow != -1 && !rowDeleted[secondRow]) {
          work++;
          HighsInt secondPos = findNonzero(secondRow, k);
          if (secondPos == -1 ||
              !checkDominationNonZero(
                  secondRow, direction * rows.coefSecondShortestRow,
                  multiplier * direction * Avalue[secondPos]))
            continue;
        }

        // check if bound is implied (computed fresh due to earlier fixings)
        bool boundImplied =
            direction > 0 ? isUpperImplied(col) : isLowerImplied(col);

        // check if variables have the same type
        bool sameVarType = varsHaveSameType(col, k);

        // skip checks if nothing to do
        if (!boundImplied && !hasCliques && !sameVarType) continue;

        // try to fix variables or strengthen bounds
        HPRESOLVE_CHECKED_CALL(checkCols(col, k, direction, multiplier,
                                         boundImplied, hasCliques,
                                         sameVarType));
        if (colDeleted[col] || rowDeleted[row] || colDeleted[k]) break;
      }
      if (colDeleted[col] || rowDeleted[row]) break;
    }
    return Result::kOk;
  };

  // compare rows by their length and use the row index as a tie breaker
  auto rowIsShorter = [&](HighsInt row1, HighsInt row2) {
    return row2 == -1 || rowsize[row1] < rowsize[row2] ||
           (rowsize[row1] == rowsize[row2] && row1 < row2);
  };

  // update the required rows of a direction with a row in which the
  // normalised coefficient (multiplied by the direction) is 'val'
  auto updateDominationRows = [&](DominationRows& rows, HighsInt row,
                                  double coef, double val) {
    if (val > 0.0) {
      // the row is required and the direction is searched
      rows.hasPositiveCoef = true;
    } else if (!isRanged(row)) {
      // the row is not required
      return;
    }
    // equality and ranged rows are required in both directions, since the
    // coefficients of dominating and dominated columns must be equal there
    if (rowIsShorter(row, rows.shortestRow)) {
      rows.secondShortestRow = rows.shortestRow;
      rows.coefSecondShortestRow = rows.coefShortestRow;
      rows.shortestRow = row;
      rows.coefShortestRow = coef;
    } else if (rowIsShorter(row, rows.secondShortestRow)) {
      rows.secondShortestRow = row;
      rows.coefSecondShortestRow = coef;
    }
  };

  // determine the required rows of a column for the given direction
  auto getDominationRows = [&](HighsInt col, HighsInt direction) {
    DominationRows rows;
    for (const HighsSliceNonzero& nonz : getColumnVector(col)) {
      HighsInt row = nonz.index();
      HighsInt scale = model->row_upper_[row] != kHighsInf ? 1 : -1;
      updateDominationRows(rows, row, nonz.value(),
                           direction * scale * nonz.value());
    }
    return rows;
  };

  // check whether two columns share their two shortest required rows
  auto haveSameRows = [&](size_t i1, size_t i2) {
    return colsToCheck[i1].shortestRow == colsToCheck[i2].shortestRow &&
           colsToCheck[i1].secondShortestRow ==
               colsToCheck[i2].secondShortestRow;
  };

  // collect the positions in 'row' whose column also appears in 'filterRow'
  // (all positions if 'filterRow' is -1)
  auto collectCandidates = [&](HighsInt row, HighsInt filterRow) {
    getRowPositions(row, collect_rowpositions);
    work += static_cast<int64_t>(collect_rowpositions.size());
    check_rowpositions.clear();
    for (HighsInt pos : collect_rowpositions) {
      if (filterRow != -1) {
        work++;
        if (findNonzero(filterRow, Acol[pos]) == -1) continue;
      }
      check_rowpositions.push_back(pos);
    }
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

  colsToCheck.reserve(2 * static_cast<size_t>(model->num_col_));
  for (HighsInt j = 0; j < model->num_col_; ++j) {
    if (colDeleted[j]) continue;
    for (HighsInt direction : {HighsInt{-1}, HighsInt{1}}) {
      DominationRows rows = getDominationRows(j, direction);
      if (rows.hasPositiveCoef)
        colsToCheck.push_back(
            {rows.shortestRow, rows.secondShortestRow, j, direction});
    }
  }
  // sort the columns by their two shortest required rows, so that
  // columns sharing these rows are processed consecutively and their candidates
  // only have to be collected once
  pdqsort(colsToCheck.begin(), colsToCheck.end(),
          [](const ColToCheck& a, const ColToCheck& b) {
            if (a.shortestRow != b.shortestRow)
              return a.shortestRow < b.shortestRow;
            if (a.secondShortestRow != b.secondShortestRow)
              return a.secondShortestRow < b.secondShortestRow;
            if (a.col != b.col) return a.col < b.col;
            return a.direction < b.direction;
          });

  for (size_t groupStart = 0;
       groupStart < colsToCheck.size() && !workLimitReached();) {
    // find the columns sharing the two shortest required rows
    size_t groupEnd = groupStart + 1;
    while (groupEnd < colsToCheck.size() && haveSameRows(groupStart, groupEnd))
      ++groupEnd;
    // the columns of a group share their candidates, so these are collected
    // only once. they are filtered by the second shortest required row only if
    // the group has more than one column: checkRow checks this row for each
    // candidate anyway, but only after the cheaper signature test
    bool filterCandidates = groupEnd - groupStart > 1;

    for (size_t i = groupStart; i < groupEnd; ++i) {
      // check for timeout according to frequency
      if (i % checkTimeFrequency == 0) HPRESOLVE_CHECKED_CALL(checkTimeLimit());

      // stop if the work limit is reached
      if (workLimitReached()) break;

      HighsInt j = colsToCheck[i].col;
      HighsInt direction = colsToCheck[i].direction;

      // skip deleted columns
      if (colDeleted[j]) continue;

      // determine the required rows again, since reductions may have changed
      // them
      DominationRows dominationRows = getDominationRows(j, direction);

      // check if there are cliques
      bool colHasCliques =
          isBinary(j) && mipsolver->mipdata_->cliquetable.numCliques(
                             j, direction > 0 ? 1 : 0) > 0;

      // directions without a positive coefficient are not searched (for
      // columns that only appear in equality rows, the dominations are also
      // found from the other column)
      if (dominationRows.hasPositiveCoef) {
        // a column whose required rows have been changed by reductions no
        // longer belongs to the group, so its candidates are used by this
        // column only and are not filtered
        HighsInt filterRow =
            filterCandidates &&
                    dominationRows.shortestRow ==
                        colsToCheck[groupStart].shortestRow &&
                    dominationRows.secondShortestRow ==
                        colsToCheck[groupStart].secondShortestRow
                ? dominationRows.secondShortestRow
                : -1;
        // the candidates only have to be collected if they were collected for
        // different rows or have become invalid
        if (dominationRows.shortestRow != candidatesRow ||
            filterRow != candidatesFilterRow) {
          collectCandidates(dominationRows.shortestRow, filterRow);
          candidatesRow = dominationRows.shortestRow;
          candidatesFilterRow = filterRow;
        }
        HPRESOLVE_CHECKED_CALL(
            checkRow(dominationRows, j, direction, colHasCliques));
      }
    }
    groupStart = groupEnd;
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Dominated columns: fixed %" HIGHSINT_FORMAT
              " columns, strengthened %" HIGHSINT_FORMAT
              " bounds, used %" PRId64 " of %" PRId64 " work units%s\n",
              numFixedCols, numModifiedBnds, work, currentWorkLimit(),
              workLimitReached() ? " (limit reached)" : "");

  return Result::kOk;
}

}  // namespace presolve
