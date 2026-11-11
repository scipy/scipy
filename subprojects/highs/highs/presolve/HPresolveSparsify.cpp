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

HPresolve::Result HPresolve::sparsify(HighsPostsolveStack& postsolve_stack) {
  assert(!may_require_basis_postsolve_);
  assert(allow_rule_[kPresolveRuleSparsify]);
  std::vector<HighsPostsolveStack::Nonzero> sparsifyRows;
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleSparsify);
  HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
  HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
  std::vector<HighsInt> tmpEquations;
  tmpEquations.reserve(equations.size());

  const double minNonzeroVal = std::sqrt(primal_feastol);

  for (const auto& eq : equations) tmpEquations.emplace_back(eq.second);
  for (HighsInt eqrow : tmpEquations) {
    if (rowDeleted[eqrow]) continue;

    assert(!rowDeleted[eqrow]);
    assert(isEquation(eqrow));

    storeRow(eqrow);

    HighsInt sparsestCol = -1;
    HighsInt secondSparsestCol = -1;
    HighsInt sparsestColLen = kHighsIInf;
    HighsInt secondSparsestColLen = kHighsIInf;
    for (HighsInt nzPos : rowpositions) {
      HighsInt col = Acol[nzPos];
      if (colsize[col] < sparsestColLen) {
        secondSparsestCol = sparsestCol;
        secondSparsestColLen = sparsestColLen;
        sparsestCol = col;
        sparsestColLen = colsize[col];
      } else if (colsize[col] < secondSparsestColLen) {
        secondSparsestCol = col;
        secondSparsestColLen = colsize[col];
      }
    }

    assert(sparsestCol != -1 && secondSparsestCol != -1);
    assert(colsize[sparsestCol] <= colsize[secondSparsestCol]);

    std::map<double, HighsInt> possibleScales;
    sparsifyRows.clear();

    for (const HighsSliceNonzero& colNz : getColumnVector(sparsestCol)) {
      HighsInt candRow = colNz.index();
      if (candRow == eqrow) continue;

      possibleScales.clear();

      HighsInt misses = 0;
      // allow no fillin if a completely continuous row is used to cancel a
      // row that has integers as there are instances where this leads to a
      // huge deterioration of cut performance
      HighsInt maxMisses = 1;
      if (rowsizeInteger[eqrow] == 0 && rowsizeInteger[candRow] != 0)
        --maxMisses;
      for (const HighsSliceNonzero& nonzero : getStoredRow()) {
        double candRowVal;
        if (nonzero.index() == sparsestCol) {
          candRowVal = colNz.value();
        } else {
          HighsInt nzPos = findNonzero(candRow, nonzero.index());
          if (nzPos == -1) {
            if (model->integrality_[nonzero.index()] ==
                    HighsVarType::kInteger &&
                model->col_upper_[nonzero.index()] -
                        model->col_lower_[nonzero.index()] >
                    1.5) {
              // do not allow fillin of general integers
              misses = 2;
              break;
            }
            ++misses;
            if (misses > maxMisses) break;
            continue;
          }
          candRowVal = Avalue[nzPos];
        }

        double scale = -candRowVal / nonzero.value();
        if (std::abs(scale) > 1e3) continue;

        double scaleTolerance = minNonzeroVal / std::abs(nonzero.value());
        auto it = possibleScales.lower_bound(scale - scaleTolerance);
        if (it != possibleScales.end() &&
            std::abs(it->first - scale) <= scaleTolerance) {
          // there already is a scale that is very close and could produce
          // a matrix value for this nonzero that is below the allowed
          // threshold. Therefore we check if the matrix value is small enough
          // for this nonzero to be deleted, in which case the number of
          // deleted nonzeros for the other scale is increased. If it is not
          // small enough we do not use this scale or the other one because
          // such small matrix values may lead to numerical troubles.

          // scale is already marked to be numerically bad
          if (it->second == -1) continue;

          if (std::abs(it->first * nonzero.value() + candRowVal) <=
              options->small_matrix_value)
            it->second += 1;
          else
            it->second = -1;
        } else
          possibleScales.emplace(scale, 1);
      }

      if (misses > maxMisses || possibleScales.empty()) continue;

      HighsInt numCancel = 0;
      double scale = 0.0;

      for (const auto& s : possibleScales) {
        if (s.second <= misses) continue;

        if (s.second > numCancel ||
            (s.second == numCancel && std::abs(s.first) < std::abs(scale))) {
          scale = s.first;
          numCancel = s.second;
        }
      }

      assert(scale != 0.0 || numCancel == 0);

      // cancels at least one nonzero if the scale cancels more than there is
      // fillin
      if (numCancel > misses) sparsifyRows.emplace_back(candRow, scale);
    }

    if (model->integrality_[sparsestCol] != HighsVarType::kInteger ||
        (model->col_upper_[sparsestCol] - model->col_lower_[sparsestCol]) <
            1.5) {
      // now check for rows which do not contain the sparsest column but all
      // other columns by scanning the second sparsest column
      for (const HighsSliceNonzero& colNz :
           getColumnVector(secondSparsestCol)) {
        HighsInt candRow = colNz.index();
        if (candRow == eqrow) continue;

        if (rowsizeInteger[eqrow] == 0 && rowsizeInteger[candRow] != 0)
          continue;

        HighsInt sparsestColPos = findNonzero(candRow, sparsestCol);

        // if the row has a nonzero for the sparsest column we have already
        // checked it
        if (sparsestColPos != -1) continue;

        possibleScales.clear();
        bool skip = false;
        for (const HighsSliceNonzero& nonzero : getStoredRow()) {
          double candRowVal;
          if (nonzero.index() == secondSparsestCol) {
            candRowVal = colNz.value();
          } else {
            HighsInt nzPos = findNonzero(candRow, nonzero.index());
            // we already have a miss for the sparsest column, so with another
            // one we want to skip the row
            skip = nzPos == -1;
            if (skip) break;

            candRowVal = Avalue[nzPos];
          }

          double scale = -candRowVal / nonzero.value();
          if (std::abs(scale) > 1e3) continue;

          double scaleTolerance = minNonzeroVal / std::abs(nonzero.value());
          auto it = possibleScales.lower_bound(scale - scaleTolerance);
          if (it != possibleScales.end() &&
              std::abs(it->first - scale) <= scaleTolerance) {
            // there already is a scale that is very close and could produce
            // a matrix value for this nonzero that is below the allowed
            // threshold. Therefore we check if the matrix value is small
            // enough for this nonzero to be deleted, in which case the number
            // of deleted nonzeros for the other scale is increased. If it is
            // not small enough we do not use this scale or the other one
            // because such small matrix values may lead to numerical
            // troubles.

            // scale is already marked to be numerically bad
            if (it->second == -1) continue;

            if (std::abs(it->first * nonzero.value() + candRowVal) <=
                options->small_matrix_value) {
              it->second += 1;
            } else {
              // mark scale to be numerically bad
              it->second = -1;
              continue;
            }
          } else
            possibleScales.emplace(scale, 1);
        }

        if (skip || possibleScales.empty()) continue;

        HighsInt numCancel = 0;
        double scale = 0.0;

        for (const auto& s : possibleScales) {
          if (s.second <= 1) continue;
          if (s.second > numCancel ||
              (s.second == numCancel && std::abs(s.first) < std::abs(scale))) {
            scale = s.first;
            numCancel = s.second;
          }
        }

        assert(scale != 0.0 || numCancel == 0);

        // cancels at least one nonzero if the scale cancels more than there
        // is fillin
        if (numCancel > 1) sparsifyRows.emplace_back(candRow, scale);
      }
    }

    if (sparsifyRows.empty()) continue;

    postsolve_stack.equalityRowAdditions(eqrow, getStoredRow(), sparsifyRows);
    double rhs = model->row_lower_[eqrow];
    for (const auto& sparsifyRow : sparsifyRows) {
      HighsInt row = sparsifyRow.index;
      double scale = sparsifyRow.value;

      addToRowLower(row, static_cast<HighsCDouble>(scale) * rhs);
      addToRowUpper(row, static_cast<HighsCDouble>(scale) * rhs);

      for (HighsInt pos : rowpositions)
        addToMatrix(row, Acol[pos], scale * Avalue[pos]);

      reinsertEquation(row);
    }

    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(removeDoubletonEquations(postsolve_stack));
  }

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleSparsify);

  return Result::kOk;
}

}  // namespace presolve
