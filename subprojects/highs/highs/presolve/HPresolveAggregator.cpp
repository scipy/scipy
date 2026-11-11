/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::aggregator(HighsPostsolveStack& postsolve_stack) {
  assert(this->allow_rule_[kPresolveRuleAggregator]);
  const bool logging_on = analysis_.logging_on_;
  if (logging_on) analysis_.startPresolveRuleLog(kPresolveRuleAggregator);
  substitutionOpportunities.erase(
      std::remove_if(substitutionOpportunities.begin(),
                     substitutionOpportunities.end(),
                     [&](const std::pair<HighsInt, HighsInt>& p) {
                       HighsInt row = p.first;
                       HighsInt col = p.second;
                       return rowDeleted[row] || colDeleted[col] ||
                              !isImpliedFree(col) || !isDualImpliedFree(row);
                     }),
      substitutionOpportunities.end());

  pdqsort(
      substitutionOpportunities.begin(), substitutionOpportunities.end(),
      [&](const std::pair<HighsInt, HighsInt>& nz1,
          const std::pair<HighsInt, HighsInt>& nz2) {
        HighsInt minLen1 = std::min(rowsize[nz1.first], colsize[nz1.second]);
        HighsInt minLen2 = std::min(rowsize[nz2.first], colsize[nz2.second]);
        if (minLen1 == 2 && minLen2 != 2) return true;
        if (minLen2 == 2 && minLen1 != 2) return false;

        int64_t sizeProd1 = int64_t(rowsize[nz1.first]) * colsize[nz1.second];
        int64_t sizeProd2 = int64_t(rowsize[nz2.first]) * colsize[nz2.second];
        if (sizeProd1 < sizeProd2) return true;
        if (sizeProd2 < sizeProd1) return false;
        if (minLen1 < minLen2) return true;
        if (minLen2 < minLen1) return false;

        return std::make_tuple(HighsHashHelpers::hash(std::make_pair(
                                   uint32_t(nz1.first), uint32_t(nz1.second))),
                               nz1.first, nz1.second) <
               std::make_tuple(HighsHashHelpers::hash(std::make_pair(
                                   uint32_t(nz2.first), uint32_t(nz2.second))),
                               nz2.first, nz2.second);
      });

  HighsInt nfail = 0;
  for (size_t i = 0; i < substitutionOpportunities.size(); ++i) {
    HighsInt row = substitutionOpportunities[i].first;
    HighsInt col = substitutionOpportunities[i].second;

    if (rowDeleted[row] || colDeleted[col] || !isImpliedFree(col) ||
        !isDualImpliedFree(row)) {
      substitutionOpportunities[i].first = -1;
      continue;
    }

    HighsInt nzPos = findNonzero(row, col);
    if (nzPos == -1) {
      substitutionOpportunities[i].first = -1;
      continue;
    }
    if (model->integrality_[col] == HighsVarType::kInteger) {
      StatusResult impliedIntegral = isImpliedIntegral(col);
      HPRESOLVE_CHECKED_CALL(static_cast<Result>(impliedIntegral));
      if (!impliedIntegral) continue;
    }

    // in the case where the row has length two or the column has length two
    // we always do the substitution since the fillin can never be problematic
    if (rowsize[row] == 2 || colsize[col] == 2) {
      storeRow(row);
      substituteFreeCol(postsolve_stack, row, col, true);
      substitutionOpportunities[i].first = -1;
      HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
      HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
      continue;
    }

    double maxVal = rowsize[row] < colsize[col] ? getMaxAbsRowVal(row)
                                                : getMaxAbsColVal(col);
    if (std::fabs(Avalue[nzPos]) < maxVal * options->presolve_pivot_threshold) {
      maxVal = rowsize[row] < colsize[col] ? getMaxAbsColVal(col)
                                           : getMaxAbsRowVal(row);
      if (std::fabs(Avalue[nzPos]) <
          maxVal * options->presolve_pivot_threshold) {
        substitutionOpportunities[i].first = -1;
        continue;
      }
    }

    storeRow(row);
    HighsInt fillin = -(rowsize[row] + colsize[col] - 1);
    for (const auto& nz : getColumnVector(col)) {
      if (nz.index() == row) continue;
      fillin += countFillin(nz.index());

      if (fillin > options->presolve_substitution_maxfillin) break;
    }

    if (fillin > options->presolve_substitution_maxfillin) {
      ++nfail;
      // if the fill in is too much for multiple tries, then we stop
      // as this indicates that the rows/columns are becoming too dense
      // for substitutions
      if (nfail == 3) break;
      continue;
    }

    nfail = 0;
    substituteFreeCol(postsolve_stack, row, col, true);
    substitutionOpportunities[i].first = -1;
    HPRESOLVE_CHECKED_CALL(removeRowSingletons(postsolve_stack));
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  substitutionOpportunities.erase(
      std::remove_if(
          substitutionOpportunities.begin(), substitutionOpportunities.end(),
          [](const std::pair<HighsInt, HighsInt>& p) { return p.first == -1; }),
      substitutionOpportunities.end());

  analysis_.logging_on_ = logging_on;
  if (logging_on) analysis_.stopPresolveRuleLog(kPresolveRuleAggregator);
  return Result::kOk;
}

HighsInt HPresolve::countFillin(HighsInt row) {
  HighsInt fillin = 0;
  for (HighsInt rowiter : rowpositions) {
    if (findNonzero(row, Acol[rowiter]) == -1) fillin += 1;
  }

  return fillin;
}

void HPresolve::substituteFreeCol(HighsPostsolveStack& postsolve_stack,
                                  HighsInt row, HighsInt col,
                                  bool relaxRowDualBounds) {
  assert(!rowDeleted[row]);
  assert(!colDeleted[col]);
  assert(isDualImpliedFree(row));
  assert(isImpliedFree(col));

  HighsPostsolveStack::RowType rowType;
  double rhs;
  dualImpliedFreeGetRhsAndRowType(row, rhs, rowType, relaxRowDualBounds);

  postsolve_stack.freeColSubstitution(row, col, rhs, model->col_cost_[col],
                                      rowType, getStoredRow(),
                                      getColumnVector(col));

  // todo, check integrality of coefficients and allow this
  substitute(row, col, rhs);
}

void HPresolve::substitute(HighsInt row, HighsInt col, double rhs) {
  assert(!rowDeleted[row]);
  assert(!colDeleted[col]);
  HighsInt pos = findNonzero(row, col);
  assert(pos != -1);

  assert(Arow[pos] == row);
  assert(Acol[pos] == col);
  double substrowscale = -1.0 / Avalue[pos];
  assert(isImpliedFree(col));

  markRowDeleted(row);
  markColDeleted(col);

  // substitute the column in each row where it occurs
  for (HighsInt coliter = colhead[col]; coliter != -1;) {
    HighsInt colrow = Arow[coliter];
    double colval = Avalue[coliter];

    // walk to the next position before doing any modifications, because
    // the current position will be deleted in the loop below
    assert(Acol[coliter] == col);
    HighsInt colpos = coliter;
    coliter = Anext[coliter];

    // skip the row that is used for substitution
    if (row == colrow) continue;

    assert(findNonzero(colrow, col) != -1);

    // cancels out and bounds of dual row for this column do not need to be
    // updated
    unlink(colpos);

    // printf("\nbefore substitution: ");
    // debugPrintRow(colrow);

    // determine the scale for the substitution row for addition to this row
    double scale = colval * substrowscale;

    // adjust the sides
    addToRowLower(colrow, static_cast<HighsCDouble>(scale) * rhs);
    addToRowUpper(colrow, static_cast<HighsCDouble>(scale) * rhs);

    for (HighsInt rowiter : rowpositions) {
      assert(Arow[rowiter] == row);

      if (Acol[rowiter] != col)
        addToMatrix(colrow, Acol[rowiter], scale * Avalue[rowiter]);
    }

    // check if this is an equation row and it now has a different size
    reinsertEquation(colrow);
    // printf("after substitution: ");
    // debugPrintRow(colrow);
  }

  assert(colsize[col] == 1);

  // substitute column in the objective function
  if (model->col_cost_[col] != 0.0) {
    HighsCDouble objscale =
        static_cast<HighsCDouble>(model->col_cost_[col]) * substrowscale;
    model->offset_ = static_cast<double>(model->offset_ - objscale * rhs);
    assert(std::isfinite(model->offset_));
    for (HighsInt rowiter : rowpositions) {
      // printf("changing col cost to %g = %g + %g * %g\n",
      // double(model->col_cost_[Acol[rowiter]] + objscale * Avalue[rowiter]),
      // model->col_cost_[Acol[rowiter]], double(objscale), Avalue[rowiter]);
      model->col_cost_[Acol[rowiter]] = static_cast<double>(
          model->col_cost_[Acol[rowiter]] + objscale * Avalue[rowiter]);
      if (std::abs(model->col_cost_[Acol[rowiter]]) <=
          options->small_matrix_value)
        model->col_cost_[Acol[rowiter]] = 0.0;
    }
    assert(std::abs(model->col_cost_[col]) <=
           std::max(options->dual_feasibility_tolerance,
                    static_cast<double>(kHighsTiny * abs(objscale))));
    model->col_cost_[col] = 0.0;
  }

  // finally remove the entries of the row that was used for substitution
  for (HighsInt rowiter : rowpositions) unlink(rowiter);
}

}  // namespace presolve
