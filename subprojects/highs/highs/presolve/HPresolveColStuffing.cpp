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

HPresolve::Result HPresolve::singletonColStuffing(
    HighsPostsolveStack& postsolve_stack, HighsInt col) {
  // singleton column stuffing
  // see Gamrath, G., Koch, T., Martin, A. et al., Progress in presolving
  // for mixed integer programming, Math. Prog. Comp. 7, 367–398 (2015).

  // count number of fixed columns
  HighsInt numFixedCols = 0;

  struct candidate {
    HighsInt col;
    double val;
    HighsInt multiplier;
  };

  auto isSingleton = [&](HighsInt col) {
    return (!colDeleted[col] && colsize[col] == 1 &&
            model->col_lower_[col] != model->col_upper_[col]);
  };

  auto sortCols = [&](std::vector<candidate>& vec) {
    pdqsort(vec.begin(), vec.end(),
            [&](const candidate& c1, const candidate& c2) {
              return model->col_cost_[c1.col] / c1.val <
                     model->col_cost_[c2.col] / c2.val;
            });
  };

  // lambda for updating row activity bounds
  auto updateActivityBounds = [&](HighsCDouble& sumLower,
                                  HighsCDouble& sumUpper, bool& sumLowerFinite,
                                  bool& sumUpperFinite, double aj,
                                  double lowerSumBound, double upperSumBound) {
    sumLowerFinite = sumLowerFinite && std::abs(lowerSumBound) != kHighsInf;
    sumUpperFinite = sumUpperFinite && std::abs(upperSumBound) != kHighsInf;
    if (sumLowerFinite)
      sumLower += aj * static_cast<HighsCDouble>(lowerSumBound);
    if (sumUpperFinite)
      sumUpper += aj * static_cast<HighsCDouble>(upperSumBound);
  };

  // lambda for storing a candidate
  auto addCandidate = [&](std::vector<candidate>& candidates, HighsInt col,
                          double val, HighsInt direction, double& minWeight,
                          double& maxWeight, size_t& numIntegerCandidates) {
    if (model->integrality_[col] == HighsVarType::kInteger)
      numIntegerCandidates++;

    minWeight = std::min(minWeight, direction * val);
    maxWeight = std::max(maxWeight, direction * val);
    candidates.push_back(candidate{col, val, direction});
  };

  // lambda for fixing a variable
  auto fixCol = [&](HighsInt col, HighsInt direction) {
    if (direction < 0)
      HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
    else
      HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
    return Result::kOk;
  };

  // lambda for computing candidates for stuffing
  auto computeCandidates = [&](HighsInt row, HighsInt direction,
                               std::vector<candidate>& candidates,
                               HighsCDouble& sumLower, HighsCDouble& sumUpper,
                               bool& sumLowerFinite, bool& sumUpperFinite,
                               size_t& numIntegerCandidates, double& minWeight,
                               double& maxWeight, bool allowIntegerCandidates) {
    // vectors for candidates and activity bounds
    candidates.clear();
    candidates.reserve(rowsize[row]);
    sumLower = 0.0;
    sumUpper = 0.0;
    sumLowerFinite = true;
    sumUpperFinite = true;
    numIntegerCandidates = 0;
    minWeight = kHighsInf;
    maxWeight = -kHighsInf;

    for (auto& nz : getRowVector(row)) {
      // get column index, coefficient, cost and bounds
      HighsInt j = nz.index();
      double aj = direction * nz.value();
      double cj = model->col_cost_[j];
      double sumLowerBound = model->col_lower_[j];
      double sumUpperBound = model->col_upper_[j];
      bool isCandidate = allowIntegerCandidates ||
                         model->integrality_[j] != HighsVarType::kInteger;
      if (isSingleton(j)) {
        // check singleton
        if (aj > 0) {
          if (cj >= 0)
            // dual fixing: fix to lower bound
            sumUpperBound = sumLowerBound;
          else if (isCandidate) {
            // candidate for stuffing
            sumUpperBound = sumLowerBound;
            addCandidate(candidates, j, aj, HighsInt{1}, minWeight, maxWeight,
                         numIntegerCandidates);
          }
        } else {
          // aj < 0
          assert(aj < 0);
          if (cj <= 0)
            // dual fixing: fix to upper bound
            sumLowerBound = sumUpperBound;
          else if (isCandidate) {
            // candidate for stuffing; multiply column with -1
            sumLowerBound = sumUpperBound;
            addCandidate(candidates, j, aj, HighsInt{-1}, minWeight, maxWeight,
                         numIntegerCandidates);
          }
        }
      }
      // update activities
      if (aj < 0) std::swap(sumLowerBound, sumUpperBound);
      updateActivityBounds(sumLower, sumUpper, sumLowerFinite, sumUpperFinite,
                           aj, sumLowerBound, sumUpperBound);
      if (!sumLowerFinite && !sumUpperFinite) return false;
    }
    return true;
  };

  // lambda for computing and checking candidates for stuffing
  auto checkCandidates = [&](HighsInt row, HighsInt direction,
                             std::vector<candidate>& candidates,
                             HighsCDouble& sumLower, HighsCDouble& sumUpper,
                             bool& sumLowerFinite, bool& sumUpperFinite) {
    // number of integer candidates and weights
    size_t numIntegerCandidates;
    double minWeight;
    double maxWeight;

    // compute candidates
    if (computeCandidates(row, direction, candidates, sumLower, sumUpper,
                          sumLowerFinite, sumUpperFinite, numIntegerCandidates,
                          minWeight, maxWeight, true)) {
      // return if there are no integer columns
      if (numIntegerCandidates == 0) return true;
      // all columns need to have same weights if we only have integer
      // columns
      if (numIntegerCandidates == candidates.size()) {
        if (minWeight != maxWeight) return false;
        return true;
      }
    }

    // recompute candidates without integer columns
    if (numIntegerCandidates > 0 &&
        computeCandidates(row, direction, candidates, sumLower, sumUpper,
                          sumLowerFinite, sumUpperFinite, numIntegerCandidates,
                          minWeight, maxWeight, false))
      return true;

    return false;
  };

  // lambda for actual stuffing
  auto checkRow = [&](HighsInt row, double rhs, HighsInt direction) {
    // skip row if rhs is not finite
    if (direction * rhs == kHighsInf) return Result::kOk;

    // vectors for candidates and activity bounds
    std::vector<candidate> candidates;
    HighsCDouble sumLower;
    HighsCDouble sumUpper;
    bool sumLowerFinite;
    bool sumUpperFinite;

    // compute candidates
    if (!checkCandidates(row, direction, candidates, sumLower, sumUpper,
                         sumLowerFinite, sumUpperFinite))
      return Result::kOk;

    // sort candidates
    sortCols(candidates);

    // check candidates
    for (const auto& t : candidates) {
      // both bounds have to be finite
      if (model->col_lower_[t.col] == -kHighsInf ||
          model->col_upper_[t.col] == kHighsInf)
        break;
      // compute delta (bound difference)
      HighsCDouble delta =
          t.multiplier * t.val *
          (static_cast<HighsCDouble>(model->col_upper_[t.col]) -
           static_cast<HighsCDouble>(model->col_lower_[t.col]));
      // check if variable can be fixed
      if (sumUpperFinite &&
          delta - (direction * rhs - sumUpper) <= primal_feastol) {
        numFixedCols++;
        HPRESOLVE_CHECKED_CALL(fixCol(t.col, t.multiplier));
      } else if (sumLowerFinite &&
                 delta - (sumLower - direction * rhs) <= primal_feastol) {
        // Previously
        //
        // direction * rhs <= sumLower + primal_feastol
        //
        // But this led to fixing at -t.multiplier until row activity
        // was at its bound, which is primal optmal but degenerate,
        // and does not yield an optimal basis if column costs are
        // positive
        numFixedCols++;
        HPRESOLVE_CHECKED_CALL(fixCol(t.col, -t.multiplier));
      }
      // update row activities
      if (sumLowerFinite) sumLower += delta;
      if (sumUpperFinite) sumUpper += delta;
    }

    return Result::kOk;
  };

  // consider only non-fixed singleton columns
  if (!isSingleton(col)) return Result::kOk;

  // get row index
  HighsInt row = Arow[colhead[col]];

  // return if we have an empty or singleton row, a large row, or the row is
  // ranged
  constexpr HighsInt maxStuffingRowSize = 1000;
  if (rowsize[row] <= 1 || rowsize[row] > maxStuffingRowSize || isRanged(row))
    return Result::kOk;

  // check row
  HPRESOLVE_CHECKED_CALL(checkRow(row, model->row_upper_[row], HighsInt{1}));
  HPRESOLVE_CHECKED_CALL(checkRow(row, model->row_lower_[row], HighsInt{-1}));

  if (numFixedCols > 0)
    highsLogDev(options->log_options, HighsLogType::kDetailed,
                "Singleton column stuffing fixed %d columns\n",
                static_cast<int>(numFixedCols));

  return Result::kOk;
}

}  // namespace presolve
