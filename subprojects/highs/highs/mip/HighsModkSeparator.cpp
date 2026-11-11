/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsModKSeparator.cpp
 */

#include "mip/HighsModkSeparator.h"

#include <unordered_set>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsGFkSolve.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"
#include "util/HighsHash.h"
#include "util/HighsIntegers.h"

template <HighsInt k, typename FoundModKCut>
static bool separateModKCuts(const std::vector<int64_t>& intSystemValue,
                             const std::vector<HighsInt>& intSystemIndex,
                             const std::vector<HighsInt>& intSystemStart,
                             const HighsCutPool& cutpool, HighsInt numCol,
                             FoundModKCut&& foundModKCut) {
  HighsGFkSolve GFkSolve;

  HighsInt numCuts = cutpool.getNumCuts();

  GFkSolve.fromCSC<k>(intSystemValue, intSystemIndex, intSystemStart,
                      numCol + 1);
  GFkSolve.setRhs<k>(numCol, 1);
  GFkSolve.solve<k>(foundModKCut);

  return cutpool.getNumCuts() != numCuts;
}

void HighsModkSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& lpAggregator,
                                            HighsTransformedLp& transLp,
                                            HighsCutPool& cutpool) {
  const HighsMipSolver& mipsolver = lpRelaxation.getMipSolver();
  const HighsLp& lp = lpRelaxation.getLp();

  std::vector<HighsBool> skipRow(lp.num_row_);

  // mark all rows that have continuous variables with a nonzero solution value
  // in the transformed LP to be skipped
  for (HighsInt col : mipsolver.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0) continue;

    const HighsInt start = lp.a_matrix_.start_[col];
    const HighsInt end = lp.a_matrix_.start_[col + 1];

    for (HighsInt i = start; i != end; ++i)
      skipRow[lp.a_matrix_.index_[i]] = true;
  }

  HighsCutGeneration cutGen(lpRelaxation, cutpool);

  std::vector<std::pair<HighsInt, double>> integralScales;
  std::vector<int64_t> intSystemValue;
  std::vector<HighsInt> intSystemIndex;
  std::vector<HighsInt> intSystemStart;

  intSystemValue.reserve(lp.a_matrix_.value_.size() + lp.num_row_);
  intSystemIndex.reserve(intSystemValue.size());
  intSystemStart.reserve(lp.num_row_ + 1);
  intSystemStart.push_back(0);

  std::vector<HighsInt> inds;
  std::vector<double> vals;
  std::vector<double> scaleVals;

  inds.reserve(lp.num_col_);
  vals.reserve(lp.num_col_);
  scaleVals.reserve(lp.num_col_);

  std::vector<double> upper;
  std::vector<double> solval;
  double rhs;

  const HighsSolution& lpSolution = lpRelaxation.getSolution();
  HighsInt numNonzeroRhs = 0;
  HighsInt maxIntRowLen = 1000 + 0.1 * lp.num_col_;

  for (HighsInt row = 0; row != lp.num_row_; ++row) {
    if (skipRow[row]) continue;

    bool leqRow;

    if (lp.row_upper_[row] - lpSolution.row_value[row] <=
        mipsolver.mipdata_->feastol)
      leqRow = true;
    else if (lpSolution.row_value[row] - lp.row_lower_[row] <=
             mipsolver.mipdata_->feastol)
      leqRow = false;
    else
      continue;

    HighsInt rowlen;
    const HighsInt* rowinds;
    const double* rowvals;

    lpRelaxation.getRow(row, rowlen, rowinds, rowvals);

    if (leqRow) {
      rhs = lp.row_upper_[row];
      inds.assign(rowinds, rowinds + rowlen);
      vals.assign(rowvals, rowvals + rowlen);
    } else {
      assert(lpSolution.row_value[row] - lp.row_lower_[row] <=
             mipsolver.mipdata_->feastol);

      rhs = -lp.row_lower_[row];
      inds.assign(rowinds, rowinds + rowlen);
      vals.resize(rowlen);
      std::transform(rowvals, rowvals + rowlen, vals.begin(),
                     [](double x) { return -x; });
    }

    bool integralPositive = false;
    if (!transLp.transform(vals, upper, solval, inds, rhs, integralPositive,
                           true))
      continue;

    rowlen = inds.size();
    if (rowlen > maxIntRowLen) {
      HighsInt intRowLen = 0;
      for (HighsInt i = 0; i < rowlen; ++i) {
        if (solval[i] <= mipsolver.mipdata_->feastol) continue;
        if (mipsolver.isColContinuous(inds[i])) continue;
        ++intRowLen;
      }

      // skip row if either too long or 0 = 0 row
      if (intRowLen > maxIntRowLen ||
          (intRowLen == 0 && fabs(rhs) <= mipsolver.mipdata_->epsilon))
        continue;
    }

    double intscale;
    int64_t intrhs;

    if (!lpRelaxation.isRowIntegral(row)) {
      scaleVals.clear();
      for (HighsInt i = 0; i != rowlen; ++i) {
        if (mipsolver.isColContinuous(inds[i])) continue;
        if (solval[i] > mipsolver.mipdata_->feastol) {
          scaleVals.push_back(vals[i]);
        }
      }

      if (fabs(rhs) > mipsolver.mipdata_->epsilon) scaleVals.push_back(-rhs);
      if (scaleVals.empty()) continue;

      intscale = HighsIntegers::integralScale(
          scaleVals, mipsolver.mipdata_->feastol, mipsolver.mipdata_->epsilon);
      if (intscale == 0.0 || intscale > 1e6) continue;

      intrhs = HighsIntegers::nearestInteger(intscale * rhs);

      for (HighsInt i = 0; i != rowlen; ++i) {
        if (mipsolver.isColContinuous(inds[i])) continue;
        if (solval[i] > mipsolver.mipdata_->feastol) {
          intSystemIndex.push_back(inds[i]);
          intSystemValue.push_back(
              HighsIntegers::nearestInteger(intscale * vals[i]));
        }
      }
    } else {
      intscale = 1.0;
      intrhs = HighsIntegers::nearestInteger(rhs);

      for (HighsInt i = 0; i != rowlen; ++i) {
        if (solval[i] > mipsolver.mipdata_->feastol) {
          intSystemIndex.push_back(inds[i]);
          intSystemValue.push_back(HighsIntegers::nearestInteger(vals[i]));
        }
      }
    }

    numNonzeroRhs += (intrhs != 0);

    intSystemIndex.push_back(lp.num_col_);
    intSystemValue.push_back(intrhs);
    intSystemStart.push_back(intSystemValue.size());
    integralScales.emplace_back(row, intscale);
  }

  if (integralScales.empty() || numNonzeroRhs == 0) return;

  std::vector<HighsInt> tmpinds;
  std::vector<double> tmpvals;

  HighsHashTable<std::vector<HighsGFkSolve::SolutionEntry>> usedWeights;
  // std::unordered_set<std::vector<HighsGFkSolve::SolutionEntry>,
  //                   HighsVectorHasher, HighsVectorEqual>
  //    usedWeights;
  HighsInt k;
  auto foundCut = [&](std::vector<HighsGFkSolve::SolutionEntry>& weights,
                      int rhsIndex) {
    // cuts which come from a single row can already be found with the
    // aggregation heuristic
    if (weights.empty()) return;

    pdqsort(weights.begin(), weights.end());
    if (!usedWeights.insert(weights)) return;

    assert(lpAggregator.isEmpty());
    for (const auto& w : weights) {
      double weight = integralScales[w.index].second *
                      (double((w.weight * (k - 1)) % k) / k);
      HighsInt row = integralScales[w.index].first;
      lpAggregator.addRow(row, weight);
    }

    lpAggregator.getCurrentAggregation(inds, vals, false);

    rhs = 0.0;
    cutGen.generateCut(transLp, inds, vals, rhs, true);

    if (k != 2) {
      lpAggregator.clear();
      for (const auto& w : weights) {
        double weight = integralScales[w.index].second * (double(w.weight) / k);
        HighsInt row = integralScales[w.index].first;
        lpAggregator.addRow(row, weight);
      }
    }

    lpAggregator.getCurrentAggregation(inds, vals, true);

    rhs = 0.0;
    cutGen.generateCut(transLp, inds, vals, rhs, true);

    lpAggregator.clear();
  };

  k = 2;
  if (separateModKCuts<2>(intSystemValue, intSystemIndex, intSystemStart,
                          cutpool, lp.num_col_, foundCut))
    return;

  usedWeights.clear();
  k = 3;
  if (separateModKCuts<3>(intSystemValue, intSystemIndex, intSystemStart,
                          cutpool, lp.num_col_, foundCut))
    return;

  usedWeights.clear();
  k = 5;
  if (separateModKCuts<5>(intSystemValue, intSystemIndex, intSystemStart,
                          cutpool, lp.num_col_, foundCut))
    return;

  usedWeights.clear();
  k = 7;
  if (separateModKCuts<7>(intSystemValue, intSystemIndex, intSystemStart,
                          cutpool, lp.num_col_, foundCut))
    return;
}
