/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTableauSeparator.cpp
 */

#include "mip/HighsTableauSeparator.h"

#include <algorithm>

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

struct FractionalInteger {
  double fractionality;
  double row_ep_norm2;
  double score;
  HighsInt basisIndex;
  std::vector<std::pair<HighsInt, double>> row_ep;

  bool operator<(const FractionalInteger& other) const {
    return score > other.score;
  }

  FractionalInteger() = default;

  FractionalInteger(HighsInt basisIndex, double fractionality)
      : fractionality(fractionality), score(-1.0), basisIndex(basisIndex) {}
};

void HighsTableauSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                               HighsLpAggregator& lpAggregator,
                                               HighsTransformedLp& transLp,
                                               HighsCutPool& cutpool) {
  Highs& lpSolver = lpRelaxation.getLpSolver();
  if (!lpSolver.hasInvert()) return;

  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  if (cutpool.getNumAvailableCuts() > mip.options_mip_->mip_pool_soft_limit)
    return;

  const HighsInt* basisinds =
      lpRelaxation.getLpSolver().getBasicVariablesArray();
  HighsInt numRow = lpRelaxation.numRows();
  HighsInt numCol = lpRelaxation.numCols();

  HighsCutGeneration cutGen(lpRelaxation, cutpool);

  std::vector<HighsInt> baseRowInds;
  std::vector<double> baseRowVals;

  const HighsSolution& lpSolution = lpRelaxation.getSolution();

  std::vector<FractionalInteger> fractionalBasisvars;
  fractionalBasisvars.reserve(numRow);
  for (HighsInt i = 0; i < numRow; ++i) {
    double my_fractionality;
    if (basisinds[i] >= numCol) {
      HighsInt row = basisinds[i] - numCol;

      if (!lpRelaxation.isRowIntegral(row)) continue;

      my_fractionality = fractionality(lpSolution.row_value[row]);
    } else {
      HighsInt col = basisinds[i];
      if (mip.isColContinuous(col)) continue;

      my_fractionality = fractionality(lpSolution.col_value[col]);
    }

    if (my_fractionality < 1000 * mip.mipdata_->feastol) continue;

    fractionalBasisvars.emplace_back(i, my_fractionality);
  }

  if (fractionalBasisvars.empty()) return;
  int64_t maxTries = 5000 + getNumCalls() * 50 +
                     (mip.mipdata_->total_lp_iterations -
                      mip.mipdata_->heuristic_lp_iterations) /
                         10;
  if (numTries >= maxTries) return;

  maxTries -= numTries;

  maxTries = std::min(
      {maxTries,
       200 + int64_t(0.1 *
                     std::min(numRow,
                              (HighsInt)mip.mipdata_->integral_cols.size()))});

  if ((HighsInt)fractionalBasisvars.size() > maxTries) {
    const double* edgeWt = lpRelaxation.getLpSolver().getDualEdgeWeights();
    if (edgeWt) {
      // printf("choosing %ld/%zu with DSE weights\n", maxTries,
      // fractionalBasisvars.size());
      pdqsort(
          fractionalBasisvars.begin(), fractionalBasisvars.end(),
          [&](const FractionalInteger& fracint1,
              const FractionalInteger& fracint2) {
            double score1 = fracint1.fractionality *
                            (1.0 - fracint1.fractionality) /
                            edgeWt[fracint1.basisIndex];
            double score2 = fracint2.fractionality *
                            (1.0 - fracint2.fractionality) /
                            edgeWt[fracint2.basisIndex];
            return std::make_pair(score1, HighsHashHelpers::hash(
                                              numTries + fracint1.basisIndex)) >
                   std::make_pair(score2, HighsHashHelpers::hash(
                                              numTries + fracint2.basisIndex));
          });
    } else {
      // printf("choosing %ld/%zu without DSE weights\n", maxTries,
      // fractionalBasisvars.size());
      pdqsort(
          fractionalBasisvars.begin(), fractionalBasisvars.end(),
          [&](const FractionalInteger& fracint1,
              const FractionalInteger& fracint2) {
            return std::make_pair(
                       fracint1.fractionality,
                       HighsHashHelpers::hash(numTries + fracint1.basisIndex)) >
                   std::make_pair(
                       fracint2.fractionality,
                       HighsHashHelpers::hash(numTries + fracint2.basisIndex));
          });
    }

    fractionalBasisvars.resize(maxTries);
  }

  HVector rowEpBuffer;
  rowEpBuffer.setup(numRow);

  numTries += fractionalBasisvars.size();

  for (auto& fracvar : fractionalBasisvars) {
    if (lpSolver.getBasisInverseRowSparse(fracvar.basisIndex, rowEpBuffer) !=
        HighsStatus::kOk)
      continue;

    // handled by other separator
    if (rowEpBuffer.count == 1) continue;

    fracvar.row_ep_norm2 = 0.0;
    double minWeight = kHighsInf;
    double maxWeight = 0.0;
    fracvar.row_ep.reserve(rowEpBuffer.count);
    for (HighsInt j = 0; j < rowEpBuffer.count; ++j) {
      HighsInt row = rowEpBuffer.index[j];
      double weight = rowEpBuffer.array[row];
      double maxAbsRowVal = lpRelaxation.getMaxAbsRowVal(row);

      double scaledWeight = maxAbsRowVal * std::abs(weight);
      if (scaledWeight <= mip.mipdata_->feastol) continue;

      minWeight = std::min(minWeight, scaledWeight);
      maxWeight = std::max(maxWeight, scaledWeight);
      fracvar.row_ep_norm2 += scaledWeight * scaledWeight;
      fracvar.row_ep.emplace_back(row, weight);
    }

    if (fracvar.row_ep.size() <= 1) continue;

    if (maxWeight / minWeight <= 1e4) {
      fracvar.score = fracvar.fractionality * (1.0 - fracvar.fractionality) /
                      fracvar.row_ep_norm2;
    }
  }

  fractionalBasisvars.erase(
      std::remove_if(fractionalBasisvars.begin(), fractionalBasisvars.end(),
                     [&](const FractionalInteger& fracInteger) {
                       return fracInteger.score <= mip.mipdata_->feastol;
                     }),
      fractionalBasisvars.end());

  if (fractionalBasisvars.empty()) return;

  pdqsort_branchless(fractionalBasisvars.begin(), fractionalBasisvars.end());
  double bestScore = -1.0;

  HighsInt numCuts = cutpool.getNumCuts();
  const double bestScoreFac[] = {0.0025, 0.01};

  for (const auto& fracvar : fractionalBasisvars) {
    if (cutpool.getNumCuts() - numCuts >= 1000) break;

    if (fracvar.score <
        bestScoreFac[cutpool.getNumCuts() - numCuts >= 50] * bestScore)
      break;

    assert(lpAggregator.isEmpty());
    for (const auto& rowWeight : fracvar.row_ep)
      lpAggregator.addRow(rowWeight.first, rowWeight.second);

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);

    if (10 * (baseRowInds.size() - fracvar.row_ep.size()) >
        10000 + static_cast<size_t>(mip.numCol())) {
      lpAggregator.clear();
      continue;
    }

    HighsInt len = baseRowInds.size();
    if (len > (HighsInt)fracvar.row_ep.size()) {
      double maxAbsVal = 0.0;
      double minAbsVal = kHighsInf;
      for (HighsInt i = 0; i < len; ++i) {
        if (baseRowInds[i] < mip.numCol()) {
          maxAbsVal = std::max(std::abs(baseRowVals[i]), maxAbsVal);
          minAbsVal = std::min(std::abs(baseRowVals[i]), minAbsVal);
        }
      }
      if (maxAbsVal / minAbsVal > 1e6) {
        lpAggregator.clear();
        continue;
      }
    }

    mip.mipdata_->debugSolution.checkRowAggregation(
        lpSolver.getLp(), baseRowInds.data(), baseRowVals.data(),
        baseRowInds.size());

    double rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);
    if (transLp.getGlobaldom().infeasible()) break;

    lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
    rhs = 0;
    cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);
    if (transLp.getGlobaldom().infeasible()) break;

    lpAggregator.clear();
    if (bestScore == -1.0 && cutpool.getNumCuts() != numCuts)
      bestScore = fracvar.score;
  }
}
