/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsPathSeparator.cpp
 */

#include "mip/HighsPathSeparator.h"

#include "mip/HighsCutGeneration.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsTransformedLp.h"

enum class RowType : int8_t {
  kUnusuable = -2,
  kGeq = -1,
  kEq = 0,
  kLeq = 1,
};

void HighsPathSeparator::separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                            HighsLpAggregator& lpAggregator,
                                            HighsTransformedLp& transLp,
                                            HighsCutPool& cutpool) {
  const HighsMipSolver& mip = lpRelaxation.getMipSolver();
  const HighsLp& lp = lpRelaxation.getLp();
  const HighsSolution& lpSolution = lpRelaxation.getSolution();

  std::vector<RowType> rowtype;
  rowtype.resize(lp.num_row_);
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    if (lp.row_lower_[i] == lp.row_upper_[i]) {
      rowtype[i] = RowType::kEq;
      continue;
    }

    double lowerslack = kHighsInf;
    double upperslack = kHighsInf;

    if (lp.row_lower_[i] != -kHighsInf)
      lowerslack = lpSolution.row_value[i] - lp.row_lower_[i];

    if (lp.row_upper_[i] != kHighsInf)
      upperslack = lp.row_upper_[i] - lpSolution.row_value[i];

    if (lowerslack > mip.mipdata_->feastol &&
        upperslack > mip.mipdata_->feastol)
      rowtype[i] = RowType::kUnusuable;
    else if (lowerslack < upperslack)
      rowtype[i] = RowType::kGeq;
    else
      rowtype[i] = RowType::kLeq;
  }

  std::vector<HighsInt> numContinuous(lp.num_row_);

  size_t maxAggrRowSize = 0;
  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;

    maxAggrRowSize += lp.a_matrix_.start_[col + 1] - lp.a_matrix_.start_[col];
    for (HighsInt i = lp.a_matrix_.start_[col];
         i != lp.a_matrix_.start_[col + 1]; ++i)
      ++numContinuous[lp.a_matrix_.index_[i]];
  }

  std::vector<std::pair<HighsInt, double>> colSubstitutions(
      lp.num_col_, std::make_pair(-1, 0.0));

  // identify equality rows where only a single continuous variable with nonzero
  // transformed solution value is present. Mark those columns and remember the
  // rows so that we can always substitute such columns away using this equation
  // and block the equation from being used as a start row
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    if (rowtype[i] != RowType::kEq || numContinuous[i] != 1) continue;

    HighsInt len;
    const HighsInt* rowinds;
    const double* rowvals;
    lpRelaxation.getRow(i, len, rowinds, rowvals);

    // find continuous variable
    HighsInt col = -1;
    double val = 0.0;
    for (HighsInt j = 0; j != len; ++j) {
      if (mip.isColIntegral(rowinds[j])) continue;
      if (transLp.boundDistance(rowinds[j]) == 0.0) continue;
      col = rowinds[j];
      val = rowvals[j];
      break;
    }

    assert(col != -1);
    assert(mip.isColContinuous(col));
    assert(transLp.boundDistance(col) > 0.0);

    if (colSubstitutions[col].first != -1) continue;

    colSubstitutions[col].first = i;
    colSubstitutions[col].second = val;
    rowtype[i] = RowType::kUnusuable;
  }

  // for each continuous variable with nonzero transformed solution value
  // remember the <= and == rows where it is present with a positive coefficient
  // in its set of in-arc rows. Treat >= rows as <= rows with reversed sign
  // The reason to only store one set of rows for one sign of the coefficients
  // is that this directs the selection to be more diverse. Consider
  // aggregations of 2 rows where we allow both directions. When one of the rows
  // is used as start row we can always select the other one. When we only
  // project out variables with negative coefficients we give the aggregation
  // path an orientation and avoid this situation. For each aggregation of rows
  // the cut generation will try the reversed orientation in any case too.

  std::vector<std::pair<HighsInt, double>> inArcRows;
  inArcRows.reserve(maxAggrRowSize);
  std::vector<std::pair<HighsInt, HighsInt>> colInArcs(lp.num_col_);

  std::vector<std::pair<HighsInt, double>> outArcRows;
  outArcRows.reserve(maxAggrRowSize);
  std::vector<std::pair<HighsInt, HighsInt>> colOutArcs(lp.num_col_);

  for (HighsInt col : mip.mipdata_->continuous_cols) {
    if (transLp.boundDistance(col) == 0.0) continue;
    if (colSubstitutions[col].first != -1) continue;

    colInArcs[col].first = inArcRows.size();
    colOutArcs[col].first = outArcRows.size();
    for (HighsInt i = lp.a_matrix_.start_[col];
         i != lp.a_matrix_.start_[col + 1]; ++i) {
      switch (rowtype[lp.a_matrix_.index_[i]]) {
        case RowType::kUnusuable:
          continue;
        case RowType::kLeq:
          if (lp.a_matrix_.value_[i] < 0)
            inArcRows.emplace_back(lp.a_matrix_.index_[i],
                                   lp.a_matrix_.value_[i]);
          else
            outArcRows.emplace_back(lp.a_matrix_.index_[i],
                                    lp.a_matrix_.value_[i]);
          break;
        case RowType::kGeq:
          if (lp.a_matrix_.value_[i] > 0)
            inArcRows.emplace_back(lp.a_matrix_.index_[i],
                                   lp.a_matrix_.value_[i]);
          else
            outArcRows.emplace_back(lp.a_matrix_.index_[i],
                                    lp.a_matrix_.value_[i]);
          break;
        case RowType::kEq:
          inArcRows.emplace_back(lp.a_matrix_.index_[i],
                                 lp.a_matrix_.value_[i]);
          outArcRows.emplace_back(lp.a_matrix_.index_[i],
                                  lp.a_matrix_.value_[i]);
      }
    }

    colInArcs[col].second = inArcRows.size();
    colOutArcs[col].second = outArcRows.size();
  }

  HighsCutGeneration cutGen(lpRelaxation, cutpool);
  std::vector<HighsInt> baseRowInds;
  std::vector<double> baseRowVals;
  constexpr HighsInt maxPathLen = 6;
  HighsInt currentPath[maxPathLen];
  std::vector<std::pair<std::vector<HighsInt>, std::vector<double>>>
      aggregatedPath;
  std::array<double, 2> scales;
  for (HighsInt i = 0; i != lp.num_row_; ++i) {
    switch (rowtype[i]) {
      case RowType::kUnusuable:
        continue;
      case RowType::kEq:
        if (lpSolution.row_dual[i] <= mip.mipdata_->epsilon) {
          scales[0] = 1.0;
          scales[1] = -1.0;
        } else {
          scales[0] = -1.0;
          scales[1] = 1.0;
        }
        break;
      case RowType::kLeq:
        scales[0] = 1.0;
        scales[1] = -1.0;
        break;
      case RowType::kGeq:
        scales[0] = -1.0;
        scales[1] = 1.0;
        break;
      default:
        assert(false);
    }

    bool success = false;

    for (HighsInt s = 0; s != 2; ++s) {
      assert(lpAggregator.isEmpty());
      lpAggregator.addRow(i, scales[s]);

      currentPath[0] = i;
      HighsInt currPathLen = 1;

      auto isRowInCurrentPath = [&](HighsInt row) {
        for (HighsInt j = 0; j < currPathLen; ++j)
          if (currentPath[j] == row) return true;
        return false;
      };

      bool tryNegatedScale = false;
      const double maxWeight = 1. / mip.mipdata_->feastol;
      const double minWeight = mip.mipdata_->feastol;

      auto checkWeight = [&](double w) {
        w = std::abs(w);
        return w <= maxWeight && w >= minWeight;
      };

      auto skipCol =
          [&](const HighsInt& col,
              const std::vector<std::pair<HighsInt, HighsInt>>& colArcs,
              const std::vector<std::pair<HighsInt, double>>& arcRows,
              const std::vector<std::pair<HighsInt, HighsInt>>& otherColArcs,
              const std::vector<std::pair<HighsInt, double>>& otherArcRows) {
            if (currPathLen == 1 && !tryNegatedScale) {
              if (colArcs[col].second - colArcs[col].first <= currPathLen) {
                for (HighsInt k = colArcs[col].first; k < colArcs[col].second;
                     ++k) {
                  if (arcRows[k].first != i) {
                    tryNegatedScale = true;
                    break;
                  }
                }
              } else
                tryNegatedScale = true;
            }

            if (otherColArcs[col].first == otherColArcs[col].second)
              return true;
            if (otherColArcs[col].second - otherColArcs[col].first <=
                currPathLen) {
              for (HighsInt k = otherColArcs[col].first;
                   k < otherColArcs[col].second; ++k) {
                if (!isRowInCurrentPath(otherArcRows[k].first)) return false;
              }
              return true;
            }
            return false;
          };

      auto findRow =
          [&](const HighsInt& bestArcCol, const double& val,
              const std::vector<std::pair<HighsInt, HighsInt>>& colArcs,
              const std::vector<std::pair<HighsInt, double>>& arcRows,
              HighsInt& row, double& weight) {
            HighsInt arcRow = randgen.integer(colArcs[bestArcCol].first,
                                              colArcs[bestArcCol].second);
            HighsInt r = arcRows[arcRow].first;
            double w = -val / arcRows[arcRow].second;
            if (!isRowInCurrentPath(r) && checkWeight(w)) {
              row = r;
              weight = w;
              return true;
            }

            for (HighsInt nextRow = arcRow + 1;
                 nextRow < colArcs[bestArcCol].second; ++nextRow) {
              r = arcRows[nextRow].first;
              w = -val / arcRows[nextRow].second;
              if (!isRowInCurrentPath(r) && checkWeight(w)) {
                row = r;
                weight = w;
                return true;
              }
            }

            for (HighsInt nextRow = colArcs[bestArcCol].first; nextRow < arcRow;
                 ++nextRow) {
              r = arcRows[nextRow].first;
              w = -val / arcRows[nextRow].second;
              if (!isRowInCurrentPath(r) && checkWeight(w)) {
                row = r;
                weight = w;
                return true;
              }
            }

            return false;
          };

      aggregatedPath.clear();

      while (currPathLen != maxPathLen) {
        lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, false);
        HighsInt baseRowLen = baseRowInds.size();
        bool addedSubstitutionRows = false;

        HighsInt bestOutArcCol = -1;
        double outArcColVal = 0.0;
        double outArcColBoundDist = 0.0;

        HighsInt bestInArcCol = -1;
        double inArcColVal = 0.0;
        double inArcColBoundDist = 0.0;

        for (HighsInt j = 0; j != baseRowLen; ++j) {
          HighsInt col = baseRowInds[j];
          if (col >= lp.num_col_ || transLp.boundDistance(col) == 0.0 ||
              lpRelaxation.isColIntegral(col))
            continue;

          if (colSubstitutions[col].first != -1) {
            addedSubstitutionRows = true;
            lpAggregator.addRow(colSubstitutions[col].first,
                                -baseRowVals[j] / colSubstitutions[col].second);
            continue;
          }

          if (addedSubstitutionRows) continue;

          if (baseRowVals[j] < 0) {
            if (skipCol(col, colOutArcs, outArcRows, colInArcs, inArcRows))
              continue;

            if (bestOutArcCol == -1 ||
                transLp.boundDistance(col) > outArcColBoundDist) {
              bestOutArcCol = col;
              outArcColVal = baseRowVals[j];
              outArcColBoundDist = transLp.boundDistance(col);
            }
          } else {
            if (skipCol(col, colInArcs, inArcRows, colOutArcs, outArcRows))
              continue;

            if (bestInArcCol == -1 ||
                transLp.boundDistance(col) > inArcColBoundDist) {
              bestInArcCol = col;
              inArcColVal = baseRowVals[j];
              inArcColBoundDist = transLp.boundDistance(col);
            }
          }
        }

        if (addedSubstitutionRows) continue;

        // generate cut
        double rhs = 0;
        success = cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

        lpAggregator.getCurrentAggregation(baseRowInds, baseRowVals, true);
        if (!aggregatedPath.empty() || bestOutArcCol != -1 ||
            bestInArcCol != -1)
          aggregatedPath.emplace_back(baseRowInds, baseRowVals);

        // generate reverse cut
        rhs = 0;
        success |= cutGen.generateCut(transLp, baseRowInds, baseRowVals, rhs);

        if (success || (bestOutArcCol == -1 && bestInArcCol == -1)) break;

        // we prefer to use an out edge if the bound distances are equal in
        // feasibility tolerance otherwise we choose an inArc. This tie
        // breaking is arbitrary, but we should direct the substitution to
        // prefer one direction to increase diversity.
        HighsInt row = -1;
        double weight = 0.0;
        if (bestInArcCol == -1 ||
            (bestOutArcCol != -1 &&
             outArcColBoundDist >= inArcColBoundDist - mip.mipdata_->feastol)) {
          if (!findRow(bestOutArcCol, outArcColVal, colInArcs, inArcRows, row,
                       weight)) {
            if (bestInArcCol == -1)
              break;
            else if (!findRow(bestInArcCol, inArcColVal, colOutArcs, outArcRows,
                              row, weight))
              break;
          }
        } else {
          if (!findRow(bestInArcCol, inArcColVal, colOutArcs, outArcRows, row,
                       weight))
            break;
        }

        currentPath[currPathLen] = row;
        lpAggregator.addRow(row, weight);
        ++currPathLen;
      }

      // if the path has length at least 2 try to separate a path mixing cut
      HighsInt pathLen = aggregatedPath.size();
      if (pathLen > 1) {
        // generate path mixing cut
        HighsHashTable<HighsInt, HighsInt> indexPos;

        std::vector<HighsInt> inds;
        std::vector<double> solval;
        std::vector<double> upper;
        std::vector<HighsBool> isIntegral;
        inds.reserve(lp.num_col_ + lp.num_row_);
        solval.reserve(lp.num_col_ + lp.num_row_);
        upper.reserve(lp.num_col_ + lp.num_row_);
        isIntegral.reserve(lp.num_col_ + lp.num_row_);

        std::vector<double> rhs(pathLen);
        std::vector<double> tmpUpper;
        std::vector<double> tmpSolval;

        double delta = 1.0;

        transLp.initMultiRowTransform();

        for (HighsInt k = 0; k < pathLen; ++k) {
          bool integralPositive = false;

          if (!transLp.transform(aggregatedPath[k].second, tmpUpper, tmpSolval,
                                 aggregatedPath[k].first, rhs[k],
                                 integralPositive, false, true)) {
            pathLen = k;
            break;
          }

          if (k == 0) {
            if (rhs[k] > kHighsTiny) {
              pathLen = k;
              break;
            }
            if (rhs[k] >= -mip.mipdata_->feastol) rhs[k] = 0;
          } else if (rhs[k] >= rhs[k - 1] - mip.mipdata_->feastol) {
            pathLen = k;
            break;
          }

          delta = std::max(std::abs(rhs[k]), delta);

          HighsInt len = aggregatedPath[k].first.size();
          for (HighsInt j = 0; j < len; ++j) {
            HighsInt index = aggregatedPath[k].first[j];
            HighsInt* pos = &indexPos[index];
            if (*pos == 0) {
              inds.push_back(index);
              solval.push_back(tmpSolval[j]);
              upper.push_back(tmpUpper[j]);
              isIntegral.push_back(lpRelaxation.isColIntegral(index));
              if (isIntegral.back())
                delta = std::max(std::abs(aggregatedPath[k].second[j]), delta);
              *pos = inds.size();
            } else {
              assert(inds[*pos - 1] == index);
              assert(solval[*pos - 1] == tmpSolval[j]);
              assert(upper[*pos - 1] == tmpUpper[j]);
              if (isIntegral[*pos - 1])
                delta = std::max(std::abs(aggregatedPath[k].second[j]), delta);
            }
          }
        }

        if (pathLen > 1) {
          delta = std::exp2(std::ceil(std::log2(delta + 1.0)));

          HighsInt numInds = inds.size();

          HighsCDouble cutRhs = 0.0;
          std::vector<double> cutVals(numInds);
          std::vector<double> maxFrac(numInds);
          std::vector<HighsCDouble> downSum(numInds);
          std::vector<HighsCDouble> fSum(numInds);

          double fLast = 0;
          double scale = -1.0 / delta;

          for (HighsInt k = 0; k < pathLen; ++k) {
            double f = rhs[k] * scale;
            HighsCDouble fDiff = HighsCDouble(f) - fLast;
            HighsInt len = aggregatedPath[k].first.size();
            cutRhs += fDiff;
            for (HighsInt j = 0; j < len; ++j) {
              HighsInt i = indexPos[aggregatedPath[k].first[j]] - 1;
              assert(i >= 0);

              double gj = aggregatedPath[k].second[j] * scale;

              switch (isIntegral[i]) {
                case 0:
                  cutVals[i] = std::max(cutVals[i], gj);
                  break;
                case 1: {
                  double gjdown = std::floor(gj);
                  double hj = gj - gjdown;
                  maxFrac[i] = std::max(maxFrac[i], hj);
                  downSum[i] += fDiff * gjdown;
                  fSum[i] += fDiff;

                  if (fSum[i] < maxFrac[i]) {
                    cutVals[i] = double(downSum[i] + fSum[i]);
                  } else {
                    cutVals[i] = double(downSum[i] + maxFrac[i]);
                  }
                  break;
                }
              }
            }

            if (k > 0) {
              double viol = double(cutRhs);
              for (HighsInt j = 0; j < numInds; ++j) {
                viol -= solval[j] * cutVals[j];
              }

              viol *= delta;

              if (viol > 10 * mip.mipdata_->feastol) {
                scale = -delta;
                double rhs = double(cutRhs) * scale;
                for (HighsInt j = 0; j < numInds; ++j) {
                  cutVals[j] *= scale;
                }

                for (HighsInt j = numInds - 1; j >= 0; --j) {
                  if (std::abs(cutVals[j]) <= mip.mipdata_->epsilon) {
                    --numInds;
                    std::swap(cutVals[j], cutVals[numInds]);
                    std::swap(inds[j], inds[numInds]);
                  }
                }

                cutVals.resize(numInds);
                inds.resize(numInds);

                if (transLp.untransform(cutVals, inds, rhs))
                  success |= cutGen.finalizeAndAddCut(transLp.getGlobaldom(),
                                                      inds, cutVals, rhs);

                // printf("cut is violated for k = %d\n", k);
                break;
              }
            }
            fLast = f;
          }
        }
      }

      lpAggregator.clear();

      if (!tryNegatedScale) break;
    }
  }
}
