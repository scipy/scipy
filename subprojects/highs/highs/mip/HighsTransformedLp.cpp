/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsTransformedLp.h"

#include "mip/HighsMipSolverData.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"

HighsTransformedLp::HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                                       HighsImplications& implications,
                                       const HighsDomain& globaldom)
    : lprelaxation(lprelaxation), globaldom_(globaldom) {
  assert(lprelaxation.scaledOptimal(lprelaxation.getStatus()));
  const HighsMipSolver& mipsolver = implications.mipsolver;
  const HighsSolution& lpSolution = lprelaxation.getLpSolver().getSolution();

  HighsInt numTransformedCol = lprelaxation.numCols() + lprelaxation.numRows();

  boundDist.resize(numTransformedCol);
  simpleLbDist.resize(numTransformedCol);
  simpleUbDist.resize(numTransformedCol);
  lbDist.resize(numTransformedCol);
  ubDist.resize(numTransformedCol);
  bestVlb.resize(numTransformedCol,
                 std::make_pair(-1, HighsImplications::VarBound()));
  bestVub.resize(numTransformedCol,
                 std::make_pair(-1, HighsImplications::VarBound()));
  boundTypes.resize(numTransformedCol);
  vectorsum.setDimension(numTransformedCol);

  for (HighsInt col : mipsolver.mipdata_->continuous_cols) {
    if (!mipsolver.mipdata_->parallelLockActive())
      mipsolver.mipdata_->implications.cleanupVarbounds(col);

    if (globaldom_.infeasible()) return;

    if (globaldom_.isFixed(col)) continue;

    double bestub = globaldom_.col_upper_[col];
    simpleUbDist[col] = bestub - lpSolution.col_value[col];
    if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
      simpleUbDist[col] = 0.0;
    bestVub[col] = implications.getBestVub(col, lpSolution, bestub, globaldom_);

    double bestlb = globaldom_.col_lower_[col];
    simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
    if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
      simpleLbDist[col] = 0.0;
    bestVlb[col] = implications.getBestVlb(col, lpSolution, bestlb, globaldom_);

    lbDist[col] = lpSolution.col_value[col] - bestlb;
    if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
    ubDist[col] = bestub - lpSolution.col_value[col];
    if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;

    boundDist[col] = std::min(lbDist[col], ubDist[col]);
  }

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    double bestub = globaldom_.col_upper_[col];
    double bestlb = globaldom_.col_lower_[col];

    if (!mipsolver.mipdata_->parallelLockActive())
      mipsolver.mipdata_->implications.cleanupVarbounds(col);

    if (globaldom_.infeasible()) return;
    simpleUbDist[col] = bestub - lpSolution.col_value[col];
    if (simpleUbDist[col] <= mipsolver.mipdata_->feastol)
      simpleUbDist[col] = 0.0;

    simpleLbDist[col] = lpSolution.col_value[col] - bestlb;
    if (simpleLbDist[col] <= mipsolver.mipdata_->feastol)
      simpleLbDist[col] = 0.0;
    double simpleBndDist = std::min(simpleLbDist[col], simpleUbDist[col]);
    if (simpleBndDist > 0 &&
        std::fabs(HighsIntegers::nearestInteger(lpSolution.col_value[col]) -
                  lpSolution.col_value[col]) < mipsolver.mipdata_->feastol) {
      bestVub[col] = mipsolver.mipdata_->implications.getBestVub(
          col, lpSolution, bestub, globaldom_);
      bestVlb[col] = mipsolver.mipdata_->implications.getBestVlb(
          col, lpSolution, bestlb, globaldom_);

      lbDist[col] = lpSolution.col_value[col] - bestlb;
      if (lbDist[col] <= mipsolver.mipdata_->feastol) lbDist[col] = 0.0;
      ubDist[col] = bestub - lpSolution.col_value[col];
      if (ubDist[col] <= mipsolver.mipdata_->feastol) ubDist[col] = 0.0;
      boundDist[col] = std::min(lbDist[col], ubDist[col]);
      if (boundDist[col] > simpleBndDist + mipsolver.mipdata_->feastol) {
        lbDist[col] = simpleLbDist[col];
        ubDist[col] = simpleUbDist[col];
        boundDist[col] = simpleBndDist;
        bestVub[col].first = -1;
        bestVlb[col].first = -1;
      }
    } else {
      lbDist[col] = simpleLbDist[col];
      ubDist[col] = simpleUbDist[col];
      boundDist[col] = simpleBndDist;
    }
  }

  // setup information of slackVariables
  HighsInt numLpRow = lprelaxation.numRows();
  HighsInt indexOffset = mipsolver.numCol();
  for (HighsInt row = 0; row != numLpRow; ++row) {
    HighsInt slackIndex = indexOffset + row;
    double bestub = lprelaxation.slackUpper(row, globaldom_);
    double bestlb = lprelaxation.slackLower(row, globaldom_);

    if (bestlb == bestub) continue;

    lbDist[slackIndex] = lpSolution.row_value[row] - bestlb;
    if (lbDist[slackIndex] <= mipsolver.mipdata_->feastol)
      lbDist[slackIndex] = 0.0;
    simpleLbDist[slackIndex] = lbDist[slackIndex];

    ubDist[slackIndex] = bestub - lpSolution.row_value[row];
    if (ubDist[slackIndex] <= mipsolver.mipdata_->feastol)
      ubDist[slackIndex] = 0.0;
    simpleUbDist[slackIndex] = ubDist[slackIndex];

    boundDist[slackIndex] = std::min(lbDist[slackIndex], ubDist[slackIndex]);
  }
}

bool HighsTransformedLp::transform(std::vector<double>& vals,
                                   std::vector<double>& upper,
                                   std::vector<double>& solval,
                                   std::vector<HighsInt>& inds, double& rhs,
                                   bool& integersPositive, bool preferVbds,
                                   bool enforceSameBds) {
  // vector sum should be empty
  assert(vectorsum.getNonzeros().empty());

  HighsCDouble tmpRhs = rhs;

  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = lprelaxation.numCols();

  HighsInt numNz = inds.size();

  auto getLb = [&](HighsInt col) {
    return (col < slackOffset
                ? globaldom_.col_lower_[col]
                : lprelaxation.slackLower(col - slackOffset, globaldom_));
  };

  auto getUb = [&](HighsInt col) {
    return (col < slackOffset
                ? globaldom_.col_upper_[col]
                : lprelaxation.slackUpper(col - slackOffset, globaldom_));
  };

  auto remove = [&](HighsInt position) {
    numNz--;
    inds[position] = inds[numNz];
    vals[position] = vals[numNz];
    inds[numNz] = 0;
    vals[numNz] = 0;
  };

  HighsInt i = 0;
  while (i < numNz) {
    HighsInt col = inds[i];

    double lb = getLb(col);
    double ub = getUb(col);

    if (ub - lb < mip.options_mip_->small_matrix_value) {
      tmpRhs -= std::min(lb, ub) * vals[i];
      remove(i);
      continue;
    }

    if (lb == -kHighsInf && ub == kHighsInf) {
      vectorsum.clear();
      return false;
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable upper bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable upper
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable upper bound.
    if (bestVub[col].first != -1 &&
        bestVub[col].second.maxValue() > ub + mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVub(col, bestVub[col].first,
                                            bestVub[col].second, ub, redundant,
                                            infeasible, false);
      if (redundant) {
        bestVub[col].first = -1;
        ubDist[col] = simpleUbDist[col];
        boundDist[col] = std::min(lbDist[col], ubDist[col]);
      }
    }

    // the code below uses the difference between the column upper and lower
    // bounds as the upper bound for the slack from the variable lower bound
    // constraint (upper[j] = ub - lb) and thus assumes that the variable lower
    // bound constraints are tight. this assumption may not be satisfied when
    // new bound changes were derived during cut generation and, therefore, we
    // tighten the best variable lower bound.
    if (bestVlb[col].first != -1 &&
        bestVlb[col].second.minValue() < lb - mip.mipdata_->feastol) {
      bool redundant = false;
      bool infeasible = false;
      mip.mipdata_->implications.cleanupVlb(col, bestVlb[col].first,
                                            bestVlb[col].second, lb, redundant,
                                            infeasible, false);
      if (redundant) {
        bestVlb[col].first = -1;
        lbDist[col] = simpleLbDist[col];
        boundDist[col] = std::min(lbDist[col], ubDist[col]);
      }
    }

    // store the old bound type so that we can restore it if the continuous
    // column is relaxed out anyways. This allows to correctly transform and
    // then untransform multiple base rows which is useful to compute cuts based
    // on several transformed base rows. It could otherwise lead to bugs if a
    // column is first transformed with a simple bound and not relaxed but for
    // another base row is transformed and relaxed with a variable bound. Should
    // the non-relaxed column now be untransformed we would wrongly use the
    // variable bound even though this is not the correct way to untransform the
    // column.
    BoundType oldBoundType = boundTypes[col];

    if (enforceSameBds && oldBoundType != BoundType::kUnused) {
      // Don't allow different bound types to be used in stacked
      // transformations unless they're going to be relaxed out
      // This ensures that a single call to untransform remains valid
      if (lprelaxation.isColIntegral(col) &&
          (ub - lb <= 1.5 || boundDist[col] != 0.0 || simpleLbDist[col] == 0 ||
           simpleUbDist[col] == 0)) {
        i++;
        continue;
      }
      if (oldBoundType == BoundType::kSimpleLb ||
          oldBoundType == BoundType::kSimpleUb) {
        if (lbDist[col] < ubDist[col] && vals[i] > 0) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else if (ubDist[col] < lbDist[col] && vals[i] < 0) {
          boundTypes[col] = BoundType::kSimpleUb;
        }
      }
    } else if (lprelaxation.isColIntegral(col)) {
      if (ub - lb <= 1.5 || boundDist[col] != 0.0 || simpleLbDist[col] == 0 ||
          simpleUbDist[col] == 0) {
        // since we skip the handling of variable bound constraints for all
        // binary and some general-integer variables here, the bound type used
        // should be a simple lower or upper bound
        if (simpleLbDist[col] < simpleUbDist[col] - mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else if (simpleUbDist[col] <
                   simpleLbDist[col] - mip.mipdata_->feastol) {
          boundTypes[col] = BoundType::kSimpleUb;
        } else if (vals[i] > 0) {
          boundTypes[col] = BoundType::kSimpleLb;
        } else {
          boundTypes[col] = BoundType::kSimpleUb;
        }
        i++;
        continue;
      }
      if (bestVlb[col].first == -1 ||
          ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        assert(bestVub[col].first != -1);
        boundTypes[col] = BoundType::kVariableUb;
      } else if (bestVub[col].first == -1 ||
                 lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        assert(bestVlb[col].first != -1);
        boundTypes[col] = BoundType::kVariableLb;
      } else if (vals[i] > 0) {
        assert(bestVub[col].first != -1);
        boundTypes[col] = BoundType::kVariableUb;
      } else {
        assert(bestVlb[col].first != -1);
        boundTypes[col] = BoundType::kVariableLb;
      }
    } else {
      if (lbDist[col] < ubDist[col] - mip.mipdata_->feastol) {
        if (bestVlb[col].first == -1)
          boundTypes[col] = BoundType::kSimpleLb;
        else if (preferVbds || vals[i] > 0 ||
                 simpleLbDist[col] > lbDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else if (ubDist[col] < lbDist[col] - mip.mipdata_->feastol) {
        if (bestVub[col].first == -1)
          boundTypes[col] = BoundType::kSimpleUb;
        else if (preferVbds || vals[i] < 0 ||
                 simpleUbDist[col] > ubDist[col] + mip.mipdata_->feastol)
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      } else if (vals[i] > 0) {
        if (bestVlb[col].first != -1)
          boundTypes[col] = BoundType::kVariableLb;
        else if (preferVbds && bestVub[col].first != -1)
          boundTypes[col] = BoundType::kVariableUb;
        else
          boundTypes[col] = BoundType::kSimpleLb;
      } else {
        if (bestVub[col].first != -1)
          boundTypes[col] = BoundType::kVariableUb;
        else if (preferVbds && bestVlb[col].first != -1)
          boundTypes[col] = BoundType::kVariableLb;
        else
          boundTypes[col] = BoundType::kSimpleUb;
      }
    }

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb:
        if (vals[i] > 0) {
          // relax away using lower bound
          tmpRhs -= lb * vals[i];
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kSimpleUb:
        if (vals[i] < 0) {
          // relax away using upper bound
          tmpRhs -= ub * vals[i];
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kVariableLb:
        tmpRhs -= bestVlb[col].second.constant * vals[i];
        vectorsum.add(bestVlb[col].first, vals[i] * bestVlb[col].second.coef);
        if (vals[i] > 0) {
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kVariableUb:
        tmpRhs -= bestVub[col].second.constant * vals[i];
        vectorsum.add(bestVub[col].first, vals[i] * bestVub[col].second.coef);
        vals[i] = -vals[i];
        if (vals[i] > 0) {
          boundTypes[col] = oldBoundType;
          remove(i);
          continue;
        }
        break;
      case BoundType::kUnused:
        assert(false);
        vectorsum.clear();
        return false;
    }
    // move to next element
    i++;
  }

  if (!vectorsum.getNonzeros().empty()) {
    for (HighsInt i = 0; i != numNz; ++i) {
      if (vals[i] != 0.0) vectorsum.add(inds[i], vals[i]);
    }

    double maxError = 0.0;
    auto IsZero = [&](HighsInt col, double val) {
      if (termIsZero(val, getLb(col), getUb(col), mip.mipdata_->epsilon))
        return true;
      return false;
    };

    vectorsum.cleanup(IsZero);
    if (maxError > mip.mipdata_->feastol) {
      vectorsum.clear();
      return false;
    }

    inds = vectorsum.getNonzeros();
    numNz = inds.size();

    vals.resize(numNz);
    for (HighsInt j = 0; j != numNz; ++j) vals[j] = vectorsum.getValue(inds[j]);

    vectorsum.clear();
  } else {
    vals.resize(numNz);
    inds.resize(numNz);
  }

  for (HighsInt j = 0; j != numNz; ++j) {
    HighsInt col = inds[j];

    if (enforceSameBds && boundTypes[col] != BoundType::kUnused) continue;

    // get bounds
    double lb = getLb(col);
    double ub = getUb(col);

    // variable should not be free
    assert(lb != -kHighsInf || ub != kHighsInf);

    // set bound type for previously unprocessed integer-constrained variables
    if (!lprelaxation.isColIntegral(col) &&
        boundTypes[col] != BoundType::kUnused)
      continue;

    // do not overwrite bound type for integral slacks from vlb / vub
    // constraints
    if (boundTypes[col] == BoundType::kVariableLb ||
        boundTypes[col] == BoundType::kVariableUb)
      continue;

    // complement integers to make coefficients positive if both bounds are
    // finite; otherwise, complement integers with closest bound.
    // take into account 'integersPositive' provided by caller.
    if (integersPositive) {
      if ((lb != -kHighsInf && vals[j] > 0) || ub == kHighsInf)
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    } else {
      if (lbDist[col] < ubDist[col])
        boundTypes[col] = BoundType::kSimpleLb;
      else
        boundTypes[col] = BoundType::kSimpleUb;
    }
  }

  upper.resize(numNz);
  solval.resize(numNz);

  for (HighsInt j = 0; j != numNz; ++j) {
    HighsInt col = inds[j];

    double lb = getLb(col);
    double ub = getUb(col);

    upper[j] = ub - lb;

    switch (boundTypes[col]) {
      case BoundType::kSimpleLb: {
        // shift (lower bound)
        assert(lb != -kHighsInf);
        tmpRhs -= lb * vals[j];
        solval[j] = lbDist[col];
        break;
      }
      case BoundType::kSimpleUb: {
        // complement (upper bound)
        assert(ub != kHighsInf);
        tmpRhs -= ub * vals[j];
        vals[j] = -vals[j];
        solval[j] = ubDist[col];
        break;
      }
      case BoundType::kVariableLb: {
        solval[j] = lbDist[col];
        break;
      }
      case BoundType::kVariableUb: {
        solval[j] = ubDist[col];
        break;
      }
      case BoundType::kUnused: {
        assert(false);
        return false;
      }
    }

    // check if all integer-constrained variables have positive coefficients
    if (lprelaxation.isColIntegral(col))
      integersPositive = integersPositive && vals[j] > 0;
  }

  rhs = double(tmpRhs);

  if (numNz == 0 && rhs >= -mip.mipdata_->feastol) return false;

  return true;
}

bool HighsTransformedLp::untransform(std::vector<double>& vals,
                                     std::vector<HighsInt>& inds, double& rhs,
                                     bool integral) {
  HighsCDouble tmpRhs = rhs;
  const HighsMipSolver& mip = lprelaxation.getMipSolver();
  const HighsInt slackOffset = mip.numCol();

  HighsInt numNz = inds.size();

  for (HighsInt i = 0; i != numNz; ++i) {
    if (vals[i] == 0.0) continue;
    HighsInt col = inds[i];

    switch (boundTypes[col]) {
      case BoundType::kVariableLb: {
        tmpRhs += bestVlb[col].second.constant * vals[i];
        vectorsum.add(bestVlb[col].first, -vals[i] * bestVlb[col].second.coef);
        vectorsum.add(col, vals[i]);
        break;
      }
      case BoundType::kVariableUb: {
        tmpRhs -= bestVub[col].second.constant * vals[i];
        vectorsum.add(bestVub[col].first, vals[i] * bestVub[col].second.coef);
        vectorsum.add(col, -vals[i]);
        break;
      }
      case BoundType::kSimpleLb: {
        if (col < slackOffset) {
          tmpRhs += vals[i] * globaldom_.col_lower_[col];
          vectorsum.add(col, vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs += vals[i] * lprelaxation.slackLower(row, globaldom_);

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
        break;
      }
      case BoundType::kSimpleUb: {
        if (col < slackOffset) {
          tmpRhs -= vals[i] * globaldom_.col_upper_[col];
          vectorsum.add(col, -vals[i]);
        } else {
          HighsInt row = col - slackOffset;
          tmpRhs -= vals[i] * lprelaxation.slackUpper(row, globaldom_);
          vals[i] = -vals[i];

          HighsInt rowlen;
          const HighsInt* rowinds;
          const double* rowvals;
          lprelaxation.getRow(row, rowlen, rowinds, rowvals);

          for (HighsInt j = 0; j != rowlen; ++j)
            vectorsum.add(rowinds[j], vals[i] * rowvals[j]);
        }
        break;
      }
      case BoundType::kUnused: {
        assert(false);
        vectorsum.clear();
        return false;
      }
    }
  }

  if (integral) {
    // if the cut is integral, we just round all coefficient values and the
    // right hand side to the nearest integral value, as small deviation
    // only come from numerical errors during resubstitution of slack variables

    auto IsZero = [&](HighsInt col, double val) {
      assert(col < mip.numCol());
      return fabs(val) < 0.5;
    };

    vectorsum.cleanup(IsZero);
    rhs = std::round(double(tmpRhs));
  } else {
    bool abort = false;
    auto IsZero = [&](HighsInt col, double val) {
      assert(col < mip.numCol());
      if (termIsZero(val, globaldom_.col_lower_[col],
                     globaldom_.col_upper_[col], mip.mipdata_->epsilon))
        return true;
      double absval = std::abs(val);

      if (absval <= mip.mipdata_->feastol) {
        if (val > 0) {
          if (globaldom_.col_lower_[col] == -kHighsInf)
            abort = true;
          else
            tmpRhs -= val * globaldom_.col_lower_[col];
        } else {
          if (globaldom_.col_upper_[col] == kHighsInf)
            abort = true;
          else
            tmpRhs -= val * globaldom_.col_upper_[col];
        }
        return true;
      }
      return false;
    };

    vectorsum.cleanup(IsZero);
    if (abort) {
      vectorsum.clear();
      return false;
    }
    rhs = double(tmpRhs);
  }

  inds = vectorsum.getNonzeros();
  numNz = inds.size();
  vals.resize(numNz);

  if (integral)
    for (HighsInt i = 0; i != numNz; ++i)
      vals[i] = std::round(vectorsum.getValue(inds[i]));
  else
    for (HighsInt i = 0; i != numNz; ++i) vals[i] = vectorsum.getValue(inds[i]);
  vectorsum.clear();

  return true;
}
