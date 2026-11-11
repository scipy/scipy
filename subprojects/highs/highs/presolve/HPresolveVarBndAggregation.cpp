/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolverData.h"
#include "presolve/HPresolve.h"

namespace presolve {

void HPresolve::extractVarBounds(HighsInt row) {
  // extract variable bound constraints from the row

  // return if row is empty or a singleton or contains no integer variables
  if (mipsolver == nullptr || rowsize[row] <= 1 || rowsizeInteger[row] == 0)
    return;

  // check if variable bounds can be derived from row
  HighsInt numInfSumLower = impliedRowBounds.getNumInfSumLower(row);
  HighsInt numInfSumUpper = impliedRowBounds.getNumInfSumUpper(row);
  bool useLhs = model->row_lower_[row] != -kHighsInf && numInfSumUpper <= 1;
  bool useRhs = model->row_upper_[row] != kHighsInf && numInfSumLower <= 1;
  if (!useLhs && !useRhs) return;

  // check if row contains a single binary variable
  HighsInt binCol = -1;
  double binCoef = 0.0;
  for (const auto& nonzero : getRowVector(row)) {
    // skip fixed variables
    if (model->col_lower_[nonzero.index()] ==
        model->col_upper_[nonzero.index()])
      continue;

    // find binary variable
    if (model->integrality_[nonzero.index()] == HighsVarType::kInteger &&
        model->col_lower_[nonzero.index()] == 0.0 &&
        model->col_upper_[nonzero.index()] == 1.0) {
      // return if there is more than one binary variable
      if (binCol != -1) return;
      binCol = nonzero.index();
      binCoef = nonzero.value();
    }
  }
  // return if there is no binary variable
  if (binCol == -1) return;

  for (const auto& nonzero : getRowVector(row)) {
    // skip fixed variables
    if (model->col_lower_[nonzero.index()] ==
        model->col_upper_[nonzero.index()])
      continue;

    // skip binary variable
    if (nonzero.index() == binCol) continue;

    // compute VLB constant
    double vlbConstant = -kHighsInf;
    if (useLhs) {
      double residual = impliedRowBounds.getResidualSumUpper(
          row, nonzero.index(), nonzero.value(), binCol, binCoef, 0.0);
      if (residual != kHighsInf) {
        vlbConstant = static_cast<double>(
            (static_cast<HighsCDouble>(model->row_lower_[row]) - residual) /
            std::abs(nonzero.value()));
        useLhs = numInfSumUpper == 0;
      }
    }

    // compute VUB constant
    double vubConstant = kHighsInf;
    if (useRhs) {
      double residual = impliedRowBounds.getResidualSumLower(
          row, nonzero.index(), nonzero.value(), binCol, binCoef, 0.0);
      if (residual != -kHighsInf) {
        vubConstant = static_cast<double>(
            (static_cast<HighsCDouble>(model->row_upper_[row]) - residual) /
            std::abs(nonzero.value()));
        useRhs = numInfSumLower == 0;
      }
    }

    // switch sign if continuous variable has a negative coefficient
    if (nonzero.value() < 0) {
      vlbConstant *= -1;
      vubConstant *= -1;
      std::swap(vlbConstant, vubConstant);
    }

    // compute coefficient for binary variable
    double vbCoef = -binCoef / nonzero.value();

    // only record origin for size-2 one-sided rows that can be safely rewritten
    HighsInt origin = rowsize[row] == 2 && !isRanged(row) ? row : -1;

    // add VLB
    if (vlbConstant != -kHighsInf)
      mipsolver->mipdata_->implications.addVLB(
          nonzero.index(), binCol, vbCoef, vlbConstant,
          model->col_lower_[nonzero.index()],
          model->integrality_[nonzero.index()] != HighsVarType::kContinuous,
          origin);

    // add VUB
    if (vubConstant != kHighsInf)
      mipsolver->mipdata_->implications.addVUB(
          nonzero.index(), binCol, vbCoef, vubConstant,
          model->col_upper_[nonzero.index()],
          model->integrality_[nonzero.index()] != HighsVarType::kContinuous,
          origin);

    // stop if no additional variable bounds can be found
    if (!useLhs && !useRhs) break;
  }
}

void HPresolve::aggregateVarBounds() {
  if (mipsolver == nullptr) return;

  HighsImplications& implications = mipsolver->mipdata_->implications;
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;

  struct colImpliedBounds {
    HighsImplications::VarBound originalBound;
    HighsImplications::VarBound standardBound;
  };

  // data structures
  HighsHashTree<HighsInt, colImpliedBounds> vlbs;
  HighsHashTree<HighsInt, colImpliedBounds> vubs;
  std::vector<HighsCliqueTable::CliqueVar> vlbsClique;
  std::vector<HighsCliqueTable::CliqueVar> vubsClique;
  std::vector<std::vector<HighsCliqueTable::CliqueVar>> vlbsCover;
  std::vector<std::vector<HighsCliqueTable::CliqueVar>> vubsCover;
  std::vector<HighsBool> rowConsumed(model->num_row_, false);

  // initialise counters
  HighsInt numRowsRemoved = 0;
  HighsInt numRowsModified = 0;
  HighsInt numVarsLifted = 0;

  for (HighsInt col = 0; col != model->num_col_; ++col) {
    // skip deleted and binary columns (aggregate VLBs/VUBs on non-binary
    // columns)
    if (colDeleted[col] || isBinary(col)) continue;

    // clear vectors
    vlbs.clear();
    vubs.clear();
    vlbsClique.clear();
    vubsClique.clear();

    // get lower bound and upper bound
    double lb = model->col_lower_[col];
    double ub = model->col_upper_[col];

    // the minimum (maximum) of every VLB (VUB) is a valid lower (upper) bound,
    // so tighten lb (ub) such that the standardized forms below imply the VLBs
    // (VUBs) exactly; the column bounds in the model are not changed
    implications.getVlbs(col).for_each(
        [&](HighsInt binaryCol, const HighsImplications::VarBound& vlb) {
          if (!colDeleted[binaryCol]) lb = std::max(lb, vlb.minValue());
        });
    implications.getVubs(col).for_each(
        [&](HighsInt binaryCol, const HighsImplications::VarBound& vub) {
          if (!colDeleted[binaryCol]) ub = std::min(ub, vub.maxValue());
        });

    if (model->integrality_[col] != HighsVarType::kContinuous) {
      lb = std::ceil(lb - mipsolver->mipdata_->feastol);
      ub = std::floor(ub + mipsolver->mipdata_->feastol);
    }

    // nothing to aggregate if the column has no VLBs and VUBs
    if (lb <= -kHighsInf && ub >= kHighsInf) continue;

    // nothing to aggregate if the tightened bounds fix the column or cross; the
    // latter indicates infeasibility, which is left to other reductions
    if (lb >= ub) continue;

    // compute range
    double range = kHighsInf;
    if (lb > -kHighsInf && ub < kHighsInf)
      range = static_cast<double>(static_cast<HighsCDouble>(ub) - lb);

    // number of VLBs (VUBs) whose origin row could be rewritten
    HighsInt numVlbOrigins = 0;
    HighsInt numVubOrigins = 0;

    // collect VLBs: y >= coef * x + constant
    // standardized form: y >= a * x + lb, where a > 0
    if (lb > -kHighsInf) {
      implications.getVlbs(col).for_each(
          [&](HighsInt binaryCol, const HighsImplications::VarBound& vlb) {
            // skip deleted cols
            if (colDeleted[binaryCol]) return;
            // skip if the VLB is dominated by the global lower bound
            if (implications.redundantVlb(vlb, lb)) return;
            // tighten so that minValue() >= lb
            HighsImplications::VarBound v = vlb;
            implications.tightenVlb(v, lb);
            // standardize: a = maxValue() - lb
            HighsCDouble newCoef = v.constant - static_cast<HighsCDouble>(lb);
            if (v.coef > 0) newCoef += v.coef;
            // skip if the standardized coefficient exceeds the variable range
            if (newCoef > range + mipsolver->mipdata_->feastol) return;
            vlbs.insert(binaryCol,
                        colImpliedBounds{v, HighsImplications::VarBound{
                                                static_cast<double>(newCoef),
                                                lb, v.origin}});
            if (v.origin >= 0) numVlbOrigins++;
          });
    }
    // collect VUBs: y <= coef * x + constant
    // standardized form: y <= ub - a * x, where a > 0
    if (ub < kHighsInf) {
      implications.getVubs(col).for_each(
          [&](HighsInt binaryCol, const HighsImplications::VarBound& vub) {
            // skip deleted cols
            if (colDeleted[binaryCol]) return;
            // skip if the VUB is dominated by the global upper bound
            if (implications.redundantVub(vub, ub)) return;
            // tighten so that maxValue() <= ub
            HighsImplications::VarBound v = vub;
            implications.tightenVub(v, ub);
            // standardize: a = ub - minValue()
            HighsCDouble newCoef = static_cast<HighsCDouble>(ub) - v.constant;
            if (v.coef < 0) newCoef -= v.coef;
            // skip if the standardized coefficient exceeds the variable range
            if (newCoef > range + mipsolver->mipdata_->feastol) return;
            vubs.insert(binaryCol,
                        colImpliedBounds{v, HighsImplications::VarBound{
                                                static_cast<double>(newCoef),
                                                ub, v.origin}});
            if (v.origin >= 0) numVubOrigins++;
          });
    }

    // only bounds with an origin row can be merged, so there is nothing to do
    // if no VLB and no VUB has one
    if (numVlbOrigins == 0 && numVubOrigins == 0) continue;

    // set up cliques for the sides with an origin row
    if (numVlbOrigins > 0)
      vlbs.for_each([&](HighsInt binaryCol, colImpliedBounds& bounds) {
        HighsInt val = bounds.originalBound.coef > 0 ? 1 : 0;
        vlbsClique.emplace_back(binaryCol, val);
      });
    if (numVubOrigins > 0)
      vubs.for_each([&](HighsInt binaryCol, colImpliedBounds& bounds) {
        HighsInt val = bounds.originalBound.coef < 0 ? 1 : 0;
        vubsClique.emplace_back(binaryCol, val);
      });

    // clique cover (clears the covers, which stay empty without cliques)
    cliquetable.cliqueCover(vlbsClique, vlbsCover);
    cliquetable.cliqueCover(vubsClique, vubsCover);

    auto mergeCliques =
        [&](std::vector<std::vector<HighsCliqueTable::CliqueVar>>& cover,
            HighsHashTree<HighsInt, colImpliedBounds>& boundsMap,
            double baseBound, HighsInt direction) {
          for (const auto& clique : cover) {
            if (clique.size() < 2) continue;

            // find an unconsumed row to reuse and remove the rest
            HighsInt row = -1;
            for (const auto& var : clique) {
              const auto* bounds = boundsMap.find(var.col);
              HighsInt currentrow = bounds->originalBound.origin;
              if (currentrow < 0) continue;
              if (!rowConsumed[currentrow]) {
                rowConsumed[currentrow] = true;
                if (row == -1)
                  row = currentrow;
                else {
                  removeRow(currentrow);
                  numRowsRemoved++;
                }
              }
            }
            if (row == -1) continue;

            // origin is invalidated when row bounds change, so the
            // row must still be one-sided
            assert(!isRanged(row));

            // rewrite the reused row with the aggregated constraint
            unlinkRow(row);
            addToMatrix(row, col, 1.0);
            numRowsModified++;

            HighsCDouble rowBound = baseBound;
            for (const auto& var : clique) {
              const auto* bounds = boundsMap.find(var.col);
              double a = bounds->standardBound.coef;
              numVarsLifted++;
              if (var.val == 1) {
                addToMatrix(row, var.col, -direction * a);
              } else {
                addToMatrix(row, var.col, direction * a);
                rowBound += direction * a;
              }
            }

            if (direction > 0) {
              changeRowLower(row, static_cast<double>(rowBound));
              changeRowUpper(row, kHighsInf);
            } else {
              changeRowLower(row, -kHighsInf);
              changeRowUpper(row, static_cast<double>(rowBound));
            }
          }
        };

    if (lb > -kHighsInf) mergeCliques(vlbsCover, vlbs, lb, HighsInt{1});
    if (ub < kHighsInf) mergeCliques(vubsCover, vubs, ub, HighsInt{-1});
  }

  if (numRowsRemoved > 0 || numRowsModified > 0)
    highsLogDev(options->log_options, HighsLogType::kInfo,
                "Implied variable bound aggregation: %" HIGHSINT_FORMAT
                " rows removed, %" HIGHSINT_FORMAT
                " rows modified, %" HIGHSINT_FORMAT " vars lifted\n",
                numRowsRemoved, numRowsModified, numVarsLifted);
}

}  // namespace presolve
