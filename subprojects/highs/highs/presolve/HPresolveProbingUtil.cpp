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
#include "util/HighsIntegers.h"

namespace presolve {

HPresolve::Result HPresolve::normaliseCliqueRows(
    HighsPostsolveStack& postsolve_stack) {
  struct nonZero {
    HighsInt index;
    double value;
    int8_t complementation;
    HighsInt position;
  };

  auto transformRow = [&](HighsInt row, std::vector<nonZero>& nzs, double scale,
                          HighsCDouble& update) {
    nzs.clear();
    update = 0.0;

    for (HighsInt rowiter : rowpositions) {
      HighsInt col = Acol[rowiter];
      double val = scale * Avalue[rowiter];

      if (val < 0) {
        nzs.push_back(nonZero{col, -val, -1, rowiter});
        update += static_cast<HighsCDouble>(val) * model->col_upper_[col];
      } else {
        nzs.push_back(nonZero{col, val, 1, rowiter});
        update += static_cast<HighsCDouble>(val) * model->col_lower_[col];
      }
    }
  };

  // maximum dynamism
  const double maxDynamism = 1e5;

  std::vector<nonZero> nzs;
  std::vector<HighsInt> perm;
  std::vector<double> rowCoefsInt;

  for (HighsInt row = 0; row < model->num_row_; row++) {
    // skip deleted and non-all-integer rows
    if (rowDeleted[row] || rowsize[row] <= 1 ||
        rowsizeInteger[row] != rowsize[row])
      continue;

    // skip redundant rows
    if (isRedundant(row)) continue;

    // store row
    storeRow(row);

    // skip rows that are not all-binary
    bool allBinary = true;
    bool isSetPpc = true;
    rowCoefsInt.clear();
    HighsInt nComp = 0;
    for (const auto& nz : getStoredRow()) {
      allBinary = isBinary(nz.index());
      if (!allBinary) break;
      isSetPpc = isSetPpc && std::abs(nz.value()) == 1.0;
      rowCoefsInt.push_back(nz.value());
      if (nz.value() < 0) nComp++;
    }
    if (!allBinary) continue;

    // skip row if it is already setppc
    if (isSetPpc && model->row_upper_[row] == 1.0 - nComp) continue;

    // skip row if dynamism is too large
    if (computeDynamism(getStoredRow()) > maxDynamism) continue;

    // compute scaling factor that makes all coefficients integral
    double intScale = HighsIntegers::integralScale(
        rowCoefsInt, options->small_matrix_value, options->small_matrix_value);
    if (intScale == 0 || intScale > 1e3) continue;

    // check if lhs / rhs are redundant
    double lhs = -kHighsInf;
    double rhs = kHighsInf;
    if (impliedRowBounds.getSumLower(row) <
        model->row_lower_[row] - primal_feastol)
      lhs = model->row_lower_[row];
    if (impliedRowBounds.getSumUpper(row) >
        model->row_upper_[row] + primal_feastol)
      rhs = model->row_upper_[row];

    // transform row
    HighsInt direction = rhs < kHighsInf ? HighsInt{1} : HighsInt{-1};
    HighsCDouble update;
    transformRow(row, nzs, direction * intScale, update);

    // update lhs / rhs
    lhs *= direction * intScale;
    rhs *= direction * intScale;
    if (direction < 0) std::swap(lhs, rhs);
    if (lhs > -kHighsInf)
      lhs = static_cast<double>(ceil(lhs - update - primal_feastol));
    if (rhs < kHighsInf)
      rhs = static_cast<double>(floor(rhs - update + primal_feastol));

    // equation?
    bool equation = rhs == lhs;

    // after rounding, a ranged row may have become an equation; skip if still
    // ranged and not an equation
    if (lhs > -kHighsInf && rhs < kHighsInf && !equation) continue;

    // sort by descending coefficient value
    HighsInt numBin = static_cast<HighsInt>(nzs.size());
    perm.resize(numBin);
    std::iota(perm.begin(), perm.end(), 0);
    pdqsort(perm.begin(), perm.end(), [&](HighsInt a, HighsInt b) {
      return nzs[a].value > nzs[b].value;
    });

    // trivial fixings: variables with coefficient > rhs must be zero
    // in the complemented space
    HighsInt start = 0;
    for (HighsInt j = 0; j < numBin; ++j) {
      if (nzs[perm[j]].value <= rhs + primal_feastol) break;
      HighsInt col = nzs[perm[j]].index;
      if (nzs[perm[j]].complementation == -1)
        HPRESOLVE_CHECKED_CALL(fixColToUpper(postsolve_stack, col));
      else
        HPRESOLVE_CHECKED_CALL(fixColToLower(postsolve_stack, col));
      ++start;
    }

    // skip row if there are less than two remaining variables or two smallest
    // remaining coefficients do not form a clique
    if (numBin - start < 2 ||
        nzs[perm[numBin - 2]].value + nzs[perm[numBin - 1]].value <=
            rhs + primal_feastol)
      continue;

    // for equations, normalisation to x1 + ... + xn = 1 is only valid
    // for set partitioning rows where all coefficients equal rhs
    if (equation && (nzs[perm[start]].value != nzs[perm[numBin - 1]].value ||
                     nzs[perm[start]].value != rhs))
      continue;

    // re-store row since trivial fixings may have removed nonzeros
    storeRow(row);

    // scale so that we have a <= inequality or equation
    intScale *= direction;
    scaleStoredRow(row, intScale, true);

    // normalize remaining coefficients to ±1
    HighsInt numComp = 0;
    for (HighsInt j = start; j < numBin; ++j) {
      HighsInt pos = nzs[perm[j]].position;
      double delta =
          static_cast<double>(nzs[perm[j]].complementation) - Avalue[pos];
      if (delta != 0.0) addToMatrix(row, nzs[perm[j]].index, delta);
      if (nzs[perm[j]].complementation == -1) numComp++;
    }

    // update row bounds
    changeRowUpper(row, 1.0 - numComp);
    if (equation) changeRowLower(row, 1.0 - numComp);
  }

  return Result::kOk;
}

HPresolve::Result HPresolve::prepareProbing(
    HighsPostsolveStack& postsolve_stack, bool& firstCall) {
  HighsDomain& domain = mipsolver->mipdata_->getDomain();
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;

  shrinkProblem(postsolve_stack);

  // first call?
  firstCall = !mipsolver->mipdata_->cliquesExtracted;

  // todo: rows that become setppc after domain propagation (below) are
  // not normalised and thus cannot be deleted by clique merging
  if (firstCall) HPRESOLVE_CHECKED_CALL(normaliseCliqueRows(postsolve_stack));

  toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
        model->a_matrix_.start_);
  okFromCSC(model->a_matrix_.value_, model->a_matrix_.index_,
            model->a_matrix_.start_);

  cliquetable.setMaxEntries(numNonzeros());

  // first tighten all bounds if they have an implied bound that is tighter
  // than their column bound before probing this is not done for continuous
  // columns since it may allow stronger dual presolve and more aggregations
  double hugeBound = primal_feastol / kHighsTiny;
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (std::abs(implColLower[i]) <= hugeBound &&
        implColLower[i] > model->col_lower_[i])
      HPRESOLVE_CHECKED_CALL(changeColLower(i, implColLower[i]));

    if (std::abs(implColUpper[i]) <= hugeBound &&
        implColUpper[i] < model->col_upper_[i])
      HPRESOLVE_CHECKED_CALL(changeColUpper(i, implColUpper[i]));
  }

  // prepare for domain propagation
  mipsolver->mipdata_->setupDomainPropagation();

  domain.propagate();
  if (domain.infeasible()) return Result::kPrimalInfeasible;

  // extract cliques that are part of the formulation every time before probing
  // after the first call we only add cliques that directly correspond to set
  // packing constraints so that the clique merging step can extend/delete them
  if (firstCall) {
    mipsolver->mipdata_->cliquesExtracted = true;

    cliquetable.extractCliques(*mipsolver);
    if (domain.infeasible()) return Result::kPrimalInfeasible;

    // during presolve we keep the objective upper bound without the current
    // offset so we need to update it

    if (mipsolver->mipdata_->upper_limit != kHighsInf) {
      double tmpLimit = mipsolver->mipdata_->upper_limit;
      mipsolver->mipdata_->upper_limit = tmpLimit - model->offset_;
      cliquetable.extractObjCliques(*mipsolver);
      mipsolver->mipdata_->upper_limit = tmpLimit;
      if (domain.infeasible()) return Result::kPrimalInfeasible;
    }

    domain.propagate();
    if (domain.infeasible()) return Result::kPrimalInfeasible;
  }

  cliquetable.cleanupFixed(domain);
  if (domain.infeasible()) return Result::kPrimalInfeasible;

  return Result::kOk;
}

HPresolve::Result HPresolve::finaliseProbing(
    HighsPostsolveStack& postsolve_stack, bool firstCall,
    HighsInt& numVarsFixed, HighsInt& numBndsTightened,
    HighsInt& numVarsSubstituted, HighsInt& liftedNonZeros) {
  HighsDomain& domain = mipsolver->mipdata_->getDomain();
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;

  cliquetable.cleanupFixed(domain);

  if (!firstCall) cliquetable.extractCliques(*mipsolver, false);
  cliquetable.runCliqueMerging(domain);

  // apply changes from probing

  // first delete redundant clique inequalities
  for (HighsInt delrow : cliquetable.getDeletedRows())
    if (!rowDeleted[delrow]) removeRow(delrow);
  cliquetable.getDeletedRows().clear();

  // add nonzeros from clique lifting before removing fixed variables, since
  // this might lead to stronger constraint sides
  auto& extensionvars = cliquetable.getCliqueExtensions();
  liftedNonZeros += static_cast<HighsInt>(extensionvars.size());
  for (const auto& cliqueextension : extensionvars) {
    if (rowDeleted[cliqueextension.first]) {
      --liftedNonZeros;
      continue;
    }
    double val = 1.0;
    if (cliqueextension.second.val == 0) {
      addToRowLower(cliqueextension.first, HighsCDouble{-1});
      addToRowUpper(cliqueextension.first, HighsCDouble{-1});
      val = -1.0;
    }
    addToMatrix(cliqueextension.first, cliqueextension.second.col, val);
  }
  extensionvars.clear();

  // now remove fixed columns and tighten domains
  for (HighsInt i = 0; i != model->num_col_; ++i) {
    if (colDeleted[i]) continue;
    bool newLowerBnd = model->col_lower_[i] < domain.col_lower_[i];
    bool newUpperBnd = model->col_upper_[i] > domain.col_upper_[i];
    if (newLowerBnd)
      HPRESOLVE_CHECKED_CALL(changeColLower(i, domain.col_lower_[i]));
    if (newUpperBnd)
      HPRESOLVE_CHECKED_CALL(changeColUpper(i, domain.col_upper_[i]));
    if (domain.isFixed(i)) {
      numVarsFixed++;
      postsolve_stack.removedFixedCol(i, model->col_lower_[i], 0.0,
                                      HighsEmptySlice());
      HPRESOLVE_CHECKED_CALL(removeFixedCol(i));
    } else {
      if (newLowerBnd) numBndsTightened++;
      if (newUpperBnd) numBndsTightened++;
    }
    // Do not check limits here: rows have already been deleted by the
    // clique table above without postsolve entries, relying on all
    // column fixings being applied to justify their redundancy.
    // HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  // finally apply substitutions
  HPRESOLVE_CHECKED_CALL(
      applyConflictGraphSubstitutions(postsolve_stack, numVarsSubstituted));

  // aggregate variable bounds now that the clique table is populated
  aggregateVarBounds();

  return checkLimits(postsolve_stack);
}

HPresolve::Result HPresolve::applyConflictGraphSubstitutions(
    HighsPostsolveStack& postsolve_stack, HighsInt& numDelCol) {
  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;
  for (const auto& substitution : implications.substitutions) {
    if (colDeleted[substitution.substcol] || colDeleted[substitution.staycol])
      continue;

    ++numDelCol;

    postsolve_stack.doubletonEquation(
        -1, substitution.substcol, substitution.staycol, 1.0,
        -substitution.scale, substitution.offset,
        model->col_lower_[substitution.substcol],
        model->col_upper_[substitution.substcol], 0.0, false, false,
        HighsPostsolveStack::RowType::kEq, HighsEmptySlice());
    HPRESOLVE_CHECKED_CALL(substitute(substitution.substcol,
                                      substitution.staycol, substitution.offset,
                                      substitution.scale));
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  implications.substitutions.clear();

  for (HighsCliqueTable::Substitution subst : cliquetable.getSubstitutions()) {
    if (colDeleted[subst.substcol] || colDeleted[subst.replace.col]) continue;

    double scale;
    double offset;

    ++numDelCol;

    if (subst.replace.val == 0) {
      scale = -1.0;
      offset = 1.0;
    } else {
      scale = 1.0;
      offset = 0.0;
    }

    postsolve_stack.doubletonEquation(
        -1, subst.substcol, subst.replace.col, 1.0, -scale, offset,
        model->col_lower_[subst.substcol], model->col_upper_[subst.substcol],
        0.0, false, false, HighsPostsolveStack::RowType::kEq,
        HighsEmptySlice());
    HPRESOLVE_CHECKED_CALL(
        substitute(subst.substcol, subst.replace.col, offset, scale));
    HPRESOLVE_CHECKED_CALL(checkLimits(postsolve_stack));
  }

  cliquetable.getSubstitutions().clear();

  return Result::kOk;
}

}  // namespace presolve
