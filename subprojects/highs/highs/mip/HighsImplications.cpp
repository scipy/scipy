/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsImplications.h"

#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"

bool HighsImplications::computeImplications(HighsInt col, bool val) {
  HighsDomain& globaldomain = mipsolver.mipdata_->getDomain();
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  globaldomain.propagate();
  if (globaldomain.infeasible() || globaldomain.isFixed(col)) return true;

  // record redundant rows for lifting
  assert(globaldomain.getRedundantRows().size() == 0);
  if (storeLiftingOpportunity != nullptr)
    globaldomain.setRecordRedundantRows(true);

  const auto& domchgstack = globaldomain.getDomainChangeStack();
  const auto& domchgreason = globaldomain.getDomainChangeReason();
  size_t changedend = globaldomain.getChangedCols().size();

  const bool dualFixProbingActive = globaldomain.getDualFixProbingActive();
  if (dualFixProbingActive) {
    globaldomain.getDualFixProbingPropagation().beginProbing();
  }

  HighsInt stackimplicstart = domchgstack.size() + 1;
  HighsInt numImplications = -stackimplicstart;
  if (val)
    globaldomain.changeBound(HighsBoundType::kLower, col, 1);
  else
    globaldomain.changeBound(HighsBoundType::kUpper, col, 0);

  auto storeLiftingOpportunities = [&](HighsInt col, bool val) {
    // use callback to store lifting opportunities
    if (storeLiftingOpportunity != nullptr) {
      for (const auto& elm : globaldomain.getRedundantRows())
        storeLiftingOpportunity(
            elm.key(), col, val ? 1 : 0,
            (val ? -1 : 1) * globaldomain.getRedundantRowValue(elm.key()));
      globaldomain.clearRedundantRows();
      globaldomain.setRecordRedundantRows(false);
    }
  };

  if (dualFixProbingActive && storeLiftingOpportunity != nullptr) {
    globaldomain.getDualFixProbingPropagation().storeLiftingOpportunity =
        [&]() { storeLiftingOpportunities(col, val); };
  }

  auto doBacktrack = [&](size_t changedend) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);
  };

  auto isInfeasible = [&](HighsInt col, bool val) {
    if (!globaldomain.infeasible()) return false;
    if (dualFixProbingActive) {
      globaldomain.getDualFixProbingPropagation().endProbing();
    }
    storeLiftingOpportunities(col, val);
    doBacktrack(changedend);
    cliquetable.vertexInfeasible(globaldomain, col, val);
    return true;
  };

  if (isInfeasible(col, val)) return true;

  globaldomain.propagate();
  if (dualFixProbingActive) {
    globaldomain.getDualFixProbingPropagation().endProbing();
  }

  if (isInfeasible(col, val)) return true;

  HighsInt stackimplicend = domchgstack.size();
  numImplications += stackimplicend;
  mipsolver.mipdata_->getPseudoCost().addInferenceObservation(
      col, numImplications, val);

  std::vector<HighsDomainChange> origStackCopy;
  if (dualFixProbingActive) {
    origStackCopy.assign(domchgstack.begin() + stackimplicstart,
                         domchgstack.begin() + stackimplicend);
  }

  std::vector<HighsDomainChange>& implics =
      val ? implicationsUp : implicationsDown;
  implics.clear();
  implics.reserve(numImplications);

  HighsInt numEntries = mipsolver.mipdata_->cliquetable.getNumEntries();
  HighsInt maxEntries = 100000 + mipsolver.numNonzero();

  HighsInt unsafeStackStart = dualFixProbingActive
                                  ? globaldomain.getDualFixProbingPropagation()
                                        .getZeroCostFixingPosition()
                                  : stackimplicend;
  HighsInt safeImplicsEnd = 0;

  for (HighsInt i = stackimplicstart; i < stackimplicend; ++i) {
    if (domchgreason[i].type == HighsDomain::Reason::kCliqueTable &&
        ((domchgreason[i].index >> 1) == col || numEntries >= maxEntries))
      continue;

    implics.push_back(domchgstack[i]);
    if (i < unsafeStackStart) safeImplicsEnd++;
  }

  // inform caller about lifting opportunities
  storeLiftingOpportunities(col, val);

  // backtrack
  doBacktrack(changedend);

  for (const HighsDomainChange& implic : origStackCopy) {
    if (globaldomain.isBinary(implic.column)) {
      recordTentativeCliques(val, implic);
    }
  }

  if (safeImplicsEnd < static_cast<HighsInt>(implics.size())) {
    auto binstart =
        std::partition(implics.begin() + safeImplicsEnd, implics.end(),
                       [&](const HighsDomainChange& a) {
                         return !globaldomain.isBinary(a.column);
                       });
    implics.erase(binstart, implics.end());
  }

  // add the implications of binary variables to the clique table
  auto binstart =
      std::partition(implics.begin(), implics.begin() + safeImplicsEnd,
                     [&](const HighsDomainChange& a) {
                       return !globaldomain.isBinary(a.column);
                     });

  std::array<HighsCliqueTable::CliqueVar, 2> clique;
  clique[0] = HighsCliqueTable::CliqueVar(col, val);

  for (auto i = binstart; i != implics.begin() + safeImplicsEnd; ++i) {
    if (i->boundtype == HighsBoundType::kLower)
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 0);
    else
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 1);

    cliquetable.addClique(mipsolver, clique.data(), 2);
    if (globaldomain.infeasible() || globaldomain.isFixed(col)) return true;
  }

  HighsInt numErasedBinaries =
      static_cast<HighsInt>(implics.begin() + safeImplicsEnd - binstart);
  implics.erase(binstart, implics.begin() + safeImplicsEnd);
  safeImplicsEnd -= numErasedBinaries;

  // store variable bounds derived from implications
  for (auto i = implics.begin(); i != implics.begin() + safeImplicsEnd; ++i) {
    if (i->boundtype == HighsBoundType::kLower) {
      if (val == 1) {
        if (globaldomain.col_lower_[i->column] != -kHighsInf)
          addVLB(i->column, col,
                 i->boundval - globaldomain.col_lower_[i->column],
                 globaldomain.col_lower_[i->column]);
      } else
        addVLB(i->column,
               col,  // in case the lower bound is infinite the varbound can
                     // still be tightened as soon as a finite upper bound is
                     // known because the offset is finite
               globaldomain.col_lower_[i->column] - i->boundval, i->boundval);
    } else {
      if (val == 1) {
        if (globaldomain.col_upper_[i->column] != kHighsInf)
          addVUB(i->column, col,
                 i->boundval - globaldomain.col_upper_[i->column],
                 globaldomain.col_upper_[i->column]);
      } else
        addVUB(i->column,
               col,  // in case the upper bound is infinite the varbound can
                     // still be tightened as soon as a finite upper bound is
                     // known because the offset is finite
               globaldomain.col_upper_[i->column] - i->boundval, i->boundval);
    }
  }

  ImplIdx idx{col, val};
  hasProbed[idx] = true;
  for (auto i = implics.begin(); i != implics.begin() + safeImplicsEnd; ++i) {
    Implication implication;
    if (i->boundtype == HighsBoundType::kLower) {
      implication.lb = i->boundval;
    } else {
      implication.ub = i->boundval;
    }
    addImplication(idx, i->column, implication);
  }

  pdqsort(implics.begin(), implics.end(),
          [](const HighsDomainChange& a, const HighsDomainChange& b) {
            return a.column < b.column;
          });

  return false;
}

static constexpr bool kSkipBadVbds = true;
static constexpr bool kUseDualsForBreakingTies = true;

std::pair<HighsInt, HighsImplications::VarBound> HighsImplications::getBestVub(
    HighsInt col, const HighsSolution& lpSolution, double& bestUb,
    const HighsDomain& globaldom) const {
  std::pair<HighsInt, VarBound> bestVub =
      std::make_pair(-1, VarBound{0.0, kHighsInf, -1});
  double minbestUb = bestUb;
  double bestUbDist = kHighsInf;
  int64_t bestvubnodes = 0;

  auto isVubBetter = [&](double ubDist, int64_t vubNodes, double minVubVal,
                         HighsInt vubCol, const VarBound& vub) {
    if (ubDist < bestUbDist - mipsolver.mipdata_->feastol) return true;
    if (vubNodes > bestvubnodes) return true;
    if (vubNodes < bestvubnodes) return false;
    if (minVubVal < minbestUb - mipsolver.mipdata_->feastol) return true;
    if (kUseDualsForBreakingTies) {
      if (minVubVal > minbestUb + mipsolver.mipdata_->feastol) return false;
      if (lpSolution.col_dual[vubCol] / vub.coef -
              lpSolution.col_dual[bestVub.first] / bestVub.second.coef >
          mipsolver.mipdata_->feastol)
        return true;
    }

    return false;
  };

  double scale = globaldom.col_upper_[col] - globaldom.col_lower_[col];
  if (scale == kHighsInf)
    scale = 1.0;
  else
    scale = 1.0 / scale;

  vubs[col].for_each([&](HighsInt vubCol, const VarBound& vub) {
    if (vub.coef == kHighsInf) return;
    if (globaldom.isFixed(vubCol)) return;
    assert(globaldom.isBinary(vubCol));
    double vubval = lpSolution.col_value[vubCol] * vub.coef + vub.constant;
    double ubDist = std::max(0.0, vubval - lpSolution.col_value[col]);

    double yDist = mipsolver.mipdata_->feastol +
                   (vub.coef > 0 ? 1 - lpSolution.col_value[vubCol]
                                 : lpSolution.col_value[vubCol]);
    // skip variable bound if the distance towards the variable bound constraint
    // is larger than the distance to the point where the binary column is
    // relaxed to it's weakest bound, i.e. 1 if its coefficient is positive.
    // The variable bound constraint has the form x <= ay + b with y binary and
    // hence the norm is sqrt(1 + a^2) and the distance of the variable bound
    // constraint is ay + b - x evaluated at the solution values of x and y
    // divided by the norm.
    double norm2 = 1.0 + vub.coef * vub.coef;
    if (kSkipBadVbds && ubDist * ubDist > yDist * yDist * norm2) return;

    assert(vubCol >= 0 && vubCol < mipsolver.numCol());
    ubDist *= scale;
    if (ubDist <= bestUbDist + mipsolver.mipdata_->feastol) {
      double minvubval = vub.minValue();
      int64_t vubnodes =
          vub.coef > 0 ? mipsolver.mipdata_->nodequeue.numNodesDown(vubCol)
                       : mipsolver.mipdata_->nodequeue.numNodesUp(vubCol);

      if (isVubBetter(ubDist, vubnodes, minvubval, vubCol, vub)) {
        bestUb = vubval;
        minbestUb = minvubval;
        bestVub = std::make_pair(vubCol, vub);
        bestvubnodes = vubnodes;
        bestUbDist = ubDist;
      }
    }
  });

  return bestVub;
}

std::pair<HighsInt, HighsImplications::VarBound> HighsImplications::getBestVlb(
    HighsInt col, const HighsSolution& lpSolution, double& bestLb,
    const HighsDomain& globaldom) const {
  std::pair<HighsInt, VarBound> bestVlb =
      std::make_pair(-1, VarBound{0.0, -kHighsInf, -1});
  double maxbestlb = bestLb;
  double bestLbDist = kHighsInf;
  int64_t bestvlbnodes = 0;

  auto isVlbBetter = [&](double lbDist, int64_t vlbNodes, double maxVlbVal,
                         HighsInt vlbCol, const VarBound& vlb) {
    if (lbDist < bestLbDist - mipsolver.mipdata_->feastol) return true;
    if (vlbNodes > bestvlbnodes) return true;
    if (vlbNodes < bestvlbnodes) return false;
    if (maxVlbVal > maxbestlb + mipsolver.mipdata_->feastol) return true;
    if (kUseDualsForBreakingTies) {
      if (maxVlbVal < maxbestlb - mipsolver.mipdata_->feastol) return false;
      if (lpSolution.col_dual[vlbCol] / vlb.coef -
              lpSolution.col_dual[bestVlb.first] / bestVlb.second.coef <
          -mipsolver.mipdata_->feastol)
        return true;
    }

    return false;
  };

  double scale = globaldom.col_upper_[col] - globaldom.col_lower_[col];
  if (scale == kHighsInf)
    scale = 1.0;
  else
    scale = 1.0 / scale;

  vlbs[col].for_each([&](HighsInt vlbCol, const VarBound& vlb) {
    if (vlb.coef == -kHighsInf) return;
    if (globaldom.isFixed(vlbCol)) return;
    assert(globaldom.isBinary(vlbCol));
    assert(vlbCol >= 0 && vlbCol < mipsolver.numCol());
    double vlbval = lpSolution.col_value[vlbCol] * vlb.coef + vlb.constant;
    double lbDist = std::max(0.0, lpSolution.col_value[col] - vlbval);

    double yDist = mipsolver.mipdata_->feastol +
                   (vlb.coef > 0 ? lpSolution.col_value[vlbCol]
                                 : 1 - lpSolution.col_value[vlbCol]);

    double norm2 = 1.0 + vlb.coef * vlb.coef;
    if (kSkipBadVbds && lbDist * lbDist > yDist * yDist * norm2) return;

    // scale the distance as if the bounded column was scaled to have ub-lb=1
    lbDist *= scale;
    if (lbDist <= bestLbDist + mipsolver.mipdata_->feastol) {
      double maxvlbval = vlb.maxValue();
      int64_t vlbnodes =
          vlb.coef > 0 ? mipsolver.mipdata_->nodequeue.numNodesUp(vlbCol)
                       : mipsolver.mipdata_->nodequeue.numNodesDown(vlbCol);

      if (isVlbBetter(lbDist, vlbnodes, maxvlbval, vlbCol, vlb)) {
        bestLb = vlbval;
        maxbestlb = maxvlbval;
        bestVlb = std::make_pair(vlbCol, vlb);
        bestvlbnodes = vlbnodes;
        bestLbDist = lbDist;
      }
    }
  });

  return bestVlb;
}

bool HighsImplications::runProbing(HighsInt col, HighsInt& numReductions) {
  HighsDomain& globaldomain = mipsolver.mipdata_->getDomain();
  if (globaldomain.isBinary(col) && !probedBefore(col, 1) &&
      !probedBefore(col, 0) &&
      mipsolver.mipdata_->cliquetable.getSubstitution(col) == nullptr) {
    const bool dualFixProbingActive = globaldomain.getDualFixProbingActive();
    if (dualFixProbingActive) {
      clearTentativeCliques();
      globaldomain.getDualFixProbingPropagation().setZeroCostFixingPosition(
          kHighsIInf);
    }

    bool infeasible = computeImplications(col, 1);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;
    if (mipsolver.mipdata_->cliquetable.getSubstitution(col) != nullptr)
      return true;

    infeasible = computeImplications(col, 0);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;
    if (mipsolver.mipdata_->cliquetable.getSubstitution(col) != nullptr)
      return true;

    if (dualFixProbingActive && !dualFixProbingBinInds_.empty()) {
      HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
      HighsCliqueTable::CliqueVar clique[2];
      for (HighsInt k : dualFixProbingBinInds_) {
        if (!globaldomain.isBinary(k) || colsubstituted[k]) continue;
        if (globaldomain.infeasible()) return true;
        const TentativeFixing& f = dualFixProbingBinFlags_[k];
        if (f.isUndecided()) continue;
        if (f.downProbe == TentativeFixing::FixLower &&
            f.upProbe == TentativeFixing::FixLower) {
          globaldomain.fixCol(k, globaldomain.col_lower_[k]);
        } else if (f.downProbe == TentativeFixing::FixUpper &&
                   f.upProbe == TentativeFixing::FixUpper) {
          globaldomain.fixCol(k, globaldomain.col_upper_[k]);
        } else if (f.downProbe == TentativeFixing::FixUpper &&
                   f.upProbe == TentativeFixing::FixLower &&
                   !cliquetable.isFull()) {
          clique[0] = HighsCliqueTable::CliqueVar(col, 1);
          clique[1] = HighsCliqueTable::CliqueVar(k, 1);
          cliquetable.addClique(mipsolver, &clique[0], 2);
          clique[0] = HighsCliqueTable::CliqueVar(col, 0);
          clique[1] = HighsCliqueTable::CliqueVar(k, 0);
          cliquetable.addClique(mipsolver, &clique[0], 2);
        } else if (f.downProbe == TentativeFixing::FixLower &&
                   f.upProbe == TentativeFixing::FixUpper &&
                   !cliquetable.isFull()) {
          clique[0] = HighsCliqueTable::CliqueVar(col, 1);
          clique[1] = HighsCliqueTable::CliqueVar(k, 0);
          cliquetable.addClique(mipsolver, &clique[0], 2);
          clique[0] = HighsCliqueTable::CliqueVar(col, 0);
          clique[1] = HighsCliqueTable::CliqueVar(k, 1);
          cliquetable.addClique(mipsolver, &clique[0], 2);
        }
        if (globaldomain.infeasible()) return true;
      }

      clearTentativeCliques();
    }

    HighsInt nimplicsdown = implicationsDown.size();
    HighsInt nimplicsup = implicationsUp.size();
    HighsInt u = 0;
    HighsInt d = 0;

    while (u < nimplicsup && d < nimplicsdown) {
      if (implicationsUp[u].column < implicationsDown[d].column)
        ++u;
      else if (implicationsDown[d].column < implicationsUp[u].column)
        ++d;
      else {
        assert(implicationsUp[u].column == implicationsDown[d].column);
        HighsInt implcol = implicationsUp[u].column;
        double lbDown = globaldomain.col_lower_[implcol];
        double ubDown = globaldomain.col_upper_[implcol];
        double lbUp = lbDown;
        double ubUp = ubDown;

        do {
          if (implicationsDown[d].boundtype == HighsBoundType::kLower) {
            lbDown = std::max(lbDown, implicationsDown[d].boundval);
          } else {
            ubDown = std::min(ubDown, implicationsDown[d].boundval);
          }
          ++d;
        } while (d < nimplicsdown && implicationsDown[d].column == implcol);

        do {
          if (implicationsUp[u].boundtype == HighsBoundType::kLower) {
            lbUp = std::max(lbUp, implicationsUp[u].boundval);
          } else {
            ubUp = std::min(ubUp, implicationsUp[u].boundval);
          }
          ++u;
        } while (u < nimplicsup && implicationsUp[u].column == implcol);

        if (colsubstituted[implcol] || globaldomain.isFixed(implcol)) continue;

        if (lbDown == ubDown && lbUp == ubUp &&
            std::abs(lbDown - lbUp) > mipsolver.mipdata_->feastol) {
          HighsSubstitution substitution;
          substitution.substcol = implcol;
          substitution.staycol = col;
          substitution.offset = lbDown;
          substitution.scale = lbUp - lbDown;
          substitutions.push_back(substitution);
          colsubstituted[implcol] = true;
          ++numReductions;
        } else if (!mipsolver.mipdata_->parallelLockActive()) {
          double lb = std::min(lbDown, lbUp);
          double ub = std::max(ubDown, ubUp);

          if (lb > globaldomain.col_lower_[implcol]) {
            globaldomain.changeBound(HighsBoundType::kLower, implcol, lb,
                                     HighsDomain::Reason::unspecified());
            ++numReductions;
            if (globaldomain.infeasible()) return true;
          }

          if (ub < globaldomain.col_upper_[implcol]) {
            globaldomain.changeBound(HighsBoundType::kUpper, implcol, ub,
                                     HighsDomain::Reason::unspecified());
            ++numReductions;
            if (globaldomain.infeasible()) return true;
          }
        }
      }
    }

    hasProbed[ImplIdx{col, 0}] = true;
    hasProbed[ImplIdx{col, 1}] = true;
    return true;
  }

  return false;
}

void HighsImplications::addImplication(ImplIdx idx, HighsInt implCol,
                                       Implication implic) {
  auto insertresult = implications[idx].insert_or_get(implCol, implic);

  if (!insertresult.second) {
    if (insertresult.first->lb == -kHighsInf && implic.lb != -kHighsInf)
      ++numImplications;
    if (insertresult.first->ub == kHighsInf && implic.ub != kHighsInf)
      ++numImplications;
    insertresult.first->lb = std::max(implic.lb, insertresult.first->lb);
    insertresult.first->ub = std::min(implic.ub, insertresult.first->ub);
  } else {
    if (implic.lb != -kHighsInf) ++numImplications;
    if (implic.ub != kHighsInf) ++numImplications;
  }

  reverseImplications[implCol].insert_or_get(idx.col);
}

void HighsImplications::strengthenVarBound(VarBound& vbnd,
                                           HighsInt multiplier) const {
  // try to strengthen variable bound constraint x + a * y <= b by computing
  // MIR cut. since x is a general-integer variable (with integral bounds) and
  // its coefficient is 1.0, it is not shifted (or complemented). similarly,
  // since y is a binary variable (i.e., its lower bound is 0.0), it does not
  // need to be shifted, and we also do not try to complement it.
  if (std::abs(vbnd.coef) == kHighsInf || std::abs(vbnd.constant) == kHighsInf)
    return;
  constexpr double f0min = 0.005;
  constexpr double f0max = 0.995;
  double downrhs = std::floor(multiplier * vbnd.constant);
  double f0 = multiplier * vbnd.constant - downrhs;
  if (f0 < f0min || f0 > f0max) return;
  double downaj = std::floor(-multiplier * vbnd.coef + kHighsTiny);
  double fj = -multiplier * vbnd.coef - downaj;
  vbnd.constant = multiplier * downrhs;
  vbnd.coef = -multiplier * (downaj + std::max(fj - f0, 0.0) / (1.0 - f0));
};

void HighsImplications::rowModified(HighsInt row) {
  rowToVarBounds[row].for_each([&](HighsInt col, HighsInt binaryCol) {
    VarBound* vb = vlbs[col].find(binaryCol);
    if (vb && vb->origin == row) vb->origin = -1;
    vb = vubs[col].find(binaryCol);
    if (vb && vb->origin == row) vb->origin = -1;
  });
  rowToVarBounds[row].clear();
}

void HighsImplications::addVUB(HighsInt col, HighsInt vubcol, double vubcoef,
                               double vubconstant, HighsInt origin) {
  addVUB(col, vubcol, vubcoef, vubconstant,
         mipsolver.mipdata_->getDomain().col_upper_[col],
         mipsolver.isColIntegral(col), origin);
}

void HighsImplications::addVUB(HighsInt col, HighsInt vubcol, double vubcoef,
                               double vubconstant, double colupperbound,
                               bool colisintegral, HighsInt origin) {
  // assume that VUBs do not have infinite coefficients and infinite constant
  // terms since such VUBs effectively evaluate to NaN.
  assert(std::abs(vubcoef) != kHighsInf || std::abs(vubconstant) != kHighsInf);
  if (tooManyVarBounds()) return;

  VarBound vub{vubcoef, vubconstant, origin};

  if (colisintegral) {
    // try to strengthen VUB
    strengthenVarBound(vub, HighsInt{1});
    if (vub.coef == 0.0) return;
  }

  mipsolver.mipdata_->debugSolution.checkVub(col, vubcol, vubcoef, vubconstant);

  double minBound = vub.minValue();
  if (minBound >= colupperbound - mipsolver.mipdata_->feastol) return;

  auto insertresult = vubs[col].insert_or_get(vubcol, vub);

  if (!insertresult.second) {
    VarBound& currentvub = *insertresult.first;
    double currentMinBound = currentvub.minValue();
    if (minBound < currentMinBound - mipsolver.mipdata_->feastol) {
      currentvub.coef = vub.coef;
      currentvub.constant = vub.constant;
      currentvub.origin = vub.origin;
    }
  } else
    numVarBounds++;

  if (origin >= 0) rowToVarBounds[origin].insert(col, vubcol);
}

void HighsImplications::addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef,
                               double vlbconstant, HighsInt origin) {
  addVLB(col, vlbcol, vlbcoef, vlbconstant,
         mipsolver.mipdata_->getDomain().col_lower_[col],
         mipsolver.isColIntegral(col), origin);
}

void HighsImplications::addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef,
                               double vlbconstant, double colllowerbound,
                               bool colisintegral, HighsInt origin) {
  // assume that VLBs do not have infinite coefficients and infinite constant
  // terms since such VLBs effectively evaluate to NaN.
  assert(std::abs(vlbcoef) != kHighsInf || std::abs(vlbconstant) != kHighsInf);
  if (tooManyVarBounds()) return;

  VarBound vlb{vlbcoef, vlbconstant, origin};

  if (colisintegral) {
    // try to strengthen VLB
    strengthenVarBound(vlb, HighsInt{-1});
    if (vlb.coef == 0.0) return;
  }

  mipsolver.mipdata_->debugSolution.checkVlb(col, vlbcol, vlbcoef, vlbconstant);

  double maxBound = vlb.maxValue();
  if (maxBound <= colllowerbound + mipsolver.mipdata_->feastol) return;

  auto insertresult = vlbs[col].insert_or_get(vlbcol, vlb);

  if (!insertresult.second) {
    VarBound& currentvlb = *insertresult.first;

    double currentMaxBound = currentvlb.maxValue();
    if (maxBound > currentMaxBound + mipsolver.mipdata_->feastol) {
      currentvlb.coef = vlb.coef;
      currentvlb.constant = vlb.constant;
      currentvlb.origin = vlb.origin;
    }
  } else
    numVarBounds++;

  if (origin >= 0) rowToVarBounds[origin].insert(col, vlbcol);
}

void HighsImplications::reindexRows(HighsInt nrows,
                                    const std::vector<HighsInt>& newRowIndex) {
  rowToVarBounds.clear();
  rowToVarBounds.resize(nrows);

  for (HighsInt col = 0; col < static_cast<HighsInt>(vlbs.size()); ++col) {
    vlbs[col].for_each([&](HighsInt vlbCol, VarBound& vlb) {
      if (vlb.origin >= 0) {
        vlb.origin = newRowIndex[vlb.origin];
        if (vlb.origin >= 0) rowToVarBounds[vlb.origin].insert(col, vlbCol);
      }
    });
    vubs[col].for_each([&](HighsInt vubCol, VarBound& vub) {
      if (vub.origin >= 0) {
        vub.origin = newRowIndex[vub.origin];
        if (vub.origin >= 0) rowToVarBounds[vub.origin].insert(col, vubCol);
      }
    });
  }
}

void HighsImplications::rebuild(HighsInt ncols,
                                const std::vector<HighsInt>& orig2reducedcol,
                                const std::vector<HighsInt>& orig2reducedrow) {
  std::vector<HighsHashTree<HighsInt, VarBound>> oldvubs;
  std::vector<HighsHashTree<HighsInt, VarBound>> oldvlbs;

  oldvlbs.swap(vlbs);
  oldvubs.swap(vubs);

  std::vector<HighsHashTree<HighsInt, Implication>> oldimplications;
  oldimplications.swap(implications);

  colsubstituted.clear();
  colsubstituted.shrink_to_fit();
  implications.clear();
  implications.shrink_to_fit();
  hasProbed.clear();
  hasProbed.shrink_to_fit();
  reverseImplications.clear();
  reverseImplications.shrink_to_fit();

  implications.resize(2 * ncols);
  hasProbed.resize(2 * ncols);
  reverseImplications.resize(ncols);
  colsubstituted.resize(ncols);
  substitutions.clear();
  vubs.clear();
  vubs.shrink_to_fit();
  vubs.resize(ncols);
  vlbs.clear();
  vlbs.shrink_to_fit();
  vlbs.resize(ncols);
  dualFixProbingBinInds_.clear();
  dualFixProbingBinInds_.reserve(ncols);
  dualFixProbingBinFlags_.assign(ncols, TentativeFixing{});
  rowToVarBounds.clear();
  rowToVarBounds.shrink_to_fit();
  rowToVarBounds.resize(mipsolver.numRow());
  numImplications = 0;
  numVarBounds = 0;
  HighsInt oldncols = oldvubs.size();

  nextCleanupCall = mipsolver.numNonzero();

  for (HighsInt i = 0; i != oldncols; ++i) {
    if (int(i) >= int(orig2reducedcol.size())) {
      printf("HighsImplications::rebuild i = %d orig2reducedcol.size = %d\n",
             int(i), int(orig2reducedcol.size()));
      assert(111 == 345);
    }
    HighsInt newi = orig2reducedcol[i];

    if (newi == -1 ||
        !mipsolver.mipdata_->postSolveStack.isColLinearlyTransformable(newi))
      continue;

    oldvubs[i].for_each([&](HighsInt vubCol, VarBound vub) {
      HighsInt newVubCol = orig2reducedcol[vubCol];
      if (newVubCol == -1) return;

      if (!mipsolver.mipdata_->getDomain().isBinary(newVubCol) ||
          !mipsolver.mipdata_->postSolveStack.isColLinearlyTransformable(
              newVubCol))
        return;

      addVUB(newi, newVubCol, vub.coef, vub.constant,
             vub.origin >= 0 ? orig2reducedrow[vub.origin] : -1);
    });

    oldvlbs[i].for_each([&](HighsInt vlbCol, VarBound vlb) {
      HighsInt newVlbCol = orig2reducedcol[vlbCol];
      if (newVlbCol == -1) return;

      if (!mipsolver.mipdata_->getDomain().isBinary(newVlbCol) ||
          !mipsolver.mipdata_->postSolveStack.isColLinearlyTransformable(
              newVlbCol))
        return;

      addVLB(newi, newVlbCol, vlb.coef, vlb.constant,
             vlb.origin >= 0 ? orig2reducedrow[vlb.origin] : -1);
    });

    if (mipsolver.mipdata_->getDomain().isBinary(newi)) {
      for (HighsInt val = 0; val != 2; val++) {
        oldimplications[ImplIdx{i, val}].for_each([&](HighsInt implCol,
                                                      Implication impl) {
          HighsInt newImplCol = orig2reducedcol[implCol];
          if (newImplCol == -1 ||
              !mipsolver.mipdata_->postSolveStack.isColLinearlyTransformable(
                  newImplCol))
            return;
          addImplication(ImplIdx{newi, val}, newImplCol, impl);
        });
      }
    }
  }
}

void HighsImplications::buildFrom(const HighsImplications& init) {
  // todo check if this should be done
  HighsInt numcol = mipsolver.numCol();

  for (HighsInt i = 0; i != numcol; ++i) {
    init.vubs[i].for_each([&](HighsInt vubCol, VarBound vub) {
      if (!mipsolver.mipdata_->getDomain().isBinary(vubCol)) return;
      addVUB(i, vubCol, vub.coef, vub.constant, vub.origin);
    });

    init.vlbs[i].for_each([&](HighsInt vlbCol, VarBound vlb) {
      if (!mipsolver.mipdata_->getDomain().isBinary(vlbCol)) return;
      addVLB(i, vlbCol, vlb.coef, vlb.constant, vlb.origin);
    });

    if (mipsolver.mipdata_->getDomain().isBinary(i)) {
      for (HighsInt val = 0; val != 2; val++) {
        init.implications[ImplIdx{i, val}].for_each(
            [&](HighsInt implCol, Implication implic) {
              addImplication(ImplIdx{i, val}, implCol, implic);
            });
      }
    }
  }
}

void HighsImplications::separateImpliedBounds(
    const HighsLpRelaxation& lpRelaxation, const std::vector<double>& sol,
    HighsCutPool& cutpool, double feastol, HighsDomain& globaldom,
    bool thread_safe) {
  std::array<HighsInt, 2> inds;
  std::array<double, 2> vals;
  double rhs;

  HighsInt numboundchgs = 0;

  // first do probing on all candidates that have not been probed yet
  if (!mipsolver.mipdata_->cliquetable.isFull() && !thread_safe) {
    auto oldNumQueries =
        mipsolver.mipdata_->cliquetable.numNeighbourhoodQueries;
    HighsInt oldNumEntries = mipsolver.mipdata_->cliquetable.getNumEntries();

    for (std::pair<HighsInt, double> fracint :
         lpRelaxation.getFractionalIntegers()) {
      HighsInt col = fracint.first;
      if (globaldom.col_lower_[col] != 0.0 ||
          globaldom.col_upper_[col] != 1.0 ||
          (probedBefore(col, 0) && probedBefore(col, 1)))
        continue;

      mipsolver.profiling_->start(kMipClockProbingImplications);
      const bool probing_result = runProbing(col, numboundchgs);
      mipsolver.profiling_->stop(kMipClockProbingImplications);
      if (probing_result) {
        if (globaldom.infeasible()) return;
      }

      if (mipsolver.mipdata_->cliquetable.isFull()) break;
    }

    // if (!mipsolver.submip)
    //   printf("numEntries: %d, beforeProbing: %d\n",
    //          mipsolver.mipdata_->cliquetable.getNumEntries(), oldNumEntries);
    HighsInt numNewEntries =
        mipsolver.mipdata_->cliquetable.getNumEntries() - oldNumEntries;

    nextCleanupCall -= std::max(HighsInt{0}, numNewEntries);

    if (nextCleanupCall < 0) {
      // HighsInt oldNumEntries =
      // mipsolver.mipdata_->cliquetable.getNumEntries();
      if (!mipsolver.mipdata_->parallelLockActive())
        mipsolver.mipdata_->cliquetable.runCliqueMerging(globaldom);

      // printf("numEntries: %d, beforeMerging: %d\n",
      //        mipsolver.mipdata_->cliquetable.getNumEntries(), oldNumEntries);
      nextCleanupCall =
          std::min(mipsolver.mipdata_->numCliqueEntriesAfterFirstPresolve,
                   mipsolver.mipdata_->cliquetable.getNumEntries());
      // printf("nextCleanupCall: %d\n", nextCleanupCall);
    }

    if (!mipsolver.mipdata_->parallelLockActive())
      mipsolver.mipdata_->cliquetable.numNeighbourhoodQueries = oldNumQueries;
  }

  auto tryAddCut = [&](HighsInt implCol) {
    double viol = sol[inds[0]] * vals[0] + sol[inds[1]] * vals[1] - rhs;
    if (viol > feastol) {
      cutpool.addCut(mipsolver, inds.data(), vals.data(), 2, rhs,
                     !mipsolver.isColContinuous(implCol), false, false, false);
    }
  };

  for (std::pair<HighsInt, double> fracint :
       lpRelaxation.getFractionalIntegers()) {
    HighsInt col = fracint.first;
    // skip non binary variables
    if (globaldom.col_lower_[col] != 0.0 || globaldom.col_upper_[col] != 1.0)
      continue;

    for (HighsInt val = 0; val != 2; val++) {
      implications[ImplIdx{col, val}].for_each(
          [&](HighsInt implCol, Implication implic) {
            if (val == 1) {
              if (implic.ub + feastol < globaldom.col_upper_[implCol]) {
                vals[0] = 1.0;
                inds[0] = implCol;
                vals[1] = globaldom.col_upper_[implCol] - implic.ub;
                inds[1] = col;
                rhs = globaldom.col_upper_[implCol];
                tryAddCut(implCol);
              }
              if (implic.lb - feastol > globaldom.col_lower_[implCol]) {
                vals[0] = -1.0;
                inds[0] = implCol;
                vals[1] = implic.lb - globaldom.col_lower_[implCol];
                inds[1] = col;
                rhs = -globaldom.col_lower_[implCol];
                tryAddCut(implCol);
              }
            } else {
              if (implic.ub + feastol < globaldom.col_upper_[implCol]) {
                vals[0] = 1.0;
                inds[0] = implCol;
                vals[1] = implic.ub - globaldom.col_upper_[implCol];
                inds[1] = col;
                rhs = implic.ub;
                tryAddCut(implCol);
              }
              if (implic.lb - feastol > globaldom.col_lower_[implCol]) {
                vals[0] = -1.0;
                inds[0] = implCol;
                vals[1] = globaldom.col_lower_[implCol] - implic.lb;
                inds[1] = col;
                rhs = -implic.lb;
                tryAddCut(implCol);
              }
            }
          });
    }
  }
}

void HighsImplications::cleanupVarbounds(HighsInt col) {
  double ub = mipsolver.mipdata_->getDomain().col_upper_[col];
  double lb = mipsolver.mipdata_->getDomain().col_lower_[col];

  if (ub == lb) {
    HighsInt numVubs = 0;
    vubs[col].for_each([&](HighsInt vubCol, VarBound& vub) { numVubs++; });
    HighsInt numVlbs = 0;
    vlbs[col].for_each([&](HighsInt vlbCol, VarBound& vlb) { numVlbs++; });
    numVarBounds -= numVubs + numVlbs;
    vlbs[col].clear();
    vubs[col].clear();
    return;
  }

  std::vector<HighsInt> delVbds;

  vubs[col].for_each([&](HighsInt vubCol, VarBound& vub) {
    bool redundant = false;
    bool infeasible = false;
    cleanupVub(col, vubCol, vub, ub, redundant, infeasible);
    if (redundant) delVbds.push_back(vubCol);
    if (infeasible) return;
  });

  for (HighsInt vubCol : delVbds) vubs[col].erase(vubCol);
  numVarBounds -= delVbds.size();
  delVbds.clear();

  vlbs[col].for_each([&](HighsInt vlbCol, VarBound& vlb) {
    bool redundant = false;
    bool infeasible = false;
    cleanupVlb(col, vlbCol, vlb, lb, redundant, infeasible);
    if (redundant) delVbds.push_back(vlbCol);
    if (infeasible) return;
  });

  for (HighsInt vlbCol : delVbds) vlbs[col].erase(vlbCol);
  numVarBounds -= delVbds.size();
}

bool HighsImplications::redundantVlb(const VarBound& vlb, double lb) const {
  return vlb.maxValue() <= lb + mipsolver.mipdata_->feastol;
}

bool HighsImplications::redundantVub(const VarBound& vub, double ub) const {
  return vub.minValue() >= ub - mipsolver.mipdata_->feastol;
}

bool HighsImplications::tightenVlb(VarBound& vlb, double lb) const {
  if (redundantVlb(vlb, lb) ||
      vlb.minValue() >= lb - mipsolver.mipdata_->epsilon)
    return false;
  double newcoef =
      static_cast<double>(lb - static_cast<HighsCDouble>(vlb.maxValue()));
  if (vlb.coef < 0)
    vlb.coef = newcoef;
  else {
    vlb.constant = lb;
    vlb.coef = -newcoef;
  }
  return true;
}

bool HighsImplications::tightenVub(VarBound& vub, double ub) const {
  if (redundantVub(vub, ub) ||
      vub.maxValue() <= ub + mipsolver.mipdata_->epsilon)
    return false;
  double newcoef =
      static_cast<double>(ub - static_cast<HighsCDouble>(vub.minValue()));
  if (vub.coef > 0)
    vub.coef = newcoef;
  else {
    vub.constant = ub;
    vub.coef = -newcoef;
  }
  return true;
}

void HighsImplications::cleanupVlb(HighsInt col, HighsInt vlbCol,
                                   HighsImplications::VarBound& vlb, double lb,
                                   bool& redundant, bool& infeasible,
                                   bool allowBoundChanges) const {
  // initialize
  redundant = false;
  infeasible = false;

  // return if there is no variable bound
  if (vlbCol == -1) return;

  // check variable lower bound
  mipsolver.mipdata_->debugSolution.checkVlb(col, vlbCol, vlb.coef,
                                             vlb.constant);

  if (redundantVlb(vlb, lb)) {
    // variable bound is redundant
    redundant = true;
  } else if (tightenVlb(vlb, lb)) {
    // coefficient can be tightened
    // check tightened variable lower bound
    mipsolver.mipdata_->debugSolution.checkVlb(col, vlbCol, vlb.coef,
                                               vlb.constant);
  } else if (allowBoundChanges &&
             vlb.minValue() > lb + mipsolver.mipdata_->epsilon) {
    mipsolver.mipdata_->getDomain().changeBound(
        HighsBoundType::kLower, col, vlb.minValue(),
        HighsDomain::Reason::unspecified());
    infeasible = mipsolver.mipdata_->getDomain().infeasible();
  }
}

void HighsImplications::cleanupVub(HighsInt col, HighsInt vubCol,
                                   HighsImplications::VarBound& vub, double ub,
                                   bool& redundant, bool& infeasible,
                                   bool allowBoundChanges) const {
  // initialize
  redundant = false;
  infeasible = false;

  // return if there is no variable bound
  if (vubCol == -1) return;

  // check variable upper bound
  mipsolver.mipdata_->debugSolution.checkVub(col, vubCol, vub.coef,
                                             vub.constant);

  if (redundantVub(vub, ub)) {
    // variable bound is redundant
    redundant = true;
  } else if (tightenVub(vub, ub)) {
    // coefficient can be tightened
    // check tightened variable upper bound
    mipsolver.mipdata_->debugSolution.checkVub(col, vubCol, vub.coef,
                                               vub.constant);
  } else if (allowBoundChanges &&
             vub.maxValue() < ub - mipsolver.mipdata_->epsilon) {
    mipsolver.mipdata_->getDomain().changeBound(
        HighsBoundType::kUpper, col, vub.maxValue(),
        HighsDomain::Reason::unspecified());
    infeasible = mipsolver.mipdata_->getDomain().infeasible();
  }
}

void HighsImplications::applyImplications(HighsDomain& domain,
                                          const HighsInt col,
                                          const HighsInt val) {
  assert(domain.isFixed(col));

  auto checkImplication = [&](const HighsInt implcol,
                              const Implication& implic) -> bool {
    assert(!domain.infeasible());
    if (domain.isFixed(implcol)) return false;
    const bool isint =
        domain.variableType(implcol) != HighsVarType::kContinuous;
    // Directly change bounds on all integer columns. Only change continuous
    // columns that fix the column, as changing their domains risks
    // suppressing further bound changes found in propagation, e.g.,
    // change [0, 100] -> [0, 50], propagation could tighten to [0, 48], but
    // such a tightening would not be applied due to min boundRange improvement.
    if ((!isint && implic.lb > domain.col_upper_[implcol] - domain.feastol()) ||
        (isint && implic.lb > domain.col_lower_[implcol] + domain.feastol())) {
      domain.changeBound(HighsBoundType::kLower, implcol, implic.lb,
                         HighsDomain::Reason::cliqueTable(col, val));
      if (domain.infeasible()) return true;
    }
    if ((!isint && implic.ub < domain.col_lower_[implcol] + domain.feastol()) ||
        (isint && implic.ub < domain.col_upper_[implcol] - domain.feastol())) {
      domain.changeBound(HighsBoundType::kUpper, implcol, implic.ub,
                         HighsDomain::Reason::cliqueTable(col, val));
    }
    return domain.infeasible();
  };

  implications[ImplIdx{col, val}].for_each(checkImplication);
}
