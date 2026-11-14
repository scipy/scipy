/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsSearch.h"

#include <numeric>
#include <tuple>

#include "lp_data/HConst.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolverData.h"

HighsSearch::HighsSearch(HighsMipWorker& mipworker, HighsPseudocost& pseudocost)
    : mipworker(mipworker),
      mipsolver(mipworker.getMipSolver()),
      lp(nullptr),
      localdom(mipworker.getGlobalDomain()),
      pseudocost(pseudocost) {
  nnodes = 0;
  nleaves = 0;
  treeweight = 0.0;
  depthoffset = 0;
  lpiterations = 0;
  heurlpiterations = 0;
  sblpiterations = 0;
  upper_limit = kHighsInf;
  inheuristic = false;
  inbranching = false;
  countTreeWeight = true;
  childselrule = mipsolver.submip ? ChildSelectionRule::kHybridInferenceCost
                                  : ChildSelectionRule::kRootSol;
  // the infeasibility flag is overwritten and lost when setDomainChangeStack is
  // called. therefore, assert that localdom is not infeasible here.
  assert(!this->localdom.infeasible());
  this->localdom.setDomainChangeStack(std::vector<HighsDomainChange>());
}

double HighsSearch::checkSol(const std::vector<double>& sol,
                             bool& integerfeasible) const {
  HighsCDouble objval = 0.0;
  integerfeasible = true;
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    objval += static_cast<HighsCDouble>(sol[i]) * mipsolver.colCost(i);
    assert(std::isfinite(sol[i]));

    if (!integerfeasible || !mipsolver.isColInteger(i)) continue;

    if (fractionality(sol[i]) > getFeasTol()) {
      integerfeasible = false;
    }
  }

  return double(objval);
}

bool HighsSearch::orbitsValidInChildNode(
    const HighsDomainChange& branchChg) const {
  HighsInt branchCol = branchChg.column;
  // if the variable is integral or we are in an up branch the stabilizer only
  // stays valid if the column has been stabilized
  const NodeData& currNode = nodestack.back();
  if (!currNode.stabilizerOrbits ||
      currNode.stabilizerOrbits->orbitCols.empty() ||
      currNode.stabilizerOrbits->isStabilized(branchCol))
    return true;

  // a down branch stays valid if the variable is binary
  if (branchChg.boundtype == HighsBoundType::kUpper &&
      localdom.isGlobalBinary(branchChg.column))
    return true;

  return false;
}

double HighsSearch::getCutoffBound() const {
  return std::min(getUpperLimit(), upper_limit);
}

void HighsSearch::setRINSNeighbourhood(const std::vector<double>& basesol,
                                       const std::vector<double>& relaxsol) {
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (!mipsolver.isColInteger(i)) continue;
    if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

    double intval = std::floor(basesol[i] + 0.5);
    if (std::abs(relaxsol[i] - intval) < getFeasTol()) {
      if (localdom.col_lower_[i] < intval)
        localdom.changeBound(HighsBoundType::kLower, i,
                             std::min(intval, localdom.col_upper_[i]),
                             HighsDomain::Reason::unspecified());
      if (localdom.col_upper_[i] > intval)
        localdom.changeBound(HighsBoundType::kUpper, i,
                             std::max(intval, localdom.col_lower_[i]),
                             HighsDomain::Reason::unspecified());
    }
  }
}

void HighsSearch::setRENSNeighbourhood(const std::vector<double>& lpsol) {
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (!mipsolver.isColInteger(i)) continue;
    if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

    double downval = std::floor(lpsol[i] + getFeasTol());
    double upval = std::ceil(lpsol[i] - getFeasTol());

    if (localdom.col_lower_[i] < downval) {
      localdom.changeBound(HighsBoundType::kLower, i,
                           std::min(downval, localdom.col_upper_[i]),
                           HighsDomain::Reason::unspecified());
      if (localdom.infeasible()) return;
    }
    if (localdom.col_upper_[i] > upval) {
      localdom.changeBound(HighsBoundType::kUpper, i,
                           std::max(upval, localdom.col_lower_[i]),
                           HighsDomain::Reason::unspecified());
      if (localdom.infeasible()) return;
    }
  }
}

void HighsSearch::createNewNode() {
  nodestack.emplace_back();
  nodestack.back().domgchgStackPos = localdom.getDomainChangeStack().size();
}

void HighsSearch::cutoffNode() { nodestack.back().opensubtrees = 0; }

void HighsSearch::setMinReliable(HighsInt minreliable) {
  pseudocost.setMinReliable(minreliable);
}

void HighsSearch::branchDownwards(HighsInt col, double newub,
                                  double branchpoint) {
  openChildNode(col, newub, branchpoint, HighsBoundType::kUpper);
}

void HighsSearch::branchUpwards(HighsInt col, double newlb,
                                double branchpoint) {
  openChildNode(col, newlb, branchpoint, HighsBoundType::kLower);
}

void HighsSearch::openChildNode(HighsInt col, double boundval,
                                double branchpoint, HighsBoundType boundtype) {
  NodeData& currnode = nodestack.back();

  assert(currnode.opensubtrees == 2);
  assert(mipsolver.isColIntegral(col));

  currnode.opensubtrees = 1;
  currnode.branching_point = branchpoint;
  currnode.branchingdecision.column = col;
  currnode.branchingdecision.boundval = boundval;
  currnode.branchingdecision.boundtype = boundtype;

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  bool passStabilizerToChildNode =
      orbitsValidInChildNode(currnode.branchingdecision);
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(
      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,
      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);
  nodestack.back().domgchgStackPos = domchgPos;
}

void HighsSearch::flipBranchingDecision(NodeData& currnode) {
  bool fallbackbranch =
      currnode.branchingdecision.boundval == currnode.branching_point;
  if (currnode.branchingdecision.boundtype == HighsBoundType::kLower) {
    currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
    currnode.branchingdecision.boundval =
        std::floor(currnode.branchingdecision.boundval - 0.5);
  } else {
    currnode.branchingdecision.boundtype = HighsBoundType::kLower;
    currnode.branchingdecision.boundval =
        std::ceil(currnode.branchingdecision.boundval + 0.5);
  }

  if (fallbackbranch)
    currnode.branching_point = currnode.branchingdecision.boundval;
}

void HighsSearch::addBoundExceedingConflict() {
  if (getUpperLimit() != kHighsInf) {
    double rhs;
    if (lp->computeDualProof(getDomain(), getUpperLimit(), inds, vals, rhs)) {
      if (getDomain().infeasible()) return;
      localdom.conflictAnalysis(inds.data(), vals.data(), inds.size(), rhs,
                                getConflictPool(), mipworker.getGlobalDomain(),
                                pseudocost);

      HighsCutGeneration cutGen(*lp, getCutPool());
      mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                                 inds.size(), rhs);
      cutGen.generateConflict(localdom, mipworker.getGlobalDomain(), inds, vals,
                              rhs);
    }
  }
}

void HighsSearch::addInfeasibleConflict() {
  double rhs;
  if (lp->getLpSolver().getModelStatus() == HighsModelStatus::kObjectiveBound)
    lp->performAging();

  if (lp->computeDualInfProof(getDomain(), inds, vals, rhs)) {
    if (getDomain().infeasible()) return;
    // double minactlocal = 0.0;
    // double minactglobal = 0.0;
    // for (HighsInt i = 0; i < int(inds.size()); ++i) {
    //  if (vals[i] > 0.0) {
    //    minactlocal += localdom.col_lower_[inds[i]] * vals[i];
    //    minactglobal += globaldom.col_lower_[inds[i]] * vals[i];
    //  } else {
    //    minactlocal += localdom.col_upper_[inds[i]] * vals[i];
    //    minactglobal += globaldom.col_upper_[inds[i]] * vals[i];
    //  }
    //}
    // HighsInt oldnumcuts = cutpool.getNumCuts();
    localdom.conflictAnalysis(inds.data(), vals.data(), inds.size(), rhs,
                              getConflictPool(), mipworker.getGlobalDomain(),
                              pseudocost);

    HighsCutGeneration cutGen(*lp, getCutPool());
    mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                               inds.size(), rhs);
    cutGen.generateConflict(localdom, mipworker.getGlobalDomain(), inds, vals,
                            rhs);

    // if (cutpool.getNumCuts() > oldnumcuts) {
    //  printf(
    //      "added cut from infeasibility proof with local min activity %g, "
    //      "global min activity %g, and rhs %g\n",
    //      minactlocal, minactglobal, rhs);
    //} else {
    //  printf(
    //      "no cut found for infeasibility proof with local min activity %g, "
    //      "global min "
    //      " activity %g, and rhs % g\n ",
    //      minactlocal, minactglobal, rhs);
    //}
    // HighsInt cutind = cutpool.addCut(inds.data(), vals.data(), inds.size(),
    // rhs); localdom.cutAdded(cutind);
  }
}

HighsInt HighsSearch::selectBranchingCandidate(int64_t maxSbIters,
                                               double& downNodeLb,
                                               double& upNodeLb) {
  assert(!lp->getFractionalIntegers().empty());

  std::vector<double> upscore;
  std::vector<double> downscore;
  std::vector<HighsBool> upscorereliable;
  std::vector<HighsBool> downscorereliable;
  std::vector<double> upbound;
  std::vector<double> downbound;

  HighsInt numfrac = lp->getFractionalIntegers().size();
  const auto& fracints = lp->getFractionalIntegers();

  upscore.resize(numfrac, kHighsInf);
  downscore.resize(numfrac, kHighsInf);
  upbound.resize(numfrac, getCurrentLowerBound());
  downbound.resize(numfrac, getCurrentLowerBound());

  upscorereliable.resize(numfrac, false);
  downscorereliable.resize(numfrac, false);

  // initialize up and down scores of variables that have a
  // reliable pseudocost so that they do not get evaluated
  for (HighsInt k = 0; k != numfrac; ++k) {
    HighsInt col = fracints[k].first;
    double fracval = fracints[k].second;

    const double lower_residual =
        (fracval - localdom.col_lower_[col]) - getFeasTol();
    const bool lower_ok = lower_residual > 0;
    if (!lower_ok)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kError,
                   "HighsSearch::selectBranchingCandidate Error fracval = %g "
                   "<= %g = %g + %g = "
                   "localdom.col_lower_[col] + getFeasTol(): "
                   "Residual %g\n",
                   fracval, localdom.col_lower_[col] + getFeasTol(),
                   localdom.col_lower_[col], getFeasTol(), lower_residual);

    const double upper_residual =
        (localdom.col_upper_[col] - fracval) - getFeasTol();
    const bool upper_ok = upper_residual > 0;
    if (!upper_ok)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kError,
                   "HighsSearch::selectBranchingCandidate Error fracval = %g "
                   ">= %g = %g - %g = "
                   "localdom.col_upper_[col] - getFeasTol(): "
                   "Residual %g\n",
                   fracval, localdom.col_upper_[col] - getFeasTol(),
                   localdom.col_upper_[col], getFeasTol(), upper_residual);

    assert(lower_residual > -1e-12 && upper_residual > -1e-12);

    //    assert(fracval > localdom.col_lower_[col] +
    //    getFeasTol()); assert(fracval <
    //    localdom.col_upper_[col] - getFeasTol());

    if (pseudocost.isReliable(col)) {
      upscore[k] = pseudocost.getPseudocostUp(col, fracval);
      downscore[k] = pseudocost.getPseudocostDown(col, fracval);
      upscorereliable[k] = true;
      downscorereliable[k] = true;
    } else {
      int flags = branchingVarReliableAtNodeFlags(col);
      if (flags & kUpReliable) {
        upscore[k] = pseudocost.getPseudocostUp(col, fracval);
        upscorereliable[k] = true;
      }

      if (flags & kDownReliable) {
        downscore[k] = pseudocost.getPseudocostDown(col, fracval);
        downscorereliable[k] = true;
      }
    }
  }

  std::vector<HighsInt> evalqueue;
  evalqueue.resize(numfrac);
  std::iota(evalqueue.begin(), evalqueue.end(), 0);

  auto numNodesUp = [&](HighsInt k) {
    return getNodeQueue().numNodesUp(fracints[k].first);
  };

  auto numNodesDown = [&](HighsInt k) {
    return getNodeQueue().numNodesDown(fracints[k].first);
  };

  double minScore = getFeasTol();

  auto selectBestScore = [&](bool finalSelection) {
    HighsInt best = -1;
    double bestscore = -1.0;
    double bestnodes = -1.0;
    int64_t bestnumnodes = 0;

    double oldminscore = minScore;
    for (HighsInt k : evalqueue) {
      double score;

      if (upscore[k] <= oldminscore) upscorereliable[k] = true;
      if (downscore[k] <= oldminscore) downscorereliable[k] = true;

      double s = 1e-3 * std::min(upscorereliable[k] ? upscore[k] : 0,
                                 downscorereliable[k] ? downscore[k] : 0);
      minScore = std::max(s, minScore);

      if (upscore[k] <= oldminscore || downscore[k] <= oldminscore)
        score = pseudocost.getScore(fracints[k].first,
                                    std::min(upscore[k], oldminscore),
                                    std::min(downscore[k], oldminscore));
      else {
        score = upscore[k] == kHighsInf || downscore[k] == kHighsInf
                    ? finalSelection ? pseudocost.getScore(fracints[k].first,
                                                           fracints[k].second)
                                     : kHighsInf
                    : pseudocost.getScore(fracints[k].first, upscore[k],
                                          downscore[k]);
      }

      assert(score >= 0.0);
      int64_t upnodes = numNodesUp(k);
      int64_t downnodes = numNodesDown(k);
      double nodes = 0;
      int64_t numnodes = upnodes + downnodes;
      if (upnodes != 0 || downnodes != 0)
        nodes =
            (downnodes / (double)(numnodes)) * (upnodes / (double)(numnodes));
      if (score > bestscore || (score > bestscore - getFeasTol() &&
                                std::make_pair(nodes, numnodes) >
                                    std::make_pair(bestnodes, bestnumnodes))) {
        bestscore = score;
        best = k;
        bestnodes = nodes;
        bestnumnodes = numnodes;
      }
    }

    return best;
  };

  HighsLpRelaxation::Playground playground = lp->playground();

  while (true) {
    bool mustStop =
        getStrongBranchingLpIterations() >= maxSbIters || checkLimits();

    HighsInt candidate = selectBestScore(mustStop);

    if ((upscorereliable[candidate] && downscorereliable[candidate]) ||
        mustStop) {
      downNodeLb = downbound[candidate];
      upNodeLb = upbound[candidate];
      return candidate;
    }

    lp->setObjectiveLimit(getUpperLimit());

    HighsInt col = fracints[candidate].first;
    double fracval = fracints[candidate].second;
    double upval = std::ceil(fracval);
    double downval = std::floor(fracval);

    auto analyzeSolution = [&](double objdelta,
                               const std::vector<double>& sol) {
      size_t numChangedCols = localdom.getChangedCols().size();
      HighsInt domchgStackSize = localdom.getDomainChangeStack().size();
      const auto& domchgstack = localdom.getDomainChangeStack();

      for (HighsInt k = 0; k != numfrac; ++k) {
        if (fracints[k].first == col) continue;
        double otherfracval = fracints[k].second;
        double otherdownval = std::floor(fracints[k].second);
        double otherupval = std::ceil(fracints[k].second);
        if (sol[fracints[k].first] <= otherdownval + getFeasTol()) {
          if (localdom.col_upper_[fracints[k].first] >
              otherdownval + getFeasTol()) {
            localdom.changeBound(HighsBoundType::kUpper, fracints[k].first,
                                 otherdownval);
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(
                  getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
              localdom.backtrack();
              localdom.clearChangedCols(numChangedCols);
              continue;
            }
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(
                  getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
              localdom.backtrack();
              localdom.clearChangedCols(numChangedCols);
              continue;
            }

            HighsInt newStackSize = localdom.getDomainChangeStack().size();

            bool solutionValid = true;
            for (HighsInt j = domchgStackSize + 1; j < newStackSize; ++j) {
              if (domchgstack[j].boundtype == HighsBoundType::kLower) {
                if (domchgstack[j].boundval >
                    sol[domchgstack[j].column] + getFeasTol()) {
                  solutionValid = false;
                  break;
                }
              } else {
                if (domchgstack[j].boundval <
                    sol[domchgstack[j].column] - getFeasTol()) {
                  solutionValid = false;
                  break;
                }
              }
            }

            localdom.backtrack();
            localdom.clearChangedCols(numChangedCols);
            if (!solutionValid) continue;
          }

          if (objdelta <= getFeasTol()) {
            pseudocost.addObservation(fracints[k].first,
                                      otherdownval - otherfracval, objdelta);
            markBranchingVarDownReliableAtNode(fracints[k].first);
          }

          downscore[k] = std::min(downscore[k], objdelta);
        } else if (sol[fracints[k].first] >= otherupval - getFeasTol()) {
          if (localdom.col_lower_[fracints[k].first] <
              otherupval - getFeasTol()) {
            localdom.changeBound(HighsBoundType::kLower, fracints[k].first,
                                 otherupval);

            if (localdom.infeasible()) {
              localdom.conflictAnalysis(
                  getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
              localdom.backtrack();
              localdom.clearChangedCols(numChangedCols);
              continue;
            }
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(
                  getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
              localdom.backtrack();
              localdom.clearChangedCols(numChangedCols);
              continue;
            }

            HighsInt newStackSize = localdom.getDomainChangeStack().size();

            bool solutionValid = true;
            for (HighsInt j = domchgStackSize + 1; j < newStackSize; ++j) {
              if (domchgstack[j].boundtype == HighsBoundType::kLower) {
                if (domchgstack[j].boundval >
                    sol[domchgstack[j].column] + getFeasTol()) {
                  solutionValid = false;
                  break;
                }
              } else {
                if (domchgstack[j].boundval <
                    sol[domchgstack[j].column] - getFeasTol()) {
                  solutionValid = false;
                  break;
                }
              }
            }

            localdom.backtrack();
            localdom.clearChangedCols(numChangedCols);

            if (!solutionValid) continue;
          }

          if (objdelta <= getFeasTol()) {
            pseudocost.addObservation(fracints[k].first,
                                      otherupval - otherfracval, objdelta);
            markBranchingVarUpReliableAtNode(fracints[k].first);
          }

          upscore[k] = std::min(upscore[k], objdelta);
        }
      }
    };

    auto strongBranch = [&](bool upbranch) -> bool {
      int64_t inferences = -(int64_t)localdom.getDomainChangeStack().size() - 1;
      HighsBoundType boundtype =
          upbranch ? HighsBoundType::kLower : HighsBoundType::kUpper;
      double boundval = upbranch ? upval : downval;
      HighsDomainChange domchg{boundval, col, boundtype};

      bool orbitalFixing =
          nodestack.back().stabilizerOrbits && orbitsValidInChildNode(domchg);
      localdom.changeBound(domchg);
      localdom.propagate();

      if (!localdom.infeasible()) {
        if (orbitalFixing)
          nodestack.back().stabilizerOrbits->orbitalFixing(localdom);
        else
          getSymmetries().propagateOrbitopes(localdom);
      }

      inferences += localdom.getDomainChangeStack().size();
      if (localdom.infeasible()) {
        localdom.conflictAnalysis(getConflictPool(),
                                  mipworker.getGlobalDomain(), pseudocost);
        pseudocost.addCutoffObservation(col, upbranch);
        localdom.backtrack();
        localdom.clearChangedCols();

        if (upbranch) {
          branchDownwards(col, downval, fracval);
        } else {
          branchUpwards(col, upval, fracval);
        }
        nodestack[nodestack.size() - 2].opensubtrees = 0;
        nodestack[nodestack.size() - 2].skipDepthCount = 1;
        depthoffset -= 1;

        return true;
      }

      pseudocost.addInferenceObservation(col, inferences, upbranch);

      int64_t numiters = lp->getNumLpIterations();
      HighsLpRelaxation::Status status = playground.solveLp(localdom);
      numiters = lp->getNumLpIterations() - numiters;
      lpiterations += numiters;
      sblpiterations += numiters;

      if (lp->scaledOptimal(status)) {
        lp->performAging();

        double delta = upbranch ? upval - fracval : downval - fracval;
        bool integerfeasible;
        const std::vector<double>& sol = lp->getSolution().col_value;
        double solobj = checkSol(sol, integerfeasible);

        double objdelta = std::max(solobj - lp->getObjective(), 0.0);
        if (objdelta <= getEpsilon()) objdelta = 0.0;

        if (upbranch) {
          upscore[candidate] = objdelta;
          upscorereliable[candidate] = true;
          markBranchingVarUpReliableAtNode(col);
        } else {
          downscore[candidate] = objdelta;
          downscorereliable[candidate] = true;
          markBranchingVarDownReliableAtNode(col);
        }
        pseudocost.addObservation(col, delta, objdelta);
        analyzeSolution(objdelta, sol);

        if (lp->unscaledPrimalFeasible(status) && integerfeasible) {
          double cutoffbnd = getCutoffBound();
          addIncumbent(lp->getLpSolver().getSolution().col_value, solobj,
                       inheuristic ? kSolutionSourceHeuristic
                                   : kSolutionSourceBranching);

          if (getUpperLimit() < cutoffbnd)
            lp->setObjectiveLimit(getUpperLimit());
        }

        if (lp->unscaledDualFeasible(status)) {
          if (upbranch) {
            upbound[candidate] = solobj;
          } else {
            downbound[candidate] = solobj;
          }
          if (solobj > mipworker.getOptimalityLimit()) {
            addBoundExceedingConflict();

            bool pruned = solobj > getCutoffBound();
            if (pruned)
              mipsolver.mipdata_->debugSolution.nodePruned(localdom,
                                                           getCutoffBound());

            localdom.backtrack();
            lp->flushDomain(localdom);

            if (upbranch) {
              branchDownwards(col, downval, fracval);
            } else {
              branchUpwards(col, upval, fracval);
            }
            nodestack[nodestack.size() - 2].opensubtrees = pruned ? 0 : 1;
            nodestack[nodestack.size() - 2].other_child_lb = solobj;
            nodestack[nodestack.size() - 2].skipDepthCount = 1;
            depthoffset -= 1;

            return true;
          }
        } else if (solobj > getCutoffBound()) {
          addBoundExceedingConflict();
          localdom.propagate();
          bool infeas = localdom.infeasible();
          if (infeas) {
            localdom.backtrack();
            lp->flushDomain(localdom);

            if (upbranch) {
              branchDownwards(col, downval, fracval);
            } else {
              branchUpwards(col, upval, fracval);
            }
            nodestack[nodestack.size() - 2].opensubtrees = 0;
            nodestack[nodestack.size() - 2].skipDepthCount = 1;
            depthoffset -= 1;

            return true;
          }
        }
      } else if (status == HighsLpRelaxation::Status::kInfeasible) {
        mipsolver.mipdata_->debugSolution.nodePruned(localdom,
                                                     getCutoffBound());
        addInfeasibleConflict();
        pseudocost.addCutoffObservation(col, upbranch);
        localdom.backtrack();
        lp->flushDomain(localdom);

        if (upbranch) {
          branchDownwards(col, downval, fracval);
        } else {
          branchUpwards(col, upval, fracval);
        }
        nodestack[nodestack.size() - 2].opensubtrees = 0;
        nodestack[nodestack.size() - 2].skipDepthCount = 1;
        depthoffset -= 1;

        return true;
      } else {
        // printf("todo2\n");
        // in case of an LP error we set the score of this variable to zero to
        // avoid choosing it as branching candidate if possible
        downscore[candidate] = 0.0;
        upscore[candidate] = 0.0;
        downscorereliable[candidate] = true;
        upscorereliable[candidate] = true;
        markBranchingVarUpReliableAtNode(col);
        markBranchingVarDownReliableAtNode(col);
      }

      localdom.backtrack();
      lp->flushDomain(localdom);
      return false;
    };

    if (!downscorereliable[candidate] &&
        (upscorereliable[candidate] ||
         std::make_pair(downscore[candidate],
                        pseudocost.getAvgInferencesDown(col)) >=
             std::make_pair(upscore[candidate],
                            pseudocost.getAvgInferencesUp(col)))) {
      // evaluate down branch
      // if (!mipsolver.submip)
      //   printf("down eval col=%d fracval=%g\n", col, fracval);
      if (strongBranch(false)) return -1;
    } else {
      // if (!mipsolver.submip)
      //  printf("up eval col=%d fracval=%g\n", col, fracval);
      // evaluate up branch
      if (strongBranch(true)) return -1;
    }
  }
}

const HighsSearch::NodeData* HighsSearch::getParentNodeData() const {
  if (nodestack.size() <= 1) return nullptr;

  return &nodestack[nodestack.size() - 2];
}

void HighsSearch::stashNodeToProcessed(double lb, double estimate,
                                       HighsInt depth) {
  std::vector<HighsInt> branchPositions;
  auto domchgStack = localdom.getReducedDomainChangeStack(branchPositions);
  mipworker.processedNodes.emplace_back(
      std::piecewise_construct,
      std::forward_as_tuple(std::move(domchgStack), std::move(branchPositions),
                            lb, estimate, depth),
      std::forward_as_tuple(countTreeWeight));
}

void HighsSearch::stashCurrentNode() {
  auto oldchangedcols = localdom.getChangedCols().size();
  bool prune = nodestack.back().lower_bound > getCutoffBound();
  if (!prune) {
    localdom.propagate();
    localdom.clearChangedCols(oldchangedcols);
    prune = localdom.infeasible();
    if (prune)
      localdom.conflictAnalysis(getConflictPool(), mipworker.getGlobalDomain(),
                                pseudocost);
  }
  if (!prune) {
    stashNodeToProcessed(std::max(nodestack.back().lower_bound,
                                  localdom.getObjectiveLowerBound()),
                         nodestack.back().estimate, getCurrentDepth());
  } else {
    mipsolver.mipdata_->debugSolution.nodePruned(localdom, getCutoffBound());
    if (countTreeWeight) treeweight += std::ldexp(1.0, 1 - getCurrentDepth());
  }
  nodestack.back().opensubtrees = 0;
}

void HighsSearch::stashOpenNodes() {
  if (nodestack.empty()) return;

  // get the basis of the node highest up in the tree
  std::shared_ptr<const HighsBasis> basis;
  for (NodeData& nodeData : nodestack) {
    if (nodeData.nodeBasis) {
      basis = std::move(nodeData.nodeBasis);
      break;
    }
  }

  if (nodestack.back().opensubtrees == 0) backtrack(false);

  while (!nodestack.empty()) {
    auto oldchangedcols = localdom.getChangedCols().size();
    bool prune = nodestack.back().lower_bound > getCutoffBound();
    if (!prune) {
      localdom.propagate();
      localdom.clearChangedCols(oldchangedcols);
      prune = localdom.infeasible();
      if (prune)
        localdom.conflictAnalysis(getConflictPool(),
                                  mipworker.getGlobalDomain(), pseudocost);
    }
    if (!prune) {
      stashNodeToProcessed(std::max(nodestack.back().lower_bound,
                                    localdom.getObjectiveLowerBound()),
                           nodestack.back().estimate, getCurrentDepth());
    } else {
      mipsolver.mipdata_->debugSolution.nodePruned(localdom, getCutoffBound());
      if (countTreeWeight) treeweight += std::ldexp(1.0, 1 - getCurrentDepth());
    }
    nodestack.back().opensubtrees = 0;
    backtrack(false);
  }

  lp->flushDomain(localdom);
  if (basis) {
    if ((HighsInt)basis->row_status.size() == lp->numRows())
      lp->setStoredBasis(std::move(basis));
    lp->recoverBasis();
  }
}

void HighsSearch::flushStatistics(HighsMipSolver& mipsolver) {
  mipsolver.mipdata_->num_nodes += nnodes;
  mipsolver.mipdata_->num_leaves += nleaves;
  mipsolver.mipdata_->pruned_treeweight += treeweight;
  mipsolver.mipdata_->total_lp_iterations += lpiterations;
  mipsolver.mipdata_->heuristic_lp_iterations += heurlpiterations;
  mipsolver.mipdata_->sb_lp_iterations += sblpiterations;
  resetStatistics();
}

void HighsSearch::resetStatistics() {
  nnodes = 0;
  nleaves = 0;
  treeweight = 0;
  lpiterations = 0;
  heurlpiterations = 0;
  sblpiterations = 0;
}

int64_t HighsSearch::getHeuristicLpIterations() const {
  return heurlpiterations + mipsolver.mipdata_->heuristic_lp_iterations;
}

int64_t HighsSearch::getTotalLpIterations() const {
  return lpiterations + mipsolver.mipdata_->total_lp_iterations;
}

int64_t HighsSearch::getLocalLpIterations() const { return lpiterations; }

int64_t& HighsSearch::getLocalNodes() { return nnodes; }

int64_t& HighsSearch::getLocalLeaves() { return nleaves; }

int64_t HighsSearch::getStrongBranchingLpIterations() const {
  return sblpiterations + mipsolver.mipdata_->sb_lp_iterations;
}

void HighsSearch::resetLocalDomain() {
  this->lp->resetToGlobalDomain(getDomain());
  localdom = getDomain();

#ifndef NDEBUG
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    assert(lp->getLpSolver().getLp().col_lower_[i] == localdom.col_lower_[i] ||
           mipsolver.isColContinuous(i));
    assert(lp->getLpSolver().getLp().col_upper_[i] == localdom.col_upper_[i] ||
           mipsolver.isColContinuous(i));
  }
#endif
}

void HighsSearch::installNode(HighsNodeQueue::OpenNode&& node) {
  localdom.setDomainChangeStack(node.domchgstack, node.branchings);
  bool globalSymmetriesValid = true;
  if (mipsolver.mipdata_->globalOrbits) {
    // if global orbits have been computed we check whether they are still valid
    // in this node
    const auto& domchgstack = localdom.getDomainChangeStack();
    for (HighsInt i : localdom.getBranchingPositions()) {
      HighsInt col = domchgstack[i].column;
      if (getSymmetries().columnPosition[col] == -1) continue;

      if (!getDomain().isBinary(col) ||
          (domchgstack[i].boundtype == HighsBoundType::kLower &&
           domchgstack[i].boundval == 1.0)) {
        globalSymmetriesValid = false;
        break;
      }
    }
  }
  nodestack.emplace_back(
      node.lower_bound, node.estimate, nullptr,
      globalSymmetriesValid ? mipsolver.mipdata_->globalOrbits : nullptr);
  subrootsol.clear();
  depthoffset = node.depth - 1;
}

HighsSearch::NodeResult HighsSearch::evaluateNode() {
  assert(!nodestack.empty());
  NodeData& currnode = nodestack.back();
  const NodeData* parent = getParentNodeData();

  const auto& domchgstack = localdom.getDomainChangeStack();

  if (!inheuristic && currnode.lower_bound > mipworker.getOptimalityLimit())
    return NodeResult::kSubOptimal;

  localdom.propagate();

  if (!inheuristic && !localdom.infeasible()) {
    if (getSymmetries().numPerms > 0 && !currnode.stabilizerOrbits &&
        (parent == nullptr || !parent->stabilizerOrbits ||
         !parent->stabilizerOrbits->orbitCols.empty())) {
      currnode.stabilizerOrbits = getSymmetries().computeStabilizerOrbits(
          localdom, stabilizerOrbitWorkspace);
    }

    if (currnode.stabilizerOrbits)
      currnode.stabilizerOrbits->orbitalFixing(localdom);
    else
      getSymmetries().propagateOrbitopes(localdom);
  }
  if (parent != nullptr) {
    int64_t inferences = domchgstack.size() - (currnode.domgchgStackPos + 1);

    pseudocost.addInferenceObservation(
        parent->branchingdecision.column, inferences,
        parent->branchingdecision.boundtype == HighsBoundType::kLower);
  }

  NodeResult result = NodeResult::kOpen;

  if (localdom.infeasible()) {
    result = NodeResult::kDomainInfeasible;
    localdom.clearChangedCols();
    if (parent != nullptr && parent->lp_objective != -kHighsInf &&
        parent->branching_point != parent->branchingdecision.boundval) {
      bool upbranch =
          parent->branchingdecision.boundtype == HighsBoundType::kLower;
      pseudocost.addCutoffObservation(parent->branchingdecision.column,
                                      upbranch);
    }

    localdom.conflictAnalysis(getConflictPool(), mipworker.getGlobalDomain(),
                              pseudocost);
  } else {
    lp->flushDomain(localdom);
    lp->setObjectiveLimit(getUpperLimit());

#ifndef NDEBUG
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      assert(lp->getLpSolver().getLp().col_lower_[i] ==
                 localdom.col_lower_[i] ||
             mipsolver.isColContinuous(i));
      assert(lp->getLpSolver().getLp().col_upper_[i] ==
                 localdom.col_upper_[i] ||
             mipsolver.isColContinuous(i));
    }
#endif
    int64_t oldnumiters = lp->getNumLpIterations();
    HighsLpRelaxation::Status status = lp->resolveLp(&localdom);
    lpiterations += lp->getNumLpIterations() - oldnumiters;

    currnode.lower_bound =
        std::max(localdom.getObjectiveLowerBound(), currnode.lower_bound);

    if (localdom.infeasible()) {
      result = NodeResult::kDomainInfeasible;
      localdom.clearChangedCols();
      if (parent != nullptr && parent->lp_objective != -kHighsInf &&
          parent->branching_point != parent->branchingdecision.boundval) {
        bool upbranch =
            parent->branchingdecision.boundtype == HighsBoundType::kLower;
        pseudocost.addCutoffObservation(parent->branchingdecision.column,
                                        upbranch);
      }

      localdom.conflictAnalysis(getConflictPool(), mipworker.getGlobalDomain(),
                                pseudocost);
    } else if (lp->scaledOptimal(status)) {
      lp->storeBasis();
      lp->performAging();

      currnode.nodeBasis = lp->getStoredBasis();
      currnode.estimate = lp->computeBestEstimate(pseudocost);
      currnode.lp_objective = lp->getObjective();

      if (parent != nullptr && parent->lp_objective != -kHighsInf &&
          parent->branching_point != parent->branchingdecision.boundval) {
        double delta =
            parent->branchingdecision.boundval - parent->branching_point;
        double objdelta =
            std::max(0.0, currnode.lp_objective - parent->lp_objective);

        pseudocost.addObservation(parent->branchingdecision.column, delta,
                                  objdelta);
      }

      if (lp->unscaledPrimalFeasible(status)) {
        if (lp->getFractionalIntegers().empty()) {
          double cutoffbnd = getCutoffBound();
          addIncumbent(lp->getLpSolver().getSolution().col_value,
                       lp->getObjective(),
                       inheuristic ? kSolutionSourceHeuristic
                                   : kSolutionSourceEvaluateNode);
          if (getUpperLimit() < cutoffbnd)
            lp->setObjectiveLimit(getUpperLimit());

          if (lp->unscaledDualFeasible(status)) {
            addBoundExceedingConflict();
            result = NodeResult::kBoundExceeding;
          }
        }
      }

      if (result == NodeResult::kOpen) {
        if (lp->unscaledDualFeasible(status)) {
          currnode.lower_bound =
              std::max(currnode.lp_objective, currnode.lower_bound);

          if (currnode.lower_bound > getCutoffBound()) {
            result = NodeResult::kBoundExceeding;
            addBoundExceedingConflict();
          } else if (getUpperLimit() != kHighsInf) {
            if (!inheuristic) {
              double gap = getUpperLimit() - lp->getObjective();
              lp->computeBasicDegenerateDuals(
                  gap + std::max(10 * getFeasTol(), getEpsilon() * gap),
                  localdom, getDomain(), getConflictPool(),
                  mipworker.getPseudocost(), true);
            }
            HighsRedcostFixing::propagateRedCost(
                mipsolver, localdom, mipworker.getGlobalDomain(), *lp,
                getConflictPool(), mipworker.getPseudocost(), getUpperLimit());
            localdom.propagate();
            if (localdom.infeasible()) {
              result = NodeResult::kDomainInfeasible;
              localdom.clearChangedCols();
              if (parent != nullptr && parent->lp_objective != -kHighsInf &&
                  parent->branching_point !=
                      parent->branchingdecision.boundval) {
                bool upbranch = parent->branchingdecision.boundtype ==
                                HighsBoundType::kLower;
                pseudocost.addCutoffObservation(
                    parent->branchingdecision.column, upbranch);
              }

              localdom.conflictAnalysis(
                  getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
            } else if (!localdom.getChangedCols().empty()) {
              return evaluateNode();
            }
          } else {
            if (!inheuristic) {
              lp->computeBasicDegenerateDuals(kHighsInf, localdom, getDomain(),
                                              getConflictPool(),
                                              mipworker.getPseudocost(), true);
              localdom.propagate();
              if (localdom.infeasible()) {
                result = NodeResult::kDomainInfeasible;
                localdom.clearChangedCols();
                if (parent != nullptr && parent->lp_objective != -kHighsInf &&
                    parent->branching_point !=
                        parent->branchingdecision.boundval) {
                  bool upbranch = parent->branchingdecision.boundtype ==
                                  HighsBoundType::kLower;
                  pseudocost.addCutoffObservation(
                      parent->branchingdecision.column, upbranch);
                }

                localdom.conflictAnalysis(
                    getConflictPool(), mipworker.getGlobalDomain(), pseudocost);
              } else if (!localdom.getChangedCols().empty()) {
                return evaluateNode();
              }
            }
          }
        } else if (lp->getObjective() > getCutoffBound()) {
          // the LP is not solved to dual feasibility due to scaling/numerics
          // therefore we compute a conflict constraint as if the LP was bound
          // exceeding and propagate the local domain again. The lp relaxation
          // class will take care to consider the dual multipliers with an
          // increased zero tolerance due to the dual infeasibility when
          // computing the proof conBoundExceedingstraint.
          addBoundExceedingConflict();
          localdom.propagate();
          if (localdom.infeasible()) {
            result = NodeResult::kBoundExceeding;
          }
        }
      }
    } else if (status == HighsLpRelaxation::Status::kInfeasible) {
      if (lp->getLpSolver().getModelStatus() ==
          HighsModelStatus::kObjectiveBound)
        result = NodeResult::kBoundExceeding;
      else
        result = NodeResult::kLpInfeasible;
      addInfeasibleConflict();
      if (parent != nullptr && parent->lp_objective != -kHighsInf &&
          parent->branching_point != parent->branchingdecision.boundval) {
        bool upbranch =
            parent->branchingdecision.boundtype == HighsBoundType::kLower;
        pseudocost.addCutoffObservation(parent->branchingdecision.column,
                                        upbranch);
      }
    }
  }

  if (result != NodeResult::kOpen) {
    mipsolver.mipdata_->debugSolution.nodePruned(localdom, getCutoffBound());
    treeweight += std::ldexp(1.0, 1 - getCurrentDepth());
    currnode.opensubtrees = 0;
  } else if (!inheuristic) {
    if (currnode.lower_bound > mipworker.getOptimalityLimit()) {
      result = NodeResult::kSubOptimal;
      addBoundExceedingConflict();
    }
  }

  return result;
}

HighsSearch::NodeResult HighsSearch::branch() {
  assert(localdom.getChangedCols().empty());

  assert(nodestack.back().opensubtrees == 2);
  nodestack.back().branchingdecision.column = -1;
  inbranching = true;

  HighsInt minrel = pseudocost.getMinReliable();
  double childLb = getCurrentLowerBound();
  NodeResult result = NodeResult::kOpen;
  while (nodestack.back().opensubtrees == 2 &&
         lp->scaledOptimal(lp->getStatus()) &&
         !lp->getFractionalIntegers().empty()) {
    int64_t sbmaxiters = 0;
    if (minrel > 0) {
      int64_t sbiters = getStrongBranchingLpIterations();
      sbmaxiters =
          100000 + ((getTotalLpIterations() - getHeuristicLpIterations() -
                     getStrongBranchingLpIterations()) >>
                    1);
      if (mipsolver.mipdata_->numRestarts <= 2) sbmaxiters = sbmaxiters >> 2;
      if (sbiters > sbmaxiters) {
        pseudocost.setMinReliable(0);
      } else if (sbiters > (sbmaxiters >> 1)) {
        double reductionratio = (sbiters - (sbmaxiters >> 1)) /
                                (double)(sbmaxiters - (sbmaxiters >> 1));

        HighsInt minrelreduced = int(minrel - reductionratio * (minrel - 1));
        pseudocost.setMinReliable(std::min(minrel, minrelreduced));
      }
    }

    double degeneracyFac = lp->computeLPDegneracy(localdom);
    pseudocost.setDegeneracyFactor(degeneracyFac);
    if (degeneracyFac >= 10.0) pseudocost.setMinReliable(0);
    // if (!mipsolver.submip)
    //  printf("selecting branching cand with minrel=%d\n",
    //         pseudocost.getMinReliable());
    double downNodeLb = getCurrentLowerBound();
    double upNodeLb = getCurrentLowerBound();
    HighsInt branchcand =
        selectBranchingCandidate(sbmaxiters, downNodeLb, upNodeLb);
    // if (!mipsolver.submip)
    //   printf("branching cand returned as %d\n", branchcand);
    NodeData& currnode = nodestack.back();
    childLb = currnode.lower_bound;
    if (branchcand != -1) {
      auto branching = lp->getFractionalIntegers()[branchcand];
      currnode.branchingdecision.column = branching.first;
      currnode.branching_point = branching.second;

      HighsInt col = branching.first;

      switch (childselrule) {
        case ChildSelectionRule::kUp:
          currnode.branchingdecision.boundtype = HighsBoundType::kLower;
          currnode.branchingdecision.boundval =
              std::ceil(currnode.branching_point);
          currnode.other_child_lb = downNodeLb;
          childLb = upNodeLb;
          break;
        case ChildSelectionRule::kDown:
          currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
          currnode.branchingdecision.boundval =
              std::floor(currnode.branching_point);
          currnode.other_child_lb = upNodeLb;
          childLb = downNodeLb;
          break;
        case ChildSelectionRule::kRootSol: {
          double downPrio = pseudocost.getAvgInferencesDown(col) + getEpsilon();
          double upPrio = pseudocost.getAvgInferencesUp(col) + getEpsilon();
          double downVal = std::floor(currnode.branching_point);
          double upVal = std::ceil(currnode.branching_point);
          if (!subrootsol.empty()) {
            double rootsol = subrootsol[col];
            if (rootsol < downVal)
              rootsol = downVal;
            else if (rootsol > upVal)
              rootsol = upVal;

            upPrio *= (1.0 + (currnode.branching_point - rootsol));
            downPrio *= (1.0 + (rootsol - currnode.branching_point));

          } else {
            if (currnode.lp_objective != -kHighsInf)
              subrootsol = lp->getSolution().col_value;
            if (!getRootLpSol().empty()) {
              double rootsol = getRootLpSol()[col];
              if (rootsol < downVal)
                rootsol = downVal;
              else if (rootsol > upVal)
                rootsol = upVal;

              upPrio *= (1.0 + (currnode.branching_point - rootsol));
              downPrio *= (1.0 + (rootsol - currnode.branching_point));
            }
          }
          if (upPrio + getEpsilon() >= downPrio) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval = upVal;
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval = downVal;
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          }
          break;
        }
        case ChildSelectionRule::kObj:
          if (mipsolver.colCost(col) >= 0) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          }
          break;
        case ChildSelectionRule::kRandom:
          if (random.bit()) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          }
          break;
        case ChildSelectionRule::kBestCost: {
          if (pseudocost.getPseudocostUp(col, currnode.branching_point,
                                         getFeasTol()) >
              pseudocost.getPseudocostDown(col, currnode.branching_point,
                                           getFeasTol())) {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          }
          break;
        }
        case ChildSelectionRule::kWorstCost:
          if (pseudocost.getPseudocostUp(col, currnode.branching_point) >=
              pseudocost.getPseudocostDown(col, currnode.branching_point)) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          }
          break;
        case ChildSelectionRule::kDisjunction: {
          int64_t numnodesup;
          int64_t numnodesdown;
          numnodesup = getNodeQueue().numNodesUp(col);
          numnodesdown = getNodeQueue().numNodesDown(col);
          if (numnodesup > numnodesdown) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else if (numnodesdown > numnodesup) {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          } else {
            if (mipsolver.colCost(col) >= 0) {
              currnode.branchingdecision.boundtype = HighsBoundType::kLower;
              currnode.branchingdecision.boundval =
                  std::ceil(currnode.branching_point);
              currnode.other_child_lb = downNodeLb;
              childLb = upNodeLb;
            } else {
              currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
              currnode.branchingdecision.boundval =
                  std::floor(currnode.branching_point);
              currnode.other_child_lb = upNodeLb;
              childLb = downNodeLb;
            }
          }
          break;
        }
        case ChildSelectionRule::kHybridInferenceCost: {
          double upVal = std::ceil(currnode.branching_point);
          double downVal = std::floor(currnode.branching_point);
          double upScore = (1 + pseudocost.getAvgInferencesUp(col)) /
                           pseudocost.getPseudocostUp(
                               col, currnode.branching_point, getFeasTol());
          double downScore = (1 + pseudocost.getAvgInferencesDown(col)) /
                             pseudocost.getPseudocostDown(
                                 col, currnode.branching_point, getFeasTol());

          if (upScore >= downScore) {
            currnode.branchingdecision.boundtype = HighsBoundType::kLower;
            currnode.branchingdecision.boundval = upVal;
            currnode.other_child_lb = downNodeLb;
            childLb = upNodeLb;
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
            currnode.branchingdecision.boundval = downVal;
            currnode.other_child_lb = upNodeLb;
            childLb = downNodeLb;
          }
        }
      }
      result = NodeResult::kBranched;
      break;
    }

    assert(!localdom.getChangedCols().empty());
    result = evaluateNode();
    if (result == NodeResult::kSubOptimal) break;
  }
  inbranching = false;
  NodeData& currnode = nodestack.back();
  pseudocost.setMinReliable(minrel);
  pseudocost.setDegeneracyFactor(1.0);

  assert(currnode.opensubtrees == 2 || currnode.opensubtrees == 0);

  if (currnode.opensubtrees != 2 || result == NodeResult::kSubOptimal)
    return result;

  if (currnode.branchingdecision.column == -1) {
    double bestscore = -1.0;
    // solution branching failed, so choose any integer variable to branch
    // on in case we have a different solution status could happen due to a
    // fail in the LP solution process
    pseudocost.setDegeneracyFactor(1e6);

    for (HighsInt i : getIntegralCols()) {
      if (localdom.col_upper_[i] - localdom.col_lower_[i] < 0.5) continue;

      double fracval;
      if (localdom.col_lower_[i] != -kHighsInf &&
          localdom.col_upper_[i] != kHighsInf)
        fracval = std::floor(0.5 * (localdom.col_lower_[i] +
                                    localdom.col_upper_[i] + 0.5)) +
                  0.5;
      else if (localdom.col_lower_[i] != -kHighsInf)
        fracval = localdom.col_lower_[i] + 0.5;
      else if (localdom.col_upper_[i] != kHighsInf)
        fracval = localdom.col_upper_[i] - 0.5;
      else
        fracval = 0.5;

      double score = pseudocost.getScore(i, fracval);
      assert(score >= 0.0);

      if (score > bestscore) {
        bestscore = score;
        bool branchUpwards;
        double cost = lp->unscaledDualFeasible(lp->getStatus())
                          ? lp->getSolution().col_dual[i]
                          : mipsolver.colCost(i);
        if (std::fabs(cost) > getFeasTol() && getCutoffBound() < kHighsInf) {
          // branch in direction of worsening cost first in case the column has
          // cost and we do have an upper bound
          branchUpwards = cost > 0;
        } else if (pseudocost.getAvgInferencesUp(i) >
                   pseudocost.getAvgInferencesDown(i) + getFeasTol()) {
          // column does not have (reduced) cost above tolerance so branch in
          // direction of more inferences
          branchUpwards = true;
        } else if (pseudocost.getAvgInferencesUp(i) <
                   pseudocost.getAvgInferencesDown(i) - getFeasTol()) {
          branchUpwards = false;
        } else {
          // number of inferences give a tie, so we branch in the direction that
          // does have a less recent domain change to avoid branching the same
          // integer column into the same direction over and over
          HighsInt colLowerPos;
          HighsInt colUpperPos;
          localdom.getColLowerPos(i, localdom.getNumDomainChanges(),
                                  colLowerPos);
          localdom.getColUpperPos(i, localdom.getNumDomainChanges(),
                                  colUpperPos);
          branchUpwards = colLowerPos <= colUpperPos;
        }
        if (branchUpwards) {
          double upval = std::ceil(fracval);
          currnode.branching_point = upval;
          currnode.branchingdecision.boundtype = HighsBoundType::kLower;
          currnode.branchingdecision.column = i;
          currnode.branchingdecision.boundval = upval;
        } else {
          double downval = std::floor(fracval);
          currnode.branching_point = downval;
          currnode.branchingdecision.boundtype = HighsBoundType::kUpper;
          currnode.branchingdecision.column = i;
          currnode.branchingdecision.boundval = downval;
        }
      }
    }

    pseudocost.setDegeneracyFactor(1);
  }

  if (currnode.branchingdecision.column == -1) {
    if (lp->getStatus() == HighsLpRelaxation::Status::kOptimal) {
      // if the LP was solved to optimality and all columns are fixed, then this
      // particular assignment is not feasible or has a worse objective in the
      // original space, otherwise the node would not be open. Hence we prune
      // this particular assignment
      currnode.opensubtrees = 0;
      result = NodeResult::kLpInfeasible;
      return result;
    }
    lp->setIterationLimit();

    // create a fresh LP only with model rows since all integer columns are
    // fixed, the cutting planes are not required and the LP could not be solved
    // so we want to make it as easy as possible
    //
    // LP relaxation instantiation
    HighsLpRelaxation lpCopy(mipsolver);
    lpCopy.setProfiling(mipsolver.profiling_);
    lpCopy.loadModel();
    lpCopy.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                          localdom.col_lower_.data(),
                                          localdom.col_upper_.data());
    // temporarily use the fresh LP for the HighsSearch class
    HighsLpRelaxation* tmpLp = &lpCopy;
    std::swap(tmpLp, lp);

    // reevaluate the node with LP presolve enabled
    lp->getLpSolver().setOptionValue("presolve", kHighsOnString);
    result = evaluateNode();

    if (result == NodeResult::kOpen) {
      // LP still not solved, reevaluate with primal simplex
      lp->getLpSolver().clearSolver();
      lp->getLpSolver().setOptionValue("simplex_strategy",
                                       kSimplexStrategyPrimal);
      result = evaluateNode();
      lp->getLpSolver().setOptionValue("simplex_strategy",
                                       kSimplexStrategyDual);
      if (result == NodeResult::kOpen) {
        // LP still not solved, reevaluate with IPM instead of simplex
        lp->getLpSolver().clearSolver();
        lp->getLpSolver().setOptionValue("solver", "ipm");
        result = evaluateNode();

        if (result == NodeResult::kOpen) {
          highsLogUser(mipsolver.options_mip_->log_options,
                       HighsLogType::kWarning,
                       "Failed to solve node with all integer columns "
                       "fixed. Declaring node infeasible.\n");
          // LP still not solved, give up and declare as infeasible
          currnode.opensubtrees = 0;
          result = NodeResult::kLpInfeasible;
        }
      }
    }

    // restore old lp relaxation
    std::swap(tmpLp, lp);

    return result;
  }

  // finally open a new node with the branching decision added
  // and remember that we have one open subtree left
  HighsInt domchgPos = localdom.getDomainChangeStack().size();

  bool passStabilizerToChildNode =
      orbitsValidInChildNode(currnode.branchingdecision);
  localdom.changeBound(currnode.branchingdecision);
  currnode.opensubtrees = 1;
  nodestack.emplace_back(
      std::max(childLb, currnode.lower_bound), currnode.estimate,
      currnode.nodeBasis,
      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);
  nodestack.back().domgchgStackPos = domchgPos;

  return NodeResult::kBranched;
}

bool HighsSearch::backtrack(bool recoverBasis) {
  if (nodestack.empty()) return false;
  assert(!nodestack.empty());
  assert(nodestack.back().opensubtrees == 0);
  while (true) {
    while (nodestack.back().opensubtrees == 0) {
      countTreeWeight = true;
      depthoffset += nodestack.back().skipDepthCount;
      if (nodestack.size() == 1) {
        if (recoverBasis && nodestack.back().nodeBasis)
          lp->setStoredBasis(std::move(nodestack.back().nodeBasis));
        nodestack.pop_back();
        localdom.backtrackToGlobal();
        lp->flushDomain(localdom);
        if (recoverBasis) lp->recoverBasis();
        return false;
      }

      nodestack.pop_back();
#ifndef NDEBUG
      HighsDomainChange branchchg =
#endif
          localdom.backtrack();

      if (nodestack.back().opensubtrees != 0) {
        countTreeWeight = nodestack.back().skipDepthCount == 0;
        // repropagate the node, as it may have become infeasible due to
        // conflicts
        HighsInt oldNumDomchgs = localdom.getNumDomainChanges();
        size_t oldNumChangedCols = localdom.getChangedCols().size();
        localdom.propagate();
        if (!localdom.infeasible() &&
            oldNumDomchgs != localdom.getNumDomainChanges()) {
          if (nodestack.back().stabilizerOrbits)
            nodestack.back().stabilizerOrbits->orbitalFixing(localdom);
          else
            getSymmetries().propagateOrbitopes(localdom);
        }
        if (localdom.infeasible()) {
          localdom.clearChangedCols(oldNumChangedCols);
          if (countTreeWeight)
            treeweight += std::ldexp(1.0, -getCurrentDepth());
          nodestack.back().opensubtrees = 0;
        }
      }

      assert(
          (branchchg.boundtype == HighsBoundType::kLower &&
           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||
          (branchchg.boundtype == HighsBoundType::kUpper &&
           branchchg.boundval <= nodestack.back().branchingdecision.boundval));
      assert(branchchg.boundtype ==
             nodestack.back().branchingdecision.boundtype);
      assert(branchchg.column == nodestack.back().branchingdecision.column);
    }

    NodeData& currnode = nodestack.back();

    assert(currnode.opensubtrees == 1);
    currnode.opensubtrees = 0;
    flipBranchingDecision(currnode);
    HighsInt domchgPos = localdom.getDomainChangeStack().size();

    size_t numChangedCols = localdom.getChangedCols().size();
    bool passStabilizerToChildNode =
        orbitsValidInChildNode(currnode.branchingdecision);
    localdom.changeBound(currnode.branchingdecision);
    double nodelb = std::max(currnode.lower_bound, currnode.other_child_lb);
    bool prune = nodelb > getCutoffBound() || localdom.infeasible();
    if (!prune) {
      localdom.propagate();
      prune = localdom.infeasible();
      if (prune)
        localdom.conflictAnalysis(getConflictPool(),
                                  mipworker.getGlobalDomain(), pseudocost);
    }
    if (!prune) {
      getSymmetries().propagateOrbitopes(localdom);
      prune = localdom.infeasible();
    }
    if (!prune && passStabilizerToChildNode && currnode.stabilizerOrbits) {
      currnode.stabilizerOrbits->orbitalFixing(localdom);
      prune = localdom.infeasible();
    }
    if (prune) {
      localdom.backtrack();
      localdom.clearChangedCols(numChangedCols);
      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());
      continue;
    }
    nodestack.emplace_back(
        nodelb, currnode.estimate, currnode.nodeBasis,
        passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);

    lp->flushDomain(localdom);
    nodestack.back().domgchgStackPos = domchgPos;
    break;
  }

  if (recoverBasis && nodestack.back().nodeBasis) {
    lp->setStoredBasis(nodestack.back().nodeBasis);
    lp->recoverBasis();
  }

  return true;
}

bool HighsSearch::backtrackPlunge() {
  const std::vector<HighsDomainChange>& domchgstack =
      localdom.getDomainChangeStack();

  if (nodestack.empty()) return false;
  assert(!nodestack.empty());
  assert(nodestack.back().opensubtrees == 0);

  while (true) {
    while (nodestack.back().opensubtrees == 0) {
      countTreeWeight = true;
      depthoffset += nodestack.back().skipDepthCount;

      if (nodestack.size() == 1) {
        if (nodestack.back().nodeBasis)
          lp->setStoredBasis(std::move(nodestack.back().nodeBasis));
        nodestack.pop_back();
        localdom.backtrackToGlobal();
        lp->flushDomain(localdom);
        lp->recoverBasis();
        return false;
      }

      nodestack.pop_back();
#ifndef NDEBUG
      HighsDomainChange branchchg =
#endif
          localdom.backtrack();

      if (nodestack.back().opensubtrees != 0) {
        countTreeWeight = nodestack.back().skipDepthCount == 0;
        // repropagate the node, as it may have become infeasible due to
        // conflicts
        HighsInt oldNumDomchgs = localdom.getNumDomainChanges();
        HighsInt oldNumChangedCols = localdom.getChangedCols().size();
        localdom.propagate();
        if (!localdom.infeasible() &&
            oldNumDomchgs != localdom.getNumDomainChanges()) {
          if (nodestack.back().stabilizerOrbits)
            nodestack.back().stabilizerOrbits->orbitalFixing(localdom);
          else
            getSymmetries().propagateOrbitopes(localdom);
        }
        if (localdom.infeasible()) {
          localdom.clearChangedCols(oldNumChangedCols);
          if (countTreeWeight)
            treeweight += std::ldexp(1.0, -getCurrentDepth());
          nodestack.back().opensubtrees = 0;
        }
      }

      assert(
          (branchchg.boundtype == HighsBoundType::kLower &&
           branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||
          (branchchg.boundtype == HighsBoundType::kUpper &&
           branchchg.boundval <= nodestack.back().branchingdecision.boundval));
      assert(branchchg.boundtype ==
             nodestack.back().branchingdecision.boundtype);
      assert(branchchg.column == nodestack.back().branchingdecision.column);
    }

    NodeData& currnode = nodestack.back();

    assert(currnode.opensubtrees == 1);
    currnode.opensubtrees = 0;
    flipBranchingDecision(currnode);

    HighsInt domchgPos = domchgstack.size();
    size_t numChangedCols = localdom.getChangedCols().size();
    bool passStabilizerToChildNode =
        orbitsValidInChildNode(currnode.branchingdecision);
    localdom.changeBound(currnode.branchingdecision);
    double nodelb = std::max(currnode.lower_bound, currnode.other_child_lb);
    bool prune = nodelb > getCutoffBound() || localdom.infeasible();
    if (!prune) {
      localdom.propagate();
      prune = localdom.infeasible();
      if (prune)
        localdom.conflictAnalysis(getConflictPool(),
                                  mipworker.getGlobalDomain(), pseudocost);
    }
    if (!prune) {
      getSymmetries().propagateOrbitopes(localdom);
      prune = localdom.infeasible();
    }
    if (!prune && passStabilizerToChildNode && currnode.stabilizerOrbits) {
      currnode.stabilizerOrbits->orbitalFixing(localdom);
      prune = localdom.infeasible();
    }
    if (prune) {
      localdom.backtrack();
      localdom.clearChangedCols(numChangedCols);
      if (countTreeWeight) treeweight += std::ldexp(1.0, -getCurrentDepth());
      continue;
    }

    nodelb = std::max(nodelb, localdom.getObjectiveLowerBound());
    bool nodeToQueue = nodelb > mipworker.getOptimalityLimit();

    if (nodeToQueue) {
      // if (!mipsolver.submip) printf("node goes to queue\n");
      stashNodeToProcessed(nodelb, nodestack.back().estimate,
                           getCurrentDepth() + 1);
      localdom.backtrack();
      localdom.clearChangedCols(numChangedCols);
      continue;
    }
    nodestack.emplace_back(
        nodelb, currnode.estimate, currnode.nodeBasis,
        passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);

    lp->flushDomain(localdom);
    nodestack.back().domgchgStackPos = domchgPos;
    break;
  }

  if (nodestack.back().nodeBasis) {
    lp->setStoredBasis(nodestack.back().nodeBasis);
    lp->recoverBasis();
  }

  return true;
}

bool HighsSearch::backtrackUntilDepth(HighsInt targetDepth) {
  if (nodestack.empty()) return false;
  assert(!nodestack.empty());
  if (getCurrentDepth() >= targetDepth) nodestack.back().opensubtrees = 0;

  while (nodestack.back().opensubtrees == 0) {
    depthoffset += nodestack.back().skipDepthCount;
    nodestack.pop_back();

#ifndef NDEBUG
    HighsDomainChange branchchg =
#endif
        localdom.backtrack();
    if (nodestack.empty()) {
      lp->flushDomain(localdom);
      return false;
    }
    assert(
        (branchchg.boundtype == HighsBoundType::kLower &&
         branchchg.boundval >= nodestack.back().branchingdecision.boundval) ||
        (branchchg.boundtype == HighsBoundType::kUpper &&
         branchchg.boundval <= nodestack.back().branchingdecision.boundval));
    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);
    assert(branchchg.column == nodestack.back().branchingdecision.column);

    if (getCurrentDepth() >= targetDepth) nodestack.back().opensubtrees = 0;
  }

  NodeData& currnode = nodestack.back();
  assert(currnode.opensubtrees == 1);
  currnode.opensubtrees = 0;
  flipBranchingDecision(currnode);

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  bool passStabilizerToChildNode =
      orbitsValidInChildNode(currnode.branchingdecision);
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(
      currnode.lower_bound, currnode.estimate, currnode.nodeBasis,
      passStabilizerToChildNode ? currnode.stabilizerOrbits : nullptr);

  lp->flushDomain(localdom);
  nodestack.back().domgchgStackPos = domchgPos;
  if (nodestack.back().nodeBasis &&
      (HighsInt)nodestack.back().nodeBasis->row_status.size() ==
          lp->getLp().num_row_)
    lp->setStoredBasis(nodestack.back().nodeBasis);
  lp->recoverBasis();

  return true;
}

HighsSearch::NodeResult HighsSearch::dive(int64_t nodeLim) {
  reliableatnode.clear();

  do {
    ++nnodes;
    NodeResult result = evaluateNode();

    if (checkLimits(nnodes)) return result;

    if (result != NodeResult::kOpen) return result;

    result = branch();
    if (result != NodeResult::kBranched) return result;
    if (nnodes >= nodeLim) return result;
  } while (true);
}

void HighsSearch::solveDepthFirst(int64_t maxbacktracks, int64_t maxnodes) {
  do {
    if (maxbacktracks == 0) break;

    NodeResult result = dive(maxnodes);
    // if a limit was reached the result might be open
    if (result == NodeResult::kOpen) break;
    if (nnodes >= maxnodes) break;

    --maxbacktracks;

  } while (backtrack());
}

double HighsSearch::getFeasTol() const { return mipsolver.mipdata_->feastol; }

double HighsSearch::getUpperLimit() const {
  if (!mipsolver.mipdata_->parallelLockActive()) {
    return mipsolver.mipdata_->upper_limit;
  } else {
    return mipworker.upper_limit;
  }
}

double HighsSearch::getEpsilon() const { return mipsolver.mipdata_->epsilon; }

const std::vector<double>& HighsSearch::getRootLpSol() const {
  return mipsolver.mipdata_->rootlpsol;
}

const std::vector<HighsInt>& HighsSearch::getIntegralCols() const {
  return mipsolver.mipdata_->integral_cols;
}

HighsDomain& HighsSearch::getDomain() const {
  return mipworker.getGlobalDomain();
}

HighsConflictPool& HighsSearch::getConflictPool() const {
  return mipworker.getConflictPool();
}

HighsCutPool& HighsSearch::getCutPool() const { return mipworker.getCutPool(); }

const HighsNodeQueue& HighsSearch::getNodeQueue() const {
  return mipsolver.mipdata_->nodequeue;
}

bool HighsSearch::checkLimits(int64_t nodeOffset) const {
  if (mipsolver.mipdata_->parallelLockActive()) {
    return checkLocalLimits();
  };
  return mipsolver.mipdata_->checkLimits(nodeOffset);
}

bool HighsSearch::checkLocalLimits() const {
  if (mipsolver.mipdata_->terminatorActive())
    if (mipsolver.mipdata_->terminatorTerminated()) return true;

  const int64_t stop = mipsolver.mipdata_->worker_lp_iterations_stop.load(
      std::memory_order_relaxed);
  if (stop <= lpiterations) {
    return true;
  }

  if (!mipsolver.submip && mipworker.upper_bound < kHighsInf &&
      mipsolver.options_mip_->objective_target > -kHighsInf) {
    const double internal_target =
        static_cast<HighsInt>(mipsolver.orig_model_->sense_) *
            mipsolver.options_mip_->objective_target -
        mipsolver.model_->offset_;
    if (mipworker.upper_bound < internal_target) {
      return true;
    }
  }

  if (mipsolver.options_mip_->mip_max_nodes != kHighsIInf &&
      mipsolver.mipdata_->num_nodes + nnodes >=
          mipsolver.options_mip_->mip_max_nodes) {
    return true;
  }

  if (mipsolver.options_mip_->mip_max_leaves != kHighsIInf &&
      mipsolver.mipdata_->num_leaves + nleaves >=
          mipsolver.options_mip_->mip_max_leaves) {
    return true;
  }

  if (mipsolver.options_mip_->time_limit < kHighsInf &&
      mipsolver.timer_.read() >= mipsolver.options_mip_->time_limit) {
    return true;
  }

  return false;
}

HighsSymmetries& HighsSearch::getSymmetries() const {
  return mipsolver.mipdata_->symmetries;
}

bool HighsSearch::addIncumbent(const std::vector<double>& sol, double solobj,
                               const int solution_source,
                               const bool print_display_line) {
  if (mipsolver.mipdata_->parallelLockActive()) {
    return mipworker.addIncumbent(sol, solobj, solution_source);
  } else {
    return mipsolver.mipdata_->addIncumbent(sol, solobj, solution_source,
                                            print_display_line);
  }
}
