/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsRedcostFixing.h"

#include "mip/HighsMipSolverData.h"

std::vector<std::pair<double, HighsDomainChange>>
HighsRedcostFixing::getLurkingBounds(const HighsMipSolver& mipsolver,
                                     const HighsDomain& globaldom) const {
  std::vector<std::pair<double, HighsDomainChange>> domchgs;
  if (lurkingColLower.empty()) return domchgs;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    for (const auto& lower : lurkingColLower[col]) {
      if (lower.second > globaldom.col_lower_[col])
        domchgs.emplace_back(
            lower.first, HighsDomainChange{static_cast<double>(lower.second),
                                           col, HighsBoundType::kLower});
    }

    for (const auto& upper : lurkingColUpper[col]) {
      if (upper.second < globaldom.col_upper_[col])
        domchgs.emplace_back(
            upper.first, HighsDomainChange{static_cast<double>(upper.second),
                                           col, HighsBoundType::kUpper});
    }
  }

  return domchgs;
}

void HighsRedcostFixing::propagateRootRedcost(const HighsMipSolver& mipsolver) {
  if (lurkingColLower.empty()) return;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    lurkingColLower[col].erase(
        lurkingColLower[col].begin(),
        lurkingColLower[col].upper_bound(mipsolver.mipdata_->lower_bound));
    lurkingColUpper[col].erase(
        lurkingColUpper[col].begin(),
        lurkingColUpper[col].upper_bound(mipsolver.mipdata_->lower_bound));

    for (auto it =
             lurkingColLower[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColLower[col].end(); ++it) {
      if (it->second > mipsolver.mipdata_->getDomain().col_lower_[col]) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kLower, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
      }
    }

    for (auto it =
             lurkingColUpper[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColUpper[col].end(); ++it) {
      if (it->second < mipsolver.mipdata_->getDomain().col_upper_[col]) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kUpper, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
      }
    }
  }

  mipsolver.mipdata_->getDomain().propagate();
}

void HighsRedcostFixing::propagateRedCost(const HighsMipSolver& mipsolver,
                                          HighsDomain& localdomain,
                                          HighsDomain& globaldom,
                                          const HighsLpRelaxation& lp,
                                          HighsConflictPool& conflictpool,
                                          HighsPseudocost& pseudocost,
                                          double upper_limit) {
  const std::vector<double>& lpredcost = lp.getSolution().col_dual;
  double lpobjective = lp.getObjective();
  HighsCDouble gap = static_cast<HighsCDouble>(upper_limit) - lpobjective;

  double tolerance = std::max(10 * mipsolver.mipdata_->feastol,
                              mipsolver.mipdata_->epsilon * double(gap));

  assert(!localdomain.infeasible());
  std::vector<HighsDomainChange> boundChanges;
  boundChanges.reserve(mipsolver.mipdata_->integral_cols.size());
  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    // lpobj + (col - bnd) * redcost <= cutoffbound
    // (col - bnd) * redcost <= gap
    // redcost * col <= gap + redcost * bnd
    //   case bnd is upper bound  => redcost < 0 :
    //      col >= gap/redcost + ub
    //      <=> redcost < gap / (lb - ub)
    //   case bnd is lower bound  => redcost > 0 :
    //      col <= gap/redcost + lb
    //      <=> redcost > gap / (ub - lb)
    if (localdomain.col_upper_[col] == localdomain.col_lower_[col]) continue;
    if (std::fabs(lpredcost[col]) <= tolerance) continue;

    double maxIncrease = lpredcost[col] * (localdomain.col_upper_[col] -
                                           localdomain.col_lower_[col]);
    if (maxIncrease > gap) {
      assert(localdomain.col_lower_[col] != -kHighsInf);
      assert(lpredcost[col] > tolerance);
      double newub =
          double(floor(gap / lpredcost[col] + localdomain.col_lower_[col] +
                       mipsolver.mipdata_->feastol));
      if (newub >= localdomain.col_upper_[col]) continue;
      assert(newub < localdomain.col_upper_[col]);

      if (globaldom.isBinary(col)) {
        boundChanges.emplace_back(
            HighsDomainChange{newub, col, HighsBoundType::kUpper});
      } else {
        localdomain.changeBound(HighsBoundType::kUpper, col, newub,
                                HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) return;
      }
    } else if (maxIncrease < -gap) {
      assert(localdomain.col_upper_[col] != kHighsInf);
      assert(lpredcost[col] < -tolerance);
      double newlb =
          double(ceil(gap / lpredcost[col] + localdomain.col_upper_[col] -
                      mipsolver.mipdata_->feastol));

      if (newlb <= localdomain.col_lower_[col]) continue;
      assert(newlb > localdomain.col_lower_[col]);

      if (globaldom.isBinary(col)) {
        boundChanges.emplace_back(
            HighsDomainChange{newlb, col, HighsBoundType::kLower});
      } else {
        localdomain.changeBound(HighsBoundType::kLower, col, newlb,
                                HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) return;
      }
    }
  }

  if (!boundChanges.empty()) {
    std::vector<HighsInt> inds;
    std::vector<double> vals;
    double rhs;

    if (boundChanges.size() <= 100 &&
        lp.computeDualProof(globaldom, upper_limit, inds, vals, rhs, false)) {
      bool addedConstraints = false;

      HighsInt oldNumConflicts = conflictpool.getNumConflicts();
      for (const HighsDomainChange& domchg : boundChanges) {
        if (localdomain.isActive(domchg)) continue;
        localdomain.conflictAnalyzeReconvergence(
            domchg, inds.data(), vals.data(),
            static_cast<HighsInt>(inds.size()), rhs, conflictpool, globaldom,
            pseudocost);
      }
      addedConstraints = conflictpool.getNumConflicts() != oldNumConflicts;

      if (addedConstraints) {
        localdomain.propagate();
        if (localdomain.infeasible()) return;

        boundChanges.erase(
            std::remove_if(boundChanges.begin(), boundChanges.end(),
                           [&](const HighsDomainChange& domchg) {
                             return localdomain.isActive(domchg);
                           }),
            boundChanges.end());
      }

      if (!boundChanges.empty()) {
        for (const HighsDomainChange& domchg : boundChanges) {
          localdomain.changeBound(domchg, HighsDomain::Reason::unspecified());
          if (localdomain.infeasible()) break;
        }

        if (!localdomain.infeasible()) localdomain.propagate();
      }
      // /printf("numConflicts: %d\n", numConflicts);
    } else {
      for (const HighsDomainChange& domchg : boundChanges) {
        localdomain.changeBound(domchg, HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) break;
      }

      if (!localdomain.infeasible()) localdomain.propagate();
    }
  }
}

void HighsRedcostFixing::addRootRedcost(const HighsMipSolver& mipsolver,
                                        const std::vector<double>& lpredcost,
                                        double lpobjective) {
  lurkingColLower.resize(mipsolver.numCol());
  lurkingColUpper.resize(mipsolver.numCol());

  // Provided domains won't be used (only used for dual proof)
  mipsolver.mipdata_->getLp().computeBasicDegenerateDuals(
      mipsolver.mipdata_->feastol, mipsolver.mipdata_->getDomain(),
      mipsolver.mipdata_->getDomain(), mipsolver.mipdata_->getConflictPool(),
      mipsolver.mipdata_->getPseudoCost(), false);

  // Compute maximum number of steps per column with large domain
  // max_steps = 2 ** k, k = max(5, min(10 ,round(log(|D| / 10)))),
  // D = {col : integral_cols | (ub - lb) >= 512}
  // This is to avoid doing 2**10 steps when there's many unbounded columns
  HighsInt numRedcostLargeDomainCols = 0;
  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    if ((mipsolver.mipdata_->getDomain().col_upper_[col] -
         mipsolver.mipdata_->getDomain().col_lower_[col]) >= 512 &&
        std::abs(lpredcost[col]) > mipsolver.mipdata_->feastol) {
      numRedcostLargeDomainCols++;
    }
  }
  HighsInt maxNumStepsExp = 10;
  int expshift = 0;
  std::frexp(numRedcostLargeDomainCols / 10, &expshift);
  if (expshift > 5) {
    expshift = std::min(expshift, static_cast<int>(maxNumStepsExp));
    maxNumStepsExp = maxNumStepsExp - expshift + 5;
  }
  HighsInt maxNumSteps = static_cast<HighsInt>(1ULL << maxNumStepsExp);

  // bounds whose absolute value is at least this are not handled: they
  // cannot be cast to HighsInt (kHighsIInf) or, for 64-bit HighsInt, leave
  // no room for the integer arithmetic in findLurkingBounds and are not
  // exact integers in double precision (kHighsInvMacheps)
  const double maxBound =
      std::min(static_cast<double>(kHighsIInf), kHighsInvMacheps);

  // lambda for finding lurking bounds
  auto findLurkingBounds =
      [&](HighsInt col, HighsInt direction, double inputBound,
          double inputOtherBound, double lpObjective, double redCost,
          HighsInt maxNumSteps, HighsInt maxNumStepsExp,
          std::multimap<double, HighsInt>& lurkingBounds,
          std::multimap<double, HighsInt>& otherLurkingBounds) {
        // bounds whose absolute value is at least maxBound cannot be handled
        if (std::abs(inputBound) >= maxBound) return;
        HighsInt bound = static_cast<HighsInt>(inputBound);

        if (direction * redCost == kHighsInf) {
          lurkingBounds.clear();
          otherLurkingBounds.clear();
          lurkingBounds.emplace(-kHighsInf, bound);
          return;
        }

        int64_t lastBound;
        if (direction * inputOtherBound >= maxBound) {
          // other bound cannot be handled, treat it as infinite
          lastBound = static_cast<int64_t>(bound) + direction * maxNumSteps;
          // keep the lurking bounds in the range of HighsInt
          lastBound = direction * std::min(direction * lastBound,
                                           static_cast<int64_t>(kHighsIInf));
        } else {
          lastBound = static_cast<int64_t>(inputOtherBound) - direction;
        }
        int64_t step = 1;
        int64_t range = direction * (lastBound - static_cast<int64_t>(bound));
        if (range > maxNumSteps)
          step = (range + maxNumSteps - 1) >> maxNumStepsExp;
        double shift = direction * (1 - 10 * mipsolver.mipdata_->feastol);
        step *= direction;

        for (int64_t lurkingBound = bound;
             direction * lurkingBound <= direction * lastBound;
             lurkingBound += step) {
          double fracBound = lurkingBound - bound + shift;
          double requiredCutoffBound = fracBound * redCost + lpObjective;
          if (requiredCutoffBound <
              mipsolver.mipdata_->lower_bound + mipsolver.mipdata_->feastol)
            continue;

          // check if we already have a better lurking bound stored
          bool useful = true;
          auto pos = lurkingBounds.lower_bound(requiredCutoffBound);
          for (auto it = pos; it != lurkingBounds.end(); ++it) {
            useful =
                direction * it->second >= direction * (lurkingBound + step);
            if (!useful) break;
          }
          if (!useful) continue;

          // we have no better lurking bound stored; store this lurking bound
          // and check if it dominates one that is already stored
          auto it = lurkingBounds.emplace_hint(
              pos, requiredCutoffBound, static_cast<HighsInt>(lurkingBound));

          auto i = lurkingBounds.begin();
          while (i != it) {
            if (direction * i->second >= direction * lurkingBound) {
              auto del = i++;
              lurkingBounds.erase(del);
            } else {
              ++i;
            }
          }
        }
      };

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    double lb = mipsolver.mipdata_->getDomain().col_lower_[col];
    double ub = mipsolver.mipdata_->getDomain().col_upper_[col];

    if (lpredcost[col] > mipsolver.mipdata_->feastol) {
      // col <= (cutoffbound - lpobj)/redcost + lb
      // so for lurkub = lb to ub - 1 we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurkub = (cutoffbound - lpobj)/redcost + lb
      //  cutoffbound = (lurkub - lb) * redcost + lpobj
      findLurkingBounds(col, HighsInt{1}, lb, ub, lpobjective, lpredcost[col],
                        maxNumSteps, maxNumStepsExp, lurkingColUpper[col],
                        lurkingColLower[col]);
    } else if (lpredcost[col] < -mipsolver.mipdata_->feastol) {
      // col >= (cutoffbound - lpobj)/redcost + ub
      // so for lurklb = lb + 1 to ub we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurklb = (cutoffbound - lpobj)/redcost + ub
      //  cutoffbound = (lurklb - ub) * redcost + lpobj
      findLurkingBounds(col, HighsInt{-1}, ub, lb, lpobjective, lpredcost[col],
                        maxNumSteps, maxNumStepsExp, lurkingColLower[col],
                        lurkingColUpper[col]);
    }
  }
}
