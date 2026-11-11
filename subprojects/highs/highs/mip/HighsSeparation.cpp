/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsSeparation.h"

#include <algorithm>
#include <cassert>
#include <queue>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpAggregator.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMachineSchedSeparator.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsModkSeparator.h"
#include "mip/HighsPathSeparator.h"
#include "mip/HighsTableauSeparator.h"
#include "mip/HighsTransformedLp.h"

HighsSeparation::HighsSeparation(HighsMipWorker& mipworker)
    : mipworker_(mipworker) {
  /*
  if (mipworker.mipsolver_.profiling_->mip_) {
    implBoundClock =
        mipworker.mipsolver_.profiling_->getSepaClockIndex(kImplboundSepaString);
    cliqueClock =
        mipworker.mipsolver_.profiling_->getSepaClockIndex(kCliqueSepaString);
  }
  */
  implBoundClock = 990;
  cliqueClock = 991;
  const HighsMipSolver& mipsolver = mipworker.getMipSolver();
  separators.emplace_back(new HighsTableauSeparator(mipsolver));
  separators.emplace_back(new HighsPathSeparator(mipsolver));
  separators.emplace_back(new HighsModkSeparator(mipsolver));
  separators.emplace_back(new HighsMachineSchedSeparator(mipsolver));
}

HighsInt HighsSeparation::separationRound(HighsDomain& propdomain,
                                          HighsLpRelaxation::Status& status) {
  const HighsSolution& sol = lp->getLpSolver().getSolution();

  HighsMipSolverData& mipdata = *lp->getMipSolver().mipdata_;

  auto propagateAndResolve = [&]() {
    if (propdomain.infeasible() || mipworker_.getGlobalDomain().infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    propdomain.propagate();
    if (propdomain.infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    // only modify cliquetable for master worker.
    if (&propdomain == &mipdata.getDomain())
      mipdata.cliquetable.cleanupFixed(mipdata.getDomain());

    if (mipworker_.getGlobalDomain().infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      propdomain.clearChangedCols();
      return -1;
    }

    int numBoundChgs = (int)propdomain.getChangedCols().size();

    while (!propdomain.getChangedCols().empty()) {
      lp->setObjectiveLimit(mipworker_.upper_limit);
      status = lp->resolveLp(&propdomain);
      if (!lp->scaledOptimal(status)) return -1;

      if (&propdomain == &mipdata.getDomain() &&
          lp->unscaledDualFeasible(status)) {
        mipdata.redcostfixing.addRootRedcost(
            mipdata.mipsolver, lp->getSolution().col_dual, lp->getObjective());
        if (mipworker_.upper_limit != kHighsInf)
          mipdata.redcostfixing.propagateRootRedcost(mipdata.mipsolver);
      }
    }

    return numBoundChgs;
  };

  if (!mipdata.parallelLockActive())
    lp->getMipSolver().profiling_->start(implBoundClock);
  mipdata.implications.separateImpliedBounds(
      *lp, lp->getSolution().col_value, mipworker_.getCutPool(),
      mipdata.feastol, mipworker_.getGlobalDomain(),
      mipdata.parallelLockActive());
  if (!mipdata.parallelLockActive())
    lp->getMipSolver().profiling_->stop(implBoundClock);

  HighsInt ncuts = 0;
  HighsInt numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  if (!mipdata.parallelLockActive())
    lp->getMipSolver().profiling_->start(cliqueClock);
  mipdata.cliquetable.separateCliques(
      lp->getMipSolver(), sol.col_value, mipworker_.getCutPool(),
      mipdata.feastol,
      mipdata.parallelLockActive() ? mipworker_.randgen
                                   : mipdata.cliquetable.getRandgen(),
      mipdata.parallelLockActive()
          ? mipworker_.getNumNeighbourhoodQueries()
          : mipdata.cliquetable.getNumNeighbourhoodQueries());
  if (!mipdata.parallelLockActive())
    lp->getMipSolver().profiling_->stop(cliqueClock);

  numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  if (&propdomain != &mipworker_.getGlobalDomain())
    lp->computeBasicDegenerateDuals(
        mipdata.feastol, propdomain, mipworker_.getGlobalDomain(),
        mipworker_.getConflictPool(), mipworker_.getPseudocost(), true);

  HighsTransformedLp transLp(*lp, mipdata.implications,
                             mipworker_.getGlobalDomain());
  if (mipworker_.getGlobalDomain().infeasible()) {
    status = HighsLpRelaxation::Status::kInfeasible;
    return 0;
  }
  HighsLpAggregator lpAggregator(*lp);

  for (const std::unique_ptr<HighsSeparator>& separator : separators) {
    separator->run(*lp, lpAggregator, transLp, mipworker_.getCutPool());
    if (mipworker_.getGlobalDomain().infeasible()) {
      status = HighsLpRelaxation::Status::kInfeasible;
      return 0;
    }
  }

  numboundchgs = propagateAndResolve();
  if (numboundchgs == -1)
    return 0;
  else
    ncuts += numboundchgs;

  mipworker_.getCutPool().separate(sol.col_value, propdomain, cutset,
                                   mipdata.feastol, mipdata.cutpools);
  // Also separate the global cut pool
  if (&mipworker_.getCutPool() != &mipdata.getCutPool()) {
    mipdata.getCutPool().separate(sol.col_value, propdomain, cutset,
                                  mipdata.feastol, mipdata.cutpools, true);
  }

  if (cutset.numCuts() > 0) {
    ncuts += cutset.numCuts();
    lp->addCuts(cutset);
    status = lp->resolveLp(&propdomain);
    lp->performAging(true);

    // only for the master domain.
    if (&propdomain == &mipdata.getDomain() &&
        lp->unscaledDualFeasible(status)) {
      mipdata.redcostfixing.addRootRedcost(
          mipdata.mipsolver, lp->getSolution().col_dual, lp->getObjective());
      if (mipdata.upper_limit != kHighsInf)
        mipdata.redcostfixing.propagateRootRedcost(mipdata.mipsolver);
    }
  }

  return ncuts;
}

void HighsSeparation::separate(HighsDomain& propdomain) {
  HighsLpRelaxation::Status status = lp->getStatus();
  const HighsMipSolver& mipsolver = lp->getMipSolver();

  if (lp->scaledOptimal(status) && !lp->getFractionalIntegers().empty()) {
    // double firstobj = lp->getObjective();
    double firstobj = mipsolver.mipdata_->rootlpsolobj;

    while (lp->getObjective() < mipworker_.optimality_limit) {
      double lastobj = lp->getObjective();

      int64_t nlpiters = -lp->getNumLpIterations();
      HighsInt ncuts = separationRound(propdomain, status);
      nlpiters += lp->getNumLpIterations();

      if (mipsolver.mipdata_->parallelLockActive()) {
        mipworker_.getSepaLpIterations() += nlpiters;
      } else {
        mipsolver.mipdata_->sepa_lp_iterations += nlpiters;
        mipsolver.mipdata_->total_lp_iterations += nlpiters;
      }

      // printf("separated %" HIGHSINT_FORMAT " cuts\n", ncuts);

      // printf(
      //     "separation round %" HIGHSINT_FORMAT " at node %" HIGHSINT_FORMAT "
      //     added %" HIGHSINT_FORMAT " cuts objective changed " "from %g to %g,
      //     first obj is %g\n", nrounds, (HighsInt)nnodes, ncuts, lastobj,
      //     lp->getObjective(), firstobj);
      if (ncuts == 0 || !lp->scaledOptimal(status) ||
          lp->getFractionalIntegers().empty())
        break;

      // if the objective improved considerably we continue
      if ((lp->getObjective() - firstobj) <=
          std::max((lastobj - firstobj), mipsolver.mipdata_->feastol) * 1.01)
        break;
    }

    // printf("done separating\n");
  } else {
    // printf("no separation, just aging. status: %" HIGHSINT_FORMAT "\n",
    //        (HighsInt)status);
    lp->performAging(true);

    mipworker_.getCutPool().performAging();
  }
}
