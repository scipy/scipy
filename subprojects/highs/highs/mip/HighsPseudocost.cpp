/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsPseudocost.h"

#include "mip/HighsMipSolverData.h"

HighsPseudocost::HighsPseudocost(const HighsMipSolver& mipsolver)
    : pseudocostup(mipsolver.numCol()),
      pseudocostdown(mipsolver.numCol()),
      nsamplesup(mipsolver.numCol()),
      nsamplesdown(mipsolver.numCol()),
      inferencesup(mipsolver.numCol()),
      inferencesdown(mipsolver.numCol()),
      ninferencesup(mipsolver.numCol()),
      ninferencesdown(mipsolver.numCol()),
      ncutoffsup(mipsolver.numCol()),
      ncutoffsdown(mipsolver.numCol()),
      conflictscoreup(mipsolver.numCol()),
      conflictscoredown(mipsolver.numCol()),
      changedpos(mipsolver.numCol(), -1),
      conflict_weight(1.0),
      conflict_avg_score(0.0),
      cost_total(0),
      inferences_total(0),
      delta_cost_sum(0.0),
      delta_inferences_sum(0.0),
      nsamplestotal(0),
      ninferencestotal(0),
      ncutoffstotal(0),
      delta_nsamplestotal(0),
      delta_ninferencestotal(0),
      minreliable(mipsolver.options_mip_->mip_pscost_minreliable),
      degeneracyFactor(1.0) {
  deltas.reserve(std::min(HighsInt{256}, mipsolver.numCol()));
  if (mipsolver.pscostinit != nullptr) {
    cost_total = mipsolver.pscostinit->cost_total;
    inferences_total = mipsolver.pscostinit->inferences_total;
    nsamplestotal = mipsolver.pscostinit->nsamplestotal;
    ninferencestotal = mipsolver.pscostinit->ninferencestotal;

    conflict_avg_score =
        mipsolver.pscostinit->conflict_avg_score * mipsolver.numCol();

    for (HighsInt i : mipsolver.mipdata_->postSolveStack.getOrigCols()) {
      HighsInt origCol = mipsolver.mipdata_->postSolveStack.getOrigColIndex(i);

      pseudocostup[i] = mipsolver.pscostinit->pseudocostup[origCol];
      nsamplesup[i] = mipsolver.pscostinit->nsamplesup[origCol];
      pseudocostdown[i] = mipsolver.pscostinit->pseudocostdown[origCol];
      nsamplesdown[i] = mipsolver.pscostinit->nsamplesdown[origCol];
      inferencesup[i] = mipsolver.pscostinit->inferencesup[origCol];
      ninferencesup[i] = mipsolver.pscostinit->ninferencesup[origCol];
      inferencesdown[i] = mipsolver.pscostinit->inferencesdown[origCol];
      ninferencesdown[i] = mipsolver.pscostinit->ninferencesdown[origCol];
      conflictscoreup[i] = mipsolver.pscostinit->conflictscoreup[origCol];
      conflictscoredown[i] = mipsolver.pscostinit->conflictscoredown[origCol];
    }
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount)
    : pseudocostup(pscost.pseudocostup),
      pseudocostdown(pscost.pseudocostdown),
      nsamplesup(pscost.nsamplesup),
      nsamplesdown(pscost.nsamplesdown),
      inferencesup(pscost.inferencesup),
      inferencesdown(pscost.inferencesdown),
      ninferencesup(pscost.ninferencesup),
      ninferencesdown(pscost.ninferencesdown),
      conflictscoreup(pscost.conflictscoreup.size()),
      conflictscoredown(pscost.conflictscoreup.size()),
      cost_total(pscost.cost_total),
      inferences_total(pscost.inferences_total),
      conflict_avg_score(pscost.conflict_avg_score),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)),
      ninferencestotal(std::min(int64_t{1}, pscost.ninferencestotal)) {
  HighsInt ncol = pseudocostup.size();
  conflict_avg_score /= ncol * pscost.conflict_weight;
  for (HighsInt i = 0; i != ncol; ++i) {
    nsamplesup[i] = std::min(nsamplesup[i], maxCount);
    nsamplesdown[i] = std::min(nsamplesdown[i], maxCount);
    ninferencesup[i] = std::min(ninferencesup[i], HighsInt{1});
    ninferencesdown[i] = std::min(ninferencesdown[i], HighsInt{1});
    conflictscoreup[i] = pscost.conflictscoreup[i] / pscost.conflict_weight;
    conflictscoredown[i] = pscost.conflictscoredown[i] / pscost.conflict_weight;
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount,
    const presolve::HighsPostsolveStack& postsolveStack)
    : cost_total(pscost.cost_total),
      inferences_total(pscost.inferences_total),
      conflict_avg_score(pscost.conflict_avg_score),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)),
      ninferencestotal(std::min(int64_t{1}, pscost.ninferencestotal)) {
  pseudocostup.resize(postsolveStack.getOrigNumCol());
  pseudocostdown.resize(postsolveStack.getOrigNumCol());
  nsamplesup.resize(postsolveStack.getOrigNumCol());
  nsamplesdown.resize(postsolveStack.getOrigNumCol());
  inferencesup.resize(postsolveStack.getOrigNumCol());
  inferencesdown.resize(postsolveStack.getOrigNumCol());
  ninferencesup.resize(postsolveStack.getOrigNumCol());
  ninferencesdown.resize(postsolveStack.getOrigNumCol());
  conflictscoreup.resize(postsolveStack.getOrigNumCol());
  conflictscoredown.resize(postsolveStack.getOrigNumCol());

  HighsInt ncols = pscost.pseudocostup.size();
  conflict_avg_score /= ncols * pscost.conflict_weight;

  for (HighsInt i : postsolveStack.getOrigCols()) {
    pseudocostup[postsolveStack.getOrigColIndex(i)] = pscost.pseudocostup[i];
    pseudocostdown[postsolveStack.getOrigColIndex(i)] =
        pscost.pseudocostdown[i];
    nsamplesup[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesup[i]);
    nsamplesdown[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesdown[i]);
    inferencesup[postsolveStack.getOrigColIndex(i)] = pscost.inferencesup[i];
    inferencesdown[postsolveStack.getOrigColIndex(i)] =
        pscost.inferencesdown[i];
    ninferencesup[postsolveStack.getOrigColIndex(i)] = 1;
    ninferencesdown[postsolveStack.getOrigColIndex(i)] = 1;
    conflictscoreup[postsolveStack.getOrigColIndex(i)] =
        pscost.conflictscoreup[i] / pscost.conflict_weight;
    conflictscoredown[postsolveStack.getOrigColIndex(i)] =
        pscost.conflictscoredown[i] / pscost.conflict_weight;
  }
}
