/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipSolver.h"

#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsMipWorker.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolve.h"
#include "presolve/HighsPostsolveStack.h"
#include "presolve/PresolveComponent.h"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"

using std::fabs;

HighsMipSolver::HighsMipSolver(HighsCallback& callback,
                               const HighsOptions& options, const HighsLp& lp,
                               const HighsSolution& solution, bool submip,
                               HighsInt submip_level)
    : callback_(&callback),
      options_mip_(&options),
      model_(&lp),
      orig_model_(&lp),
      solution_objective_(kHighsInf),
      submip(submip),
      submip_level(submip_level),
      rootbasis(nullptr),
      pscostinit(nullptr),
      clqtableinit(nullptr),
      implicinit(nullptr) {
  timer_.setPrintfFlag(options_mip_->output_flag, options_mip_->log_to_console);
  assert(!submip || submip_level > 0);
  max_submip_level = 0;
  // Initialise empty terminator
  initialiseTerminator();
  assert(termination_status_ == HighsModelStatus::kNotset);
  if (solution.value_valid) {
#ifndef NDEBUG
    // MIP solver doesn't check row residuals, but they should be OK
    // so validate using assert
    bool valid, integral, feasible;
    assessLpPrimalSolution("For debugging: ", options, lp, solution, valid,
                           integral, feasible);
    assert(valid);
#endif
    // Original comment: "Initial solution can be infeasible, but need
    // to set values for violation and objective"
    //
    // The initial solution shouldn't be infeasible in the original
    // space, since it's tested externally, but logging of
    // infeasibility has been added to
    // HighsMipSolverData::checkAddSolution()
    MipViolation violation;
    HighsCDouble quad_solution_objective_;
    solutionFeasible(orig_model_, solution.col_value, &solution.row_value,
                     violation, quad_solution_objective_);
    violation.copy(bound_violation_, row_violation_, integrality_violation_);
    solution_objective_ = double(quad_solution_objective_);
    solution_ = solution.col_value;
  }
}

HighsMipSolver::~HighsMipSolver() = default;

template <class F>
void HighsMipSolver::runTask(F&& f, highs::parallel::TaskGroup& tg,
                             bool parallel_lock, bool force_serial,
                             const std::vector<HighsInt>& indices) {
  if (indices.empty()) return;
  setParallelLock(parallel_lock);
  const bool spawn_tasks = !force_serial && indices.size() > 1 &&
                           !options_mip_->mip_search_simulate_concurrency;
  for (HighsInt i : indices) {
    if (spawn_tasks) {
      tg.spawn([&f, i] { f(i); });
    } else {
      f(i);
    }
  }
  if (spawn_tasks) {
    tg.taskWait();
  }
  setParallelLock(false);
}

void HighsMipSolver::run() {
  modelstatus_ = HighsModelStatus::kNotset;
  max_submip_level = std::max(submip_level, max_submip_level);
  // Start the timer local to HighsMipSolver - independent of the
  // timer passed from Highs as a pointer that's used in
  // HighsProfiling
  this->timer_.start();
  improving_solution_file_ = nullptr;
  if (!submip && options_mip_->mip_improving_solution_file != "")
    improving_solution_file_ =
        fopen(options_mip_->mip_improving_solution_file.c_str(), "w");

  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  for (HighsInt iLp = 0; iLp < static_cast<HighsInt>(mipdata_->lps.size());
       iLp++)
    mipdata_->lps[iLp].setProfiling(this->profiling_);
  assert(profiling_);
  // The solve time clock shouldn't be running on entry
  assert(!profiling_->running(kSolveTime));
  profiling_->start(kPresolveTime);
  profiling_->start(kMipClockInit);
  mipdata_->init();
#ifdef HIGHS_DEBUGSOL
  mipdata_->debugSolution.activate();
  bool debugSolActive = false;
  std::swap(mipdata_->debugSolution.debugSolActive, debugSolActive);
#endif
  profiling_->stop(kMipClockInit);
  profiling_->start(kMipClockRunPresolve);
  mipdata_->runMipPresolve(options_mip_->presolve_reduction_limit);
  profiling_->stop(kMipClockRunPresolve);
  profiling_->stop(kPresolveTime);
  if (profiling_->mip_ && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed presolve\n", timer_.read());
  if (modelstatus_ == HighsModelStatus::kNotset) {
    // Time limit has been reached (in presolve)
    if (mipdata_->presolve_status == HighsPresolveStatus::kTimeout
        // I think this is a legacy timeout check dataing from before there was
        // a timeout in presolve
        //
        // || timer_.read() >= options_mip_->time_limit
    ) {
      modelstatus_ = HighsModelStatus::kTimeLimit;
    } else if (mipdata_->presolve_status == HighsPresolveStatus::kOutOfMemory ||
               mipdata_->presolve_status == HighsPresolveStatus::kException ||
               mipdata_->presolve_status ==
                   HighsPresolveStatus::kMatrixValueError) {
      // Handle failures in presolve
      if (mipdata_->presolve_status == HighsPresolveStatus::kOutOfMemory) {
        modelstatus_ = HighsModelStatus::kMemoryLimit;
      } else {
        modelstatus_ = HighsModelStatus::kPresolveError;
      }
    }
  }
  if (modelstatus_ != HighsModelStatus::kNotset) {
    HighsModelStatus presolveStatus = modelstatus_;
    if (presolveStatus == HighsModelStatus::kInfeasible &&
        solution_objective_ != kHighsInf &&
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance) {
      presolveStatus = HighsModelStatus::kOptimal;
    }
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "Presolve: %s\n",
                 utilModelStatusToString(presolveStatus).c_str());
    if (modelstatus_ == HighsModelStatus::kOptimal) {
      // Model has been reduced to empty in presolve, so optimality
      // has been deduced
      mipdata_->lower_bound = 0;
      mipdata_->upper_bound = 0;
      double solution_objective =
          mipdata_->transformNewIntegerFeasibleSolution(std::vector<double>());
      if (solution_objective == kHighsInf) {
        // Solution is not feasible (due to an error in presolve) so
        // don't return a model status of HighsModelStatus::kOptimal!
        modelstatus_ = HighsModelStatus::kSolveError;
      } else {
        // Set solution_objective = -kHighsInf, which ensures that the
        // incumbent is viewed as an improving solution in
        // saveReportMipSolution
        solution_objective = -kHighsInf;
      }
      // Pass solution_objective as new_upper_limit
      mipdata_->saveReportMipSolution(solution_objective);
    }
    cleanupSolve();
    return;
  }

  profiling_->start(kSolveTime);
  if (profiling_->mip_ && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - starting  setup\n", timer_.read());
  profiling_->start(kMipClockRunSetup);
#ifdef HIGHS_DEBUGSOL
  mipdata_->debugSolution.debugSolActive = debugSolActive;
#endif
  mipdata_->runSetup();
  profiling_->stop(kMipClockRunSetup);
  if (profiling_->mip_ && !submip)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "MIP-Timing: %11.2g - completed setup\n", timer_.read());

  if (mipdata_->getDomain().infeasible()) {
    cleanupSolve();
    return;
  }
  // Initialise master worker.
  mipdata_->workers.emplace_back(
      *this, &mipdata_->getLp(), &mipdata_->getDomain(),
      &mipdata_->getCutPool(), &mipdata_->getConflictPool(),
      &mipdata_->getPseudoCost());

  HighsMipWorker& master_worker = mipdata_->workers[0];

restart:
  if (modelstatus_ == HighsModelStatus::kNotset) {
    // Check limits have not been reached before evaluating root node
    if (mipdata_->checkLimits()) {
      cleanupSolve();
      return;
    }
    // Possibly query existence of an external solution
    if (!submip)
      mipdata_->queryExternalSolution(
          solution_objective_, kExternalMipSolutionQueryOriginAfterSetup);

    // Apply the trivial heuristics
    profiling_->start(kMipClockTrivialHeuristics);
    HighsModelStatus returned_model_status = mipdata_->trivialHeuristics();
    profiling_->stop(kMipClockTrivialHeuristics);
    if (modelstatus_ == HighsModelStatus::kNotset &&
        returned_model_status == HighsModelStatus::kInfeasible) {
      // trivialHeuristics can spot trivial infeasibility, so act on it
      modelstatus_ = returned_model_status;
      cleanupSolve();
      return;
    }
    // Apply the feasibility jump heuristic (if enabled)
    if (options_mip_->mip_heuristic_run_feasibility_jump) {
      profiling_->start(kMipClockFeasibilityJump);
      HighsModelStatus returned_model_status = mipdata_->feasibilityJump();
      profiling_->stop(kMipClockFeasibilityJump);
      if (modelstatus_ == HighsModelStatus::kNotset &&
          returned_model_status == HighsModelStatus::kInfeasible) {
        // feasibilityJump can spot trivial infeasibility, so act on it
        modelstatus_ = returned_model_status;
        cleanupSolve();
        return;
      }
    }
    // End of pre-root-node heuristics
    if (profiling_->mip_ && !submip)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "MIP-Timing: %11.2g - starting evaluate root node\n",
                   timer_.read());
    profiling_->start(kMipClockEvaluateRootNode);
    mipdata_->evaluateRootNode(master_worker);
    profiling_->stop(kMipClockEvaluateRootNode);
    if (this->terminate()) {
      modelstatus_ = this->terminationStatus();
      cleanupSolve();
      return;
    }
    if (profiling_->mip_ && !submip)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "MIP-Timing: %11.2g - completed evaluate root node\n",
                   timer_.read());
    // age 5 times to remove stored but never violated cuts after root
    // separation
    profiling_->start(kMipClockPerformAging0);
    mipdata_->getCutPool().performAging();
    mipdata_->getCutPool().performAging();
    mipdata_->getCutPool().performAging();
    mipdata_->getCutPool().performAging();
    mipdata_->getCutPool().performAging();
    profiling_->stop(kMipClockPerformAging0);
  }
  if (mipdata_->nodequeue.empty() || mipdata_->checkLimits()) {
    cleanupSolve();
    return;
  }

  mipdata_->updateLowerBound(mipdata_->nodequeue.getBestLowerBound());
  mipdata_->printDisplayLine();

  const HighsInt max_num_workers = getMaxNumWorkers();
  HighsInt num_workers = 1;
  highs::parallel::TaskGroup tg;

  auto destroyOldWorkers = [&]() {
    if (mipdata_->workers.size() <= 1) return;
    while (mipdata_->domains.size() > 1) {
      mipdata_->domains.pop_back();
    }
    while (mipdata_->lps.size() > 1) {
      mipdata_->lps.pop_back();
    }
    while (mipdata_->pseudocosts.size() > 1) {
      mipdata_->pseudocosts.pop_back();
    }
    while (mipdata_->workers.size() > 1) {
      mipdata_->workers.pop_back();
    }
    while (mipdata_->cutpools.size() > 1) {
      mipdata_->cutpools.pop_back();
    }
    while (mipdata_->conflictpools.size() > 1) {
      mipdata_->conflictpools.pop_back();
    }
  };

  auto createNewWorkers = [&](HighsInt num_new_workers) {
    if (num_new_workers <= 0) return;
    // Remove all cuts from non-global pool for copied LP
    mipdata_->lps.emplace_back(mipdata_->getLp());
    // Have to set the global sub-solver call time pointer here for
    // new worker 0, since
    // HighsLpRelaxation::removeWorkerSpecificRows(); solves an LP
    mipdata_->lps.back().setProfiling(this->profiling_);
    mipdata_->lps.back().removeWorkerSpecificRows();
    for (HighsInt i = 0; i != num_new_workers; ++i) {
      if (i != 0) {
        mipdata_->lps.emplace_back(mipdata_->lps.back());
        mipdata_->lps.back().setProfiling(this->profiling_);
      }
      mipdata_->domains.emplace_back(mipdata_->getDomain());
      mipdata_->cutpools.emplace_back(
          numCol(), options_mip_->mip_pool_age_limit,
          options_mip_->mip_pool_soft_limit, mipdata_->cutpools.size());
      mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                           options_mip_->mip_pool_soft_limit);
      mipdata_->domains.back().addCutpool(mipdata_->cutpools.back());
      assert(mipdata_->domains.back().getDomainChangeStack().empty());
      mipdata_->domains.back().addConflictPool(mipdata_->conflictpools.back());
      mipdata_->pseudocosts.emplace_back(*this);
      mipdata_->workers.emplace_back(
          *this, &mipdata_->lps.back(), &mipdata_->domains.back(),
          &mipdata_->cutpools.back(), &mipdata_->conflictpools.back(),
          &mipdata_->pseudocosts.back());
      mipdata_->lps.back().setMipWorker(mipdata_->workers.back());
      mipdata_->lps.back().notifyCutPoolsLpCopied(1);
      mipdata_->workers.back().randgen.initialise(options_mip_->random_seed +
                                                  mipdata_->workers.size() - 1);
      mipdata_->debugSolution.registerDomain(
          mipdata_->workers.back().search_ptr_->getLocalDomain());
    }
  };

  // Use case: Change pointers in master worker to local copies of global info
  auto constructAdditionalWorkerData = [&](HighsMipWorker& worker) {
    assert(mipdata_->cutpools.size() == 1 &&
           mipdata_->conflictpools.size() == 1);
    assert(&worker == &mipdata_->workers[0]);
    mipdata_->cutpools.emplace_back(numCol(), options_mip_->mip_pool_age_limit,
                                    options_mip_->mip_pool_soft_limit, 1);
    worker.setCutPool(&mipdata_->cutpools.back());
    mipdata_->conflictpools.emplace_back(5 * options_mip_->mip_pool_age_limit,
                                         options_mip_->mip_pool_soft_limit);
    worker.setConflictPool(&mipdata_->conflictpools.back());
    mipdata_->domains.emplace_back(mipdata_->getDomain());
    worker.setGlobalDomain(&mipdata_->domains.back());
    worker.getGlobalDomain().addCutpool(worker.getCutPool());
    assert(worker.getGlobalDomain().getDomainChangeStack().empty());
    worker.getGlobalDomain().addConflictPool(worker.getConflictPool());
    mipdata_->pseudocosts.emplace_back(*this);
    worker.setPseudocost(&mipdata_->pseudocosts.back());
    worker.getLpRelaxation().setMipWorker(worker);
    worker.resetSearch();
    worker.resetSepa();
  };

  auto syncSolutions = [&]() -> void {
    // Note: Upper bound / limit of workers updated via addIncumbent
    for (HighsMipWorker& worker : mipdata_->workers) {
      for (auto& sol : worker.solutions_) {
        mipdata_->addIncumbent(std::get<0>(sol), std::get<1>(sol),
                               std::get<2>(sol));
      }
      worker.solutions_.clear();
    }
  };

  auto syncPools = [&](std::vector<HighsInt>& indices) -> void {
    if (!mipdata_->hasMultipleWorkers() || mipdata_->parallelLockActive())
      return;
    for (const HighsInt i : indices) {
      mipdata_->workers[i].getConflictPool().syncConflictPool(
          mipdata_->getConflictPool());
      mipdata_->workers[i].getCutPool().syncCutPool(*this,
                                                    mipdata_->getCutPool());
    }
    mipdata_->getCutPool().performAging();
    mipdata_->getConflictPool().performAging();
  };

  auto syncGlobalDomain = [&](std::vector<HighsInt>& indices) -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    for (const HighsInt i : indices) {
      HighsMipWorker& worker = mipdata_->workers[i];
      const auto& domchgstack = worker.getGlobalDomain().getDomainChangeStack();
      for (const HighsDomainChange& domchg : domchgstack) {
        if ((domchg.boundtype == HighsBoundType::kLower &&
             domchg.boundval >
                 mipdata_->getDomain().col_lower_[domchg.column]) ||
            (domchg.boundtype == HighsBoundType::kUpper &&
             domchg.boundval <
                 mipdata_->getDomain().col_upper_[domchg.column])) {
          mipdata_->getDomain().changeBound(domchg,
                                            HighsDomain::Reason::unspecified());
        }
      }
    }
  };

  auto doResetWorkerDomain = [&](HighsInt i) {
    HighsMipWorker& worker = mipdata_->workers[i];
    for (const HighsDomainChange& domchg :
         mipdata_->getDomain().getDomainChangeStack()) {
      worker.getGlobalDomain().changeBound(domchg,
                                           HighsDomain::Reason::unspecified());
      if (worker.getGlobalDomain().infeasible()) {
        return;
      }
    }
#ifndef NDEBUG
    for (HighsInt col = 0; col < numCol(); ++col) {
      assert(mipdata_->getDomain().col_lower_[col] ==
             worker.getGlobalDomain().col_lower_[col]);
      assert(mipdata_->getDomain().col_upper_[col] ==
             worker.getGlobalDomain().col_upper_[col]);
    }
#endif
    worker.getGlobalDomain().setDomainChangeStack(
        std::vector<HighsDomainChange>());
    // Warning: Resetting local domain cannot be done in parallel (changes
    // propagationDomains of main pool)
    worker.search_ptr_->resetLocalDomain();
    worker.getGlobalDomain().clearChangedCols();
  };

  auto resetGlobalDomain = [&](bool force, bool resetWorkers) -> void {
    // if global propagation found bound changes, we update the domain
    if (!mipdata_->getDomain().getChangedCols().empty() || force) {
      profiling_->start(kMipClockUpdateLocalDomain);
      highsLogDev(options_mip_->log_options, HighsLogType::kInfo,
                  "added %" HIGHSINT_FORMAT " global bound changes\n",
                  (HighsInt)mipdata_->getDomain().getChangedCols().size());
      mipdata_->cliquetable.cleanupFixed(mipdata_->getDomain());
      if (mipdata_->hasMultipleWorkers() && resetWorkers) {
        // Sync worker domains here. cleanupFixed might have found extra changes
        std::vector<HighsInt> indices(num_workers);
        std::iota(indices.begin(), indices.end(), 0);
        runTask(doResetWorkerDomain, tg, false, true, indices);
      }
      for (const HighsInt col : mipdata_->getDomain().getChangedCols())
        mipdata_->implications.cleanupVarbounds(col);

      mipdata_->getDomain().setDomainChangeStack(
          std::vector<HighsDomainChange>());
      if (!mipdata_->hasMultipleWorkers())
        master_worker.search_ptr_->resetLocalDomain();
      mipdata_->getDomain().clearChangedCols();
      mipdata_->removeFixedIndices();
      profiling_->stop(kMipClockUpdateLocalDomain);
    }
  };

  auto syncGlobalPseudoCost = [&]() -> void {
    if (!mipdata_->hasMultipleWorkers()) return;
    for (HighsMipWorker& worker : mipdata_->workers) {
      mipdata_->getPseudoCost().flushPseudoCost(worker.getPseudocost());
    }
  };

  auto resetWorkerPseudoCosts = [&](std::vector<HighsInt>& indices) {
    if (!mipdata_->hasMultipleWorkers()) return;
    auto doResetWorkerPseudoCost = [&](HighsInt i) -> void {
      mipdata_->getPseudoCost().syncPseudoCost(
          mipdata_->workers[i].getPseudocost());
    };
    runTask(doResetWorkerPseudoCost, tg, false, false, indices);
  };

  master_worker.resetSearch();
  master_worker.resetSepa();
  master_worker.upper_bound = mipdata_->upper_bound;
  master_worker.upper_limit = mipdata_->upper_limit;
  master_worker.optimality_limit = mipdata_->optimality_limit;
  destroyOldWorkers();
  mipdata_->debugSolution.registerDomain(
      master_worker.search_ptr_->getLocalDomain());

  profiling_->start(kMipClockSearch);
  // Maximum node processing budget for each worker in a batch during parallel
  constexpr HighsInt kMaxNodesPerWorker = 100;
  HighsInt maxNodesPerWorkerLim =
      max_num_workers > 1 ? kMaxNodesPerWorker : 10000;
  int64_t numStallNodes = 0;
  int64_t lastLbLeave = 0;
  int64_t numQueueLeaves = 0;
  HighsInt numHugeTreeEstim = 0;
  int64_t numNodesLastCheck = mipdata_->num_nodes;
  int64_t nextCheck = mipdata_->num_nodes;
  double treeweightLastCheck = 0.0;
  double upperLimLastCheck = mipdata_->upper_limit;
  double lowerBoundLastCheck = mipdata_->lower_bound;

  enum class RestartVote { kNoCheck, kNoHugeTree, kHugeTree, kWouldRestart };

  auto nodesInstalled = [&]() -> bool {
    for (const HighsMipWorker& worker : mipdata_->workers) {
      if (worker.search_ptr_->hasNode()) return true;
    }
    return false;
  };

  auto checkRestart = [&](const HighsMipWorker& worker, HighsInt numWorkerVotes,
                          HighsCDouble currTreeWeight) -> RestartVote {
    int64_t nNodes = worker.search_ptr_->nnodes + mipdata_->num_nodes;
    if (!submip && nNodes >= nextCheck && options_mip_->mip_allow_restart) {
      const HighsInt nTreeRestarts =
          mipdata_->numRestarts - mipdata_->numRestartsRoot;
      const HighsCDouble treeWeight =
          currTreeWeight + mipdata_->pruned_treeweight;
      double currNodeEstim =
          numNodesLastCheck - mipdata_->num_nodes_before_run +
          (nNodes - numNodesLastCheck) * double(1.0 - treeWeight) /
              std::max(double(treeWeight - treeweightLastCheck),
                       mipdata_->epsilon);
      // printf(
      //     "nTreeRestarts: %d, numNodesThisRun: %ld, numNodesLastCheck: %ld,
      //     " "currNodeEstim: %g, " "prunedTreeWeightDelta: %g,
      //     numHugeTreeEstim: %d, numLeavesThisRun:
      //     "
      //     "%ld\n",
      //     nTreeRestarts, mipdata_->num_nodes -
      //     mipdata_->num_nodes_before_run, numNodesLastCheck -
      //     mipdata_->num_nodes_before_run, currNodeEstim, 100.0 *
      //     double(mipdata_->pruned_treeweight - treeweightLastCheck),
      //     numHugeTreeEstim,
      //     mipdata_->num_leaves - mipdata_->num_leaves_before_run);

      bool doRestart = false;

      double activeIntegerRatio =
          1.0 - mipdata_->percentageInactiveIntegers() / 100.0;
      activeIntegerRatio *= activeIntegerRatio;
      bool hugeTree = false;

      if (!doRestart) {
        double gapReduction = 1.0;
        if (worker.upper_limit != kHighsInf) {
          double oldGap = upperLimLastCheck - lowerBoundLastCheck;
          double newGap = worker.upper_limit - mipdata_->lower_bound;
          gapReduction = oldGap / newGap;
        }

        if (gapReduction < 1.0 + (0.05 / activeIntegerRatio) &&
            currNodeEstim >= activeIntegerRatio * 20 *
                                 (nNodes - mipdata_->num_nodes_before_run)) {
          hugeTree = true;
        } else {
          numWorkerVotes = -kHighsIInf;
        }

        int64_t minHugeTreeOffset =
            (mipdata_->num_leaves + worker.search_ptr_->nleaves -
             mipdata_->num_leaves_before_run) /
            1000;
        int64_t minHugeTreeEstim = HighsIntegers::nearestInteger(
            activeIntegerRatio * (10 + minHugeTreeOffset) *
            std::pow(1.5, nTreeRestarts));

        doRestart = numHugeTreeEstim + numWorkerVotes >= minHugeTreeEstim;
      }

      if (doRestart) return RestartVote::kWouldRestart;
      if (hugeTree) return RestartVote::kHugeTree;
      return RestartVote::kNoHugeTree;
    }
    return RestartVote::kNoCheck;
  };

  auto performRestart = [&] {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "\nRestarting search from the root node\n");
    mipdata_->performRestart();
    profiling_->stop(kMipClockSearch);
  };

  auto checkWorkerRestartVotes =
      [&](const std::vector<RestartVote>& votes) -> bool {
    HighsInt numVoters = 0;
    HighsInt numHugeTreeVotes = 0;
    HighsInt numRestartVotes = 0;
    for (const RestartVote vote : votes) {
      if (vote != RestartVote::kNoCheck) {
        ++numVoters;
        if (vote == RestartVote::kWouldRestart) {
          numRestartVotes++;
        } else if (vote == RestartVote::kHugeTree) {
          numHugeTreeVotes++;
        }
      }
    }
    bool forcedRestart = false;
    if (options_mip_->mip_allow_restart && numVoters >= 2 &&
        numRestartVotes / static_cast<double>(numVoters) >= 0.25) {
      forcedRestart = true;
    }
    // Force a restart if enough individual workers vote for it
    if (forcedRestart) {
      performRestart();
      return true;
    }
    // Using joint information after workers are synced, query a restart
    const RestartVote vote = checkRestart(
        master_worker, (numRestartVotes + numHugeTreeVotes + 2) / 2,
        master_worker.search_ptr_->treeweight);
    if (vote == RestartVote::kWouldRestart) {
      performRestart();
      return true;
    }

    // Update the appropriate values in case of no restart
    if (vote == RestartVote::kHugeTree) {
      nextCheck = mipdata_->num_nodes + 100;
      numHugeTreeEstim += (numRestartVotes + numHugeTreeVotes + 2) / 2;
    } else if (vote == RestartVote::kNoHugeTree) {
      numHugeTreeEstim = 0;
      treeweightLastCheck = static_cast<double>(mipdata_->pruned_treeweight);
      numNodesLastCheck = mipdata_->num_nodes;
      upperLimLastCheck = mipdata_->upper_limit;
      lowerBoundLastCheck = mipdata_->lower_bound;
    }
    return false;
  };

  auto getSearchIndicesWithNoNodes = [&](std::vector<HighsInt>& indices,
                                         std::vector<RestartVote>& restarts) {
    indices.clear();
    std::fill(restarts.begin(), restarts.end(), RestartVote::kNoCheck);
    const HighsInt num_search_workers =
        std::min(num_workers,
                 static_cast<HighsInt>(mipdata_->nodequeue.numActiveNodes()));

    for (HighsInt i = 0; i < num_search_workers; i++) {
      assert(!mipdata_->workers[i].search_ptr_->hasNode());
      indices.emplace_back(i);
      mipdata_->workers[i].setAllowHeuristics(true);
    }
  };

  auto prepareNodes = [&](const std::vector<HighsInt>& indices) {
    const HighsInt numNodesPerWorker =
        indices.size() > 1 &&
                static_cast<size_t>(mipdata_->nodequeue.numActiveNodes()) >=
                    2 * indices.size()
            ? 2
            : 1;
    for (HighsInt j = 0; j != numNodesPerWorker; j++) {
      for (const HighsInt i : indices) {
        if ((indices.size() == 1 && numQueueLeaves - lastLbLeave >= 10) ||
            (indices.size() > 1 && i == 0 && j == 0)) {
          mipdata_->workers[i].preparedNodes.push_back(
              std::move(mipdata_->nodequeue.popBestBoundNode()));
          lastLbLeave = numQueueLeaves;
        } else {
          HighsInt bestBoundNodeStackSize =
              mipdata_->nodequeue.getBestBoundDomchgStackSize();
          double bestBoundNodeLb = mipdata_->nodequeue.getBestLowerBound();
          mipdata_->workers[i].preparedNodes.push_back(
              std::move(mipdata_->nodequeue.popBestNode()));
          HighsNodeQueue::OpenNode& nextNode =
              mipdata_->workers[i].preparedNodes.back();
          if (nextNode.lower_bound == bestBoundNodeLb &&
              static_cast<HighsInt>(nextNode.domchgstack.size()) ==
                  bestBoundNodeStackSize)
            lastLbLeave = numQueueLeaves;
        }
        ++numQueueLeaves;
      }
    }
  };

  auto evaluateNode = [&](HighsInt i) -> bool {
    if (!mipdata_->parallelLockActive())
      profiling_->start(kMipClockEvaluateNode1);
    HighsMipWorker& worker = mipdata_->workers[i];
    if (worker.search_ptr_->evaluateNode() ==
        HighsSearch::NodeResult::kSubOptimal) {
      worker.search_ptr_->stashCurrentNode();
      if (!mipdata_->parallelLockActive())
        profiling_->stop(kMipClockEvaluateNode1);
      return true;
    }
    if (!mipdata_->parallelLockActive())
      profiling_->stop(kMipClockEvaluateNode1);
    return false;
  };

  auto assignEarlyTermination = [&](HighsMipWorker& worker) {
    if (worker.getGlobalDomain().infeasible() ||
        mipdata_->lower_bound > worker.getOptimalityLimit()) {
      mipdata_->updateWorkerEarlyTermination(worker);
    }
  };

  auto pruneNode = [&](const HighsInt i, bool& infeasible) -> bool {
    assert(infeasible == false);
    if (!mipdata_->parallelLockActive())
      profiling_->start(kMipClockNodePrunedLoop);
    bool pruned = false;
    if (mipdata_->workers[i].search_ptr_->currentNodePruned()) {
      mipdata_->workers[i].search_ptr_->backtrack();
      mipdata_->workers[i].getGlobalDomain().propagate();
      pruned = true;
      ++mipdata_->workers[i].search_ptr_->getLocalNodes();
      ++mipdata_->workers[i].search_ptr_->getLocalLeaves();
      if (mipdata_->workers[i].getGlobalDomain().infeasible()) {
        assignEarlyTermination(mipdata_->workers[i]);
        infeasible = true;
      }
    }
    if (!mipdata_->parallelLockActive())
      profiling_->stop(kMipClockNodePrunedLoop);
    return pruned;
  };

  auto separateAndStoreBasis = [&](const HighsInt i, bool& infeasible) -> void {
    HighsMipWorker& worker = mipdata_->workers[i];
    if (options_mip_->mip_allow_cut_separation_at_nodes) {
      if (!mipdata_->parallelLockActive())
        profiling_->start(kMipClockNodeSearchSeparation);
      worker.sepa_ptr_->separate(worker.search_ptr_->getLocalDomain());
      if (!mipdata_->parallelLockActive())
        profiling_->stop(kMipClockNodeSearchSeparation);
    } else {
      worker.getCutPool().performAging();
    }

    if (worker.getGlobalDomain().infeasible()) {
      worker.search_ptr_->cutoffNode();
      worker.search_ptr_->stashOpenNodes();
      assignEarlyTermination(worker);
      infeasible = true;
      return;
    }

    if (worker.getLpRelaxation().getStatus() !=
            HighsLpRelaxation::Status::kError &&
        worker.getLpRelaxation().getStatus() !=
            HighsLpRelaxation::Status::kNotSet)
      worker.getLpRelaxation().storeBasis();

    std::shared_ptr<const HighsBasis> basis =
        worker.getLpRelaxation().getStoredBasis();
    if (!basis ||
        !isBasisConsistent(worker.getLpRelaxation().getLp(), *basis)) {
      HighsBasis b = mipdata_->firstrootbasis;
      b.row_status.resize(worker.getLpRelaxation().numRows(),
                          HighsBasisStatus::kBasic);
      basis = std::make_shared<const HighsBasis>(std::move(b));
      worker.getLpRelaxation().setStoredBasis(basis);
    }
  };

  auto backtrackPlunge = [&](HighsInt i) {
    if (!mipdata_->parallelLockActive())
      profiling_->start(kMipClockBacktrackPlunge);
    HighsMipWorker& worker = mipdata_->workers[i];
    const bool backtrack_plunge = worker.search_ptr_->backtrackPlunge();
    if (!mipdata_->parallelLockActive())
      profiling_->stop(kMipClockBacktrackPlunge);

    if (!backtrack_plunge) return true;

    assert(worker.search_ptr_->hasNode());

    if (worker.getConflictPool().getNumConflicts() >
        options_mip_->mip_pool_soft_limit) {
      worker.getConflictPool().performAging();
    }
    return false;
  };

  auto runHeuristics = [&](const HighsInt i, bool& infeasible) -> bool {
    HighsMipWorker& worker = mipdata_->workers[i];
    if (!mipdata_->parallelLockActive())
      profiling_->start(kMipClockDiveEvaluateNode);
    const HighsSearch::NodeResult evaluate_node_result =
        worker.search_ptr_->evaluateNode();
    if (!mipdata_->parallelLockActive())
      profiling_->stop(kMipClockDiveEvaluateNode);

    if (evaluate_node_result == HighsSearch::NodeResult::kSubOptimal) {
      return true;
    }

    if (worker.search_ptr_->currentNodePruned()) {
      ++worker.search_ptr_->getLocalLeaves();
      return false;
    }

    if (!mipdata_->parallelLockActive())
      profiling_->start(kMipClockDivePrimalHeuristics);
    if (mipdata_->incumbent.empty()) {
      if (!mipdata_->parallelLockActive())
        profiling_->start(kMipClockDiveRandomizedRounding);
      mipdata_->heuristics.randomizedRounding(
          worker,
          worker.getLpRelaxation().getLpSolver().getSolution().col_value);
      if (!mipdata_->parallelLockActive())
        profiling_->stop(kMipClockDiveRandomizedRounding);
    }
    if (mipdata_->incumbent.empty()) {
      if (options_mip_->mip_heuristic_run_rens) {
        if (!mipdata_->parallelLockActive())
          profiling_->start(kMipClockDiveRens);
        mipdata_->heuristics.RENS(
            worker,
            worker.getLpRelaxation().getLpSolver().getSolution().col_value);
        if (!mipdata_->parallelLockActive())
          profiling_->stop(kMipClockDiveRens);
      }
    } else {
      if (options_mip_->mip_heuristic_run_rins) {
        if (!mipdata_->parallelLockActive())
          profiling_->start(kMipClockDiveRins);
        mipdata_->heuristics.RINS(
            worker,
            worker.getLpRelaxation().getLpSolver().getSolution().col_value);
        if (!mipdata_->parallelLockActive())
          profiling_->stop(kMipClockDiveRins);
      }
    }

    if (!mipdata_->parallelLockActive())
      profiling_->stop(kMipClockDivePrimalHeuristics);

    if (worker.getGlobalDomain().infeasible()) infeasible = true;
    return false;
  };

  auto dive = [&](HighsInt i, HighsInt nodeLim) -> bool {
    HighsMipWorker& worker = mipdata_->workers[i];
    const int64_t diveNodeLim = nodeLim == kHighsIInf
                                    ? std::numeric_limits<int64_t>::max()
                                    : worker.search_ptr_->nnodes + nodeLim;
    if (!worker.search_ptr_->currentNodePruned()) {
      if (!mipdata_->parallelLockActive()) profiling_->start(kMipClockTheDive);
      const HighsSearch::NodeResult search_dive_result =
          worker.search_ptr_->dive(diveNodeLim);
      if (!mipdata_->parallelLockActive()) profiling_->stop(kMipClockTheDive);
      if (search_dive_result == HighsSearch::NodeResult::kSubOptimal) {
        return true;
      }
      worker.search_ptr_->getLocalLeaves()++;
    }
    return nodeLim != kHighsIInf && worker.search_ptr_->nnodes >= diveNodeLim;
  };

  auto processNodes = [&](const std::vector<HighsInt>& indices,
                          std::vector<RestartVote>& restarts,
                          std::vector<HighsInt>& stallNodes,
                          const bool skip_separation, const HighsInt nodeLim,
                          const HighsInt plungeLimit, double avgiter) {
    auto processWorkerNodes = [&](const HighsInt i) {
      HighsMipWorker& worker = mipdata_->workers[i];
      worker.prepNodeIdx = 0;
      int64_t total_nodes_explored = 0;
      bool considerHeuristics = true;
      bool stop = false;
      while (worker.prepNodeIdx != worker.preparedNodes.size()) {
        int64_t nodes_explored = 0;
        int64_t lastLocalNodeCount = worker.search_ptr_->getLocalNodes();
        worker.search_ptr_->installNode(
            std::move(worker.preparedNodes[worker.prepNodeIdx]));
        ++worker.prepNodeIdx;
        if (worker.search_ptr_->getCurrentEstimate() >= worker.upper_limit) {
          stallNodes[i]++;
        }
        if (!skip_separation) {
          evaluateNode(i);
          assignEarlyTermination(worker);
          if (pruneNode(i, stop)) {
            total_nodes_explored++;
            if (worker.early_termination) stop = true;
            if (stop || total_nodes_explored >= nodeLim) break;
            continue;
          }
          if (stop || worker.search_ptr_->checkLimits()) {
            stop = true;
            break;
          };
          if (!mipdata_->parallelLockActive()) mipdata_->printDisplayLine();
          separateAndStoreBasis(i, stop);
          if (stop) break;
        }
        worker.getConflictPool().performAging();
        HighsInt iterlimit = 10 * std::max(avgiter, mipdata_->avgrootlpiters);
        iterlimit = std::max(
            {HighsInt{10000}, iterlimit,
             static_cast<HighsInt>((3 * mipdata_->firstrootlpiters) / 2)});
        worker.getLpRelaxation().setIterationLimit(iterlimit);
        while (true) {
          if (considerHeuristics &&
              (skip_separation || nodeLim == maxNodesPerWorkerLim) &&
              worker.getAllowHeuristics() &&
              mipdata_->moreHeuristicsAllowed()) {
            if (runHeuristics(i, stop)) break;
            if (stop) break;
          }
          considerHeuristics = false;
          if (worker.getGlobalDomain().infeasible()) {
            stop = true;
            break;
          }
          if (dive(i, nodeLim)) break;
          if (worker.search_ptr_->checkLimits(
                  worker.search_ptr_->getLocalNodes())) {
            stop = true;
            break;
          }

          nodes_explored +=
              worker.search_ptr_->getLocalNodes() - lastLocalNodeCount;
          lastLocalNodeCount = worker.search_ptr_->getLocalNodes();
          if (nodes_explored >= plungeLimit) break;
          if (backtrackPlunge(i)) break;
          if (!mipdata_->parallelLockActive()) {
            worker.search_ptr_->flushStatistics(*this);
            mipdata_->printDisplayLine();
            lastLocalNodeCount = 0;
          }
        }
        nodes_explored +=
            worker.search_ptr_->getLocalNodes() - lastLocalNodeCount;
        lastLocalNodeCount = worker.search_ptr_->getLocalNodes();
        total_nodes_explored += nodes_explored;
        assignEarlyTermination(worker);
        if (worker.early_termination) stop = true;
        worker.search_ptr_->stashOpenNodes();
        if (total_nodes_explored >= nodeLim || stop) {
          break;
        }
      }
      if (worker.search_ptr_->hasNode()) {
        worker.search_ptr_->stashOpenNodes();
      }
      for (size_t j = worker.prepNodeIdx; j != worker.preparedNodes.size();
           j++) {
        worker.processedNodes.emplace_back(std::move(worker.preparedNodes[j]),
                                           worker.search_ptr_->countTreeWeight);
      }
      worker.preparedNodes.clear();
      worker.prepNodeIdx = 0;
      if (nodeLim == maxNodesPerWorkerLim) {
        // Need to spoof tree weight
        HighsCDouble currTreeWeight = worker.search_ptr_->treeweight;
        for (auto& node_treeweight_pair : worker.processedNodes) {
          if (node_treeweight_pair.second &&
              node_treeweight_pair.first.lower_bound >
                  worker.optimality_limit) {
            currTreeWeight +=
                std::ldexp(1.0, 1 - node_treeweight_pair.first.depth);
          }
        }
        restarts[i] = checkRestart(worker, 1, currTreeWeight);
      }
    };
    runTask(processWorkerNodes, tg, true, false, indices);
  };

  auto syncSepaStats = [&](HighsMipWorker& worker) {
    mipdata_->cliquetable.getNumNeighbourhoodQueries() +=
        worker.getNumNeighbourhoodQueries();
    mipdata_->sepa_lp_iterations += worker.getSepaLpIterations();
    mipdata_->total_lp_iterations += worker.getSepaLpIterations();
    worker.resetSepaStats();
  };

  // Main solve loop
  std::vector<HighsInt> search_indices(1, 0);
  search_indices.reserve(max_num_workers);
  std::vector<RestartVote> workerRestartVotes(getMaxNumWorkers(),
                                              RestartVote::kNoCheck);
  std::vector<HighsInt> stallNodesPerWorker(getMaxNumWorkers());
  bool root_node = true;  // Don't separate the root node again
  HighsInt nodeLim = max_num_workers > 1 ? 1 : maxNodesPerWorkerLim;  // ramp-up
  while (!mipdata_->nodequeue.empty()) {
    // Possibly query existence of an external solution
    if (!submip)
      mipdata_->queryExternalSolution(
          solution_objective_, kExternalMipSolutionQueryOriginBeforeDive);

    // Update global pseudo-cost with worker information
    syncGlobalPseudoCost();

    // Get new candidate worker search indices
    getSearchIndicesWithNoNodes(search_indices, workerRestartVotes);

    // Only update worker's pseudo-costs that have been assigned a node
    resetWorkerPseudoCosts(search_indices);

    // Assign nodes to workers
    bool limit_reached = false;
    if (root_node) {
      master_worker.preparedNodes.push_back(
          std::move(mipdata_->nodequeue.popBestBoundNode()));
    } else {
      prepareNodes(search_indices);
    }

    if (nodeLim != maxNodesPerWorkerLim) {
      if (num_workers >= max_num_workers) {
        nodeLim = std::max(HighsInt{20}, 2 * nodeLim);
      }
      if (nodeLim >= maxNodesPerWorkerLim) {
        nodeLim = maxNodesPerWorkerLim;
      }
    }

    // Process nodes (separation / heuristics / dives)
    processNodes(search_indices, workerRestartVotes, stallNodesPerWorker,
                 root_node, nodeLim, 100, mipdata_->getLp().getAvgSolveIters());

    root_node = false;

    // Check if one worker sent an early global termination signal
    const int64_t early_termination_lp_iterations =
        mipdata_->worker_lp_iterations_stop.load(std::memory_order_relaxed);
    HighsInt early_terminated_worker = -1;
    if (early_termination_lp_iterations !=
        std::numeric_limits<int64_t>::max()) {
      for (const HighsInt i : search_indices) {
        if (mipdata_->workers[i].early_termination &&
            mipdata_->workers[i].search_ptr_->lpiterations ==
                early_termination_lp_iterations) {
          early_terminated_worker = i;
          break;
        }
      }
      assert(early_terminated_worker != -1);
      if (early_terminated_worker == -1) {
        mipdata_->worker_lp_iterations_stop.store(
            std::numeric_limits<int64_t>::max(), std::memory_order_relaxed);
        for (const HighsInt i : search_indices) {
          mipdata_->workers[i].early_termination = false;
        }
      } else {
        // Clear all buffered solve information from workers that are not
        // the earliest to terminate
        for (const HighsInt i : search_indices) {
          if (i == early_terminated_worker) continue;
          HighsMipWorker& worker = mipdata_->workers[i];
          worker.processedNodes.clear();
          worker.search_ptr_->resetStatistics();
          worker.resetHeurStats();
          worker.resetSepaStats();
          worker.solutions_.clear();
        }
        search_indices.clear();
        search_indices.emplace_back(early_terminated_worker);
      }
    }

    // Sync statistics, check infeasibility, and flush nodes from worker queues
    bool infeasible = false;
    for (const HighsInt i : search_indices) {
      HighsMipWorker& worker = mipdata_->workers[i];
      if (worker.getGlobalDomain().infeasible()) {
        infeasible = true;
      }
      profiling_->start(kMipClockOpenNodesToQueue0);
      assert(!worker.search_ptr_->hasNode());
      for (auto& node_treeweight_pair : worker.processedNodes) {
        HighsNodeQueue::OpenNode& node = node_treeweight_pair.first;
        double tmpTreeWeight = mipdata_->nodequeue.emplaceNode(
            std::move(node.domchgstack), std::move(node.branchings),
            node.lower_bound, node.estimate, node.depth);
        if (node_treeweight_pair.second) {
          worker.search_ptr_->treeweight += tmpTreeWeight;
        }
      }
      worker.processedNodes.clear();
      profiling_->stop(kMipClockOpenNodesToQueue0);
      worker.search_ptr_->flushStatistics(*this);
      syncSepaStats(worker);
      mipdata_->heuristics.flushStatistics(*this, worker);
    }

    if (infeasible) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
      break;
    }

    syncSolutions();

    if (early_terminated_worker == -1) {
      mipdata_->updateLowerBound(std::min(
          mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
    }

    limit_reached = mipdata_->checkLimits();
    if (limit_reached || early_terminated_worker != -1) {
      mipdata_->printDisplayLine();
      break;
    }

    for (const HighsInt i : search_indices) {
      if (stallNodesPerWorker[i] == 0) {
        numStallNodes = 0;
      } else {
        numStallNodes += stallNodesPerWorker[i];
      }
      stallNodesPerWorker[i] = 0;
      if (options_mip_->mip_max_stall_nodes != kHighsIInf &&
          numStallNodes >=
              std::max(HighsInt{1}, options_mip_->mip_max_stall_nodes)) {
        limit_reached = true;
        modelstatus_ = HighsModelStatus::kSolutionLimit;
        break;
      }
    }
    if (limit_reached) break;

    assert(!nodesInstalled());

    // Sync global information
    profiling_->start(kMipClockDomainPropagate);
    syncPools(search_indices);
    syncGlobalDomain(search_indices);
    mipdata_->getDomain().propagate();
    profiling_->stop(kMipClockDomainPropagate);

    profiling_->start(kMipClockPruneInfeasibleNodes);
    mipdata_->pruned_treeweight += mipdata_->nodequeue.pruneInfeasibleNodes(
        mipdata_->getDomain(), mipdata_->feastol);
    profiling_->stop(kMipClockPruneInfeasibleNodes);

    if (mipdata_->getDomain().infeasible()) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
      mipdata_->printDisplayLine();
      break;
    }

    mipdata_->updateLowerBound(std::min(
        mipdata_->upper_bound, mipdata_->nodequeue.getBestLowerBound()));
    mipdata_->printDisplayLine();

    if (mipdata_->nodequeue.empty()) break;

    // Reset global domain and sync worker's global domains.
    const bool spawn_more_workers =
        num_workers < max_num_workers &&
        mipdata_->nodequeue.numNodes() > num_workers;
    resetGlobalDomain(spawn_more_workers, mipdata_->hasMultipleWorkers());

    for (HighsMipWorker& worker : mipdata_->workers) {
      if (worker.getGlobalDomain().infeasible()) {
        infeasible = true;
        break;
      }
    }
    if (infeasible) {
      mipdata_->nodequeue.clear();
      mipdata_->pruned_treeweight = 1.0;
      mipdata_->updateLowerBound(std::min(kHighsInf, mipdata_->upper_bound));
      mipdata_->printDisplayLine();
      break;
    }

    if (nodeLim == maxNodesPerWorkerLim &&
        checkWorkerRestartVotes(workerRestartVotes))
      goto restart;

    if (spawn_more_workers) {
      HighsInt new_max_num_workers =
          std::min(static_cast<HighsInt>(mipdata_->nodequeue.numNodes()),
                   max_num_workers);
      mipdata_->getPseudoCost().removeChanged();
      if (num_workers == 1) {
        constructAdditionalWorkerData(master_worker);
      }
      createNewWorkers(new_max_num_workers - num_workers);
      num_workers = new_max_num_workers;
    }
  }
  syncSolutions();
  profiling_->stop(kMipClockSearch);

  cleanupSolve();
}

void HighsMipSolver::cleanupSolve() {
  for (HighsMipWorker& worker : mipdata_->workers) {
    assert(worker.solutions_.empty());
  }
  if (mipdata_->terminatorActive()) {
    if (mipdata_->terminatorTerminated()) {
      // Indicate that this instance has been interrupted
      modelstatus_ = HighsModelStatus::kHighsInterrupt;
    } else if (!submip) {
      // When sub-MIPs call cleanupSolve(), they generally don't have
      // a termination criterion for the whole MIP solver
      //
      // Possibly allow sub-MIPs to terminate if the time limit is
      // reached
      //
      // No other instance has terminated, so terminate
      mipdata_->terminatorTerminate();
    }
  }

  // Force a final logging line
  mipdata_->printDisplayLine(kSolutionSourceCleanup);
  // Stop the solve clock - which won't be running if presolve
  // determines the model status
  if (profiling_->running(kSolveTime)) profiling_->stop(kSolveTime);
  // Need to complete the calculation of P-D integral, checking for NO
  // gap change
  mipdata_->updatePrimalDualIntegral(
      mipdata_->lower_bound, mipdata_->lower_bound, mipdata_->upper_bound,
      mipdata_->upper_bound, false);
  profiling_->start(kPostsolveTime);

  bool havesolution = solution_objective_ != kHighsInf;
  bool feasible;
  if (havesolution) {
    feasible =
        bound_violation_ <= options_mip_->mip_feasibility_tolerance &&
        integrality_violation_ <= options_mip_->mip_feasibility_tolerance &&
        row_violation_ <= options_mip_->mip_feasibility_tolerance;
    if (feasible && mipdata_->upper_bound == kHighsInf) {
      mipdata_->checkAddSolution();
      if (modelstatus_ == HighsModelStatus::kInfeasible &&
          mipdata_->upper_bound != kHighsInf) {
        mipdata_->lower_bound = mipdata_->upper_bound;
      }
    }
  } else {
    feasible = false;
  }

  dual_bound_ = mipdata_->lower_bound;
  if (mipdata_->objectiveFunction.isIntegral()) {
    double rounded_lower_bound =
        std::ceil(mipdata_->lower_bound *
                      mipdata_->objectiveFunction.integralScale() -
                  mipdata_->feastol) /
        mipdata_->objectiveFunction.integralScale();
    dual_bound_ = std::max(dual_bound_, rounded_lower_bound);
  }
  dual_bound_ += model_->offset_;
  primal_bound_ = mipdata_->upper_bound + model_->offset_;
  node_count_ = mipdata_->num_nodes;
  total_lp_iterations_ = mipdata_->total_lp_iterations;
  dual_bound_ = std::min(dual_bound_, primal_bound_);
  primal_dual_integral_ = mipdata_->primal_dual_integral.value;

  // adjust objective sense in case of maximization problem
  if (orig_model_->sense_ == ObjSense::kMaximize) {
    dual_bound_ = -dual_bound_;
    primal_bound_ = -primal_bound_;
  }

  if (modelstatus_ == HighsModelStatus::kNotset ||
      modelstatus_ == HighsModelStatus::kInfeasible) {
    if (feasible && havesolution)
      modelstatus_ = HighsModelStatus::kOptimal;
    else
      modelstatus_ = HighsModelStatus::kInfeasible;
  }

  profiling_->stop(kPostsolveTime);
  timer_.stop();

  gap_ = fabs(primal_bound_ - dual_bound_);
  if (primal_bound_ == 0.0)
    gap_ = dual_bound_ == 0.0 ? 0.0 : kHighsInf;
  else if (primal_bound_ != kHighsInf)
    gap_ = fabs(primal_bound_ - dual_bound_) / fabs(primal_bound_);
  else
    gap_ = kHighsInf;

  //
  if (options_mip_->output_flag) {
    // solutionstatus is only used in solvingReport, but depends on
    // values of havesolution and feasible that are local to
    // HighsMipSolver::cleanupSolve
    const std::string solutionstatus =
        havesolution ? (feasible ? "feasible" : "infeasible") : "-";
    solvingReport(solutionstatus);
  }

  assert(modelstatus_ != HighsModelStatus::kNotset);

  if (improving_solution_file_ != nullptr) fclose(improving_solution_file_);
}

void HighsMipSolver::solvingReport(const std::string& solutionstatus) const {
  std::array<char, 128> gapString =
      getGapString(gap_, primal_bound_, options_mip_);

  bool timeless_log = options_mip_->timeless_log;
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "\nSolving report\n");
  if (this->orig_model_->model_name_.length())
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Model             %s\n",
                 this->orig_model_->model_name_.c_str());
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Status            %s\n",
               utilModelStatusToString(modelstatus_).c_str());
  if (modelstatus_ != HighsModelStatus::kSolveError &&
      modelstatus_ != HighsModelStatus::kMemoryLimit &&
      modelstatus_ != HighsModelStatus::kPresolveError) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Primal bound      %.12g\n"
                 "  Dual bound        %.12g\n"
                 "  Gap               %s\n",
                 primal_bound_, dual_bound_, gapString.data());
    if (!timeless_log)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "  P-D integral      %.2f\n",
                   mipdata_->primal_dual_integral.value);
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Solution status   %s\n", solutionstatus.c_str());
  if (solutionstatus != "-") {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "                    %.12g (objective)\n",
                 solution_objective_);
    if (bound_violation_)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "                    %.12g (variable bound violation)\n",
                   bound_violation_);
    if (row_violation_)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "                    %.12g (constraint bound violation)\n",
                   row_violation_);
    if (integrality_violation_)
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "                    %.12g (integrality violation)\n",
                   integrality_violation_);
  }
  if (!timeless_log) {
    auto callRecord = [&](HighsInt clock) {
      double mip_time = profiling_->read(clock, kMipRecord);
      double submip_time = profiling_->read(clock, kSubMipRecord);
      HighsInt mip_calls = profiling_->numCall(clock, kMipRecord);
      HighsInt submip_calls = profiling_->numCall(clock, kSubMipRecord);
      double total_time = mip_time + submip_time;
      // Only log postsolve if it's written as nonzero
      highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                   "                    %.2f (%s)\n", total_time,
                   profiling_->name[clock].c_str());
      if (mip_calls > 1 || submip_calls > 0) {
        highsLogUser(
            options_mip_->log_options, HighsLogType::kInfo,
            "                        MIP    time [calls] = %.2f [%d]\n",
            mip_time, int(mip_calls));
        if (submip_calls > 0)
          highsLogUser(
              options_mip_->log_options, HighsLogType::kInfo,
              "                        subMIP time [calls] = %.2f [%d]\n",
              submip_time, int(submip_calls));
      }
    };
    double total = timer_.read();
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Timing            %.2f\n", total);
    callRecord(kPresolveTime);
    callRecord(kSolveTime);
    callRecord(kPostsolveTime);
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  Max sub-MIP depth %d\n"
               "  Nodes             %llu\n",
               int(max_submip_level), (long long unsigned)mipdata_->num_nodes);
  if (mipdata_->total_repair_lp) {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Repair LPs        %llu (%llu feasible; %llu iterations)\n",
                 (long long unsigned)mipdata_->total_repair_lp,
                 (long long unsigned)mipdata_->total_repair_lp_feasible,
                 (long long unsigned)mipdata_->total_repair_lp_iterations);
  } else {
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "  Repair LPs        0\n");
  }
  highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
               "  LP iterations     %llu\n",
               (long long unsigned)mipdata_->total_lp_iterations);
  if (mipdata_->total_lp_iterations)
    highsLogUser(options_mip_->log_options, HighsLogType::kInfo,
                 "                    %llu (strong br.)\n"
                 "                    %llu (separation)\n"
                 "                    %llu (heuristics)\n",
                 (long long unsigned)mipdata_->sb_lp_iterations,
                 (long long unsigned)mipdata_->sepa_lp_iterations,
                 (long long unsigned)mipdata_->heuristic_lp_iterations);
}

// Only called in Highs::runPresolve
void HighsMipSolver::runMipPresolve(const HighsInt presolve_reduction_limit) {
  mipdata_ = decltype(mipdata_)(new HighsMipSolverData(*this));
  mipdata_->init();
  mipdata_->runMipPresolve(presolve_reduction_limit);
}

const HighsLp& HighsMipSolver::getPresolvedModel() const {
  return mipdata_->presolvedModel;
}

HighsPresolveStatus HighsMipSolver::getPresolveStatus() const {
  return mipdata_->presolve_status;
}

presolve::HighsPostsolveStack HighsMipSolver::getPostsolveStack() const {
  return mipdata_->postSolveStack;
}

void HighsMipSolver::callbackGetCutPool() const {
  assert(callback_->user_callback);
  assert(callback_->callbackActive(kCallbackMipGetCutPool));
  callback_->clearHighsCallbackOutput();
  HighsCallbackOutput& data_out = callback_->data_out;

  HighsSparseMatrix cut_matrix;
  mipdata_->getLp().getCutPool(data_out.cutpool_num_col,
                               data_out.cutpool_num_cut, data_out.cutpool_lower,
                               data_out.cutpool_upper, cut_matrix);

  // take ownership
  data_out.cutpool_start = std::move(cut_matrix.start_);
  data_out.cutpool_index = std::move(cut_matrix.index_);
  data_out.cutpool_value = std::move(cut_matrix.value_);

  const bool interrupt = mipdata_->interruptFromCallbackWithData(
      kCallbackMipGetCutPool, solution_objective_, "MIP cut pool");
  assert(!interrupt);
}

std::array<char, 128> getGapString(const double gap_,
                                   const double primal_bound_,
                                   const HighsOptions* options_mip_) {
  std::array<char, 128> gapString = {};
  if (gap_ == kHighsInf)
    std::strcpy(gapString.data(), "inf");
  else {
    auto gapValString = highsGapToString(gap_);
    double gapTol = options_mip_->mip_rel_gap;

    if (options_mip_->mip_abs_gap > options_mip_->mip_feasibility_tolerance) {
      gapTol = primal_bound_ == 0.0
                   ? kHighsInf
                   : std::max(gapTol,
                              options_mip_->mip_abs_gap / fabs(primal_bound_));
    }

    if (gapTol == 0.0)
      std::snprintf(gapString.data(), gapString.size(), "%s%%",
                    gapValString.data());
    else if (gapTol != kHighsInf) {
      auto gapTolString = highsGapToString(gapTol);
      std::snprintf(gapString.data(), gapString.size(),
                    "%s%% (tolerance: %s%%)", gapValString.data(),
                    gapTolString.data());
    } else
      std::snprintf(gapString.data(), gapString.size(), "%s%% (tolerance: inf)",
                    gapValString.data());
  }

  return gapString;
}

bool HighsMipSolver::solutionFeasible(const HighsLp* lp,
                                      const std::vector<double>& col_value,
                                      const std::vector<double>* pass_row_value,
                                      MipViolation& violation,
                                      HighsCDouble& obj) const {
  violation.clear();
  double& bound_violation = violation.bound_violation;
  double& integrality_violation = violation.integrality_violation;
  double& row_violation = violation.row_violation;
  HighsInt& num_bound_violations = violation.num_bound_violations;
  HighsInt& num_integrality_violations = violation.num_integrality_violations;
  HighsInt& num_row_violations = violation.num_row_violations;
  HighsInt& col_of_max_bound_violation = violation.col_of_max_bound_violation;
  HighsInt& col_of_max_integrality_violation =
      violation.col_of_max_integrality_violation;
  HighsInt& row_of_max_row_violation = violation.row_of_max_row_violation;
  assert(bound_violation == 0);
  assert(integrality_violation == 0);
  assert(row_violation == 0);
  assert(num_bound_violations == 0);
  assert(num_integrality_violations == 0);
  assert(num_row_violations == 0);
  assert(col_of_max_bound_violation == -1);
  assert(col_of_max_integrality_violation == -1);
  assert(row_of_max_row_violation == -1);
  const double mip_feasibility_tolerance =
      options_mip_->mip_feasibility_tolerance;

  obj = lp->offset_;

  if (kAllowDeveloperAssert)
    assert(col_value.size() == static_cast<size_t>(lp->num_col_));

  // Lambda to update violation record
  auto updateViolation = [&](double v, HighsInt index, HighsInt& num,
                             double& max_v, HighsInt& max_index) {
    if (v <= 0) return;
    if (v > mip_feasibility_tolerance) num++;
    if (v > max_v) {
      max_v = v;
      max_index = index;
    }
  };

  auto updatePrimalViolation = [&](double value, double lower, double upper,
                                   HighsInt index, HighsInt& num, double& max_v,
                                   HighsInt& max_index) {
    double primal_infeasibility;
    if (value < lower - mip_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + mip_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    } else {
      return;
    }
    updateViolation(primal_infeasibility, index, num, max_v, max_index);
  };

  auto updatePrimalViolationCol = [&](double value, HighsInt index) {
    updatePrimalViolation(value, lp->col_lower_[index], lp->col_upper_[index],
                          index, num_bound_violations, bound_violation,
                          col_of_max_bound_violation);
  };

  auto updatePrimalViolationRow = [&](double value, HighsInt index) {
    updatePrimalViolation(value, lp->row_lower_[index], lp->row_upper_[index],
                          index, num_row_violations, row_violation,
                          row_of_max_row_violation);
  };

  // Check column integrality and feasibility
  for (HighsInt i = 0; i != lp->num_col_; ++i) {
    const double value = col_value[i];
    obj += static_cast<HighsCDouble>(lp->col_cost_[i]) * value;
    if (lp->integrality_[i] == HighsVarType::kInteger)
      updateViolation(fractionality(value), i, num_integrality_violations,
                      integrality_violation, col_of_max_integrality_violation);
    updatePrimalViolationCol(value, i);
  }

  // Check row feasibility if there is a positive number of rows
  //
  // If there are no rows and pass_row_value is nullptr, then
  // row_value_p is also nullptr since row_value is not resized
  if (lp->num_row_ > 0) {
    std::vector<double> row_value;
    if (pass_row_value) {
      if (kAllowDeveloperAssert)
        assert((*pass_row_value).size() == static_cast<size_t>(lp->num_row_));
    } else {
      calculateRowValuesQuad(*lp, col_value, row_value);
    }
    const double* row_value_p =
        pass_row_value ? (*pass_row_value).data() : row_value.data();
    assert(row_value_p);

    for (HighsInt i = 0; i != lp->num_row_; ++i)
      updatePrimalViolationRow(row_value_p[i], i);
  }

  const bool feasible = bound_violation <= mip_feasibility_tolerance &&
                        integrality_violation <= mip_feasibility_tolerance &&
                        row_violation <= mip_feasibility_tolerance;
  return feasible;
}

std::vector<HighsModelStatus> HighsMipSolver::initialiseTerminatorRecord(
    HighsInt num_instance) const {
  std::vector<HighsModelStatus> record(num_instance, HighsModelStatus::kNotset);
  return record;
}

void HighsMipSolver::initialiseTerminator(HighsInt num_instance_,
                                          HighsInt my_instance_,
                                          HighsModelStatus* record_) {
  this->termination_status_ = HighsModelStatus::kNotset;
  this->terminator_.clear();
  this->terminator_.initialise(num_instance_, my_instance_, record_);
}

void HighsMipSolver::initialiseTerminator(const HighsMipSolver& mip_solver) {
  this->terminator_.clear();
  if (!mip_solver.mipdata_->terminatorActive()) return;
  assert(mip_solver.mipdata_->terminatorConcurrency() > 0);
  this->initialiseTerminator(mip_solver.mipdata_->terminatorConcurrency(),
                             mip_solver.mipdata_->terminatorMyInstance(),
                             mip_solver.terminator_.record);
}

void HighsMipSolver::setParallelLock(bool lock) const {
  if (!mipdata_->hasMultipleWorkers()) return;
  mipdata_->parallel_lock = lock;
  for (HighsConflictPool& conflictpool : mipdata_->conflictpools) {
    conflictpool.setAgeLock(lock);
  }
  mipdata_->cliquetable.setAllowParallel(!lock && !submip);
}

void HighsMipSolver::setProfiling(HighsProfiling* profiling) {
  assert(profiling);
  this->profiling_ = profiling;
}

void MipViolation::clear() {
  this->bound_violation = 0;
  this->integrality_violation = 0;
  this->row_violation = 0;
  this->num_bound_violations = 0;
  this->num_integrality_violations = 0;
  this->num_row_violations = 0;
  this->col_of_max_bound_violation = -1;
  this->col_of_max_integrality_violation = -1;
  this->row_of_max_row_violation = -1;
}

void MipViolation::copy(double& bound_violation_, double& row_violation_,
                        double& integrality_violation_) const {
  bound_violation_ = this->bound_violation;
  row_violation_ = this->row_violation;
  integrality_violation_ = this->integrality_violation;
}

void MipViolation::log(const HighsLogOptions& log_options,
                       const double objective_value,
                       const std::string& source) const {
  if (*(log_options.log_dev_level) > 0) {
    highsLogDev(log_options, HighsLogType::kWarning,
                "%s with objective %g has untransformed violations:\n",
                source.c_str(), objective_value);
    if (this->num_bound_violations)
      highsLogDev(log_options, HighsLogType::kWarning,
                  "   variable bound   (%d) with max of %.4g in column %d\n",
                  int(this->num_bound_violations), this->bound_violation,
                  int(this->col_of_max_bound_violation));
    if (this->num_row_violations)
      highsLogDev(log_options, HighsLogType::kWarning,
                  "   constraint bound (%d) with max of %.4g in row %d\n",
                  int(this->num_row_violations), this->row_violation,
                  int(this->row_of_max_row_violation));
    if (this->num_integrality_violations)
      highsLogDev(log_options, HighsLogType::kWarning,
                  "   integrality      (%d) with max of %.4g in column %d\n",
                  int(this->num_integrality_violations),
                  this->integrality_violation,
                  int(this->col_of_max_integrality_violation));
  } else {
    highsLogUser(
        log_options, HighsLogType::kWarning,
        "%s with objective %g has untransformed violations: "
        "variable bound = %.4g; constraint bound = %.4g; integrality = %.4g\n",
        source.c_str(), objective_value, this->bound_violation,
        this->row_violation, this->integrality_violation);
  }
}
