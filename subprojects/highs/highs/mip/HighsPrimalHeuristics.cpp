/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsPrimalHeuristics.h"

#include <numeric>
#include <unordered_set>

#include "../extern/pdqsort/pdqsort.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLpUtils.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"
#include "util/HighsHash.h"
#include "util/HighsIntegers.h"

// GCC floating point errors are well-known for 32-bit architectures;
// see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=323.
// An easy workaround is to add the "volatile" keyword to avoid
// problematic GCC optimizations that impact precision.
#ifdef __i386__
#define FP_32BIT_VOLATILE volatile
#else
#define FP_32BIT_VOLATILE
#endif

HighsPrimalHeuristics::HighsPrimalHeuristics(HighsMipSolver& mipsolver)
    : mipsolver(mipsolver),
      successObservations(0.0),
      numSuccessObservations(0),
      infeasObservations(0.0),
      numInfeasObservations(0),
      randgen(mipsolver.options_mip_->random_seed) {}

bool HighsPrimalHeuristics::subMipShouldRetry(HighsMipWorker& worker,
                                              const HighsSearch& heur,
                                              double fixingrate,
                                              double& maxfixingrate,
                                              HighsInt& targetdepth) {
  int64_t new_lp_iterations =
      worker.getHeurLpIterations() + heur.getLocalLpIterations();
  if (new_lp_iterations + mipsolver.mipdata_->heuristic_lp_iterations >
      100000 + ((mipsolver.mipdata_->total_lp_iterations -
                 mipsolver.mipdata_->heuristic_lp_iterations -
                 mipsolver.mipdata_->sb_lp_iterations) >>
                1)) {
    worker.getHeurLpIterations() = new_lp_iterations;
    return false;
  }

  targetdepth = heur.getCurrentDepth() / 2;
  if (targetdepth <= 1 || worker.search_ptr_->checkLimits()) {
    worker.getHeurLpIterations() = new_lp_iterations;
    return false;
  }
  maxfixingrate = fixingrate * 0.5;
  return true;
}

void HighsPrimalHeuristics::setupIntCols() {
  intcols = mipsolver.mipdata_->integer_cols;

  pdqsort(intcols.begin(), intcols.end(), [&](HighsInt c1, HighsInt c2) {
    const FP_32BIT_VOLATILE double lockScore1 =
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->uplocks[c1]) *
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->downlocks[c1]);

    const FP_32BIT_VOLATILE double lockScore2 =
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->uplocks[c2]) *
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->downlocks[c2]);

    if (lockScore1 > lockScore2) return true;
    if (lockScore2 > lockScore1) return false;

    const FP_32BIT_VOLATILE double cliqueScore1 =
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c1, 1)) *
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c1, 0));

    const FP_32BIT_VOLATILE double cliqueScore2 =
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c2, 1)) *
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c2, 0));

    return std::make_tuple(cliqueScore1, HighsHashHelpers::hash(uint64_t(c1)),
                           c1) >
           std::make_tuple(cliqueScore2, HighsHashHelpers::hash(uint64_t(c2)),
                           c2);
  });
}

bool HighsPrimalHeuristics::solveSubMip(
    HighsMipWorker& worker, const HighsLp& lp, const HighsBasis& basis,
    double fixingRate, std::vector<double> colLower,
    std::vector<double> colUpper, HighsInt maxleaves, HighsInt maxnodes,
    HighsInt stallnodes) {
  HighsOptions submipoptions = *mipsolver.options_mip_;
  HighsLp submip = lp;

  // set bounds and restore integrality of the lp relaxation copy
  submip.col_lower_ = std::move(colLower);
  submip.col_upper_ = std::move(colUpper);
  submip.integrality_ = mipsolver.model_->integrality_;
  submip.offset_ = 0;

  // set limits
  submipoptions.mip_max_leaves = maxleaves;
  submipoptions.output_flag = false;

  const bool allow_submip_log = true;
  if (allow_submip_log && lp.num_col_ == -54 && lp.num_row_ == -172) {
    submipoptions.output_flag = true;
    if (mipsolver.profiling_->sub_solver_)
      printf(
          "HighsPrimalHeuristics::solveSubMip (%d, %d) with output_flag = %s\n",
          int(lp.num_col_), int(lp.num_row_),
          highsBoolToString(submipoptions.output_flag).c_str());
  }

  submipoptions.mip_max_nodes = maxnodes;
  submipoptions.mip_max_stall_nodes = stallnodes;
  submipoptions.mip_pscost_minreliable = 0;
  submipoptions.time_limit -= mipsolver.timer_.read();
  submipoptions.objective_bound = worker.upper_limit;

  if (!mipsolver.submip) {
    double curr_abs_gap = worker.upper_limit - mipsolver.mipdata_->lower_bound;

    if (curr_abs_gap == kHighsInf) {
      curr_abs_gap = fabs(mipsolver.mipdata_->lower_bound);
      if (curr_abs_gap == kHighsInf) curr_abs_gap = 0.0;
    }

    submipoptions.mip_rel_gap = 0.0;
    submipoptions.mip_abs_gap =
        mipsolver.mipdata_->feastol * std::max(curr_abs_gap, 1000.0);
  }

  // check if only root presolve is allowed
  if (submipoptions.mip_root_presolve_only)
    submipoptions.presolve = kHighsOffString;
  else
    submipoptions.presolve = kHighsOnString;
  submipoptions.mip_detect_symmetry = false;
  submipoptions.mip_heuristic_effort = 0.8;
  // setup solver and run it

  HighsSolution solution;
  solution.value_valid = false;
  solution.dual_valid = false;
  if (!mipsolver.submip && !mipsolver.mipdata_->parallelLockActive()) {
    mipsolver.profiling_->start(kMipClockSubMipSolve);
  }
  mipsolver.profiling_->solveCall("MIP", mipsolver.submip);
  // Create HighsMipSolver instance for sub-MIP
  HighsMipSolver submipsolver(*mipsolver.callback_, submipoptions, submip,
                              solution, true, mipsolver.submip_level + 1);
  // Initialise termination_status_ and propagate any terminator to
  // the sub-MIP
  submipsolver.initialiseTerminator(mipsolver);
  submipsolver.rootbasis = &basis;
  HighsPseudocostInitialization pscostinit(worker.getPseudocost(), 1);
  submipsolver.pscostinit = &pscostinit;
  submipsolver.clqtableinit = &mipsolver.mipdata_->cliquetable;
  submipsolver.implicinit = &mipsolver.mipdata_->implications;
  // Solve the sub-MIP
  //
  // Copy the pointer to global sub-solver data into the sub-MIP
  // solver
  submipsolver.setProfiling(mipsolver.profiling_);
  // Stop the solve timer so that presolve/solve/postsolve for the
  // sub-MIP are timed independently
  const bool was_running_solve = mipsolver.profiling_->running(kSolveTime);
  if (was_running_solve) mipsolver.profiling_->stop(kSolveTime);
  // Only start timing the submip if the calling MIP isn't a sub-MIP
  if (mipsolver.profiling_->sub_solver_)
    printf(
        "\nHighsPrimalHeuristics::solveSubMip Before run() for %sMIP at depth "
        "%2d on thread %2d\n",
        mipsolver.submip ? "sub-" : "    ", int(mipsolver.submip_level),
        int(mipsolver.profiling_->myThread()));
  if (!mipsolver.submip) mipsolver.profiling_->start(kSubSolverSubMip);
  // Ensure that sub-solver call time data accumulate in the sub-MIP record
  mipsolver.profiling_->setSubMip(true);
  submipsolver.run();
  if (mipsolver.profiling_->sub_solver_)
    printf(
        "HighsPrimalHeuristics::solveSubMip After  run() for %sMIP at depth "
        "%2d "
        "on thread %2d\n\n",
        mipsolver.submip ? "sub-" : "    ", int(mipsolver.submip_level),
        int(mipsolver.profiling_->myThread()));
  // Ensure that further sub-solver call time data accumulate in the
  // MIP or sub-MIP record, according to whether the calling MIP is a
  // sub-MIP
  mipsolver.profiling_->setSubMip(mipsolver.submip);
  if (!mipsolver.submip) mipsolver.profiling_->stop(kSubSolverSubMip);
  worker.updateHeurStatsMaxSubMipLevel(submipsolver.max_submip_level + 1);
  if (!mipsolver.submip && !mipsolver.mipdata_->parallelLockActive()) {
    // Only stop timing the submip if the calling MIP isn't a sub-MIP
    mipsolver.profiling_->stop(kMipClockSubMipSolve);
  }
  if (was_running_solve) {
    // Re-start the solve timer now that presolve/solve/postsolve for the
    // sub-MIP have been timed independently
    const bool restart = true;
    mipsolver.profiling_->start(kSolveTime, restart);
  }
  // 22/07/25: Seems impossible for submipsolver.mipdata_ to be a null
  // pointer after calling HighsMipSolver::run(), and assert isn't
  // triggered for anything in ctest, but use direct test of
  // submipsolver.termination_status_, rather than
  // submipsolver.mipdata_.terminatorTerminated()
  if (!submipsolver.mipdata_) {
    printf(
        "HighsPrimalHeuristics::solveSubMip: submipsolver.mipdata_ is "
        "nullptr\n");
    assert(submipsolver.mipdata_);
  }
  if (submipsolver.termination_status_ != HighsModelStatus::kNotset) {
    worker.setHeurTerminationStatus(submipsolver.termination_status_);
    return false;
  }
  if (submipsolver.mipdata_) {
    double numUnfixed = mipsolver.mipdata_->integral_cols.size() +
                        mipsolver.mipdata_->continuous_cols.size();
    double adjustmentfactor =
        ((1 - fixingRate) * mipsolver.mipdata_->integral_cols.size() +
         mipsolver.mipdata_->continuous_cols.size()) /
        std::max(1.0, numUnfixed);
    // (double)mipsolver.orig_model_->a_matrix_.value_.size();
    int64_t adjusted_lp_iterations =
        (size_t)(adjustmentfactor * submipsolver.mipdata_->total_lp_iterations);
    worker.updateHeurStatsLpIters(
        adjusted_lp_iterations, submipsolver.mipdata_->total_repair_lp,
        submipsolver.mipdata_->total_repair_lp_feasible,
        submipsolver.mipdata_->total_repair_lp_iterations);
    // Warning: This will not be deterministic if sub-mips are run in parallel
    if (mipsolver.submip)
      mipsolver.mipdata_->num_nodes += std::max(
          int64_t{1}, int64_t(adjustmentfactor * submipsolver.node_count_));
  }

  if (submipsolver.modelstatus_ == HighsModelStatus::kInfeasible) {
    worker.updateHeurStatsInfeasObservations(fixingRate);
  }
  if (submipsolver.node_count_ <= 1 &&
      submipsolver.modelstatus_ == HighsModelStatus::kInfeasible)
    return false;
  double oldUpperLimit = worker.upper_limit;
  if (submipsolver.modelstatus_ != HighsModelStatus::kInfeasible &&
      !submipsolver.solution_.empty()) {
    trySolution(submipsolver.solution_, kSolutionSourceSubMip, worker);
  }

  if (worker.upper_limit < oldUpperLimit) {
    // remember fixing rate as good
    worker.updateHeurStatsSuccessObservations(fixingRate);
  }

  return true;
}

double HighsPrimalHeuristics::determineTargetFixingRate(
    HighsMipWorker& worker) {
  double lowFixingRate = 0.6;
  double highFixingRate = 0.6;

  HighsRandom& randgen =
      mipsolver.mipdata_->parallelLockActive() ? worker.randgen : this->randgen;

  if (getNumInfeasObservations(worker) != 0) {
    double infeasRate =
        getInfeasObservations(worker) / getNumInfeasObservations(worker);
    highFixingRate = 0.9 * infeasRate;
    lowFixingRate = std::min(lowFixingRate, highFixingRate);
  }

  if (getNumSuccessObservations(worker) != 0) {
    double successFixingRate =
        getSuccessObservations(worker) / getNumSuccessObservations(worker);
    lowFixingRate = std::min(lowFixingRate, 0.9 * successFixingRate);
    highFixingRate = std::max(successFixingRate * 1.1, highFixingRate);
  }

  double fixingRate = randgen.real(lowFixingRate, highFixingRate);
  // if (!mipsolver.submip) printf("fixing rate: %.2f\n", 100.0 * fixingRate);
  return fixingRate;
}

class HeuristicNeighbourhood {
  HighsDomain& localdom;
  HighsInt numFixed;
  HighsHashTable<HighsInt> fixedCols;
  size_t startCheckedChanges;
  size_t nCheckedChanges;
  HighsInt numTotal;

 public:
  HeuristicNeighbourhood(const HighsMipSolver& mipsolver, HighsDomain& localdom)
      : localdom(localdom),
        numFixed(0),
        startCheckedChanges(localdom.getDomainChangeStack().size()),
        nCheckedChanges(startCheckedChanges) {
    for (HighsInt i : mipsolver.mipdata_->integral_cols)
      if (localdom.col_lower_[i] == localdom.col_upper_[i]) ++numFixed;

    numTotal = mipsolver.mipdata_->integral_cols.size() - numFixed;
  }

  double getFixingRate() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      HighsInt col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (localdom.variableType(col) == HighsVarType::kContinuous) continue;
      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return numTotal ? static_cast<double>(fixedCols.size()) /
                          static_cast<double>(numTotal)
                    : 0.0;
  }

  void backtracked() {
    nCheckedChanges = startCheckedChanges;
    if (fixedCols.size()) fixedCols.clear();
  }
};

void HighsPrimalHeuristics::rootReducedCost(HighsMipWorker& worker) {
  std::vector<std::pair<double, HighsDomainChange>> lurkingBounds =
      mipsolver.mipdata_->redcostfixing.getLurkingBounds(
          mipsolver, worker.getGlobalDomain());
  if (10 * lurkingBounds.size() < mipsolver.mipdata_->integral_cols.size())
    return;
  pdqsort(lurkingBounds.begin(), lurkingBounds.end(),
          [](const std::pair<double, HighsDomainChange>& a,
             const std::pair<double, HighsDomainChange>& b) {
            return a.first > b.first;
          });

  HighsDomain localdom = worker.getGlobalDomain();

  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);

  double currCutoff = kHighsInf;
  double lower_bound =
      mipsolver.mipdata_->lower_bound + mipsolver.mipdata_->feastol;

  for (const std::pair<double, HighsDomainChange>& domchg : lurkingBounds) {
    currCutoff = domchg.first;

    if (currCutoff <= lower_bound) break;

    if (localdom.isActive(domchg.second)) continue;
    localdom.changeBound(domchg.second);

    while (true) {
      localdom.propagate();
      if (localdom.infeasible()) {
        localdom.conflictAnalysis(worker.getConflictPool(),
                                  worker.getGlobalDomain(),
                                  worker.getPseudocost());

        mipsolver.mipdata_->updateLowerBound(
            std::max(mipsolver.mipdata_->lower_bound, currCutoff));

        localdom.backtrack();
        if (localdom.getBranchDepth() == 0) break;
        neighbourhood.backtracked();
        continue;
      }
      break;
    }
    double fixingRate = neighbourhood.getFixingRate();
    if (fixingRate >= 0.5) break;
    // double gap = (currCutoff - mipsolver.mipdata_->lower_bound) /
    //             std::max(std::abs(mipsolver.mipdata_->lower_bound), 1.0);
    // if (gap < 0.001) break;
  }

  double fixingRate = neighbourhood.getFixingRate();
  if (fixingRate < 0.3) return;

  solveSubMip(worker, *mipsolver.model_, mipsolver.mipdata_->firstrootbasis,
              fixingRate, localdom.col_lower_, localdom.col_upper_,
              500,  // std::max(50, int(0.05 *
                    // (mipsolver.mipdata_->num_leaves))),
              200 + static_cast<HighsInt>(mipsolver.mipdata_->num_nodes / 20),
              12);
}

static double calcFixVal(double rootchange, double fracval, double cost) {
  // reinforce direction of this solution away from root
  // solution if the change is at least 0.4
  // otherwise take the direction where the objective gets worse
  // if objective is zero round to nearest integer
  if (rootchange >= 0.4)
    return std::ceil(fracval);
  else if (rootchange <= -0.4)
    return std::floor(fracval);
  else if (cost > 0.0)
    return std::ceil(fracval);
  else if (cost < 0.0)
    return std::floor(fracval);
  else
    return std::floor(fracval + 0.5);
}

static constexpr int64_t kDepthFirstNodeLimit = 1000;

void HighsPrimalHeuristics::RENS(HighsMipWorker& worker,
                                 const std::vector<double>& tmp) {
  // return if domain is infeasible
  if (worker.getGlobalDomain().infeasible()) return;

  HighsPseudocost pscost(worker.getPseudocost());
  HighsSearch heur(worker, pscost);

  HighsDomain& localdom = heur.getLocalDomain();
  heur.setHeuristic(true);

  std::vector<HighsInt> intcols_;
  if (mipsolver.mipdata_->parallelLockActive()) {
    intcols_ = intcols;
  }
  std::vector<HighsInt>& intcols =
      mipsolver.mipdata_->parallelLockActive() ? intcols_ : this->intcols;
  intcols.erase(std::remove_if(intcols.begin(), intcols.end(),
                               [&](HighsInt i) {
                                 return worker.getGlobalDomain().isFixed(i);
                               }),
                intcols.end());

  // LP relaxation instantiation
  HighsLpRelaxation heurlp(worker.getLpRelaxation());
  heurlp.setMipWorker(worker);
  heurlp.setProfiling(mipsolver.profiling_);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(worker.upper_limit);
  heurlp.setAdjustSymmetricBranchingCol(false);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.col_lower_.data(),
                                        localdom.col_upper_.data());
  localdom.clearChangedCols();
  heur.createNewNode();

  // determine the initial number of unfixed variables fixing rate to decide if
  // the problem is restricted enough to be considered for solving a submip
  double maxfixingrate = determineTargetFixingRate(worker);
  double fixingrate = 0.0;
  bool stop = false;
  // heurlp.setIterationLimit(2 * mipsolver.mipdata_->maxrootlpiters);
  // printf("iterlimit: %" HIGHSINT_FORMAT "\n",
  //       heurlp.getLpSolver().getOptions().simplex_iteration_limit);
  HighsInt targetdepth = 1;
  HighsInt nbacktracks = -1;
  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);
retry:
  ++nbacktracks;
  neighbourhood.backtracked();
  // printf("current depth : %" HIGHSINT_FORMAT
  //        "   target depth : %" HIGHSINT_FORMAT "\n",
  //        heur.getCurrentDepth(), targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) {
      worker.getHeurLpIterations() += heur.getLocalLpIterations();
      return;
    }
  }

  // printf("fixingrate before loop is %g\n", fixingrate);
  assert(heur.hasNode());
  while (true) {
    // printf("evaluating node\n");
    heur.evaluateNode();
    // printf("done evaluating node\n");
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      if (worker.getGlobalDomain().infeasible()) {
        worker.getHeurLpIterations() += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      neighbourhood.backtracked();
      continue;
    }

    fixingrate = neighbourhood.getFixingRate();
    // printf("after evaluating node current fixingrate is %g\n", fixingrate);
    if (fixingrate >= maxfixingrate) break;
    if (stop) break;
    if (nbacktracks >= 10) break;

    HighsInt numBranched = 0;
    double stopFixingRate = std::min(
        1.0 - (1.0 - neighbourhood.getFixingRate()) * 0.9, maxfixingrate);
    const auto& relaxationsol = heurlp.getSolution().col_value;
    for (HighsInt i : intcols) {
      if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

      double downval =
          std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
      double upval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);

      downval = std::min(downval, localdom.col_upper_[i]);
      upval = std::max(upval, localdom.col_lower_[i]);
      if (localdom.col_lower_[i] < downval) {
        ++numBranched;
        heur.branchUpwards(i, downval, downval - 0.5);
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          break;
        }
      }
      if (localdom.col_upper_[i] > upval) {
        ++numBranched;
        heur.branchDownwards(i, upval, upval + 0.5);
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          break;
        }
      }

      if (neighbourhood.getFixingRate() >= stopFixingRate) break;
    }

    if (numBranched == 0) {
      auto getFixVal = [&](HighsInt col, double fracval) {
        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objective is zero round to nearest integer
        double fixval =
            calcFixVal(mipsolver.mipdata_->rootlpsol.empty()
                           ? 0.0
                           : fracval - mipsolver.mipdata_->rootlpsol[col],
                       fracval, mipsolver.model_->col_cost_[col]);
        // make sure we do not set an infeasible domain
        fixval = std::min(localdom.col_upper_[col], fixval);
        fixval = std::max(localdom.col_lower_[col], fixval);
        return fixval;
      };

      pdqsort(heurlp.getFractionalIntegers().begin(),
              heurlp.getFractionalIntegers().end(),
              [&](const std::pair<HighsInt, double>& a,
                  const std::pair<HighsInt, double>& b) {
                return std::make_pair(
                           std::abs(getFixVal(a.first, a.second) - a.second),
                           HighsHashHelpers::hash(
                               (uint64_t(a.first) << 32) +
                               heurlp.getFractionalIntegers().size())) <
                       std::make_pair(
                           std::abs(getFixVal(b.first, b.second) - b.second),
                           HighsHashHelpers::hash(
                               (uint64_t(b.first) << 32) +
                               heurlp.getFractionalIntegers().size()));
              });

      double change = 0.0;
      // select a set of fractional variables to fix
      for (auto fracint : heurlp.getFractionalIntegers()) {
        double fixval = getFixVal(fracint.first, fracint.second);

        if (localdom.col_lower_[fracint.first] < fixval) {
          ++numBranched;
          heur.branchUpwards(fracint.first, fixval, fracint.second);
          localdom.propagate();
          if (localdom.infeasible()) {
            localdom.conflictAnalysis(worker.getConflictPool(),
                                      worker.getGlobalDomain(),
                                      worker.getPseudocost());
            break;
          }

          fixingrate = neighbourhood.getFixingRate();
        }

        if (localdom.col_upper_[fracint.first] > fixval) {
          ++numBranched;
          heur.branchDownwards(fracint.first, fixval, fracint.second);
          localdom.propagate();
          if (localdom.infeasible()) {
            localdom.conflictAnalysis(worker.getConflictPool(),
                                      worker.getGlobalDomain(),
                                      worker.getPseudocost());
            break;
          }

          fixingrate = neighbourhood.getFixingRate();
        }

        if (fixingrate >= maxfixingrate) break;

        change += std::abs(fixval - fracint.second);
        if (change >= 0.5) break;
      }
    }

    if (numBranched == 0) break;
    heurlp.flushDomain(localdom);
  }

  // printf("stopped heur dive with fixing rate %g\n", fixingrate);
  // if there is no node left it means we backtracked to the global domain and
  // the subproblem was solved with the dive
  if (!heur.hasNode()) {
    worker.getHeurLpIterations() += heur.getLocalLpIterations();
    return;
  }
  // determine the fixing rate to decide if the problem is restricted enough to
  // be considered for solving a submip

  fixingrate = neighbourhood.getFixingRate();
  // printf("fixing rate is %g\n", fixingrate);
  if (fixingrate < 0.1 ||
      (mipsolver.submip && mipsolver.mipdata_->numImprovingSols != 0)) {
    // heur.childselrule = ChildSelectionRule::kBestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10, kDepthFirstNodeLimit);
    worker.getHeurLpIterations() += heur.getLocalLpIterations();
    if (mipsolver.submip) mipsolver.mipdata_->num_nodes += heur.getLocalNodes();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  HighsInt node_reduction_factor =
      mipsolver.mipdata_->parallelLockActive()
          ? std::max(
                HighsInt{1},
                static_cast<HighsInt>(mipsolver.mipdata_->workers.size()) / 4)
          : 1;
  const bool solve_sub_mip_return = solveSubMip(
      worker, heurlp.getLp(), heurlp.getLpSolver().getBasis(), fixingrate,
      localdom.col_lower_, localdom.col_upper_,
      500,  // std::max(50, int(0.05 *
      // (mipsolver.mipdata_->num_leaves))),
      200 + mipsolver.mipdata_->num_nodes / (node_reduction_factor * 20), 12);
  if (worker.terminatorTerminated()) return;
  if (!solve_sub_mip_return) {
    if (subMipShouldRetry(worker, heur, fixingrate, maxfixingrate, targetdepth))
      goto retry;
    return;
  }

  worker.getHeurLpIterations() += heur.getLocalLpIterations();
}

void HighsPrimalHeuristics::RINS(HighsMipWorker& worker,
                                 const std::vector<double>& relaxationsol) {
  // return if domain is infeasible
  if (worker.getGlobalDomain().infeasible()) return;

  if (relaxationsol.size() != static_cast<size_t>(mipsolver.numCol())) return;

  std::vector<HighsInt> intcols_;
  if (mipsolver.mipdata_->parallelLockActive()) {
    intcols_ = intcols;
  }
  std::vector<HighsInt>& intcols =
      mipsolver.mipdata_->parallelLockActive() ? intcols_ : this->intcols;
  intcols.erase(std::remove_if(intcols.begin(), intcols.end(),
                               [&](HighsInt i) {
                                 return worker.getGlobalDomain().isFixed(i);
                               }),
                intcols.end());

  HighsPseudocost pscost(worker.getPseudocost());
  HighsSearch heur(worker, pscost);

  HighsDomain& localdom = heur.getLocalDomain();
  heur.setHeuristic(true);

  // LP relaxation instantiation
  HighsLpRelaxation heurlp(worker.getLpRelaxation());
  heurlp.setMipWorker(worker);
  heurlp.setProfiling(mipsolver.profiling_);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(worker.upper_limit);
  heurlp.setAdjustSymmetricBranchingCol(false);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.col_lower_.data(),
                                        localdom.col_upper_.data());
  localdom.clearChangedCols();
  heur.createNewNode();

  // determine the initial number of unfixed variables fixing rate to decide if
  // the problem is restricted enough to be considered for solving a submip
  double maxfixingrate = determineTargetFixingRate(worker);
  double minfixingrate = 0.25;
  double fixingrate = 0.0;
  bool stop = false;
  HighsInt nbacktracks = -1;
  HighsInt targetdepth = 1;
  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);
retry:
  ++nbacktracks;
  neighbourhood.backtracked();
  // printf("current depth : %" HIGHSINT_FORMAT "   target depth : %"
  // HIGHSINT_FORMAT "\n", heur.getCurrentDepth(),
  //       targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) {
      worker.getHeurLpIterations() += heur.getLocalLpIterations();
      return;
    }
  }

  assert(heur.hasNode());

  while (true) {
    heur.evaluateNode();
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      // printf("backtrack1\n");
      if (worker.getGlobalDomain().infeasible()) {
        worker.getHeurLpIterations() += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      neighbourhood.backtracked();
      continue;
    }

    fixingrate = neighbourhood.getFixingRate();

    if (stop) break;
    if (fixingrate >= maxfixingrate) break;
    if (nbacktracks >= 10) break;

    std::vector<std::pair<HighsInt, double>>::iterator fixcandend;

    // partition the fractional variables to consider which ones should we fix
    // in this dive first if there is an incumbent, we dive towards the RINS
    // neighbourhood
    fixcandend = std::partition(
        heurlp.getFractionalIntegers().begin(),
        heurlp.getFractionalIntegers().end(),
        [&](const std::pair<HighsInt, double>& fracvar) {
          return std::abs(relaxationsol[fracvar.first] -
                          mipsolver.mipdata_->incumbent[fracvar.first]) <=
                 mipsolver.mipdata_->feastol;
        });

    bool fixtolpsol = true;

    auto getFixVal = [&](HighsInt col, double fracval) {
      double fixval;
      if (fixtolpsol) {
        // RINS neighbourhood (with extension)
        fixval = std::floor(relaxationsol[col] + 0.5);
      } else {
        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objective is zero round to nearest integer
        fixval = calcFixVal(fracval - mipsolver.mipdata_->rootlpsol[col],
                            fracval, mipsolver.model_->col_cost_[col]);
      }
      // make sure we do not set an infeasible domain
      fixval = std::min(localdom.col_upper_[col], fixval);
      fixval = std::max(localdom.col_lower_[col], fixval);
      return fixval;
    };

    // no candidates left to fix for getting to the neighbourhood, therefore we
    // switch to a different diving strategy until the minimal fixing rate is
    // reached
    HighsInt numBranched = 0;
    if (heurlp.getFractionalIntegers().begin() == fixcandend) {
      fixingrate = neighbourhood.getFixingRate();
      double stopFixingRate =
          std::min(maxfixingrate, 1.0 - (1.0 - fixingrate) * 0.9);
      const auto& currlpsol = heurlp.getSolution().col_value;
      for (HighsInt i : intcols) {
        if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

        if (std::abs(currlpsol[i] - mipsolver.mipdata_->incumbent[i]) <=
            mipsolver.mipdata_->feastol) {
          double fixval = HighsIntegers::nearestInteger(currlpsol[i]);
          if (localdom.col_lower_[i] < fixval) {
            ++numBranched;
            heur.branchUpwards(i, fixval, fixval - 0.5);
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(worker.getConflictPool(),
                                        worker.getGlobalDomain(),
                                        worker.getPseudocost());
              break;
            }

            fixingrate = neighbourhood.getFixingRate();
          }
          if (localdom.col_upper_[i] > fixval) {
            ++numBranched;
            heur.branchDownwards(i, fixval, fixval + 0.5);
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(worker.getConflictPool(),
                                        worker.getGlobalDomain(),
                                        worker.getPseudocost());
              break;
            }

            fixingrate = neighbourhood.getFixingRate();
          }

          if (fixingrate >= stopFixingRate) break;
        }
      }

      if (numBranched != 0) {
        // printf(
        //    "fixed %" HIGHSINT_FORMAT " additional cols, old fixing rate:
        //    %.2f%%, new fixing " "rate: %.2f%%\n", numBranched, fixingrate,
        //    getFixingRate());
        heurlp.flushDomain(localdom);
        continue;
      }

      if (fixingrate >= minfixingrate)
        break;  // if the RINS neighbourhood achieved a high enough fixing rate
                // by itself we stop here
      fixcandend = heurlp.getFractionalIntegers().end();
      // now sort the variables by their distance towards the value they will
      // be fixed to
      fixtolpsol = false;
    }

    // now sort the variables by their distance towards the value they will be
    // fixed to
    pdqsort(heurlp.getFractionalIntegers().begin(), fixcandend,
            [&](const std::pair<HighsInt, double>& a,
                const std::pair<HighsInt, double>& b) {
              return std::make_pair(
                         std::abs(getFixVal(a.first, a.second) - a.second),
                         HighsHashHelpers::hash(
                             (uint64_t(a.first) << 32) +
                             heurlp.getFractionalIntegers().size())) <
                     std::make_pair(
                         std::abs(getFixVal(b.first, b.second) - b.second),
                         HighsHashHelpers::hash(
                             (uint64_t(b.first) << 32) +
                             heurlp.getFractionalIntegers().size()));
            });

    double change = 0.0;
    // select a set of fractional variables to fix
    for (auto fracint = heurlp.getFractionalIntegers().begin();
         fracint != fixcandend; ++fracint) {
      double fixval = getFixVal(fracint->first, fracint->second);

      if (localdom.col_lower_[fracint->first] < fixval) {
        ++numBranched;
        heur.branchUpwards(fracint->first, fixval, fracint->second);
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          break;
        }

        fixingrate = neighbourhood.getFixingRate();
      }

      if (localdom.col_upper_[fracint->first] > fixval) {
        ++numBranched;
        heur.branchDownwards(fracint->first, fixval, fracint->second);
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          break;
        }

        fixingrate = neighbourhood.getFixingRate();
      }

      if (fixingrate >= maxfixingrate) break;

      change += std::abs(fixval - fracint->second);
      if (change >= 0.5) break;
    }

    if (numBranched == 0) break;

    heurlp.flushDomain(localdom);

    // printf("%" HIGHSINT_FORMAT "/%" HIGHSINT_FORMAT " fixed, fixingrate is
    // %g\n", nfixed, ntotal, fixingrate);
  }

  // if there is no node left it means we backtracked to the global domain and
  // the subproblem was solved with the dive
  if (!heur.hasNode()) {
    worker.getHeurLpIterations() += heur.getLocalLpIterations();
    return;
  }
  // determine the fixing rate to decide if the problem is restricted enough
  // to be considered for solving a submip

  // printf("fixing rate is %g\n", fixingrate);
  fixingrate = neighbourhood.getFixingRate();
  if (fixingrate < 0.1 ||
      (mipsolver.submip && mipsolver.mipdata_->numImprovingSols != 0)) {
    // heur.childselrule = ChildSelectionRule::kBestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10, kDepthFirstNodeLimit);
    worker.getHeurLpIterations() += heur.getLocalLpIterations();
    if (mipsolver.submip) mipsolver.mipdata_->num_nodes += heur.getLocalNodes();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  HighsInt node_reduction_factor =
      mipsolver.mipdata_->parallelLockActive()
          ? std::max(
                HighsInt{1},
                static_cast<HighsInt>(mipsolver.mipdata_->workers.size()) / 4)
          : 1;
  const bool solve_sub_mip_return = solveSubMip(
      worker, heurlp.getLp(), heurlp.getLpSolver().getBasis(), fixingrate,
      localdom.col_lower_, localdom.col_upper_,
      500,  // std::max(50, int(0.05 *
      // (mipsolver.mipdata_->num_leaves))),
      200 + mipsolver.mipdata_->num_nodes / (node_reduction_factor * 20), 12);
  if (worker.terminatorTerminated()) return;
  if (!solve_sub_mip_return) {
    if (subMipShouldRetry(worker, heur, fixingrate, maxfixingrate, targetdepth))
      goto retry;
    return;
  }

  worker.getHeurLpIterations() += heur.getLocalLpIterations();
}

bool HighsPrimalHeuristics::tryRoundedPoint(HighsMipWorker& worker,
                                            const std::vector<double>& point,
                                            const int solution_source) {
  HighsDomain localdom = worker.getGlobalDomain();
  bool integerFeasible = true;

  HighsInt numintcols = intcols.size();
  for (HighsInt i = 0; i != numintcols; ++i) {
    HighsInt col = intcols[i];
    double intval = point[col];
    double rounded;
    // check if solution value of integer-constrained variable is actually
    // integral
    bool feasible =
        fractionality(intval, &rounded) <= mipsolver.mipdata_->feastol;
    integerFeasible = integerFeasible && feasible;
    if (!feasible) continue;
    // use rounded solution value and check against bounds
    intval = rounded;
    intval = std::min(localdom.col_upper_[col], intval);
    intval = std::max(localdom.col_lower_[col], intval);

    localdom.fixCol(col, intval, HighsDomain::Reason::branching());
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(worker.getConflictPool(),
                                worker.getGlobalDomain(),
                                worker.getPseudocost());
      return false;
    }
    localdom.propagate();
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(worker.getConflictPool(),
                                worker.getGlobalDomain(),
                                worker.getPseudocost());
      return false;
    }
  }

  if (numintcols != mipsolver.numCol()) {
    // LP relaxation instantiation
    HighsLpRelaxation lprelax(mipsolver);
    lprelax.setMipWorker(worker);
    lprelax.setProfiling(mipsolver.profiling_);
    lprelax.loadModel();
    lprelax.setIterationLimit(
        std::max(int64_t{10000}, 2 * mipsolver.mipdata_->firstrootlpiters));
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.col_lower_.data(),
                                           localdom.col_upper_.data());

    // check if only root presolve is allowed
    if (mipsolver.options_mip_->mip_root_presolve_only)
      lprelax.getLpSolver().setOptionValue("presolve", kHighsOffString);
    if (!mipsolver.options_mip_->mip_root_presolve_only &&
        (5 * numintcols) / mipsolver.numCol() >= 1)
      lprelax.getLpSolver().setOptionValue("presolve", kHighsOnString);
    else
      lprelax.getLpSolver().setBasis(mipsolver.mipdata_->firstrootbasis,
                                     "HighsPrimalHeuristics::tryRoundedPoint");

    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::kInfeasible) {
      std::vector<HighsInt> inds;
      std::vector<double> vals;
      double rhs;
      if (lprelax.computeDualInfProof(worker.getGlobalDomain(), inds, vals,
                                      rhs)) {
        HighsCutGeneration cutGen(lprelax, worker.getCutPool());
        cutGen.generateConflict(localdom, worker.getGlobalDomain(), inds, vals,
                                rhs);
      }
      return false;
    } else if (lprelax.unscaledPrimalFeasible(st)) {
      const auto& lpsol = lprelax.getLpSolver().getSolution().col_value;
      if (!integerFeasible) {
        // there may be fractional integer variables -> try ziRound heuristic
        ziRound(worker, lpsol);
        return trySolution(lpsol, solution_source, worker);
      } else {
        // all integer variables are fixed -> add incumbent
        addIncumbent(lpsol, lprelax.getObjective(), solution_source, worker);
        return true;
      }
    }
  }

  return trySolution(localdom.col_lower_, solution_source, worker);
}

bool HighsPrimalHeuristics::linesearchRounding(
    HighsMipWorker& worker, const std::vector<double>& point1,
    const std::vector<double>& point2, const int solution_source) {
  std::vector<double> roundedpoint;

  HighsInt numintcols = intcols.size();
  roundedpoint.resize(mipsolver.numCol());

  double alpha = 0.0;
  assert(int(mipsolver.mipdata_->uplocks.size()) == mipsolver.numCol());
  assert(int(point1.size()) == mipsolver.numCol());
  assert(int(point2.size()) == mipsolver.numCol());

  while (alpha < 1.0) {
    double nextalpha = 1.0;
    bool reachedpoint2 = true;
    // printf("trying alpha = %g\n", alpha);
    for (HighsInt i = 0; i != numintcols; ++i) {
      HighsInt col = intcols[i];
      assert(col >= 0);
      assert(col < mipsolver.numCol());
      if (mipsolver.mipdata_->uplocks[col] == 0) {
        roundedpoint[col] = std::ceil(std::max(point1[col], point2[col]) -
                                      mipsolver.mipdata_->feastol);
        continue;
      }

      if (mipsolver.mipdata_->downlocks[col] == 0) {
        roundedpoint[col] = std::floor(std::min(point1[col], point2[col]) +
                                       mipsolver.mipdata_->feastol);
        continue;
      }

      double convexcomb = (1.0 - alpha) * point1[col] + alpha * point2[col];
      double intpoint2 = std::floor(point2[col] + 0.5);
      roundedpoint[col] = std::floor(convexcomb + 0.5);

      if (roundedpoint[col] == intpoint2) continue;

      reachedpoint2 = false;
      double tmpalpha = (roundedpoint[col] + 0.5 + mipsolver.mipdata_->feastol -
                         point1[col]) /
                        std::abs(point2[col] - point1[col]);
      if (tmpalpha < nextalpha && tmpalpha > alpha + 1e-2) nextalpha = tmpalpha;
    }

    if (tryRoundedPoint(worker, roundedpoint, solution_source)) return true;

    if (reachedpoint2) return false;

    alpha = nextalpha;
  }

  return false;
}

void HighsPrimalHeuristics::randomizedRounding(
    HighsMipWorker& worker, const std::vector<double>& relaxationsol) {
  if (relaxationsol.size() != static_cast<size_t>(mipsolver.numCol())) return;

  HighsDomain localdom = worker.getGlobalDomain();
  HighsRandom& randgen =
      mipsolver.mipdata_->parallelLockActive() ? worker.randgen : this->randgen;

  for (HighsInt i : intcols) {
    double intval;
    if (mipsolver.mipdata_->uplocks[i] == 0)
      intval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);
    else if (mipsolver.mipdata_->downlocks[i] == 0)
      intval = std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
    else
      intval = std::floor(relaxationsol[i] + randgen.real(0.1, 0.9));

    intval = std::min(localdom.col_upper_[i], intval);
    intval = std::max(localdom.col_lower_[i], intval);

    localdom.fixCol(i, intval, HighsDomain::Reason::branching());
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(worker.getConflictPool(),
                                worker.getGlobalDomain(),
                                worker.getPseudocost());
      return;
    }
    localdom.propagate();
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(worker.getConflictPool(),
                                worker.getGlobalDomain(),
                                worker.getPseudocost());
      return;
    }
  }

  if (mipsolver.mipdata_->integer_cols.size() !=
      static_cast<size_t>(mipsolver.numCol())) {
    // LP relaxation instantiation
    HighsLpRelaxation lprelax(mipsolver);
    lprelax.setMipWorker(worker);
    lprelax.setProfiling(mipsolver.profiling_);
    lprelax.loadModel();
    lprelax.setIterationLimit(
        std::max(int64_t{10000}, 2 * mipsolver.mipdata_->firstrootlpiters));
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.col_lower_.data(),
                                           localdom.col_upper_.data());

    // check if only root presolve is allowed
    if (mipsolver.options_mip_->mip_root_presolve_only)
      lprelax.getLpSolver().setOptionValue("presolve", kHighsOffString);

    if (!mipsolver.options_mip_->mip_root_presolve_only &&
        (5 * intcols.size()) / mipsolver.numCol() >= 1) {
      // LP to solve is very much smaller, so use presolve rather than
      // the root basis
      lprelax.getLpSolver().setOptionValue("presolve", kHighsOnString);
    } else {
      lprelax.getLpSolver().setBasis(
          mipsolver.mipdata_->firstrootbasis,
          "HighsPrimalHeuristics::randomizedRounding");
    }

    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::kInfeasible) {
      std::vector<HighsInt> inds;
      std::vector<double> vals;
      double rhs;
      if (lprelax.computeDualInfProof(worker.getGlobalDomain(), inds, vals,
                                      rhs)) {
        HighsCutGeneration cutGen(lprelax, worker.getCutPool());
        cutGen.generateConflict(localdom, worker.getGlobalDomain(), inds, vals,
                                rhs);
      }

    } else if (HighsLpRelaxation::unscaledPrimalFeasible(st)) {
      addIncumbent(lprelax.getLpSolver().getSolution().col_value,
                   lprelax.getObjective(), kSolutionSourceRandomizedRounding,
                   worker);
    }
  } else {
    trySolution(localdom.col_lower_, kSolutionSourceRandomizedRounding, worker);
  }
}

void HighsPrimalHeuristics::shifting(HighsMipWorker& worker,
                                     const std::vector<double>& relaxationsol) {
  if (relaxationsol.size() != static_cast<size_t>(mipsolver.numCol())) return;

  std::vector<double> current_relax_solution = relaxationsol;
  HighsInt t = 0;
  const HighsLp& currentLp = *mipsolver.model_;
  // LP relaxation instantiation
  HighsLpRelaxation lprelax(worker.getLpRelaxation());
  lprelax.setMipWorker(worker);
  lprelax.setProfiling(mipsolver.profiling_);
  HighsRandom& randgen =
      mipsolver.mipdata_->parallelLockActive() ? worker.randgen : this->randgen;
  std::vector<std::pair<HighsInt, double>> current_fractional_integers =
      lprelax.getFractionalIntegers();
  std::vector<std::tuple<HighsInt, HighsInt, double>> current_infeasible_rows =
      mipsolver.mipdata_->getInfeasibleRows(current_relax_solution);
  size_t previous_infeasible_rows_size = current_infeasible_rows.size();
  bool hasInfeasibleConstraints = current_infeasible_rows.size() != 0;
  HighsInt iterationsWithoutReductions = 0;
  HighsInt maxIterationsWithoutReductions = 5;
  std::unordered_map<HighsInt, std::vector<HighsInt>> shift_iterations_set;
  std::vector<HighsInt> shifts;

  auto findPairByIndex = [](std::vector<std::pair<HighsInt, double>>& vec,
                            HighsInt k) {
    return std::find_if(
        vec.begin(), vec.end(),
        [k](const std::pair<HighsInt, double>& p) { return p.first == k; });
  };

  auto findShiftsByIndex =
      [](const std::unordered_map<HighsInt, std::vector<HighsInt>>& shifts,
         HighsInt k) -> std::vector<HighsInt> {
    auto it = shifts.find(k);
    if (it != shifts.end()) {
      return it->second;
    }
    return {};
  };

  while ((current_fractional_integers.size() > 0 || hasInfeasibleConstraints) &&
         iterationsWithoutReductions <= maxIterationsWithoutReductions &&
         t <= static_cast<HighsInt>(mipsolver.mipdata_->integer_cols.size())) {
    t++;
    bool fractionalIntegersReduced = false;
    iterationsWithoutReductions++;
    if (hasInfeasibleConstraints) {
      // find an infeasible row that has a non-zero coefficient a(row,col) where
      // col in currentFractInt
      bool fractionalIntegerFound = false;
      HighsInt rIndex = 0;
      while (!fractionalIntegerFound &&
             rIndex != static_cast<HighsInt>(current_infeasible_rows.size())) {
        HighsInt r = std::get<0>(current_infeasible_rows[rIndex]);
        HighsInt start = mipsolver.mipdata_->ARstart_[r];
        HighsInt end = mipsolver.mipdata_->ARstart_[r + 1];
        for (HighsInt jInd = start; jInd != end; ++jInd) {
          auto it = findPairByIndex(current_fractional_integers,
                                    mipsolver.mipdata_->ARindex_[jInd]);
          fractionalIntegerFound = it != current_fractional_integers.end();
          if (fractionalIntegerFound) break;
        }
        rIndex++;
      }
      HighsInt row_index = rIndex - 1;
      if (!fractionalIntegerFound) {
        // otherwise select a random infeasible row
        row_index = randgen.integer(current_infeasible_rows.size());
      }

      HighsInt row = std::get<0>(current_infeasible_rows[row_index]);
      HighsInt row_sense = std::get<1>(current_infeasible_rows[row_index]);
      double infeasibility = std::get<2>(current_infeasible_rows[row_index]);
      double score_min = kHighsInf;
      HighsInt j_min = std::numeric_limits<HighsInt>::max();
      double x_j_min = kHighsInf;
      double aij_min = 0.0;
      bool moveValueUp = false;
      HighsInt start = mipsolver.mipdata_->ARstart_[row];
      HighsInt end = mipsolver.mipdata_->ARstart_[row + 1];
      for (HighsInt jInd = start; jInd != end; ++jInd) {
        HighsInt j = mipsolver.mipdata_->ARindex_[jInd];

        // skip fixed variables
        if (currentLp.col_lower_[j] == currentLp.col_upper_[j]) continue;

        // lambda for finding best shift
        auto repair = [&findPairByIndex, &current_fractional_integers,
                       &findShiftsByIndex, &shift_iterations_set, &t,
                       &score_min, &j_min, &aij_min, &x_j_min,
                       &current_relax_solution, &moveValueUp](
                          HighsInt col, HighsInt direction, HighsInt row_sense,
                          HighsInt numLocks, double coef, double cost,
                          bool isMaximization, bool isInteger, bool isAtBound) {
          if ((row_sense < 0 || direction * coef > 0) &&
              (row_sense > 0 || direction * coef < 0))
            return;

          // skip variables at bounds
          if (isAtBound) return;

          // search for column
          auto it = findPairByIndex(current_fractional_integers, col);

          // add data
          bool found = it != current_fractional_integers.end();

          double score = kHighsInf;
          if (found) {
            score = -1.0 + 1.0 / static_cast<double>(numLocks + 1);
          } else {
            const auto& shifts = findShiftsByIndex(shift_iterations_set, col);
            if (shifts.empty())
              score = direction * (isMaximization ? -cost : cost);
            else {
              score = 0.0;
              for (double shift : shifts) {
                if (direction * shift > 0)
                  score += pow(1.1, direction * shift - t);
              }
            }
            if (isInteger) score += 1;
          }

          if (score < score_min) {
            score_min = score;
            j_min = col;
            aij_min = coef;
            x_j_min = current_relax_solution[col];
            moveValueUp = direction > 0;
          }
        };

        // repair with up-rounding
        repair(j, HighsInt{1}, row_sense, mipsolver.mipdata_->downlocks[j],
               mipsolver.mipdata_->ARvalue_[jInd], currentLp.col_cost_[j],
               mipsolver.orig_model_->sense_ == ObjSense::kMaximize,
               currentLp.integrality_[j] == HighsVarType::kInteger,
               std::abs(currentLp.col_upper_[j] - current_relax_solution[j]) <=
                   mipsolver.mipdata_->feastol);

        // repair with down-rounding
        repair(j, HighsInt{-1}, row_sense, mipsolver.mipdata_->uplocks[j],
               mipsolver.mipdata_->ARvalue_[jInd], currentLp.col_cost_[j],
               mipsolver.orig_model_->sense_ == ObjSense::kMaximize,
               currentLp.integrality_[j] == HighsVarType::kInteger,
               std::abs(current_relax_solution[j] - currentLp.col_lower_[j]) <=
                   mipsolver.mipdata_->feastol);
      }

      if (j_min != std::numeric_limits<HighsInt>::max()) {
        // Update current_fractional_integers
        auto it = findPairByIndex(current_fractional_integers, j_min);
        if (it != current_fractional_integers.end()) {
          current_fractional_integers.erase(it);
          fractionalIntegersReduced = true;
        }
        // Update current_relax_solution and shift_iterations_set (for not
        // fractional integers)
        if (moveValueUp) {
          if (fractionalIntegersReduced) {
            current_relax_solution[j_min] =
                std::ceil(x_j_min - mipsolver.mipdata_->feastol);
          } else {
            if (currentLp.integrality_[j_min] == HighsVarType::kInteger) {
              // variable is integer and not at the upper bound, so increment
              // by 1.
              current_relax_solution[j_min] = x_j_min + 1.0;
            } else {
              current_relax_solution[j_min] = std::min(
                  x_j_min + infeasibility / std::abs(aij_min),
                  currentLp.col_upper_[j_min] + mipsolver.mipdata_->feastol);
            }
            shift_iterations_set[j_min].push_back(t);
          }
        } else {
          if (fractionalIntegersReduced) {
            current_relax_solution[j_min] =
                std::floor(x_j_min + mipsolver.mipdata_->feastol);
          } else {
            if (currentLp.integrality_[j_min] == HighsVarType::kInteger) {
              // variable is integer and not at the lower bound, so decrement
              // by 1.
              current_relax_solution[j_min] = x_j_min - 1.0;
            } else {
              current_relax_solution[j_min] = std::max(
                  x_j_min - infeasibility / std::abs(aij_min),
                  currentLp.col_lower_[j_min] - mipsolver.mipdata_->feastol);
            }

            shift_iterations_set[j_min].push_back(-t);
          }
        }
      }
    } else {
      double xi_max = -1;
      double delta_c_min = kHighsInf;
      HighsInt pind_j_min = std::numeric_limits<HighsInt>::max();
      HighsInt j_min = std::numeric_limits<HighsInt>::max();
      double x_j_min = kHighsInf;
      HighsInt sigma = 0;
      for (HighsInt i = 0;
           i != static_cast<HighsInt>(current_fractional_integers.size());
           ++i) {
        std::pair<HighsInt, double> it = current_fractional_integers[i];
        HighsInt col = it.first;
        assert(col >= 0);
        assert(col < mipsolver.numCol());

        auto isBetter = [&currentLp, &it, &xi_max, &delta_c_min, &pind_j_min,
                         &j_min, &x_j_min, &sigma,
                         &i](double col, double xi, double roundedval,
                             HighsInt direction) {
          double c_min = currentLp.col_cost_[col] * (roundedval - it.second);
          if (xi > xi_max || (xi == xi_max && c_min < delta_c_min)) {
            xi_max = xi;
            delta_c_min = c_min;
            pind_j_min = i;
            j_min = col;
            x_j_min = roundedval;
            sigma = direction;
          }
        };

        isBetter(col, mipsolver.mipdata_->uplocks[col],
                 std::floor(it.second + mipsolver.mipdata_->feastol),
                 HighsInt{-1});
        isBetter(col, mipsolver.mipdata_->downlocks[col],
                 std::ceil(it.second - mipsolver.mipdata_->feastol),
                 HighsInt{1});
      }
      if (sigma != 0) {
        current_relax_solution[j_min] = x_j_min;
      }
      if (pind_j_min != std::numeric_limits<HighsInt>::max()) {
        current_fractional_integers.erase(current_fractional_integers.begin() +
                                          pind_j_min);
        fractionalIntegersReduced = true;
      }
    }
    current_infeasible_rows =
        mipsolver.mipdata_->getInfeasibleRows(current_relax_solution);
    hasInfeasibleConstraints = current_infeasible_rows.size() != 0;
    if (current_infeasible_rows.size() < previous_infeasible_rows_size ||
        fractionalIntegersReduced)
      iterationsWithoutReductions = 0;
    previous_infeasible_rows_size = current_infeasible_rows.size();
  }
  // re-check for feasibility and add incumbent
  if (hasInfeasibleConstraints) {
    tryRoundedPoint(worker, current_relax_solution, kSolutionSourceShifting);
  } else {
    if (current_fractional_integers.size() > 0) {
      ziRound(worker, current_relax_solution);
    } else {
      trySolution(current_relax_solution, kSolutionSourceShifting, worker);
    }
  }
}

void HighsPrimalHeuristics::ziRound(HighsMipWorker& worker,
                                    const std::vector<double>& relaxationsol) {
  // if (mipsolver.submip) return;
  if (relaxationsol.size() != static_cast<size_t>(mipsolver.numCol())) return;

  std::vector<double> current_relax_solution = relaxationsol;

  auto zi = [this](double x) {
    return std::min(std::ceil(x - mipsolver.mipdata_->feastol) - x,
                    x - std::floor(x + mipsolver.mipdata_->feastol));
  };

  // auto localdom = mipsolver.mipdata_->getDomain();

  HighsCDouble zi_total = 0.0;
  for (HighsInt i : intcols) {
    zi_total += zi(current_relax_solution[i]);
  }

  if (zi_total <= mipsolver.mipdata_->feastol) return;

  const HighsLp& currentLp = *mipsolver.model_;

  std::vector<double> rowActivities;
  std::vector<double> XrowLower;
  std::vector<double> XrowUpper;
  rowActivities.resize(currentLp.num_row_);
  XrowLower.resize(currentLp.num_row_);
  XrowUpper.resize(currentLp.num_row_);

  HighsInt loop_count = 0;
  HighsInt max_loop_count = 5;
  HighsCDouble previous_zi_total;
  HighsCDouble improvement_in_feasibility = kHighsInf;

  while (zi_total > mipsolver.mipdata_->feastol &&
         improvement_in_feasibility > mipsolver.mipdata_->feastol &&
         loop_count <= max_loop_count) {
    previous_zi_total = zi_total;
    loop_count++;

    if (currentLp.num_row_ > 0)
      getLpRowBounds(currentLp, 0, currentLp.num_row_ - 1, XrowLower.data(),
                     XrowUpper.data());

    for (HighsInt j : intcols) {
      double relax_solution = current_relax_solution[j];
      if (fractionality(relax_solution) <= mipsolver.mipdata_->feastol)
        continue;

      rowActivities.assign(currentLp.num_row_, 0.0);
      calculateRowValuesQuad(currentLp, current_relax_solution, rowActivities);

      double min_row_ratio_for_upper = kHighsInf;
      double min_row_ratio_for_lower = kHighsInf;

      for (HighsInt el = currentLp.a_matrix_.start_[j];
           el < currentLp.a_matrix_.start_[j + 1]; el++) {
        HighsInt i = currentLp.a_matrix_.index_[el];
        double aij = currentLp.a_matrix_.value_[el];

        double slack_upper = XrowUpper[i] - rowActivities[i];
        double slack_lower = rowActivities[i] - XrowLower[i];
        min_row_ratio_for_upper =
            std::min(min_row_ratio_for_upper,
                     (aij > 0 ? slack_upper : -slack_lower) / aij);
        min_row_ratio_for_lower =
            std::min(min_row_ratio_for_lower,
                     (aij > 0 ? slack_lower : -slack_upper) / aij);
      }

      double upper_bound = std::min(currentLp.col_upper_[j] - relax_solution,
                                    min_row_ratio_for_upper);
      double lower_bound = std::min(relax_solution - currentLp.col_lower_[j],
                                    min_row_ratio_for_lower);

      auto performUpdates = [&](HighsInt col, double change) {
        double old_relax_solution = current_relax_solution[col];
        current_relax_solution[col] += change;
        zi_total =
            zi_total - zi(old_relax_solution) + zi(current_relax_solution[col]);
      };

      if (std::abs(zi(relax_solution + upper_bound) -
                   zi(relax_solution - lower_bound)) <=
              mipsolver.mipdata_->feastol &&
          zi(relax_solution + upper_bound) < zi(relax_solution)) {
        double XcolCost = currentLp.col_cost_[j];
        bool ubObjChangeSmaller = XcolCost * (relax_solution + upper_bound) <=
                                  XcolCost * (relax_solution - lower_bound);
        bool isMinimization = currentLp.sense_ == ObjSense::kMinimize;
        if ((isMinimization && ubObjChangeSmaller) ||
            (!isMinimization && !ubObjChangeSmaller)) {
          performUpdates(j, upper_bound);
        } else {
          performUpdates(j, -lower_bound);
        }
      } else if (zi(relax_solution + upper_bound) <
                     zi(relax_solution - lower_bound) &&
                 zi(relax_solution + upper_bound) < zi(relax_solution)) {
        performUpdates(j, upper_bound);
      } else if (zi(relax_solution + upper_bound) >
                     zi(relax_solution - lower_bound) &&
                 zi(relax_solution - lower_bound) < zi(relax_solution)) {
        performUpdates(j, -lower_bound);
      }
    }
    improvement_in_feasibility = previous_zi_total - zi_total;
  }
  // re-check for feasibility and add incumbent
  trySolution(current_relax_solution, kSolutionSourceZiRound, worker);
}

void HighsPrimalHeuristics::feasibilityPump(HighsMipWorker& worker) {
  // LP relaxation instantiation
  HighsLpRelaxation lprelax(worker.getLpRelaxation());
  lprelax.setMipWorker(worker);
  lprelax.setProfiling(mipsolver.profiling_);
  std::unordered_set<std::vector<HighsInt>, HighsVectorHasher, HighsVectorEqual>
      referencepoints;
  std::vector<double> roundedsol;
  HighsLpRelaxation::Status status = lprelax.resolveLp();
  worker.getHeurLpIterations() += lprelax.getNumLpIterations();

  HighsRandom& randgen =
      mipsolver.mipdata_->parallelLockActive() ? worker.randgen : this->randgen;

  std::vector<double> fracintcost;
  std::vector<HighsInt> fracintset;

  std::vector<HighsInt> mask(mipsolver.numCol(), 1);
  std::vector<double> cost(mipsolver.numCol(), 0.0);

  lprelax.getLpSolver().setOptionValue("simplex_strategy",
                                       kSimplexStrategyPrimal);
  lprelax.setObjectiveLimit();
  lprelax.getLpSolver().setOptionValue(
      "primal_simplex_bound_perturbation_multiplier", 0.0);

  lprelax.setIterationLimit(5 * mipsolver.mipdata_->avgrootlpiters);

  while (!lprelax.getFractionalIntegers().empty()) {
    const auto& lpsol = lprelax.getLpSolver().getSolution().col_value;
    roundedsol = lprelax.getLpSolver().getSolution().col_value;

    std::vector<HighsInt> referencepoint;
    referencepoint.reserve(mipsolver.mipdata_->integer_cols.size());

    HighsDomain localdom = worker.getGlobalDomain();
    for (HighsInt i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.isColInteger(i));
      double intval = std::floor(roundedsol[i] + randgen.real(0.4, 0.6));
      intval = std::max(intval, localdom.col_lower_[i]);
      intval = std::min(intval, localdom.col_upper_[i]);
      roundedsol[i] = intval;
      referencepoint.push_back((HighsInt)intval);
      if (!localdom.infeasible()) {
        localdom.fixCol(i, intval, HighsDomain::Reason::branching());
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          continue;
        }
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(worker.getConflictPool(),
                                    worker.getGlobalDomain(),
                                    worker.getPseudocost());
          continue;
        }
      }
    }

    bool havecycle = !referencepoints.emplace(referencepoint).second;
    for (HighsInt k = 0; havecycle && k < 2; ++k) {
      for (HighsInt i = 0; i != 10; ++i) {
        HighsInt flippos =
            randgen.integer(mipsolver.mipdata_->integer_cols.size());
        HighsInt col = mipsolver.mipdata_->integer_cols[flippos];
        if (roundedsol[col] > lpsol[col])
          roundedsol[col] = (HighsInt)std::floor(lpsol[col]);
        else if (roundedsol[col] < lpsol[col])
          roundedsol[col] = (HighsInt)std::ceil(lpsol[col]);
        else if (roundedsol[col] < worker.getGlobalDomain().col_upper_[col])
          roundedsol[col] = worker.getGlobalDomain().col_upper_[col];
        else
          roundedsol[col] = worker.getGlobalDomain().col_lower_[col];

        referencepoint[flippos] = (HighsInt)roundedsol[col];
      }
      havecycle = !referencepoints.emplace(referencepoint).second;
    }

    if (havecycle) return;

    if (linesearchRounding(worker, lpsol, roundedsol,
                           kSolutionSourceFeasibilityPump))
      return;

    if (lprelax.getNumLpIterations() >=
        1000 + mipsolver.mipdata_->avgrootlpiters * 5)
      break;

    for (HighsInt i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.isColInteger(i));

      if (mipsolver.mipdata_->uplocks[i] == 0 ||
          mipsolver.mipdata_->downlocks[i] == 0)
        cost[i] = 0.0;
      else if (lpsol[i] > roundedsol[i] - mipsolver.mipdata_->feastol)
        cost[i] = -1.0 + randgen.real(-1e-4, 1e-4);
      else
        cost[i] = 1.0 + randgen.real(-1e-4, 1e-4);
    }

    lprelax.getLpSolver().changeColsCost(mask.data(), cost.data());
    int64_t niters = -lprelax.getNumLpIterations();
    status = lprelax.resolveLp();
    niters += lprelax.getNumLpIterations();
    if (niters == 0) break;
    worker.getHeurLpIterations() += niters;
  }

  if (lprelax.getFractionalIntegers().empty() &&
      HighsLpRelaxation::unscaledPrimalFeasible(status)) {
    addIncumbent(lprelax.getLpSolver().getSolution().col_value,
                 lprelax.getObjective(), kSolutionSourceFeasibilityPump,
                 worker);
  }
}

void HighsPrimalHeuristics::centralRounding(HighsMipWorker& worker) {
  if (mipsolver.mipdata_->analyticCenter.size() !=
      static_cast<size_t>(mipsolver.numCol()))
    return;

  if (!mipsolver.mipdata_->firstlpsol.empty())
    linesearchRounding(worker, mipsolver.mipdata_->firstlpsol,
                       mipsolver.mipdata_->analyticCenter,
                       kSolutionSourceCentralRounding);
  else if (!mipsolver.mipdata_->rootlpsol.empty())
    linesearchRounding(worker, mipsolver.mipdata_->rootlpsol,
                       mipsolver.mipdata_->analyticCenter,
                       kSolutionSourceCentralRounding);
  else
    linesearchRounding(worker, mipsolver.mipdata_->analyticCenter,
                       mipsolver.mipdata_->analyticCenter,
                       kSolutionSourceCentralRounding);
}

#if 0
void HighsPrimalHeuristics::clique() {
  HighsHashTable<HighsInt, double> entries;
  double offset = 0.0;

  HighsDomain& globaldom = mipsolver.mipdata_->getDomain();
  for (HighsInt j = 0; j != mipsolver.numCol(); ++j) {
    HighsInt col = j;
    double val = mipsolver.colCost(col);
    if (val == 0.0) continue;

    if (!globaldom.isBinary(col)) {
      offset += val * globaldom.col_lower_[col];
      continue;
    }

    mipsolver.mipdata_->cliquetable.resolveSubstitution(col, val, offset);
    entries[col] += val;
  }

  std::vector<double> profits;
  std::vector<HighsCliqueTable::CliqueVar> objvars;

  for (const auto& entry : entries) {
    double objprofit = -entry.value();
    if (objprofit < 0) {
      offset += objprofit;
      profits.push_back(-objprofit);
      objvars.emplace_back(entry.key(), 0);
    } else {
      profits.push_back(objprofit);
      objvars.emplace_back(entry.key(), 1);
    }
  }

  std::vector<double> solution(mipsolver.numCol());

  HighsInt nobjvars = profits.size();
  for (HighsInt i = 0; i != nobjvars; ++i) solution[objvars[i].col] = objvars[i].val;

  std::vector<std::vector<HighsCliqueTable::CliqueVar>> cliques;
  double bestviol;
  HighsInt bestviolpos;
  HighsInt numcliques;

  cliques = mipsolver.mipdata_->cliquetable.separateCliques(
      solution, mipsolver.mipdata_->getDomain(), mipsolver.mipdata_->feastol);
  numcliques = cliques.size();
  while (numcliques != 0) {
    bestviol = 0.5;
    bestviolpos = -1;

    for (HighsInt c = 0; c != numcliques; ++c) {
      double viol = -1.0;
      for (HighsCliqueTable::CliqueVar clqvar : cliques[c])
        viol += clqvar.weight(solution);

      if (viol > bestviolpos) {
        bestviolpos = c;
        bestviol = viol;
      }
    }

    cliques = mipsolver.mipdata_->cliquetable.separateCliques(
        solution, mipsolver.mipdata_->getDomain(), mipsolver.mipdata_->feastol);
    numcliques = cliques.size();
  }
}
#endif

bool HighsPrimalHeuristics::addIncumbent(const std::vector<double>& sol,
                                         double solobj,
                                         const int solution_source,
                                         HighsMipWorker& worker) {
  if (mipsolver.mipdata_->parallelLockActive()) {
    return worker.addIncumbent(sol, solobj, solution_source);
  } else {
    return mipsolver.mipdata_->addIncumbent(sol, solobj, solution_source);
  }
}

bool HighsPrimalHeuristics::trySolution(const std::vector<double>& solution,
                                        const int solution_source,
                                        HighsMipWorker& worker) {
  if (mipsolver.mipdata_->parallelLockActive()) {
    return worker.trySolution(solution, solution_source);
  } else {
    return mipsolver.mipdata_->trySolution(solution, solution_source);
  }
}

HighsInt HighsPrimalHeuristics::getNumSuccessObservations(
    HighsMipWorker& worker) const {
  return numSuccessObservations + worker.getHeurNumSuccessObservations();
}

HighsInt HighsPrimalHeuristics::getNumInfeasObservations(
    HighsMipWorker& worker) const {
  return numInfeasObservations + worker.getHeurNumInfeasObservations();
}

double HighsPrimalHeuristics::getSuccessObservations(
    HighsMipWorker& worker) const {
  return successObservations + worker.getHeurSuccessObservations();
}

double HighsPrimalHeuristics::getInfeasObservations(
    HighsMipWorker& worker) const {
  return infeasObservations + worker.getHeurInfeasObservations();
}

void HighsPrimalHeuristics::flushStatistics(HighsMipSolver& mipsolver,
                                            HighsMipWorker& worker) {
  int64_t total_repair_lp;
  int64_t total_repair_lp_feasible;
  int64_t total_repair_lp_iterations;
  int64_t lp_iterations;
  double successObservations;
  HighsInt numSuccessObservations;
  double infeasObservations;
  HighsInt numInfeasObservations;
  HighsInt max_submip_level;
  HighsModelStatus termination_status;
  worker.getHeurStatsValues(total_repair_lp, total_repair_lp_feasible,
                            total_repair_lp_iterations, lp_iterations,
                            successObservations, numSuccessObservations,
                            infeasObservations, numInfeasObservations,
                            max_submip_level, termination_status);

  mipsolver.mipdata_->total_repair_lp += total_repair_lp;
  mipsolver.mipdata_->total_repair_lp_feasible += total_repair_lp_feasible;
  mipsolver.mipdata_->total_repair_lp_iterations += total_repair_lp_iterations;
  mipsolver.mipdata_->heuristic_lp_iterations += lp_iterations;
  mipsolver.mipdata_->total_lp_iterations += lp_iterations;
  mipsolver.max_submip_level =
      std::max(mipsolver.max_submip_level, max_submip_level);
  if (termination_status != HighsModelStatus::kNotset &&
      mipsolver.termination_status_ == HighsModelStatus::kNotset) {
    mipsolver.termination_status_ = termination_status;
  }
  this->successObservations += successObservations;
  this->numSuccessObservations += numSuccessObservations;
  this->infeasObservations += infeasObservations;
  this->numInfeasObservations += numInfeasObservations;
  worker.resetHeurStats();
}
