/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_MIP_SOLVER_DATA_H_
#define HIGHS_MIP_SOLVER_DATA_H_

#include <atomic>
#include <vector>

#include "mip/HighsCliqueTable.h"
#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsDebugSol.h"
#include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipWorker.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsObjectiveFunction.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsRedcostFixing.h"
#include "mip/HighsSearch.h"
#include "mip/HighsSeparation.h"
#include "parallel/HighsParallel.h"
#include "presolve/HighsPostsolveStack.h"
#include "presolve/HighsSymmetry.h"
#include "util/HighsTimer.h"

struct HighsPrimaDualIntegral {
  double value;
  double prev_lb;
  double prev_ub;
  double prev_gap;
  double prev_time;
  void initialise();
};

enum MipSolutionSource : int {
  kSolutionSourceNone = -1,
  kSolutionSourceMin = kSolutionSourceNone,
  //  kSolutionSourceInitial, // 0
  kSolutionSourceBranching,           // B
  kSolutionSourceCentralRounding,     // C
  kSolutionSourceFeasibilityPump,     // F
  kSolutionSourceHeuristic,           // H
  kSolutionSourceShifting,            // I
  kSolutionSourceFeasibilityJump,     // J
  kSolutionSourceSubMip,              // L
  kSolutionSourceEmptyMip,            // P
  kSolutionSourceRandomizedRounding,  // R
  kSolutionSourceSolveLp,             // S
  kSolutionSourceEvaluateNode,        // T
  kSolutionSourceUnbounded,           // U
  kSolutionSourceUserSolution,        // X
  kSolutionSourceHighsSolution,       // Y
  kSolutionSourceZiRound,             // Z
  kSolutionSourceTrivialL,            // l
  kSolutionSourceTrivialP,            // p
  kSolutionSourceTrivialU,            // u
  kSolutionSourceTrivialZ,            // z
  kSolutionSourceCleanup,
  kSolutionSourceCount
};

struct HighsMipSolverData {
  HighsMipSolver& mipsolver;

  std::deque<HighsLpRelaxation> lps;
  std::deque<HighsCutPool> cutpools;
  std::deque<HighsConflictPool> conflictpools;
  std::deque<HighsDomain> domains;
  std::deque<HighsPseudocost> pseudocosts;
  std::deque<HighsMipWorker> workers;
  bool parallel_lock;
  std::atomic<int64_t> worker_lp_iterations_stop;

  HighsPrimalHeuristics heuristics;
  HighsCliqueTable cliquetable;
  HighsImplications implications;
  HighsRedcostFixing redcostfixing;
  HighsObjectiveFunction objectiveFunction;
  presolve::HighsPostsolveStack postSolveStack;
  HighsPresolveStatus presolve_status;
  HighsLp presolvedModel;
  bool cliquesExtracted;
  bool rowMatrixSet;
  bool analyticCenterComputed;
  HighsModelStatus analyticCenterStatus;
  bool detectSymmetries;
  HighsInt numRestarts;
  HighsInt numRestartsRoot;
  HighsInt numCliqueEntriesAfterPresolve;
  HighsInt numCliqueEntriesAfterFirstPresolve;

  std::vector<HighsInt> ARstart_;
  std::vector<HighsInt> ARindex_;
  std::vector<double> ARvalue_;
  std::vector<double> maxAbsRowCoef;
  std::vector<HighsBool> rowintegral;
  std::vector<HighsInt> uplocks;
  std::vector<HighsInt> downlocks;
  std::vector<HighsInt> integer_cols;
  std::vector<HighsInt> implint_cols;
  std::vector<HighsInt> integral_cols;
  std::vector<HighsInt> continuous_cols;

  HighsSymmetries symmetries;
  std::shared_ptr<const StabilizerOrbits> globalOrbits;

  double feastol;
  double epsilon;
  double heuristic_effort;
  int64_t dispfreq;
  std::vector<double> analyticCenter;
  std::vector<double> firstlpsol;
  std::vector<double> rootlpsol;
  double firstlpsolobj;
  HighsBasis firstrootbasis;
  double rootlpsolobj;
  HighsInt numintegercols;
  HighsInt maxTreeSizeLog2;

  HighsCDouble pruned_treeweight;
  double avgrootlpiters;
  double disptime;
  double last_disptime;
  int64_t firstrootlpiters;
  int64_t num_nodes;
  int64_t num_leaves;
  int64_t num_leaves_before_run;
  int64_t num_nodes_before_run;
  int64_t total_repair_lp;
  int64_t total_repair_lp_feasible;
  int64_t total_repair_lp_iterations;
  int64_t total_lp_iterations;
  int64_t heuristic_lp_iterations;
  int64_t sepa_lp_iterations;
  int64_t sb_lp_iterations;
  int64_t total_lp_iterations_before_run;
  int64_t heuristic_lp_iterations_before_run;
  int64_t sepa_lp_iterations_before_run;
  int64_t sb_lp_iterations_before_run;
  int64_t num_disp_lines;

  HighsInt numImprovingSols;
  double lower_bound;
  double upper_bound;
  double upper_limit;
  double optimality_limit;
  std::vector<double> incumbent;

  HighsNodeQueue nodequeue;

  HighsPrimaDualIntegral primal_dual_integral;

  HighsDebugSol debugSolution;

  HighsMipSolverData(HighsMipSolver& mipsolver);

  bool solutionRowFeasible(const std::vector<double>& solution) const;
  HighsModelStatus feasibilityJump();
  HighsModelStatus trivialHeuristics();

  void startAnalyticCenterComputation(
      const highs::parallel::TaskGroup& taskGroup);
  void finishAnalyticCenterComputation(
      const highs::parallel::TaskGroup& taskGroup);

  struct SymmetryDetectionData {
    HighsSymmetryDetection symDetection;
    HighsSymmetries symmetries;
    double detectionTime = 0.0;
  };

  void startSymmetryDetection(const highs::parallel::TaskGroup& taskGroup,
                              std::unique_ptr<SymmetryDetectionData>& symData);
  void finishSymmetryDetection(const highs::parallel::TaskGroup& taskGroup,
                               std::unique_ptr<SymmetryDetectionData>& symData);

  void updatePrimalDualIntegral(const double from_lower_bound,
                                const double to_lower_bound,
                                const double from_upper_bound,
                                const double to_upper_bound,
                                const bool check_bound_change = true,
                                const bool check_prev_data = true);
  double limitsToGap(const double use_lower_bound, const double use_upper_bound,
                     double& lb, double& ub) const;

  double computeNewUpperLimit(double upper_bound, double mip_abs_gap,
                              double mip_rel_gap) const;
  bool moreHeuristicsAllowed() const;
  void removeFixedIndices();
  void init();
  void basisTransfer();
  void checkObjIntegrality();
  void runMipPresolve(const HighsInt presolve_reduction_limit);
  void setupDomainPropagation();
  void saveReportMipSolution(const double new_upper_limit);
  void checkAddSolution();
  void runSetup();
  double transformNewIntegerFeasibleSolution(
      const std::vector<double>& sol,
      const bool possibly_store_as_new_incumbent = true);
  double percentageInactiveIntegers() const;
  void performRestart();
  bool checkSolution(const std::vector<double>& solution) const;
  std::vector<std::tuple<HighsInt, HighsInt, double>> getInfeasibleRows(
      const std::vector<double>& solution) const;
  bool trySolution(const std::vector<double>& solution,
                   const int solution_source = kSolutionSourceNone);
  bool rootSeparationRound(HighsMipWorker& worker, HighsSeparation& sepa,
                           HighsInt& ncuts, HighsLpRelaxation::Status& status);
  HighsLpRelaxation::Status evaluateRootLp(HighsMipWorker& worker);

  void evaluateRootNode(HighsMipWorker& worker);

  bool addIncumbent(const std::vector<double>& sol, double solobj,
                    const int solution_source,
                    const bool print_display_line = true,
                    const bool is_user_solution = false);

  const std::vector<double>& getSolution() const;

  std::string solutionSourceToString(const int solution_source,
                                     const bool code = true) const;
  void printSolutionSourceKey() const;
  void printDisplayLine(const int solution_source = kSolutionSourceNone);

  void getRow(HighsInt row, HighsInt& rowlen, const HighsInt*& rowinds,
              const double*& rowvals) const {
    HighsInt start = ARstart_[row];
    rowlen = ARstart_[row + 1] - start;
    rowinds = ARindex_.data() + start;
    rowvals = ARvalue_.data() + start;
  }

  bool checkLimits(int64_t nodeOffset = 0) const;
  void limitsToBounds(double& dual_bound, double& primal_bound,
                      double& mip_rel_gap) const;
  void updateLowerBound(double new_lower_bound,
                        const bool check_bound_change = true,
                        const bool check_prev_data = true);
  void setCallbackDataOut(const double mipsolver_objective_value) const;
  bool interruptFromCallbackWithData(const int callback_type,
                                     const double mipsolver_objective_value,
                                     const std::string& message = "") const;
  void queryExternalSolution(
      const double mipsolver_objective_value,
      const ExternalMipSolutionQueryOrigin external_solution_query_origin);

  HighsInt terminatorConcurrency() const;
  bool terminatorActive() const { return terminatorConcurrency() > 0; }
  HighsInt terminatorMyInstance() const;
  void terminatorTerminate();
  bool terminatorTerminated() const;
  void terminatorReport() const;

  bool parallelLockActive() const {
    return (parallel_lock && hasMultipleWorkers());
  }

  void updateWorkerEarlyTermination(HighsMipWorker& worker) {
    int64_t current = worker_lp_iterations_stop.load(std::memory_order_relaxed);
    const int64_t candidate = worker.search_ptr_->lpiterations;
    while (candidate < current &&
           !worker_lp_iterations_stop.compare_exchange_weak(
               current, candidate, std::memory_order_relaxed,
               std::memory_order_relaxed)) {
    }
    worker.early_termination = true;
  }

  bool hasMultipleWorkers() const { return workers.size() > 1; }

  HighsDomain& getDomain() { return domains[0]; }
  HighsConflictPool& getConflictPool() { return conflictpools[0]; }
  HighsCutPool& getCutPool() { return cutpools[0]; }
  HighsLpRelaxation& getLp() { return lps[0]; }
  HighsPseudocost& getPseudoCost() { return pseudocosts[0]; }
  const HighsDomain& getDomain() const { return domains[0]; }
  const HighsConflictPool& getConflictPool() const { return conflictpools[0]; }
  const HighsCutPool& getCutPool() const { return cutpools[0]; }
  const HighsLpRelaxation& getLp() const { return lps[0]; }
  const HighsPseudocost& getPseudoCost() const { return pseudocosts[0]; }
  void reportOriginalPresolvedCol(
      const HighsInt original_col,
      const std::vector<double>& presolved_solution);
};

#endif
