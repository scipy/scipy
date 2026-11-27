/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available Hias open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_MIP_WORKER_H_
#define HIGHS_MIP_WORKER_H_

#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSeparation.h"

class HighsSearch;

class HighsMipWorker {
 private:
  struct SepaStatistics {
    SepaStatistics() : numNeighbourhoodQueries(0), sepa_lp_iterations(0) {}

    int64_t numNeighbourhoodQueries;
    int64_t sepa_lp_iterations;
  };
  struct HeurStatistics {
    HeurStatistics()
        : total_repair_lp(0),
          total_repair_lp_feasible(0),
          total_repair_lp_iterations(0),
          lp_iterations(0) {
      successObservations = 0;
      numSuccessObservations = 0;
      infeasObservations = 0;
      numInfeasObservations = 0;
      max_submip_level = 0;
      termination_status_ = HighsModelStatus::kNotset;
    }

    int64_t total_repair_lp;
    int64_t total_repair_lp_feasible;
    int64_t total_repair_lp_iterations;
    int64_t lp_iterations;

    double successObservations;
    HighsInt numSuccessObservations;
    double infeasObservations;
    HighsInt numInfeasObservations;
    HighsInt max_submip_level;
    HighsModelStatus termination_status_;
  };
  const HighsMipSolver& mipsolver_;
  const HighsMipSolverData& mipdata_;

  HighsLpRelaxation* lp_;
  HighsDomain* globaldom_;
  HighsCutPool* cutpool_;
  HighsConflictPool* conflictpool_;
  HighsPseudocost* pseudocost_;

  bool heuristics_allowed;
  HeurStatistics heur_stats;
  SepaStatistics sepa_stats;

 public:
  std::unique_ptr<HighsSearch> search_ptr_;
  std::unique_ptr<HighsSeparation> sepa_ptr_;
  std::vector<std::pair<HighsNodeQueue::OpenNode, bool>> processedNodes;
  std::vector<HighsNodeQueue::OpenNode> preparedNodes;
  size_t prepNodeIdx;

  double upper_bound;
  double upper_limit;
  double optimality_limit;
  bool early_termination;

  std::vector<std::tuple<std::vector<double>, double, int>> solutions_;

  HighsRandom randgen;

  const HighsMipSolver& getMipSolver() const;

  HighsMipWorker(const HighsMipSolver& mipsolver, HighsLpRelaxation* lp,
                 HighsDomain* domain, HighsCutPool* cutpool,
                 HighsConflictPool* conflictpool, HighsPseudocost* pseudocost);

  ~HighsMipWorker();

  void resetSearch();

  void resetSepa();

  HighsDomain& getGlobalDomain() const { return *globaldom_; }

  void setGlobalDomain(HighsDomain* globaldom) { globaldom_ = globaldom; }

  HighsPseudocost& getPseudocost() const { return *pseudocost_; }

  void setPseudocost(HighsPseudocost* pseudocost) { pseudocost_ = pseudocost; }

  HighsConflictPool& getConflictPool() const { return *conflictpool_; }

  void setConflictPool(HighsConflictPool* conflictpool) {
    conflictpool_ = conflictpool;
  }

  HighsCutPool& getCutPool() const { return *cutpool_; }

  void setCutPool(HighsCutPool* cutpool) { cutpool_ = cutpool; }

  HighsLpRelaxation& getLpRelaxation() const { return *lp_; }

  void setLpRelaxation(HighsLpRelaxation* lp) { lp_ = lp; }

  bool addIncumbent(const std::vector<double>& sol, double solobj,
                    int solution_source);

  std::pair<bool, double> transformNewIntegerFeasibleSolution(
      const std::vector<double>& sol) const;

  bool trySolution(const std::vector<double>& solution, int solution_source);

  void setAllowHeuristics(const bool allowed) { heuristics_allowed = allowed; }

  bool getAllowHeuristics() const { return heuristics_allowed; }

  double getOptimalityLimit() const;

  int64_t& getNumNeighbourhoodQueries() {
    return sepa_stats.numNeighbourhoodQueries;
  }

  int64_t& getSepaLpIterations() { return sepa_stats.sepa_lp_iterations; }

  void resetSepaStats();

  void updateHeurStatsLpIters(int64_t lp_iters, int64_t total_repair_lp,
                              int64_t total_repair_lp_feasible,
                              int64_t total_repair_lp_iters);

  void updateHeurStatsInfeasObservations(double fixingRate);

  void updateHeurStatsSuccessObservations(double fixingRate);

  int64_t& getHeurLpIterations() { return heur_stats.lp_iterations; }

  HighsInt& getHeurNumSuccessObservations() {
    return heur_stats.numSuccessObservations;
  }

  HighsInt& getHeurNumInfeasObservations() {
    return heur_stats.numInfeasObservations;
  }

  double& getHeurSuccessObservations() {
    return heur_stats.successObservations;
  }

  double& getHeurInfeasObservations() { return heur_stats.infeasObservations; }

  void setHeurTerminationStatus(HighsModelStatus status) {
    heur_stats.termination_status_ = status;
  }

  HighsModelStatus getHeurTerminationStatus() const {
    return heur_stats.termination_status_;
  }

  bool terminatorTerminated() const {
    return heur_stats.termination_status_ != HighsModelStatus::kNotset;
  }

  void updateHeurStatsMaxSubMipLevel(HighsInt new_level) {
    heur_stats.max_submip_level =
        std::max(new_level, heur_stats.max_submip_level);
  }

  void getHeurStatsValues(
      int64_t& total_repair_lp, int64_t& total_repair_lp_feasible,
      int64_t& total_repair_lp_iterations, int64_t& lp_iterations,
      double& successObservations, HighsInt& numSuccessObservations,
      double& infeasObservations, HighsInt& numInfeasObservations,
      HighsInt& max_submip_level, HighsModelStatus& termination_status) const;

  void resetHeurStats();
};

#endif
