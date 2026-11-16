/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/restart.hpp
 * @brief
 */
// restart.hpp
#ifndef PDLP_HIPDLP_RESTART_HPP
#define PDLP_HIPDLP_RESTART_HPP

#include <cmath>
#include <iostream>
#include <vector>

#include "defs.hpp"
#include "scaling.hpp"
#include "solver_results.hpp"

// Struct to communicate restart decisions
struct RestartInfo {
  bool should_restart = false;
  bool restart_to_average =
      true;  // If true, restart to average; otherwise, to current

  RestartInfo(bool should_restart_, bool restart_to_average_)
      : should_restart(should_restart_),
        restart_to_average(restart_to_average_){};

  RestartInfo() = default;
};

class RestartScheme {
 public:
  RestartScheme() = default;
  void initialize(const SolverResults& results);

  // Checks if a restart should be performed based on the chosen strategy
  RestartInfo check(HighsInt current_iter, const SolverResults& current_results,
                    const SolverResults& average_results);

  HighsInt getLastRestartIter() const { return last_restart_iter_; }
  HighsInt last_restart_iter_ = 0;
  void passParams(PrimalDualParams* params) { params_ = params; };
  void passLogOptions(const HighsLogOptions* log_options) {
    log_options_ = log_options;
  };
  void passDebugPdlpLogFile(FILE* debug_pdlp_log_file) {
    debug_pdlp_log_file_ = debug_pdlp_log_file;
  };
  void passDebugPdlpData(DebugPdlpData* debug_pdlp_data) {
    debug_pdlp_data_ = debug_pdlp_data;
  };
  double getBeta() const { return beta_; };
  void updateBeta(double beta) {
    beta_ = beta;
    params_->omega = std::sqrt(beta);
  };
  void setLastRestartIter(int iter) { last_restart_iter_ = iter; };
  double getSufficientDecayFactor() const { return sufficient_decay_factor_; }
  double getNecessaryDecayFactor() const { return necessary_decay_factor_; }
  double getArtificialRestartThreshold() const {
    return artificial_restart_threshold_;
  }

  // State for adaptive restart
  // dPrimalFeasLastRestart = primal_feas_last_restart_
  // dDualFeasLastRestart = dual_feas_last_restart_;
  // dDualityGapLastRestart = duality_gap_last_restart_;
  double primal_feas_last_restart_ = 0.0;
  double dual_feas_last_restart_ = 0.0;
  double duality_gap_last_restart_ = 0.0;

 private:
  PrimalDualParams* params_;
  const HighsLogOptions* log_options_;
  FILE* debug_pdlp_log_file_ = nullptr;
  DebugPdlpData* debug_pdlp_data_;

  RestartStrategy strategy_ = RestartStrategy::NO_RESTART;
  HighsInt fixed_restart_interval_ = 100;
  double beta_;

  double primal_feas_last_candidate_ = 0.0;
  double dual_feas_last_candidate_ = 0.0;
  double duality_gap_last_candidate_ = 0.0;

  double sufficient_decay_factor_ = 0.2;
  double necessary_decay_factor_ = 0.8;
  double artificial_restart_threshold_ = 0.36;
};

#endif  // PDLP_HIPDLP_RESTART_HPP
