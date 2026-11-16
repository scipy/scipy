/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/restart.cc
 * @brief
 */
// restart.cc
#include "restart.hpp"

#include <algorithm>
#include <cmath>

#include "io/HighsIO.h"          // For pdlpLogging
#include "pdlp/cupdlp/cupdlp.h"  // For pdlpLogging

// Helper function to compute the score
static double compute_score(double beta, double primal_feas, double dual_feas,
                            double duality_gap) {
  return std::sqrt(beta * primal_feas * primal_feas +
                   dual_feas * dual_feas / beta + duality_gap * duality_gap);
}

// Initializes the restart scheme with parameters and initial results
void RestartScheme::initialize(const SolverResults& results) {
  strategy_ = params_->restart_strategy;
  fixed_restart_interval_ = params_->fixed_restart_interval;
  last_restart_iter_ = 0;

  // Stepsize ratio,
  //  β = dBeta = dDualStep / dPrimalStep,
  //    in the paper, primal weight is the ω:
  //    ω = √β
  beta_ = params_->omega * params_->omega;
}

RestartInfo RestartScheme::check(HighsInt current_iter,
                                 const SolverResults& current_results,
                                 const SolverResults& average_results) {
  if (strategy_ != RestartStrategy::ADAPTIVE_RESTART) {
    if (strategy_ == RestartStrategy::FIXED_RESTART &&
        current_iter - last_restart_iter_ >= fixed_restart_interval_) {
      return RestartInfo(true, true);
    }
    return RestartInfo(false, false);  // ✨ Corrected return
  }

  // Base Case: First check after a restart. Initialize state and exit.
  if (current_iter == last_restart_iter_) {
    primal_feas_last_restart_ = current_results.primal_feasibility;
    dual_feas_last_restart_ = current_results.dual_feasibility;
    duality_gap_last_restart_ = current_results.duality_gap;

    primal_feas_last_candidate_ = current_results.primal_feasibility;
    dual_feas_last_candidate_ = current_results.dual_feasibility;
    duality_gap_last_candidate_ = current_results.duality_gap;

    return {false, false};
  }

  // 1. Calculate scores and determine the best candidate
  double mu_current = compute_score(beta_, current_results.primal_feasibility,
                                    current_results.dual_feasibility,
                                    current_results.duality_gap);
  double mu_average = compute_score(beta_, average_results.primal_feasibility,
                                    average_results.dual_feasibility,
                                    average_results.duality_gap);
  debugPdlpRestartLog(debug_pdlp_log_file_, current_iter, mu_current,
                      mu_average);
  double candidate_score;
  bool restart_to_average = false;
  if (mu_current < mu_average) {
    candidate_score = mu_current;
    restart_to_average = false;
  } else {
    candidate_score = mu_average;
    restart_to_average = true;
  }

  // 2. Check the three restart conditions in order
  bool should_restart = false;
  PDHG_restart_choice restart_choice_for_logic = PDHG_NO_RESTART;

  if ((current_iter - last_restart_iter_) >=
      artificial_restart_threshold_ * current_iter) {
    // Condition 1: Artificial Restart
    should_restart = true;
  } else {
    double mu_last_restart =
        compute_score(beta_, primal_feas_last_restart_, dual_feas_last_restart_,
                      duality_gap_last_restart_);

    if (candidate_score < sufficient_decay_factor_ * mu_last_restart) {
      // Condition 2: Sufficient Decay
      should_restart = true;
    } else {
      double mu_last_candidate =
          compute_score(beta_, primal_feas_last_candidate_,
                        dual_feas_last_candidate_, duality_gap_last_candidate_);

      if (candidate_score < necessary_decay_factor_ * mu_last_restart &&
          candidate_score > mu_last_candidate) {
        // Condition 3: Necessary Decay
        should_restart = true;
      }
    }
  }

  // 3. ALWAYS update the "last candidate" metrics for the next check
  if (restart_to_average) {
    primal_feas_last_candidate_ = average_results.primal_feasibility;
    dual_feas_last_candidate_ = average_results.dual_feasibility;
    duality_gap_last_candidate_ = average_results.duality_gap;
  } else {
    primal_feas_last_candidate_ = current_results.primal_feasibility;
    dual_feas_last_candidate_ = current_results.dual_feasibility;
    duality_gap_last_candidate_ = current_results.duality_gap;
  }
  /*
    if (should_restart ) {
      if (restart_to_average){
        std::cout << "Last restart was iter " << last_restart_iter_ << ":
    average\n"; } else { std::cout << "Last restart was iter " <<
    last_restart_iter_ << ": current\n";
      }
    }
  */
  return RestartInfo(should_restart, restart_to_average);
}
