/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/solver_results.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_SOLVER_RESULTS_HPP
#define PDLP_HIPDLP_SOLVER_RESULTS_HPP

#include <vector>

enum class TerminationStatus {
  NOTSET,
  OPTIMAL,
  INFEASIBLE,
  UNBOUNDED,
  MAXITER,
  TIMEOUT,
  ERROR
};

enum class TerminationIterate { LAST_ITERATE, AVERAGE_ITERATE };

struct SolverResults {
  TerminationStatus term_code = TerminationStatus::TIMEOUT;
  TerminationIterate term_iterate = TerminationIterate::LAST_ITERATE;

  double primal_obj = 0.0;
  double dual_obj = 0.0;
  double duality_gap = 0.0;
  double complementarity = 0.0;
  double primal_feasibility = 0.0;
  double dual_feasibility = 0.0;
  double relative_obj_gap = 0.0;

  // --- Averaged Metrics ---
  double primal_obj_average = 0.0;
  double dual_obj_average = 0.0;
  double duality_gap_average = 0.0;
  double primal_feasibility_average = 0.0;
  double dual_feasibility_average = 0.0;

  // --- Residual Vectors ---
  std::vector<double> primal_residual;
  std::vector<double> dual_residual;
  std::vector<double> primal_residual_average;
  std::vector<double> dual_residual_average;

  // --- Infeasibility Detection Metrics ---
  double primal_infeasibility_obj = 0.0;
  double dual_infeasibility_obj = 0.0;
  double primal_infeasibility_res = 0.0;
  double dual_infeasibility_res = 0.0;

  // --- Restart Statistics ---
  double primal_feasibility_last_restart = 0.0;
  double dual_feasibility_last_restart = 0.0;
  double duality_gap_last_restart = 0.0;
};

#endif  // PDLP_HIPDLP_SOLVER_RESULTS_HPP
