/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplex.h"
#include "simplex/HSimplexDebug.h"

void reportSimplexPhaseIterations(const HighsModelObject& highs_model_object,
                                  const SimplexAlgorithm algorithm,
                                  const bool initialise) {
  if (highs_model_object.simplex_info_.run_quiet) return;
  static int iteration_count0 = 0;
  static int dual_phase1_iteration_count0 = 0;
  static int dual_phase2_iteration_count0 = 0;
  static int primal_phase1_iteration_count0 = 0;
  static int primal_phase2_iteration_count0 = 0;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const HighsOptions& options = highs_model_object.options_;
  if (initialise) {
    iteration_count0 = highs_model_object.iteration_counts_.simplex;
    dual_phase1_iteration_count0 = simplex_info.dual_phase1_iteration_count;
    dual_phase2_iteration_count0 = simplex_info.dual_phase2_iteration_count;
    primal_phase1_iteration_count0 = simplex_info.primal_phase1_iteration_count;
    primal_phase2_iteration_count0 = simplex_info.primal_phase2_iteration_count;
    return;
  }
  const int delta_iteration_count =
      highs_model_object.iteration_counts_.simplex - iteration_count0;
  const int delta_dual_phase1_iteration_count =
      simplex_info.dual_phase1_iteration_count - dual_phase1_iteration_count0;
  const int delta_dual_phase2_iteration_count =
      simplex_info.dual_phase2_iteration_count - dual_phase2_iteration_count0;
  const int delta_primal_phase1_iteration_count =
      simplex_info.primal_phase1_iteration_count -
      primal_phase1_iteration_count0;
  const int delta_primal_phase2_iteration_count =
      simplex_info.primal_phase2_iteration_count -
      primal_phase2_iteration_count0;

  if (delta_dual_phase1_iteration_count + delta_dual_phase2_iteration_count +
          delta_primal_phase1_iteration_count +
          delta_primal_phase2_iteration_count !=
      delta_iteration_count) {
    printf("Iteration total error %d + %d + %d + %d != %d\n",
           delta_dual_phase1_iteration_count, delta_dual_phase2_iteration_count,
           delta_primal_phase1_iteration_count,
           delta_primal_phase2_iteration_count, delta_iteration_count);
  }
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "Primal simplex iterations [Ph1 %d; Ph2 %d] Total %d",
                    delta_primal_phase1_iteration_count,
                    delta_primal_phase2_iteration_count, delta_iteration_count);
  } else {
    HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "Dual simplex iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d",
                    delta_dual_phase1_iteration_count,
                    delta_dual_phase2_iteration_count,
                    delta_primal_phase2_iteration_count, delta_iteration_count);
  }
}
