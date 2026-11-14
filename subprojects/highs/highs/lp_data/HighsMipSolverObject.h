/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsMipSolverObject.h
 * @brief Collection of class instances required to solve a MIP
 */
#ifndef LP_DATA_HIGHS_MIP_SOLVER_OBJECT_H_
#define LP_DATA_HIGHS_MIP_SOLVER_OBJECT_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsOptions.h"

class HighsMipSolverObject {
 public:
  HighsMipSolverObject(
      HighsLp& lp, HighsSolution& solution,
      std::vector<HighsObjectiveSolution>& saved_objective_and_solution,
      HighsInfo& highs_info, HighsCallback& callback, HighsOptions& options,
      HighsTimer& timer)
      : lp_(lp),
        solution_(solution),
        saved_objective_and_solution_(saved_objective_and_solution),
        highs_info_(highs_info),
        callback_(callback),
        options_(options),
        timer_(timer) {}

  HighsLp& lp_;
  HighsSolution& solution_;
  std::vector<HighsObjectiveSolution>& saved_objective_and_solution_;
  HighsInfo& highs_info_;
  HighsCallback& callback_;
  HighsOptions& options_;
  HighsTimer& timer_;
  HighsProfiling* profiling_ = nullptr;
  HighsModelStatus model_status_ = HighsModelStatus::kNotset;
  void setProfiling(HighsProfiling* profiling) {
    assert(profiling);
    profiling_ = profiling;
  }
};

#endif  // LP_DATA_HIGHS_MIP_SOLVER_OBJECT_H_
