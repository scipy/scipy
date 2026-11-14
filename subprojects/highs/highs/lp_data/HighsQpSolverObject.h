/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsQpSolverObject.h
 * @brief Collection of class instances required to solve a QP
 */
#ifndef LP_DATA_HIGHS_QP_SOLVER_OBJECT_H_
#define LP_DATA_HIGHS_QP_SOLVER_OBJECT_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsOptions.h"
#include "model/HighsModel.h"

class HighsQpSolverObject {
 public:
  HighsQpSolverObject(HighsModel& model, HighsBasis& basis,
                      HighsSolution& solution, HighsInfo& highs_info,
                      HighsCallback& callback, HighsOptions& options,
                      HighsTimer& timer)
      : model_(model),
        basis_(basis),
        solution_(solution),
        highs_info_(highs_info),
        callback_(callback),
        options_(options),
        timer_(timer) {}

  HighsModel& model_;
  HighsBasis& basis_;
  HighsSolution& solution_;
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

#endif  // LP_DATA_HIGHS_QP_SOLVER_OBJECT_H_
