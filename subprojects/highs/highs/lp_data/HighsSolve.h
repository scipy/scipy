/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLVE_H_
#define LP_DATA_HIGHSSOLVE_H_

#include "lp_data/HighsMipSolverObject.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsQpSolverObject.h"

HighsStatus solveLp(HighsLpSolverObject& solver_object,
                    const std::string& message);
HighsStatus solveMip(HighsMipSolverObject& solver_object,
                     const std::string& message);
HighsStatus solveQp(HighsQpSolverObject& solver_object,
                    const std::string& message);

HighsStatus solveUnconstrainedLp(HighsLpSolverObject& solver_object);
HighsStatus solveUnconstrainedLp(const HighsOptions& options, const HighsLp& lp,
                                 HighsModelStatus& model_status,
                                 HighsInfo& highs_info, HighsSolution& solution,
                                 HighsBasis& basis);
void assessExcessiveObjectiveBoundScaling(const HighsOptions& options,
                                          const HighsModel& model,
                                          HighsUserScaleData& user_scale_data);
bool useHipo(const HighsOptions& options,
             const std::string& specific_solver_option, const HighsLp& lp,
             const bool logging = false);

HighsStatus checkOptimality(const std::string& solver_type,
                            const HighsOptions& options_,
                            const HighsInfo& info_,
                            HighsModelStatus& model_status_);
#endif  // LP_DATA_HIGHSSOLVE_H_
