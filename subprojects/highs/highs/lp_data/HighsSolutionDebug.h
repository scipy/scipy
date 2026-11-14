/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.h
 * @brief
 */
#ifndef SIMPLEX_HIGHSSOLUTIONDEBUG_H_
#define SIMPLEX_HIGHSSOLUTIONDEBUG_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"

HighsDebugStatus debugHighsLpSolution(const std::string& message,
                                      const HighsLpSolverObject& solver_object);

HighsDebugStatus debugHighsSolution(const string& message,
                                    const HighsOptions& options,
                                    const HighsModel& model,
                                    const HighsSolution& solution,
                                    const HighsBasis& basis);

HighsDebugStatus debugHighsSolution(
    const string& message, const HighsOptions& options, const HighsModel& model,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status, const HighsInfo& info);

HighsDebugStatus debugHighsSolution(
    const std::string& message, const HighsOptions& options, const HighsLp& lp,
    const HighsHessian& hessian, const HighsSolution& solution,
    const HighsBasis& basis, const HighsModelStatus model_status,
    const HighsInfo& highs_info, const bool check_model_status_and_highs_info);

void debugReportHighsSolution(const string& message,
                              const HighsLogOptions& log_options,
                              const HighsInfo& highs_info,
                              const HighsModelStatus model_status);

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp& lp,
                                     const HighsBasis& basis);

HighsDebugStatus debugPrimalSolutionRightSize(const HighsOptions& options,
                                              const HighsLp& lp,
                                              const HighsSolution& solution);

HighsDebugStatus debugDualSolutionRightSize(const HighsOptions& options,
                                            const HighsLp& lp,
                                            const HighsSolution& solution);

HighsDebugStatus debugHighsBasisConsistent(const HighsOptions& options,
                                           const HighsLp& lp,
                                           const HighsBasis& basis);

// Methods below are not called externally

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors);

HighsDebugStatus debugCompareHighsInfo(const HighsOptions& options,
                                       const HighsInfo& highs_info0,
                                       const HighsInfo& highs_info1);
HighsDebugStatus debugCompareHighsInfoObjective(const HighsOptions& options,
                                                const HighsInfo& highs_info0,
                                                const HighsInfo& highs_info1);
HighsDebugStatus debugCompareHighsInfoStatus(const HighsOptions& options,
                                             const HighsInfo& highs_info0,
                                             const HighsInfo& highs_info1);
HighsDebugStatus debugCompareHighsInfoInfeasibility(
    const HighsOptions& options, const HighsInfo& highs_info0,
    const HighsInfo& highs_info1);

HighsDebugStatus debugCompareHighsInfoDouble(const string& name,
                                             const HighsOptions& options,
                                             const double v0, const double v1);

HighsDebugStatus debugCompareHighsInfoInteger(const string& name,
                                              const HighsOptions& options,
                                              const HighsInt v0,
                                              const HighsInt v1);

#endif  // SIMPLEX_HIGHSSOLUTIONDEBUG_H_
