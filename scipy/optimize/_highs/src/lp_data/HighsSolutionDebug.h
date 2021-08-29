/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolutionDebug.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSSOLUTIONDEBUG_H_
#define SIMPLEX_HIGHSSOLUTIONDEBUG_H_

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"

struct HighsPrimalDualErrors {
  int num_nonzero_basic_duals;
  int num_large_nonzero_basic_duals;
  double max_nonzero_basic_dual;
  double sum_nonzero_basic_duals;
  int num_off_bound_nonbasic;
  double max_off_bound_nonbasic;
  double sum_off_bound_nonbasic;
  int num_primal_residual;
  double max_primal_residual;
  double sum_primal_residual;
  int num_dual_residual;
  double max_dual_residual;
  double sum_dual_residual;
};

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp lp, const HighsBasis& basis);

HighsDebugStatus debugSolutionRightSize(const HighsOptions& options,
                                        const HighsLp lp,
                                        const HighsSolution& solution);

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp lp,
                                      const HighsBasis& basis);

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsModelObject& highs_model_object);

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsBasis& basis, const HighsSolution& solution,
    const HighsInfo& info, const HighsModelStatus model_status);

HighsDebugStatus debugHighsBasicSolution(const string message,
                                         const HighsOptions& options,
                                         const HighsLp& lp,
                                         const HighsBasis& basis,
                                         const HighsSolution& solution);

HighsDebugStatus debugHighsBasicSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsBasis& basis, const HighsSolution& solution,
    const HighsSolutionParams& solution_params,
    const HighsModelStatus model_status);

// Methods below are not called externally

HighsDebugStatus debugHaveBasisAndSolutionData(const HighsLp& lp,
                                               const HighsBasis& basis,
                                               const HighsSolution& solution);

void debugHighsBasicSolutionPrimalDualInfeasibilitiesAndErrors(
    const HighsOptions& options, const HighsLp& lp, const HighsBasis& basis,
    const HighsSolution& solution, double& primal_objective_value,
    double& dual_objective_value, HighsSolutionParams& solution_params,
    HighsPrimalDualErrors& primal_dual_errors);

bool debugBasicSolutionVariable(
    bool report, const double primal_feasibility_tolerance,
    const double dual_feasibility_tolerance, const HighsBasisStatus status,
    const double lower, const double upper, const double value,
    const double dual, int& num_non_basic_var, int& num_basic_var,
    double& off_bound_nonbasic, double& primal_infeasibility,
    double& dual_infeasibility);

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors);

HighsDebugStatus debugCompareSolutionParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionObjectiveParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionStatusParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);
HighsDebugStatus debugCompareSolutionInfeasibilityParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1);

HighsDebugStatus debugCompareSolutionParamValue(const string name,
                                                const HighsOptions& options,
                                                const double v0,
                                                const double v1);

HighsDebugStatus debugCompareSolutionParamInteger(const string name,
                                                  const HighsOptions& options,
                                                  const int v0, const int v1);

void debugReportHighsBasicSolution(const string message,
                                   const HighsOptions& options,
                                   const HighsSolutionParams& solution_params,
                                   const HighsModelStatus model_status);

#endif  // SIMPLEX_HIGHSSOLUTIONDEBUG_H_
