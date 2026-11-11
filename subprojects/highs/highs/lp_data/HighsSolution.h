/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolution.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSSOLUTION_H_
#define LP_DATA_HIGHSSOLUTION_H_

#include <string>
#include <vector>

#include "io/HighsIO.h"
#include "lp_data/HStruct.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsStatus.h"
#include "model/HighsModel.h"

class HighsLp;
struct IpxSolution;
class HighsOptions;

using std::string;

// Collecting absolute and relative errors, and the corresponding
// indices is required for Glpsol output
struct HighsError {
  double absolute_value;
  HighsInt absolute_index;
  double relative_value;
  HighsInt relative_index;
  void print(const std::string& message);
  void reset();
  void invalidate();
};

// These are values used for HighsSolutionDebug, or for Glpsol output,
// so not worthy of being in HighsInfo
struct HighsPrimalDualErrors {
  HighsInt num_nonzero_basic_duals;
  double max_nonzero_basic_dual;
  double sum_nonzero_basic_duals;
  HighsInt num_off_bound_nonbasic;
  double max_off_bound_nonbasic;
  double sum_off_bound_nonbasic;
  HighsInt glpsol_num_primal_residual_errors;
  HighsInt glpsol_num_dual_residual_errors;
  HighsError glpsol_max_primal_residual;
  HighsError glpsol_max_primal_infeasibility;
  HighsError glpsol_max_dual_residual;
  HighsError glpsol_max_dual_infeasibility;
};

void getKktFailures(const HighsOptions& options, const HighsModel& model,
                    const HighsSolution& solution, const HighsBasis& basis,
                    HighsInfo& highs_info,
                    HighsPrimalDualErrors* primal_dual_errors = nullptr,
                    const bool get_residuals = false);

void getQpKktFailures(const HighsOptions& options, const HighsModel& model,
                      const HighsSolution& solution, HighsInfo& highs_info);

void getLpKktFailures(const HighsOptions& options, const HighsLp& lp,
                      const HighsSolution& solution, const HighsBasis& basis,
                      HighsInfo& highs_info,
                      HighsPrimalDualErrors* primal_dual_errors = nullptr,
                      const bool get_residuals = false);

// Inner getKktFailures
void getKktFailures(const HighsOptions& options, const bool is_qp,
                    const HighsLp& lp, const std::vector<double>& gradient,
                    const HighsSolution& solution, HighsInfo& highs_info,
                    const bool get_residuals = false);

void getVariableKktFailures(const double primal_feasibility_tolerance,
                            const double dual_feasibility_tolerance,
                            const double mip_feasibility_tolerance,
                            const double lower, const double upper,
                            const double value, const double dual,
                            const HighsVarType integrality,
                            double& primal_infeasibility,
                            double& dual_infeasibility,
                            double& semi_infeasibility, uint8_t& at_status,
                            uint8_t& mid_status, const HighsInt index = 0);

void lpNoBasisKktCheck(HighsModelStatus& model_status, HighsInfo& info,
                       const HighsLp& lp, const HighsSolution& solution,
                       const HighsOptions& options, const std::string& message);

void lpKktCheck(HighsModelStatus& model_status, HighsInfo& info,
                const HighsLp& lp, const HighsSolution& solution,
                const HighsBasis& basis, const HighsOptions& options,
                const std::string& message);

std::vector<double> getEffectiveCosts(const HighsLp& lp,
                                      const HighsOptions& options);

void getPrimalDualGlpsolErrors(const HighsOptions& options, const HighsLp& lp,
                               const std::vector<double>& gradient,
                               const HighsSolution& solution,
                               HighsPrimalDualErrors& primal_dual_errors);

void getPrimalDualBasisErrors(const HighsOptions& options, const HighsLp& lp,
                              const HighsSolution& solution,
                              const HighsBasis& basis,
                              HighsPrimalDualErrors& primal_dual_errors);

bool getComplementarityViolations(const HighsLp& lp,
                                  const HighsSolution& solution,
                                  const double optimality_tolerance,
                                  HighsInt& num_complementarity_violations,
                                  double& max_complementarity_violation);

bool computeDualObjectiveValue(const HighsModel& model,
                               const HighsSolution& solution,
                               double& dual_objective_value);

bool computeDualObjectiveValue(const double* gradient, const HighsLp& lp,
                               const HighsSolution& solution,
                               double& dual_objective_value);

double computeObjectiveValue(const HighsLp& lp, const HighsSolution& solution);

void refineBasis(const HighsLp& lp, const HighsSolution& solution,
                 HighsBasis& basis);

HighsStatus ipxSolutionToHighsSolution(
    const HighsOptions& options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const HighsInt ipx_num_col, const HighsInt ipx_num_row,
    const std::vector<double>& ipx_x, const std::vector<double>& ipx_slack_vars,
    const std::vector<double>& ipx_y, const std::vector<double>& ipx_zl,
    const std::vector<double>& ipx_zu, HighsSolution& highs_solution);

HighsStatus ipxBasicSolutionToHighsBasicSolution(
    const HighsLogOptions& log_options, const HighsLp& lp,
    const std::vector<double>& rhs, const std::vector<char>& constraint_type,
    const IpxSolution& ipx_solution, HighsBasis& highs_basis,
    HighsSolution& highs_solution);

HighsStatus formSimplexLpBasisAndFactorReturn(
    const HighsStatus return_status, HighsLpSolverObject& solver_object);
HighsStatus formSimplexLpBasisAndFactor(
    HighsLpSolverObject& solver_object,
    const bool only_from_known_basis = false);

void accommodateAlienBasis(HighsLpSolverObject& solver_object);

void resetModelStatusAndHighsInfo(HighsLpSolverObject& solver_object);
void resetModelStatusAndHighsInfo(HighsModelStatus& model_status,
                                  HighsInfo& highs_info);
bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis);

bool isColPrimalSolutionRightSize(const HighsLp& lp,
                                  const HighsSolution& solution);
bool isRowPrimalSolutionRightSize(const HighsLp& lp,
                                  const HighsSolution& solution);
bool isPrimalSolutionRightSize(const HighsLp& lp,
                               const HighsSolution& solution);

bool isColDualSolutionRightSize(const HighsLp& lp,
                                const HighsSolution& solution);
bool isRowDualSolutionRightSize(const HighsLp& lp,
                                const HighsSolution& solution);
bool isDualSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);

bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis);

bool reportKktFailures(const HighsLp& lp, const HighsOptions& options,
                       const HighsInfo& highs_info,
                       const std::string& message = "");

HighsStatus assessUserColValues(const HighsLp& lp, const HighsOptions& options,
                                const HighsInt num_entries,
                                const HighsInt* index, const double* value,
                                const bool perform_feasibility_check = true);
void assignUserColValues(const HighsInt num_col, std::vector<double>& col_value,
                         const HighsInt num_entries, const HighsInt* index,
                         const double* value);
#endif  // LP_DATA_HIGHSSOLUTION_H_
