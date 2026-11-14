/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/IpxWrapper.h
 * @brief
 */
#ifndef IPM_IPX_WRAPPER_H_
#define IPM_IPX_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "ipm/IpxSolution.h"
#include "ipm/hipo/ipm/Solver.h"
#include "ipm/ipx/ipx_status.h"
#include "ipm/ipx/lp_solver.h"
#include "lp_data/HighsQpSolverObject.h"
#include "lp_data/HighsSolution.h"

HighsStatus solveLpIpx(HighsLpSolverObject& solver_object);

HighsStatus solveLpIpx(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, HighsBasis& highs_basis,
                       HighsSolution& highs_solution,
                       HighsModelStatus& model_status, HighsInfo& highs_info,
                       HighsCallback& callback);

HighsStatus solveLpHipo(HighsLpSolverObject& solver_object);

HighsStatus solveQpHipo(HighsQpSolverObject& solver_object);

HighsStatus solveHipo(const HighsOptions& options, HighsTimer& timer,
                      const HighsLp& lp, const HighsHessian& H,
                      HighsBasis& highs_basis, HighsSolution& highs_solution,
                      HighsModelStatus& model_status, HighsInfo& highs_info,
                      HighsCallback& callback);

HighsStatus reportHipoStatus(const HighsOptions& options,
                             const hipo::Int status, const hipo::Solver& hipo);

HighsStatus reportHipoCrossoverStatus(const HighsOptions& options,
                                      const ipx::Int status, bool is_qp);

void reportHipoNoProgress(const HighsOptions& options,
                          const hipo::Info& hipo_info);

void getHipoNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
                              const hipo::Int num_col, const hipo::Int num_row,
                              const std::vector<double>& rhs,
                              const std::vector<char>& constraint_type,
                              const hipo::Solver& hipo,
                              const HighsModelStatus model_status,
                              HighsSolution& highs_solution);

void fillInIpxData(const HighsLp& lp, ipx::Int& num_col, ipx::Int& num_row,
                   double& offset, std::vector<double>& obj,
                   std::vector<double>& col_lb, std::vector<double>& col_ub,
                   std::vector<ipx::Int>& Ap, std::vector<ipx::Int>& Ai,
                   std::vector<double>& Ax, std::vector<double>& rhs,
                   std::vector<char>& constraint_type);

void fillInRhsAndConstraints(const HighsLp& lp, std::vector<double>& rhs,
                             std::vector<char>& constraint_type);

HighsStatus reportIpxSolveStatus(const HighsOptions& options,
                                 const ipx::Int solve_status,
                                 const ipx::Int error_flag);

HighsStatus reportIpxIpmCrossoverStatus(const HighsOptions& options,
                                        const ipx::Int status,
                                        const bool ipm_status);

bool ipxStatusError(const bool status_error, const HighsOptions& options,
                    const std::string& solver, const std::string& message,
                    const int value = -1);

bool illegalIpxSolvedStatus(const ipx::Info& ipx_info,
                            const HighsOptions& options);

bool illegalIpxStoppedIpmStatus(const ipx::Info& ipx_info,
                                const HighsOptions& options);

bool illegalIpxStoppedCrossoverStatus(const ipx::Info& ipx_info,
                                      const HighsOptions& options);

void reportIpmNoProgress(const HighsOptions& options,
                         const ipx::Info& ipx_info);

void getHighsNonVertexSolution(const HighsOptions& options, const HighsLp& lp,
                               const ipx::Int num_col, const ipx::Int num_row,
                               const std::vector<double>& rhs,
                               const std::vector<char>& constraint_type,
                               const ipx::LpSolver& lps,
                               const HighsModelStatus model_status,
                               HighsSolution& highs_solution);

void reportSolveData(const HighsLogOptions& log_options,
                     const ipx::Info& ipx_info);
#endif
