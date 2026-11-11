/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file qpsolver/QpAsmWrapper.h
 * @brief
 */
#ifndef QPSOLVER_QPASM_WRAPPER_H_
#define QPSOLVER_QPASM_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "lp_data/HighsQpSolverObject.h"
#include "lp_data/HighsSolution.h"

HighsStatus solveQpAsm(HighsQpSolverObject& solver_object);

HighsStatus solveQpAsm(const HighsOptions& options, HighsTimer& timer,
                       const HighsLp& lp, const HighsHessian& hessian,
                       HighsBasis& basis, HighsSolution& solution,
                       HighsModelStatus& model_status, HighsInfo& info,
                       HighsCallback& callback);

#endif
