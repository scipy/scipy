//
// Created by chuwen on 23-11-27.
//

#ifndef CUPDLP_CUPDLP_SOLVER_H
#define CUPDLP_CUPDLP_SOLVER_H

#include "pdlp/cupdlp/cupdlp_defs.h"
#include "pdlp/cupdlp/glbopts.h"

#ifdef __cplusplus
extern "C" {
#endif
#define CUPDLP_CHECK_TIMEOUT(pdhg)                               \
  {                                                              \
    PDHG_Compute_SolvingTime(pdhg);                              \
    if (pdhg->timers->dSolvingTime > pdhg->settings->dTimeLim) { \
      retcode = RETCODE_FAILED;                                  \
      goto exit_cleanup;                                         \
    }                                                            \
  }

void PDHG_Compute_Primal_Feasibility(CUPDLPwork *work, cupdlp_float *primalResidual,
                                     const cupdlp_float *ax, const cupdlp_float *x,
                                     cupdlp_float *dPrimalFeasibility,
                                     cupdlp_float *dPrimalObj);

void PDHG_Compute_Dual_Feasibility(CUPDLPwork *work, cupdlp_float *dualResidual,
                                   const cupdlp_float *aty, const cupdlp_float *x,
                                   const cupdlp_float *y, cupdlp_float *dDualFeasibility,
                                   cupdlp_float *dDualObj, cupdlp_float *dComplementarity,
                                   cupdlp_float *dSlackPos, cupdlp_float *dSlackNeg);

void PDHG_Compute_Residuals(CUPDLPwork *work);

void PDHG_Compute_Primal_Infeasibility(CUPDLPwork *work, const cupdlp_float *y,
                                       const cupdlp_float *dSlackPos,
                                       const cupdlp_float *dSlackNeg,
                                       const cupdlp_float *aty,
                                       const cupdlp_float dualObj,
                                       cupdlp_float *dPrimalInfeasObj,
                                       cupdlp_float *dPrimalInfeasRes);

void PDHG_Compute_Dual_Infeasibility(CUPDLPwork *work, const cupdlp_float *x,
                                     const cupdlp_float *ax,
                                     const cupdlp_float primalObj,
                                     cupdlp_float *dDualInfeasObj,
                                     cupdlp_float *dDualInfeasRes);

void PDHG_Compute_Infeas_Residuals(CUPDLPwork *work);

void PDHG_Init_Variables(const cupdlp_int* has_variables, CUPDLPwork *work);

void PDHG_Check_Data(CUPDLPwork *work);

cupdlp_bool PDHG_Check_Termination(CUPDLPwork *pdhg, int bool_print);

cupdlp_bool PDHG_Check_Termination_Average(CUPDLPwork *pdhg, int bool_print);

termination_code PDHG_Check_Infeasibility(CUPDLPwork *pdhg, int bool_print);

termination_code PDHG_Check_Primal_Infeasibility(CUPDLPwork *pdhg,
                                                 cupdlp_float dPrimalInfeasObj,
                                                 cupdlp_float dPrimalInfeasRes);
termination_code PDHG_Check_Dual_Infeasibility(CUPDLPwork *pdhg,
                                               cupdlp_float dDualInfeasObj,
                                               cupdlp_float dDualInfeasRes);

void PDHG_Print_Header(CUPDLPwork *pdhg);

void PDHG_Print_Iter(CUPDLPwork *pdhg);

void PDHG_Print_Iter_Average(CUPDLPwork *pdhg);

void PDHG_Compute_SolvingTime(CUPDLPwork *pdhg);

cupdlp_retcode PDHG_Solve(const cupdlp_int* has_variables, CUPDLPwork *pdhg);

cupdlp_retcode PDHG_PostSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
                              cupdlp_int *constraint_new_idx,
                              cupdlp_int *constraint_type,
                              cupdlp_float *col_value, cupdlp_float *col_dual,
                              cupdlp_float *row_value, cupdlp_float *row_dual,
                              cupdlp_int *value_valid, cupdlp_int *dual_valid);

cupdlp_retcode PDHG_PreSolve(CUPDLPwork *pdhg, cupdlp_int nCols_origin,
			     cupdlp_int *constraint_new_idx,
			     cupdlp_int *constraint_type,
			     cupdlp_float *col_value, cupdlp_float *col_dual,
			     cupdlp_float *row_value, cupdlp_float *row_dual,
			     cupdlp_int *value_valid, cupdlp_int *dual_valid);

cupdlp_retcode LP_SolvePDHG(
    CUPDLPwork *pdhg, cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
    cupdlp_bool *ifChangeFloatParam, cupdlp_float *floatParam, char *fp,
    cupdlp_int nCols_origin, cupdlp_float *col_value, cupdlp_float *col_dual,
    cupdlp_float *row_value, cupdlp_float *row_dual, cupdlp_int *value_valid,
    cupdlp_int *dual_valid, cupdlp_bool ifSaveSol, char *fp_sol,
    cupdlp_int *constraint_new_idx, cupdlp_int *constraint_type,
    cupdlp_int *model_status, cupdlp_int* num_iter);

#ifdef __cplusplus
}
#endif
#endif  // CUPDLP_CUPDLP_SOLVER_H
