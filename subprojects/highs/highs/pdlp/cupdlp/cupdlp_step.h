//
// Created by chuwen on 23-11-28.
//

#ifndef CUPDLP_CUPDLP_STEP_H
#define CUPDLP_CUPDLP_STEP_H

#include "pdlp/cupdlp/cupdlp_defs.h"
// #include "cupdlp_scaling.h"
#include "pdlp/cupdlp/glbopts.h"

cupdlp_retcode PDHG_Power_Method(CUPDLPwork *work, double *lambda);

void PDHG_Compute_Step_Size_Ratio(CUPDLPwork *pdhg);

void PDHG_Update_Iterate_Constant_Step_Size(CUPDLPwork *pdhg);

void PDHG_Update_Iterate_Malitsky_Pock(CUPDLPwork *pdhg);

cupdlp_retcode PDHG_Update_Iterate_Adaptive_Step_Size(CUPDLPwork *pdhg);

cupdlp_retcode PDHG_Init_Step_Sizes(CUPDLPwork *pdhg);

void PDHG_Compute_Average_Iterate(CUPDLPwork *work);

void PDHG_Update_Average(CUPDLPwork *work);

cupdlp_retcode PDHG_Update_Iterate(CUPDLPwork *pdhg);

void PDHG_primalGradientStep(CUPDLPwork *work, CUPDLPvec *xUpdate,
                             const CUPDLPvec *x, const CUPDLPvec *ATy,
                             cupdlp_float dPrimalStepSize);
void PDHG_dualGradientStep(CUPDLPwork *work, CUPDLPvec *yUpdate,
                           const CUPDLPvec *y, const CUPDLPvec *Ax,
                           const CUPDLPvec *AxUpdate, cupdlp_float dDualStepSize);

#endif  // CUPDLP_CUPDLP_STEP_H
