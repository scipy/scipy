//
// Created by chuwen on 23-11-28.
//

#ifndef CUPDLP_CUPDLP_RESTART_H
#define CUPDLP_CUPDLP_RESTART_H

#include "pdlp/cupdlp/cupdlp_defs.h"
#include "pdlp/cupdlp/cupdlp_linalg.h"
#include "pdlp/cupdlp/cupdlp_proj.h"
// #include "cupdlp_scaling.h"
#include "pdlp/cupdlp/cupdlp_step.h"
#include "pdlp/cupdlp/cupdlp_utils.h"
#include "pdlp/cupdlp/glbopts.h"

typedef enum {
  PDHG_NO_RESTART = 0,
  PDHG_RESTART_TO_CURRENT,
  PDHG_RESTART_TO_AVERAGE
} PDHG_restart_choice;

cupdlp_bool PDHG_Check_Restart_Merit_Function(CUPDLPwork *work);

PDHG_restart_choice PDHG_Check_Restart_GPU(CUPDLPwork *work);

cupdlp_float PDHG_Restart_Score_GPU(cupdlp_float weightSquared,
                                    cupdlp_float dPrimalFeas,
                                    cupdlp_float dDualFeas,
                                    cupdlp_float dDualityGap);

#endif  // CUPDLP_CUPDLP_RESTART_H
