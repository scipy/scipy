//
// Created by chuwen on 23-11-28.
//

#ifndef CUPDLP_CUPDLP_PROJ_H
#define CUPDLP_CUPDLP_PROJ_H

#include "pdlp/cupdlp/cupdlp_defs.h"
#include "pdlp/cupdlp/glbopts.h"

void PDHG_Project_Bounds(CUPDLPwork *work, double *r);

// void PDHG_Project_Row_Duals(CUPDLPwork *work, double *r);

void PDHG_Restart_Iterate(CUPDLPwork *pdhg);

void PDHG_Restart_Iterate_GPU(CUPDLPwork *pdhg);

#endif  // CUPDLP_CUPDLP_PROJ_H
