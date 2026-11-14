//
// Created by chuwen on 23-11-26.
//

#ifndef CUPDLP_CUPDLP_UTILS_H
#define CUPDLP_CUPDLP_UTILS_H

#include <stdio.h>
#ifdef CUPDLP_TIMER
#include <time.h>

#if defined(_WIN32) || defined(_WIN64)
#include <time.h>
#else
#include <sys/time.h>
#endif

#endif
#include "pdlp/cupdlp/cupdlp_defs.h"

#ifdef __cplusplus
extern "C" {
#endif

void data_clear(CUPDLPdata *data);

void problem_clear(CUPDLPproblem *problem);

cupdlp_int vec_clear(CUPDLPvec *vec);

void settings_clear(CUPDLPsettings *settings);

void iterates_clear(CUPDLPiterates *iterates);

void resobj_clear(CUPDLPresobj *resobj);

void stepsize_clear(CUPDLPstepsize *stepsize);

void timers_clear(int log_level, CUPDLPtimers *timers);

void scaling_clear(CUPDLPscaling *scaling);

cupdlp_int PDHG_Clear(CUPDLPwork *w);

void PDHG_PrintPDHGParam(CUPDLPwork *w);

void PDHG_PrintHugeCUPDHG(void);

void PDHG_PrintUserParamHelper(void);

cupdlp_retcode getUserParam(int argc, char **argv,
                            cupdlp_bool *ifChangeIntParam, cupdlp_int *intParam,
                            cupdlp_bool *ifChangeFloatParam,
                            cupdlp_float *floatParam);

cupdlp_retcode settings_SetUserParam(CUPDLPsettings *settings,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam);

cupdlp_retcode resobj_SetUserParam(CUPDLPresobj *resobj,
                                   cupdlp_bool *ifChangeIntParam,
                                   cupdlp_int *intParam,
                                   cupdlp_bool *ifChangeFloatParam,
                                   cupdlp_float *floatParam);

cupdlp_retcode iterates_SetUserParam(CUPDLPiterates *iterates,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam);

cupdlp_retcode stepsize_SetUserParam(CUPDLPstepsize *stepsize,
                                     cupdlp_bool *ifChangeIntParam,
                                     cupdlp_int *intParam,
                                     cupdlp_bool *ifChangeFloatParam,
                                     cupdlp_float *floatParam);

cupdlp_retcode scaling_SetUserParam(CUPDLPscaling *scaling,
                                    cupdlp_bool *ifChangeIntParam,
                                    cupdlp_int *intParam,
                                    cupdlp_bool *ifChangeFloatParam,
                                    cupdlp_float *floatParam);

cupdlp_retcode timers_SetUserParam(CUPDLPtimers *timers,
                                   cupdlp_bool *ifChangeIntParam,
                                   cupdlp_int *intParam,
                                   cupdlp_bool *ifChangeFloatParam,
                                   cupdlp_float *floatParam);

cupdlp_retcode PDHG_SetUserParam(CUPDLPwork *w, cupdlp_bool *ifChangeIntParam,
                                 cupdlp_int *intParam,
                                 cupdlp_bool *ifChangeFloatParam,
                                 cupdlp_float *floatParam);

cupdlp_retcode vec_Alloc(CUPDLPvec *vec, cupdlp_int n);

cupdlp_retcode settings_Alloc(CUPDLPsettings *settings);

cupdlp_retcode resobj_Alloc(CUPDLPresobj *resobj, CUPDLPproblem *problem,
                            cupdlp_int ncols, cupdlp_int nrows);

cupdlp_retcode iterates_Alloc(CUPDLPiterates *iterates, cupdlp_int ncols,
                              cupdlp_int nrows);

cupdlp_retcode stepsize_Alloc(CUPDLPstepsize *stepsize);

cupdlp_retcode timers_Alloc(CUPDLPtimers *timers);

void PDHG_Init_Data(CUPDLPwork *w);

cupdlp_retcode PDHG_Alloc(CUPDLPwork *w);

cupdlp_retcode PDHG_Create(CUPDLPwork **ww, CUPDLPproblem *lp,
                           CUPDLPscaling *scaling);

void PDHG_Destroy(CUPDLPwork **w);

void vecPrint(const char *s, const cupdlp_float *a, cupdlp_int n);

void vecIntPrint(const char *s, const cupdlp_int *a, cupdlp_int n);

void PDHG_Dump_Stats(CUPDLPwork *w);

/* TODO: Add compatibility for Windows platform */
double my_clock(void);

extern double getTimeStamp(void);

cupdlp_retcode dense_create(CUPDLPdense **dense);

cupdlp_retcode csr_create(CUPDLPcsr **csr);

cupdlp_retcode csc_create(CUPDLPcsc **csc);

cupdlp_retcode dense_alloc(CUPDLPdense *dense, cupdlp_int nRows,
                           cupdlp_int nCols, cupdlp_float *val);

cupdlp_retcode csr_alloc(CUPDLPcsr *csr, cupdlp_int nRows, cupdlp_int nCols,
                         cupdlp_int nnz, cupdlp_int *col_ptr,
                         cupdlp_int *row_ind, cupdlp_float *val);

cupdlp_retcode csc_alloc(CUPDLPcsc *csc, cupdlp_int nRows, cupdlp_int nCols,
                         cupdlp_int nnz, cupdlp_int *row_ptr,
                         cupdlp_int *col_ind, cupdlp_float *val);

cupdlp_retcode dense_alloc_matrix(CUPDLPdense *dense, cupdlp_int nRows,
                                  cupdlp_int nCols, void *src,
                                  CUPDLP_MATRIX_FORMAT src_matrix_format);

cupdlp_retcode csr_alloc_matrix(CUPDLPcsr *csr, cupdlp_int nRows,
                                cupdlp_int nCols, void *src,
                                CUPDLP_MATRIX_FORMAT src_matrix_format);

cupdlp_retcode csc_alloc_matrix(CUPDLPcsc *csc, cupdlp_int nRows,
                                cupdlp_int nCols, void *src,
                                CUPDLP_MATRIX_FORMAT src_matrix_format);

void dense2csr(CUPDLPcsr *csr, CUPDLPdense *dense);

void dense2csc(CUPDLPcsc *csc, CUPDLPdense *dense);

void csr2csc(CUPDLPcsc *csc, CUPDLPcsr *csr);

cupdlp_int csc2csr(CUPDLPcsr *csr, CUPDLPcsc *csc);

void csr2dense(CUPDLPdense *dense, CUPDLPcsr *csr);

void csc2dense(CUPDLPdense *dense, CUPDLPcsc *csc);

cupdlp_int csc_clear(CUPDLPcsc *csc);

cupdlp_int csc_clear_host(CUPDLPcsc *csc);

cupdlp_int csr_clear(CUPDLPcsr *csr);

void dense_clear(CUPDLPdense *dense);

void dense_copy(CUPDLPdense *dst, CUPDLPdense *src);

void csr_copy(CUPDLPcsr *dst, CUPDLPcsr *src);

cupdlp_int csc_copy(CUPDLPcsc *dst, CUPDLPcsc *src);

void csrPrintDense(const char *s, CUPDLPcsr *csr);

void cscPrintDense(const char *s, CUPDLPcsc *csc);

void writeJson(const char *fout, CUPDLPwork *work);

void writeSol(const char *fout, cupdlp_int nCols, cupdlp_int nRows,
              cupdlp_float *col_value, cupdlp_float *col_dual,
              cupdlp_float *row_value, cupdlp_float *row_dual);

  void debugPdlpIterHeaderLog(FILE* file);
  void debugPdlpDataInitialise(struct DebugPdlpData* debug_pdlp);
  void debugPdlpIterLog(FILE* file, const int iter_num, const struct DebugPdlpData* debug_pdlp, const double beta, const double primal_step, const double dual_step);
  void debugPdlpFeasOptLog(FILE* file,
			   const int iter_num,
			   const double primal_obj, const double dual_obj,
			   const double gap, const double primal_feas, const double dual_feas,
			   const char* type);
  void debugPdlpRestartLog(FILE* file, const int iter_num, const double current_score, const double average_score);
  void debugPdlpRestarScoretLog(FILE* file, const double weight_squared, const double primal_feas, const double dual_feas, const double obj_gap);
  void debugPdlpFinalSolutionLog(FILE* file, const double* x, int nCols, const double* y, int nRows);


#ifdef __cplusplus
}
#endif
#endif  // CUPDLP_CUPDLP_UTILS_H
