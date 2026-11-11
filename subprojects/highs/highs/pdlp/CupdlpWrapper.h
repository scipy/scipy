/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/CupdlpWrapper.h
 * @brief
 */
#ifndef PDLP_CUPDLP_WRAPPER_H_
#define PDLP_CUPDLP_WRAPPER_H_

#include <algorithm>
#include <cassert>

#include "lp_data/HighsSolution.h"
#include "pdlp/cupdlp/cupdlp.h"

#ifdef CUPDLP_GPU
#include <cuda_runtime.h>
#endif

// #define CUPDLP_CPP_INIT(var, type, size)                                \
//   {                                                                     \
//     (var) = (type *)malloc((size) * sizeof(type));                      \
//     if ((var) == NULL) {                                                \
//       retcode = 1;                                                      \
//       goto exit_cleanup;                                                \
//     }                                                                   \
//   }

#define cupdlp_init_int(var, size) \
  { (var) = (int*)malloc((size) * sizeof(int)); }

#define cupdlp_init_double(var, size) \
  { (var) = (double*)malloc((size) * sizeof(double)); }

#define cupdlp_init_work(var, size) \
  { (var) = (CUPDLPwork*)malloc((size) * sizeof(CUPDLPwork)); }

#define cupdlp_init_problem(var, size) \
  { (var) = (CUPDLPproblem*)malloc((size) * sizeof(CUPDLPproblem)); }

#define cupdlp_init_data(var, size) \
  { (var) = (CUPDLPdata*)malloc((size) * sizeof(CUPDLPdata)); }

#ifdef CUPDLP_CPU

#define cupdlp_init_vec_double(var, size) \
  { (var) = (double*)malloc((size) * sizeof(double)); }

#define cupdlp_copy_vec(dst, src, type, size) \
  memcpy(dst, src, sizeof(type) * (size))

#endif

//#define cupdlp_init_csc_cpu(var, size)	\
//   {\
//     (var) = (CUPDLPcsc*)malloc((size) * sizeof(CUPDLPcsc));\
//   }

cupdlp_retcode problem_create(CUPDLPproblem** prob);
// cupdlp_retcode csc_create(CUPDLPcsc **csc_cpu);

cupdlp_retcode problem_alloc(
    CUPDLPproblem* prob, cupdlp_int nRows, cupdlp_int nCols, cupdlp_int nEqs,
    cupdlp_float* cost, cupdlp_float offset, cupdlp_float sign_origin,
    void* matrix, CUPDLP_MATRIX_FORMAT src_matrix_format,
    CUPDLP_MATRIX_FORMAT dst_matrix_format, cupdlp_float* rhs,
    cupdlp_float* lower, cupdlp_float* upper, cupdlp_float* alloc_matrix_time,
    cupdlp_float* copy_vec_time);

cupdlp_retcode data_alloc(CUPDLPdata* data, cupdlp_int nRows, cupdlp_int nCols,
                          void* matrix, CUPDLP_MATRIX_FORMAT src_matrix_format,
                          CUPDLP_MATRIX_FORMAT dst_matrix_format);

double infNorm(double* x, cupdlp_int n);

void cupdlp_hasLower(cupdlp_float* haslb, const cupdlp_float* lb,
                     const cupdlp_float bound, const cupdlp_int len);

void cupdlp_hasUpper(cupdlp_float* hasub, const cupdlp_float* ub,
                     const cupdlp_float bound, const cupdlp_int len);

void cupdlp_haslb(cupdlp_float* haslb, const cupdlp_float* lb,
                  const cupdlp_float bound, const cupdlp_int len);

void cupdlp_hasub(cupdlp_float* hasub, const cupdlp_float* ub,
                  const cupdlp_float bound, const cupdlp_int len);

HighsStatus solveLpCupdlp(HighsLpSolverObject& solver_object);

HighsStatus solveLpCupdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback);
int formulateLP_highs(const cupdlp_int local_log_level, const HighsLp& lp,
                      double** cost, int* nCols, int* nRows, int* nnz,
                      int* nEqs, int** csc_beg, int** csc_idx, double** csc_val,
                      double** rhs, double** lower, double** upper,
                      double* offset, double* sign_origin, int* nCols_origin,
                      int** constraint_new_idx, int** constraint_type);

cupdlp_int getCupdlpLogLevel(const HighsOptions& options);

#endif
