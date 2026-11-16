#ifndef CUPDLP_CUDA_LINALG_H
#define CUPDLP_CUDA_LINALG_H

#include <cublas_v2.h>         // cublas
#include <cuda_runtime_api.h>  // cudaMalloc, cudaMemcpy, etc.
#include <cusparse.h>          // cusparseSpMV

#include "cupdlp_cuda_kernels.cuh"

#define PRINT_CUDA_INFO (1)
#define PRINT_DETAILED_CUDA_INFO (0)

#ifdef __cplusplus
extern "C" {
#endif

cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer_csc_ATy, void **dBuffer_csr_Ax);

/*
cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                       cusparseSpMatDescr_t cuda_csc,
                       cusparseDnVecDescr_t vecX,
                       cusparseDnVecDescr_t vecAx, void *dBuffer,
                       cupdlp_float alpha, cupdlp_float beta);
*/
cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                       cusparseSpMatDescr_t cuda_csr,
                       cusparseDnVecDescr_t vecX,
                       cusparseDnVecDescr_t vecAx, void *dBuffer,
                       cupdlp_float alpha, cupdlp_float beta);
cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                        cusparseSpMatDescr_t cuda_csc,
                        cusparseDnVecDescr_t vecY,
                        cusparseDnVecDescr_t vecATy, void *dBuffer,
                        cupdlp_float alpha, cupdlp_float beta);
/*
cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                        cusparseSpMatDescr_t cuda_csr,
                        cusparseDnVecDescr_t vecY,
                        cusparseDnVecDescr_t vecATy, void *dBuffer,
                        cupdlp_float alpha, cupdlp_float beta);
*/

// z[i] = min(ax[i] - rhs[i], i >= nEqs ? 0.0 : ax[i] - rhs[i]) * (ifScaled ? rowScale[i] : 1.0)
void cupdlp_primal_feasibility_kernel_cuda(cupdlp_float *z,
                                           const cupdlp_float *ax,
                                           const cupdlp_float *rhs,
                                           const cupdlp_float *rowScale,
                                           int ifScaled,
                                           int nEqs, int nRows);

// z[i] = cost[i] - aty[i]
void cupdlp_dual_feasibility_kernel_1_cuda(cupdlp_float *z,
                                           const cupdlp_float *aty,
                                           const cupdlp_float *cost,
                                           int nCols);

// z[i] = max(dualResidual, 0) * hasLower[i]
void cupdlp_dual_feasibility_kernel_2_cuda(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasLower,
                                           int nCols);

// z[i] = - min(dualResidual, 0) * hasUpper[i]
void cupdlp_dual_feasibility_kernel_3_cuda(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasUpper,
                                           int nCols);

// z[i] = alpha * (aty[i] + dSlackPos[i] - dSlackNeg[i]) * (ifScaled ? colScale[i] : 1.0)
void cupdlp_primal_infeasibility_kernel_cuda(cupdlp_float *z, const cupdlp_float *aty,
                                             const cupdlp_float *dSlackPos,
                                             const cupdlp_float *dSlackNeg,
                                             const cupdlp_float *colScale,
                                             cupdlp_float alpha, int ifScaled, int nCols);

// z[i] = min(alpha * x[i], 0.0) * hasLower[i] / (ifScaled ? colScale[i] : 1.0)
void cupdlp_dual_infeasibility_kernel_lb_cuda(cupdlp_float *z, const cupdlp_float *x,
                                              const cupdlp_float *hasLower,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols);

// z[i] = max(alpha * x[i], 0.0) * hasUpper[i] / (ifScaled ? colScale[i] : 1.0)
void cupdlp_dual_infeasibility_kernel_ub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                              const cupdlp_float *hasUpper,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols);

// z[i] = min(alpha * ax[i], i >= nEqs ? 0.0 : alpha * ax[i]) * (ifScaled ? rowScale[i] : 1.0)
void cupdlp_dual_infeasibility_kernel_constr_cuda(cupdlp_float *z, const cupdlp_float *ax,
                                                  const cupdlp_float *rowScale,
                                                  cupdlp_float alpha, int ifScaled,
                                                  int nEqs, int nRows);

void cupdlp_projSameub_cuda(cupdlp_float *x, const cupdlp_float ub, int n);
void cupdlp_projSamelb_cuda(cupdlp_float *x, const cupdlp_float lb, int n);

void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub, int n);
void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb, int n);

void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y, int n);

void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y, int n);

void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                       cupdlp_float bound, int n);
void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                       cupdlp_float bound, int n);

void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                          cupdlp_float bound, int n);
void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                          cupdlp_float bound, int n);

void cupdlp_initvec_cuda(cupdlp_float *x, cupdlp_float val, int n);

void cupdlp_pgrad_cuda(cupdlp_float *xUpdate, const cupdlp_float *x,
                       const cupdlp_float *cost, const cupdlp_float *ATy,
                       const cupdlp_float *lb, const cupdlp_float *ub,
                       cupdlp_float dPrimalStep, int nCols);

void cupdlp_dgrad_cuda(cupdlp_float *yUpdate,
                       const cupdlp_float *y, const cupdlp_float *b,
                       const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
                       cupdlp_float dDualStep, int nRows, int nEqs);

/*
void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                const cupdlp_float *y, const cupdlp_int len);
*/

void cupdlp_movement_interaction_cuda(
    cupdlp_float *dX2, cupdlp_float *dY2, cupdlp_float *dInter, cupdlp_float *buffer,
    const cupdlp_float *xUpdate, const cupdlp_float *x,
    const cupdlp_float *yUpdate, const cupdlp_float *y,
    const cupdlp_float *atyUpdate, const cupdlp_float *aty,
    int nRows, int nCols);

cupdlp_int print_cuda_info(cusparseHandle_t handle);

#ifdef __cplusplus
}
#endif

#endif
