#include <stdio.h>   // printf
#include <stdlib.h>  // EXIT_FAILURE

#include "cupdlp_cudalinalg.cuh"

inline int nBlocks256(int n) {
  constexpr int BLOCKS_PER_SM = 32;
  int numSMs;
  CHECK_CUDA_IGNORE(cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0))
  return min((n + 256 - 1) / 256, BLOCKS_PER_SM * numSMs);
}

extern "C" {

cupdlp_int cuda_alloc_MVbuffer(
    cusparseHandle_t handle, cusparseSpMatDescr_t cuda_csc,
    cusparseDnVecDescr_t vecX, cusparseDnVecDescr_t vecAx,
    cusparseSpMatDescr_t cuda_csr, cusparseDnVecDescr_t vecY,
    cusparseDnVecDescr_t vecATy, void **dBuffer_csc_ATy, void **dBuffer_csr_Ax) {

  size_t AxBufferSize = 0;
  size_t ATyBufferSize = 0;
  cupdlp_float alpha = 1.0;
  cupdlp_float beta = 0.0;
  // cusparseSpSVAlg_t alg = CUSPARSE_SPSV_ALG_DEFAULT;
  cusparseSpMVAlg_t alg = CUSPARSE_SPMV_CSR_ALG2; //deterministic

  // get the buffer size needed by csr Ax
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, cuda_csr, vecX, &beta,
      vecAx, CudaComputeType, alg, &AxBufferSize))

  // allocate an external buffer if needed
  CHECK_CUDA(cudaMalloc(dBuffer_csr_Ax, AxBufferSize))

  // get the buffer size needed by csc ATy
  CHECK_CUSPARSE(cusparseSpMV_bufferSize(
      handle, CUSPARSE_OPERATION_TRANSPOSE, &alpha, cuda_csc, vecY, &beta,
      vecATy, CudaComputeType, alg, &ATyBufferSize))

  // allocate an external buffer if needed
  CHECK_CUDA(cudaMalloc(dBuffer_csc_ATy, ATyBufferSize))

  return EXIT_SUCCESS;
}

/*
cupdlp_int cuda_csc_Ax(cusparseHandle_t handle,
                       cusparseSpMatDescr_t cuda_csc,
                       cusparseDnVecDescr_t vecX,
                       cusparseDnVecDescr_t vecAx, void *dBuffer,
                       cupdlp_float alpha, cupdlp_float beta) {
  // Ax = alpha * Acsc * X + beta * Ax

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecX, &beta, vecAx,
                              // CudaComputeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              CudaComputeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}
*/

cupdlp_int cuda_csr_Ax(cusparseHandle_t handle,
                       cusparseSpMatDescr_t cuda_csr,
                       cusparseDnVecDescr_t vecX,
                       cusparseDnVecDescr_t vecAx, void *dBuffer,
                       cupdlp_float alpha, cupdlp_float beta) {
  // Ax = alpha * Acsr * X + beta * Ax

  cusparseOperation_t op = CUSPARSE_OPERATION_NON_TRANSPOSE;

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecX, &beta, vecAx,
                              // CudaComputeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              CudaComputeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

cupdlp_int cuda_csc_ATy(cusparseHandle_t handle,
                        cusparseSpMatDescr_t cuda_csc,
                        cusparseDnVecDescr_t vecY,
                        cusparseDnVecDescr_t vecATy, void *dBuffer,
                        cupdlp_float alpha, cupdlp_float beta) {
  // ATy = alpha * Acsc^T * Y + beta * ATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csc, vecY, &beta, vecATy,
                              // CudaComputeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              CudaComputeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}

/*
cupdlp_int cuda_csr_ATy(cusparseHandle_t handle,
                        cusparseSpMatDescr_t cuda_csr,
                        cusparseDnVecDescr_t vecY,
                        cusparseDnVecDescr_t vecATy, void *dBuffer,
                        cupdlp_float alpha, cupdlp_float beta) {
  // ATy = alpha * Acsr^T * Y + beta * ATy
  cusparseOperation_t op = CUSPARSE_OPERATION_TRANSPOSE;

  CHECK_CUSPARSE(cusparseSpMV(handle, op, &alpha, cuda_csr, vecY, &beta, vecATy,
                              // CudaComputeType, CUSPARSE_SPMV_ALG_DEFAULT, dBuffer))
                              CudaComputeType, CUSPARSE_SPMV_CSR_ALG2, dBuffer))

  return EXIT_SUCCESS;
}
*/

void cupdlp_primal_feasibility_kernel_cuda(cupdlp_float *z,
                                           const cupdlp_float *ax,
                                           const cupdlp_float *rhs,
                                           const cupdlp_float *rowScale,
                                           int ifScaled,
                                           int nEqs, int nRows) {
  element_primal_feas_kernel<<<nBlocks256(nRows), 256>>>(z, ax, rhs, rowScale,
                                                         ifScaled, nEqs, nRows);
}

void cupdlp_dual_feasibility_kernel_1_cuda(cupdlp_float *z,
                                           const cupdlp_float *aty,
                                           const cupdlp_float *cost,
                                           int nCols) {
  element_dual_feas_kernel_1<<<nBlocks256(nCols), 256>>>(z, aty, cost, nCols);
}

void cupdlp_dual_feasibility_kernel_2_cuda(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasLower,
                                           int nCols) {
  element_dual_feas_kernel_2<<<nBlocks256(nCols), 256>>>(z, dualResidual,
                                                         hasLower, nCols);
}

void cupdlp_dual_feasibility_kernel_3_cuda(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasUpper,
                                           int nCols) {
  element_dual_feas_kernel_3<<<nBlocks256(nCols), 256>>>(z, dualResidual,
                                                         hasUpper, nCols);
}

void cupdlp_primal_infeasibility_kernel_cuda(cupdlp_float *z, const cupdlp_float *aty,
                                             const cupdlp_float *dSlackPos,
                                             const cupdlp_float *dSlackNeg,
                                             const cupdlp_float *colScale,
                                             cupdlp_float alpha, int ifScaled, int nCols) {
  element_primal_infeas_kernel<<<nBlocks256(nCols), 256>>>(z, aty, dSlackPos,
                                                           dSlackNeg, colScale, alpha,
                                                           ifScaled, nCols);
}

void cupdlp_dual_infeasibility_kernel_lb_cuda(cupdlp_float *z, const cupdlp_float *x,
                                              const cupdlp_float *hasLower,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols) {
  element_dual_infeas_kernel_lb<<<nBlocks256(nCols), 256>>>(z, x, hasLower,
                                                            colScale, alpha,
                                                            ifScaled, nCols);
}

void cupdlp_dual_infeasibility_kernel_ub_cuda(cupdlp_float *z, const cupdlp_float *x,
                                              const cupdlp_float *hasUpper,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols) {
  element_dual_infeas_kernel_ub<<<nBlocks256(nCols), 256>>>(z, x, hasUpper,
                                                            colScale, alpha,
                                                            ifScaled, nCols);
}

void cupdlp_dual_infeasibility_kernel_constr_cuda(cupdlp_float *z, const cupdlp_float *ax,
                                                  const cupdlp_float *rowScale,
                                                  cupdlp_float alpha, int ifScaled,
                                                  int nEqs, int nRows) {
  element_dual_infeas_kernel_constr<<<nBlocks256(nRows), 256>>>(z, ax, rowScale,
                                                                alpha, ifScaled, nEqs, nRows);
}

void cupdlp_projSameub_cuda(cupdlp_float *x, cupdlp_float ub, int n) {
  element_wise_projSameub_kernel<<<nBlocks256(n), 256>>>(x, ub, n);
}

void cupdlp_projSamelb_cuda(cupdlp_float *x, cupdlp_float lb, int n) {
  element_wise_projSamelb_kernel<<<nBlocks256(n), 256>>>(x, lb, n);
}

void cupdlp_projub_cuda(cupdlp_float *x, const cupdlp_float *ub, int n) {
  element_wise_projub_kernel<<<nBlocks256(n), 256>>>(x, ub, n);
}

void cupdlp_projlb_cuda(cupdlp_float *x, const cupdlp_float *lb, int n) {
  element_wise_projlb_kernel<<<nBlocks256(n), 256>>>(x, lb, n);
}

void cupdlp_ediv_cuda(cupdlp_float *x, const cupdlp_float *y, int n) {
  element_wise_div_kernel<<<nBlocks256(n), 256>>>(x, y, n);
}

void cupdlp_edot_cuda(cupdlp_float *x, const cupdlp_float *y, int n) {
  element_wise_dot_kernel<<<nBlocks256(n), 256>>>(x, y, n);
}

void cupdlp_haslb_cuda(cupdlp_float *haslb, const cupdlp_float *lb,
                       cupdlp_float bound, int n) {
  element_wise_initHaslb_kernel<<<nBlocks256(n), 256>>>(haslb, lb, bound, n);
}

void cupdlp_hasub_cuda(cupdlp_float *hasub, const cupdlp_float *ub,
                       cupdlp_float bound, int n) {
  element_wise_initHasub_kernel<<<nBlocks256(n), 256>>>(hasub, ub, bound, n);
}

void cupdlp_filterlb_cuda(cupdlp_float *x, const cupdlp_float *lb,
                          cupdlp_float bound, int n) {
  element_wise_filterlb_kernel<<<nBlocks256(n), 256>>>(x, lb, bound, n);
}

void cupdlp_filterub_cuda(cupdlp_float *x, const cupdlp_float *ub,
                          cupdlp_float bound, int n) {
  element_wise_filterub_kernel<<<nBlocks256(n), 256>>>(x, ub, bound, n);
}

void cupdlp_initvec_cuda(cupdlp_float *x, cupdlp_float val, int n) {
  init_cuda_vec_kernel<<<nBlocks256(n), 256>>>(x, val, n);
}

void cupdlp_pgrad_cuda(cupdlp_float *xUpdate, const cupdlp_float *x,
                       const cupdlp_float *cost, const cupdlp_float *ATy,
                       const cupdlp_float *lb, const cupdlp_float *ub,
                       cupdlp_float dPrimalStep, int nCols) {
  primal_grad_step_kernel<<<nBlocks256(nCols), 256>>>(xUpdate, x, cost, ATy, lb, ub, dPrimalStep, nCols);
}

void cupdlp_dgrad_cuda(cupdlp_float *yUpdate,
                       const cupdlp_float *y, const cupdlp_float *b,
                       const cupdlp_float *Ax, const cupdlp_float *AxUpdate,
                       cupdlp_float dDualStep, int nRows, int nEqs) {
  dual_grad_step_kernel<<<nBlocks256(nRows), 256>>>(yUpdate, y, b, Ax, AxUpdate, dDualStep, nRows, nEqs);
}

/*
void cupdlp_sub_cuda(cupdlp_float *z, const cupdlp_float *x,
                     const cupdlp_float *y, int n)
{
  naive_sub_kernel<<<nBlocks256(n), 256>>>(z, x, y, n);
}
*/


void cupdlp_movement_interaction_cuda(
    cupdlp_float *dX2, cupdlp_float *dY2, cupdlp_float *dInter, cupdlp_float *buffer,
    const cupdlp_float *xUpdate, const cupdlp_float *x,
    const cupdlp_float *yUpdate, const cupdlp_float *y,
    const cupdlp_float *atyUpdate, const cupdlp_float *aty,
    int nRows, int nCols)
{
  int warpSize;
  CHECK_CUDA_IGNORE(cudaDeviceGetAttribute(&warpSize, cudaDevAttrWarpSize, 0))
  if (warpSize != 32) {
    printf("warpSize\n");
    exit(1);
  }

  constexpr int RED_BLOCK_SIZE = 256;
  constexpr int RED_ELS_PER_THREAD = 4;
  constexpr int RED_ELS_PER_BLOCK = RED_BLOCK_SIZE * RED_ELS_PER_THREAD;

  constexpr int SUM_BLOCK_SIZE = 512;
  constexpr int SUM_ELS_PER_THREAD = 2;
  constexpr int SUM_ELS_PER_BLOCK = SUM_BLOCK_SIZE * SUM_ELS_PER_THREAD;

  int nBlocksCols = (nCols + RED_ELS_PER_BLOCK - 1) / RED_ELS_PER_BLOCK;
  int nBlocksRows = (nRows + RED_ELS_PER_BLOCK - 1) / RED_ELS_PER_BLOCK;

  int buf_size = 256 * ((max(nBlocksCols,nBlocksRows) + 256 - 1) / 256);

  // buffer has size max(nCols, nRows, 2048)

  // we have buf_size <= max(nBlocksCols, nBlocksRows) + 256 - 1
  //                  <= max(nCols, nRows) / (4 * 256) + 1 + 256 - 1,
  // since RED_ELS_PER_BLOCK = 4 * 256

  // assume wlog that nRows <= nCols

  // if nCols <= 2 * 4 * 256, then
  // buf_size <= 2 + 256 <= 258,
  // so 5 * buf_size <= 1290 <= max(nCols, nRows, 2048)

  // if nCols > 2 * 4 * 256, then
  // buf_size <= nCols / (4 * 256) + nCols / 8 <= nCols * 3/16
  // so 5 * buf_size <= 15/16 * nCols <= max(nCols, nRows, 2048)

  cupdlp_float *buf_1 = buffer + 0 * buf_size;
  cupdlp_float *buf_2 = buffer + 1 * buf_size;
  cupdlp_float *buf_3 = buffer + 2 * buf_size;
  cupdlp_float *buf_4 = buffer + 3 * buf_size;
  cupdlp_float *buf_5 = buffer + 4 * buf_size;

  int nBlocks = nBlocksCols;
  movement_1_kernel<<<nBlocks, RED_BLOCK_SIZE>>>(buf_1, buf_2, xUpdate, x, atyUpdate, aty, nCols);

  while (nBlocks > 1) {
    int nBlocks2 = (nBlocks + SUM_ELS_PER_BLOCK - 1) / SUM_ELS_PER_BLOCK;
    sum_kernel<<<nBlocks2, SUM_BLOCK_SIZE>>>(buf_3, buf_1, nBlocks);
    sum_kernel<<<nBlocks2, SUM_BLOCK_SIZE>>>(buf_4, buf_2, nBlocks);
    nBlocks = nBlocks2;
    cupdlp_float *tmp = buf_1;
    buf_1 = buf_3;
    buf_3 = tmp;
    tmp = buf_2;
    buf_2 = buf_4;
    buf_4 = tmp;
  }

  CHECK_CUDA_STRICT(cudaMemcpyAsync(buf_5 + 0, buf_1, sizeof(cupdlp_float), cudaMemcpyDeviceToDevice))
  CHECK_CUDA_STRICT(cudaMemcpyAsync(buf_5 + 1, buf_2, sizeof(cupdlp_float), cudaMemcpyDeviceToDevice))

  nBlocks = nBlocksRows;
  movement_2_kernel<<<nBlocks, RED_BLOCK_SIZE>>>(buf_1, yUpdate, y, nRows);

  while (nBlocks > 1) {
    int nBlocks2 = (nBlocks + SUM_ELS_PER_BLOCK - 1) / SUM_ELS_PER_BLOCK;
    sum_kernel<<<nBlocks2, SUM_BLOCK_SIZE>>>(buf_2, buf_1, nBlocks);
    nBlocks = nBlocks2;
    cupdlp_float *tmp = buf_1;
    buf_1 = buf_2;
    buf_2 = tmp;
  }

  cupdlp_float res[3];
  CHECK_CUDA_STRICT(cudaMemcpyAsync(buf_5 + 2, buf_1, sizeof(cupdlp_float), cudaMemcpyDeviceToDevice))
  CHECK_CUDA_STRICT(cudaDeviceSynchronize())
  CHECK_CUDA_STRICT(cudaMemcpy(res, buf_5, 3 * sizeof(cupdlp_float), cudaMemcpyDeviceToHost))
  CHECK_CUDA_LAST();

  *dX2 = res[0];
  *dY2 = res[2];
  *dInter = res[1];
}

cupdlp_int print_cuda_info(cusparseHandle_t handle)
{
#if PRINT_CUDA_INFO

  int v_cuda_runtime = 0;
  int v_cuda_driver = 0;
  int v_cusparse = 0;
  CHECK_CUDA(cudaRuntimeGetVersion(&v_cuda_runtime))
  CHECK_CUDA(cudaDriverGetVersion(&v_cuda_driver))
  CHECK_CUSPARSE(cusparseGetVersion(handle, &v_cusparse))

  printf("Cuda runtime version %d\n", v_cuda_runtime);
  printf("Cuda driver  version %d\n", v_cuda_driver);
  printf("cuSparse     version %d\n", v_cusparse);

  int n_devices = 0;
  CHECK_CUDA(cudaGetDeviceCount(&n_devices))

  for (int i = 0; i < n_devices; i++) {
    cudaDeviceProp prop;
    CHECK_CUDA(cudaGetDeviceProperties(&prop, i));

    printf("Cuda device %d: %s\n", i, prop.name);
#if PRINT_DETAILED_CUDA_INFO
    printf("  Clock rate (KHz): %d\n", prop.clockRate);
    printf("  Memory clock rate (KHz): %d\n", prop.memoryClockRate);
    printf("  Memory bus width (bits): %d\n", prop.memoryBusWidth);
    printf("  Peak memory bandwidth (GB/s): %f\n",
            2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
    printf("  Global memory available on device (GB): %f\n", prop.totalGlobalMem / 1.0e9);
    printf("  Shared memory available per block (B): %zu\n", prop.sharedMemPerBlock);
    printf("  Warp size in threads: %d\n", prop.warpSize);
    printf("  Maximum number of threads per block: %d\n", prop.maxThreadsPerBlock);
    printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
    printf("  Number of multiprocessors on device: %d\n", prop.multiProcessorCount);
#endif
  }
#endif

  return EXIT_SUCCESS;
}

}
