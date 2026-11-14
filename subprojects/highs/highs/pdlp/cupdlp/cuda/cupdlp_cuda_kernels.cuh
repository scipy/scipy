#ifndef CUPDLP_CUDA_KERNALS_H
#define CUPDLP_CUDA_KERNALS_H

#include <stdio.h>
#include <stdlib.h>  /* EXIT_FAILURE */
#include <cublas_v2.h>
#include <cusparse.h>
#include <cuda_runtime.h>

//#define CUPDLP_BLOCK_SIZE 512

#ifndef SFLOAT

#ifdef DLONG
typedef long long cupdlp_int;
#else
typedef int cupdlp_int;
#endif
typedef double cupdlp_float;
#define CudaComputeType CUDA_R_64F
#define cupdlp_fma_rn __fma_rn

#else

typedef float cupdlp_float;
#define CudaComputeType CUDA_R_32F
#define cupdlp_fma_rn __fmaf_rn

#endif

static inline cudaError_t check_cuda_call(cudaError_t status,
                                          const char *filename, int line)
{
  if (status != cudaSuccess) {
    printf("CUDA API failed at line %d of %s with error: %s (%d)\n",
      line, filename, cudaGetErrorString(status), status);
  }
  return status;
}

static inline cusparseStatus_t check_cusparse_call(cusparseStatus_t status,
                                                   const char *filename, int line)
{
  if (status != CUSPARSE_STATUS_SUCCESS) {
    printf("CUSPARSE API failed at line %d of %s with error: %s (%d)\n",
      line, filename, cusparseGetErrorString(status), status);
  }
  return status;
}

static inline cublasStatus_t check_cublas_call(cublasStatus_t status,
                                               const char *filename, int line)
{
  if (status != CUBLAS_STATUS_SUCCESS) {
    printf("CUBLAS API failed at line %d of %s with error: %s (%d)\n",
      line, filename, cublasGetStatusString(status), status);
  }
  return status;
}

static inline cudaError_t check_cuda_last(const char *filename, int line)
{
  cudaError_t status = cudaGetLastError();
  if (status != cudaSuccess) {
    printf("CUDA API failed at line %d of %s with error: %s (%d)\n",
      line, filename, cudaGetErrorString(status), status);
  }
  return status;
}

#define CHECK_CUDA(res) \
  { if (check_cuda_call(res, __FILE__, __LINE__) != cudaSuccess) \
      return EXIT_FAILURE; }
#define CHECK_CUDA_STRICT(res) \
  { if (check_cuda_call(res, __FILE__, __LINE__) != cudaSuccess) \
      exit(EXIT_FAILURE); }
#define CHECK_CUDA_IGNORE(res) \
  { check_cuda_call(res, __FILE__, __LINE__); }

#define CHECK_CUSPARSE(res) \
  { if (check_cusparse_call(res, __FILE__, __LINE__) != CUSPARSE_STATUS_SUCCESS) \
      return EXIT_FAILURE; }
#define CHECK_CUSPARSE_STRICT(res) \
  { if (check_cusparse_call(res, __FILE__, __LINE__) != CUSPARSE_STATUS_SUCCESS) \
      exit(EXIT_FAILURE); }
#define CHECK_CUSPARSE_IGNORE(res) \
  { check_cusparse_call(res, __FILE__, __LINE__); }

#define CHECK_CUBLAS(res) \
  { if (check_cublas_call(res, __FILE__, __LINE__) != CUBLAS_STATUS_SUCCESS) \
      return EXIT_FAILURE; }
#define CHECK_CUBLAS_STRICT(res) \
  { if (check_cublas_call(res, __FILE__, __LINE__) != CUBLAS_STATUS_SUCCESS) \
      exit(EXIT_FAILURE); }
#define CHECK_CUBLAS_IGNORE(res) \
  { check_cublas_call(res, __FILE__, __LINE__); }

#define CHECK_CUDA_LAST() check_cuda_last(__FILE__, __LINE__)


#define CUPDLP_FREE_VEC(x) \
  { check_cuda_call(cudaFree(x), __FILE__, __LINE__); x = cupdlp_NULL; }

#define CUPDLP_COPY_VEC(dst, src, type, size) \
  check_cuda_call( \
    cudaMemcpy(dst, src, sizeof(type) * (size), cudaMemcpyDefault), \
    __FILE__, __LINE__)

#define CUPDLP_ZERO_VEC(var, type, size) \
  check_cuda_call( \
    cudaMemset(var, 0, sizeof(type) * (size)), __FILE__, __LINE__)

#define CUPDLP_INIT_VEC(var, size)                                             \
  {                                                                            \
    cudaError_t status = cudaMalloc((void **)&var, (size) * sizeof(__typeof__(*var))); \
    check_cuda_call(status, __FILE__, __LINE__);                               \
    if (status != cudaSuccess) goto exit_cleanup;                              \
  }

#define CUPDLP_INIT_ZERO_VEC(var, size)                                         \
  {                                                                            \
    cudaError_t status = cudaMalloc((void **)&var, (size) * sizeof(__typeof__(*var))); \
    check_cuda_call(status, __FILE__, __LINE__);                               \
    if (status != cudaSuccess) goto exit_cleanup;                              \
    status = cudaMemset(var, 0, (size) * sizeof(__typeof__(*var)));            \
    if (status != cudaSuccess) goto exit_cleanup;                              \
  }

// IG 


#define cupdlp_init_vec_int(var, size)                                   \
  {                                                                      \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(int)); \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
  }

#define cupdlp_init_vec_double(var, size)                                   \
  {                                                                         \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(double)); \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
  }

#define cupdlp_init_zero_vec_int(var, size)                              \
  {                                                                      \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(int)); \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
    status = cudaMemset(var, 0, (size) * sizeof(int));                   \
    if (status != cudaSuccess) {                                         \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",   \
             __LINE__, __FILE__, cudaGetErrorString(status), status);    \
      goto exit_cleanup;                                                 \
    }                                                                    \
  }
#define cupdlp_init_zero_vec_double(var, size)                              \
  {                                                                         \
    cudaError_t status = cudaMalloc((void**)&var, (size) * sizeof(double)); \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
    status = cudaMemset(var, 0, (size) * sizeof(double));                   \
    if (status != cudaSuccess) {                                            \
      printf("CUDA API failed at line %d of %s with error: %s (%d)\n",      \
             __LINE__, __FILE__, cudaGetErrorString(status), status);       \
      goto exit_cleanup;                                                    \
    }                                                                       \
  }

#define get_gpu_vec(vec_ptr, index, result, success)       \
  {                                                                \
    cublasStatus_t stat =                                          \
      cublasGetVector(N, sizeof(*vec_ptr), vec_ptr, 1, result, 1); \
    if (stat != CUBLAS_STATUS_SUCCESS) {                           \
      printf("data upload failed");                                \
      *success = 0;                                                 \
    } else {                                                       \
      *success = 1;                                                 \
    }                                                              \
  }

#define get_gpu_vec_element(vec_ptr, index, result, success)       \
  {                                                                \
    cublasStatus_t stat =                                          \
        cublasGetVector(1, sizeof(double), vec_ptr, 1, result, 1); \
    if (stat != CUBLAS_STATUS_SUCCESS) {                           \
      printf("data upload failed");                                \
      *success = 0;                                                 \
    } else {                                                       \
      *success = 1;                                                 \
    }                                                              \
  }

__global__ void element_primal_feas_kernel(cupdlp_float *z,
                                           const cupdlp_float *ax,
                                           const cupdlp_float *rhs,
                                           const cupdlp_float *rowScale,
                                           int ifScaled,
                                           int nEqs, int nRows);

__global__ void element_dual_feas_kernel_1(cupdlp_float *z,
                                           const cupdlp_float *aty,
                                           const cupdlp_float *cost,
                                           int nCols);

__global__ void element_dual_feas_kernel_2(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasLower,
                                           int nCols);

__global__ void element_dual_feas_kernel_3(cupdlp_float *z,
                                           const cupdlp_float *dualResidual,
                                           const cupdlp_float *hasUpper,
                                           int nCols);

__global__ void element_primal_infeas_kernel(cupdlp_float *z, const cupdlp_float *aty,
                                             const cupdlp_float *dSlackPos,
                                             const cupdlp_float *dSlackNeg,
                                             const cupdlp_float *colScale,
                                             cupdlp_float alpha, int ifScaled, int nCols);

__global__ void element_dual_infeas_kernel_lb(cupdlp_float *z,
                                              const cupdlp_float *x,
                                              const cupdlp_float *hasLower,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols);

__global__ void element_dual_infeas_kernel_ub(cupdlp_float *z,
                                              const cupdlp_float *x,
                                              const cupdlp_float *hasUpper,
                                              const cupdlp_float *colScale,
                                              cupdlp_float alpha, int ifScaled, int nCols);

__global__ void element_dual_infeas_kernel_constr(cupdlp_float *z,
                                                  const cupdlp_float *ax,
                                                  const cupdlp_float *rowScale,
                                                  cupdlp_float alpha, int ifScaled,
                                                  int nEqs, int nRows);

__global__ void element_wise_dot_kernel(cupdlp_float *x, const cupdlp_float *y, int n);

__global__ void element_wise_div_kernel(cupdlp_float *x, const cupdlp_float *y, int n);

__global__ void element_wise_projlb_kernel(cupdlp_float *x,
                                           const cupdlp_float *lb, int n);

__global__ void element_wise_projub_kernel(cupdlp_float *x,
                                           const cupdlp_float *ub, int n);

__global__ void element_wise_projSamelb_kernel(cupdlp_float *x,
                                               cupdlp_float lb, int n);

__global__ void element_wise_projSameub_kernel(cupdlp_float *x,
                                               cupdlp_float ub, int n);

__global__ void element_wise_initHaslb_kernel(cupdlp_float *haslb,
                                              const cupdlp_float *lb,
                                              cupdlp_float bound, int n);

__global__ void element_wise_initHasub_kernel(cupdlp_float *hasub,
                                              const cupdlp_float *ub,
                                              cupdlp_float bound, int n);

__global__ void element_wise_filterlb_kernel(cupdlp_float *x,
                                             const cupdlp_float *lb,
                                             cupdlp_float bound, int n);

__global__ void element_wise_filterub_kernel(cupdlp_float *x,
                                             const cupdlp_float *ub,
                                             cupdlp_float bound, int n);

__global__ void init_cuda_vec_kernel(cupdlp_float *x, cupdlp_float val, int n);

__global__ void primal_grad_step_kernel(cupdlp_float * __restrict__ xUpdate,
                                        const cupdlp_float * __restrict__ x,
                                        const cupdlp_float * __restrict__ cost,
                                        const cupdlp_float * __restrict__ ATy,
                                        const cupdlp_float * __restrict__ lb,
                                        const cupdlp_float * __restrict__ ub,
                                        cupdlp_float dPrimalStep, int nCols);

__global__ void dual_grad_step_kernel(cupdlp_float * __restrict__ yUpdate,
                                      const cupdlp_float * __restrict__ y,
                                      const cupdlp_float * __restrict__ b,
                                      const cupdlp_float * __restrict__ Ax,
                                      const cupdlp_float * __restrict__ AxUpdate,
                                      cupdlp_float dDualStep, int nRows, int nEqs);

/*
__global__ void naive_sub_kernel(cupdlp_float *z, const cupdlp_float *x,
                                 const cupdlp_float *y, int n);
*/


__global__ void movement_1_kernel(cupdlp_float * __restrict__ res_x, cupdlp_float * __restrict__ res_y,
                                  const cupdlp_float * __restrict__ xUpdate, const cupdlp_float * __restrict__ x,
                                  const cupdlp_float * __restrict__ atyUpdate, const cupdlp_float * __restrict__ aty,
                                  int nCols);

__global__ void movement_2_kernel(cupdlp_float * __restrict__ res,
                                  const cupdlp_float * __restrict__ yUpdate, const cupdlp_float * __restrict__ y,
                                  int nRows);

__global__ void sum_kernel(cupdlp_float * __restrict__ res, const cupdlp_float * __restrict__ x, int n);

#endif
