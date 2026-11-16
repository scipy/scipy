#pragma once

// GPU portability layer for HiPDLP.
// Provides unified gpu* types and macros that map to either the CUDA or HIP
// runtime depending on whether USE_HIP is defined at compile time.

#if defined(CUPDLP_GPU) || defined(HIPDLP_GPU)

#ifdef USE_HIP
#include <hip/hip_runtime.h>
#include <hipblas/hipblas.h>
#include <hipsparse/hipsparse.h>

#define GPU_BACKEND_NAME "HIP"

using gpuStream_t = hipStream_t;
using gpuSparseHandle_t = hipsparseHandle_t;
using gpuBlasHandle_t = hipblasHandle_t;
using gpuSpMatDescr_t = hipsparseSpMatDescr_t;
using gpuDnVecDescr_t = hipsparseDnVecDescr_t;
using gpuError_t = hipError_t;
using gpuBlasStatus_t = hipblasStatus_t;
using gpuSparseStatus_t = hipsparseStatus_t;

#define gpuMalloc(ptr, size) hipMalloc(ptr, size)
#define gpuFree(ptr) hipFree(ptr)
#define gpuMemcpy(dst, src, size, kind) hipMemcpy(dst, src, size, kind)
#define gpuMemcpyAsync(dst, src, size, kind, stream) \
  hipMemcpyAsync(dst, src, size, kind, stream)
#define gpuMemset(ptr, val, size) hipMemset(ptr, val, size)
#define gpuMemsetAsync(ptr, val, size, stream) \
  hipMemsetAsync(ptr, val, size, stream)
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice

#define gpuStreamCreate(s) hipStreamCreate(s)
#define gpuStreamDestroy(s) hipStreamDestroy(s)
#define gpuStreamSynchronize(s) hipStreamSynchronize(s)
#define gpuDeviceSynchronize() hipDeviceSynchronize()
#define gpuGetLastError() hipGetLastError()
#define gpuGetErrorString(e) hipGetErrorString(e)
#define gpuSuccess hipSuccess

#define gpuBlasCreate(h) hipblasCreate(h)
#define gpuBlasDestroy(h) hipblasDestroy(h)
#define gpuBlasSetStream(h, s) hipblasSetStream(h, s)
#define gpuBlasDdot(h, n, x, ix, y, iy, r) hipblasDdot(h, n, x, ix, y, iy, r)
#define gpuBlasDnrm2(h, n, x, ix, r) hipblasDnrm2(h, n, x, ix, r)
#define gpuBlasDaxpy(h, n, a, x, ix, y, iy) hipblasDaxpy(h, n, a, x, ix, y, iy)
#define gpuBlasDscal(h, n, a, x, ix) hipblasDscal(h, n, a, x, ix)
#define GPU_BLAS_STATUS_SUCCESS HIPBLAS_STATUS_SUCCESS

#define gpuSparseCreate(h) hipsparseCreate(h)
#define gpuSparseDestroy(h) hipsparseDestroy(h)
#define gpuSparseSetStream(h, s) hipsparseSetStream(h, s)
#define gpuSparseCreateCsr(mat, r, c, nnz, rp, ci, cv, rpt, cit, vt, dt) \
  hipsparseCreateCsr(mat, r, c, nnz, rp, ci, cv, rpt, cit, vt, dt)
#define gpuSparseDestroySpMat(m) hipsparseDestroySpMat(m)
#define gpuSparseCreateDnVec(v, s, d, t) hipsparseCreateDnVec(v, s, d, t)
#define gpuSparseDestroyDnVec(v) hipsparseDestroyDnVec(v)
#define gpuSparseSpMV_bufferSize(h, op, alpha, mat, x, beta, y, dt, alg, buf) \
  hipsparseSpMV_bufferSize(h, op, alpha, mat, x, beta, y, dt, alg, buf)
#define gpuSparseSpMV(h, op, alpha, mat, x, beta, y, dt, alg, buf) \
  hipsparseSpMV(h, op, alpha, mat, x, beta, y, dt, alg, buf)
#define GPU_SPARSE_STATUS_SUCCESS HIPSPARSE_STATUS_SUCCESS
#define gpuSparseGetErrorString(s) hipsparseGetErrorString(s)

#define GPU_OPERATION_NON_TRANSPOSE HIPSPARSE_OPERATION_NON_TRANSPOSE
#define GPU_OPERATION_TRANSPOSE HIPSPARSE_OPERATION_TRANSPOSE
#define GPU_INDEX_32I HIPSPARSE_INDEX_32I
#define GPU_INDEX_BASE_ZERO HIPSPARSE_INDEX_BASE_ZERO
#define GPU_R_64F HIP_R_64F
#define GPU_SPMV_ALG_DEFAULT HIPSPARSE_SPMV_ALG_DEFAULT
#define GPU_SPMV_CSR_ALG2 HIPSPARSE_SPMV_CSR_ALG2
#define GPU_SPMV_ALG_MAIN HIPSPARSE_SPMV_CSR_ALG1

// Device info
#define gpuGetDeviceCount(n) hipGetDeviceCount(n)
#define gpuGetDeviceProperties(p, d) hipGetDeviceProperties(p, d)
using gpuDeviceProp_t = hipDeviceProp_t;

// Sparse vector update
#define gpuSparseSetDnVecValues(v, p) hipsparseDnVecSetValues(v, p)

#define gpuSparseSpMV_preprocess(h, op, a, mat, x, b, y, dt, alg, buf) \
  hipsparseSpMV_preprocess(h, op, a, mat, x, b, y, dt, alg, buf)

// Graph API
using gpuGraph_t = hipGraph_t;
using gpuGraphExec_t = hipGraphExec_t;
#define gpuStreamBeginCapture(s, m) hipStreamBeginCapture(s, m)
#define gpuStreamEndCapture(s, g) hipStreamEndCapture(s, g)
#define gpuStreamCaptureModeGlobal hipStreamCaptureModeGlobal
#define gpuGraphInstantiate(ge, g, n, e, f) hipGraphInstantiate(ge, g, n, e, f)
#define gpuGraphDestroy(g) hipGraphDestroy(g)
#define gpuGraphLaunch(ge, s) hipGraphLaunch(ge, s)
#define gpuGraphExecDestroy(ge) hipGraphExecDestroy(ge)

#else
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#define GPU_BACKEND_NAME "CUDA"

using gpuStream_t = cudaStream_t;
using gpuSparseHandle_t = cusparseHandle_t;
using gpuBlasHandle_t = cublasHandle_t;
using gpuSpMatDescr_t = cusparseSpMatDescr_t;
using gpuDnVecDescr_t = cusparseDnVecDescr_t;
using gpuError_t = cudaError_t;
using gpuBlasStatus_t = cublasStatus_t;
using gpuSparseStatus_t = cusparseStatus_t;

#define gpuMalloc(ptr, size) cudaMalloc(ptr, size)
#define gpuFree(ptr) cudaFree(ptr)
#define gpuMemcpy(dst, src, size, kind) cudaMemcpy(dst, src, size, kind)
#define gpuMemcpyAsync(dst, src, size, kind, stream) \
  cudaMemcpyAsync(dst, src, size, kind, stream)
#define gpuMemset(ptr, val, size) cudaMemset(ptr, val, size)
#define gpuMemsetAsync(ptr, val, size, stream) \
  cudaMemsetAsync(ptr, val, size, stream)
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice

#define gpuStreamCreate(s) cudaStreamCreate(s)
#define gpuStreamDestroy(s) cudaStreamDestroy(s)
#define gpuStreamSynchronize(s) cudaStreamSynchronize(s)
#define gpuDeviceSynchronize() cudaDeviceSynchronize()
#define gpuGetLastError() cudaGetLastError()
#define gpuGetErrorString(e) cudaGetErrorString(e)
#define gpuSuccess cudaSuccess

#define gpuBlasCreate(h) cublasCreate(h)
#define gpuBlasDestroy(h) cublasDestroy(h)
#define gpuBlasSetStream(h, s) cublasSetStream(h, s)
#define gpuBlasDdot(h, n, x, ix, y, iy, r) cublasDdot(h, n, x, ix, y, iy, r)
#define gpuBlasDnrm2(h, n, x, ix, r) cublasDnrm2(h, n, x, ix, r)
#define gpuBlasDaxpy(h, n, a, x, ix, y, iy) cublasDaxpy(h, n, a, x, ix, y, iy)
#define gpuBlasDscal(h, n, a, x, ix) cublasDscal(h, n, a, x, ix)
#define GPU_BLAS_STATUS_SUCCESS CUBLAS_STATUS_SUCCESS

#define gpuSparseCreate(h) cusparseCreate(h)
#define gpuSparseDestroy(h) cusparseDestroy(h)
#define gpuSparseSetStream(h, s) cusparseSetStream(h, s)
#define gpuSparseCreateCsr(mat, r, c, nnz, rp, ci, cv, rpt, cit, vt, dt) \
  cusparseCreateCsr(mat, r, c, nnz, rp, ci, cv, rpt, cit, vt, dt)
#define gpuSparseDestroySpMat(m) cusparseDestroySpMat(m)
#define gpuSparseCreateDnVec(v, s, d, t) cusparseCreateDnVec(v, s, d, t)
#define gpuSparseDestroyDnVec(v) cusparseDestroyDnVec(v)
#define gpuSparseSpMV_bufferSize(h, op, alpha, mat, x, beta, y, dt, alg, buf) \
  cusparseSpMV_bufferSize(h, op, alpha, mat, x, beta, y, dt, alg, buf)
#define gpuSparseSpMV(h, op, alpha, mat, x, beta, y, dt, alg, buf) \
  cusparseSpMV(h, op, alpha, mat, x, beta, y, dt, alg, buf)
#define GPU_SPARSE_STATUS_SUCCESS CUSPARSE_STATUS_SUCCESS
#define gpuSparseGetErrorString(s) cusparseGetErrorString(s)

#define GPU_OPERATION_NON_TRANSPOSE CUSPARSE_OPERATION_NON_TRANSPOSE
#define GPU_OPERATION_TRANSPOSE CUSPARSE_OPERATION_TRANSPOSE
#define GPU_INDEX_32I CUSPARSE_INDEX_32I
#define GPU_INDEX_BASE_ZERO CUSPARSE_INDEX_BASE_ZERO
#define GPU_R_64F CUDA_R_64F
#define GPU_SPMV_ALG_DEFAULT CUSPARSE_SPMV_ALG_DEFAULT
#define GPU_SPMV_CSR_ALG2 CUSPARSE_SPMV_CSR_ALG2
#define GPU_SPMV_ALG_MAIN CUSPARSE_SPMV_ALG_DEFAULT

// Device info
#define gpuGetDeviceCount(n) cudaGetDeviceCount(n)
#define gpuGetDeviceProperties(p, d) cudaGetDeviceProperties(p, d)
using gpuDeviceProp_t = cudaDeviceProp;

// Sparse vector update
#define gpuSparseSetDnVecValues(v, p) cusparseDnVecSetValues(v, p)

// cuSPARSE preprocess for SpMV (no-op equivalent in hipSPARSE)
#define gpuSparseSpMV_preprocess(h, op, a, mat, x, b, y, dt, alg, buf) \
  cusparseSpMV_preprocess(h, op, a, mat, x, b, y, dt, alg, buf)

// Graph API
using gpuGraph_t = cudaGraph_t;
using gpuGraphExec_t = cudaGraphExec_t;
#define gpuStreamBeginCapture(s, m) cudaStreamBeginCapture(s, m)
#define gpuStreamEndCapture(s, g) cudaStreamEndCapture(s, g)
#define gpuStreamCaptureModeGlobal cudaStreamCaptureModeGlobal
#define gpuGraphInstantiate(ge, g, n, e, f) cudaGraphInstantiate(ge, g, n, e, f)
#define gpuGraphDestroy(g) cudaGraphDestroy(g)
#define gpuGraphLaunch(ge, s) cudaGraphLaunch(ge, s)
#define gpuGraphExecDestroy(ge) cudaGraphExecDestroy(ge)

#endif  // USE_HIP

#define GPU_CHECK(call)                                               \
  do {                                                                \
    gpuError_t err = (call);                                          \
    if (err != gpuSuccess) {                                          \
      fprintf(stderr, "GPU Error at %s:%d: %s\n", __FILE__, __LINE__, \
              gpuGetErrorString(err));                                \
      exit(EXIT_FAILURE);                                             \
    }                                                                 \
  } while (0)

#define GPU_BLAS_CHECK(call)                                               \
  do {                                                                     \
    gpuBlasStatus_t status = (call);                                       \
    if (status != GPU_BLAS_STATUS_SUCCESS) {                               \
      fprintf(stderr, "GPU BLAS Error at %s:%d: %d\n", __FILE__, __LINE__, \
              (int)status);                                                \
      exit(EXIT_FAILURE);                                                  \
    }                                                                      \
  } while (0)

#define GPU_SPARSE_CHECK(call)                                               \
  do {                                                                       \
    gpuSparseStatus_t status = (call);                                       \
    if (status != GPU_SPARSE_STATUS_SUCCESS) {                               \
      fprintf(stderr, "GPU Sparse Error at %s:%d: %s\n", __FILE__, __LINE__, \
              gpuSparseGetErrorString(status));                              \
      exit(EXIT_FAILURE);                                                    \
    }                                                                        \
  } while (0)

#endif  // CUPDLP_GPU or HIPDLP_GPU
