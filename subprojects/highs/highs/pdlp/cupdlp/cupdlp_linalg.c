
#include "pdlp/cupdlp/cupdlp_linalg.h"

/**
 * The function `ScatterCol` performs a scatter operation on a specific
 * column of a matrix.
 *
 * @param pdhg A pointer to a structure of type pdhg.
 * @param iCol The parameter "iCol" represents the index of the column in the
 * matrix that we want to scatter.
 * @param multiplier The multiplier is a scalar value that is multiplied with
 * each element in the column of the matrix. It is used to scale the values
 * before adding them to the target array.
 * @param target The "target" parameter is a pointer to a cupdlp_float array
 * where the scattered column values will be added.
 */
void ScatterCol(CUPDLPwork *w, cupdlp_int iCol, cupdlp_float multiplier,
                cupdlp_float *target) {
  CUPDLPcsc *matrix = w->problem->data->csc_matrix;

  for (cupdlp_int p = matrix->colMatBeg[iCol]; p < matrix->colMatBeg[iCol + 1];
       ++p)
    target[matrix->colMatIdx[p]] += matrix->colMatElem[p] * multiplier;
}

void ScatterRow(CUPDLPwork *w, cupdlp_int iRow, cupdlp_float multiplier,
                cupdlp_float *target) {
  CUPDLPcsr *matrix = w->problem->data->csr_matrix;

  for (cupdlp_int p = matrix->rowMatBeg[iRow]; p < matrix->rowMatBeg[iRow + 1];
       ++p)
    target[matrix->rowMatIdx[p]] += matrix->rowMatElem[p] * multiplier;
}

void AxCPU(CUPDLPwork *w, cupdlp_float *ax, const cupdlp_float *x) {
  // #if PDHG_USE_TIMERS
  //     ++w->timers->nAxCalls;
  //     cupdlp_float dStartTime = getTimeStamp();
  // #endif

  CUPDLPproblem *lp = w->problem;

  /* no indentity currently

  FILL_ZERO(ax, lp->data->nRows);

  // [A I]*x
  for (cupdlp_int iSeq = ncols, iRow = 0; iSeq < lp->nSeq; ++iSeq, ++iRow)
  {
      if ((pdhg->lower[iSeq] > -INFINITY) && (pdhg->upper[iSeq] < INFINITY))
      {
          ax[iRow] = scaling->rowScale ? scaling->rowScale[iRow] * x[iSeq] :
  x[iSeq];
      }
      else
      {
          ax[iRow] = 0.0;
      }
  }
  */

  memset(ax, 0, sizeof(cupdlp_float) * lp->data->nRows);

  for (cupdlp_int iCol = 0; iCol < lp->data->nCols; ++iCol) {
    ScatterCol(w, iCol, x[iCol], ax);
  }

  // #if PDHG_USE_TIMERS
  //     w->timers->dAxTime += getTimeStamp() - dStartTime;
  // #endif
}

void ATyCPU(CUPDLPwork *w, cupdlp_float *aty, const cupdlp_float *y) {
  // #if PDHG_USE_TIMERS
  //     ++w->timers->nAtyCalls;
  //     cupdlp_float dStartTime = getTimeStamp();
  // #endif

  CUPDLPproblem *lp = w->problem;

  /* no indentity currently
  FILL_ZERO(aty, lp->nSeq);

  // [A I]'*y
  for (cupdlp_int iSeq = ncols, iRow = 0; iSeq < lp->nSeq; ++iSeq, ++iRow)
  {
      ScatterRow(pdhg, iRow, y[iRow], aty);

      if ((pdhg->lower[iSeq] > -INFINITY) && (pdhg->upper[iSeq] < INFINITY))
      {
          aty[iSeq] = (scaling->rowScale ? scaling->rowScale[iRow] * y[iRow] :
  y[iRow]);
      }
      else
      {
          aty[iSeq] = 0.0;
      }
  }
  */

  memset(aty, 0, sizeof(cupdlp_float) * lp->data->nCols);
  for (cupdlp_int iRow = 0; iRow < lp->data->nRows; ++iRow) {
    ScatterRow(w, iRow, y[iRow], aty);
  }

  // #if PDHG_USE_TIMERS
  //     w->timers->dAtyTime += getTimeStamp() - dStartTime;
  // #endif
}

double nrm2(cupdlp_int n, const double *x, cupdlp_int incx) {
#ifdef USE_MY_BLAS
  assert(incx == 1);

  double nrm = 0.0;

  for (int i = 0; i < n; ++i) {
    nrm += x[i] * x[i];
  }

  return sqrt(nrm);
#else
  return dnrm2(n, x, incx);
#endif
}

double nrminf(cupdlp_int n, const double *x, cupdlp_int incx) {
#ifdef USE_MY_BLAS
  assert(incx == 1);

  double nrm = 0.0;

  for (int i = 0; i < n; ++i) {
    double tmp = fabs(x[i]);
    if (tmp > nrm) nrm = tmp;
  }

  return nrm;
#else
  return dnrm2(n, x, incx);
#endif
}

static cupdlp_int nrminfindex(cupdlp_int n, const double *x, cupdlp_int incx) {
#ifdef USE_MY_BLAS
  assert(incx == 1);

  double nrm = 0.0;
  cupdlp_int index = 0;

  for (int i = 0; i < n; ++i) {
    double tmp = fabs(x[i]);
    if (tmp > nrm) {
      nrm = tmp;
      index = i;
    }
  }

  return index;
#else
  return dnrminfindex(n, x, incx);
#endif
}

double twoNorm(double *x, cupdlp_int n) { return nrm2(n, x, 1); }

double twoNormSquared(double *x, cupdlp_int n) { return pow(twoNorm(x, n), 2); }

double infNorm(double *x, cupdlp_int n) { return nrminf(n, x, 1); }

/*------------------------ new added --------------------*/
double GenNorm(double *x, cupdlp_int n, cupdlp_float p) {
  if (p == 2.0) {
    return twoNorm(x, n);
  } else if (p == INFINITY) {
    return infNorm(x, n);
  } else {
    double nrm = 0.0;

    for (int i = 0; i < n; ++i) {
      nrm += pow(fabs(x[i]), p);
    }

    return pow(nrm, 1.0 / p);
  }
}

/* x = x .* y*/
void cupdlp_cdot(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] *= y[i];
  }
}

/* x = x ./ y*/
void cupdlp_cdiv(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] /= y[i];
  }
}

/* xout = weight * x */
// void cupdlp_scaleVector(cupdlp_float *xout, cupdlp_float *x, cupdlp_float
// weight, const cupdlp_int len)
// {
//     for (int i = 0; i < len; i++)
//     {
//         xout[i] = weight * x[i];
//     }
// }

/* xout = max(x, lb), lb is vector */
void cupdlp_projLowerBound(cupdlp_float *x, const cupdlp_float *lb,
                           const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = x[i] > lb[i] ? x[i] : lb[i];
  }
}

/* xout = min(x, ub), ub is vector */
void cupdlp_projUpperBound(cupdlp_float *x, const cupdlp_float *ub,
                           const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = x[i] < ub[i] ? x[i] : ub[i];
  }
}

/* xout = max(x, lb), lb is number */
void cupdlp_projSameLowerBound(cupdlp_float *x, const cupdlp_float lb,
                               const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = x[i] > lb ? x[i] : lb;
  }
}

/* xout = min(x, ub), ub is number */
void cupdlp_projSameUpperBound(cupdlp_float *x, const cupdlp_float ub,
                               const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = x[i] < ub ? x[i] : ub;
  }
}

/* xout = max(x, 0) */
void cupdlp_projPositive(cupdlp_float *x, const cupdlp_int len) {
  cupdlp_projSameLowerBound(x, 0.0, len);
}

/* xout = min(x, 0) */
void cupdlp_projNegative(cupdlp_float *x, const cupdlp_int len) {
  cupdlp_projSameUpperBound(x, 0.0, len);
}

// Unused warning
/* ||x - y||_2^2 */
// cupdlp_float cupdlp_diffTwoNormSquared(cupdlp_float *x, cupdlp_float *y,
// const cupdlp_int len)
// static cupdlp_float diffTwoNormSquared(cupdlp_float *x, cupdlp_float *y,
//                                 const cupdlp_int len) {
//   cupdlp_float res = 0.0;
//   for (int i = 0; i < len; i++) {
//     cupdlp_float tmp = x[i] - y[i];
//     res += tmp * tmp;
//   }
//   return res;
// }

// Unused warning
/* ||x - y||_2 */
// cupdlp_float cupdlp_diffTwoNorm(cupdlp_float *x, cupdlp_float *y, const
// cupdlp_int len)
// static cupdlp_float diffTwoNorm(cupdlp_float *x, cupdlp_float *y,
//                          const cupdlp_int len) {
//   // return sqrt(cupdlp_diffTwoNormSquared(x, y, len));
//   return sqrt(diffTwoNormSquared(x, y, len));
// }

// Unused warning
/* ||x - y||_inf */
// cupdlp_float cupdlp_diffInfNorm(cupdlp_float *x, cupdlp_float *y, const
// cupdlp_int len)
// static cupdlp_float diffInfNorm(cupdlp_float *x, cupdlp_float *y,
//                          const cupdlp_int len) {
//   cupdlp_float res = 0.0;
//   for (int i = 0; i < len; i++) {
//     cupdlp_float tmp = fabs(x[i] - y[i]);
//     if (tmp > res) res = tmp;
//   }
//   return res;
// }

// Unused warning
/* (x1 - x2)' (y1 - y2) */
// cupdlp_float cupdlp_diffDotDiff(cupdlp_float *x1, cupdlp_float *x2,
// cupdlp_float *y1, cupdlp_float *y2, const cupdlp_int len)
// static cupdlp_float diffDotDiff(cupdlp_float *x1, cupdlp_float *x2, cupdlp_float *y1,
//                          cupdlp_float *y2, const cupdlp_int len) {
//   cupdlp_float x1y1 = dot(len, x1, 1, y1, 1);
//   cupdlp_float x2y2 = dot(len, x2, 1, y2, 1);
//   cupdlp_float x1y2 = dot(len, x1, 1, y2, 1);
//   cupdlp_float x2y1 = dot(len, x2, 1, y1, 1);

//   return x1y1 - x1y2 - x2y1 + x2y2;
// }

/* x = x .* y */
// void cupdlp_cdot_fb(cupdlp_float *x, const cupdlp_bool *y, const cupdlp_int
// len)
// {
//     for (int i = 0; i < len; i++)
//     {
//         x[i] *= y[i];
//     }
// }

/*------------------------ new added --------------------*/

double dot(cupdlp_int n, const cupdlp_float *x, cupdlp_int incx, const cupdlp_float *y,
           cupdlp_int incy) {
#ifdef USE_MY_BLAS
  assert(incx == 1 && incy == 1);

  double dres = 0.0;

  for (int i = 0; i < n; ++i) {
    dres += x[i] * y[i];
  }

  return dres;
#else
  return ddot(n, x, incx, y, incy);
#endif
}

double Dotprod(const cupdlp_float *x, const cupdlp_float *y, cupdlp_int n) {
  return dot(n, x, 1, y, 1);
}

double Dotprod_Neumaier(const cupdlp_float *x, const cupdlp_float *y, cupdlp_int n) {
  return dot(n, x, 1, y, 1);
}

/* x = x + weight * y */
void AddToVector(cupdlp_float *x, const cupdlp_float weight,
                 const cupdlp_float *y, const cupdlp_int n) {
#ifdef USE_MY_BLAS

  for (int i = 0; i < n; ++i) {
    x[i] += weight * y[i];
  }

#else
  return ddot(n, x, incx, y, incy);
#endif
}

/* x = weight * x */
void ScaleVector(cupdlp_float weight, cupdlp_float *x, cupdlp_int n) {
#ifdef USE_MY_BLAS

  for (int i = 0; i < n; ++i) {
    x[i] *= weight;
  }

#else
  return ddot(n, x, incx, y, incy);
#endif
}


static void cupdlp_filter_lower_bound(cupdlp_float *x, const cupdlp_float *lb,
                               const cupdlp_float bound, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = lb[i] > bound ? lb[i] : 0.0;
  }
}

static void cupdlp_filter_upper_bound(cupdlp_float *x, const cupdlp_float *ub,
                               const cupdlp_float bound, const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = ub[i] < bound ? ub[i] : 0.0;
  }
}

static void cupdlp_init_vector(cupdlp_float *x, const cupdlp_float val,
                        const cupdlp_int len) {
  for (int i = 0; i < len; i++) {
    x[i] = val;
  }
}

#ifndef CUPDLP_CPU

static void Ax_single_gpu(CUPDLPwork *w, cusparseDnVecDescr_t vecX,
                   cusparseDnVecDescr_t vecAx) {
  cupdlp_float begin = getTimeStamp();
  cupdlp_float alpha = 1.0;
  cupdlp_float beta = 0.0;

  switch (w->problem->data->matrix_format) {
    case CSR_CSC:
      cuda_csr_Ax(w->cusparsehandle, w->problem->data->csr_matrix->cuda_csr,
                  vecX, vecAx, w->dBuffer_csr_Ax, alpha, beta);
      break;
    case CSC:
      cupdlp_printf("Error: Ax_single_gpu requires CSR matrix\n");
      exit(1);
      break;
    case CSR:
      cuda_csr_Ax(w->cusparsehandle, w->problem->data->csr_matrix->cuda_csr,
                  vecX, vecAx, w->dBuffer_csr_Ax, alpha, beta);
      break;
    default:
      cupdlp_printf("Error: Unknown matrix format in Ax_single_gpu\n");
      exit(1);
  }
  w->timers->DeviceMatVecProdTime += getTimeStamp() - begin;
}

static void Ax_multi_gpu(CUPDLPdata *d, cupdlp_float *ax, const cupdlp_float *x) {
  cupdlp_printf("Error: Ax_multi_gpu not implemented\n");
  exit(1);
}

static void ATy_single_gpu(CUPDLPwork *w, cusparseDnVecDescr_t vecY,
                    cusparseDnVecDescr_t vecATy) {
  cupdlp_float begin = getTimeStamp();

  cupdlp_float alpha = 1.0;
  cupdlp_float beta = 0.0;

  switch (w->problem->data->matrix_format) {
    case CSR_CSC:
      cuda_csc_ATy(w->cusparsehandle, w->problem->data->csc_matrix->cuda_csc,
                   vecY, vecATy, w->dBuffer_csc_ATy, alpha, beta);
      break;
    case CSC:
      cuda_csc_ATy(w->cusparsehandle, w->problem->data->csc_matrix->cuda_csc,
                   vecY, vecATy, w->dBuffer_csc_ATy, alpha, beta);
      break;
    case CSR:
      cupdlp_printf("Error: ATy_single_gpu requires CSC matrix\n");
      exit(1);
      break;
    default:
      printf("Error: Unknown matrix format in Ax_single_gpu\n");
      exit(1);
  }

  w->timers->DeviceMatVecProdTime += getTimeStamp() - begin;
}

static void ATy_multi_gpu(CUPDLPdata *d, cupdlp_float *aty, const cupdlp_float *y) {
  cupdlp_printf("Error: ATy_multi_gpu not implemented\n");
  exit(1);
}

#endif

void Ax(CUPDLPwork *w, CUPDLPvec *ax, const CUPDLPvec *x) {
  cupdlp_float begin = getTimeStamp();

  CUPDLPdata *d = w->problem->data;
  switch (d->device) {
    case CPU:
      AxCPU(w, ax->data, x->data);
      break;
    case SINGLE_GPU:

#ifndef CUPDLP_CPU
      Ax_single_gpu(w, x->cuda_vec, ax->cuda_vec);
#else
      printf("GPU not supported in CPU build\n");
      exit(1);
#endif
      break;
    case MULTI_GPU:
#ifndef CUPDLP_CPU
      Ax_multi_gpu(d, ax->data, x->data);
#else
      printf("GPU not supported in CPU build\n");
      exit(1);
#endif
      break;
    default:
      printf("Error: Unknown device type in Ax\n");
      exit(1);
  }

#if PDHG_USE_TIMERS
  w->timers->dAxTime += getTimeStamp() - begin;
  w->timers->nAxCalls++;
#endif
}

void ATy(CUPDLPwork *w, CUPDLPvec *aty, const CUPDLPvec *y)

{
  cupdlp_float begin = getTimeStamp();

  CUPDLPdata *d = w->problem->data;
  switch (d->device) {
    case CPU:
      ATyCPU(w, aty->data, y->data);
      break;
    case SINGLE_GPU:
#ifndef CUPDLP_CPU
      ATy_single_gpu(w, y->cuda_vec, aty->cuda_vec);
#else
      printf("GPU not supported in CPU build\n");
      exit(1);
#endif
      break;
    case MULTI_GPU:
#ifndef CUPDLP_CPU
      ATy_multi_gpu(d, aty->data, y->data);
#else
      printf("GPU not supported in CPU build\n");
      exit(1);
#endif
      break;
    default:
      printf("Error: Unknown device type in ATy\n");
      exit(1);
  }
#if PDHG_USE_TIMERS
  w->timers->dAtyTime += getTimeStamp() - begin;
  w->timers->nAtyCalls++;
#endif
}

/*-------------- Apis compatible with both CPU and GPU -------------------*/
// only implemented the APis need to be used on GPU

// functions in cublas

cupdlp_int cupdlp_axpy(CUPDLPwork *w, const cupdlp_int n,
                       const cupdlp_float *alpha, const cupdlp_float *x,
                       cupdlp_float *y) {
#ifndef CUPDLP_CPU
#ifndef SFLOAT
  CHECK_CUBLAS(cublasDaxpy(w->cublashandle, n, alpha, x, 1, y, 1));
#else
  CHECK_CUBLAS(cublasSaxpy(w->cublashandle, n, alpha, x, 1, y, 1));
#endif
#else
  // AddToVector(x, *alpha, y, n);
  AddToVector(y, *alpha, x, n);
#endif
  return 0;
}

cupdlp_int cupdlp_dot(CUPDLPwork *w, const cupdlp_int n, const cupdlp_float *x,
                      const cupdlp_float *y, cupdlp_float *res) {
#ifndef CUPDLP_CPU
#ifndef SFLOAT
  CHECK_CUBLAS(cublasDdot(w->cublashandle, n, x, 1, y, 1, res));
#else
  CHECK_CUBLAS(cublasSdot(w->cublashandle, n, x, 1, y, 1 res));
#endif
#else
  *res = dot(n, x, 1, y, 1);
#endif
  return 0;
}

cupdlp_int cupdlp_twoNorm(CUPDLPwork *w, const cupdlp_int n,
                          const cupdlp_float *x, cupdlp_float *res) {
#ifndef CUPDLP_CPU
#ifndef SFLOAT
  CHECK_CUBLAS(cublasDnrm2(w->cublashandle, n, x, 1, res));
#else
  CHECK_CUBLAS(cublasSnrm2(w->cublashandle, n, x, 1, res));
#endif
#else
  *res = nrm2(n, x, 1);
#endif
  return 0;
}

cupdlp_int cupdlp_infNormIndex(CUPDLPwork *w, const cupdlp_int n,
                               const cupdlp_float *x, cupdlp_int *res) {
#ifndef CUPDLP_CPU
#ifndef SFLOAT
  CHECK_CUBLAS(cublasIdamax(w->cublashandle, n, x, 1, res));
#else
  CHECK_CUBLAS(cublasIsamax(w->cublashandle, n, x, 1, res));
#endif
#else
  *res = nrminfindex(n, x, 1);
#endif
  return 0;
}

cupdlp_int cupdlp_scaleVector(CUPDLPwork *w, const cupdlp_float weight,
                              cupdlp_float *x, const cupdlp_int n) {
#ifndef CUPDLP_CPU
#ifndef SFLOAT
  CHECK_CUBLAS(cublasDscal(w->cublashandle, n, &weight, x, 1));
#else
  CHECK_CUBLAS(cublasSscal(w->cublashandle, n, &weight, x, 1));
#endif
#else
  ScaleVector(weight, x, n);
#endif
  return 0;
}

void cupdlp_twoNormSquared(CUPDLPwork *w, const cupdlp_int n,
                           const cupdlp_float *x, cupdlp_float *res) {
  cupdlp_dot(w, n, x, x, res);
}

/* ||x - y||_2^2 */
void cupdlp_diffTwoNormSquared(CUPDLPwork *w, const cupdlp_float *x,
                               const cupdlp_float *y, const cupdlp_int len,
                               cupdlp_float *res) {
  CUPDLP_COPY_VEC(w->buffer2, x, cupdlp_float, len);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(w, len, &alpha, y, w->buffer2);
  cupdlp_twoNormSquared(w, len, w->buffer2, res);
}

/* ||x - y||_2 */
void cupdlp_diffTwoNorm(CUPDLPwork *w, const cupdlp_float *x,
                        const cupdlp_float *y, const cupdlp_int len,
                        cupdlp_float *res) {
  CUPDLP_COPY_VEC(w->buffer2, x, cupdlp_float, len);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(w, len, &alpha, y, w->buffer2);
  cupdlp_twoNorm(w, len, w->buffer2, res);
}

/* (x1 - x2)' (y1 - y2) */
void cupdlp_diffDotDiff(CUPDLPwork *w, const cupdlp_float *x1,
                        const cupdlp_float *x2, const cupdlp_float *y1,
                        const cupdlp_float *y2, const cupdlp_int len,
                        cupdlp_float *res) {
  CUPDLP_COPY_VEC(w->buffer2, x1, cupdlp_float, len);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(w, len, &alpha, x2, w->buffer2);
  CUPDLP_COPY_VEC(w->buffer3, y1, cupdlp_float, len);
  cupdlp_axpy(w, len, &alpha, y2, w->buffer3);
  // reduce step
  cupdlp_dot(w, len, w->buffer2, w->buffer3, res);
}

// functions not in cublas

/* element wise dot: x = x .* y*/
void cupdlp_edot(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_edot_cuda(x, y, len);
#else
  cupdlp_cdot(x, y, len);
#endif
}

/* element wise div: x = x ./ y*/
void cupdlp_ediv(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_ediv_cuda(x, y, len);
#else
  cupdlp_cdiv(x, y, len);
#endif
}

void cupdlp_projlb(cupdlp_float *x, const cupdlp_float *lb,
                   const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_projlb_cuda(x, lb, len);
#else
  cupdlp_projLowerBound(x, lb, len);
#endif
}

void cupdlp_projub(cupdlp_float *x, const cupdlp_float *ub,
                   const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_projub_cuda(x, ub, len);
#else
  cupdlp_projUpperBound(x, ub, len);
#endif
}

void cupdlp_projSamelb(cupdlp_float *x, const cupdlp_float lb,
                       const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_projSamelb_cuda(x, lb, len);
#else
  cupdlp_projSameLowerBound(x, lb, len);
#endif
}

void cupdlp_projSameub(cupdlp_float *x, const cupdlp_float ub,
                       const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_projSameub_cuda(x, ub, len);
#else
  cupdlp_projSameUpperBound(x, ub, len);
#endif
}

/* xout = max(x, 0) */
void cupdlp_projPos(cupdlp_float *x, const cupdlp_int len) {
  cupdlp_projSamelb(x, 0.0, len);
}

/* xout = min(x, 0) */
void cupdlp_projNeg(cupdlp_float *x, const cupdlp_int len) {
  cupdlp_projSameub(x, 0.0, len);
}

// void cupdlp_haslb(cupdlp_float *haslb, const cupdlp_float *lb,
//                   const cupdlp_float bound, const cupdlp_int len) {
// #ifndef CUPDLP_CPU
//   cupdlp_haslb_cuda(haslb, lb, bound, len);
// #else
//   cupdlp_hasLower(haslb, lb, bound, len);
// #endif
// }

// void cupdlp_hasub(cupdlp_float *hasub, const cupdlp_float *ub,
//                   const cupdlp_float bound, const cupdlp_int len) {
// #ifndef CUPDLP_CPU
//   cupdlp_hasub_cuda(hasub, ub, bound, len);
// #else
//   cupdlp_hasUpper(hasub, ub, bound, len);
// #endif
// }

void cupdlp_filterlb(cupdlp_float *x, const cupdlp_float *lb,
                     const cupdlp_float bound, const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_filterlb_cuda(x, lb, bound, len);
#else
  cupdlp_filter_lower_bound(x, lb, bound, len);
#endif
}

void cupdlp_filterub(cupdlp_float *x, const cupdlp_float *ub,
                     const cupdlp_float bound, const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_filterub_cuda(x, ub, bound, len);
#else
  cupdlp_filter_upper_bound(x, ub, bound, len);
#endif
}

void cupdlp_initvec(cupdlp_float *x, const cupdlp_float val,
                    const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_initvec_cuda(x, val, len);
#else
  cupdlp_init_vector(x, val, len);
#endif
}

/*
void cupdlp_sub(cupdlp_float *xout, const cupdlp_float *x1,
                const cupdlp_float *x2, const cupdlp_int len) {
#ifndef CUPDLP_CPU
  cupdlp_sub_cuda(xout, x1, x2, len);
#else
  CUPDLP_COPY_VEC(xout, x1, cupdlp_float, len);
  cupdlp_float alpha = -1.0;
  cupdlp_axpy(NULL, len, &alpha, x2, xout);
#endif
}
*/

void cupdlp_compute_interaction_and_movement(CUPDLPwork *w,
                                             cupdlp_float *dMovement,
                                             cupdlp_float *dInteraction) {
  CUPDLPiterates *iterates = w->iterates;
  cupdlp_int nCols = w->problem->nCols;
  cupdlp_int nRows = w->problem->nRows;
  cupdlp_float beta = sqrt(w->stepsize->dBeta);
  cupdlp_float dX = 0.0;
  cupdlp_float dY = 0.0;

  cupdlp_int iter = w->timers->nIter;
  CUPDLPvec *x = iterates->x[iter % 2];
  CUPDLPvec *y = iterates->y[iter % 2];
  CUPDLPvec *aty = iterates->aty[iter % 2];
  CUPDLPvec *xUpdate = iterates->x[(iter + 1) % 2];
  CUPDLPvec *yUpdate = iterates->y[(iter + 1) % 2];
  CUPDLPvec *atyUpdate = iterates->aty[(iter + 1) % 2];

#if !defined(CUPDLP_CPU) && USE_KERNELS
  cupdlp_movement_interaction_cuda(&dX, &dY, dInteraction, w->buffer2,
      xUpdate->data, x->data, yUpdate->data, y->data, atyUpdate->data, aty->data, nRows, nCols);
#else
    cupdlp_diffTwoNormSquared(w, x->data, xUpdate->data, nCols, &dX);
    cupdlp_diffTwoNormSquared(w, y->data, yUpdate->data, nRows, &dY);
    //      Δx' (AΔy)
    cupdlp_diffDotDiff(w, x->data, xUpdate->data, aty->data, atyUpdate->data,
                       nCols, dInteraction);
#endif
  *dMovement = dX * 0.5 * beta + dY / (2.0 * beta);
}

