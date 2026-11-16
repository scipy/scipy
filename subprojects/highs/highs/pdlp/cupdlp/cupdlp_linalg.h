#ifndef CUPDLP_CUPDLP_LINALG_H
#define CUPDLP_CUPDLP_LINALG_H

#include "cupdlp_defs.h"
#include "cupdlp_utils.h"

#ifndef CUPDLP_CPU
#include "cuda/cupdlp_cudalinalg.cuh"
#endif

void ScatterCol(CUPDLPwork *w, cupdlp_int iCol, cupdlp_float multiplier,
                cupdlp_float *target);

void ScatterRow(CUPDLPwork *w, cupdlp_int iRow, cupdlp_float multiplier,
                cupdlp_float *target);

void AxCPU(CUPDLPwork *w, cupdlp_float *ax, const cupdlp_float *x);

void ATyCPU(CUPDLPwork *w, cupdlp_float *aty, const cupdlp_float *y);

extern double nrm2(cupdlp_int n, const double *x, cupdlp_int incx);

extern double nrminf(cupdlp_int n, const double *x, cupdlp_int incx);

double twoNorm(double *x, cupdlp_int n);

double twoNormSquared(double *x, cupdlp_int n);

double infNorm(double *x, cupdlp_int n);

cupdlp_int infNormIndex(double *x, cupdlp_int n);

/*------------------------ new added --------------------*/

double GenNorm(double *x, cupdlp_int n, cupdlp_float p);

void cupdlp_cdot(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len);

void cupdlp_cdiv(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len);

// void cupdlp_scaleVector(cupdlp_float *xout, cupdlp_float *x, cupdlp_float
// weight, const cupdlp_int len);
void cupdlp_projLowerBound(cupdlp_float *x, const cupdlp_float *lb,
                           const cupdlp_int len);
void cupdlp_projUpperBound(cupdlp_float *x, const cupdlp_float *ub,
                           const cupdlp_int len);
void cupdlp_projSameLowerBound(cupdlp_float *x, const cupdlp_float lb,
                               const cupdlp_int len);
void cupdlp_projSameUpperBound(cupdlp_float *x, const cupdlp_float ub,
                               const cupdlp_int len);
void cupdlp_projPositive(cupdlp_float *x, const cupdlp_int len);
void cupdlp_projNegative(cupdlp_float *x, const cupdlp_int len);

// void cupdlp_projLowerBound(cupdlp_float *xout, cupdlp_float *x, cupdlp_float
// *lb, const cupdlp_int len); void cupdlp_projUpperBound(cupdlp_float *xout,
// cupdlp_float *x, cupdlp_float *ub, const cupdlp_int len); void
// cupdlp_projSameLowerBound(cupdlp_float *xout, cupdlp_float *x, cupdlp_float
// lb, const cupdlp_int len); void cupdlp_projSameUpperBound(cupdlp_float *xout,
// cupdlp_float *x, cupdlp_float ub, const cupdlp_int len); void
// cupdlp_projPositive(cupdlp_float *xout, cupdlp_float *x, const cupdlp_int
// len); void cupdlp_projNegative(cupdlp_float *xout, cupdlp_float *x, const
// cupdlp_int len); cupdlp_float cupdlp_diffTwoNormSquared(cupdlp_float *x,
// cupdlp_float *y, const cupdlp_int len); cupdlp_float
// cupdlp_diffTwoNorm(cupdlp_float *x, cupdlp_float *y, const cupdlp_int len);
// cupdlp_float cupdlp_diffInfNorm(cupdlp_float *x, cupdlp_float *y, const
// cupdlp_int len); cupdlp_float cupdlp_diffDotDiff(cupdlp_float *x1,
// cupdlp_float *x2, cupdlp_float *y1, cupdlp_float *y2, const cupdlp_int len);
// void cupdlp_cdot_fb(cupdlp_float *x, const cupdlp_bool *y, const cupdlp_int
// len);

/*------------------------ new added --------------------*/

extern double dot(cupdlp_int n, const cupdlp_float *x, cupdlp_int incx,
                  const cupdlp_float *y, cupdlp_int incy);

extern double Dotprod(const cupdlp_float *x, const cupdlp_float *y, cupdlp_int n);

// todo, add this
extern double Dotprod_Neumaier(const cupdlp_float *x, const cupdlp_float *y, cupdlp_int n);

/* x = x + weight * y */
void AddToVector(cupdlp_float *x, const cupdlp_float weight,
                 const cupdlp_float *y, const cupdlp_int n);

void ScaleVector(cupdlp_float weight, cupdlp_float *x, cupdlp_int n);

// The main matrix-vector multiplication routines
// #ifndef CUPDLP_CPU
// Ax currently only works for CSC matrix multiply dense vector
// void Ax(CUPDLPwork *w, cupdlp_float *ax, const cupdlp_float *x);
// void Ax(CUPDLPwork *w, cupdlp_float *ax, void* vecAx, const cupdlp_float *x,
// void *vecX);
void Ax(CUPDLPwork *w, CUPDLPvec *ax, const CUPDLPvec *x);

// ATy currently only works for CSR matrix multiply dense vector
// void ATy(CUPDLPwork *w, cupdlp_float *aty, const cupdlp_float *y);
// void ATy(CUPDLPwork *w, cupdlp_float *aty, void *vecATy, const cupdlp_float
// *y, void *vecY);
void ATy(CUPDLPwork *w, CUPDLPvec *aty, const CUPDLPvec *y);

// #endif

/*-------------- Apis compatible with both CPU and GPU -------------------*/
// only implemented the APis need to be used on GPU

// functions in cublas
cupdlp_int cupdlp_axpy(CUPDLPwork *w, const cupdlp_int n,
                       const cupdlp_float *alpha, const cupdlp_float *x,
                       cupdlp_float *y);

cupdlp_int cupdlp_dot(CUPDLPwork *w, const cupdlp_int n, const cupdlp_float *x,
                      const cupdlp_float *y, cupdlp_float *res);

cupdlp_int cupdlp_twoNorm(CUPDLPwork *w, const cupdlp_int n,
                          const cupdlp_float *x, cupdlp_float *res);

cupdlp_int cupdlp_infNorm(CUPDLPwork *w, const cupdlp_int n,
                          const cupdlp_float *x, cupdlp_float *res);

cupdlp_int cupdlp_infNormIndex(CUPDLPwork *w, const cupdlp_int n,
                               const cupdlp_float *x, cupdlp_int *res);

cupdlp_int cupdlp_scaleVector(CUPDLPwork *w, const cupdlp_float weight,
                              cupdlp_float *x, const cupdlp_int n);

void cupdlp_twoNormSquared(CUPDLPwork *w, const cupdlp_int n,
                           const cupdlp_float *x, cupdlp_float *res);

void cupdlp_diffTwoNormSquared(CUPDLPwork *w, const cupdlp_float *x,
                               const cupdlp_float *y, const cupdlp_int len,
                               cupdlp_float *res);

void cupdlp_diffTwoNorm(CUPDLPwork *w, const cupdlp_float *x,
                        const cupdlp_float *y, const cupdlp_int len,
                        cupdlp_float *res);

void cupdlp_diffDotDiff(CUPDLPwork *w, const cupdlp_float *x1,
                        const cupdlp_float *x2, const cupdlp_float *y1,
                        const cupdlp_float *y2, const cupdlp_int len,
                        cupdlp_float *res);

// functions not in cublas
/* element wise dot: x = x .* y*/
void cupdlp_edot(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len);
/* element wise div: x = x ./ y*/
void cupdlp_ediv(cupdlp_float *x, const cupdlp_float *y, const cupdlp_int len);

void cupdlp_projlb(cupdlp_float *x, const cupdlp_float *lb,
                   const cupdlp_int len);

void cupdlp_projub(cupdlp_float *x, const cupdlp_float *ub,
                   const cupdlp_int len);

void cupdlp_projSamelb(cupdlp_float *x, const cupdlp_float lb,
                       const cupdlp_int len);

void cupdlp_projSameub(cupdlp_float *x, const cupdlp_float ub,
                       const cupdlp_int len);

/* xout = max(x, 0) */
void cupdlp_projPos(cupdlp_float *x, const cupdlp_int len);

/* xout = min(x, 0) */
void cupdlp_projNeg(cupdlp_float *x, const cupdlp_int len);

// void cupdlp_haslb(cupdlp_float *haslb, const cupdlp_float *lb,
//                   const cupdlp_float bound, const cupdlp_int len);

// void cupdlp_hasub(cupdlp_float *hasub, const cupdlp_float *ub,
//                   const cupdlp_float bound, const cupdlp_int len);

void cupdlp_filterlb(cupdlp_float *x, const cupdlp_float *lb,
                     const cupdlp_float bound, const cupdlp_int len);

void cupdlp_filterub(cupdlp_float *x, const cupdlp_float *ub,
                     const cupdlp_float bound, const cupdlp_int len);

void cupdlp_initvec(cupdlp_float *x, const cupdlp_float val,
                    const cupdlp_int len);

void cupdlp_compute_interaction_and_movement(CUPDLPwork *w,
                                             cupdlp_float *dMovement,
                                             cupdlp_float *dIteraction);

// WIP
// double get_fabs_value(double* vec, int index);
// double get_fabs_value(double* vec, int index, int N);

#endif  // CUPDLP_CUPDLP_LINALG_H
