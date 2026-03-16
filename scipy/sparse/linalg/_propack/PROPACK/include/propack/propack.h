#ifndef PROPACK_H
#define PROPACK_H

#include "types.h"


/**
 * @brief Compute a partial SVD of an m-by-n matrix (single precision) using Lanczos bidiagonalization.
 *
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of A.
 * @param n        Number of columns of A.
 * @param k        Number of singular triplets requested.
 * @param kmax     Maximum dimension of the Lanczos subspace.
 * @param aprod    Matrix–vector callback: y := A*x or A^T*x depending on transa.
 * @param U        Output array for U (size m-by-k).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array (length k) of singular values.
 * @param bnd      Output array (length k) of error bounds for sigma.
 * @param V        Output array for V (size n-by-k).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param dparm    Real parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void slansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_s aprod,
             float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
             float tolin, float* work, int lwork, int* iwork,
             float* doption, int* ioption, int* info, float* dparm, int* iparm,
             uint64_t* rng_state);

/**
 * @brief Compute a partial SVD (double precision) using Lanczos bidiagonalization.
 *
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of A.
 * @param n        Number of columns of A.
 * @param k        Number of singular triplets requested.
 * @param kmax     Maximum dimension of the Lanczos subspace.
 * @param aprod    Matrix–vector callback: y := A*x or A^T*x depending on transa.
 * @param U        Output array for U (size m-by-k).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array (length k) of singular values.
 * @param bnd      Output array (length k) of error bounds for sigma.
 * @param V        Output array for V (size n-by-k).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param dparm    Real parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void dlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_d aprod,
             double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
             double tolin, double* work, int lwork, int* iwork,
             double* doption, int* ioption, int* info, double* dparm, int* iparm,
             uint64_t* rng_state);

/**
 * @brief Compute a partial SVD (single-precision complex) using Lanczos bidiagonalization.
 *
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of A.
 * @param n        Number of columns of A.
 * @param k        Number of singular triplets requested.
 * @param kmax     Maximum dimension of the Lanczos subspace.
 * @param aprod    Matrix–vector callback: y := A*x or A^H*x depending on transa.
 * @param U        Output array for U (size m-by-k).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array (length k) of singular values (real).
 * @param bnd      Output array (length k) of error bounds (real).
 * @param V        Output array for V (size n-by-k).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param cwork    Complex workspace array.
 * @param lcwork   Length of cwork.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param cparm    Complex parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void clansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_c aprod,
             PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd, PROPACK_CPLXF_TYPE* V, int ldv,
             float tolin, float* work, int lwork, PROPACK_CPLXF_TYPE* cwork, int lcwork,
             int* iwork, float* doption, int* ioption, int* info,
             PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state);

/**
 * @brief Compute a partial SVD (double-precision complex) using Lanczos bidiagonalization.
 *
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of A.
 * @param n        Number of columns of A.
 * @param k        Number of singular triplets requested.
 * @param kmax     Maximum dimension of the Lanczos subspace.
 * @param aprod    Matrix–vector callback: y := A*x or A^H*x depending on transa.
 * @param U        Output array for U (size m-by-k).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array (length k) of singular values (real).
 * @param bnd      Output array (length k) of error bounds (real).
 * @param V        Output array for V (size n-by-k).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param zwork    Complex workspace array.
 * @param lzwork   Length of zwork.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param zparm    Complex parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void zlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_z aprod,
             PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd, PROPACK_CPLX_TYPE* V, int ldv,
             double tolin, double* work, int lwork, PROPACK_CPLX_TYPE* zwork, int lzwork,
             int* iwork, double* doption, int* ioption, int* info,
             PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state);


/**
 * @brief Implicitly restarted Lanczos SVD (single precision).
 *
 * @param which    Selector for target singular values (1/0 -> largest/smallest).
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of the implicit matrix A.
 * @param n        Number of columns of the implicit matrix A.
 * @param dim      Subspace dimension (typically p + neig).
 * @param p        Number of shifts (restart/augmentation parameter).
 * @param neig     In/out: number of requested/converged singular values.
 * @param maxiter  Maximum number of iterations.
 * @param aprod    Matrix–vector callback: y := A*x or A^T*x depending on transa.
 * @param U        Output array for U (size m-by-neig).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array of singular values (length >= neig).
 * @param bnd      Output array of error bounds (length >= neig).
 * @param V        Output array for V (size n-by-neig).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param dparm    Real parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void slansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_s aprod, float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
                 float tolin, float* work, int lwork, int* iwork, float* doption, int* ioption,
                 int* info, float* dparm, int* iparm, uint64_t* rng_state);

/**
 * @brief Implicitly restarted Lanczos SVD (double precision).
 *
 * @param which    Selector for target singular values (1/0 -> largest/smallest).
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of the implicit matrix A.
 * @param n        Number of columns of the implicit matrix A.
 * @param dim      Subspace dimension (typically p + neig).
 * @param p        Number of shifts (restart/augmentation parameter).
 * @param neig     In/out: number of requested/converged singular values.
 * @param maxiter  Maximum number of iterations.
 * @param aprod    Matrix–vector callback: y := A*x or A^T*x depending on transa.
 * @param U        Output array for U (size m-by-neig).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array of singular values (length >= neig).
 * @param bnd      Output array of error bounds (length >= neig).
 * @param V        Output array for V (size n-by-neig).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param dparm    Real parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void dlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_d aprod, double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
                 double tolin, double* work, int lwork, int* iwork, double* doption, int* ioption,
                 int* info, double* dparm, int* iparm, uint64_t* rng_state);

/**
 * @brief Implicitly restarted Lanczos SVD (single-precision complex).
 *
 * @param which    Selector for target singular values (1/0 -> largest/smallest).
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of the implicit matrix A.
 * @param n        Number of columns of the implicit matrix A.
 * @param dim      Subspace dimension (typically p + neig).
 * @param p        Number of shifts (restart/augmentation parameter).
 * @param neig     In/out: number of requested/converged singular values.
 * @param maxiter  Maximum number of iterations.
 * @param aprod    Matrix–vector callback: y := A*x or A^H*x depending on transa.
 * @param U        Output array for U (size m-by-neig).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array of singular values (real; length >= neig).
 * @param bnd      Output array of error bounds (real; length >= neig).
 * @param V        Output array for V (size n-by-neig).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param cwork    Complex workspace array.
 * @param lcwork   Length of cwork.
 * @param iwork    Integer workspace array.
 * @param soption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param cparm    Complex parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void clansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd,
                 PROPACK_CPLXF_TYPE* V, int ldv, float tolin, float* work, int lwork,
                 PROPACK_CPLXF_TYPE* cwork, int lcwork, int* iwork, float* soption,
                 int* ioption, int* info, PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state);

/**
 * @brief Implicitly restarted Lanczos SVD (double-precision complex).
 *
 * @param which    Selector for target singular values (1/0 -> largest/smallest).
 * @param jobu     Flag to compute left singular vectors U (0 no / 1 yes).
 * @param jobv     Flag to compute right singular vectors V (0 no / 1 yes).
 * @param m        Number of rows of the implicit matrix A.
 * @param n        Number of columns of the implicit matrix A.
 * @param dim      Subspace dimension (typically p + neig).
 * @param p        Number of shifts (restart/augmentation parameter).
 * @param neig     In/out: number of requested/converged singular values.
 * @param maxiter  Maximum number of iterations.
 * @param aprod    Matrix–vector callback: y := A*x or A^H*x depending on transa.
 * @param U        Output array for U (size m-by-neig).
 * @param ldu      Leading dimension of U (>= m).
 * @param sigma    Output array of singular values (real; length >= neig).
 * @param bnd      Output array of error bounds (real; length >= neig).
 * @param V        Output array for V (size n-by-neig).
 * @param ldv      Leading dimension of V (>= n).
 * @param tolin    Convergence tolerance.
 * @param work     Real workspace array.
 * @param lwork    Length of work.
 * @param zwork    Complex workspace array.
 * @param lzwork   Length of zwork.
 * @param iwork    Integer workspace array.
 * @param doption  Real options array (algorithm controls).
 * @param ioption  Integer options array.
 * @param info     Output status/return code.
 * @param zparm    Complex parameter array.
 * @param iparm    Integer parameter array.
 * @param rng_state Random number generator state (input/output).
 */
void zlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd,
                 PROPACK_CPLX_TYPE* V, int ldv, double tolin, double* work, int lwork,
                 PROPACK_CPLX_TYPE* zwork, int lzwork, int* iwork, double* doption,
                 int* ioption, int* info, PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state);

#endif
