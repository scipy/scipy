#ifndef ARNAUD_H
#define ARNAUD_H

#include "types.h"

enum ARNAUD_which {
    which_LM = 0,    /** want the NEV eigenvalues of largest magnitude.                      */
    which_SM = 1,    /** want the NEV eigenvalues of smallest magnitude.                     */
    which_LR = 2,    /** want the NEV eigenvalues of largest real part.                      */
    which_SR = 3,    /** want the NEV eigenvalues of smallest real part.                     */
    which_LI = 4,    /** want the NEV eigenvalues of largest imaginary part.                 */
    which_SI = 5,    /** want the NEV eigenvalues of smallest imaginary part.                */
    which_LA = 6,    /** compute the NEV largest (algebraic) eigenvalues. (sym)              */
    which_SA = 7,    /** compute the NEV smallest (algebraic) eigenvalues. (sym)             */
    which_BE = 8     /** compute NEV eigenvalues, half from each end of the spectrum. (sym)  */
};

/**
 * @enum ARNAUD_ido
 * @brief Flow control for reverse communication interface.
 *
 * Used to indicate the action required by the calling program.
 */
enum ARNAUD_ido {
    ido_FIRST      = 0,  /** First call                                        */
    ido_OPX        = 1,  /** OP*x needed                                       */
    ido_BX         = 2,  /** B*x needed                                        */
    ido_USER_SHIFT = 3,  /** User shifts are needed                            */
    ido_RANDOM     = 4,  /** A random vector is needed to be written in resid  */
    ido_RANDOM_OPX = 5,  /** Force random vector to be in the range of OP      */
    ido_DONE       = 99  /** Done                                              */
};


/**
 * @struct ARNAUD_state_s
 * @brief State structure for single precision ARPACK-like solvers.
 *
 * This structure holds all persistent state and workspace variables
 * required for the single precision (float) nonsymmetric and symmetric
 * eigenvalue routines. It is designed to mimic the "SAVE" variables
 * in the original Fortran ARPACK code, ensuring reentrancy and
 * thread-safety.
 *
 * Members are grouped by their usage in different parts of the program.
 */
struct ARNAUD_state_s {
    float tol;               /** problem parameter                input parameter       */
    float getv0_rnorm0;      /** getv0 internal compute           internal              */
    float aitr_betaj;        /** naitr internal compute           internal              */
    float aitr_rnorm1;       /** naitr internal compute           internal              */
    float aitr_wnorm;        /** naitr internal compute           internal              */
    float aup2_rnorm;        /** naup2 internal compute           internal              */
    enum ARNAUD_which which; /** naupd flow control               input                 */
    enum ARNAUD_ido ido;     /** naupd flow control               input/output          */
    int info;                /** problem outcome,                 input/output          */
    int bmat;                /** problem parameter, boolean       input                 */
    int mode;                /** problem parameter,               input                 */
    int n;                   /** problem parameter,               input                 */
    int ncv;                 /** problem parameter,               input                 */
    int nev;                 /** problem parameter,               input                 */
    int shift;               /** problem parameter, boolean       input                 */
    int maxiter;             /** problem parameter,               input                 */
    int nconv;               /** problem outcome,                 output                */
    int iter;                /** problem intermediate,            internal              */
    int np;                  /** problem intermediate,            internal              */
    int getv0_first;         /** getv0 flow control               internal              */
    int getv0_iter;          /** getv0 flow control               internal              */
    int getv0_itry;          /** getv0 flow control               internal              */
    int getv0_orth;          /** getv0 flow control               internal              */
    int aitr_iter;           /** naitr flow control               internal              */
    int aitr_j;              /** naitr flow control               internal              */
    int aitr_orth1;          /** naitr flow control               internal              */
    int aitr_orth2;          /** naitr flow control               internal              */
    int aitr_restart;        /** naitr flow control               internal              */
    int aitr_step3;          /** naitr flow control               internal              */
    int aitr_step4;          /** naitr flow control               internal              */
    int aitr_ierr;           /** naitr flow control               internal              */
    int aup2_initv;          /** naupd2 flow control              internal              */
    int aup2_iter;           /** naupd2 flow control              internal              */
    int aup2_getv0;          /** naupd2 flow control              internal              */
    int aup2_cnorm;          /** naupd2 flow control              internal              */
    int aup2_kplusp;         /** naupd2 flow control              internal              */
    int aup2_nev0;           /** naupd2 internal compute          internal              */
    int aup2_np0;            /** naupd2 internal compute          internal              */
    int aup2_numcnv;         /** naupd2 internal compute          internal              */
    int aup2_update;         /** naupd2 flow control              internal              */
    int aup2_ushift;         /** naupd2 flow control              internal              */
};

/**
 * @struct ARNAUD_state_d
 * @brief State structure for double precision ARPACK-like solvers.
 *
 * This structure holds all persistent state and workspace variables
 * required for the double precision (double) nonsymmetric and symmetric
 * eigenvalue routines. It is designed to mimic the "SAVE" variables
 * in the original Fortran ARPACK code, ensuring reentrancy and
 * thread-safety.
 *
 * Members are grouped by their usage in different parts of the program.
 */
struct ARNAUD_state_d {
    double tol;              /** problem parameter                input parameter       */
    double getv0_rnorm0;     /** getv0 internal compute           internal              */
    double aitr_betaj;       /** naitr internal compute           internal              */
    double aitr_rnorm1;      /** naitr internal compute           internal              */
    double aitr_wnorm;       /** naitr internal compute           internal              */
    double aup2_rnorm;       /** naup2 internal compute           internal              */
    enum ARNAUD_which which; /** naupd flow control               input                 */
    enum ARNAUD_ido ido;     /** naupd flow control               input/output          */
    int info;                /** problem outcome,                 input/output          */
    int bmat;                /** problem parameter, boolean       input                 */
    int mode;                /** problem parameter,               input                 */
    int n;                   /** problem parameter,               input                 */
    int ncv;                 /** problem parameter,               input                 */
    int nev;                 /** problem parameter,               input                 */
    int shift;               /** problem parameter, boolean       input                 */
    int maxiter;             /** problem parameter,               input                 */
    int nconv;               /** problem outcome,                 output                */
    int iter;                /** problem intermediate,            internal              */
    int np;                  /** problem intermediate,            internal              */
    int getv0_first;         /** getv0 flow control               internal              */
    int getv0_iter;          /** getv0 flow control               internal              */
    int getv0_itry;          /** getv0 flow control               internal              */
    int getv0_orth;          /** getv0 flow control               internal              */
    int aitr_iter;           /** naitr flow control               internal              */
    int aitr_j;              /** naitr flow control               internal              */
    int aitr_orth1;          /** naitr flow control               internal              */
    int aitr_orth2;          /** naitr flow control               internal              */
    int aitr_restart;        /** naitr flow control               internal              */
    int aitr_step3;          /** naitr flow control               internal              */
    int aitr_step4;          /** naitr flow control               internal              */
    int aitr_ierr;           /** naitr flow control               internal              */
    int aup2_initv;          /** naupd2 flow control              internal              */
    int aup2_iter;           /** naupd2 flow control              internal              */
    int aup2_getv0;          /** naupd2 flow control              internal              */
    int aup2_cnorm;          /** naupd2 flow control              internal              */
    int aup2_kplusp;         /** naupd2 flow control              internal              */
    int aup2_nev0;           /** naupd2 internal compute          internal              */
    int aup2_np0;            /** naupd2 internal compute          internal              */
    int aup2_numcnv;         /** naupd2 internal compute          internal              */
    int aup2_update;         /** naupd2 flow control              internal              */
    int aup2_ushift;         /** naupd2 flow control              internal              */
};



/**
 * @brief Main reverse communication loop for solving the standard or generalized
 *        eigenvalue problem in single precision real arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param resid    Array, float, n; (input/output). If info = 1 in the first call,
 *                 it should contain the initial starting vector. Otherwise holds the
 *                 residual and also used in the matvec operations.
 * @param v        Work array, float, (ldv, ncv), holds Arnoldi basis vectors (intput/output).
 * @param ldv      Integer, leading dimension of v (typically n), (input).
 * @param ipntr    Integer array, 14; For communicating pointer addresses etc (input/output).
 * @param workd    Work array, float, >= 3*n; workspace for user-supplied operations (input/output).
 * @param workl    Work array, float, >= 3*ncv*(ncv+2); internal solver workspace (input/output).
 */
void ARNAUD_snaupd(struct ARNAUD_state_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);

/**
 * @brief Main reverse communication loop for solving the standard or generalized
 *        eigenvalue problem in double precision real arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_d struct for the solver state (input/output).
 * @param resid    Array, double, n; (input/output). If info = 1 in the first call,
 *                 it should contain the initial starting vector. Otherwise holds the
 *                 residual and also used in the matvec operations.
 * @param v        Work array, double, (ldv, ncv), holds Arnoldi basis vectors (intput/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; For communicating pointer addresses etc (input/output).
 * @param workd    Work array, double, >= 3*n; workspace for user-supplied operations (input/output).
 * @param workl    Work array, double, >= 3*ncv*(ncv+2); internal solver workspace (input/output).
 */
void ARNAUD_dnaupd(struct ARNAUD_state_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

/**
 * @brief Main reverse communication loop for solving the standard or generalized
 *        eigenvalue problem in single precision complex arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param resid    Array, float complex, n; (input/output). If info = 1 in the first call,
 *                 it should contain the initial starting vector. Otherwise holds the
 *                 residual and also used in the matvec operations.
 * @param v        Work array, float complex, (ldv, ncv), holds Arnoldi basis vectors (intput/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; For communicating pointer addresses etc (input/output).
 * @param workd    Work array, float complex, >= 3*n; workspace for user-supplied operations (input/output).
 * @param workl    Work array, float complex, >= 3*ncv*(ncv+2); internal solver workspace (input/output).
 * @param rwork    Work array, float, >= 3*ncv; additional workspace for real valued computations (input/output).
 */
void ARNAUD_cnaupd(struct ARNAUD_state_s *V, ARNAUD_CPLXF_TYPE* resid, ARNAUD_CPLXF_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLXF_TYPE* workd, ARNAUD_CPLXF_TYPE* workl, float* rwork);

/**
 * @brief Main reverse communication loop for solving the standard or generalized
 *        eigenvalue problem in double precision complex arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param resid    Array, double complex, n; (input/output). If info = 1 in the first call,
 *                 it should contain the initial starting vector. Otherwise holds the
 *                 residual and also used in the matvec operations.
 * @param v        Work array, double complex, (ldv, ncv), holds Arnoldi basis vectors (intput/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; For communicating pointer addresses etc (input/output).
 * @param workd    Work array, double complex, >= 3*n; workspace for user-supplied operations (input/output).
 * @param workl    Work array, double complex, >= 3*ncv*(ncv+2); internal solver workspace (input/output).
 * @param rwork    Work array, double, >= 3*ncv; additional workspace for real valued computations (input/output).
 */
void ARNAUD_znaupd(struct ARNAUD_state_d *V, ARNAUD_CPLX_TYPE* resid, ARNAUD_CPLX_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLX_TYPE* workd, ARNAUD_CPLX_TYPE* workl, double* rwork);

/**
 * @brief Computes the eigenvalues and eigenvectors of the matrix in single precision.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer array, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param dr       Array, float, ncv; real parts of the computed eigenvalues (output).
 * @param di       Array, float, ncv; imaginary parts of the computed eigenvalues (output).
 * @param z        Work array, float, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigmar   Float, real part of the shift value (input).
 * @param sigmai   Float, imaginary part of the shift value (input).
 * @param workev   Work array, float, 3*ncv; used for storing intermediate results.
 * @param resid    Array, float, n; residual vector obtained from snaupd (input).
 * @param v        Work array, float, (ldv, ncv); holds Arnoldi basis vectors obtained from snaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from snaupd (input/output).
 * @param workd    Work array, float, >= 3*n; workspace (input).
 * @param workl    Work array, float, >= 3*ncv*(ncv+2); holds the results of snaupd (input/output).
 */
void ARNAUD_sneupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);

/**
 * @brief Computes the eigenvalues and eigenvectors of the matrix in double precision.
 *
 * @param V        Pointer to ARNAUD_state_d struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer array, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param dr       Array, double, ncv; real parts of the computed eigenvalues (output).
 * @param di       Array, double, ncv; imaginary parts of the computed eigenvalues (output).
 * @param z        Work array, double, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigmar   Float, real part of the shift value (input).
 * @param sigmai   Float, imaginary part of the shift value (input).
 * @param workev   Work array, double, 3*ncv; used for storing intermediate results.
 * @param resid    Array, double, n; residual vector obtained from dnaupd (input).
 * @param v        Work array, double, (ldv, ncv); holds Arnoldi basis vectors obtained from dnaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from dnaupd (input/output).
 * @param workd    Work array, double, >= 3*n; workspace (input).
 * @param workl    Work array, double, >= 3*ncv*(ncv+2); holds the results of dnaupd (input/output).
 */
void ARNAUD_dneupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

/**
 * @brief Computes the eigenvalues and eigenvectors of the matrix in single precision complex arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer array, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param d        Array, float complex, ncv; computed eigenvalues (output).
 * @param z        Work array, float complex, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigma    Float complex, shift value (input).
 * @param workev   Work array, float complex, 3*ncv; used for storing intermediate results (input/output).
 * @param resid    Array, float complex, n; residual vector obtained from cnaupd (input).
 * @param v        Work array, float complex, (ldv, ncv); holds Arnoldi basis vectors obtained from cnaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from cnaupd (input/output).
 * @param workd    Work array, float complex, >= 3*n; workspace (input).
 * @param workl    Work array, float complex, >= 3*ncv*(ncv+2); holds the results of cnaupd (input/output).
 * @param rwork    Work array, float, >= 3*ncv; additional workspace for real valued computations (input/output).
 */
void ARNAUD_cneupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, ARNAUD_CPLXF_TYPE* d, ARNAUD_CPLXF_TYPE* z, int ldz, ARNAUD_CPLXF_TYPE sigma, ARNAUD_CPLXF_TYPE* workev, ARNAUD_CPLXF_TYPE* resid, ARNAUD_CPLXF_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLXF_TYPE* workd, ARNAUD_CPLXF_TYPE* workl, float* rwork);

/**
 * @brief Computes the eigenvalues and eigenvectors of the matrix in double precision complex arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_d struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer array, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param d        Array, double complex, ncv; computed eigenvalues (output).
 * @param z        Work array, double complex, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigma    Double complex, shift value (input).
 * @param workev   Work array, double complex, 3*ncv; used for storing intermediate results (input/output).
 * @param resid    Array, double complex, n; residual vector obtained from znaupd (input).
 * @param v        Work array, double complex, (ldv, ncv); holds Arnoldi basis vectors obtained from znaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from znaupd (input/output).
 * @param workd    Work array, double complex, >= 3*n; workspace (input).
 * @param workl    Work array, double complex, >= 3*ncv*(ncv+2); holds the results of znaupd (input/output).
 * @param rwork    Work array, float, >= 3*ncv; additional workspace for real valued computations (input/output).
 */
void ARNAUD_zneupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, ARNAUD_CPLX_TYPE* d, ARNAUD_CPLX_TYPE* z, int ldz, ARNAUD_CPLX_TYPE sigma, ARNAUD_CPLX_TYPE* workev, ARNAUD_CPLX_TYPE* resid, ARNAUD_CPLX_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLX_TYPE* workd, ARNAUD_CPLX_TYPE* workl, double* rwork);

/**
 * @brief Computes the eigenvalues and eigenvectors of the symmetric matrix in single precision.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param resid    Array, float, n; residual vector obtained from snaupd (input).
 * @param v        Work array, float, (ldv, ncv); holds Arnoldi basis vectors obtained from snaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from snaupd (input/output).
 * @param workd    Work array, float, >= 3*n; workspace (input/output).
 * @param workl    Work array, float, >= ncv*(ncv+8); holds the results of snaupd (input/output).
 */
void ARNAUD_ssaupd(struct ARNAUD_state_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);

/**
 * @brief Computes the eigenvalues and eigenvectors of the symmetric matrix in double precision.
 *
 * @param V        Pointer to ARNAUD_state_d struct for the solver state (input/output).
 * @param resid    Array, double, n; residual vector obtained from dnaupd (input).
 * @param v        Work array, double, (ldv, ncv); holds Arnoldi basis vectors obtained from dnaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from dnaupd (input/output).
 * @param workd    Work array, double, >= 3*n; workspace (input/output).
 * @param workl    Work array, double, >= ncv*(ncv+8); holds the results of dnaupd (input/output).
 */
void ARNAUD_dsaupd(struct ARNAUD_state_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

/**
 * @brief Computes the eigenvalues and eigenvectors of the symmetric matrix in single precision arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_s struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer array, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param d        Array, float complex, ncv; computed eigenvalues (output).
 * @param z        Work array, float complex, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigma    Float complex, shift value (input).
 * @param resid    Array, float complex, n; residual vector obtained from ssaupd (input).
 * @param v        Work array, float complex, (ldv, ncv); holds Arnoldi basis vectors obtained from ssaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from ssaupd (input/output).
 * @param workd    Work array, float complex, >= 3*n; workspace (input).
 * @param workl    Work array, float complex, >= 3*ncv*(ncv+2); holds the results of ssaupd (input/output).
 */
void ARNAUD_sseupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, float* d, float* z, int ldz, float sigma, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);

/**
 * @brief Computes the eigenvalues and eigenvectors of the symmetric matrix in double precision arithmetic.
 *
 * @param V        Pointer to ARNAUD_state_d struct for the solver state (input/output).
 * @param rvec     Integer, whether to compute eigenvectors (1) or not (0).
 * @param howmny   Integer, which eigenvalues to compute, All (0), Select the ones with select array parameter (1) (input).
 * @param select   Integer, ncv; used as a boolean array for selecting which eigenvalues to compute (input).
 * @param d        Array, double complex, ncv; computed eigenvalues (output).
 * @param z        Work array, double complex, (ldz, ncv); holds the computed eigenvectors (output).
 * @param ldz      Integer, leading dimension of z (must be at least n), (input).
 * @param sigma    Double complex, shift value (input).
 * @param resid    Array, double complex, n; residual vector obtained from dsaupd (input).
 * @param v        Work array, double complex, (ldv, ncv); holds Arnoldi basis vectors obtained from dsaupd (input/output).
 * @param ldv      Integer, leading dimension of v (must be at least n), (input).
 * @param ipntr    Integer array, 14; must hold the result from dsaupd (input/output).
 * @param workd    Work array, double complex, >= 3*n; workspace (input).
 * @param workl    Work array, double complex, >= 3*ncv*(ncv+2); holds the results of dsaupd (input/output).
 */
void ARNAUD_dseupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, double* d, double* z, int ldz, double sigma, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

#endif
