#ifndef _ARPACK_H
#define _ARPACK_H

#include <math.h>
#include <complex.h>

#if defined(_MSC_VER)
    // MSVC definitions
    typedef _Dcomplex ARPACK_CPLX_TYPE;
    typedef _Fcomplex ARPACK_CPLXF_TYPE;

#else
    // C99 compliant compilers
    typedef double complex ARPACK_CPLX_TYPE;
    typedef float complex ARPACK_CPLXF_TYPE;

#endif



/*
 * ARPACK uses the so-called reverse-communication style that typically exits
 * the program with its in/out arguments to signal, in what stage the algorithm
 * is and what it needs. Then user modifies the arguments and calls again with
 * the necessary information. Thus the state of the whole program is sent back
 * and forth through in-place modified arguments. On top of this, ARPACK also
 * uses lots of variables through the Fortran's dreadful SAVE attribute that
 * persists the variable values across runs. Instead we move all those variables
 * into the reverse communication layer.
 *
 * For scalar arguments we use a C struct <-> Python dict bridge and for array
 * arguments we use NumPy array pointers that are originated from the user side
 * to modify things in-place and not to deal with ref counting or alloc/free.
 *
 * To generate random vectors that are used in ARPACK algorithm, we also expect
 * a NumPy generator object from the user for reproducible runs. However this
 * can be replaced with a different number generation routine.
*/

enum ARPACK_which {
    which_LM = 0,    // want the NEV eigenvalues of largest magnitude.
    which_SM = 1,    // want the NEV eigenvalues of smallest magnitude.
    which_LR = 2,    // want the NEV eigenvalues of largest real part.
    which_SR = 3,    // want the NEV eigenvalues of smallest real part.
    which_LI = 4,    // want the NEV eigenvalues of largest imaginary part.
    which_SI = 5,    // want the NEV eigenvalues of smallest imaginary part.
    which_LA = 6,    // compute the NEV largest (algebraic) eigenvalues. (sym)
    which_SA = 7,    // compute the NEV smallest (algebraic) eigenvalues. (sym)
    which_BE = 8     // compute NEV eigenvalues, half from each end of the spectrum. (sym)
};

enum ARPACK_ido {
    ido_FIRST      = 0,  // First call
    ido_OPX        = 1,  // OP*x needed
    ido_BX         = 2,  // B*x needed
    ido_USER_SHIFT = 3,  // User shifts are needed
    ido_RANDOM     = 4,  // A random vector is needed to be written in resid
    ido_RANDOM_OPX = 5,  // Force random vector to be in the range of OP
    ido_DONE       = 99  // Done
};

/**
 * With the following structs, we collect all "SAVE"d Fortran variables to track
 * the problem and avoid reentry issues. It is not the cleanest and is laborious
 * but otherwise reentrancy is compromised. There are additional variables in the
 * original Fortran code that are also "SAVE"d however upon inspection, they
 * are assigned and then used in the same call and thus used without saving.
**/

struct ARPACK_arnoldi_update_vars_s {
    float tol;               // problem parameter
    float getv0_rnorm0;      // getv0 internal compute
    float aitr_betaj;        // naitr internal compute
    float aitr_rnorm1;       // naitr internal compute
    float aitr_wnorm;        // naitr internal compute
    float aup2_rnorm;        // naup2 internal compute
    enum ARPACK_ido ido;     // naupd flow control
    enum ARPACK_which which; // naupd flow control
    int bmat;                // problem parameter,          boolean
    int info;                // problem outcome,            integer
    int iter;                // problem intermediate,       integer
    int maxiter;             // problem parameter,          integer
    int mode;                // problem parameter,          integer
    int n;                   // problem parameter,          integer
    int nconv;               // problem outcome,            integer
    int ncv;                 // problem parameter,          integer
    int nev;                 // problem parameter,          integer
    int np;                  // problem intermediate,       integer
    int numop;               // problem intermediate,       integer
    int numpb;               // problem intermediate,       integer
    int numreo;              // problem intermediate,       integer
    int shift;               // problem parameter,          boolean
    int getv0_first;         // getv0 flow control
    int getv0_iter;          // getv0 flow control
    int getv0_itry;          // getv0 flow control
    int getv0_orth;          // getv0 flow control
    int aitr_iter;           // naitr flow control
    int aitr_j;              // naitr flow control
    int aitr_orth1;          // naitr flow control
    int aitr_orth2;          // naitr flow control
    int aitr_restart;        // naitr flow control
    int aitr_step3;          // naitr flow control
    int aitr_step4;          // naitr flow control
    int aitr_ierr;           // naitr flow control
    int aup2_initv;          // naupd2 flow control
    int aup2_iter;           // naupd2 flow control
    int aup2_getv0;          // naupd2 flow control
    int aup2_cnorm;          // naupd2 flow control
    int aup2_kplusp;         // naupd2 flow control
    int aup2_nev0;           // naupd2 internal compute
    int aup2_np0;            // naupd2 internal compute
    int aup2_numcnv;         // naupd2 internal compute
    int aup2_update;         // naupd2 flow control
    int aup2_ushift;         // naupd2 flow control
};


struct ARPACK_arnoldi_update_vars_d {
    double tol;              // problem parameter
    double getv0_rnorm0;     // getv0 internal compute
    double aitr_betaj;       // naitr internal compute
    double aitr_rnorm1;      // naitr internal compute
    double aitr_wnorm;       // naitr internal compute
    double aup2_rnorm;       // naup2 internal compute
    enum ARPACK_ido ido;     // naupd flow control
    enum ARPACK_which which; // naupd flow control
    int bmat;                // problem parameter,    boolean
    int info;                // problem outcome,      integer
    int iter;                // problem intermediate, integer
    int maxiter;             // problem parameter,    integer
    int mode;                // problem parameter,    integer
    int n;                   // problem parameter,    integer
    int nconv;               // problem outcome,      integer
    int ncv;                 // problem parameter,    integer
    int nev;                 // problem parameter,    integer
    int np;                  // problem intermediate, integer
    int numop;               // problem intermediate, integer
    int numpb;               // problem intermediate, integer
    int numreo;              // problem intermediate, integer
    int shift;               // problem parameter,    boolean
    int getv0_first;         // getv0 flow control
    int getv0_iter;          // getv0 flow control
    int getv0_itry;          // getv0 flow control
    int getv0_orth;          // getv0 flow control
    int aitr_iter;           // naitr flow control
    int aitr_j;              // naitr flow control
    int aitr_orth1;          // naitr flow control
    int aitr_orth2;          // naitr flow control
    int aitr_restart;        // naitr flow control
    int aitr_step3;          // naitr flow control
    int aitr_step4;          // naitr flow control
    int aitr_ierr;           // naitr flow control
    int aup2_initv;          // naupd2 flow control
    int aup2_iter;           // naupd2 flow control
    int aup2_getv0;          // naupd2 flow control
    int aup2_cnorm;          // naupd2 flow control
    int aup2_kplusp;         // naupd2 flow control
    int aup2_nev0;           // naupd2 internal compute
    int aup2_np0;            // naupd2 internal compute
    int aup2_numcnv;         // naupd2 internal compute
    int aup2_update;         // naupd2 flow control
    int aup2_ushift;         // naupd2 flow control
};


void snaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dnaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void cnaupd(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork);
void znaupd(struct ARPACK_arnoldi_update_vars_d *V, ARPACK_CPLX_TYPE* resid, ARPACK_CPLX_TYPE* v, int ldv, int* ipntr, ARPACK_CPLX_TYPE* workd, ARPACK_CPLX_TYPE* workl, double* rwork);

void sneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dneupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void cneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, ARPACK_CPLXF_TYPE* d, ARPACK_CPLXF_TYPE* z, int ldz, ARPACK_CPLXF_TYPE sigma, ARPACK_CPLXF_TYPE* workev, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork);
void zneupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, ARPACK_CPLX_TYPE* d, ARPACK_CPLX_TYPE* z, int ldz, ARPACK_CPLX_TYPE sigma, ARPACK_CPLX_TYPE* workev, ARPACK_CPLX_TYPE* resid, ARPACK_CPLX_TYPE* v, int ldv, int* ipntr, ARPACK_CPLX_TYPE* workd, ARPACK_CPLX_TYPE* workl, double* rwork);

void ssaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dsaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

void sseupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, float* d, float* z, int ldz, float sigma, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void dseupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, double* d, double* z, int ldz, double sigma, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

#endif /* ifndef */
