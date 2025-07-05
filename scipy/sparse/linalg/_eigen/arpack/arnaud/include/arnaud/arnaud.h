#ifndef ARNAUD_H
#define ARNAUD_H

#include <math.h>
#include <complex.h>

#if defined(_MSC_VER)
    // MSVC definitions
    typedef _Dcomplex ARNAUD_CPLX_TYPE;
    typedef _Fcomplex ARNAUD_CPLXF_TYPE;

#else
    // C99 compliant compilers
    typedef double complex ARNAUD_CPLX_TYPE;
    typedef float complex ARNAUD_CPLXF_TYPE;

#endif


enum ARNAUD_which {
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

enum ARNAUD_ido {
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

struct ARNAUD_state_s {
    float tol;               // problem parameter                input parameter
    float getv0_rnorm0;      // getv0 internal compute           internal
    float aitr_betaj;        // naitr internal compute           internal
    float aitr_rnorm1;       // naitr internal compute           internal
    float aitr_wnorm;        // naitr internal compute           internal
    float aup2_rnorm;        // naup2 internal compute           internal
    enum ARNAUD_which which; // naupd flow control               input
    enum ARNAUD_ido ido;     // naupd flow control               input/output
    int info;                // problem outcome,                 input/output
    int bmat;                // problem parameter, boolean       input
    int mode;                // problem parameter,               input
    int n;                   // problem parameter,               input
    int ncv;                 // problem parameter,               input
    int nev;                 // problem parameter,               input
    int shift;               // problem parameter, boolean       input
    int maxiter;             // problem parameter,               input
    int nconv;               // problem outcome,                 output
    int iter;                // problem intermediate,            internal
    int np;                  // problem intermediate,            internal
    int getv0_first;         // getv0 flow control               internal
    int getv0_iter;          // getv0 flow control               internal
    int getv0_itry;          // getv0 flow control               internal
    int getv0_orth;          // getv0 flow control               internal
    int aitr_iter;           // naitr flow control               internal
    int aitr_j;              // naitr flow control               internal
    int aitr_orth1;          // naitr flow control               internal
    int aitr_orth2;          // naitr flow control               internal
    int aitr_restart;        // naitr flow control               internal
    int aitr_step3;          // naitr flow control               internal
    int aitr_step4;          // naitr flow control               internal
    int aitr_ierr;           // naitr flow control               internal
    int aup2_initv;          // naupd2 flow control              internal
    int aup2_iter;           // naupd2 flow control              internal
    int aup2_getv0;          // naupd2 flow control              internal
    int aup2_cnorm;          // naupd2 flow control              internal
    int aup2_kplusp;         // naupd2 flow control              internal
    int aup2_nev0;           // naupd2 internal compute          internal
    int aup2_np0;            // naupd2 internal compute          internal
    int aup2_numcnv;         // naupd2 internal compute          internal
    int aup2_update;         // naupd2 flow control              internal
    int aup2_ushift;         // naupd2 flow control              internal
};


struct ARNAUD_state_d {
    double tol;              // problem parameter                input parameter
    double getv0_rnorm0;     // getv0 internal compute           internal
    double aitr_betaj;       // naitr internal compute           internal
    double aitr_rnorm1;      // naitr internal compute           internal
    double aitr_wnorm;       // naitr internal compute           internal
    double aup2_rnorm;       // naup2 internal compute           internal
    enum ARNAUD_which which; // naupd flow control               input
    enum ARNAUD_ido ido;     // naupd flow control               input/output
    int info;                // problem outcome,                 input/output
    int bmat;                // problem parameter, boolean       input
    int mode;                // problem parameter,               input
    int n;                   // problem parameter,               input
    int ncv;                 // problem parameter,               input
    int nev;                 // problem parameter,               input
    int shift;               // problem parameter, boolean       input
    int maxiter;             // problem parameter,               input
    int nconv;               // problem outcome,                 output
    int iter;                // problem intermediate,            internal
    int np;                  // problem intermediate,            internal
    int getv0_first;         // getv0 flow control               internal
    int getv0_iter;          // getv0 flow control               internal
    int getv0_itry;          // getv0 flow control               internal
    int getv0_orth;          // getv0 flow control               internal
    int aitr_iter;           // naitr flow control               internal
    int aitr_j;              // naitr flow control               internal
    int aitr_orth1;          // naitr flow control               internal
    int aitr_orth2;          // naitr flow control               internal
    int aitr_restart;        // naitr flow control               internal
    int aitr_step3;          // naitr flow control               internal
    int aitr_step4;          // naitr flow control               internal
    int aitr_ierr;           // naitr flow control               internal
    int aup2_initv;          // naupd2 flow control              internal
    int aup2_iter;           // naupd2 flow control              internal
    int aup2_getv0;          // naupd2 flow control              internal
    int aup2_cnorm;          // naupd2 flow control              internal
    int aup2_kplusp;         // naupd2 flow control              internal
    int aup2_nev0;           // naupd2 internal compute          internal
    int aup2_np0;            // naupd2 internal compute          internal
    int aup2_numcnv;         // naupd2 internal compute          internal
    int aup2_update;         // naupd2 flow control              internal
    int aup2_ushift;         // naupd2 flow control              internal
};


void ARNAUD_snaupd(struct ARNAUD_state_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARNAUD_dnaupd(struct ARNAUD_state_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void ARNAUD_cnaupd(struct ARNAUD_state_s *V, ARNAUD_CPLXF_TYPE* resid, ARNAUD_CPLXF_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLXF_TYPE* workd, ARNAUD_CPLXF_TYPE* workl, float* rwork);
void ARNAUD_znaupd(struct ARNAUD_state_d *V, ARNAUD_CPLX_TYPE* resid, ARNAUD_CPLX_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLX_TYPE* workd, ARNAUD_CPLX_TYPE* workl, double* rwork);

void ARNAUD_sneupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARNAUD_dneupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void ARNAUD_cneupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, ARNAUD_CPLXF_TYPE* d, ARNAUD_CPLXF_TYPE* z, int ldz, ARNAUD_CPLXF_TYPE sigma, ARNAUD_CPLXF_TYPE* workev, ARNAUD_CPLXF_TYPE* resid, ARNAUD_CPLXF_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLXF_TYPE* workd, ARNAUD_CPLXF_TYPE* workl, float* rwork);
void ARNAUD_zneupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, ARNAUD_CPLX_TYPE* d, ARNAUD_CPLX_TYPE* z, int ldz, ARNAUD_CPLX_TYPE sigma, ARNAUD_CPLX_TYPE* workev, ARNAUD_CPLX_TYPE* resid, ARNAUD_CPLX_TYPE* v, int ldv, int* ipntr, ARNAUD_CPLX_TYPE* workd, ARNAUD_CPLX_TYPE* workl, double* rwork);

void ARNAUD_ssaupd(struct ARNAUD_state_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARNAUD_dsaupd(struct ARNAUD_state_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

void ARNAUD_sseupd(struct ARNAUD_state_s *V, int rvec, int howmny, int* select, float* d, float* z, int ldz, float sigma, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARNAUD_dseupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select, double* d, double* z, int ldz, double sigma, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

#endif /* ifndef */
