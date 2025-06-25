/*

 This is a C adaptation and rewrite of the well-known Fortran77 ARPACK large-scale
 eigenvalue problem solver, which is widely used in scientific computing, authored
 by Richard Lehoucq, Kristi Maschhoff, Danny Sorensen, and Chao Yang.

 The source is based on the original Fortran77 and a few of the patches collected
 over the years. The patched Fortran code can be found at arpack-ng repository
 on GitHub, at the time of writing version 3.9.1:

 https://github.com/opencollab/arpack-ng/

 While the translation is done mostly, in a straightforward fashion, however,
 still there are significant changes, namely, XYapps.f and Xstqrb.f are rewritten
 to avoid the goto-based flow. This version also includes API breaking changes to
 make it more flexible to be included in other projects.

 ARPACK uses the so-called reverse-communication style that typically exits
 the program with its in/out arguments to signal, in what stage the algorithm
 is and what it needs. Then user modifies the arguments and calls again with
 the necessary information. Thus the state of the whole program is sent back
 and forth through in-place modified arguments. On top of this, ARPACK also
 uses lots of variables through the Fortran's dreadful use of SAVE attribute
 (similar to that of C language STATIC keyword inside a function body) that
 persists the variable values across consecutive calls. Instead we move all
 those variables into the reverse communication layer by a C-struct bridge and
 for array arguments pointers that are provided by the user to make modifications
 in-place without any alloc/free. This struct bridge also allows for reentrancy
 and avoids the issues that come with thread safety.

 Compared to the original Fortran code, random number generation is now delegated
 to the user side to allow for seed control, custom generators and replicable runs.
 Hence, the ido_RANDOM and ido_RANDOM_OPX codes are used to signal that the user
 input is needed. In turn the ido mode -1 is removed.


 ==============================================================================

 Author: Ilhan Polat
 Copyright (C) 2025 SciPy developers

  Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
  a. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
 b. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution.
 c. Names of the SciPy Developers may not be used to endorse or promote
    products derived from this software without specific prior written
    permission.
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 THE POSSIBILITY OF SUCH DAMAGE.


 Original Fortran77 ARPACK code license;

-------------

 The ARPACK license is the BSD 3-clause license ("New BSD License")

 BSD Software License

 Pertains to ARPACK and P_ARPACK

 Copyright (c) 1996-2008 Rice University.
 Developed by D.C. Sorensen, R.B. Lehoucq, C. Yang, and K. Maschhoff.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer listed
   in this license in the documentation and/or other materials
   provided with the distribution.

 - Neither the name of the copyright holders nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/



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
    float tol;               // problem parameter                input parameter
    float getv0_rnorm0;      // getv0 internal compute           internal
    float aitr_betaj;        // naitr internal compute           internal
    float aitr_rnorm1;       // naitr internal compute           internal
    float aitr_wnorm;        // naitr internal compute           internal
    float aup2_rnorm;        // naup2 internal compute           internal
    enum ARPACK_which which; // naupd flow control               input
    enum ARPACK_ido ido;     // naupd flow control               input/output
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


struct ARPACK_arnoldi_update_vars_d {
    double tol;              // problem parameter                input parameter
    double getv0_rnorm0;     // getv0 internal compute           internal
    double aitr_betaj;       // naitr internal compute           internal
    double aitr_rnorm1;      // naitr internal compute           internal
    double aitr_wnorm;       // naitr internal compute           internal
    double aup2_rnorm;       // naup2 internal compute           internal
    enum ARPACK_which which; // naupd flow control               input
    enum ARPACK_ido ido;     // naupd flow control               input/output
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


void ARPACK_snaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARPACK_dnaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void ARPACK_cnaupd(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork);
void ARPACK_znaupd(struct ARPACK_arnoldi_update_vars_d *V, ARPACK_CPLX_TYPE* resid, ARPACK_CPLX_TYPE* v, int ldv, int* ipntr, ARPACK_CPLX_TYPE* workd, ARPACK_CPLX_TYPE* workl, double* rwork);

void ARPACK_sneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, float* dr, float* di, float* z, int ldz, float sigmar, float sigmai, float* workev, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARPACK_dneupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, double* dr, double* di, double* z, int ldz, double sigmar, double sigmai, double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);
void ARPACK_cneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, ARPACK_CPLXF_TYPE* d, ARPACK_CPLXF_TYPE* z, int ldz, ARPACK_CPLXF_TYPE sigma, ARPACK_CPLXF_TYPE* workev, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork);
void ARPACK_zneupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, ARPACK_CPLX_TYPE* d, ARPACK_CPLX_TYPE* z, int ldz, ARPACK_CPLX_TYPE sigma, ARPACK_CPLX_TYPE* workev, ARPACK_CPLX_TYPE* resid, ARPACK_CPLX_TYPE* v, int ldv, int* ipntr, ARPACK_CPLX_TYPE* workd, ARPACK_CPLX_TYPE* workl, double* rwork);

void ARPACK_ssaupd(struct ARPACK_arnoldi_update_vars_s *V, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARPACK_dsaupd(struct ARPACK_arnoldi_update_vars_d *V, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

void ARPACK_sseupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select, float* d, float* z, int ldz, float sigma, float* resid, float* v, int ldv, int* ipntr, float* workd, float* workl);
void ARPACK_dseupd(struct ARPACK_arnoldi_update_vars_d *V, int rvec, int howmny, int* select, double* d, double* z, int ldz, double sigma, double* resid, double* v, int ldv, int* ipntr, double* workd, double* workl);

#endif /* ifndef */
