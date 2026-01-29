#ifndef _ZVODE_H
#define _ZVODE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "blaslapack_declarations.h"


/**
 * @brief Struct to hold the ZVODE common block variables.
 *
 * This struct serves as a C representation of the Fortran common blocks used in ZVODE.
 * Moreover, original Fortran ZVODE code, type puns doubles and ints in the same common
 * block making it impossible to decipher which variables are used in which way. Hence,
 * those punned variables are replicated as separate variables. While this slightly
 * increases the memory usage, it greatly improves code clarity and maintainability.
 *
 * NOTE: The struct is organized with all doubles first, then all ints, to enable
 * efficient serialization via memcpy for state persistence between Python calls.
 */
typedef struct {
    /**
     * All double precision variables (51 total)
     * Combining double common blocks ZVOD01 and ZVOD02
     */
    double acnrm, ccmxj, conp, crate, drc, el[13], eta, etamax, etaq, etaqm1, h, hmin, hmxi, hnew, hrl1, hscal, prl1, rc, rl1,
           srur, tau[13], tq[5], tn, uround, hu;

    /**
     * All integer variables (41 total)
     * Combining integer common blocks ZVOD01 and ZVOD02
     */
    int icf, init, ipup, jcur, jstart, jsv, kflag, kuth, l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm, locjs, maxord, meth, miter, msbj,
        mxhnil, mxstep, n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh, /* DVOD02 */ ncfn, netf, nfe, nje, nlu, nni, nqu, nst;

} zvode_common_struct_t;

#define ZVODE_STATE_DOUBLE_SIZE 53
#define ZVODE_STATE_INT_SIZE 41

// Helper functions for state serialization
static inline void pack_zvode_state(const zvode_common_struct_t* S, double* d_state, int* i_state) {
    // All doubles are contiguous at the start of the struct
    memcpy(d_state, &S->acnrm, ZVODE_STATE_DOUBLE_SIZE * sizeof(double));
    // All ints are contiguous after the doubles
    memcpy(i_state, &S->icf, ZVODE_STATE_INT_SIZE * sizeof(int));
}

static inline void unpack_zvode_state(const double* d_state, const int* i_state, zvode_common_struct_t* S) {
    // Restore doubles
    memcpy(&S->acnrm, d_state, ZVODE_STATE_DOUBLE_SIZE * sizeof(double));
    // Restore ints
    memcpy(&S->icf, i_state, ZVODE_STATE_INT_SIZE * sizeof(int));
}


typedef void (*zvode_func_t)(int neq, double t, ZVODE_CPLX_TYPE* y, ZVODE_CPLX_TYPE* ydot, ZVODE_CPLX_TYPE* rpar, int* ipar);
typedef void (*zvode_jac_t)(int neq, double t, ZVODE_CPLX_TYPE* y, int ml, int mu, ZVODE_CPLX_TYPE* pd, int nrowpd, ZVODE_CPLX_TYPE* rpar, int* ipar);


#if defined(_MSC_VER)
    // Complex arithmetic helpers for MSVC (no _Cadd, _Csub, _Cdiv in MSVC)
    static inline ZVODE_CPLX_TYPE CMPLX_ADD(ZVODE_CPLX_TYPE a, ZVODE_CPLX_TYPE b) {
        return ZVODE_cplx(
            a._Val[0] + b._Val[0],
            a._Val[1] + b._Val[1]
        );
    }

    static inline ZVODE_CPLX_TYPE CMPLX_SUB(ZVODE_CPLX_TYPE a, ZVODE_CPLX_TYPE b) {
        return ZVODE_cplx(
            a._Val[0] - b._Val[0],
            a._Val[1] - b._Val[1]
        );
    }

    static inline ZVODE_CPLX_TYPE CMPLX_DIV(ZVODE_CPLX_TYPE a, ZVODE_CPLX_TYPE b) {
        double denom = b._Val[0] * b._Val[0] + b._Val[1] * b._Val[1];
        return ZVODE_cplx(
            (a._Val[0] * b._Val[0] + a._Val[1] * b._Val[1]) / denom,
            (a._Val[1] * b._Val[0] - a._Val[0] * b._Val[1]) / denom
        );
    }
#endif

// Main ZVODE interface
void
zvode(
    zvode_common_struct_t* S,
    zvode_func_t f,
    int neq,
    ZVODE_CPLX_TYPE* restrict y,
    double* t,
    double* tout,
    int itol,
    double* rtol,
    double* atol,
    int* itask,
    int* istate,
    int* iopt,
    ZVODE_CPLX_TYPE* restrict zwork,
    int lzw,
    double* restrict rwork,
    int lrw,
    int* restrict iwork,
    int liw,
    zvode_jac_t jac,
    int mf,
    ZVODE_CPLX_TYPE* rpar,
    int* ipar
);

#endif
