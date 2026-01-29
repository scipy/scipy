#ifndef _VODE_H
#define _VODE_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "blaslapack_declarations.h"

/**
 * @brief Struct to hold the VODE common block variables.
 *
 * This struct serves as a C representation of the Fortran common blocks used in VODE.
 * Moreover, original Fortran VODE code, type puns doubles and ints in the same common
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
     * Combining double common blocks DVOD01 and DVOD02
     */
    double acnrm, ccmxj, conp, crate, drc, el[13], eta, etamax, etaq, etaqm1, h, hmin, hmxi, hnew, hscal, prl1, rc, rl1,
           tau[13], tq[5], tn, uround, hu;

    /**
     * All integer variables (41 total)
     * Combining integer common blocks DVOD01 and DVOD02
     */
    int icf, init, ipup, jcur, jstart, jsv, kflag, kuth, l, lmax, lyh, lewt, lacor, lsavf, lwm, liwm, locjs, maxord, meth, miter, msbj,
        mxhnil, mxstep, n, newh, newq, nhnil, nq, nqnyh, nqwait, nslj, nslp, nyh, /* DVOD02 */ ncfn, netf, nfe, nje, nlu, nni, nqu, nst;

} vode_common_struct_t;

#define VODE_STATE_DOUBLE_SIZE 51
#define VODE_STATE_INT_SIZE 41

// Helper functions for state serialization
static inline void pack_vode_state(const vode_common_struct_t* S, double* d_state, int* i_state) {
    // All doubles are contiguous at the start of the struct
    memcpy(d_state, &S->acnrm, VODE_STATE_DOUBLE_SIZE * sizeof(double));
    // All ints are contiguous after the doubles
    memcpy(i_state, &S->icf, VODE_STATE_INT_SIZE * sizeof(int));
}

static inline void unpack_vode_state(const double* d_state, const int* i_state, vode_common_struct_t* S) {
    // Restore doubles
    memcpy(&S->acnrm, d_state, VODE_STATE_DOUBLE_SIZE * sizeof(double));
    // Restore ints
    memcpy(&S->icf, i_state, VODE_STATE_INT_SIZE * sizeof(int));
}


typedef void (*vode_func_t)(int neq, double t, double* y, double* ydot, double* rpar, int* ipar);
typedef void (*vode_jac_t)(int neq, double t, double* y, int ml, int mu, double* pd, int nrowpd, double* rpar, int* ipar);

// Main DVODE interface
void
dvode(
    vode_common_struct_t* S,
    vode_func_t f,
    int neq,
    double* restrict y,
    double* t,
    double* tout,
    int itol,
    double* rtol,
    double* atol,
    int* itask,
    int* istate,
    int* iopt,
    double* restrict rwork,
    int lrw,
    int* restrict iwork,
    int liw,
    vode_jac_t jac,
    int mf,
    double* rpar,
    int* ipar
);

#endif
