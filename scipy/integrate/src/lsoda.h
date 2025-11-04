#ifndef __LSODA_H
#define __LSODA_H

#include <math.h>
#include <string.h>
#include "blaslapack_declarations.h"


/**
 * @brief Struct to hold the LSODA common block variables.
 *
 * This struct serves as a C representation of the Fortran common blocks used in LSODA.
 * Moreover, original Fortran LSODA code, type puns doubles and ints in the same common
 * block making it impossible to decipher which variables are used in which way. Hence,
 * those punned variables are replicated as separate variables. While this slightly
 * increases the memory usage, it greatly improves code clarity and maintainability.
 *
 * NOTE: The struct is organized with all doubles first, then all ints, to enable
 * efficient serialization via memcpy for state persistence between Python calls.
 */
typedef struct {
    /*
     * All double precision variables (240 total)
     * Combining double common blocks ls0001 and lsa001
     */
    double conit, crate, el[13], elco[156], hold, rmax, tesco[36], ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround;
    double tsw, pdest, pdlast, ratio, cm1[12], cm2[5], pdnorm;

    /*
     * All integer variables (48 total)
     * Combining integer common blocks ls0001 and lsa001
     */
    int illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil,
        nhnil, ntrep, nslast, nyh, ialth, ipup, lmax, meo, nqnyh, nslp,
        icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord,
        maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu, /* lsa001 part*/ insufr,
        insufi, ixpr, icount, irflag, jtyp, mused, mxordn, mxords;
} lsoda_common_struct_t;

#define LSODA_STATE_DOUBLE_SIZE 240
#define LSODA_STATE_INT_SIZE 48

// Helper functions for state serialization
static inline void pack_lsoda_state(const lsoda_common_struct_t* S, double* d_state, int* i_state) {
    // All doubles are contiguous at the start of the struct
    memcpy(d_state, &S->conit, LSODA_STATE_DOUBLE_SIZE * sizeof(double));
    // All ints are contiguous after the doubles
    memcpy(i_state, &S->illin, LSODA_STATE_INT_SIZE * sizeof(int));
}

static inline void unpack_lsoda_state(const double* d_state, const int* i_state, lsoda_common_struct_t* S) {
    // Restore doubles
    memcpy(&S->conit, d_state, LSODA_STATE_DOUBLE_SIZE * sizeof(double));
    // Restore ints
    memcpy(&S->illin, i_state, LSODA_STATE_INT_SIZE * sizeof(int));
}


typedef void (*lsoda_func_t)(int* neq, double* t, double* y, double* ydot);
typedef void (*lsoda_jac_t)(int* neq, double* t, double* y, int* ml, int* mu, double* pd, int* nrowpd);


void
lsoda(
    lsoda_func_t f,
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
    lsoda_jac_t jac,
    const int jt,
    lsoda_common_struct_t* S
);

#endif
