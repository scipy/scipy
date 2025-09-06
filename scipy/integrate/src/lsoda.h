#ifndef __LSODA_H
#define __LSODA_H

#include <math.h>
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
 */
typedef struct {
    /*
     * Common block ls0001
     */
    double conit, crate, el[13], elco[156], hold, rmax, tesco[36],                                  // rowns[209]
        ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround;
    int illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, // iownd[14]
        ialth, ipup, lmax, meo, nqnyh, nslp,                                                        // iowns[6]
        icf, ierpj, iersl, jcur, jstart, kflag, l,
        meth, miter, maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu;
    /*
     * Common block lsa001 - We potentially disturb the data alignment
     * due to double, int, double, int variable grouping, in order to have
     * consistent mapping of the common blocks ls0001 and lsa001.
     */
    double tsw, pdest, pdlast, ratio, cm1[12], cm2[5], pdnorm;
    int insufr, insufi, ixpr, icount, irflag, jtyp, mused, mxordn, mxords;
} lsoda_common_struct_t;


typedef void (*lsoda_func_t)(int* neq, double* t, double* y, double* ydot);
typedef void (*lsoda_jac_t)(int* neq, double* t, double* y, int* ml, int* mu, double* pd, int* nrowpd);


void lsoda(
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
