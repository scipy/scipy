#ifndef PROPACK_H
#define PROPACK_H

#include "types.h"

void slansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_s aprod,
             float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
             float tolin, float* work, int lwork, int* iwork,
             float* doption, int* ioption, int* info, float* dparm, int* iparm,
             uint64_t* rng_state);

void dlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_d aprod,
             double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
             double tolin, double* work, int lwork, int* iwork,
             double* doption, int* ioption, int* info, double* dparm, int* iparm,
             uint64_t* rng_state);

void clansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_c aprod,
             PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd, PROPACK_CPLXF_TYPE* V, int ldv,
             float tolin, float* work, int lwork, PROPACK_CPLXF_TYPE* cwork, int lcwork,
             int* iwork, float* doption, int* ioption, int* info,
             PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state);

void zlansvd(int jobu, int jobv, int m, int n, int k, int kmax, PROPACK_aprod_z aprod,
             PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd, PROPACK_CPLX_TYPE* V, int ldv,
             double tolin, double* work, int lwork, PROPACK_CPLX_TYPE* zwork, int lzwork,
             int* iwork, double* doption, int* ioption, int* info,
             PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state);


void slansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_s aprod, float* U, int ldu, float* sigma, float* bnd, float* V, int ldv,
                 float tolin, float* work, int lwork, int* iwork, float* doption, int* ioption,
                 int* info, float* dparm, int* iparm, uint64_t* rng_state);

void dlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_d aprod, double* U, int ldu, double* sigma, double* bnd, double* V, int ldv,
                 double tolin, double* work, int lwork, int* iwork, double* doption, int* ioption,
                 int* info, double* dparm, int* iparm, uint64_t* rng_state);

void clansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_c aprod, PROPACK_CPLXF_TYPE* U, int ldu, float* sigma, float* bnd,
                 PROPACK_CPLXF_TYPE* V, int ldv, float tolin, float* work, int lwork,
                 PROPACK_CPLXF_TYPE* cwork, int lcwork, int* iwork, float* soption,
                 int* ioption, int* info, PROPACK_CPLXF_TYPE* cparm, int* iparm, uint64_t* rng_state);

void zlansvd_irl(int which, int jobu, int jobv, int m, int n, int dim, int p, int *neig, int maxiter,
                 PROPACK_aprod_z aprod, PROPACK_CPLX_TYPE* U, int ldu, double* sigma, double* bnd,
                 PROPACK_CPLX_TYPE* V, int ldv, double tolin, double* work, int lwork,
                 PROPACK_CPLX_TYPE* zwork, int lzwork, int* iwork, double* doption,
                 int* ioption, int* info, PROPACK_CPLX_TYPE* zparm, int* iparm, uint64_t* rng_state);

#endif
