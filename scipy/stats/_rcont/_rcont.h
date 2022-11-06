#ifndef RCONT_H
#define RCONT_H

#include <numpy/random/distributions.h>

void rcont1_init(int *work, int nc, const double *c);

void rcont1(double *matrix, int nr, const double *r, int nc, const double *c,
            double ntot, int *work, bitgen_t *rstate);

void rcont2(double *matrix, int nr, const double *r, int nc, const double *c,
            double ntot, bitgen_t *rstate);

#endif