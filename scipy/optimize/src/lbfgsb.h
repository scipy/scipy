#ifndef __LBFGSB_H
#define __LBFGSB_H

#include "blaslapack_declarations.h"
#include <math.h>

void setulb(
    CBLAS_INT n, CBLAS_INT m, double* x, double* l, double* u, CBLAS_INT* nbd, double* f, double* g, double factr,
    double pgtol, double* wa, CBLAS_INT* iwa, CBLAS_INT* task, CBLAS_INT* lsave, CBLAS_INT* isave, double* dsave, CBLAS_INT maxls, CBLAS_INT* ln_task
);


#endif /* ifndef */
