#ifndef __LBFGSB_H
#define __LBFGSB_H

#include "blaslapack_declarations.h"
#include <math.h>

void setulb(
    int n, int m, double* x, double* l, double* u, int* nbd, double* f, double* g, double factr,
    double pgtol, double* wa, int* iwa, int* task, int* lsave, int* isave, double* dsave, int maxls, int* ln_task
);


#endif /* ifndef */
