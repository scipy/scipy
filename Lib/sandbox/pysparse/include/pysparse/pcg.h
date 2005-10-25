#ifndef PCG_H
#define PCG_H

#include "Python.h"

int Itsolvers_pcg_kernel(int n, 
			 double *x, 
			 double *b,
			 double tol, 
			 int maxit,
			 int clvl,
			 int *iter, 
			 double *relres, 
			 int *flag,
			 double *work,
			 PyObject *mat_obj,
			 PyObject *prec_obj);
#endif
