#ifndef CGS_H
#define CGS_H

#include <Python.h>

int Itsolvers_cgs_kernel(int n, 
			 double *b, 
			 double *x, 
			 int maxit, 
			 double tol,
			 double *work,
			 int *iter,
			 double *res,
			 PyObject *mat_obj,
			 PyObject *prec_obj);

#endif
