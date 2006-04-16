#ifndef QMRS_H
#define QMRS_H

#include "fortran.h"
#include "Python.h"

int Itsolvers_qmrs_kernel(int n, 
			  double *b, 
			  double *x, 
			  double *work, 
			  double tol, 
			  int maxitera, 
			  int *itera, 
			  double *err,
			  PyObject *mat_obj,
			  PyObject *prec_obj);
#endif
