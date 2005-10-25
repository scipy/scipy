#ifndef CSR_MAT_H
#define CSR_MAT_H

#include "Python.h"

typedef struct {
  PyObject_VAR_HEAD
  int dim[2];			/* array dimension */
  int nnz;			/* number of stored items */
  double *val;			/* pointer to array of values */
  int *col;			/* pointer to array of indices */
  int *ind;			/* pointer to array of indices */
} CSRMatObject;

#ifdef SPMATRIX_MODULE
/* forward declarations */
static PyTypeObject CSRMatType;	
static PyObject *newCSRMatObject(int dim[], int nnz);
#endif

#endif
