#ifndef SSS_MAT_H
#define SSS_MAT_H

#include "Python.h"

typedef struct {
  PyObject_VAR_HEAD
  int n;			/* array dimension */
  int nnz;			/* number of stored items */
  double *val;			/* pointer to array of values */
  double *diag;			/* pointer to diagonal elements */
  int *col;			/* pointer to array of indices */
  int *ind;			/* pointer to array of indices */
} SSSMatObject;

#ifdef SPMATRIX_MODULE
/* forward declarations */
static PyTypeObject SSSMatType;	
static PyObject *newSSSMatObject(int n, int nnz);
#endif

#endif
