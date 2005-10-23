#ifndef LL_MAT_H
#define LL_MAT_H

#include "Python.h"

typedef struct {
  PyObject_VAR_HEAD
  int dim[2];			/* array dimension */
  int issym;			/* non-zero, if obj represents a symmetric matrix */
  int nnz;			/* number of stored items */
  int nalloc;			/* allocated size of value and index arrays */
  int free;			/* index to first element in free chain */
  double *val;			/* pointer to array of values */
  int *col;			/* pointer to array of indices */
  int *link;			/* pointer to array of indices */
  int *root;			/* pointer to array of indices */
} LLMatObject;

/******************************************************************************
 *                                                                            *
 * llColIndexlinked -- list data structure which links the entries of         *
 *                     each column of a llmat matrix                          *
 *                                                                            *
 ******************************************************************************/

struct llColIndex{
  int *root;			/* ptr to array storing first element of each column */
  int *row;			/* ptr to array of row indices */
  int *link;			/* ptr to array storing index of next element in column */
  int nzLo;			/* number of non-zero entries in lower triangle */
  int nzDiag;			/* number of non-zero entries on diagonal */
  int nzUp;			/* number of non-zero entries in upper triangle */
};


#ifdef SPMATRIX_MODULE
/* forward declarations */
static PyTypeObject LLMatType;	/* forward declaration */
#endif

#endif
