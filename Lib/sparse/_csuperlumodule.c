
/* Copyright 1999 Travis Oliphant
   Permision to copy and modified this file is granted under the LGPL.
   No warranty is expressed or IMPLIED
*/

/* 
   This file implements glue between the SuperLU library for 
   sparse matrix inversion and Python.
*/


/* We want a low-level interface to:
   xGSSV

   These will be done in separate files due to the include structure of
   SuperLU.

   Define a user abort and a user malloc and free (to keep pointers 
     that will be released on errors)
*/

#include <setjmp.h>
#include "Python.h"
#include "arrayobject.h"
#include "csp_defs.h"
#include "util.h"

extern jmp_buf _superlu_py_jmpbuf;
extern PyObject *_superlumodule_memory_dict;


/* Natively handles Compressed Sparse Row */

static int NRFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz, PyArrayObject *nzvals, PyArrayObject *colind, PyArrayObject *rowptr)
{
  int retval = -1, err=0;

  err = (nzvals->descr->type_num != PyArray_CFLOAT);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fourth argument must be a 1-D complex array at least as big as third argument.");
    return retval;
  }

  cCreate_CompRow_Matrix(A, m, n, nnz, (complex *)nzvals->data, (int *)colind->data, (int *)rowptr->data, NR, _C, GE);
  retval = 0;
  return retval;
}

static int Dense_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
  int m, n, ldx, nd;
  PyArrayObject *aX;
 
  if (!PyArray_Check(PyX)) {
    PyErr_SetString(PyExc_TypeError, "cgssv: Second argument is not an array.");
    return -1;
  }
  aX = (PyArrayObject *)PyX;

  nd = aX->nd;
  if (nd == 1) {
    m = aX->dimensions[0];
    n = 1;
    ldx = m;
  }
  else {  /* nd == 2 */
    m = aX->dimensions[1];
    n = aX->dimensions[0];
    ldx = m;
  }

  cCreate_Dense_Matrix(X, m, n, (complex *)aX->data, ldx, DN, _C, GE);

  return 0;
}

static PyObject *Numeric_from_Dense(SuperMatrix *B)
{
  PyArrayObject *aB;
  DNformat *Bstore;
  int nd, dims[2];

  nd = 1 + (B->ncol > 1);
  
  if (nd == 1) {
    dims[0] = B->nrow;
  }
  else {  /* nd == 2 */
    dims[0] = B->ncol;
    dims[1] = B->nrow;
  }
  Bstore = (DNformat *) B->Store;
  aB = (PyArrayObject *)PyArray_FromDimsAndData(nd, dims, PyArray_CFLOAT, (char *)Bstore->nzval);
  
  return (PyObject *)aB;
}

static char doc_cgssv[] = "Direct inversion of sparse matrix.\n\nX = cgssv(A,B) solves A*X = B for X.";

static PyObject *Py_cgssv (PyObject *self, PyObject *args, PyObject *kwdict)
{
  PyObject *Py_B=NULL, *perm_c=NULL, *Py_X=NULL;
  PyArrayObject *aperm_c=NULL, *aperm_r=NULL, *nzvals=NULL;
  PyArrayObject *colind=NULL, *rowptr=NULL;
  int M, N, nnz;
  int info, dims[1], full_output=0;
  SuperMatrix A, B, L, U;
  
  static char *kwlist[] = {"M","N","nnz","nzvals","colind","rowptr","B", "perm_c","full_output"};

  /* Get input arguments */
  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiiO!O!O!O|Oi", kwlist, &M, &N, &nnz, &PyArray_Type, &nzvals, &PyArray_Type, &colind, &PyArray_Type, &rowptr, &Py_B, &perm_c, &full_output))
    return NULL;

  /* Create Space for output */
  Py_X = PyArray_CopyFromObject(Py_B,PyArray_CFLOAT,1,2);
  if (Py_X == NULL) goto fail;
  if (NRFormat_from_spMatrix(&A, M, N, nnz, nzvals, colind, rowptr)) goto fail; 
  if (Dense_from_Numeric(&B, Py_X)) goto fail;

  if (perm_c == NULL) {
    dims[0] = A.nrow;
    aperm_c = (PyArrayObject *)PyArray_FromDims(1,dims,PyArray_INT);
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     */                                                                  
    if (aperm_c == NULL) goto fail;
    get_perm_c(0, &A, (int *)aperm_c->data);
  }
  else {
    aperm_c = (PyArrayObject *)PyArray_ContiguousFromObject(perm_c,PyArray_INT,1,1);
    if (aperm_c == NULL) goto fail;
    dims[0] = aperm_c -> dimensions[0];
    if (dims[0] != A.nrow) {
      char msg[256];
      sprintf(msg, "Permuation vector should be of size %d", A.nrow);
      PyErr_SetString(PyExc_ValueError, msg);
      goto fail;
    }
  }
  dims[0] = A.ncol;
  aperm_r = (PyArrayObject *)PyArray_FromDims(1,dims,PyArray_INT);
  if (aperm_r == NULL) goto fail;  

  /* Compute direct inverse of sparse Matrix */
  if (setjmp(_superlu_py_jmpbuf)) goto fail;
  else cgssv(&A, (int *)aperm_c->data, (int *)aperm_r->data, &L, &U, &B, &info);

  Py_DECREF(aperm_c);
  Py_DECREF(aperm_r);
  /* Output results (jsut solution for now)*/
  if (full_output)
    Py_X = Py_BuildValue("Ni",Numeric_from_Dense(&B),info);
  else
    Py_X = Numeric_from_Dense(&B);

  /* Delete all memory allocated in calls to SuperLU routines */
  superlu_delete_allkeys();

  return Py_X;
    
 fail:
  superlu_delete_allkeys();
  Py_XDECREF(Py_X);
  Py_XDECREF(aperm_c);
  Py_XDECREF(aperm_r);
  return NULL;
}

   
static PyMethodDef cSuperLU_Methods[] = {
   {"cgssv", (PyCFunction) Py_cgssv, METH_VARARGS|METH_KEYWORDS, doc_cgssv},  
   /*   {"_cgstrf", Py_cgstrf, METH_VARARGS, doc_cgstrf},
   {"_cgstrs", Py_cgstrs, METH_VARARGS, doc_cgstrs},
   {"_cgscon", Py_cgscon, METH_VARARGS, doc_cgscon},
   {"_cgsequ", Py_cgsequ, METH_VARARGS, doc_cgsequ},
   {"_claqgs", Py_claqgs, METH_VARARGS, doc_claqgs},
   {"_cgsrfs", Py_cgsrfs, METH_VARARGS, doc_cgsrfs}, */
  {NULL, NULL}
};


void init_csuperlu()
{
  Py_InitModule("_csuperlu", cSuperLU_Methods);
  import_array();

}




