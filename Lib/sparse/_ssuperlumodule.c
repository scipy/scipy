
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
#include "Numeric/arrayobject.h"
#include "SuperLU/SRC/ssp_defs.h"
#include "SuperLU/SRC/util.h"

extern jmp_buf _superlu_py_jmpbuf;

/* Natively handles Compressed Sparse Row */

static int NRFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz, PyArrayObject *nzvals, PyArrayObject *colind, PyArrayObject *rowptr)
{
  int retval = -1, err=0;

  err = (nzvals->descr->type_num != PyArray_FLOAT);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D float array at least as big as fourth argument.");
    return retval;
  }

  if (setjmp(_superlu_py_jmpbuf)) return retval;
  else sCreate_CompRow_Matrix(A, m, n, nnz, (float *)nzvals->data, (int *)colind->data, (int *)rowptr->data, SLU_NR, SLU_S, SLU_GE);
  retval = 0;
  return retval;
}

static int NCFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz, PyArrayObject *nzvals, PyArrayObject *rowind, PyArrayObject *colptr)
{
  int retval = -1, err=0;

  err = (nzvals->descr->type_num != PyArray_FLOAT);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D float array at least as big as fourth argument.");
    return retval;
  }

  if (setjmp(_superlu_py_jmpbuf)) return retval;
  else sCreate_CompCol_Matrix(A, m, n, nnz, (float *)nzvals->data, (int *)rowind->data, (int *)colptr->data, SLU_NC, SLU_S, SLU_GE);
  retval = 0;
  return retval;
}


static int Dense_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
  int m, n, ldx, nd;
  PyArrayObject *aX;
 
  if (!PyArray_Check(PyX)) {
    PyErr_SetString(PyExc_TypeError, "sgssv: Second argument is not an array.");
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

  if (setjmp(_superlu_py_jmpbuf)) return -1;
  else sCreate_Dense_Matrix(X, m, n, (float *)aX->data, ldx, SLU_DN, SLU_S, SLU_GE);

  return 0;
}

static colperm_t superlu_module_getpermc(int permc_spec)
{
  switch(permc_spec) {
  case 0:
    return NATURAL;
  case 1:
    return MMD_ATA;
  case 2:
    return MMD_AT_PLUS_A;
  case 3:
    return COLAMD;
  }
  ABORT("Invalid input for permc_spec.");
}

static char doc_sgssv[] = "Direct inversion of sparse matrix.\n\nX = sgssv(A,B) solves A*X = B for X.";

static PyObject *Py_sgssv (PyObject *self, PyObject *args, PyObject *kwdict)
{
  PyObject *Py_B=NULL, *Py_X=NULL;
  PyArrayObject *nzvals=NULL;
  PyArrayObject *colind=NULL, *rowptr=NULL;
  int M, N, nnz;
  int info, full_output=0;
  int csc=0, permc_spec=0;
  int *perm_r=NULL, *perm_c=NULL;
  SuperMatrix A, B, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;

  static char *kwlist[] = {"M","N","nnz","nzvals","colind","rowptr","B", "csc", "permc_spec", "full_output",NULL};

  /* Get input arguments */
  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiiO!O!O!O|iii", kwlist, &M, &N, &nnz, &PyArray_Type, &nzvals, &PyArray_Type, &colind, &PyArray_Type, &rowptr, &Py_B, &csc, &permc_spec, &full_output))
    return NULL;

  /* Create Space for output */
  Py_X = PyArray_CopyFromObject(Py_B,PyArray_FLOAT,1,2);
  if (Py_X == NULL) goto fail;
  if (csc) {
      if (NCFormat_from_spMatrix(&A, M, N, nnz, nzvals, colind, rowptr)) goto fail;
  }
  else {
      if (NRFormat_from_spMatrix(&A, M, N, nnz, nzvals, colind, rowptr)) goto fail; 
  }

  if (Dense_from_Numeric(&B, Py_X)) goto fail;

  /* B and Py_X  share same data now but Py_X "owns" it */
    
  /* Setup options */
  
  if (setjmp(_superlu_py_jmpbuf)) goto fail;
  else {
      perm_c = intMalloc(N);
      perm_r = intMalloc(M);
      set_default_options(&options);
      options.ColPerm=superlu_module_getpermc(permc_spec);
      StatInit(&stat);

  /* Compute direct inverse of sparse Matrix */
      sgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
  }

  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);

  if (full_output)
      return Py_BuildValue("Ni", Py_X, info);
  else
      return Py_X;

 fail:
  SUPERLU_FREE(perm_r);
  SUPERLU_FREE(perm_c);
  Destroy_SuperMatrix_Store(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  StatFree(&stat);

  Py_XDECREF(Py_X);
  return NULL;
}

   
static PyMethodDef sSuperLU_Methods[] = {
   {"sgssv", (PyCFunction) Py_sgssv, METH_VARARGS|METH_KEYWORDS, doc_sgssv},  
   /*   {"_sgstrf", Py_sgstrf, METH_VARARGS, doc_sgstrf},
   {"_sgstrs", Py_sgstrs, METH_VARARGS, doc_sgstrs},
   {"_sgscon", Py_sgscon, METH_VARARGS, doc_sgscon},
   {"_sgsequ", Py_sgsequ, METH_VARARGS, doc_sgsequ},
   {"_slaqgs", Py_slaqgs, METH_VARARGS, doc_slaqgs},
   {"_sgsrfs", Py_sgsrfs, METH_VARARGS, doc_sgsrfs}, */
  {NULL, NULL}
};


void init_ssuperlu()
{
  Py_InitModule("_ssuperlu", sSuperLU_Methods);
  import_array();

}




