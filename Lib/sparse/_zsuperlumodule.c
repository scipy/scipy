
/* Copyright 1999 Travis Oliphant
   Permision to copy and modified this file is granted under the LGPL.
   No warranty is expressed or IMPLIED
 
   Changes:  2004 converted to SuperLU_3.0 and added factor and solve routines for
             more flexible handling. 

             Also added NC (compressed sparse column handling -- best to use CSC)
*/

/* 
   This file implements glue between the SuperLU library for 
   sparse matrix inversion and Python.
*/


/* We want a low-level interface to:
   xGSSV
   xgstrf  -- factor
   xgstrs  -- solve

   These will be done in separate files due to the include structure of
   SuperLU.

   Define a user abort and a user malloc and free (to keep pointers 
     that will be released on errors)
*/

#include <setjmp.h>
#include "Python.h"
#include "Numeric/arrayobject.h"
#include "SuperLU/SRC/zsp_defs.h"
#include "SuperLU/SRC/util.h"

extern jmp_buf _superlu_py_jmpbuf;

/* Natively handles Compressed Sparse Row or Compressed Sparse Column */

static int NRFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz, PyArrayObject *nzvals, PyArrayObject *colind, PyArrayObject *rowptr)
{
  int retval = -1, err=0;

  err = (nzvals->descr->type_num != PyArray_CDOUBLE);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D complex double array at least as big as fourth argument.");
    return retval;
  }

  if (setjmp(_superlu_py_jmpbuf)) return retval;
  else zCreate_CompRow_Matrix(A, m, n, nnz, (doublecomplex *)nzvals->data, (int *)colind->data, (int *)rowptr->data, SLU_NR, SLU_Z, SLU_GE);
  retval = 0;
  return retval;
}

static int NCFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz, PyArrayObject *nzvals, PyArrayObject *rowind, PyArrayObject *colptr)
{
  int retval = -1, err=0;

  err = (nzvals->descr->type_num != PyArray_CDOUBLE);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D complex double array at least as big as fourth argument.");
    return retval;
  }

  if (setjmp(_superlu_py_jmpbuf)) return retval;
  else zCreate_CompCol_Matrix(A, m, n, nnz, (doublecomplex *)nzvals->data, (int *)rowind->data, (int *)colptr->data, SLU_NC, SLU_Z, SLU_GE);
  retval = 0;
  return retval;
}


static int Dense_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
  int m, n, ldx, nd;
  PyArrayObject *aX;
 
  if (!PyArray_Check(PyX)) {
    PyErr_SetString(PyExc_TypeError, "zgssv: Second argument is not an array.");
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
  else zCreate_Dense_Matrix(X, m, n, (doublecomplex *)aX->data, ldx, SLU_DN, SLU_Z, SLU_GE);

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


static char doc_zgssv[] = "Direct inversion of sparse matrix.\n\nX = zgssv(A,B) solves A*X = B for X.";

static PyObject *Py_zgssv (PyObject *self, PyObject *args, PyObject *kwdict)
{
  PyObject *Py_B=NULL, *Py_X=NULL;
  PyArrayObject *nzvals=NULL;
  PyArrayObject *colind=NULL, *rowptr=NULL;
  int M, N, nnz;
  int csc=0, permc_spec=0;
  int info, full_output=0;
  int *perm_r=NULL, *perm_c=NULL;
  SuperMatrix A, B, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;
  
  static char *kwlist[] = {"M","N","nnz","nzvals","colind","rowptr","B", "csc", "permc_spec", "full_output",NULL};

  /* Get input arguments */
  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiiO!O!O!O|iii", kwlist, &M, &N, &nnz, &PyArray_Type, &nzvals, &PyArray_Type, &colind, &PyArray_Type, &rowptr, &Py_B, &csc, &permc_spec, &full_output))
    return NULL;


  /* Create Space for output */
  Py_X = PyArray_CopyFromObject(Py_B,PyArray_CDOUBLE,1,2);
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
      zgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
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



static PyMethodDef zSuperLU_Methods[] = {
   {"zgssv", (PyCFunction) Py_zgssv, METH_VARARGS|METH_KEYWORDS, doc_zgssv},  
   /* {"zgstrf", (PyCFunction) Py_zgstrf, METH_VARARGS|METH_KEYWORDS, doc_zgstrf},*/
   /* {"zgstrs", (PyCFunction) Py_zgstrs, METH_VARARGS|METH_KEYWORDS, doc_zgstrs},
      {"_zgscon", Py_zgscon, METH_VARARGS, doc_zgscon},
   {"_zgsequ", Py_zgsequ, METH_VARARGS, doc_zgsequ},
   {"_zlaqgs", Py_zlaqgs, METH_VARARGS, doc_zlaqgs},
   {"_zgsrfs", Py_zgsrfs, METH_VARARGS, doc_zgsrfs}, */
  {NULL, NULL}
};


/* This should be imported first */
void init_zsuperlu()
{

  PyObject *m=NULL, *d=NULL;

  /* SuperLUType.ob_type = &PyType_Type;*/

  m = Py_InitModule("_zsuperlu", zSuperLU_Methods);
  d = PyModule_GetDict(m);

  /* PyDict_SetItemString(d, "ZSuperLUFactoredType", (PyObject *)&ZSuperLUType);*/

  import_array();

  if (PyErr_Occurred())
    Py_FatalError("can't initialize module zsuperlu");
}




