
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
#include "SuperLU/SRC/zsp_defs.h"
#include "_superluobject.h"

extern jmp_buf _superlu_py_jmpbuf;


static char doc_zgssv[] = "Direct inversion of sparse matrix.\n\nX = zgssv(A,B) solves A*X = B for X.";

static PyObject *Py_zgssv (PyObject *self, PyObject *args, PyObject *kwdict)
{
  PyObject *Py_B=NULL, *Py_X=NULL;
  PyArrayObject *nzvals=NULL;
  PyArrayObject *colind=NULL, *rowptr=NULL;
  int N, nnz;
  int csc=0, permc_spec=2;
  int info;
  int *perm_r=NULL, *perm_c=NULL;
  SuperMatrix A, B, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;
  
  static char *kwlist[] = {"N","nnz","nzvals","colind","rowptr","B", "csc", "permc_spec",NULL};

  /* Get input arguments */
  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiO!O!O!O|ii", kwlist, &N, &nnz, &PyArray_Type, &nzvals, &PyArray_Type, &colind, &PyArray_Type, &rowptr, &Py_B, &csc, &permc_spec))
      return NULL;


  /* Create Space for output */
  Py_X = PyArray_CopyFromObject(Py_B,PyArray_CDOUBLE,1,2);
  if (Py_X == NULL) goto fail;
  if (csc) {
      if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr, PyArray_CDOUBLE)) goto fail;
  }
  else {
      if (NRFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr, PyArray_CDOUBLE)) goto fail; 
  }
  if (DenseSuper_from_Numeric(&B, Py_X)) goto fail;  

  /* B and Py_X  share same data now but Py_X "owns" it */
    
  /* Setup options */
  
  if (setjmp(_superlu_py_jmpbuf)) goto fail;
  else {
      perm_c = intMalloc(N);
      perm_r = intMalloc(N);
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

  return Py_BuildValue("Ni", Py_X, info);

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

/*******************************Begin Code Adapted from PySparse *****************/

static char doc_zgstrf[] = "zgstrf(A, ...)\n\
\n\
performs a factorization of the sparse matrix A=*(N,nnz,nzvals,rowind,colptr) and \n\
returns a factored_lu object.\n\
\n\
see dgstrf for more information.";

static PyObject *
Py_zgstrf(PyObject *self, PyObject *args, PyObject *keywds) {

  /* default value for SuperLU parameters*/
  double diag_pivot_thresh = 1.0;
  double drop_tol = 0.0;
  int relax = 1;
  int panel_size = 10;
  int permc_spec = 2;
  int N, nnz;
  PyArrayObject *rowind, *colptr, *nzvals;
  SuperMatrix A;
  PyObject *result;
  
  static char *kwlist[] = {"N","nnz","nzvals","rowind","colptr","permc_spec","diag_pivot_thresh", "drop_tol", "relax", "panel_size", NULL};

  int res = PyArg_ParseTupleAndKeywords(args, keywds, "iiO!O!O!|iddii", kwlist, 
                                        &N, &nnz,
					&PyArray_Type, &nzvals,
                                        &PyArray_Type, &rowind,
                                        &PyArray_Type, &colptr,
					&permc_spec,
					&diag_pivot_thresh,
					&drop_tol,
					&relax,
					&panel_size);
  if (!res)
    return NULL;

  if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, rowind, colptr, PyArray_CDOUBLE)) goto fail;
 
  result = newSciPyLUObject(&A, diag_pivot_thresh, drop_tol, relax, panel_size,\
                            permc_spec, PyArray_CDOUBLE);
  if (result == NULL) goto fail;

  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */  
  return result;

 fail:
  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */
  return NULL;
}


/*******************************End Code Adapted from PySparse *****************/


static PyMethodDef zSuperLU_Methods[] = {
   {"zgssv", (PyCFunction) Py_zgssv, METH_VARARGS|METH_KEYWORDS, doc_zgssv},  
   {"zgstrf", (PyCFunction) Py_zgstrf, METH_VARARGS|METH_KEYWORDS, doc_zgstrf},
   /* {"zgstrs", (PyCFunction) Py_zgstrs, METH_VARARGS|METH_KEYWORDS, doc_zgstrs},
      {"_zgscon", Py_zgscon, METH_VARARGS, doc_zgscon},
   {"_zgsequ", Py_zgsequ, METH_VARARGS, doc_zgsequ},
   {"_zlaqgs", Py_zlaqgs, METH_VARARGS, doc_zlaqgs},
   {"_zgsrfs", Py_zgsrfs, METH_VARARGS, doc_zgsrfs}, */
  {NULL, NULL}
};


/* This should be imported first */
DL_EXPORT(void)
init_zsuperlu()
{

  Py_InitModule("_zsuperlu", zSuperLU_Methods);

  import_array();

  if (PyErr_Occurred())
    Py_FatalError("can't initialize module zsuperlu");
}




