
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
#include "SuperLU/SRC/ssp_defs.h"
#include "_superluobject.h"

extern jmp_buf _superlu_py_jmpbuf;


static char doc_sgssv[] = "Direct inversion of sparse matrix.\n\nX = sgssv(A,B) solves A*X = B for X.";

static PyObject *Py_sgssv (PyObject *self, PyObject *args, PyObject *kwdict)
{
  PyObject *Py_B=NULL, *Py_X=NULL;
  PyArrayObject *nzvals=NULL;
  PyArrayObject *colind=NULL, *rowptr=NULL;
  int N, nnz;
  int info;
  int csc=0, permc_spec=2;
  int *perm_r=NULL, *perm_c=NULL;
  SuperMatrix A, B, L, U;
  superlu_options_t options;
  SuperLUStat_t stat;

  static char *kwlist[] = {"N","nnz","nzvals","colind","rowptr","B", "csc", "permc_spec",NULL};

  /* Get input arguments */
  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiO!O!O!O|ii", kwlist, &N, &nnz, &PyArray_Type, &nzvals, &PyArray_Type, &colind, &PyArray_Type, &rowptr, &Py_B, &csc, &permc_spec))
    return NULL;

  /* Create Space for output */
  Py_X = PyArray_CopyFromObject(Py_B,PyArray_FLOAT,1,2);
  if (Py_X == NULL) goto fail;
  if (csc) {
      if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr, PyArray_FLOAT)) goto fail;
  }
  else {
      if (NRFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr, PyArray_FLOAT)) goto fail; 
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
      sgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
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


static char doc_sgstrf[] = "sgstrf(A, ...)\n\
\n\
performs a factorization of the sparse matrix A=*(N,nnz,nzvals,rowind,colptr) and \n\
returns a factored_lu object.\n\
\n\
see dgstrf for more information.";

static PyObject *
Py_sgstrf(PyObject *self, PyObject *args, PyObject *keywds) {

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

  if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, rowind, colptr, PyArray_FLOAT)) goto fail;
 
  result = newSciPyLUObject(&A, diag_pivot_thresh, drop_tol, relax, panel_size,\
                            permc_spec, PyArray_FLOAT);
  if (result == NULL) goto fail;

  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */  
  return result;

 fail:
  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */
  return NULL;
}


/*******************************End Code Adapted from PySparse *****************/

   
static PyMethodDef sSuperLU_Methods[] = {
   {"sgssv", (PyCFunction) Py_sgssv, METH_VARARGS|METH_KEYWORDS, doc_sgssv},  
   {"sgstrf", (PyCFunction) Py_sgstrf, METH_VARARGS|METH_KEYWORDS, doc_sgstrf},
   /*   {"_sgstrs", Py_sgstrs, METH_VARARGS, doc_sgstrs},
   {"_sgscon", Py_sgscon, METH_VARARGS, doc_sgscon},
   {"_sgsequ", Py_sgsequ, METH_VARARGS, doc_sgsequ},
   {"_slaqgs", Py_slaqgs, METH_VARARGS, doc_slaqgs},
   {"_sgsrfs", Py_sgsrfs, METH_VARARGS, doc_sgsrfs}, */
  {NULL, NULL}
};

DL_EXPORT(void)
init_ssuperlu()
{
  Py_InitModule("_ssuperlu", sSuperLU_Methods);
  import_array();

}




