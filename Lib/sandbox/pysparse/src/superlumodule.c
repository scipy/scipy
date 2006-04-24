#include <stdio.h>
#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL superlu
#include "numpy/arrayobject.h"
#include "pysparse/dsp_defs.h"
#include "pysparse/util.h"
#include "pysparse/spmatrix.h"

/*********************************************************************** 
 * static variable
 */

static int StatInit_done = 0;	/* flag showing whether StatInit was already called */

/*********************************************************************** 
 * SuperLUObject definition
 */

typedef struct SuperLUObject {
  PyObject_VAR_HEAD
  int n;
  SuperMatrix L;
  SuperMatrix U;
  int *perm_r;
  int *perm_c;
} SuperLUObject;

/*********************************************************************** 
 * SuperLUObject methods
 */

static char solve_doc[] = "self.solve(b, x, trans)\n\
\n\
solves linear system of equations with one or sereral right hand sides.\n\
\n\
parameters\n\
----------\n\
\n\
b        array, right hand side(s) of equation\n\
x        array, solution vector(s)\n\
trans    'N': solve A   * x == b\n\
         'T': solve A^T * x == b\n\
         (optional, default value 'N')\n\
";

static PyObject *
SuperLU_solve(SuperLUObject *self, PyObject *args) {
  PyArrayObject *b, *x;
  SuperMatrix B;
  char trans = 'N';
  int i, info;

  if (!PyArg_ParseTuple(args, "O!O!|c", 
			&PyArray_Type, &b, 
			&PyArray_Type, &x,
			&trans))
    return NULL;
  
  SPMATRIX_CHECK_ARR_DIM_SIZE(b, 1, self->n);
  SPMATRIX_CHECK_ARR_DIM_SIZE(x, 1, self->n);

  /* solve transposed system: matrix was passed row-wise instead of column-wise */
  if (trans == 'n' || trans == 'N')
    trans = 'T';
  else if (trans == 't' || trans == 'T')
    trans = 'N';
  else {
    PyErr_SetString(PyExc_ValueError, "trans");
    return NULL;
  }

  /* copy b to x */
  for (i = 0; i < self->n; i ++)
    ((double *)x->data)[i] = ((double *)b->data)[i];

  /* Create data structure for right hand side */
#ifdef OLD_SUPERLU
  dCreate_Dense_Matrix(&B, self->n, 1, (double *)x->data, self->n, DN, _D, GE);
#else
  dCreate_Dense_Matrix(&B, self->n, 1, (double *)x->data, self->n, DN, D_, GE);
#endif

  /* Solve the system, overwriting vector x. */
  dgstrs(&trans, &self->L, &self->U, self->perm_r, self->perm_c, &B, &info);

  /* free memory */
  Destroy_SuperMatrix_Store(&B);

  if (info) {
    PyErr_SetString(PyExc_SystemError, "dgstrs was called with invalid arguments");
    return NULL;
  } else {
    Py_INCREF(Py_None); 
    return Py_None;
  }
}

static char precon_doc[] = "";

static PyObject *
SuperLU_precon(LLMatObject *self, PyObject *args) {}

/** table of object methods
 */
PyMethodDef SuperLU_methods[] = {
  {"solve", (PyCFunction)SuperLU_solve, METH_VARARGS, solve_doc},
  {NULL, NULL}			/* sentinel */
};


/*********************************************************************** 
 * SuperLUType methods
 */

static void
SuperLU_dealloc(SuperLUObject *self)
{
  SUPERLU_FREE(self->perm_r);
  SUPERLU_FREE(self->perm_c);
  Destroy_SuperNode_Matrix(&self->L);
  Destroy_CompCol_Matrix(&self->U);
  PyObject_Del(self);
}

static PyObject *
SuperLU_getattr(SuperLUObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "nnz") == 0)
    return Py_BuildValue("i", ((SCformat *)self->L.Store)->nnz + ((SCformat *)self->U.Store)->nnz);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
	PyList_SetItem(list, i, PyString_FromString(members[i]));
      if (PyErr_Occurred()) {
	Py_DECREF(list);
	list = NULL;
      }
    }
    return list;
  }
  return Py_FindMethod(SuperLU_methods, (PyObject *)self, name);
}


/***********************************************************************
 * SuperLUType structure
 */

PyTypeObject SuperLUType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "superlu_context",
  sizeof(SuperLUObject),
  0,
  (destructor)SuperLU_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)SuperLU_getattr,  /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
};

/*********************************************************************** 
 * Object construction functions
 */

/*********************************************************************** 
 * Module functions
 */

static PyObject *
newSuperLUObject(int n, CSRMatObject *matrix, 
		 double diag_pivot_thresh, double drop_tol, int relax, int panel_size, int permc_spec)
{
  SuperLUObject *self;
  char refact[1];
  SuperMatrix A;      /* A in NC format used by the factorization routine. */
  SuperMatrix AC;     /* Matrix postmultiplied by Pc */
  mem_usage_t mem_usage;
  int lwork = 0;
  int *etree;
  int info;
  
  *refact = 'N';

  /* make sure StatInit is only called once */
  if (!StatInit_done) {
    StatInit(panel_size, relax);
    StatInit_done = 1;
  }

  /* Create SuperLUObject */
  self = PyObject_New(SuperLUObject, &SuperLUType);
  if (self == NULL)
    return PyErr_NoMemory();
  self->n = n;
  self->perm_r = NULL;
  self->perm_c = NULL;

  /* Create matrix structure */
#ifdef OLD_SUPERLU
  dCreate_CompCol_Matrix(&A, n, n, matrix->nnz, matrix->val, matrix->col, matrix->ind, NC, _D, GE);
#else
  dCreate_CompCol_Matrix(&A, n, n, matrix->nnz, matrix->val, matrix->col, matrix->ind, NC, D_, GE);
#endif

  /* Calculate and apply minimum degree ordering*/
  etree = intMalloc(n);
  self->perm_r = intMalloc(n);
  self->perm_c = intMalloc(n);
  if (self->perm_r == NULL || self->perm_c == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  get_perm_c(permc_spec, &A, self->perm_c); /* calc column permutation */
  sp_preorder(refact, &A, self->perm_c, etree, &AC); /* apply column permutation */

  /* Perform factorization */
  dgstrf(refact, &AC, diag_pivot_thresh, drop_tol, relax, panel_size,
         etree, NULL, lwork, self->perm_r, self->perm_c,
         &self->L, &self->U, &info);

  /* free memory */
  SUPERLU_FREE(etree);
  Destroy_SuperMatrix_Store(&A); /* arrays of input matrix will not be freed */
  Destroy_CompCol_Permuted(&AC);

  if (info) {
    if (info < 0)
      PyErr_SetString(PyExc_SystemError, "dgstrf was called with invalid arguments");
    else {
      if (info <= n) 
	PyErr_SetString(PyExc_RuntimeError, "Factor is exactly singular");
      else
	PyErr_NoMemory();
    }
    goto fail;
  }
  
  return (PyObject *)self;

 fail:
  PyMem_Del(self->perm_r);
  PyMem_Del(self->perm_c);
  PyObject_Del(self);
  return NULL;
}

static char factorize_doc[] = "factorize(A, ...)\n\
\n\
performs a factorization of the sparse matrix A (which is a csr_mat object) and \n\
return a superlu_context object.\n\
\n\
arguments\n\
---------\n\
\n\
A    spmatrix.csr_mat object.\n\
     Matrix to be factorized\n\
\n\
additional keyword arguments:\n\
-----------------------------\n\
\n\
diag_pivot_thresh   threshhold for partial pivoting.\n\
                    0.0 <= diag_pivot_thresh <= 1.0\n\
                    0.0 corresponds to no pivoting\n\
                    1.0 corresponds to partial pivoting\n\
                    (default: 1.0)\n\
\n\
drop_tol            drop tolerance parameter\n\
                    0.0 <= drop_tol <= 1.0\n\
                    0.0 corresponds to exact factorization\n\
                    CAUTION: the drop_tol is not implemented in SuperLU 2.0\n\
                    (default: 0.0)\n\
\n\
relax               to control degree of relaxing supernodes\n\
                    (default: 1)\n\
\n\
panel_size          a panel consist of at most panel_size consecutive columns.\n\
                    (default: 10)\n\
\n\
permc_spec          specifies the matrix ordering used for the factorization\n\
                    0: natural ordering\n\
                    1: MMD applied to the structure of A^T * A\n\
                    2: MMD applied to the structure of A^T + A\n\
                    3: COLAMD, approximate minimum degree column ordering\n\
                    (default: 2)\n\
";

static PyObject *
factorize(PyObject *self, PyObject *args, PyObject *keywds) {
  int n;			/* dimension of matrix */

  /* default value for SuperLU parameters*/
  double diag_pivot_thresh = 1.0;
  double drop_tol = 0.0;
  int relax = 1;
  int panel_size = 10;
  int permc_spec = 2;
  
  CSRMatObject *matrix;
  static char *kwlist[] = {"", "diag_pivot_thresh", "drop_tol", "relax", "panel_size", "permc_spec", NULL};

  int res = PyArg_ParseTupleAndKeywords(args, keywds, "O!|ddiii", kwlist, 
					&CSRMatType, &matrix,
					&diag_pivot_thresh,
					&drop_tol,
					&relax,
					&panel_size,
					&permc_spec);
  if (!res)
    return NULL;

  if (drop_tol != 0.0)
    printf("CAUTION: the drop_tol is not implemented in SuperLU 2.0 and earlier\n");
 
  /* check shape of matrix object */
  if (SpMatrix_GetOrder((PyObject *)matrix, &n))
    return NULL;

  return newSuperLUObject(n, matrix, diag_pivot_thresh, drop_tol, relax, panel_size, permc_spec);
}


/** table of module functions
 */
static PyMethodDef precon_methods[] = {
  {"factorize", (PyCFunction)factorize, METH_VARARGS|METH_KEYWORDS, factorize_doc},
  {NULL, NULL}	/* sentinel */
};


PyMODINIT_FUNC
initsuperlu(void)
{
  PyObject *m, *d;
  
  SuperLUType.ob_type = &PyType_Type;

  m = Py_InitModule("superlu", precon_methods);
  d = PyModule_GetDict(m);

  PyDict_SetItemString(d, "SuperLUType", (PyObject *)&SuperLUType);

  /* initialize scipy array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();

  /* No need to check the error here, the caller will do that */
}
