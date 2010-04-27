#include <Python.h>

#define NO_IMPORT_ARRAY
#include "_superluobject.h"
#include <setjmp.h>

extern jmp_buf _superlu_py_jmpbuf;

/*********************************************************************** 
 * SciPyLUObject methods
 */

static char solve_doc[] = "x = self.solve(b, trans)\n\
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
         'H': solve A^H * x == b (not yet implemented)\n\
         (optional, default value 'N')\n\
";

static PyObject *
SciPyLU_solve(SciPyLUObject *self, PyObject *args, PyObject *kwds) {
  PyArrayObject *b, *x=NULL;
  SuperMatrix B;
  char itrans = 'N';
  int info;
  trans_t trans;
  SuperLUStat_t stat;

  static char *kwlist[] = {"rhs","trans",NULL};

  if (!CHECK_SLU_TYPE(self->type)) {
      PyErr_SetString(PyExc_ValueError, "unsupported data type");
      return NULL;
  }

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|c", kwlist,
                                   &PyArray_Type, &b, 
                                   &itrans))
    return NULL;

  /* solve transposed system: matrix was passed row-wise instead of
   * column-wise */
  if (itrans == 'n' || itrans == 'N')
      trans = NOTRANS;
  else if (itrans == 't' || itrans == 'T')
      trans = TRANS;
  else if (itrans == 'h' || itrans == 'H')
      trans = CONJ;
  else {
    PyErr_SetString(PyExc_ValueError, "trans must be N, T, or H");
    return NULL;
  }

  if ((x = (PyArrayObject *) \
       PyArray_CopyFromObject((PyObject *)b,self->type,1,2))==NULL) return NULL;

  if (b->dimensions[0] != self->n) goto fail;


  if (setjmp(_superlu_py_jmpbuf)) goto fail; 

  if (DenseSuper_from_Numeric(&B, (PyObject *)x)) goto fail;

  StatInit(&stat);

  /* Solve the system, overwriting vector x. */
  gstrs(self->type,
        trans, &self->L, &self->U, self->perm_c, self->perm_r, &B,
        &stat, &info);

  if (info) { 
      PyErr_SetString(PyExc_SystemError,
                      "gstrs was called with invalid arguments");
      goto fail;
  }
  
  /* free memory */
  Destroy_SuperMatrix_Store(&B);
  StatFree(&stat);
  return (PyObject *)x;

fail:
  Destroy_SuperMatrix_Store(&B);  
  StatFree(&stat);
  Py_XDECREF(x);
  return NULL;
}

/** table of object methods
 */
PyMethodDef SciPyLU_methods[] = {
  {"solve", (PyCFunction)SciPyLU_solve, METH_VARARGS|METH_KEYWORDS, solve_doc},
  {NULL, NULL}			/* sentinel */
};


/*********************************************************************** 
 * SciPySuperLUType methods
 */

static void
SciPyLU_dealloc(SciPyLUObject *self)
{
  SUPERLU_FREE(self->perm_r);
  SUPERLU_FREE(self->perm_c);
  Destroy_SuperNode_Matrix(&self->L);
  Destroy_CompCol_Matrix(&self->U);
  PyObject_Del(self);
}

static PyObject *
SciPyLU_getattr(SciPyLUObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->m, self->n);
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
  return Py_FindMethod(SciPyLU_methods, (PyObject *)self, name);
}


/***********************************************************************
 * SciPySuperLUType structure
 */

PyTypeObject SciPySuperLUType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "factored_lu",
  sizeof(SciPyLUObject),
  0,
  (destructor)SciPyLU_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)SciPyLU_getattr,  /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
};


int DenseSuper_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
  int m, n, ldx, nd;
  PyArrayObject *aX;
  
  if (!PyArray_Check(PyX)) {
    PyErr_SetString(PyExc_TypeError, "dgssv: Second argument is not an array.");
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
  
  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(aX->descr->type_num)) {
          PyErr_SetString(PyExc_ValueError, "unsupported data type");
          return -1;
      }
      Create_Dense_Matrix(aX->descr->type_num, X, m, n,
                          aX->data, ldx, SLU_DN,
                          NPY_TYPECODE_TO_SLU(aX->descr->type_num), SLU_GE);
  }
  return 0;
}

/* Natively handles Compressed Sparse Row and CSC */

int NRFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz,
                           PyArrayObject *nzvals, PyArrayObject *colind,
                           PyArrayObject *rowptr, int typenum)
{
  int err = 0;
    
  err = (nzvals->descr->type_num != typenum);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fourth argument must be a 1-D array at least as big as third argument.");
    return -1;
  }

  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(nzvals->descr->type_num)) {
          PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
          return -1;  
      }
      Create_CompRow_Matrix(nzvals->descr->type_num,
                            A, m, n, nnz, nzvals->data, (int *)colind->data,
                            (int *)rowptr->data, SLU_NR,
                            NPY_TYPECODE_TO_SLU(nzvals->descr->type_num),
                            SLU_GE);
  }

  return 0;
}

int NCFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz,
                           PyArrayObject *nzvals, PyArrayObject *rowind,
                           PyArrayObject *colptr, int typenum)
{
  int err=0;

  err = (nzvals->descr->type_num != typenum);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D array at least as big as fourth argument.");
    return -1;
  }


  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(nzvals->descr->type_num)) {
          PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
          return -1;  
      }
      Create_CompCol_Matrix(nzvals->descr->type_num,
                            A, m, n, nnz, nzvals->data, (int *)rowind->data,
                            (int *)colptr->data, SLU_NC,
                            NPY_TYPECODE_TO_SLU(nzvals->descr->type_num),
                            SLU_GE);
  }

  return 0;
}

colperm_t superlu_module_getpermc(int permc_spec)
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
  return NATURAL; /* compiler complains... */
}

PyObject *
newSciPyLUObject(SuperMatrix *A, double diag_pivot_thresh,
		 int relax, int panel_size, int permc_spec,
                 int intype)
{

   /* A must be in SLU_NC format used by the factorization routine. */
  SciPyLUObject *self;
  SuperMatrix AC;     /* Matrix postmultiplied by Pc */
  int lwork = 0;
  int *etree=NULL;
  int info;
  int n;
  superlu_options_t options;
  SuperLUStat_t stat;
  
  n = A->ncol;

  /* Create SciPyLUObject */
  self = PyObject_New(SciPyLUObject, &SciPySuperLUType);
  if (self == NULL)
    return PyErr_NoMemory();
  self->m = A->nrow;
  self->n = n;
  self->perm_r = NULL;
  self->perm_c = NULL;
  self->type = intype;

  if (setjmp(_superlu_py_jmpbuf)) goto fail;
  
  /* Calculate and apply minimum degree ordering*/
  etree = intMalloc(n);
  self->perm_r = intMalloc(n);
  self->perm_c = intMalloc(n);

  set_default_options(&options);
  options.ColPerm=superlu_module_getpermc(permc_spec);
  options.DiagPivotThresh = diag_pivot_thresh;
  StatInit(&stat);
  
  get_perm_c(permc_spec, A, self->perm_c); /* calc column permutation */
  sp_preorder(&options, A, self->perm_c, etree, &AC); /* apply column
                                                       * permutation */

  /* Perform factorization */
  if (!CHECK_SLU_TYPE(SLU_TYPECODE_TO_NPY(A->Dtype))) {
      PyErr_SetString(PyExc_ValueError, "Invalid type in SuperMatrix.");
      goto fail;
  }
  gstrf(SLU_TYPECODE_TO_NPY(A->Dtype),
        &options, &AC, relax, panel_size,
        etree, NULL, lwork, self->perm_c, self->perm_r,
        &self->L, &self->U, &stat, &info);

  if (info) {
    if (info < 0)
        PyErr_SetString(PyExc_SystemError,
                        "gstrf was called with invalid arguments");
    else {
        if (info <= n) 
            PyErr_SetString(PyExc_RuntimeError, "Factor is exactly singular");
        else
            PyErr_NoMemory();
    }
    goto fail;
  }

  /* free memory */
  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&AC);
  StatFree(&stat);
  
  return (PyObject *)self;

fail:
  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&AC);
  StatFree(&stat);
  SciPyLU_dealloc(self);
  return NULL;
}
