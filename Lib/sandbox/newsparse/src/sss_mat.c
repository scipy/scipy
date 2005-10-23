#include "Python.h"

#define SPMATRIX_MODULE
#include "newsparse/spmatrix.h"

#define PY_ARRAY_UNIQUE_SYMBOL spmatrix
/* Was: #include "Numeric/arrayobject.h" */
#include "scipy/arrayobject.h"

/** getitem
 *
 */
static PyObject *
getitem(SSSMatObject *self, int i, int j) {
  int k, t;

  if (i == j)
    return PyFloat_FromDouble(self->diag[i]);
  if (i < j) {
    t = i; i = j; j = t;
  }
  for (k = self->ind[i]; k < self->ind[i+1]; k ++) {
    if (self->col[k] == j)
      return PyFloat_FromDouble(self->val[k]);
  }
  return PyFloat_FromDouble(0.0);
}

/*********************************************************************** 
 * SSSMatObject methods
 */

static char SSSMat_matvec_doc[] = "a.matvec(x, y)\n\
\n\
compute the sparse matrix-vector product y := a * x. \n\
a is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numeric arrays of appropriate size.";

void sss_matvec(int n, double *x, double *y,
		double *va, double *da, int *ja, int *ia) {
  double s, v, xi;
  int i, j, k;
 
  for (i = 0; i < n; i ++) {
    xi = x[i];
    s = 0.0;
    for (k = ia[i]; k < ia[i+1]; k ++) {
      j = ja[k];
      v = va[k];
      s += v * x[j];
      y[j] += v * xi;
    }
    y[i] = s + da[i]*xi;
  }
}

void sss_matvec_stride(int n, 
		       double *x, int incx, 
		       double *y, int incy,
		       double *va, double *da, int *ja, int *ia) {
  double s, v, xi;
  int i, j, k;
 
  for (i = 0; i < n; i ++) {
    xi = x[i * incx];
    s = 0.0;
    for (k = ia[i]; k < ia[i+1]; k ++) {
      j = ja[k];
      v = va[k];
      s += v * x[j * incx];
      y[j * incy] += v * xi;
    }
    y[i * incy] = s + da[i]*xi;
  }
}

static PyObject *
SSSMat_matvec(SSSMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->n, self->n);

  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS) {
    sss_matvec(self->n, (double *)(xp->data), (double *)(yp->data), 
	       self->val, self->diag, self->col, self->ind);
  } else {
    sss_matvec_stride(self->n, 
		      (double *)(xp->data), xp->strides[0] / sizeof(double),
		      (double *)(yp->data), yp->strides[0] / sizeof(double),
		      self->val, self->diag, self->col, self->ind);
  }
  Py_INCREF(Py_None); 
  return Py_None;
}

static char SSSMat_matvec_transp_doc[] = "a.matvec_transp(x, y)\n\
\n\
compute the sparse matrix-vector product y := a^T * x. \n\
a^T is the transpose of a, which is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional Numeric arrays of appropriate size.";

/** table of object methods
 */
PyMethodDef SSSMat_methods[] = {
  {"matvec", (PyCFunction)SSSMat_matvec, METH_VARARGS, SSSMat_matvec_doc},
  {"matvec_transp", (PyCFunction)SSSMat_matvec, METH_VARARGS, SSSMat_matvec_transp_doc},
  {NULL, NULL}			/* sentinel */
};

/*********************************************************************** 
 * SSSMatType methods
 */

static void
SSSMatType_dealloc(SSSMatObject *self)
{
  PyMem_DEL(self->ind);
  PyMem_DEL(self->val);
  PyMem_DEL(self->col);
  PyMem_DEL(self->diag);  
  PyObject_Del(self);
}

static int
SSSMatType_print(SSSMatObject *a, FILE *fp, int flags)
{
  int i, k, first = 1;

  if (a->nnz == 0) {
    fprintf(fp, "sss_mat([%d,%d])", a->n, a->n);
    return 0;
  }
  fprintf(fp, "sss_mat([%d,%d], [", a->n, a->n);
  for (i = 0; i < a->n; i ++) {
    for (k = a->ind[i]; k < a->ind[i+1]; k ++) {
      if (!first)
	fprintf(fp, ", ");
      first = 0;
      fprintf(fp, "(%d,%d): %g", i, a->col[k], a->val[k]);
    }
    fprintf(fp, "(%d,%d): %g", i, i, a->diag[i]);
  }
  fprintf(fp, "])");
  return 0;
}

static PyObject *
SSSMatType_getattr(SSSMatObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "nnz") == 0)
    return PyInt_FromLong(self->nnz + self->n);
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
  return Py_FindMethod(SSSMat_methods, (PyObject *)self, name);
}

/***********************************************************************
 * mapping functions
 */

/** SSSMat_length - number of items in mapping
 *    == number of matrix entries
 */
static int SSSMat_length(SSSMatObject *self) {
  return self->n * self->n;
}

/** SSSMat_subscript
 *    Called when treating array object like a mapping. This is used
 *    implement two-dimensional idices, e.g. A[i,j] or A[i1:i2,j1:j2]
 */
static PyObject *
SSSMat_subscript(SSSMatObject *self, PyObject *idx) {
  int type, start0, stop0, start1, stop1;
  int dim[2];

  dim[0] = dim[1] = self->n;
  if ((type = LLMat_parse_index(idx, dim, &start0, &stop0, &start1, &stop1)) == -1)
    return NULL;
  if (type == 1)
    return getitem(self, start0, start1);
  else {
    PyErr_SetString(PyExc_IndexError, "slices not supported");
    return NULL;
  }
}

static PyMappingMethods SSSMat_as_mapping = {
    (inquiry)SSSMat_length,	/*mp_length*/
    (binaryfunc)SSSMat_subscript, /*mp_subscript*/
    (objobjargproc)0,
};

/***********************************************************************
 * SSSMatType structure
 */

static PyTypeObject SSSMatType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "sss_mat",
  sizeof(SSSMatObject),
  0,
  (destructor)SSSMatType_dealloc, /* tp_dealloc */
  (printfunc)SSSMatType_print,	/* tp_print */
  (getattrfunc)SSSMatType_getattr, /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  &SSSMat_as_mapping,		/* tp_as_mapping*/
  0,				/* tp_hash */
};

/** newSSSMatObject -- allocate a new SSSMatObject instance
 *
 *    a newly allocated, uninitialized SSSMatObject is returned
 */
static PyObject *
newSSSMatObject(int n, int nnz) {
  SSSMatObject *op;

  /* create new SparseArrayt object */
  op = PyObject_New(SSSMatObject, &SSSMatType);
  if (op == NULL)
    PyErr_NoMemory();

  op->val = NULL;
  op->diag = NULL;
  op->ind = NULL;
  op->col = NULL;

  /* allocate arrays */
  op->ind = PyMem_New(int, n + 1);
  if (op->ind == NULL)
    goto fail;
  op->diag = PyMem_New(double, n);
  if (op->diag == NULL)
    goto fail;
  op->val = PyMem_New(double, nnz);
  if (op->val == NULL)
    goto fail;
  op->col = PyMem_New(int, nnz);
  if (op->col == NULL)
    goto fail;

  /* initialize rest of fields */
  op->n = n;
  op->nnz = nnz;

  return (PyObject *) op;

 fail:
    PyMem_Del(op->ind);    
    PyMem_Del(op->diag);    
    PyMem_Del(op->val);    
    PyMem_Del(op->col);    
    PyObject_Del(op);
    return PyErr_NoMemory();
}
