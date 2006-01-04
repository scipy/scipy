#include "Python.h"

#define SPMATRIX_MODULE
#include "pysparse/spmatrix.h"

#define PY_ARRAY_UNIQUE_SYMBOL spmatrix
#include "numpy/arrayobject.h"

#define UNROLL_LOOPS 0

/*********************************************************************** 
 * CSRMatObject methods
 */

static void
csr_matvec_kernel(int m, double *x, double *y,
		  double *va, int *ja, int *ia) {
  double s;
  int i, k;
  
#if UNROLL_LOOPS
  int k1, *ja_ptr;

  ja_ptr = ja;
  k = ia[0];
  for (i = 0; i < m; i ++) {
    s = 0.0;
    k1 = ia[i+1];
    for (; k+10 < k1; k += 10) {
      s += va[k+9]*x[ja_ptr[9]] +
        va[k+8]*x[ja_ptr[8]] +
        va[k+7]*x[ja_ptr[7]] +
        va[k+6]*x[ja_ptr[6]] +
        va[k+5]*x[ja_ptr[5]] +
        va[k+4]*x[ja_ptr[4]] +
        va[k+3]*x[ja_ptr[3]] +
        va[k+2]*x[ja_ptr[2]] +
        va[k+1]*x[ja_ptr[1]] +
        va[k]*x[*ja_ptr];
      ja_ptr += 10;
    }
    for (; k < ia[i+1]; k ++)
      s += va[k]*x[*ja_ptr++];
    y[i] = s;
  }
#else
  for (i = 0; i < m; i ++) {
    s = 0.0;
    for (k = ia[i]; k < ia[i+1]; k ++)
      s += va[k] * x[ja[k]];
    y[i] = s;
  }
#endif
}

static void
csr_matvec_kernel_stride(int m, 
			 double *x, int incx, 
			 double *y, int incy,
			 double *va, int *ja, int *ia) {
  double s;
  int i, k;
  
  for (i = 0; i < m; i ++) {
    s = 0.0;
    for (k = ia[i]; k < ia[i+1]; k ++)
      s += va[k] * x[ja[k]*incx];
    y[i*incy] = s;
  }
}

static void
csr_matvec_transp_kernel(int m, int n, double *x, double *y,
			 double *va, int *ja, int *ia) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i];
    for (k = ia[i];  k < ia[i+1]; k ++)
      y[ja[k]] += va[k] * xi;
  }
}

static void
csr_matvec_transp_kernel_stride(int m, int n, 
				double *x, int incx, 
				double *y, int incy,
				double *va, int *ja, int *ia) {
  double xi;
  int i, k;
  
  for (i = 0; i < n; i ++)
    y[i*incy] = 0.0;
  
  for (i = 0; i < m; i ++) {
    xi = x[i*incx];
    for (k = ia[i];  k < ia[i+1]; k ++)
      y[ja[k]*incy] += va[k] * xi;
  }
}

static char matvec_transp_doc[] = "a.matvec_transp(x, y)\n\
\n\
compute the sparse matrix-vector product y := a^T * x. \n\
a^T is the transpose of a, which is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional scipy arrays of appropriate size.";

static PyObject *
CSRMat_matvec_transp(CSRMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[0], self->dim[1]);
  
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    csr_matvec_transp_kernel(self->dim[0], self->dim[1], (double *)(xp->data), (double *)(yp->data), 
			     self->val, self->col, self->ind);
  else {
    csr_matvec_transp_kernel_stride(self->dim[0], self->dim[1], 
				    (double *)(xp->data), xp->strides[0] / sizeof(double),
				    (double *)(yp->data), yp->strides[0] / sizeof(double),
				    self->val, self->col, self->ind);
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

static char matvec_doc[] = "a.matvec(x, y)\n\
\n\
compute the sparse matrix-vector product y := a * x. \n\
a is a d1 by d2 sparse matrix.\n\
x and y are two 1-dimensional scipy arrays of appropriate size.";

static PyObject *
CSRMat_matvec(CSRMatObject *self, PyObject *args)
{
  PyArrayObject *xp, *yp;

  SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, xp, yp, self->dim[1], self->dim[0]);
     
  if (xp->flags & CONTIGUOUS &&  yp->flags & CONTIGUOUS)
    csr_matvec_kernel(self->dim[0], (double *)(xp->data), (double *)(yp->data), 
		      self->val, self->col, self->ind);
  else {
    int incx = xp->strides[0] / sizeof(double);
    int incy = yp->strides[0] / sizeof(double);
    
    csr_matvec_kernel_stride(self->dim[0], 
			     (double *)(xp->data), incx, 
			     (double *)(yp->data), incy,
			     self->val, self->col, self->ind);
  }

  Py_INCREF(Py_None); 
  return Py_None;
}

/** table of object methods
 */
PyMethodDef CSRMat_methods[] = {
  {"matvec", (PyCFunction)CSRMat_matvec, METH_VARARGS, matvec_doc},
  {"matvec_transp", (PyCFunction)CSRMat_matvec_transp, METH_VARARGS, matvec_transp_doc},
  {NULL, NULL}			/* sentinel */
};

/*********************************************************************** 
 * CSRMatType methods
 */

static void
CSRMatType_dealloc(CSRMatObject *a)
{
  PyMem_DEL(a->ind);
  PyMem_DEL(a->val);
  PyMem_DEL(a->col);
  PyObject_Del(a);
}

static int
CSRMatType_print(CSRMatObject *a, FILE *fp, int flags)
{
  int i, k, first = 1;

  if (a->nnz == 0) {
    fprintf(fp, "csr_mat([%d,%d])", a->dim[0], a->dim[1]);
    return 0;
  }
  fprintf(fp, "csr_mat([%d,%d], [", a->dim[0], a->dim[1]);
  for (i = 0; i < a->dim[0]; i ++) {
    for (k = a->ind[i]; k < a->ind[i+1]; k ++) {
      if (!first)
	fprintf(fp, ", ");
      first = 0;
      fprintf(fp, "(%d,%d): %g", i, a->col[k], a->val[k]);
    }
  }
  fprintf(fp, "])");
  return 0;
}

static PyObject *
CSRMatType_getattr(CSRMatObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->dim[0], self->dim[1]);
  if (strcmp(name, "nnz") == 0)
    return PyInt_FromLong(self->nnz);
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
  return Py_FindMethod(CSRMat_methods, (PyObject *)self, name);
}

/***********************************************************************
 * CSRMatType structure
 */

static PyTypeObject CSRMatType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "csr_mat",
  sizeof(CSRMatObject),
  0,
  (destructor)CSRMatType_dealloc, /* tp_dealloc */
  (printfunc)CSRMatType_print,	/* tp_print */
  (getattrfunc)CSRMatType_getattr, /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
};

/** newCSRMatObject -- allocate a new CSRMatObject instance
 *
 *    a newly allocated, uninitialized CSRMatObject is returned
 */
static PyObject *
newCSRMatObject(int dim[], int nnz) {
  CSRMatObject *op;

  /* create new SparseArrayt object */
  op = PyObject_New(CSRMatObject, &CSRMatType);
  if (op == NULL)
    PyErr_NoMemory();

  op->val = NULL;
  op->ind = NULL;
  op->col = NULL;

  /* allocate arrays */
  op->ind = PyMem_New(int, dim[0] + 1);
  if (op->ind == NULL)
    goto fail;
  op->val = PyMem_New(double, nnz);
  if (op->val == NULL)
    goto fail;
  op->col = PyMem_New(int, nnz);
  if (op->col == NULL)
    goto fail;

  /* initialize rest of fields */
  op->dim[0] = dim[0];
  op->dim[1] = dim[1];
  op->nnz = nnz;

  return (PyObject *) op;

 fail:
    PyMem_Del(op->ind);    
    PyMem_Del(op->val);    
    PyMem_Del(op->col);    
    PyObject_Del(op);
    return PyErr_NoMemory();
}
