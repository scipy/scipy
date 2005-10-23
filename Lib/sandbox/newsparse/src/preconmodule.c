#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL precon
#include "Numeric/arrayobject.h"
#include "pysparse/blas.h"
#include "pysparse/fortran.h"
#include "pysparse/spmatrix.h"

typedef struct JacobiObject {
  PyObject_VAR_HEAD
  int n;
  PyObject *matrix;
  double *dinv;
  double *temp;
  double omega;
  int steps;
} JacobiObject;

typedef struct SSORObject {
  PyObject_VAR_HEAD
  int n;
  SSSMatObject *matrix;
  double *temp;
  double *temp2;
  double omega;
  int steps;
} SSORObject;

/*********************************************************************** 
 * JacobiObject methods
 */

static int
jacobi(PyObject *matrix, int n, double *dinv, int steps, double *x, double *y, double *temp) {
  int ONE = 1;
  int i, step, res;
  
  /* 1st step */
  for (i = 0; i < n; i ++)
    y[i] = x[i]*dinv[i];
 
  /* following steps */
  for(step = 1; step < steps; step ++) {
    F77(dcopy)(&n, y, &ONE, temp, &ONE);
    res = SpMatrix_Matvec(matrix, n, temp, n, y);
    if (res == -1)
      return res;
    for (i = 0; i < n; i ++)
      y[i] = (x[i] - y[i])*dinv[i] + temp[i];
  }
  return 0;
}

static char Jacobi_precon_doc[] = "self.precon(x, y)\n\
\n\
apply preconditioner self on x, store result in y. x is unchanged.";

static PyObject *
Jacobi_precon(JacobiObject *self, PyObject *args) {
  PyArrayObject *xp, *yp;
  int res;
  
  /* parse input arguments */
  SPMATRIX_PARSE_ARGS_ARR_ARR(args, xp, yp, self->n, self->n);
  
  res = jacobi(self->matrix, 
	       self->n, self->dinv, self->steps, 
	       (double *)(xp->data), (double *)(yp->data), 
	       self->temp);

  if (res) {
    PyErr_SetString(PyExc_RuntimeError, "unknown error in Jacobi iteration");
    return NULL;
  } else {
    Py_INCREF(Py_None); 
    return Py_None;
  }
}

/** table of object methods
 */
PyMethodDef Jacobi_methods[] = {
  {"precon", (PyCFunction)Jacobi_precon, METH_VARARGS, Jacobi_precon_doc},
  {NULL, NULL}			/* sentinel */
};

/*********************************************************************** 
 * SSORObject methods
 */

/** SSOR_KERNEL -- SSOR iteration kernel
 */
static void 
ssor_kernel(int n, double *b, double *x, double *h, double *temp,
	     double *va, double *da, int *ja, int *ia,
	     double omega, int steps) {

  double s;                     /* runnung sum in innermost loop */
  int step;                     /* current iteration step */
  int i;                        /* current row of L */
  int j;                        /* column index */
  int k;                        /* loop variable */
 
  /* main iteration loop */
  for (step = 0; step < steps; step ++) {
 
    /* 1st half-step: x = (omega*L + D) \ (D*x - omega*(L^T + D)*x) and h = -omega*L*x */
    if (step == 0)
      for (i = 0; i < n; i ++)
	temp[i] = omega*b[i];
    else
      for (i = 0; i < n; i ++)
	temp[i] = (1.0 - omega)*x[i]*da[i] + h[i] + omega*b[i];

    for (i = 0; i < n; i ++) {
      s = 0.0;
      for (k = ia[i]; k < ia[i+1]; k ++) {
        j = ja[k];
        s -= va[k] * x[j];
      }
      h[i] = omega*s;
      x[i] = (temp[i] + h[i])/da[i];
    }
 
    /* 2nd half-step: x = (omega*L^T + D) \ (D*x - omega*(L+D)*x) and h = -omega*L^T*x */
    for (i = 0; i < n; i ++) {
      temp[i] = (1.0 - omega)*x[i]*da[i] + h[i] + omega*b[i];
      h[i] = 0.0;
    }
    for (i = n-1; i >= 0; i --) {
      h[i] = omega*h[i];
      x[i] = (temp[i] + h[i]) / da[i];
      s = x[i];
      for (k = ia[i]; k < ia[i+1]; k ++) {
        j = ja[k];
        h[j] -= va[k] * s;
      }
    }
 
  } /* end of main iteration loop */
}


/** SYMGS_KERNEL -- SYMGS iteration kernel, optimized version of
 *   SSOR_KERNEL for omega==1.0
 */
static void 
symgs_kernel(int n, double *b, double *x, double *y,
	     double *va, double *da, int *ja, int *ia,
	     int steps) {

  double s;                     /* runnung sum in innermost loop */
  int step;                     /* current iteration step */
  int i;                        /* current row of L */
  int j;                        /* column index */
  int k;                        /* loop variable */
 
  /* initialize iteration: y = 0 */
  for (k = 0; k < n; k ++)
    y[k] = 0.0;
 
  /* main iteration loop */
  for (step = 0; step < steps; step ++) {
 
    /* 1st half-step: x = (L + D) \ (b - y) and y = L*x */
    for (i = 0; i < n; i ++) {
      s = 0.0;
      for (k = ia[i]; k < ia[i+1]; k ++) {
        j = ja[k];
        s += va[k] * x[j];
      }
      x[i] = (b[i] - y[i] - s)/da[i];
      y[i] = s;
    }
 
    /* 2nd half-step: x = (L^T + D) \ (b - y) and y = L^T*x */
    for (k = 0; k < n; k ++) {
      x[k] = y[k];
      y[k] = 0.0;
    }
    for (i = n-1; i >= 0; i --) {
      x[i] = (b[i] - x[i] - y[i]) / da[i];
      s = x[i];
      for (k = ia[i]; k < ia[i+1]; k ++) {
        j = ja[k];
        y[j] += va[k] * s;
      }
    }
 
  } /* end of main iteration loop */
}

static char SSOR_precon_doc[] = "self.precon(x, y)\n\
\n\
apply preconditioner self on x, store result in y. x is unchanged.";

static PyObject *
SSOR_precon(SSORObject *self, PyObject *args) {
  PyArrayObject *xp, *yp;
  
  /* parse input arguments */
  SPMATRIX_PARSE_ARGS_ARR_ARR(args, xp, yp, self->n, self->n);

  if (self->omega == 1.0)
    /* call SSOR iteration kernel */
    symgs_kernel(self->n, (double *)(xp->data), (double *)(yp->data), 
		 self->temp,
		 self->matrix->val, self->matrix->diag,
		 self->matrix->col, self->matrix->ind,
		 self->steps);
  else
    /* call SSOR iteration kernel */
    ssor_kernel(self->n, (double *)(xp->data), (double *)(yp->data), 
		self->temp, self->temp2,
		self->matrix->val, self->matrix->diag,
		self->matrix->col, self->matrix->ind,
		self->omega, self->steps);
  
  Py_INCREF(Py_None); 
  return Py_None;
}

/** table of object methods
 */
PyMethodDef SSOR_methods[] = {
  {"precon", (PyCFunction)SSOR_precon, METH_VARARGS, SSOR_precon_doc},
  {NULL, NULL}			/* sentinel */
};

/*********************************************************************** 
 * JacobiType methods
 */

static void
Jacobi_dealloc(JacobiObject *self)
{
  Py_DECREF(self->matrix);
  PyMem_DEL(self->dinv);
  PyMem_DEL(self->temp);
  PyObject_Del(self);
}

static PyObject *
Jacobi_getattr(JacobiObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape"};
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
  return Py_FindMethod(Jacobi_methods, (PyObject *)self, name);
}

/*********************************************************************** 
 * SSORType methods
 */

static void
SSOR_dealloc(SSORObject *self)
{
  Py_DECREF(self->matrix);
  PyMem_DEL(self->temp);
  PyMem_DEL(self->temp2);
  PyObject_Del(self);
}

static PyObject *
SSOR_getattr(SSORObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->n, self->n);
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape"};
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
  return Py_FindMethod(SSOR_methods, (PyObject *)self, name);
}

/***********************************************************************
 * JacobiType structure
 */

PyTypeObject JacobiType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "jacobi_prec",
  sizeof(JacobiObject),
  0,
  (destructor)Jacobi_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)Jacobi_getattr,  /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
};

/***********************************************************************
 * SSORType structure
 */

PyTypeObject SSORType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "ssor_prec",
  sizeof(SSORObject),
  0,
  (destructor)SSOR_dealloc,     /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)SSOR_getattr,    /* tp_getattr */
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

static PyObject *
newJacobiObject(PyObject *matrix, double omega, int steps) {
  JacobiObject *op;
  int n;			/* order of matrix */
  double d;
  int i;

  /* check shape of matrix object */
  if (SpMatrix_GetOrder(matrix, &n))
    return NULL;
  
  /* create new JacobiObject */
  op = PyObject_New(JacobiObject, &JacobiType);
  if (op == NULL)
    return PyErr_NoMemory();
  op->matrix = matrix;
  op->n = n;
  op->omega = omega;
  op->steps = steps;
  op->dinv = NULL;
  op->temp = NULL;

  /* allocate temp array if steps > 1 */
  if (steps > 1) {
    op->temp = PyMem_New(double, n);
    if (op->temp == NULL) {
      PyErr_NoMemory();
      goto fail;
    }
  }
    
  /* allocate and fill dinv array */
  op->dinv = PyMem_New(double, n);
  if (op->dinv == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  for (i = 0; i < n; i ++) {
    /* get diagonal element (i,i) */
    d = SpMatrix_GetItem(matrix, i, i);
    if (PyErr_Occurred())
      goto fail;
    /* check for singularity */
    if (1.0 + d == 1.0) {
      PyErr_SetString(PyExc_ValueError, "diagonal element close to zero");
      goto fail;
    }
    /* store inverse */
    op->dinv[i] = omega / d;
  }

  /* increment ref count of matrix */
  Py_INCREF(matrix);
  return (PyObject *)op;

 fail:
  PyMem_Del(op->dinv);
  PyMem_Del(op->temp);
  PyObject_Del(op);
  return NULL;
}

static PyObject *
newSSORObject(SSSMatObject *matrix, double omega, int steps) {
  SSORObject *self;
  int n;			/* order of matrix */

  /* check shape of matrix object */
  if (SpMatrix_GetOrder((PyObject *)matrix, &n))
    return NULL;
  
  /* create new SSORObject */
  self = PyObject_New(SSORObject, &SSORType);
  if (self == NULL)
    return PyErr_NoMemory();
  self->matrix = matrix;
  self->n = n;
  self->omega = omega;
  self->steps = steps;
  self->temp = NULL;
  self->temp2 = NULL;

  /* allocate temp array */
  self->temp = PyMem_New(double, n);
  if (self->temp == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
    
  /* allocate temp2 array if omega != 0.0 */
  if (omega != 1.0) {
    self->temp2 = PyMem_New(double, n);
    if (self->temp2 == NULL) {
      PyErr_NoMemory();
      goto fail;
    }
  }
    
  /* increment ref count of matrix */
  Py_INCREF(matrix);
  return (PyObject *)self;

 fail:
  PyMem_Del(self->temp);
  PyMem_Del(self->temp2);
  PyObject_Del(self);
  return NULL;
}

/*********************************************************************** 
 * Module functions
 */

static char jacobi_doc[] = 
"jacobi(A, omega=1.0, steps=1)\n\
\n\
return Jacobi preconditioner object.";

static PyObject *
jacobi_prec(PyObject *self, PyObject *args)
{
  PyObject *matrix;
  double omega;
  int steps;

  /* parse input arguments */
  omega = 1.0;
  steps = 1;
  if (!PyArg_ParseTuple(args, "O|di", &matrix, &omega, &steps))
    return NULL;

  /* construct SSORObject */
  return newJacobiObject(matrix, omega, steps);
}

static char ssor_doc[] = 
"ssor(A, omega, steps) -- return SSOR preconditioner object\n\
\n\
This preconditioner executes 'steps' SSOR steps with a zero initial guess.\n\
\n\
parameters\n\
----------\n\
\n\
A      'sss_mat' object, symmetric sparse matrix\n\
omega  relaxation parameter (default value: 1.0)\n\
steps  number of SSOR steps";

static PyObject *
ssor(PyObject *self, PyObject *args)
{
  SSSMatObject *matrix;
  double omega;
  int steps;

  /* parse input arguments */
  omega = 1.0;
  steps = 1;
  if (!PyArg_ParseTuple(args, "O!|di", &SSSMatType, &matrix, &omega, &steps))
    return NULL;

  /* construct JacobiObject */
  return newSSORObject(matrix, omega, steps);
}


/** table of module functions
 */
static PyMethodDef precon_methods[] = {
  {"jacobi", (PyCFunction)jacobi_prec, METH_VARARGS, jacobi_doc},
  {"ssor", (PyCFunction)ssor, METH_VARARGS, ssor_doc},
  {NULL, NULL}	/* sentinel */
};


DL_EXPORT(void)
initprecon(void)
{
  PyObject *m, *d;
  
  JacobiType.ob_type = &PyType_Type;
  SSORType.ob_type = &PyType_Type;

  m = Py_InitModule("precon", precon_methods);
  d = PyModule_GetDict(m);

  PyDict_SetItemString(d, "JacobiType", (PyObject *)&JacobiType);
  PyDict_SetItemString(d, "SSORType", (PyObject *)&SSORType);

  /* initialize Numeric array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();

  /* No need to check the error here, the caller will do that */
}
