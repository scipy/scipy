#include "Python.h"

#define SPMATRIX_MODULE
#include "pysparse/spmatrix.h"

#define PY_ARRAY_UNIQUE_SYMBOL spmatrix
#include "numpy/arrayobject.h"

PyObject *SpMatrix_ErrorObject;

/******************************************************************************
 *                                                                            *
 * Include C code from other C files                                          *
 *                                                                            *
 ******************************************************************************/
#include "mmio_patched.c"
#include "ll_mat.c"
#include "csr_mat.c"
#include "sss_mat.c"

/******************************************************************************
 ******************************************************************************
 *                                                                            *
 * C API functions of the spmatrix package                                    *
 *                                                                            *
 ******************************************************************************
 ******************************************************************************/


/******************************************************************************
 *                                                                            *
 * SpMatrix_ParseVecOpArgs -- parse arguments                                 *
 *                                                                            *
 *   parses arguments of Python functions, that expect two one dimensional    *
 *   Py_Array objects.                                                        *
 *                                                                            *
 *   This function should not be used, since it doesn't allow free the        *
 *   created objects.                                                         *
 *                                                                            *
 *   Use the macro SPMATRIX_PARSE_ARGS_ARR_ARR defined in spmatrix.h instead. *
 *                                                                            *
 ******************************************************************************/
static int
SpMatrix_ParseVecOpArgs(PyObject *args, double **x_data, double **y_data, int n) {
  PyObject *x_obj, *y_obj;
  int nx, ny, res;
  
  /* parse input arguments */
  if (!PyArg_ParseTuple(args, "OO", &x_obj, &y_obj))
    return -1;

  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x_obj, (char **)x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert first argument to double array");
    return -1;
  }
  res = PyArray_As1D((PyObject **)&y_obj, (char **)y_data, &ny, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert second argument to double array");
    return -1;
  }
  
  /* check operand shapes */
  if (nx != n || ny != n) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return -1;
  }

  return 0;
}

/******************************************************************************
 *                                                                            *
 * SpMatrix_GetShape -- get object shape                                      *
 *                                                                            *
 *   Query the shape attribute of an object.                                  *
 *                                                                            *
 *   Used to get the number of row and columns of a sparse matrix or a        *
 *   preconditioner.                                                          *
 *                                                                            *
 ******************************************************************************/
static int
SpMatrix_GetShape(PyObject *op, int dim[]) {
  PyObject *sh, *elem;
  
  if ((sh = PyObject_GetAttrString(op, "shape")) == NULL)
    return -1;
  if (PySequence_Size(sh) != 2) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return -1;
  }
  elem = PySequence_GetItem(sh, 0); dim[0] = PyInt_AsLong(elem); Py_DECREF(elem);
  elem = PySequence_GetItem(sh, 1); dim[1] = PyInt_AsLong(elem); Py_DECREF(elem);
  Py_DECREF(sh);
  if (PyErr_Occurred() != NULL) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return -1;
  }
  return 0;
}

/******************************************************************************
 *                                                                            *
 * SpMatrix_GetOrder -- get object order                                      *
 *                                                                            *
 *   Queries the shape attribute of an object.                                *
 *                                                                            *
 *   Used to get the number of row and columns of a sparse matrix or a        *
 *   preconditioner.                                                          *
 *                                                                            *
 *   Fails if the shape is not square.                                        *
 *                                                                            *
 ******************************************************************************/
static int
SpMatrix_GetOrder(PyObject *op, int *n) {
  int dim[2];

  if (SpMatrix_GetShape(op, dim) == -1)
    return -1;
  
  if (dim[0] != dim[1]) {
    PyErr_SetString(PyExc_ValueError, "matrix is not square");
    return -1;
  }

  *n = dim[0];
  return 0;
}

/******************************************************************************
 *                                                                            *
 * SpMatrix_GetItem -- access matrix entry (i,j)                              *
 *                                                                            *
 *      returns matrix entry op[i,j] as a double                              *
 *                                                                            *
 ******************************************************************************/
static double 
SpMatrix_GetItem(PyObject *op, int i, int j) {
  PyObject *index;
  PyObject *fo;
  double d;

  index = Py_BuildValue("(ii)", i, i);
  assert(index);
  fo = PyObject_GetItem(op, index);
  Py_DECREF(index);
  if (fo == NULL)
    return 0.0;
  d = PyFloat_AsDouble(fo);
  Py_DECREF(fo);
  return d;
}

/******************************************************************************
 *                                                                            *
 * SpMatrix_Matvec -- invoke matrix-vector multiplication                     *
 *                                                                            *
 *   Invokes matrix-vector multiplication, by calling the 'matvec'-method of  *
 *   'matrix',  y = matrix*x. The vectors 'x' and 'y' are given as arrays of  *
 *   double.                                                                  *
 *                                                                            *
 *   Returns 0 if the operation was successful, or -1 if an error occured.    *
 *                                                                            *
 ******************************************************************************/
static int
SpMatrix_Matvec(PyObject *matrix, int nx, double *x, int ny, double *y) {
  PyObject *x_arr = NULL;
  PyObject *y_arr = NULL;
  PyObject *res;

  /* create array objects from x and y */
  x_arr = PyArray_FromDimsAndData(1, &nx, PyArray_DOUBLE, (char *)x);
  if (x_arr == NULL)
    goto fail;
  y_arr = PyArray_FromDimsAndData(1, &ny, PyArray_DOUBLE, (char *)y);
  if (y_arr == NULL)
    goto fail;

  /* Call matvec method of matrix object */
  res = PyObject_CallMethod(matrix, "matvec", "OO", x_arr, y_arr);
  if (res == NULL)
    goto fail;
  Py_DECREF(res);
  
  /* free array objects */
  Py_DECREF(x_arr);
  Py_DECREF(y_arr);
  return 0;

 fail:
  Py_XDECREF(x_arr);
  Py_XDECREF(y_arr);
  return -1;
}


/******************************************************************************
 *                                                                            *
 * SpMatrix_Precon -- invoke preconditioner                                   *
 *                                                                            *
 *  Applies the preconditioner 'prec' on the vector 'x' and stores the result *
 *  in vector 'y'. This is done by calling the 'precon' method of 'prec'. The *
 *  vectors 'x' and 'y' are given as arrays of double.                        *
 *                                                                            *
 *  Returns 0 if the operation was successful, or -1 if an error occured.     *
 *                                                                            *
 ******************************************************************************/
static int
SpMatrix_Precon(PyObject *prec, int n, double *x, double *y) {
  PyObject *x_arr = NULL;
  PyObject *y_arr = NULL;
  int dimensions[1];
  PyObject *res;

  dimensions[0] = n;

  /* create array objects from x and y */
  x_arr = PyArray_FromDimsAndData(1, dimensions, PyArray_DOUBLE, (char *)x);
  if (x_arr == NULL)
    goto fail;
  y_arr = PyArray_FromDimsAndData(1, dimensions, PyArray_DOUBLE, (char *)y);
  if (y_arr == NULL)
    goto fail;

  /* Call matvec method of matrix object */
  res = PyObject_CallMethod(prec, "precon", "OO", x_arr, y_arr);
  if (res == NULL)
    goto fail;

  Py_DECREF(res);
  
  /* free array objects */
  Py_DECREF(x_arr);
  Py_DECREF(y_arr);
  return 0;

 fail:
  Py_XDECREF(x_arr);
  Py_XDECREF(y_arr);
  return -1;
}

/******************************************************************************
 *                                                                            *
 * Itsolvers_Solve -- Invoke linear solver                                    *
 *                                                                            *
 * Invoke iterative linear solver 'linsolver', to (approximately) solve the   *
 * linear system                                                              *
 *                                                                            *
 *     A * x = b                                                              *
 *                                                                            *
 * to an accuracy of 'tol'. The maximum number of iteration steps taken is    *
 * 'itmax'. The vectors 'x' and 'y' are given as arrays of double of length   *
 * 'n'.                                                                       *
 *                                                                            *
 * Returns 0 if the operation was successful, or -1 if an error occured.      *
 *                                                                            *
 ******************************************************************************/
static int
ItSolvers_Solve(PyObject *linsolver, PyObject *A, int n, 
		double *b, double *x, double tol, int itmax, PyObject *K,
		int *info, int *iter, double *relres) {
  PyObject *b_arr = NULL;
  PyObject *x_arr = NULL;
  int dimensions[1];
  PyObject *res;

  dimensions[0] = n;

  /* create array objects from x and y */
  b_arr = PyArray_FromDimsAndData(1, dimensions, PyArray_DOUBLE, (char *)b);
  if (b_arr == NULL)
    goto fail;
  x_arr = PyArray_FromDimsAndData(1, dimensions, PyArray_DOUBLE, (char *)x);
  if (x_arr == NULL)
    goto fail;

  /* Call iterative solver */
  if (K == NULL)
    res = PyObject_CallFunction(linsolver, "OOOdi", A, b_arr, x_arr, tol, itmax);
  else
    res = PyObject_CallFunction(linsolver, "OOOdiO", A, b_arr, x_arr, tol, itmax, K);
  if (res == NULL)
    goto fail;

  /* Parse result 
     Abuse PyArg_ParseTuple to parse res tuple (is this safe?) */
  PyArg_ParseTuple(res, "iid", info, iter, relres);
  Py_DECREF(res);
  
  /* free array objects */
  Py_DECREF(b_arr);
  Py_DECREF(x_arr);
  return 0;

 fail:
  Py_XDECREF(b_arr);
  Py_XDECREF(x_arr);
  return -1;
}


/******************************************************************************
 *                                                                            *
 * LLMatType_alloc                                                            *
 *                                                                            *
 *   construct an ll_mat object that represents an empty sparse matrix.       *
 *                                                                            *
 ******************************************************************************/
static PyObject *
LLMat_zeros(PyObject *self, PyObject *args)
{
  int dim[2], sizeHint;
  
  sizeHint = 1000;
  if (!PyArg_ParseTuple(args, "ii|i", dim, dim + 1, &sizeHint))
    return NULL;

  return SpMatrix_NewLLMatObject(dim, 0, sizeHint);
}

static PyObject *
LLMat_sym_zeros(PyObject *self, PyObject *args)
{
  int dim[2], n, sizeHint;
  
  sizeHint = 1000;
  if (!PyArg_ParseTuple(args, "i|i", &n, &sizeHint))
    return NULL;
  dim[0] = dim[1] = n;

  return SpMatrix_NewLLMatObject(dim, 1, sizeHint);
}

/******************************************************************************
 *                                                                            *
 * table of module functions                                                  *
 *                                                                            *
 ******************************************************************************/
static PyMethodDef spmatrix_methods[] = {
  {"ll_mat",          LLMat_zeros,          METH_VARARGS, 0},
  {"ll_mat_sym",      LLMat_sym_zeros,      METH_VARARGS, 0},
  {"ll_mat_from_mtx", LLMat_from_mtx,       METH_VARARGS, 0},
  {"matrixmultiply",  LLMat_matrixmultiply, METH_VARARGS, LLMat_matrixmultiply_doc},
  {"dot",             LLMat_dot,            METH_VARARGS, LLMat_dot_doc},
  {NULL, NULL}	/* sentinel */
};

DL_EXPORT(void)
initspmatrix(void)
{
  PyObject *m, *d;

  LLMatType.ob_type = &PyType_Type;
  CSRMatType.ob_type = &PyType_Type;
  SSSMatType.ob_type = &PyType_Type;
  m = Py_InitModule("spmatrix", spmatrix_methods);
  if (m == NULL)
    goto fail;
  d = PyModule_GetDict(m);
  if (d == NULL)
    goto fail;
  PyDict_SetItemString(d, "LLMatType", (PyObject *)&LLMatType);
  PyDict_SetItemString(d, "CSRMatType", (PyObject *)&CSRMatType);
  PyDict_SetItemString(d, "SSSMatType", (PyObject *)&SSSMatType);
  SpMatrix_ErrorObject = PyString_FromString("spmatrix.error");
  PyDict_SetItemString(d, "error", SpMatrix_ErrorObject);
 
  /* initialise C API */
  init_c_api(d);

  /* initialise scipy array module */
  import_array();

  /* Check for errors */
  if (PyErr_Occurred())
    goto fail;

  return;

fail:
  Py_FatalError("can't initialize module spmatrix");
}
