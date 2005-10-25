/** MODULE jdsym
 *
 */

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL jdsym_API
#include "scipy/arrayobject.h"
#include "pysparse/fortran.h"
#include "pysparse/blas.h"
#include "pysparse/lapack.h"
#include "pysparse/spmatrix.h"


/*********************************************************************** 
 * Static constant variables
 */

static double DONE  =  1.0;
static double DZER  =  0.0;
static double DMONE = -1.0;
static int    MONE  = -1;
static int    ONE   =  1;


/*********************************************************************** 
 * Helper code for jdsym: Calls the 'proj' method of a python object
 * passing a double array.
 */

static int
Jdsym_Proj(PyObject *proj, int n, double *x) {
  PyObject *x_arr = NULL;
  int dimensions[1];
  PyObject *res;

  dimensions[0] = n;

  /* create array objects from x and y */
  x_arr = PyArray_FromDimsAndData(1, dimensions, PyArray_DOUBLE, (char *)x);
  if (x_arr == NULL)
    goto fail;

  /* Call "project" method of "proj" object */
  res = PyObject_CallMethod(proj, "project", "O", x_arr);
  if (res != Py_None) {
    if (!PyErr_Occurred())
      PyErr_SetString(PyExc_RuntimeError, "error when calling projector");
    goto fail;
  }
  Py_DECREF(res);
  
  /* free array objects */
  Py_DECREF(x_arr);
  return 0;

 fail:
  Py_XDECREF(x_arr);
  return -1;
}


/*********************************************************************** 
 * Include additional C code
 */

#include "orthopack.c"
#include "correq.c"
#include "jdsym.c"


/*********************************************************************** 
 * Module functions
 */

static char jdsym_doc[] = "\
Jacobi davidson eigenvalue solver for symmetric matrices\n\
\n\
kconv, lmbd, Q, it, it_inner = jdsym(amat, mmat, prec, kmax, tau, jdtol, itmax, linsolver, ...)\n\
\n\
jdsym computes kmax eigenpairs in the vincinity of target tau of the generalized \n\
eigenvalue problem\n\
\n\
   amat * x = lambda * mmat * x\n\
\n\
If mmat is None, jdsym computes eigenpairs of the standart eigenvalue problem\n\
\n\
   amat * x = lambda * x\n\
\n\
Required arguments:\n\
\n\
   amat: matrix object. Must be symmetric.\n\
\n\
   mmat: matrix object or None. Must be symmetric positive definite.\n\
\n\
   prec: preconditioner object or None. The preconditioner must be symmetric and\n\
         should approximate amat - tau mmat.\n\
\n\
   kmax: number of eigenpairs to be computed\n\
\n\
   tau: float value, target value\n\
\n\
   jdtol: required error tolerance for computed eigensolutions\n\
\n\
   itmax: maximum number of JD iterations to be performed\n\
\n\
   linsolver: callable object which implements a linear solver such as\n\
              itsolvers.cgs or itsolvers.qmrs\n\
\n\
Additional keyword arguments:\n\
\n\
   jmax: integer, maximal dimension of search subspace (default value: 25)\n\
\n\
   jmin: integer, dimension of search subspace after restart (default value: 10)\n\
\n\
Result values:\n\
\n\
   kconv: integer, number of converged eigenpairs\n\
\n\
   lmbd: double array of dimension kconv, eigenvalues\n\
\n\
   Q: double array of dimension n x kconv, eigenvectors, \n\
    lmbd[k], Q[:,k] is the kth eigenpair\n\
\n\
   it: integer, number of JD iterations\n\
\n\
   it_inner: number of inner iterations performed";

static PyObject *
JDSym_jdsym(PyObject *self, PyObject *args, PyObject *keywds) {
  /* required arguments */
  PyObject *amat;
  PyObject *mmat;
  PyObject *prec;
  int kmax;
  double tau;
  double jdtol;
  int itmax;
  PyObject* linsolver;
  /* optional arguments (default values) */
  int jmax = 25;
  int jmin = 10;
  int blksize = 1;
  int blkwise = 0;
  PyArrayObject *V0 = NULL;
  int optype = 2;
  int linitmax = 200;
  double eps_tr = 1.0e-3;
  double toldecay = 1.5;
  int clvl = 0;
  int strategy = 0;
  PyObject *proj = NULL;
  /* result */
  int kconv;
  int it_outer;
  int it_inner;
  double *Q = NULL;
  double *lambda = NULL;
  PyArrayObject *Q_obj = NULL;
  PyArrayObject *lambda_obj = NULL;
  PyObject *result;
  /* other local variables */
  int n, na, nm, nk, np;	/* order of eigensystem */
  int i, k, ret;
  int dimensions[2];
  
  static char *kwlist[] = {
    "A", "M", "K", "kmax", "tau", "jdtol", "itmax", "linsolver", 
    "jmax", "jmin", "blksize", "blkwise", "V0", "optype", "linitmax", 
    "eps_tr", "toldecay", "clvl", "strategy", "projector", NULL};

  ret = PyArg_ParseTupleAndKeywords(args, keywds, "OOOiddiO|iiiiO!iiddiiO", kwlist,
				    /* required args */
				    &amat,
				    &mmat,
				    &prec,
				    &kmax,
				    &tau,
				    &jdtol,
				    &itmax,
				    &linsolver,
				    /* optional args */
				    &jmax,
				    &jmin,
				    &blksize,
				    &blkwise,
				    &PyArray_Type,
				    &V0,
				    &optype,
				    &linitmax,
				    &eps_tr,
				    &toldecay,
				    &clvl,
				    &strategy,
				    &proj);
  if (!ret)
    return NULL;
  
  /* 
   *  Check shapes of matrix, preconditioner and projector objects 
   */
  if (SpMatrix_GetOrder((PyObject *)amat, &na))
    return NULL;
  if (mmat != Py_None) {
    if (SpMatrix_GetOrder((PyObject *)mmat, &nm))
      return NULL;
  } else
    nm = na;
  if (prec != Py_None) {
    if (SpMatrix_GetOrder((PyObject *)prec, &nk))
      return NULL;
  } else
    nk = na;
  if (proj != NULL) {
    if (SpMatrix_GetOrder((PyObject *)proj, &np))
      return NULL;
  } else
    np = na;
  if (na != nm || na != nk || na != np) {
    PyErr_SetString(PyExc_ValueError, "matrix, preconditioner or projector shapes differ");
    return NULL;
  }
  n = na;

  /* Check V0 argument 
   *
   * Here the V0 argument should be converted to Fortran ordering
   */
  if (V0 != NULL) {
    if (!((V0->nd == 1  || V0->nd == 2) && 
	  V0->descr->type_num == PyArray_DOUBLE && 
	  V0->dimensions[0] == n)) {
      PyErr_SetString(PyExc_ValueError, "V0 is not of correct type or shape");
      goto fail;
    }
  }
  
  /* allocate (max.) space for eigenvectors and eigenvalues */
  Q = PyMem_New(double, n*kmax);
  if (Q == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  lambda = PyMem_New(double, kmax);
  if (lambda == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  
  /* call Jacobi-Davidson eigensolver */
  ret = jdsym(n, tau, jdtol, 
	      kmax, jmax, jmin, itmax,
	      blksize, blkwise,
	      V0, 
	      linsolver,
	      optype, 
	      linitmax, eps_tr, toldecay,
	      strategy, clvl,
	      &kconv, Q, lambda, &it_outer, &it_inner,
	      amat,
	      mmat == Py_None ? NULL : mmat,
	      prec == Py_None ? NULL : prec,
	      proj);
  if (ret == -1)
    goto fail;
  
  /* 
   * prepare and return results 
   */
  /* allocate array objects for *converged* eigenvectors */
  dimensions[0] = n; dimensions[1] = kconv;
  Q_obj = (PyArrayObject *)PyArray_FromDims(2, dimensions, PyArray_DOUBLE);
  if (Q_obj == NULL)
    goto fail;
  /* copy converged eigenvectors, convert from Fortran to C ordering */
  for (k = 0; k < kconv; k ++)
    for (i = 0; i < n; i ++)
      ((double *)Q_obj->data)[i*kconv + k] = Q[k*n + i];
  /* allocate array objects for *converged* eigenvalues */
  dimensions[0] = kconv;
  lambda_obj = (PyArrayObject *)PyArray_FromDims(1, dimensions, PyArray_DOUBLE);
  if (lambda_obj == NULL)
    goto fail;
  for (k = 0; k < kconv; k ++)
    ((double *)lambda_obj->data)[k] = lambda[k];
  /* delete original arrays */
  PyMem_DEL(Q);
  PyMem_DEL(lambda);
  result =  Py_BuildValue("iOOii", kconv, lambda_obj, Q_obj, it_outer, it_inner);
  /* Py_BuildValue increased refcount, so decrease it or it would never be freed */
  Py_DECREF(lambda_obj);
  Py_DECREF(Q_obj);
  return result;

 fail:
  PyMem_DEL(Q);
  PyMem_DEL(lambda);
  Py_XDECREF(Q_obj);
  Py_XDECREF(lambda_obj);
  return NULL;
}


/** table of module functions
 */
static PyMethodDef jdsym_methods[] = {
  {"jdsym",    (PyCFunction)JDSym_jdsym,  METH_VARARGS|METH_KEYWORDS, jdsym_doc},
  {NULL, NULL}	/* sentinel */
};

static char module_doc[] = "This module ...\n\
";

DL_EXPORT(void)
initjdsym(void) {
  PyObject *m;
  
  m = Py_InitModule3("jdsym", jdsym_methods, module_doc);

  /* initialize scipy array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();

  /* No need to check the error here, the caller will do that */
}
