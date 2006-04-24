/** MODULE itsolvers
 *
 */

#include "Python.h"

#define PY_ARRAY_UNIQUE_SYMBOL itsolvers_pyarray
#include "numpy/arrayobject.h"

#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"

#include "pysparse/pcg.h"
#include "pysparse/bicgstab.h"
#include "pysparse/minres.h"
#include "pysparse/gmres.h"
#include "pysparse/qmrs.h"
#include "pysparse/cgs.h"

/*********************************************************************** 
 * Module functions
 */

static char pcg_doc[] = 
"info, iter, relres = pcg(A, b, x, tol, maxit)\n\
\n\
Preconditioned Conjugate Gradient method.";

static PyObject *
ItSolvers_pcg(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  int clvl = 0;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|O",
			&amat, 
			&b,
			&x,
			&tol,
			&maxit,
			&precon))
    return NULL;

  /* check shape of matrix object */
  res = SpMatrix_GetOrder((PyObject *)amat, &n);
  if (res)
    return NULL;
  
  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, 4*n);

  /* call pcg routine */
  Itsolvers_pcg_kernel(n,
		       x_data,
		       b_data,
		       tol,
		       maxit,
		       clvl,
		       &iter,
		       &relres,
		       &info,
		       work,
		       amat,
		       precon == Py_None ? NULL : precon);  
  
  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

static char bicgstab_doc[] = 
"info, iter, relres = bicgstab(A, b, x, tol, maxit)\n\
\n\
Preconditioned Stabilized BiConjugate Gradient method.";

static PyObject *
ItSolvers_bicgstab(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  int clvl = 0;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|O",
			&amat, 
			&b,
			&x,
			&tol,
			&maxit,
			&precon))
    return NULL;

  /* check shape of matrix object */
  res = SpMatrix_GetOrder((PyObject *)amat, &n);
  if (res)
    return NULL;
  
  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, 8*n);

  /* call bicgstab routine */
  Itsolvers_bicgstab_kernel(n,
			    x_data,
			    b_data,
			    tol,
			    maxit,
			    clvl,
			    &iter,
			    &relres,
			    &info,
			    work,
			    amat,
			    precon == Py_None ? NULL : precon);  
  
  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

static char minres_doc[] = 
"info, iter, relres = minres(A, b, x, tol, maxit, K)\n\
\n\
Minimal Residual method.";

static PyObject *
ItSolvers_minres(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  int clvl = 0;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  int dim[2];			/* shape of amat */
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|O",
			&amat, 
			&b,
			&x,
			&tol,
			&maxit,
			&precon))
    return NULL;

  /* check shape of matrix object */
  SpMatrix_GetShape((PyObject *)amat, dim);
  if (dim[0] != dim[1] || dim[0] <= 0) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return NULL;
  }
  n = dim[0];

  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, 7*n);

  /* call minres routine */
  info = Itsolvers_minres_kernel(n, 
				 tol, 
				 maxit,
				 &iter, 
				 &relres, 
				 clvl,
				 x_data, 
				 b_data, 
				 work,
				 amat,
				 precon == Py_None ? NULL : precon);
  
  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

static char gmres_doc[] =
"info, iter, relres = gmres(A, b, x, tol, maxit, K, dim)\n\
\n\
General Minimal Residual method GMRES(m) (m=dim) of Saad and Schultz.";

static PyObject *
ItSolvers_gmres(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  int dim_gmres = 20;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  int dim[2];			/* shape of amat */
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|Oi",
			&amat,
			&b,
			&x,
			&tol,
			&maxit,
			&precon,&dim_gmres))
    return NULL;

  /* check shape of matrix object */
  SpMatrix_GetShape((PyObject *)amat, dim);
  if (dim[0] != dim[1] || dim[0] <= 0) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return NULL;
  }
  n = dim[0];

  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, n);

  /* call gmres routine */
  info = Itsolvers_gmres_kernel(n,
				 tol,
				 maxit,
				 &iter,
				 &relres,
				 dim_gmres,
				 x_data,
				 b_data,
				 work,
				 amat,
				 precon == Py_None ? NULL : precon);

  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

static char qmrs_doc[] = 
"info, iter, relres = qmrs(A, b, x, tol, maxit, K)\n\
\n\
Minimal Residual method.";

static PyObject *
ItSolvers_qmrs(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  int dim[2];			/* shape of amat */
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|O",
			&amat, 
			&b,
			&x,
			&tol,
			&maxit,
			&precon))
    return NULL;

  /* check shape of matrix object */
  SpMatrix_GetShape((PyObject *)amat, dim);
  if (dim[0] != dim[1] || dim[0] <= 0) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return NULL;
  }
  n = dim[0];

  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, 6*n);

  /* call qmrs routine */
  info = Itsolvers_qmrs_kernel(n, 
			       b_data, 
			       x_data, 
			       work, 
			       tol, 
			       maxit, 
			       &iter, 
			       &relres,
			       amat,
			       precon == Py_None ? NULL : precon);
  
  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

static char cgs_doc[] = 
"info, iter, relres = cgs(A, b, x, tol, maxit, K)\n\
\n\
Conjugate Gradient Square method.";

static PyObject *
ItSolvers_cgs(PyObject *self, PyObject *args)
{
  /* input arguments */
  PyObject *amat;
  PyArrayObject *b;
  PyArrayObject *x;
  double tol;
  int maxit;
  /* output arguments */
  int info;
  int iter;
  double relres;
  /* other variables */
  int n;
  double *work;
  int res;
  int nx, nb;
  double *x_data, *b_data;
  int dim[2];			/* shape of amat */
  PyObject *precon = Py_None;	/* default value for precon */

  /* Parse input arguments */
  if (!PyArg_ParseTuple(args, "OOOdi|O",
			&amat, 
			&b,
			&x,
			&tol,
			&maxit,
			&precon))
    return NULL;

  /* check shape of matrix object */
  SpMatrix_GetShape((PyObject *)amat, dim);
  if (dim[0] != dim[1] || dim[0] <= 0) {
    PyErr_SetString(PyExc_ValueError, "invalid matrix shape");
    return NULL;
  }
  n = dim[0];

  /* Make sure that x and b are continous double arrays */
  res = PyArray_As1D((PyObject **)&x, (char **)&x_data, &nx, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert x to double array");
    return NULL;
  }
  res = PyArray_As1D((PyObject **)&b, (char **)&b_data, &nb, PyArray_DOUBLE);
  if (res == -1) {
    PyErr_SetString(PyExc_ValueError, "Unable to convert b to double array");
    return NULL;
  }

  if (nx != nb || n != nx) {
    PyErr_SetString(PyExc_ValueError, "incompatible operand shapes");
    return NULL;
  }

  /* allocate workspace for temporary vectors */
  work = PyMem_New(double, 8*n);

  /* call cgs routine */
  info = Itsolvers_cgs_kernel(n, 
			      b_data, 
			      x_data, 
			      maxit, 
			      tol,
			      work,
			      &iter,
			      &relres,
			      amat,
			      precon == Py_None ? NULL : precon);
  
  /* free workspace */
  PyMem_DEL(work);
  res = PyArray_Free((PyObject *)x, (char *)x_data);
  assert(res != -1);
  res = PyArray_Free((PyObject *)b, (char *)b_data);
  assert(res != -1);

  /* return result tuple */
  if (PyErr_Occurred())
    return NULL;
  else
    return Py_BuildValue("iid", info, iter, relres);
}

/** table of module functions
 */
static PyMethodDef itsolvers_methods[] = {
  {"pcg",     (PyCFunction)ItSolvers_pcg,      METH_VARARGS, pcg_doc},
  {"bicgstab",(PyCFunction)ItSolvers_bicgstab, METH_VARARGS, bicgstab_doc},
  {"minres",  (PyCFunction)ItSolvers_minres,   METH_VARARGS, minres_doc},
  {"gmres",   (PyCFunction)ItSolvers_gmres,    METH_VARARGS, gmres_doc},
  {"qmrs",    (PyCFunction)ItSolvers_qmrs,     METH_VARARGS, qmrs_doc},
  {"cgs",     (PyCFunction)ItSolvers_cgs,      METH_VARARGS, cgs_doc}, 
  {NULL, NULL}	/* sentinel */
};

static char module_doc[] = "This module ...\n\
\n\
All iterative solvers in this module provide the following interface\n\
\n\
    info, iter, relres = pcg(A, b, x, tol, maxit, K)\n\
\n\
parameters\n\
----------\n\
\n\
A     system matrix\n\
b     right hand side\n\
x     initial guess on input, solution vector on output\n\
tol   requested error tolerance\n\
maxit maximum number of iteration to be executed\n\
K     preconditioner\n\
      (optional parameter)\n\
\n\
the iterative solvers defined in this module may accept additional\n\
parameters, which are passed as keyword arguments.\n\
\n\
return value\n\
------------\n\
\n\
The result is tuple with 3 elements:\n\
\n\
 info    return code with the following meaning\n\
\n\
          2  iteration converged, residual is as small as seems reasonable on this machine.\n\
\n\
          1  iteration converged, b = 0,  so the exact solution is  x = 0.\n\
\n\
          0  iteration converged, relative error appears to be less than tol\n\
\n\
         -1  iteration not converged, maximum number of iterations was reached\n\
\n\
         -2  iteration not converged, the system involving the preconditioner was ill conditioned\n\
\n\
         -3  iteration not converged, an inner product of the form  x(t)*K^(-1)*x\n\
             was not positive, so the preconditioning matrix K does not appear to be positive \n\
             definite.\n\
\n\
         -4  iteration not converged, the matrix A appears to be very ill-conditioned\n\
\n\
         -5  iteration not converged, the method stagnated\n\
\n\
         -6  iteration not converged, a scalar quantity became too small or too large to \n\
             continue computing\n\
\n\
         Note that not all iterative solvers check for all above error conditions.\n\
\n\
 iter    number of iterations executed\n\
\n\
 relres  relative error of the solution\n\
\n\
";

PyMODINIT_FUNC
inititsolvers(void)
{
  PyObject *m;
  
  m = Py_InitModule3("itsolvers", itsolvers_methods, module_doc);

  /* initialize scipy array module */
  import_array();
  /* initialize spmatrix module */
  import_spmatrix();

  /* No need to check the error here, the caller will do that */
}
