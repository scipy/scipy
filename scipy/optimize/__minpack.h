/* This file is used to make _multipackmodule.c */
/* $Revision$ */
/* module_methods:
  {"_hybrd", minpack_hybrd, METH_VARARGS, doc_hybrd},
  {"_hybrj", minpack_hybrj, METH_VARARGS, doc_hybrj},
  {"_lmdif", minpack_lmdif, METH_VARARGS, doc_lmdif},
  {"_lmder", minpack_lmder, METH_VARARGS, doc_lmder},
  {"_chkder", minpack_chkder, METH_VARARGS, doc_chkder},
 */
#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"
#include "ccallback.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}

#define STORE_VARS() ccallback_t callback; int callback_inited = 0; jac_callback_info_t jac_callback_info;
#define STORE_VARS_NO_INFO() ccallback_t callback; int callback_inited = 0;

#define INIT_FUNC(fun,arg,errobj) do { /* Get extra arguments or set to zero length tuple */ \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail_free; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun)) \
    PYERR(errobj,"First argument must be a callable function."); \
  if (init_callback(&callback, fun, arg) != 0) \
    PYERR(errobj,"Could not init callback");\
  callback_inited = 1; \
  } while(0)

#define INIT_JAC_FUNC(fun,Dfun,arg,col_deriv,errobj) do { \
  if (arg == NULL) { \
    if ((arg = PyTuple_New(0)) == NULL) goto fail_free; \
  } \
  else \
    Py_INCREF(arg);   /* We decrement on exit. */ \
  if (!PyTuple_Check(arg))  \
    PYERR(errobj,"Extra Arguments must be in a tuple"); \
  /* Set up callback functions */ \
  if (!PyCallable_Check(fun) || (Dfun != Py_None && !PyCallable_Check(Dfun))) \
    PYERR(errobj,"The function and its Jacobian must be callable functions."); \
  if (init_jac_callback(&callback, &jac_callback_info, fun, Dfun, arg, col_deriv) != 0) \
    PYERR(errobj,"Could not init callback");\
  callback_inited = 1; \
} while(0)

#define RESTORE_JAC_FUNC() do { \
  if (callback_inited && release_callback(&callback) != 0) { \
    goto fail_free; \
  }} while(0)

#define RESTORE_FUNC() do { \
  if (callback_inited && release_callback(&callback) != 0) { \
    goto fail_free; \
  }} while(0)

#define SET_DIAG(ap_diag,o_diag,mode) { /* Set the diag vector from input */ \
  if (o_diag == NULL || o_diag == Py_None) { \
    ap_diag = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)PyArray_DATA(ap_diag); \
    mode = 1; \
  } \
  else { \
    ap_diag = (PyArrayObject *)PyArray_ContiguousFromObject(o_diag, NPY_DOUBLE, 1, 1); \
    if (ap_diag == NULL) goto fail; \
    diag = (double *)PyArray_DATA(ap_diag); \
    mode = 2; } }

#define MATRIXC2F(jac,data,m,n) {double *p1=(double *)(jac), *p2, *p3=(double *)(data);\
int i,j;\
for (j=0;j<(m);p3++,j++) \
  for (p2=p3,i=0;i<(n);p2+=(m),i++,p1++) \
    *p1 = *p2; }

static PyObject *minpack_error;

typedef struct {
  PyObject *Dfun;
  PyObject *extra_args;
  int jac_transpose;
} jac_callback_info_t;

static PyObject *call_python_function(PyObject *func, npy_intp n, double *x, PyObject *args, int dim, PyObject *error_obj, npy_intp out_size)
{
  /*
    This is a generic function to call a python function that takes a 1-D
    sequence as a first argument and optional extra_arguments (should be a
    zero-length tuple if none desired).  The result of the function is
    returned in a multiarray object.
        -- build sequence object from values in x.
	-- add extra arguments (if any) to an argument list.
	-- call Python callable object
 	-- check if error occurred:
	         if so return NULL
	-- if no error, place result of Python code into multiarray object.
  */

  PyArrayObject *sequence = NULL;
  PyObject *arglist = NULL;
  PyObject *arg1 = NULL;
  PyObject *result = NULL;
  PyArrayObject *result_array = NULL;
  npy_intp fvec_sz = 0;

  /* Build sequence argument from inputs */
  sequence = (PyArrayObject *)PyArray_SimpleNewFromData(1, &n, NPY_DOUBLE, (char *)x);
  if (sequence == NULL) PYERR2(error_obj,"Internal failure to make an array of doubles out of first\n                 argument to function call.");

  /* Build argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    Py_DECREF(sequence);
    return NULL;
  }
  PyTuple_SET_ITEM(arg1, 0, (PyObject *)sequence);
                /* arg1 now owns sequence reference */
  if ((arglist = PySequence_Concat( arg1, args)) == NULL)
    PYERR2(error_obj,"Internal error constructing argument list.");

  Py_DECREF(arg1);    /* arglist has a reference to sequence, now. */
  arg1 = NULL;

  /* Call function object --- variable passed to routine.  Extra
          arguments are in another passed variable.
   */
  if ((result = PyObject_CallObject(func, arglist))==NULL) {
      goto fail;
  }

  if ((result_array = (PyArrayObject *)PyArray_ContiguousFromObject(result, NPY_DOUBLE, dim-1, dim))==NULL)
    PYERR2(error_obj,"Result from function call is not a proper array of floats.");

  fvec_sz = PyArray_SIZE(result_array);
  if(out_size != -1 && fvec_sz != out_size){
      PyErr_SetString(PyExc_ValueError, "The array returned by a function changed size between calls");
      Py_DECREF(result_array);
      goto fail;
  }

  Py_DECREF(result);
  Py_DECREF(arglist);
  return (PyObject *)result_array;

 fail:
  Py_XDECREF(arglist);
  Py_XDECREF(result);
  Py_XDECREF(arg1);
  return NULL;
}

void CHKDER(const int,const int,double*,const double*,const double*,const int,double*,const double*,const int,double*);
void HYBRD(int(*)(int*,double*,double*,int*),const int,double*,double*,const double,const int,const int,const int,const double,double*,const int,const double,const int,int*,int*,double*,const int,double*,const int,double*,double*,double*,double*,double*);
void HYBRJ(int(*)(int*,double*,double*,double*,int*,int*),const int,double*,double*,double*,const int,const double,const int,double*,const int,const double,const int,int*,int*,int*,double*,const int,double*,double*,double*,double*,double*);
void LMDIF(int(*)(int*,int*,double*,double*,int*),const int,const int,double*,double*,const double,const double,const double,const int,const double,double*,const int,const double,const int,int*,int*,double*,const int,int*,double*,double*,double*,double*,double*);
void LMDER(int(*)(int*,int*,double*,double*,double*,int*,int*),const int,const int,double*,double*,double*,const int,const double,const double,const double,const int,double*,const int,const double,const int,int*,int*,int*,int*,double*,double*,double*,double*,double*);
void LMSTR(int(*)(int*,int*,double*,double*,double*,int*),const int,const int,double*,double*,double*,const int,const double,const double,const double,const int,double*,const int,const double,const int,int*,int*,int*,int*,double*,double*,double*,double*,double*);

/* We only use ccallback with Python functions right now */
static ccallback_signature_t call_signatures[] = {
  {NULL}
};

static int init_callback(ccallback_t *callback, PyObject *fcn, PyObject *extra_args)
{
  int ret;
  int flags = CCALLBACK_OBTAIN;

  ret = ccallback_prepare(callback, call_signatures, fcn, flags);
  if (ret == -1) {
    return -1;
  }

  callback->info_p = (void *)extra_args;

  return 0;
}

static int release_callback(ccallback_t *callback)
{
  return ccallback_release(callback) != 0;
}

static int init_jac_callback(ccallback_t *callback, jac_callback_info_t *jac_callback_info, PyObject *fcn, PyObject *Dfun, PyObject *extra_args, int col_deriv)
{
  int ret;
  int flags = CCALLBACK_OBTAIN;

  ret = ccallback_prepare(callback, call_signatures, fcn, flags);
  if (ret == -1) {
    return -1;
  }

  jac_callback_info->Dfun = Dfun;
  jac_callback_info->extra_args = extra_args;
  jac_callback_info->jac_transpose = !col_deriv;

  callback->info_p = (void *)jac_callback_info;

  return 0;
}


int raw_multipack_calling_function(int *n, double *x, double *fvec, int *iflag)
{
  /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any
	-- otherwise place result of calculation in *fvec
  */

  ccallback_t *callback = ccallback_obtain();
  PyObject *multipack_python_function = callback->py_function;
  PyObject *multipack_extra_arguments = (PyObject *)callback->info_p;

  PyArrayObject *result_array = NULL;

  result_array = (PyArrayObject *)call_python_function(multipack_python_function, *n, x, multipack_extra_arguments, 1, minpack_error, *n);
  if (result_array == NULL) {
    *iflag = -1;
    return -1;
  }
  memcpy(fvec, PyArray_DATA(result_array), (*n)*sizeof(double));
  Py_DECREF(result_array);
  return 0;

}


int jac_multipack_calling_function(int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag)
{
  /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any
	-- otherwise place result of calculation in *fvec or *fjac.

     If iflag = 1 this should compute the function.
     If iflag = 2 this should compute the jacobian (derivative matrix)
  */

  ccallback_t *callback = ccallback_obtain();
  PyObject *multipack_python_function = callback->py_function,
           *multipack_python_jacobian = ((jac_callback_info_t *)callback->info_p)->Dfun;
  PyObject *multipack_extra_arguments = ((jac_callback_info_t *)callback->info_p)->extra_args;
  int multipack_jac_transpose = ((jac_callback_info_t *)callback->info_p)->jac_transpose;

  PyArrayObject *result_array;

  if (*iflag == 1) {
    result_array = (PyArrayObject *)call_python_function(multipack_python_function, *n, x, multipack_extra_arguments, 1, minpack_error, *n);
    if (result_array == NULL) {
      *iflag = -1;
      return -1;
    }
    memcpy(fvec, PyArray_DATA(result_array), (*n)*sizeof(double));
  }
  else {         /* iflag == 2 */
    result_array = (PyArrayObject *)call_python_function(multipack_python_jacobian, *n, x, multipack_extra_arguments, 2, minpack_error, (*n)*(*ldfjac));
    if (result_array == NULL) {
      *iflag = -1;
      return -1;
    }
    if (multipack_jac_transpose == 1)
      MATRIXC2F(fjac, PyArray_DATA(result_array), *n, *ldfjac)
    else
      memcpy(fjac, PyArray_DATA(result_array), (*n)*(*ldfjac)*sizeof(double));
  }

  Py_DECREF(result_array);
  return 0;
}

int raw_multipack_lm_function(int *m, int *n, double *x, double *fvec, int *iflag)
{
  /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any
	-- otherwise place result of calculation in *fvec
  */

  ccallback_t *callback = ccallback_obtain();
  PyObject *multipack_python_function = callback->py_function;
  PyObject *multipack_extra_arguments = (PyObject *)callback->info_p;

  PyArrayObject *result_array = NULL;

  result_array = (PyArrayObject *)call_python_function(multipack_python_function,*n, x, multipack_extra_arguments, 1, minpack_error, *m);
  if (result_array == NULL) {
    *iflag = -1;
    return -1;
  }
  memcpy(fvec, PyArray_DATA(result_array), (*m)*sizeof(double));
  Py_DECREF(result_array);
  return 0;
}

int jac_multipack_lm_function(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag)
{
  /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any
	-- otherwise place result of calculation in *fvec or *fjac.

     If iflag = 1 this should compute the function.
     If iflag = 2 this should compute the jacobian (derivative matrix)
  */

  ccallback_t *callback = ccallback_obtain();
  PyObject *multipack_python_function = callback->py_function,
           *multipack_python_jacobian = ((jac_callback_info_t *)callback->info_p)->Dfun;
  PyObject *multipack_extra_arguments = ((jac_callback_info_t *)callback->info_p)->extra_args;
  int multipack_jac_transpose = ((jac_callback_info_t *)callback->info_p)->jac_transpose;

  PyArrayObject *result_array;

  if (*iflag == 1) {
    result_array = (PyArrayObject *)call_python_function(multipack_python_function, *n, x, multipack_extra_arguments, 1, minpack_error, *m);
    if (result_array == NULL) {
      *iflag = -1;
      return -1;
    }
    memcpy(fvec, PyArray_DATA(result_array), (*m)*sizeof(double));
  }
  else {         /* iflag == 2 */
    result_array = (PyArrayObject *)call_python_function(multipack_python_jacobian, *n, x, multipack_extra_arguments, 2, minpack_error, (*n)*(*ldfjac));
    if (result_array == NULL) {
      *iflag = -1;
      return -1;
    }
    if (multipack_jac_transpose == 1)
      MATRIXC2F(fjac, PyArray_DATA(result_array), *n, *ldfjac)
    else
      memcpy(fjac, PyArray_DATA(result_array), (*n)*(*ldfjac)*sizeof(double));
  }

  Py_DECREF(result_array);
  return 0;
}


static char doc_hybrd[] = "[x,infodict,info] = _hybrd(fun, x0, args, full_output, xtol, maxfev, ml, mu, epsfcn, factor, diag)";

static PyObject *minpack_hybrd(PyObject *dummy, PyObject *args) {
  PyObject *fcn, *x0, *extra_args = NULL, *o_diag = NULL;
  int      full_output = 0, maxfev = -10, ml = -10, mu = -10;
  double   xtol = 1.49012e-8, epsfcn = 0.0, factor = 1.0e2;
  int      mode = 2, nprint = 0, info, nfev, ldfjac;
  npy_intp n,lr;
  int      n_int, lr_int;  /* for casted storage to pass int into HYBRD */
  double   *x, *fvec, *diag, *fjac, *r, *qtf;

  PyArrayObject *ap_x = NULL, *ap_fvec = NULL;
  PyArrayObject *ap_fjac = NULL, *ap_r = NULL, *ap_qtf = NULL;
  PyArrayObject *ap_diag = NULL;

  npy_intp dims[2];
  int      allocated = 0;
  double   *wa = NULL;

  STORE_VARS_NO_INFO();    /* Define storage variables for global variables. */

  if (!PyArg_ParseTuple(args, "OO|OidiiiddO", &fcn, &x0, &extra_args, &full_output, &xtol, &maxfev, &ml, &mu, &epsfcn, &factor, &o_diag)) return NULL;

  INIT_FUNC(fcn,extra_args,minpack_error);

  /* Initial input vector */
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x0, NPY_DOUBLE, 1, 1);
  if (ap_x == NULL) goto fail;
  x = (double *) PyArray_DATA(ap_x);
  n = PyArray_DIMS(ap_x)[0];

  lr = n * (n + 1) / 2;
  if (ml < 0) ml = n-1;
  if (mu < 0) mu = n-1;
  if (maxfev < 0) maxfev = 200*(n+1);

  /* Setup array to hold the function evaluations */
  ap_fvec = (PyArrayObject *)call_python_function(fcn, n, x, extra_args, 1, minpack_error, -1);
  if (ap_fvec == NULL) goto fail;
  fvec = (double *) PyArray_DATA(ap_fvec);
  if (PyArray_NDIM(ap_fvec) == 0)
    n = 1;
  else if (PyArray_DIMS(ap_fvec)[0] < n)
    n = PyArray_DIMS(ap_fvec)[0];

  SET_DIAG(ap_diag,o_diag,mode);

  dims[0] = n; dims[1] = n;
  ap_r = (PyArrayObject *)PyArray_SimpleNew(1,&lr,NPY_DOUBLE);
  ap_qtf = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE);
  ap_fjac = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  if (ap_r == NULL || ap_qtf == NULL || ap_fjac ==NULL) goto fail;

  r = (double *) PyArray_DATA(ap_r);
  qtf = (double *) PyArray_DATA(ap_qtf);
  fjac = (double *) PyArray_DATA(ap_fjac);
  ldfjac = dims[1];

  if ((wa = malloc(4*n * sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  allocated = 1;

  /* Call the underlying FORTRAN routines. */
  n_int = n; lr_int = lr; /* cast/store/pass into HYBRD */
  HYBRD(raw_multipack_calling_function, n_int, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode, factor, nprint, &info, &nfev, fjac, ldfjac, r, lr_int, qtf, wa, wa+n, wa+2*n, wa+3*n);

  RESTORE_FUNC();

  if (info < 0) goto fail;            /* Python Terminated */


  free(wa);
  Py_DECREF(extra_args);
  Py_DECREF(ap_diag);

  if (full_output) {
    return Py_BuildValue("N{s:N,s:i,s:N,s:N,s:N}i",PyArray_Return(ap_x),"fvec",PyArray_Return(ap_fvec),"nfev",nfev,"fjac",PyArray_Return(ap_fjac),"r",PyArray_Return(ap_r),"qtf",PyArray_Return(ap_qtf),info);
  }
  else {
    Py_DECREF(ap_fvec);
    Py_DECREF(ap_fjac);
    Py_DECREF(ap_r);
    Py_DECREF(ap_qtf);
    return Py_BuildValue("Ni",PyArray_Return(ap_x),info);
  }

 fail:
  RESTORE_FUNC();
 fail_free:
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_fvec);
  Py_XDECREF(ap_diag);
  Py_XDECREF(ap_fjac);
  Py_XDECREF(ap_r);
  Py_XDECREF(ap_qtf);
  if (allocated) free(wa);
  return NULL;
}


static char doc_hybrj[] = "[x,infodict,info] = _hybrj(fun, Dfun, x0, args, full_output, col_deriv, xtol, maxfev, factor, diag)";

static PyObject *minpack_hybrj(PyObject *dummy, PyObject *args) {
  PyObject *fcn, *Dfun, *x0, *extra_args = NULL, *o_diag = NULL;
  int      full_output = 0, maxfev = -10, col_deriv = 1;
  double   xtol = 1.49012e-8, factor = 1.0e2;
  int      mode = 2, nprint = 0, info, nfev, njev, ldfjac;
  npy_intp n, lr;
  int n_int, lr_int;
  double   *x, *fvec, *diag, *fjac, *r, *qtf;

  PyArrayObject *ap_x = NULL, *ap_fvec = NULL;
  PyArrayObject *ap_fjac = NULL, *ap_r = NULL, *ap_qtf = NULL;
  PyArrayObject *ap_diag = NULL;

  npy_intp dims[2];
  int      allocated = 0;
  double   *wa = NULL;

  STORE_VARS();

  if (!PyArg_ParseTuple(args, "OOO|OiididO", &fcn, &Dfun, &x0, &extra_args, &full_output, &col_deriv, &xtol, &maxfev, &factor, &o_diag)) return NULL;

  INIT_JAC_FUNC(fcn,Dfun,extra_args,col_deriv,minpack_error);

  /* Initial input vector */
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x0, NPY_DOUBLE, 1, 1);
  if (ap_x == NULL) goto fail;
  x = (double *) PyArray_DATA(ap_x);
  n = PyArray_DIMS(ap_x)[0];
  lr = n * (n + 1) / 2;

  if (maxfev < 0) maxfev = 100*(n+1);

  /* Setup array to hold the function evaluations */
  ap_fvec = (PyArrayObject *)call_python_function(fcn, n, x, extra_args, 1, minpack_error, -1);
  if (ap_fvec == NULL) goto fail;
  fvec = (double *) PyArray_DATA(ap_fvec);
  if (PyArray_NDIM(ap_fvec) == 0)
    n = 1;
  else if (PyArray_DIMS(ap_fvec)[0] < n)
    n = PyArray_DIMS(ap_fvec)[0];

  SET_DIAG(ap_diag,o_diag,mode);

  dims[0] = n; dims[1] = n;
  ap_r = (PyArrayObject *)PyArray_SimpleNew(1,&lr,NPY_DOUBLE);
  ap_qtf = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE);
  ap_fjac = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  if (ap_r == NULL || ap_qtf == NULL || ap_fjac ==NULL) goto fail;

  r = (double *) PyArray_DATA(ap_r);
  qtf = (double *) PyArray_DATA(ap_qtf);
  fjac = (double *) PyArray_DATA(ap_fjac);

  ldfjac = dims[1];

  if ((wa = malloc(4*n * sizeof(double)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  allocated = 1;

  /* Call the underlying FORTRAN routines. */
  n_int = n; lr_int = lr; /* cast/store/pass into HYBRJ */
  HYBRJ(jac_multipack_calling_function, n_int, x, fvec, fjac, ldfjac, xtol, maxfev, diag, mode, factor, nprint, &info, &nfev, &njev, r, lr_int, qtf, wa, wa+n, wa+2*n, wa+3*n);

  RESTORE_JAC_FUNC();

  if (info < 0) goto fail;            /* Python Terminated */

  free(wa);
  Py_DECREF(extra_args);
  Py_DECREF(ap_diag);

  if (full_output) {
    return Py_BuildValue("N{s:N,s:i,s:i,s:N,s:N,s:N}i",PyArray_Return(ap_x),"fvec",PyArray_Return(ap_fvec),"nfev",nfev,"njev",njev,"fjac",PyArray_Return(ap_fjac),"r",PyArray_Return(ap_r),"qtf",PyArray_Return(ap_qtf),info);
  }
  else {
    Py_DECREF(ap_fvec);
    Py_DECREF(ap_fjac);
    Py_DECREF(ap_r);
    Py_DECREF(ap_qtf);
    return Py_BuildValue("Ni",PyArray_Return(ap_x),info);
  }

 fail:
  RESTORE_JAC_FUNC();
 fail_free:
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_fvec);
  Py_XDECREF(ap_fjac);
  Py_XDECREF(ap_diag);
  Py_XDECREF(ap_r);
  Py_XDECREF(ap_qtf);
  if (allocated) free(wa);
  return NULL;

}

/************************ Levenberg-Marquardt *******************/

static char doc_lmdif[] = "[x,infodict,info] = _lmdif(fun, x0, args, full_output, ftol, xtol, gtol, maxfev, epsfcn, factor, diag)";

static PyObject *minpack_lmdif(PyObject *dummy, PyObject *args) {
  PyObject *fcn, *x0, *extra_args = NULL, *o_diag = NULL;
  int      full_output = 0, maxfev = -10;
  double   xtol = 1.49012e-8, ftol = 1.49012e-8;
  double   gtol = 0.0, epsfcn = 0.0, factor = 1.0e2;
  int      m, mode = 2, nprint = 0, info = 0, nfev, ldfjac, *ipvt;
  npy_intp n;
  int      n_int;  /* for casted storage to pass int into LMDIF */
  double   *x, *fvec, *diag, *fjac, *qtf;

  PyArrayObject *ap_x = NULL, *ap_fvec = NULL;
  PyArrayObject *ap_fjac = NULL, *ap_ipvt = NULL, *ap_qtf = NULL;
  PyArrayObject *ap_diag = NULL;

  npy_intp dims[2];
  int      allocated = 0;
  double   *wa = NULL;

  STORE_VARS_NO_INFO();

  if (!PyArg_ParseTuple(args, "OO|OidddiddO", &fcn, &x0, &extra_args, &full_output, &ftol, &xtol, &gtol, &maxfev, &epsfcn, &factor, &o_diag)) return NULL;

  INIT_FUNC(fcn,extra_args,minpack_error);

  /* Initial input vector */
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x0, NPY_DOUBLE, 1, 1);
  if (ap_x == NULL) goto fail;
  x = (double *) PyArray_DATA(ap_x);
  n = PyArray_DIMS(ap_x)[0];
  dims[0] = n;

  SET_DIAG(ap_diag,o_diag,mode);

  if (maxfev < 0) maxfev = 200*(n+1);

  /* Setup array to hold the function evaluations and find it's size*/
  ap_fvec = (PyArrayObject *)call_python_function(fcn, n, x, extra_args, 1, minpack_error, -1);
  if (ap_fvec == NULL) goto fail;
  fvec = (double *) PyArray_DATA(ap_fvec);
  m = (PyArray_NDIM(ap_fvec) > 0 ? PyArray_DIMS(ap_fvec)[0] : 1);

  dims[0] = n; dims[1] = m;
  ap_ipvt = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_INT);
  ap_qtf = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE);
  ap_fjac = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  if (ap_ipvt == NULL || ap_qtf == NULL || ap_fjac ==NULL) goto fail;

  ipvt = (int *) PyArray_DATA(ap_ipvt);
  qtf = (double *) PyArray_DATA(ap_qtf);
  fjac = (double *) PyArray_DATA(ap_fjac);
  ldfjac = dims[1];
  wa = (double *)malloc((3*n + m)* sizeof(double));
  if (wa == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  allocated = 1;

  /* Call the underlying FORTRAN routines. */
  n_int = n; /* to provide int*-pointed storage for int argument of LMDIF */
  LMDIF(raw_multipack_lm_function, m, n_int, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, diag, mode, factor, nprint, &info, &nfev, fjac, ldfjac, ipvt, qtf, wa, wa+n, wa+2*n, wa+3*n);

  RESTORE_FUNC();

  if (info < 0) goto fail;           /* Python error */

  free(wa);
  Py_DECREF(extra_args);
  Py_DECREF(ap_diag);

  if (full_output) {
    return Py_BuildValue("N{s:N,s:i,s:N,s:N,s:N}i",PyArray_Return(ap_x),"fvec",PyArray_Return(ap_fvec),"nfev",nfev,"fjac",PyArray_Return(ap_fjac),"ipvt",PyArray_Return(ap_ipvt),"qtf",PyArray_Return(ap_qtf),info);
  }
  else {
    Py_DECREF(ap_fvec);
    Py_DECREF(ap_fjac);
    Py_DECREF(ap_ipvt);
    Py_DECREF(ap_qtf);
    return Py_BuildValue("Ni",PyArray_Return(ap_x),info);
  }

 fail:
  RESTORE_FUNC();
 fail_free:
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_fvec);
  Py_XDECREF(ap_fjac);
  Py_XDECREF(ap_diag);
  Py_XDECREF(ap_ipvt);
  Py_XDECREF(ap_qtf);
  if (allocated) free(wa);
  return NULL;
}


static char doc_lmder[] = "[x,infodict,info] = _lmder(fun, Dfun, x0, args, full_output, col_deriv, ftol, xtol, gtol, maxfev, factor, diag)";

static PyObject *minpack_lmder(PyObject *dummy, PyObject *args) {
  PyObject *fcn, *x0, *Dfun, *extra_args = NULL, *o_diag = NULL;
  int      full_output = 0, maxfev = -10, col_deriv = 1;
  double   xtol = 1.49012e-8, ftol = 1.49012e-8;
  double   gtol = 0.0, factor = 1.0e2;
  int      m, mode = 2, nprint = 0, info, nfev, njev, ldfjac, *ipvt;
  npy_intp n;
  int n_int;
  double   *x, *fvec, *diag, *fjac, *qtf;

  PyArrayObject *ap_x = NULL, *ap_fvec = NULL;
  PyArrayObject *ap_fjac = NULL, *ap_ipvt = NULL, *ap_qtf = NULL;
  PyArrayObject *ap_diag = NULL;

  npy_intp dims[2];
  int      allocated = 0;
  double   *wa = NULL;

  STORE_VARS();

  if (!PyArg_ParseTuple(args, "OOO|OiidddidO", &fcn, &Dfun, &x0, &extra_args, &full_output, &col_deriv, &ftol, &xtol, &gtol, &maxfev, &factor, &o_diag)) return NULL;

  INIT_JAC_FUNC(fcn,Dfun,extra_args,col_deriv,minpack_error);

  /* Initial input vector */
  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(x0, NPY_DOUBLE, 1, 1);
  if (ap_x == NULL) goto fail;
  x = (double *) PyArray_DATA(ap_x);
  n = PyArray_DIMS(ap_x)[0];

  if (maxfev < 0) maxfev = 100*(n+1);

  /* Setup array to hold the function evaluations */
  ap_fvec = (PyArrayObject *)call_python_function(fcn, n, x, extra_args, 1, minpack_error, -1);
  if (ap_fvec == NULL) goto fail;
  fvec = (double *) PyArray_DATA(ap_fvec);

  SET_DIAG(ap_diag,o_diag,mode);

  m = (PyArray_NDIM(ap_fvec) > 0 ? PyArray_DIMS(ap_fvec)[0] : 1);

  dims[0] = n; dims[1] = m;
  ap_ipvt = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_INT);
  ap_qtf = (PyArrayObject *)PyArray_SimpleNew(1,&n,NPY_DOUBLE);
  ap_fjac = (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);

  if (ap_ipvt == NULL || ap_qtf == NULL || ap_fjac ==NULL) goto fail;

  ipvt = (int *) PyArray_DATA(ap_ipvt);
  qtf = (double *) PyArray_DATA(ap_qtf);
  fjac = (double *) PyArray_DATA(ap_fjac);
  ldfjac = dims[1];
  wa = (double *)malloc((3*n + m)* sizeof(double));
  if (wa == NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  allocated = 1;

  /* Call the underlying FORTRAN routines. */
  n_int = n;
  LMDER(jac_multipack_lm_function, m, n_int, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, diag, mode, factor, nprint, &info, &nfev, &njev, ipvt, qtf, wa, wa+n, wa+2*n, wa+3*n);

  RESTORE_JAC_FUNC();

  if (info < 0) goto fail;           /* Python error */

  free(wa);
  Py_DECREF(extra_args);
  Py_DECREF(ap_diag);

  if (full_output) {
    return Py_BuildValue("N{s:N,s:i,s:i,s:N,s:N,s:N}i",PyArray_Return(ap_x),"fvec",PyArray_Return(ap_fvec),"nfev",nfev,"njev",njev,"fjac",PyArray_Return(ap_fjac),"ipvt",PyArray_Return(ap_ipvt),"qtf",PyArray_Return(ap_qtf),info);
  }
  else {
    Py_DECREF(ap_fvec);
    Py_DECREF(ap_fjac);
    Py_DECREF(ap_ipvt);
    Py_DECREF(ap_qtf);
    return Py_BuildValue("Ni",PyArray_Return(ap_x),info);
  }

 fail:
  RESTORE_JAC_FUNC();
 fail_free:
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_x);
  Py_XDECREF(ap_fvec);
  Py_XDECREF(ap_fjac);
  Py_XDECREF(ap_diag);
  Py_XDECREF(ap_ipvt);
  Py_XDECREF(ap_qtf);
  if (allocated) free(wa);
  return NULL;
}


/** Check gradient function **/

static char doc_chkder[] = "_chkder(m,n,x,fvec,fjac,ldfjac,xp,fvecp,mode,err)";

static PyObject *minpack_chkder(PyObject *self, PyObject *args)
{
  PyArrayObject *ap_fvecp = NULL, *ap_fjac = NULL, *ap_err = NULL;
  PyArrayObject *ap_x = NULL, *ap_fvec = NULL, *ap_xp = NULL;
  PyObject *o_x, *o_fvec, *o_fjac, *o_fvecp;
  double *xp, *fvecp, *fjac, *fvec, *x;
  double *err;
  int mode, m, n, ldfjac;

  if (!PyArg_ParseTuple(args,"iiOOOiO!OiO!",&m, &n, &o_x, &o_fvec, &o_fjac, &ldfjac, &PyArray_Type, (PyObject **)&ap_xp, &o_fvecp, &mode, &PyArray_Type, (PyObject **)&ap_err)) return NULL;

  ap_x = (PyArrayObject *)PyArray_ContiguousFromObject(o_x,NPY_DOUBLE,1,1);
  if (ap_x == NULL) goto fail;
  if (n != PyArray_DIMS(ap_x)[0])
     PYERR(minpack_error,"Input data array (x) must have length n");
  x = (double *) PyArray_DATA(ap_x);
  if (!PyArray_IS_C_CONTIGUOUS(ap_xp) || (PyArray_TYPE(ap_xp) != NPY_DOUBLE))
     PYERR(minpack_error,"Seventh argument (xp) must be contiguous array of type Float64.");

  if (mode == 1) {
    fvec = NULL;
    fjac = NULL;
    xp = (double *)PyArray_DATA(ap_xp);
    fvecp = NULL;
    err = NULL;
    CHKDER(m, n, x, fvec, fjac, ldfjac, xp, fvecp, mode, err);
  }
  else if (mode == 2) {
    if (!PyArray_IS_C_CONTIGUOUS(ap_err) || (PyArray_TYPE(ap_err) != NPY_DOUBLE))
       PYERR(minpack_error,"Last argument (err) must be contiguous array of type Float64.");
    ap_fvec = (PyArrayObject *)PyArray_ContiguousFromObject(o_fvec,NPY_DOUBLE,1,1);
    ap_fjac = (PyArrayObject *)PyArray_ContiguousFromObject(o_fjac,NPY_DOUBLE,2,2);
    ap_fvecp = (PyArrayObject *)PyArray_ContiguousFromObject(o_fvecp,NPY_DOUBLE,1,1);
    if (ap_fvec == NULL || ap_fjac == NULL || ap_fvecp == NULL) goto fail;

    fvec = (double *)PyArray_DATA(ap_fvec);
    fjac = (double *)PyArray_DATA(ap_fjac);
    xp = (double *)PyArray_DATA(ap_xp);
    fvecp = (double *)PyArray_DATA(ap_fvecp);
    err = (double *)PyArray_DATA(ap_err);

    CHKDER(m, n, x, fvec, fjac, m, xp, fvecp, mode, err);

    Py_DECREF(ap_fvec);
    Py_DECREF(ap_fjac);
    Py_DECREF(ap_fvecp);
  }
  else
    PYERR(minpack_error,"Invalid mode, must be 1 or 2.");

  Py_DECREF(ap_x);

  Py_INCREF(Py_None);
  return Py_None;

 fail:
  Py_XDECREF(ap_fvec);
  Py_XDECREF(ap_fjac);
  Py_XDECREF(ap_fvecp);
  Py_XDECREF(ap_x);
  return NULL;
}

static struct PyMethodDef minpack_module_methods[] = {
{"_hybrd", minpack_hybrd, METH_VARARGS, doc_hybrd},
{"_hybrj", minpack_hybrj, METH_VARARGS, doc_hybrj},
{"_lmdif", minpack_lmdif, METH_VARARGS, doc_lmdif},
{"_lmder", minpack_lmder, METH_VARARGS, doc_lmder},
{"_chkder", minpack_chkder, METH_VARARGS, doc_chkder},
{NULL,		NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_minpack",
    NULL,
    -1,
    minpack_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__minpack(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }
    minpack_error = PyErr_NewException ("_minpack.error", NULL, NULL);
    if (minpack_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", minpack_error)) {
        return NULL;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

    return module;
}


