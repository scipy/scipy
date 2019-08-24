/* This file should be included into the _multipackmodule file */
/* $Revision$ */
/* module_methods:
  {"_qagse", quadpack_qagse, METH_VARARGS, doc_qagse},
  {"_qagie", quadpack_qagie, METH_VARARGS, doc_qagie},
  {"_qagpe", quadpack_qagpe, METH_VARARGS, doc_qagpe},
  {"_qawoe", quadpack_qawoe, METH_VARARGS, doc_qawoe},
  {"_qawfe", quadpack_qawfe, METH_VARARGS, doc_qawfe},
  {"_qawse", quadpack_qawse, METH_VARARGS, doc_qawse},
  {"_qawce", quadpack_qawce, METH_VARARGS, doc_qawce},
 */
/* link libraries: (should be listed in separate lines)
   quadpack
   linpack_lite
   blas
   mach
 */
/* python files: (to be imported to Multipack.py)
   quadpack.py
 */


#include <Python.h>
#include <setjmp.h>

#include "ccallback.h"

#include "numpy/arrayobject.h"

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

static PyObject *quadpack_error;


#if defined(NO_APPEND_FORTRAN)
  #if defined(UPPERCASE_FORTRAN)
  /* nothing to do here */
  #else
    #define DQAGSE dqagse
    #define DQAGIE dqagie
    #define DQAGPE dqagpe
    #define DQAWOE dqawoe
    #define DQAWFE dqawfe
    #define DQAWSE dqawse
    #define DQAWCE dqawce
  #endif
#else
  #if defined(UPPERCASE_FORTRAN)
    #define DQAGSE DQAGSE_
    #define DQAGIE DQAGIE_
    #define DQAGPE DQAGPE_
    #define DQAWOE DQAWOE_
    #define DQAWFE DQAWFE_
    #define DQAWSE DQAWSE_
    #define DQAWCE DQAWCE_
#else
    #define DQAGSE dqagse_
    #define DQAGIE dqagie_
    #define DQAGPE dqagpe_
    #define DQAWOE dqawoe_
    #define DQAWFE dqawfe_
    #define DQAWSE dqawse_
    #define DQAWCE dqawce_ 
  #endif
#endif

void DQAGSE();
void DQAGIE();
void DQAGPE();
void DQAWOE();
void DQAWFE();
void DQAWSE();
void DQAWCE();



typedef enum {
    CB_1D_USER = 0,
    CB_ND_USER = 1,
    CB_1D = 2,
    CB_ND = 3
} quadpack_signature_t;


static ccallback_signature_t quadpack_call_signatures[] = {
    {"double (double, void *)", CB_1D_USER},
    {"double (int, double *, void *)", CB_ND_USER},
    {"double (double)", CB_1D},
    {"double (int, double *)", CB_ND},
#if NPY_SIZEOF_SHORT == NPY_SIZEOF_INT
    {"double (short, double *)", CB_ND},
    {"double (short, double *, void *)", CB_ND_USER},
#endif
#if NPY_SIZEOF_LONG == NPY_SIZEOF_INT
    {"double (long, double *)", CB_ND},
    {"double (long, double *, void *)", CB_ND_USER},
#endif
    {NULL}
};

static ccallback_signature_t quadpack_call_legacy_signatures[] = {
    {"double (double)", CB_1D},
    {"double (int, double)", CB_ND}, /* sic -- for backward compat only */
#if NPY_SIZEOF_SHORT == NPY_SIZEOF_INT
    {"double (short, double)", CB_ND}, /* sic -- for backward compat only */
#endif
#if NPY_SIZEOF_LONG == NPY_SIZEOF_INT
    {"double (long, double)", CB_ND}, /* sic -- for backward compat only */
#endif
    {NULL}
};


static int
init_multivariate_data(ccallback_t *callback, int ndim, PyObject *extra_arguments)
{
    double *p;
    Py_ssize_t i, size;

    callback->info_p = NULL;

    p = (double *)malloc(sizeof(double) * ndim);
    if (p == NULL) {
        free(p);
        PyErr_SetString(PyExc_MemoryError, "failed to allocate memory");
        return -1;
    }

    size = PyTuple_Size(extra_arguments);
    if (size != ndim - 1) {
        free(p);
        PyErr_SetString(PyExc_ValueError, "extra arguments don't match ndim");
        return -1;
    }

    p[0] = 0;
    for (i = 0; i < size; ++i) {
        PyObject *item;

        item = PyTuple_GET_ITEM(extra_arguments, i);
        p[i+1] = PyFloat_AsDouble(item);
        if (PyErr_Occurred()) {
            free(p);
            return -1;
        }
    }

    callback->info_p = (void *)p;
    return 0;
}


static int
init_callback(ccallback_t *callback, PyObject *func, PyObject *extra_arguments)
{
    static PyObject *cfuncptr_type = NULL;

    int ret;
    int ndim;
    int flags = CCALLBACK_OBTAIN;
    int legacy = 0;
    ccallback_signature_t *signatures = quadpack_call_signatures;

    if (cfuncptr_type == NULL) {
        PyObject *module;

        module = PyImport_ImportModule("ctypes");
        if (module == NULL) {
            return -1;
        }

        cfuncptr_type = PyObject_GetAttrString(module, "_CFuncPtr");
        Py_DECREF(module);
        if (cfuncptr_type == NULL) {
            return -1;
        }
    }

    if (PyObject_TypeCheck(func, (PyTypeObject *)cfuncptr_type)) {
        /* Legacy support --- ctypes objects can be passed in as-is */
        flags |= CCALLBACK_PARSE;
        signatures = quadpack_call_legacy_signatures;
        legacy = 1;
    }

    ret = ccallback_prepare(callback, signatures, func, flags);
    if (ret == -1) {
        return -1;
    }

    if (callback->signature == NULL) {
        /* pure-Python */
        callback->info_p = (void *)extra_arguments;
    }
    else if (callback->signature->value == CB_1D || callback->signature->value == CB_1D_USER) {
        /* extra_arguments is just ignored */
        callback->info_p = NULL;
    }
    else {
        if (!PyTuple_Check(extra_arguments)) {
            PyErr_SetString(PyExc_ValueError, "multidimensional integrand but invalid extra args");
            return -1;
        }

        ndim = PyTuple_GET_SIZE(extra_arguments) + 1;

        callback->info = ndim;

        if (init_multivariate_data(callback, ndim, extra_arguments) == -1) {
            return -1;
        }
    }

    return 0;
}


static int
free_callback(ccallback_t *callback)
{
    if (callback->signature && (callback->signature->value == CB_ND ||
                                callback->signature->value == CB_ND_USER)) {
        free(callback->info_p);
        callback->info_p = NULL;
    }

    if (ccallback_release(callback) != 0) {
        return -1;
    }

    return 0;
}


double quad_thunk(double *x)
{
    ccallback_t *callback = ccallback_obtain();
    double result = 0;
    int error = 0;

    if (callback->py_function) {
        PyObject *arg1 = NULL, *argobj = NULL, *arglist = NULL, *res = NULL;
        PyObject *extra_arguments = (PyObject *)callback->info_p;

        argobj = PyFloat_FromDouble(*x);
        if (argobj == NULL) {
            error = 1;
            goto done;
        }

        arg1 = PyTuple_New(1);
        if (arg1 == NULL) {
            error = 1;
            goto done;
        }

        PyTuple_SET_ITEM(arg1, 0, argobj);
        argobj = NULL;

        arglist = PySequence_Concat(arg1, extra_arguments);
        if (arglist == NULL) {
            error = 1;
            goto done;
        }

        res = PyEval_CallObject(callback->py_function, arglist);
        if (res == NULL) {
            error = 1;
            goto done;
        }

        result = PyFloat_AsDouble(res);
        if (PyErr_Occurred()) {
            error = 1;
            goto done;
        }

    done:
        Py_XDECREF(arg1);
        Py_XDECREF(argobj);
        Py_XDECREF(arglist);
        Py_XDECREF(res);
    }
    else {
        switch (callback->signature->value) {
        case CB_1D_USER:
            result = ((double(*)(double, void *))callback->c_function)(*x, callback->user_data);
            break;
        case CB_1D:
            result = ((double(*)(double))callback->c_function)(*x);
            break;
        case CB_ND_USER:
            ((double *)callback->info_p)[0] = *x;
            result = ((double(*)(int, double *, void *))callback->c_function)(
                (int)callback->info, (double *)callback->info_p, callback->user_data);
            break;
        case CB_ND:
            ((double *)callback->info_p)[0] = *x;
            result = ((double(*)(int, double *))callback->c_function)(
                (int)callback->info, (double *)callback->info_p);
            break;
        default:
            error = 1;
            Py_FatalError("scipy.integrate.quad: internal error (this is a bug!): invalid callback type");
            break;
        }
    }

    if (error) {
        longjmp(callback->error_buf, 1);
    }

    return result;
}


static char doc_qagse[] = "[result,abserr,infodict,ier] = _qagse(fun, a, b, | args, full_output, epsabs, epsrel, limit)";

static PyObject *quadpack_qagse(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn;

  int      limit=50;
  npy_intp limit_shape[1];
  int      full_output = 0;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  int      ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Odd|Oiddi", &fcn, &a, &b, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAGSE(quad_thunk, &a, &b, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, 
         blist, rlist, elist, iord, &last);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_iord);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}

static char doc_qagie[] = "[result,abserr,infodict,ier] = _qagie(fun, bound, inf, | args, full_output, epsabs, epsrel, limit)";

static PyObject *quadpack_qagie(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn;

  int      limit=50;
  npy_intp limit_shape[1];
  int      full_output = 0;
  double   bound, epsabs=1.49e-8, epsrel=1.49e-8;
  int      inf, neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Odi|Oiddi", &fcn, &bound, &inf, &extra_args, 
                        &full_output, &epsabs, &epsrel, &limit)) 
    return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL
      || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAGIE(quad_thunk, &bound, &inf, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_iord);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}


static char doc_qagpe[] = "[result,abserr,infodict,ier] = _qagpe(fun, a, b, points, | args, full_output, epsabs, epsrel, limit)";

static PyObject *quadpack_qagpe(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL, *ap_points = NULL;
  PyArrayObject *ap_pts = NULL, *ap_level = NULL;
  PyArrayObject *ap_ndin = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn, *o_points;

  int      limit=50, npts2;
  npy_intp limit_shape[1], npts2_shape[1];
  int      full_output = 0;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  int      *level, *ndin;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  double   *pts, *points;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "OddO|Oiddi", &fcn, &a, &b, &o_points, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  ap_points = (PyArrayObject *)PyArray_ContiguousFromObject(o_points, NPY_DOUBLE, 1, 1);
  if (ap_points == NULL) goto fail;
  npts2 = ap_points->dimensions[0];
  npts2_shape[0] = npts2;
  points = (double *)ap_points->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_pts = (PyArrayObject *)PyArray_SimpleNew(1,npts2_shape,NPY_DOUBLE);
  ap_level = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_ndin = (PyArrayObject *)PyArray_SimpleNew(1,npts2_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL || ap_pts == NULL || ap_level == NULL || ap_ndin == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;
  pts = (double *)ap_pts->data;
  level = (int *)ap_level->data;
  ndin = (int *)ap_level->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAGPE(quad_thunk, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, pts, iord, level, ndin, &last);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  Py_DECREF(ap_points);

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist), "pts", PyArray_Return(ap_pts), "level", PyArray_Return(ap_level), "ndin", PyArray_Return(ap_ndin),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_pts);
    Py_DECREF(ap_iord);
    Py_DECREF(ap_ndin);
    Py_DECREF(ap_level);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  Py_XDECREF(ap_pts);
  Py_XDECREF(ap_points);
  Py_XDECREF(ap_ndin);
  Py_XDECREF(ap_level);
  return NULL;
}


static char doc_qawoe[] = "[result,abserr,infodict,ier] = _qawoe(fun, a, b, omega, integr, | args, full_output, epsabs, epsrel, limit, maxp1, icall, momcom, chebmo)";

static PyObject *quadpack_qawoe(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL, *ap_nnlog = NULL;
  PyArrayObject *ap_chebmo = NULL;

  PyObject *extra_args = NULL, *o_chebmo = NULL;
  PyObject *fcn;

  int      limit=50;
  npy_intp limit_shape[1], sz[2];
  int      full_output = 0, maxp1=50, icall=1;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, integr=1, last=0, momcom=0, *iord;
  int      *nnlog;
  double   result=0.0, abserr=0.0, omega=0.0;
  double   *chebmo;
  double   *alist, *blist, *rlist, *elist;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Odddi|OiddiiiiO", &fcn, &a, &b, &omega, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit, &maxp1, &icall, &momcom, &o_chebmo)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  if (o_chebmo != NULL) {
    ap_chebmo = (PyArrayObject *)PyArray_ContiguousFromObject(o_chebmo, NPY_DOUBLE, 2, 2);
    if (ap_chebmo == NULL) goto fail;
    if (ap_chebmo->dimensions[1] != maxp1 || ap_chebmo->dimensions[0] != 25)
      PYERR(quadpack_error,"Chebyshev moment array has the wrong size.");
  }
  else {
    sz[0] = 25;
    sz[1] = maxp1;
    ap_chebmo = (PyArrayObject *)PyArray_SimpleNew(2,sz,NPY_DOUBLE);
    if (ap_chebmo == NULL) goto fail;
  }
  chebmo = (double *) ap_chebmo->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_nnlog = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_nnlog == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  nnlog = (int *)ap_nnlog->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAWOE(quad_thunk, &a, &b, &omega, &integr, &epsabs, &epsrel, &limit, &icall, &maxp1, &result, &abserr, &neval, &ier, &last, alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N,s:N,s:i,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist), "nnlog", PyArray_Return(ap_nnlog), "momcom", momcom, "chebmo", PyArray_Return(ap_chebmo),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_iord);
    Py_DECREF(ap_nnlog);
    Py_DECREF(ap_chebmo);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  Py_XDECREF(ap_nnlog);
  Py_XDECREF(ap_chebmo);
  return NULL;
}


static char doc_qawfe[] = "[result,abserr,infodict,ier] = _qawfe(fun, a, omega, integr, | args, full_output, epsabs, limlst, limit, maxp1)";

static PyObject *quadpack_qawfe(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL, *ap_nnlog = NULL;
  PyArrayObject *ap_chebmo = NULL, *ap_rslst = NULL;
  PyArrayObject *ap_erlst = NULL, *ap_ierlst = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn;

  int      limlst = 50, limit=50;
  npy_intp limlst_shape[1], limit_shape[1], sz[2];
  int      full_output = 0, maxp1=50;
  double   a, epsabs=1.49e-8;
  int      neval=0, ier=6, integr=1, *iord;
  int      lst, *nnlog, *ierlst;
  double   *chebmo, *rslst, *erlst;
  double   result=0.0, abserr=0.0, omega=0.0;
  double   *alist, *blist, *rlist, *elist;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Oddi|Oidiii", &fcn, &a, &omega, &integr, &extra_args, &full_output, &epsabs, &limlst, &limit, &maxp1)) return NULL;
  limit_shape[0] = limit;
  limlst_shape[0] = limlst;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  sz[0] = 25;
  sz[1] = maxp1;
  ap_chebmo = (PyArrayObject *)PyArray_SimpleNew(2,sz,NPY_DOUBLE);
  if (ap_chebmo == NULL) goto fail;
  chebmo = (double *) ap_chebmo->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_nnlog = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rslst = (PyArrayObject *)PyArray_SimpleNew(1,limlst_shape,NPY_DOUBLE);
  ap_erlst = (PyArrayObject *)PyArray_SimpleNew(1,limlst_shape,NPY_DOUBLE);
  ap_ierlst = (PyArrayObject *)PyArray_SimpleNew(1,limlst_shape,NPY_INT);
  if (ap_iord == NULL || ap_nnlog == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL || ap_rslst == NULL || ap_erlst == NULL || ap_ierlst == NULL) goto fail;
  iord = (int *)ap_iord->data;
  nnlog = (int *)ap_nnlog->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;
  rslst = (double *)ap_rslst->data;
  erlst = (double *)ap_erlst->data;
  ierlst = (int *)ap_ierlst->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAWFE(quad_thunk, &a, &omega, &integr, &epsabs, &limlst, &limit, &maxp1, &result, &abserr, &neval, &ier, rslst, erlst, ierlst, &lst, alist, blist, rlist, elist, iord, nnlog, chebmo);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  Py_DECREF(ap_nnlog);
  Py_DECREF(ap_alist);
  Py_DECREF(ap_blist);
  Py_DECREF(ap_rlist);
  Py_DECREF(ap_elist);
  Py_DECREF(ap_iord);
  Py_DECREF(ap_chebmo);

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N}i", result, abserr, "neval", neval, "lst", lst, "rslst", PyArray_Return(ap_rslst), "erlst", PyArray_Return(ap_erlst), "ierlst", PyArray_Return(ap_ierlst), ier);
  }
  else {
    Py_DECREF(ap_rslst);
    Py_DECREF(ap_erlst);
    Py_DECREF(ap_ierlst);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  Py_XDECREF(ap_nnlog);
  Py_XDECREF(ap_chebmo);
  Py_XDECREF(ap_rslst);
  Py_XDECREF(ap_erlst);
  Py_XDECREF(ap_ierlst);
  return NULL;
}


static char doc_qawce[] = "[result,abserr,infodict,ier] = _qawce(fun, a, b, c, | args, full_output, epsabs, epsrel, limit)";

static PyObject *quadpack_qawce(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn;

  int      limit;
  npy_intp limit_shape[1];
  int      full_output = 0;
  double   a, b, c, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Oddd|Oiddi", &fcn, &a, &b, &c, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAWCE(quad_thunk, &a, &b, &c, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_iord);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}


static char doc_qawse[] = "[result,abserr,infodict,ier] = _qawse(fun, a, b, (alfa, beta), integr, | args, full_output, epsabs, epsrel, limit)";

static PyObject *quadpack_qawse(PyObject *dummy, PyObject *args) {

  PyArrayObject *ap_alist = NULL, *ap_iord = NULL;
  PyArrayObject *ap_blist = NULL, *ap_elist = NULL;
  PyArrayObject *ap_rlist = NULL;

  PyObject *extra_args = NULL;
  PyObject *fcn;

  int      full_output = 0, integr;
  int      limit=50;
  npy_intp limit_shape[1];
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  double   alfa, beta;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  int ret;
  ccallback_t callback;

  if (!PyArg_ParseTuple(args, "Odd(dd)i|Oiddi", &fcn, &a, &b, &alfa, &beta, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  ret = init_callback(&callback, fcn, extra_args);
  if (ret == -1) {
      return NULL;
  }

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,limit_shape,NPY_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(callback.error_buf) != 0) {
      goto fail;
  }

  DQAWSE(quad_thunk, &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);

  if (free_callback(&callback) != 0) {
      goto fail_free;
  }

  if (full_output) {
    return Py_BuildValue("dd{s:i,s:i,s:N,s:N,s:N,s:N,s:N}i", result, abserr, "neval", neval, "last", last, "iord", PyArray_Return(ap_iord), "alist", PyArray_Return(ap_alist), "blist", PyArray_Return(ap_blist), "rlist", PyArray_Return(ap_rlist), "elist", PyArray_Return(ap_elist),ier);
  }
  else {
    Py_DECREF(ap_alist);
    Py_DECREF(ap_blist);
    Py_DECREF(ap_rlist);
    Py_DECREF(ap_elist);
    Py_DECREF(ap_iord);
    return Py_BuildValue("ddi",result,abserr,ier);
  }

 fail:
  free_callback(&callback);
 fail_free:
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}
