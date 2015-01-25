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
  Error=-3, 
  Not_Callable=-2,
  Invalid_Ctype=-1,
  /* Acceptable returns below */
  Callable=1,
  Valid_Ctype=2,
  Valid_Multivariate_Ctype=3
} FuncType;


/* Checks a callable object:
   Returns Valid_Multivariate_Ctype if the Python Object is a CType Function of the type (int, double*) -> (double)
   Returns Valid_Ctype if the Python Object is a CType Function of the type (double) -> (double)
   Return  Callable if it is a Python callable (including a Ctypes function with no stored types)
   Returns Invalid_Ctype if it is a CType Function but of the incorrect type
   Returns Not_Callable if it is not a Python callable
   Returns Error if other error occurs.
*/

static FuncType
get_func_type(PyObject *func) {
  PyObject *ctypes_module, *CFuncPtr, *check, *c_double, *c_int;
  int result;
  
  if (!PyCallable_Check(func)) {
    PyErr_SetString(quadpack_error, "quad: first argument is not callable");
    return Not_Callable;
  }

  ctypes_module = PyImport_ImportModule("ctypes");
  if (ctypes_module == NULL) { /* We don't have ctypes... just return Callable */
      PyErr_Clear();
      return Callable;
  }
  CFuncPtr = PyObject_GetAttrString(ctypes_module, "_CFuncPtr");
  if (CFuncPtr == NULL) {
    Py_DECREF(ctypes_module);
    return Error;
  }
  result = PyObject_TypeCheck(func, (PyTypeObject *)CFuncPtr);
  Py_DECREF(CFuncPtr);
  if (!result) {
    Py_DECREF(ctypes_module);
    return Callable;
  }
  
  /* We are a ctypes function */
  /* Check for restype and argtypes */
  if (!PyObject_HasAttrString(func, "restype") ||
      !PyObject_HasAttrString(func, "argtypes")) {
    Py_DECREF(ctypes_module);
    return Callable;
  }
  
  c_double = PyObject_GetAttrString(ctypes_module, "c_double");
  c_int = PyObject_GetAttrString(ctypes_module, "c_int");
  Py_DECREF(ctypes_module);
  
  check = PyObject_GetAttrString(func, "restype");
  if (check != c_double) {  /* Checking for pointer identity */
    goto fail;
  }
  Py_DECREF(check);
  check = PyObject_GetAttrString(func, "argtypes");
  if (!PyTuple_Check(check) || (PyTuple_GET_SIZE(check) != 1) || 
      (PyTuple_GET_ITEM(check, 0) != c_double)) {
    if ((PyTuple_GET_ITEM(check, 0) == c_int) && (PyTuple_GET_ITEM(check, 1) == c_double)){   
      /*If first param is int and second is of double type, we have valid_multivariate */
      Py_DECREF(check);
      Py_DECREF(c_double);
      Py_DECREF(c_int);
      return Valid_Multivariate_Ctype;
    }
    else 
      goto fail;
  }
  Py_DECREF(check);
  Py_DECREF(c_double);
  Py_DECREF(c_int);
  return Valid_Ctype;
  
fail:
  Py_DECREF(check);
  Py_XDECREF(c_double);
  Py_XDECREF(c_int);
  PyErr_SetString(quadpack_error,
                  "quad: first argument is a ctypes function pointer with incorrect signature");
  return Invalid_Ctype;
}

/*
  This is the equivalent to the above which only uses ctypes calls... However it is about 8x slower in a tight loop
  than the above code: 
 
static _sp_double_func
get_ctypes_function_pointer(PyObject *obj) {
  PyObject *ctypes_module=NULL, *c_void_p=NULL, *castfunc=NULL, *value=NULL, *result=NULL;
  void *final;

  ctypes_module = PyImport_ImportModule("ctypes");
  if (ctypes_module == NULL) return NULL;
  
  castfunc = PyObject_GetAttrString(ctypes_module, "cast");
  if (castfunc == NULL) goto fail;
  c_void_p = PyObject_GetAttrString(ctypes_module, "c_void_p");
  if (c_void_p == NULL) goto fail;
  result = PyObject_CallFunctionObjArgs(castfunc, obj, c_void_p);
  if (result == NULL) goto fail;
  value = PyObject_GetAttrString(result, "value");
  if (value == NULL) goto fail;
  final = PyLong_AsVoidPtr(value);
  if (PyErr_Occurred()) goto fail;
  
  Py_DECREF(ctypes_module);
  Py_DECREF(castfunc);
  Py_DECREF(c_void_p);
  Py_DECREF(result);
  Py_DECREF(value);
  
  return (_sp_double_func) final;
  
fail:
  Py_XDECREF(ctypes_module);
  Py_XDECREF(castfunc);
  Py_XDECREF(c_void_p);
  Py_XDECREF(result);
  Py_XDECREF(value);
  return NULL;
}
*/

double quad_function2(double *x) {
  return quadpack_ctypes_function(*x);
}

double quad_function(double *x) {

  double d_result;
  PyObject *arg1 = NULL, *arglist=NULL, *result=NULL;

  /* Build argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) goto fail;

  PyTuple_SET_ITEM(arg1, 0, PyFloat_FromDouble(*x));
                /* arg1 now owns reference to Float object*/
  if ((arglist = PySequence_Concat( arg1, quadpack_extra_arguments)) == NULL) goto fail;

  /* Call function object --- stored as a global variable.  Extra
          arguments are in another global variable.
   */
  if ((result = PyEval_CallObject(quadpack_python_function, arglist))==NULL) goto fail;

  /* Have to do own error checking because PyFloat_AsDouble returns -1 on
     error -- making that return value from the function unusable.
     No; Solution is to test for Python Error Occurrence if -1 is return of PyFloat_AsDouble.
  */

  d_result = PyFloat_AsDouble(result);
  if (PyErr_Occurred())
    PYERR(quadpack_error, "Supplied function does not return a valid float.")

  Py_DECREF(arg1); /* arglist has the reference to Float object. */
  Py_DECREF(arglist);
  Py_DECREF(result);

  return d_result;

 fail:
  Py_XDECREF(arg1);
  Py_XDECREF(arglist);
  Py_XDECREF(result);
  longjmp(quadpack_jmpbuf, 1);
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Odd|Oiddi", &fcn, &a, &b, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable) 
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;

    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAGSE(quad_function, &a, &b, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, 
    blist, rlist, elist, iord, &last);
    }

    quad_restore_func(&storevar, &ier);
  }
  else {
      /* Can't allow another thread to run because of the global variables
         quadpack_raw_function and quad_function2 being used */
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAGSE(call_c_multivariate, &a, &b, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_c_multivariate(&storevar);
    }
    else{ /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAGSE(quad_function2, &a, &b, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Odi|Oiddi", &fcn, &bound, &inf, &extra_args, 
                        &full_output, &epsabs, &epsrel, &limit)) 
    return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable) 
      return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;

    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAGIE(quad_function, &bound, &inf, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
    }
    quad_restore_func(&storevar, &ier);
  }
  else {
      /* Can't allow another thread to run because of the global variables
         quadpack_raw_function and quad_function2 being used */
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAGIE(call_c_multivariate, &bound, &inf, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAGIE(quad_function2, &bound, &inf, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "OddO|Oiddi", &fcn, &a, &b, &o_points, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable)
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;

    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAGPE(quad_function, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, pts, iord, level, ndin, &last);
    }
    
    quad_restore_func(&storevar, &ier);
  }
  else {
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAGPE(call_c_multivariate, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, pts, iord, level, ndin, &last);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */ 
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAGPE(quad_function2, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, pts, iord, level, ndin, &last);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Odddi|OiddiiiiO", &fcn, &a, &b, &omega, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit, &maxp1, &icall, &momcom, &o_chebmo)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable)
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;
    
    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAWOE(quad_function, &a, &b, &omega, &integr, &epsabs, &epsrel, &limit, &icall, &maxp1, &result, &abserr, &neval, &ier, &last, alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);
    }

    quad_restore_func(&storevar, &ier);

  }
  else {
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAWOE(call_c_multivariate, &a, &b, &omega, &integr, &epsabs, &epsrel, &limit, &icall, &maxp1, &result, &abserr, &neval, &ier, &last, alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAWOE(quad_function2, &a, &b, &omega, &integr, &epsabs, &epsrel, &limit, &icall, &maxp1, &result, &abserr, &neval, &ier, &last, alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Oddi|Oidiii", &fcn, &a, &omega, &integr, &extra_args, &full_output, &epsabs, &limlst, &limit, &maxp1)) return NULL;
  limit_shape[0] = limit;
  limlst_shape[0] = limlst;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable)
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;

    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAWFE(quad_function, &a, &omega, &integr, &epsabs, &limlst, &limit, &maxp1, &result, &abserr, &neval, &ier, rslst, erlst, ierlst, &lst, alist, blist, rlist, elist, iord, nnlog, chebmo);
    }

    quad_restore_func(&storevar, &ier);
  }
  else {
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAWFE(call_c_multivariate, &a, &omega, &integr, &epsabs, &limlst, &limit, &maxp1, &result, &abserr, &neval, &ier, rslst, erlst, ierlst, &lst, alist, blist, rlist, elist, iord, nnlog, chebmo);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAWFE(quad_function2, &a, &omega, &integr, &epsabs, &limlst, &limit, &maxp1, &result, &abserr, &neval, &ier, rslst, erlst, ierlst, &lst, alist, blist, rlist, elist, iord, nnlog, chebmo);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Oddd|Oiddi", &fcn, &a, &b, &c, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable)
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;

    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAWCE(quad_function, &a, &b, &c, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
    }

    quad_restore_func(&storevar, &ier);

  } 
  else {
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAWCE(call_c_multivariate, &a, &b, &c, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAWCE(quad_function2, &a, &b, &c, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_ctypes_func(&storevar);
    }
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
  FuncType func_type;
  QStorage storevar;

  if (!PyArg_ParseTuple(args, "Odd(dd)i|Oiddi", &fcn, &a, &b, &alfa, &beta, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;
  limit_shape[0] = limit;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  if ((func_type = get_func_type(fcn)) < Callable)
    return NULL;

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

  if (func_type == Callable) {
    if (quad_init_func(&storevar, fcn, extra_args) == NPY_FAIL)
      goto fail;
   
    if (setjmp(quadpack_jmpbuf)) {
      quad_restore_func(&storevar, NULL);
      goto fail;
    }
    else {
      DQAWSE(quad_function, &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
    }

    quad_restore_func(&storevar, &ier);
  }
  else {
    if (func_type != Valid_Ctype) { /* func_type == VALID_MULTIVARIATE_CTYPE */
      if (init_c_multivariate(&storevar, fcn, extra_args) == NPY_FAIL) goto fail;
      DQAWSE(call_c_multivariate, &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_c_multivariate(&storevar);
    }
    else { /* func_type == VALID_CTYPE */
      if (init_ctypes_func(&storevar, fcn) == NPY_FAIL) goto fail;
      DQAWSE(quad_function2, &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
      restore_ctypes_func(&storevar);
    }
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
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}
