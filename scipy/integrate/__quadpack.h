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


static int already_printed_python_error = 0;

#define QUAD_INIT_FUNC(fun,arg) {\
  INIT_FUNC(fun,arg,quadpack_error); \
  already_printed_python_error = 0;\
}

double quad_function(double *x) {

  double d_result;
  PyObject *arg1 = NULL, *arglist=NULL, *result=NULL;
  PyNumberMethods *nb;

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

  Py_DECREF(arg1);    /* arglist has the reference to Float object. */
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

  npy_intp limit=50;
  int      full_output = 0;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;

  STORE_VARS();

  if (!PyArg_ParseTuple(args, "Odd|Oiddi", &fcn, &a, &b, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1) 
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAGSE(quad_function, &a, &b, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);

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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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

  npy_intp limit=50;
  int      full_output = 0;
  double   bound, epsabs=1.49e-8, epsrel=1.49e-8;
  int      inf, neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  
  STORE_VARS();
  
  if (!PyArg_ParseTuple(args, "Odi|Oiddi", &fcn, &bound, &inf, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args);

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAGIE(quad_function, &bound, &inf, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }

  Py_DECREF(extra_args);

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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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

  npy_intp limit=50, npts2;
  int      full_output = 0;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  int      *level, *ndin;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;
  double   *pts, *points;

  STORE_VARS();

  if (!PyArg_ParseTuple(args, "OddO|Oiddi", &fcn, &a, &b, &o_points, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  ap_points = (PyArrayObject *)PyArray_ContiguousFromObject(o_points, PyArray_DOUBLE, 1, 1);
  if (ap_points == NULL) goto fail;
  npts2 = ap_points->dimensions[0];
  points = (double *)ap_points->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_pts = (PyArrayObject *)PyArray_SimpleNew(1,&npts2,PyArray_DOUBLE);
  ap_level = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_ndin = (PyArrayObject *)PyArray_SimpleNew(1,&npts2,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL || ap_pts == NULL || ap_level == NULL || ap_ndin == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;
  pts = (double *)ap_pts->data;
  level = (int *)ap_level->data;
  ndin = (int *)ap_level->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {    
    DQAGPE(quad_function, &a, &b, &npts2, points, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, pts, iord, level, ndin, &last);
  }

  RESTORE_FUNC()

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);
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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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

  npy_intp limit=50, sz[2];
  int      full_output = 0, maxp1=50, icall=1;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, integr=1, last=0, momcom=0, *iord;
  int      *nnlog;
  double   result=0.0, abserr=0.0, omega=0.0;
  double   *chebmo;
  double   *alist, *blist, *rlist, *elist;
  
  STORE_VARS();

  if (!PyArg_ParseTuple(args, "Odddi|OiddiiiiO", &fcn, &a, &b, &omega, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit, &maxp1, &icall, &momcom, &o_chebmo)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  if (o_chebmo != NULL) {
    ap_chebmo = (PyArrayObject *)PyArray_ContiguousFromObject(o_chebmo, PyArray_DOUBLE, 2, 2);
    if (ap_chebmo == NULL) goto fail;
    if (ap_chebmo->dimensions[1] != maxp1 || ap_chebmo->dimensions[0] != 25)
      PYERR(quadpack_error,"Chebyshev moment array has the wrong size.");
  }
  else {
    sz[0] = 25;
    sz[1] = maxp1;
    ap_chebmo = (PyArrayObject *)PyArray_SimpleNew(2,sz,PyArray_DOUBLE);
    if (ap_chebmo == NULL) goto fail;
  }
  chebmo = (double *) ap_chebmo->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_nnlog = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_nnlog == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  nnlog = (int *)ap_nnlog->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAWOE(quad_function, &a, &b, &omega, &integr, &epsabs, &epsrel, &limit, &icall, &maxp1, &result, &abserr, &neval, &ier, &last, alist, blist, rlist, elist, iord, nnlog, &momcom, chebmo);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);

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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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

  npy_intp limlst = 50, limit=50, sz[2];
  int      full_output = 0, maxp1=50;
  double   a, epsabs=1.49e-8;
  int      neval=0, ier=6, integr=1, *iord;
  int      lst, *nnlog, *ierlst;
  double   *chebmo, *rslst, *erlst;
  double   result=0.0, abserr=0.0, omega=0.0;
  double   *alist, *blist, *rlist, *elist;

  STORE_VARS();

  if (!PyArg_ParseTuple(args, "Oddi|Oidiii", &fcn, &a, &omega, &integr, &extra_args, &full_output, &epsabs, &limlst, &limit, &maxp1)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  sz[0] = 25;
  sz[1] = maxp1;
  ap_chebmo = (PyArrayObject *)PyArray_SimpleNew(2,sz,PyArray_DOUBLE);
  if (ap_chebmo == NULL) goto fail;
  chebmo = (double *) ap_chebmo->data;

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_nnlog = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rslst = (PyArrayObject *)PyArray_SimpleNew(1,&limlst,PyArray_DOUBLE);
  ap_erlst = (PyArrayObject *)PyArray_SimpleNew(1,&limlst,PyArray_DOUBLE);
  ap_ierlst = (PyArrayObject *)PyArray_SimpleNew(1,&limlst,PyArray_INT);
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

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAWFE(quad_function, &a, &omega, &integr, &epsabs, &limlst, &limit, &maxp1, &result, &abserr, &neval, &ier, rslst, erlst, ierlst, &lst, alist, blist, rlist, elist, iord, nnlog, chebmo);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);
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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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

  npy_intp limit=50;
  int      full_output = 0;
  double   a, b, c, epsabs=1.49e-8, epsrel=1.49e-8;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;

  STORE_VARS();
  
  if (!PyArg_ParseTuple(args, "Oddd|Oiddi", &fcn, &a, &b, &c, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAWCE(quad_function, &a, &b, &c, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);

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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
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
  npy_intp limit=50;
  double   a, b, epsabs=1.49e-8, epsrel=1.49e-8;
  double   alfa, beta;
  int      neval=0, ier=6, last=0, *iord;
  double   result=0.0, abserr=0.0;
  double   *alist, *blist, *rlist, *elist;

  STORE_VARS();
  
  if (!PyArg_ParseTuple(args, "Odd(dd)i|Oiddi", &fcn, &a, &b, &alfa, &beta, &integr, &extra_args, &full_output, &epsabs, &epsrel, &limit)) return NULL;

  /* Need to check that limit is bigger than 1 */
  if (limit < 1)
    return Py_BuildValue("ddi",result,abserr,ier);

  QUAD_INIT_FUNC(fcn,extra_args)

  /* Setup iwork and work arrays */
  ap_iord = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_INT);
  ap_alist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_blist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_rlist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  ap_elist = (PyArrayObject *)PyArray_SimpleNew(1,&limit,PyArray_DOUBLE);
  if (ap_iord == NULL || ap_alist == NULL || ap_blist == NULL || ap_rlist == NULL || ap_elist == NULL) goto fail;
  iord = (int *)ap_iord->data;
  alist = (double *)ap_alist->data;
  blist = (double *)ap_blist->data;
  rlist = (double *)ap_rlist->data;
  elist = (double *)ap_elist->data;

  if (setjmp(quadpack_jmpbuf)) {
    goto fail;
  }
  else {
    DQAWSE(quad_function, &a, &b, &alfa, &beta, &integr, &epsabs, &epsrel, &limit, &result, &abserr, &neval, &ier, alist, blist, rlist, elist, iord, &last);
  }

  RESTORE_FUNC();

  if (PyErr_Occurred()) {
    ier = 80;             /* Python error */
    PyErr_Clear();
  }
  Py_DECREF(extra_args);

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
  RESTORE_FUNC();
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_alist);
  Py_XDECREF(ap_blist);
  Py_XDECREF(ap_rlist);
  Py_XDECREF(ap_elist);
  Py_XDECREF(ap_iord);
  return NULL;
}




