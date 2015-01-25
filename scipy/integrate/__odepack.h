/* This file should be included in the multipack module */
/* $Revision$ */
/*  module_methods:
 {"odeint", (PyCFunction) odepack_odeint, METH_VARARGS|METH_KEYWORDS, doc_odeint},
 */
/* link libraries: (should be listed in separate lines)
   odepack
   linpack_lite
   blas
   mach
 */
/* python files: (to be appended to Multipack.py)
   odepack.py
 */


#if defined(UPPERCASE_FORTRAN)
	#if defined(NO_APPEND_FORTRAN)
		/* nothing to do here */
	#else
		#define LSODA  LSODA_
	#endif
#else
	#if defined(NO_APPEND_FORTRAN)
		#define LSODA  lsoda
	#else
		#define LSODA  lsoda_
	#endif
#endif

void LSODA();

/*
void ode_function(int *n, double *t, double *y, double *ydot)
{
  ydot[0] = -0.04*y[0] + 1e4*y[1]*y[2];
  ydot[2] = 3e7*y[1]*y[1];
  ydot[1] = -ydot[0] - ydot[2];
  return;
}
*/

void ode_function(int *n, double *t, double *y, double *ydot)
{
 /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and set *n to -1 if any
 	-- otherwise place result of calculation in ydot
  */

  PyArrayObject *result_array = NULL;
  PyObject *arg1, *arglist;

  /* Append t to argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    *n = -1;
    return;
  }
  PyTuple_SET_ITEM(arg1, 0, PyFloat_FromDouble(*t)); 
                /* arg1 now owns newly created reference */
  if ((arglist = PySequence_Concat(arg1, multipack_extra_arguments)) == NULL) {
    *n = -1;
    Py_DECREF(arg1);
    return;
  }
  Py_DECREF(arg1);    /* arglist has reference */

  result_array = (PyArrayObject *)call_python_function(multipack_python_function,
                                                       *n, y, arglist, odepack_error);
  if (result_array == NULL) {
    *n = -1;
    Py_DECREF(arglist);
    return;
  }

  if (PyArray_NDIM(result_array) > 1) {
      *n = -1;
      PyErr_Format(PyExc_RuntimeError,
                   "The array return by func must be one-dimensional, but got ndim=%d.",
                   PyArray_NDIM(result_array));
      Py_DECREF(arglist);
      Py_DECREF(result_array);
      return;
  }

  if (PyArray_Size((PyObject *)result_array) != *n) {
      PyErr_Format(PyExc_RuntimeError,
                   "The size of the array returned by func (%ld) does not match the size of y0 (%d).",
                   PyArray_Size((PyObject *)result_array), *n);
      *n = -1;
      Py_DECREF(arglist);
      Py_DECREF(result_array);
      return;
  }

  memcpy(ydot, result_array->data, (*n)*sizeof(double));
  Py_DECREF(result_array);
  Py_DECREF(arglist);
  return;
}


/*
 *  Copy a contiguous matrix at `c` to a Fortran-ordered matrix at `f`.
 *  `ldf` is the leading dimension of the Fortran array at `f`.
 *  `nrows` and `ncols` are the number of rows and columns of the matrix, resp.
 *  If `transposed` is 0, c[i, j] is *(c + ncols*i + j).
 *  If `transposed` is nonzero, c[i, j] is *(c + i + nrows*j)  (i.e. `c` is
 *  stored in F-contiguous order).
 */

static void
copy_array_to_fortran(double *f, int ldf, int nrows, int ncols,
                      double *c, int transposed)
{
    int i, j;
    int row_stride, col_stride;

    /* The strides count multiples of sizeof(double), not bytes. */
    if (transposed) {
        row_stride = 1;
        col_stride = nrows;
    }
    else {
        row_stride = ncols;
        col_stride = 1;
    }
    for (i = 0; i < nrows; ++i) {
        for (j = 0; j < ncols; ++j) {
            double value;
            /* value = c[i,j] */
            value = *(c + row_stride*i + col_stride*j);
            /* f[i,j] = value */
            *(f + ldf*j + i) = value;
        }
    }
}


int ode_jacobian_function(int *n, double *t, double *y, int *ml, int *mu, double *pd, int *nrowpd)
{
  /* This is the function called from the Fortran code it should
        -- use call_python_function to get a multiarrayobject result
	-- check for errors and return -1 if any (though this is ignored 
	                         by calling program).
	-- otherwise place result of calculation in pd
  */

  PyArrayObject *result_array;
  PyObject *arglist, *arg1;

  int ndim, nrows, ncols, dim_error;
  npy_intp *dims;

  /* Append t to argument list */
  if ((arg1 = PyTuple_New(1)) == NULL) {
    *n = -1;
    return -1;
  }
  PyTuple_SET_ITEM(arg1, 0, PyFloat_FromDouble(*t)); 
                /* arg1 now owns newly created reference */
  if ((arglist = PySequence_Concat( arg1, multipack_extra_arguments)) == NULL) {
    *n = -1;
    Py_DECREF(arg1);
    return -1;
  }
  Py_DECREF(arg1);    /* arglist has reference */

  result_array = (PyArrayObject *)call_python_function(multipack_python_jacobian, *n, y,
                                                       arglist, odepack_error);
  if (result_array == NULL) {
    *n = -1;
    Py_DECREF(arglist);
    return -1;
  }

  ncols = *n;
  if (multipack_jac_type == 4) {
      nrows = *ml + *mu + 1;
  }
  else {
      nrows = *n;
  }

  if (!multipack_jac_transpose) {
      int tmp;
      tmp = nrows;
      nrows = ncols;
      ncols = tmp;
  }

  ndim = PyArray_NDIM(result_array);
  if (ndim > 2) {
      PyErr_Format(PyExc_RuntimeError,
                   "The Jacobian array must be two dimensional, but got ndim=%d.",
                   ndim);
      Py_DECREF(arglist);
      Py_DECREF(result_array);
      return -1;
  }

  dims = PyArray_DIMS(result_array);
  dim_error = 0;
  if (ndim == 0) {
      if ((nrows != 1) || (ncols != 1)) {
          dim_error = 1;
      }
  }
  if (ndim == 1) {
      if ((nrows != 1) || (dims[0] != ncols)) {
          dim_error = 1;
      }
  }
  if (ndim == 2) {
      if ((dims[0] != nrows) || (dims[1] != ncols)) {
          dim_error = 1;
      }
  }
  if (dim_error) {
      char *b = "";
      if (multipack_jac_type == 4) {
          b = "banded ";
      }
      PyErr_Format(PyExc_RuntimeError,
                   "Expected a %sJacobian array with shape (%d, %d)",
                   b, nrows, ncols);
      Py_DECREF(arglist);
      Py_DECREF(result_array);
      return -1;
  }

  /*
   *  multipack_jac_type is either 1 (full Jacobian) or 4 (banded Jacobian).
   *  multipack_jac_transpose is !col_deriv, so if multipack_jac_transpose is 0,
   *  the array created by the user is already in Fortran order, and a transpose is
   *  not needed when it is copied to pd.
   */

  if ((multipack_jac_type == 1) && !multipack_jac_transpose) {
      /* Full Jacobian, no transpose needed, so we can use memcpy. */
      memcpy(pd, result_array->data, (*n)*(*nrowpd)*sizeof(double));
  }
  else {
      /*
       *  multipack_jac_type == 4 (banded Jacobian), or
       *  multipack_jac_type == 1 and multipack_jac_transpose == 1.
       *
       *  We can't use memcpy when multipack_jac_type is 4 because the leading
       *  dimension of pd doesn't necessarily equal the number of rows of the
       *  matrix.
       */
      int m;  /* Number of rows in the (full or packed banded) Jacobian. */
      if (multipack_jac_type == 4) {
          m = *ml + *mu + 1;
      }
      else {
          m = *n;
      }
      copy_array_to_fortran(pd, *nrowpd, m, *n, (double *) result_array->data,
                            !multipack_jac_transpose);
  }

  Py_DECREF(arglist);
  Py_DECREF(result_array);
  return 0;
}


int setup_extra_inputs(PyArrayObject **ap_rtol, PyObject *o_rtol, PyArrayObject **ap_atol, PyObject *o_atol, PyArrayObject **ap_tcrit, PyObject *o_tcrit, int *numcrit, int neq)
{
  int itol = 0;
  double tol=1.49012e-8;
  npy_intp one = 1;

  /* Setup tolerances */
  if (o_rtol == NULL) {
    *ap_rtol = (PyArrayObject *)PyArray_SimpleNew(1, &one, NPY_DOUBLE);
    if (*ap_rtol == NULL) PYERR2(odepack_error,"Error constructing relative tolerance.");
    *(double *)(*ap_rtol)->data = tol;                /* Default */
  }
  else {
    *ap_rtol = (PyArrayObject *)PyArray_ContiguousFromObject(o_rtol,NPY_DOUBLE,0,1);
    if (*ap_rtol == NULL) PYERR2(odepack_error,"Error converting relative tolerance.");
    if ((*ap_rtol)->nd == 0); /* rtol is scalar */
    else if ((*ap_rtol)->dimensions[0] == neq)
      itol |= 2;      /* Set rtol array flag */
    else
      PYERR(odepack_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
  }

  if (o_atol == NULL) {
    *ap_atol = (PyArrayObject *)PyArray_SimpleNew(1,&one,NPY_DOUBLE);
    if (*ap_atol == NULL) PYERR2(odepack_error,"Error constructing absolute tolerance");
    *(double *)(*ap_atol)->data = tol;
  }
  else {
    *ap_atol = (PyArrayObject *)PyArray_ContiguousFromObject(o_atol,NPY_DOUBLE,0,1);
    if (*ap_atol == NULL) PYERR2(odepack_error,"Error converting absolute tolerance.");
    if ((*ap_atol)->nd == 0); /* atol is scalar */
    else if ((*ap_atol)->dimensions[0] == neq) 
      itol |= 1;        /* Set atol array flag */
    else
      PYERR(odepack_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
  }
  itol++;             /* increment to get correct value */


  /* Setup t-critical */
  if (o_tcrit != NULL) {
    *ap_tcrit = (PyArrayObject *)PyArray_ContiguousFromObject(o_tcrit,NPY_DOUBLE,0,1);
    if (*ap_tcrit == NULL) PYERR2(odepack_error,"Error constructing critical times.");
    *numcrit = PyArray_Size((PyObject *)(*ap_tcrit));
  }
  return itol;

 fail:       /* Needed for use of PYERR */
  return -1;
}


int compute_lrw_liw(int *lrw, int *liw, int neq, int jt, int ml, int mu, int mxordn, int mxords)
{
  int lrn, lrs, nyh, lmat;

  if (jt == 1 || jt == 2)
    lmat = neq*neq + 2;
  else if (jt == 4 || jt == 5)
    lmat = (2*ml + mu + 1)*neq + 2;
  else PYERR(odepack_error,"Incorrect value for jt");

  if (mxordn < 0) PYERR(odepack_error,"Incorrect value for mxordn");
  if (mxords < 0) PYERR(odepack_error,"Incorrect value for mxords");
  nyh = neq;

  lrn = 20 + nyh*(mxordn+1) + 3*neq;
  lrs = 20 + nyh*(mxords+1) + 3*neq + lmat;

  *lrw = PyArray_MAX(lrn,lrs);
  *liw = 20 + neq;
  return 0;

 fail:
  return -1;
}

static char doc_odeint[] = "[y,{infodict,}istate] = odeint(fun, y0, t, args=(), Dfun=None, col_deriv=0, ml=, mu=, full_output=0, rtol=, atol=, tcrit=, h0=0.0, hmax=0.0, hmin=0.0, ixpr=0.0, mxstep=0.0, mxhnil=0, mxordn=0, mxords=0)\n  yprime = fun(y,t,...)";

static PyObject *odepack_odeint(PyObject *dummy, PyObject *args, PyObject *kwdict) {
  PyObject *fcn, *y0, *p_tout, *o_rtol=NULL, *o_atol=NULL;
  PyArrayObject *ap_y = NULL, *ap_yout= NULL;
  PyArrayObject *ap_rtol=NULL, *ap_atol=NULL;
  PyArrayObject *ap_tout = NULL;
  PyObject *extra_args = NULL;
  PyObject *Dfun = Py_None;
  int      neq, itol=1, itask=1, istate=1, iopt=0, lrw, *iwork, liw, jt=4;
  double   *y, t, *tout, *rtol, *atol, *rwork;
  double   h0=0.0, hmax=0.0, hmin=0.0;
  int      ixpr=0, mxstep=0, mxhnil=0, mxordn=12, mxords=5, ml=-1, mu=-1;
  PyObject *o_tcrit=NULL;
  PyArrayObject *ap_tcrit=NULL;
  PyArrayObject *ap_hu=NULL, *ap_tcur=NULL, *ap_tolsf=NULL, *ap_tsw=NULL;
  PyArrayObject *ap_nst=NULL, *ap_nfe=NULL, *ap_nje=NULL, *ap_nqu=NULL;
  PyArrayObject *ap_mused=NULL;
  int      imxer=0, lenrw=0, leniw=0, col_deriv = 0;
  npy_intp out_sz=0,dims[2];
  int      k, ntimes, crit_ind=0;
  int      allocated = 0, full_output = 0, numcrit=0;
  double   *yout, *yout_ptr, *tout_ptr, *tcrit;
  double   *wa;
  static char *kwlist[] = {"fun","y0","t","args","Dfun","col_deriv","ml","mu","full_output","rtol","atol","tcrit","h0","hmax","hmin","ixpr","mxstep","mxhnil","mxordn","mxords",NULL};

  STORE_VARS();

  if (!PyArg_ParseTupleAndKeywords(args, kwdict, "OOO|OOiiiiOOOdddiiiii", kwlist,
                                   &fcn, &y0, &p_tout, &extra_args, &Dfun,
                                   &col_deriv, &ml, &mu, &full_output, &o_rtol, &o_atol,
                                   &o_tcrit, &h0, &hmax, &hmin, &ixpr, &mxstep, &mxhnil,
                                   &mxordn, &mxords)) {
    return NULL;
  }

  if (o_tcrit == Py_None) {
    o_tcrit = NULL;
  }
  if (o_rtol == Py_None) {
    o_rtol = NULL;
  }
  if (o_atol == Py_None) {
    o_atol = NULL;
  }

  /* Set up jt, ml, and mu */
  if (Dfun == Py_None) jt++;    /* set jt for internally generated */
  if (ml < 0 && mu < 0) jt -= 3;     /* neither ml nor mu given, 
                                        mark jt for full jacobian */
  if (ml < 0) ml = 0;    /* if one but not both are given */
  if (mu < 0) mu = 0;

  INIT_JAC_FUNC(fcn, Dfun, extra_args, col_deriv, odepack_error, jt);

  /* Initial input vector */
  ap_y = (PyArrayObject *)PyArray_ContiguousFromObject(y0, NPY_DOUBLE, 0, 0);
  if (ap_y == NULL) goto fail;
  if (PyArray_NDIM(ap_y) > 1) {
      PyErr_SetString(PyExc_ValueError, "Initial condition y0 must be one-dimensional.");
      goto fail;
  }
  y = (double *) ap_y->data;
  neq = PyArray_Size((PyObject *)ap_y);
  dims[1] = neq;

  /* Set of output times for integration */
  ap_tout = (PyArrayObject *)PyArray_ContiguousFromObject(p_tout, NPY_DOUBLE, 0, 0);
  if (ap_tout == NULL) goto fail;
  if (PyArray_NDIM(ap_tout) > 1) {
      PyErr_SetString(PyExc_ValueError, "Output times t must be one-dimensional.");
      goto fail;
  }
  tout = (double *)ap_tout->data;
  ntimes = PyArray_Size((PyObject *)ap_tout);
  dims[0] = ntimes;
  t = tout[0];

  /* Setup array to hold the output evaluations*/
  ap_yout= (PyArrayObject *)PyArray_SimpleNew(2,dims,NPY_DOUBLE);
  if (ap_yout== NULL) goto fail;
  yout = (double *) ap_yout->data;
  /* Copy initial vector into first row of output */
  memcpy(yout, y, neq*sizeof(double));  /* copy initial value to output */
  yout_ptr = yout + neq;    /* set output pointer to next position */

  itol = setup_extra_inputs(&ap_rtol, o_rtol, &ap_atol, o_atol, &ap_tcrit, o_tcrit, &numcrit, neq);
  if (itol < 0 ) goto fail;  /* Something didn't work */
  rtol = (double *) ap_rtol->data;
  atol = (double *) ap_atol->data;
  if (o_tcrit != NULL) tcrit = (double *)(ap_tcrit->data);

  /* Find size of working arrays*/
  if (compute_lrw_liw(&lrw, &liw, neq, jt, ml, mu, mxordn, mxords) < 0) goto fail;

  if ((wa = (double *)malloc(lrw*sizeof(double) + liw*sizeof(int)))==NULL) {
    PyErr_NoMemory();
    goto fail;
  }
  allocated = 1;
  rwork = wa;
  iwork = (int *)(wa + lrw);

  iwork[0] = ml; iwork[1] = mu;      /* ignored if not needed */

  if (h0 != 0.0 || hmax != 0.0 || hmin != 0.0 || ixpr != 0 || mxstep != 0 || mxhnil != 0 || mxordn != 0 || mxords != 0) {
  rwork[4] = h0; rwork[5] = hmax; rwork[6] = hmin;
  iwork[4] = ixpr; iwork[5] = mxstep; iwork[6] = mxhnil;
  iwork[7] = mxordn; iwork[8] = mxords;
  iopt = 1;
  }
  istate = 1;
  k = 1;

  /* If full output make some useful output arrays */
  if (full_output) {
    out_sz = ntimes-1;
    ap_hu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_DOUBLE);
    ap_tcur = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_DOUBLE);
    ap_tolsf = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_DOUBLE);
    ap_tsw = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_DOUBLE);
    ap_nst = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_INT);
    ap_nfe = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_INT);
    ap_nje = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_INT);
    ap_nqu = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_INT);
    ap_mused = (PyArrayObject *)PyArray_SimpleNew(1,&out_sz,NPY_INT);
    if (ap_hu == NULL || ap_tcur == NULL || ap_tolsf == NULL || ap_tsw == NULL || ap_nst == NULL || ap_nfe == NULL || ap_nje == NULL || ap_nqu == NULL || ap_mused == NULL) goto fail;
  }

  if (o_tcrit != NULL) {itask = 4; rwork[0] = *tcrit;}  /* There are critical points */
  while (k < ntimes && istate > 0) {    /* loop over desired times */

    tout_ptr = tout + k;
    /* Use tcrit if relevant */
    if (itask == 4 && *tout_ptr > *(tcrit + crit_ind)) {crit_ind++; rwork[0] = *(tcrit+crit_ind);}
    if (crit_ind >= numcrit) itask = 1;  /* No more critical values */

    LSODA(ode_function, &neq, y, &t, tout_ptr, &itol, rtol, atol, &itask, &istate, &iopt, rwork, &lrw, iwork, &liw, ode_jacobian_function, &jt);
    if (full_output) {
      *((double *)ap_hu->data + (k-1)) = rwork[10];
      *((double *)ap_tcur->data + (k-1)) = rwork[12];
      *((double *)ap_tolsf->data + (k-1)) = rwork[13];
      *((double *)ap_tsw->data + (k-1)) = rwork[14];
      *((int *)ap_nst->data + (k-1)) = iwork[10];
      *((int *)ap_nfe->data + (k-1)) = iwork[11];
      *((int *)ap_nje->data + (k-1)) = iwork[12];
      *((int *)ap_nqu->data + (k-1)) = iwork[13];
      if (istate == -5 || istate == -4) {
        imxer = iwork[15];
      } else {
        imxer = -1;
      }
      lenrw = iwork[16];
      leniw = iwork[17];
      *((int *)ap_mused->data + (k-1)) = iwork[18];
    }
    if (PyErr_Occurred()) goto fail;
    memcpy(yout_ptr, y, neq*sizeof(double));  /* copy integration result to output*/
    yout_ptr += neq;  k++;
  }

  RESTORE_JAC_FUNC();

  Py_DECREF(extra_args);
  Py_DECREF(ap_atol);
  Py_DECREF(ap_rtol);
  Py_XDECREF(ap_tcrit);
  Py_DECREF(ap_y);
  Py_DECREF(ap_tout);
  free(wa);

  /* Do Full output */
    if (full_output) {
      return Py_BuildValue("N{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:i,s:i,s:i,s:N}i",PyArray_Return(ap_yout),
                      "hu",PyArray_Return(ap_hu),
                      "tcur",PyArray_Return(ap_tcur),
                      "tolsf",PyArray_Return(ap_tolsf),
                      "tsw",PyArray_Return(ap_tsw),
                      "nst",PyArray_Return(ap_nst),
                      "nfe",PyArray_Return(ap_nfe),
                      "nje",PyArray_Return(ap_nje),
                      "nqu",PyArray_Return(ap_nqu),
                      "imxer",imxer,
                      "lenrw",lenrw,
                      "leniw",leniw,
                      "mused",PyArray_Return(ap_mused),
                      istate);
    }
    else {
      return Py_BuildValue("Ni",PyArray_Return(ap_yout),istate);
    }
    
 fail:
  RESTORE_JAC_FUNC();
  Py_XDECREF(extra_args);
  Py_XDECREF(ap_y);
  Py_XDECREF(ap_rtol);
  Py_XDECREF(ap_atol);
  Py_XDECREF(ap_tcrit);
  Py_XDECREF(ap_tout);
  Py_XDECREF(ap_yout);
  if (allocated) free(wa);
  if (full_output) {
    Py_XDECREF(ap_hu);
    Py_XDECREF(ap_tcur);
    Py_XDECREF(ap_tolsf);
    Py_XDECREF(ap_tsw);
    Py_XDECREF(ap_nst);
    Py_XDECREF(ap_nfe);
    Py_XDECREF(ap_nje);
    Py_XDECREF(ap_nqu);
    Py_XDECREF(ap_mused);
  }
  return NULL;
  
}
