/* Anti-Copyright
 *
 * I hereby release this code into the PUBLIC DOMAIN AS IS.  There is no
 * support, warranty, or guarantee.  I will gladly accept comments, bug 
 * reports, and patches, however.
 *
 * Robert Kern
 * kern@caltech.edu
 *
 */

#include "odrpack.h"


void F_FUNC(dodrc,DODRC)(void (*fcn)(int *n, int *m, int *np, int *nq, int *ldn, int *ldm, 
            int *ldnp, double *beta, double *xplusd, int *ifixb, int *ifixx, 
            int *ldifx, int *ideval, double *f, double *fjacb, double *fjacd, 
            int *istop), 
           int *n, int *m, int *np, int *nq, double *beta, double *y, int *ldy,
           double *x, int *ldx, double *we, int *ldwe, int *ld2we, double *wd,
           int *ldwd, int *ld2wd, int *ifixb, int *ifixx, int *ldifx, int *job,
           int *ndigit, double *taufac, double *sstol, double *partol, 
           int *maxit, int *iprint, int *lunerr, int *lunrpt, double *stpb,
           double *stpd, int *ldstpd, double *sclb, double *scld, int *ldscld,
           double *work, int *lwork, int *iwork, int *liwork, int *info);
void F_FUNC(dwinf,DWINF)(int *n, int *m, int *np, int *nq, int *ldwe, int *ld2we, int *isodr,
        int *delta, int *eps, int *xplus, int *fn, int *sd, int *vcv, int *rvar,
        int *wss, int *wssde, int *wssep, int *rcond, int *eta, int *olmav, 
        int *tau, int *alpha, int *actrs, int *pnorm, int *rnors, int *prers,
        int *partl, int *sstol, int *taufc, int *apsma, int *betao, int *betac,
        int *betas, int *betan, int *s, int *ss, int *ssf, int *qraux, int *u,
        int *fs, int *fjacb, int *we1, int *diff, int *delts, int *deltn, 
        int *t, int *tt, int *omega, int *fjacd, int *wrk1, int *wrk2, 
        int *wrk3, int *wrk4, int *wrk5, int *wrk6, int *wrk7, int *lwkmn);
void F_FUNC(dluno,DLUNO)(int *lun, char *fn, int fnlen);
void F_FUNC(dlunc,DLUNC)(int *lun);



/* callback to pass to DODRC; calls the Python function in the global structure |odr_global| */
void fcn_callback(int *n, int *m, int *np, int *nq, int *ldn, int *ldm,
                  int *ldnp, double *beta, double *xplusd, int *ifixb,
                  int *ifixx, int *ldfix, int *ideval, double *f,
                  double *fjacb, double *fjacd, int *istop)
{
  PyObject *arg01, *arglist;
  PyObject *result;
  PyArrayObject *result_array = NULL;
  PyArrayObject *pyXplusD;

  arg01 = PyTuple_New(2);

  if (*m != 1)
    {
      npy_intp dim2[2];
      dim2[0] = *m;
      dim2[1] = *n;
      pyXplusD = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
      memcpy(pyXplusD->data, (void *)xplusd, (*m) * (*n) * sizeof(double));
    }
  else
    {
      npy_intp dim1[1];
      dim1[0] = *n;
      pyXplusD = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      memcpy(pyXplusD->data, (void *)xplusd, (*n) * sizeof(double));
    }

  PyTuple_SetItem(arg01, 0, odr_global.pyBeta);
  Py_INCREF(odr_global.pyBeta);
  PyTuple_SetItem(arg01, 1, (PyObject *) pyXplusD);
  Py_INCREF((PyObject *) pyXplusD);

  if (odr_global.extra_args != NULL)
    {
      arglist = PySequence_Concat(arg01, odr_global.extra_args);
    }
  else
    {
      arglist = PySequence_Tuple(arg01);        /* make a copy */
    }

  Py_DECREF(arg01);
  *istop = 0;

  memcpy(((PyArrayObject *) (odr_global.pyBeta))->data, (void *)beta,
         (*np) * sizeof(double));

  if ((*ideval % 10) >= 1)
    {
      /* compute f with odr_global.fcn */

      if (odr_global.fcn == NULL)
        {
          /* we don't have a function to call */
          PYERR2(odr_error, "Function has not been initialized");
        }

      if ((result = PyEval_CallObject(odr_global.fcn, arglist)) == NULL)
        {
          PyObject *tmpobj, *str1;

          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }

          PyErr_Print();
          tmpobj = PyObject_GetAttrString(odr_global.fcn, "__name__");
          if (tmpobj == NULL)
            goto fail;

          str1 =
            PyString_FromString
            ("Error occurred while calling the Python function named ");
          if (str1 == NULL)
            {
              Py_DECREF(tmpobj);
              goto fail;
            }
          PyString_ConcatAndDel(&str1, tmpobj);
          PyErr_SetString(odr_error, PyString_AsString(str1));
          Py_DECREF(str1);
          goto fail;
        }

      if ((result_array =
           (PyArrayObject *) PyArray_ContiguousFromObject(result,
                                                          NPY_DOUBLE, 0,
                                                          2)) == NULL)
        {
          PYERR2(odr_error,
                 "Result from function call is not a proper array of floats.");
        }

      memcpy((void *)f, result_array->data, (*n) * (*nq) * sizeof(double));
      Py_DECREF(result_array);
    }

  if (((*ideval) / 10) % 10 >= 1)
    {
      /* compute fjacb with odr_global.fjacb */

      if (odr_global.fjacb == NULL)
        {
          /* we don't have a function to call */
          PYERR2(odr_error, "Function has not been initialized");
        }

      if ((result = PyEval_CallObject(odr_global.fjacb, arglist)) == NULL)
        {
          PyObject *tmpobj, *str1;

          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }

          PyErr_Print();
          tmpobj = PyObject_GetAttrString(odr_global.fjacb, "__name__");
          if (tmpobj == NULL)
            goto fail;

          str1 =
            PyString_FromString
            ("Error occurred while calling the Python function named ");
          if (str1 == NULL)
            {
              Py_DECREF(tmpobj);
              goto fail;
            }
          PyString_ConcatAndDel(&str1, tmpobj);
          PyErr_SetString(odr_error, PyString_AsString(str1));
          Py_DECREF(str1);
          goto fail;
        }

      if ((result_array =
           (PyArrayObject *) PyArray_ContiguousFromObject(result,
                                                          NPY_DOUBLE, 0,
                                                          2)) == NULL)
        {
          PYERR2(odr_error,
                 "Result from function call is not a proper array of floats.");
        }

      if (*nq != 1 && *np != 1)
        {
          /* result_array should be rank-3 */

          if (result_array->nd != 3)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "Beta Jacobian is not rank-3");
            }
        }
      else if (*nq == 1)
        {
          /* result_array should be rank-2 */

          if (result_array->nd != 2)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "Beta Jacobian is not rank-2");
            }
        }

      memcpy((void *)fjacb, result_array->data,
             (*n) * (*nq) * (*np) * sizeof(double));
      Py_DECREF(result_array);

    }

  if (((*ideval) / 100) % 10 >= 1)
    {
      /* compute fjacd with odr_global.fjacd */

      if (odr_global.fjacd == NULL)
        {
          /* we don't have a function to call */
          PYERR2(odr_error, "fjcad has not been initialized");
        }

      if ((result = PyEval_CallObject(odr_global.fjacd, arglist)) == NULL)
        {
          PyObject *tmpobj, *str1;

          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }

          PyErr_Print();
          tmpobj = PyObject_GetAttrString(odr_global.fjacd, "__name__");
          if (tmpobj == NULL)
            goto fail;

          str1 =
            PyString_FromString
            ("Error occurred while calling the Python function named ");
          if (str1 == NULL)
            {
              Py_DECREF(tmpobj);
              goto fail;
            }
          PyString_ConcatAndDel(&str1, tmpobj);
          PyErr_SetString(odr_error, PyString_AsString(str1));
          Py_DECREF(str1);
          goto fail;
        }

      if ((result_array =
           (PyArrayObject *) PyArray_ContiguousFromObject(result,
                                                          NPY_DOUBLE, 0,
                                                          2)) == NULL)
        {
          PYERR2(odr_error,
                 "Result from function call is not a proper array of floats.");
        }

      if (*nq != 1 && *m != 1)
        {
          /* result_array should be rank-3 */

          if (result_array->nd != 3)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-3");
            }
        }
      else if (*nq == 1 && *m != 1)
        {
          /* result_array should be rank-2 */

          if (result_array->nd != 2)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-2");
            }
        }
      else if (*nq == 1 && *m == 1)
        {
          /* result_array should be rank-1 */

          if (result_array->nd != 1)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-1");
            }
        }

      memcpy((void *)fjacd, result_array->data,
             (*n) * (*nq) * (*m) * sizeof(double));
      Py_DECREF(result_array);
    }

  Py_DECREF(result);
  Py_DECREF(arglist);
  Py_DECREF(pyXplusD);

  return;

fail:
  Py_XDECREF(result);
  Py_XDECREF(arglist);
  Py_XDECREF(pyXplusD);
  *istop = -1;
  return;
}


/* generates Python output from the raw output from DODRC */
PyObject *gen_output(int n, int m, int np, int nq, int ldwe, int ld2we,
                     PyArrayObject * beta, PyArrayObject * work,
                     PyArrayObject * iwork, int isodr, int info,
                     int full_output)
{
  PyArrayObject *sd_beta, *cov_beta;

  int delta, eps, xplus, fn, sd, vcv, rvar, wss, wssde, wssep, rcond;
  int eta, olmav, tau, alpha, actrs, pnorm, rnors, prers, partl, sstol;
  int taufc, apsma, betao, betac, betas, betan, s, ss, ssf, qraux, u;
  int fs, fjacb, we1, diff, delts, deltn, t, tt, omega, fjacd;
  int wrk1, wrk2, wrk3, wrk4, wrk5, wrk6, wrk7, lwkmn;

  PyObject *retobj;

  npy_intp dim1[1], dim2[2];

  if (info == 50005) {
      /* fatal error in fcn call; return NULL to propogate the exception */

      return NULL;
  }

  lwkmn = work->dimensions[0];

  F_FUNC(dwinf,DWINF)(&n, &m, &np, &nq, &ldwe, &ld2we, &isodr,
        &delta, &eps, &xplus, &fn, &sd, &vcv, &rvar, &wss, &wssde,
        &wssep, &rcond, &eta, &olmav, &tau, &alpha, &actrs, &pnorm,
        &rnors, &prers, &partl, &sstol, &taufc, &apsma, &betao, &betac,
        &betas, &betan, &s, &ss, &ssf, &qraux, &u, &fs, &fjacb, &we1,
        &diff, &delts, &deltn, &t, &tt, &omega, &fjacd, &wrk1, &wrk2,
        &wrk3, &wrk4, &wrk5, &wrk6, &wrk7, &lwkmn);

  /* convert FORTRAN indices to C indices */
  delta--;
  eps--;
  xplus--;
  fn--;
  sd--;
  vcv--;
  rvar--;
  wss--;
  wssde--;
  wssep--;
  rcond--;
  eta--;
  olmav--;
  tau--;
  alpha--;
  actrs--;
  pnorm--;
  rnors--;
  prers--;
  partl--;
  sstol--;
  taufc--;
  apsma--;
  betao--;
  betac--;
  betas--;
  betan--;
  s--;
  ss--;
  ssf--;
  qraux--;
  u--;
  fs--;
  fjacb--;
  we1--;
  diff--;
  delts--;
  deltn--;
  t--;
  tt--;
  omega--;
  fjacd--;
  wrk1--;
  wrk2--;
  wrk3--;
  wrk4--;
  wrk5--;
  wrk6--;
  wrk7--;

  dim1[0] = beta->dimensions[0];
  sd_beta = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
  dim2[0] = beta->dimensions[0];
  dim2[1] = beta->dimensions[0];
  cov_beta = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);

  memcpy(sd_beta->data, (void *)((double *)(work->data) + sd),
         np * sizeof(double));
  memcpy(cov_beta->data, (void *)((double *)(work->data) + vcv),
         np * np * sizeof(double));

  if (!full_output)
    {
      retobj = Py_BuildValue("OOO", PyArray_Return(beta),
                             PyArray_Return(sd_beta),
                             PyArray_Return(cov_beta));
      Py_DECREF((PyObject *) sd_beta);
      Py_DECREF((PyObject *) cov_beta);

      return retobj;
    }
  else
    {
      PyArrayObject *deltaA, *epsA, *xplusA, *fnA;
      double res_var, sum_square, sum_square_delta, sum_square_eps;
      double inv_condnum, rel_error;
      PyObject *work_ind;

      work_ind =
        Py_BuildValue
        ("{s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i}",
         "delta", delta, "eps", eps, "xplus", xplus, "fn", fn, "sd", sd, "sd",
         vcv, "rvar", rvar, "wss", wss, "wssde", wssde, "wssep", wssep,
         "rcond", rcond, "eta", eta, "olmav", olmav, "tau", tau, "alpha",
         alpha, "actrs", actrs, "pnorm", pnorm, "rnors", rnors, "prers",
         prers, "partl", partl, "sstol", sstol, "taufc", taufc, "apsma",
         apsma, "betao", betao, "betac", betac, "betas", betas, "betan",
         betan, "s", s, "ss", ss, "ssf", ssf, "qraux", qraux, "u", u, "fs",
         fs, "fjacb", fjacb, "we1", we1, "diff", diff, "delts", delts,
         "deltn", deltn, "t", t, "tt", tt, "omega", omega, "fjacd", fjacd,
         "wrk1", wrk1, "wrk2", wrk2, "wrk3", wrk3, "wrk4", wrk4, "wrk5", wrk5,
         "wrk6", wrk6, "wrk7", wrk7);

      if (m == 1)
        {
          dim1[0] = n;
          deltaA =
            (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
          xplusA =
            (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
        }
      else
        {
          dim2[0] = m;
          dim2[1] = n;
          deltaA =
            (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
          xplusA =
            (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
        }

      if (nq == 1)
        {
          dim1[0] = n;
          epsA = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
          fnA = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
        }
      else
        {
          dim2[0] = nq;
          dim2[1] = n;
          epsA = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
          fnA = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
        }

      memcpy(deltaA->data, (void *)((double *)(work->data) + delta),
             m * n * sizeof(double));
      memcpy(epsA->data, (void *)((double *)(work->data) + eps),
             nq * n * sizeof(double));
      memcpy(xplusA->data, (void *)((double *)(work->data) + xplus),
             m * n * sizeof(double));
      memcpy(fnA->data, (void *)((double *)(work->data) + fn),
             nq * n * sizeof(double));

      res_var = *((double *)(work->data) + rvar);
      sum_square = *((double *)(work->data) + wss);
      sum_square_delta = *((double *)(work->data) + wssde);
      sum_square_eps = *((double *)(work->data) + wssep);
      inv_condnum = *((double *)(work->data) + rcond);
      rel_error = *((double *)(work->data) + eta);

      retobj =
        Py_BuildValue
        ("OOO{s:O,s:O,s:O,s:O,s:d,s:d,s:d,s:d,s:d,s:d,s:O,s:O,s:O,s:i}",
         PyArray_Return(beta), PyArray_Return(sd_beta),
         PyArray_Return(cov_beta), "delta", PyArray_Return(deltaA), "eps",
         PyArray_Return(epsA), "xplus", PyArray_Return(xplusA), "y",
         PyArray_Return(fnA), "res_var", res_var, "sum_square", sum_square,
         "sum_square_delta", sum_square_delta, "sum_square_eps",
         sum_square_eps, "inv_condnum", inv_condnum, "rel_error", rel_error,
         "work", PyArray_Return(work), "work_ind", work_ind, "iwork",
         PyArray_Return(iwork), "info", info);
      Py_DECREF((PyObject *) sd_beta);
      Py_DECREF((PyObject *) cov_beta);
      Py_DECREF((PyObject *) deltaA);
      Py_DECREF((PyObject *) epsA);
      Py_DECREF((PyObject *) xplusA);
      Py_DECREF((PyObject *) fnA);
      Py_DECREF((PyObject *) work_ind);

      return retobj;
    }
}

PyObject *odr(PyObject * self, PyObject * args, PyObject * kwds)
{
  PyObject *fcn, *initbeta, *py, *px, *pwe = NULL, *pwd = NULL, *fjacb = NULL;
  PyObject *fjacd = NULL, *pifixb = NULL, *pifixx = NULL;
  PyObject *pstpb = NULL, *pstpd = NULL, *psclb = NULL, *pscld = NULL;
  PyObject *pwork = NULL, *piwork = NULL, *extra_args = NULL;
  int job = 0, ndigit = 0, maxit = -1, iprint = 0;
  int full_output = 0;
  double taufac = 0.0, sstol = -1.0, partol = -1.0;
  char *errfile = NULL, *rptfile = NULL;
  int lerrfile = 0, lrptfile = 0;
  PyArrayObject *beta = NULL, *y = NULL, *x = NULL, *we = NULL, *wd = NULL;
  PyArrayObject *ifixb = NULL, *ifixx = NULL;
  PyArrayObject *stpb = NULL, *stpd = NULL, *sclb = NULL, *scld = NULL;
  PyArrayObject *work = NULL, *iwork = NULL;
  int n, m, np, nq, ldy, ldx, ldwe, ld2we, ldwd, ld2wd, ldifx;
  int lunerr = -1, lunrpt = -1, ldstpd, ldscld, lwork, liwork, info = 0;
  static char *kw_list[] = { "fcn", "initbeta", "y", "x", "we", "wd", "fjacb",
    "fjacd", "extra_args", "ifixb", "ifixx", "job", "iprint", "errfile",
    "rptfile", "ndigit", "taufac", "sstol", "partol",
    "maxit", "stpb", "stpd", "sclb", "scld", "work",
    "iwork", "full_output", NULL
  };
  int isodr = 1;
  PyObject *result;
  npy_intp dim1[1], dim2[2], dim3[3];
  int implicit;                 /* flag for implicit model */


  if (kwds == NULL)
    {
      if (!PyArg_ParseTuple(args, "OOOO|OOOOOOOiiz#z#idddiOOOOOOi:odr",
                            &fcn, &initbeta, &py, &px, &pwe, &pwd,
                            &fjacb, &fjacd, &extra_args, &pifixb, &pifixx,
                            &job, &iprint, &errfile, &lerrfile, &rptfile,
                            &lrptfile, &ndigit, &taufac, &sstol, &partol,
                            &maxit, &pstpb, &pstpd, &psclb, &pscld, &pwork,
                            &piwork, &full_output))
        {
          return NULL;
        }
    }
  else
    {
      if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                       "OOOO|OOOOOOOiiz#z#idddiOOOOOOi:odr",
                                       kw_list, &fcn, &initbeta, &py, &px,
                                       &pwe, &pwd, &fjacb, &fjacd,
                                       &extra_args, &pifixb, &pifixx, &job,
                                       &iprint, &errfile, &lerrfile, &rptfile,
                                       &lrptfile, &ndigit, &taufac, &sstol,
                                       &partol, &maxit, &pstpb, &pstpd,
                                       &psclb, &pscld, &pwork, &piwork,
                                       &full_output))
        {
          return NULL;
        }
    }

  /* Check the validity of all arguments */

  if (!PyCallable_Check(fcn))
    {
      PYERR(PyExc_TypeError, "fcn must be callable");
    }
  if (!PySequence_Check(initbeta))
    {
      PYERR(PyExc_TypeError, "initbeta must be a sequence");
    }
  if (!PySequence_Check(py))
    {
      /* Checking whether py is an int 
       *
       * XXX: PyInt_Check for np.int32 instances does not work on python 2.6 -
       * we should fix this in numpy, workaround by trying to cast to an int
       * for now */
      long val;

      PyErr_Clear();
      val = PyInt_AsLong(py);
      if (val == -1 && PyErr_Occurred()) {
        PYERR(PyExc_TypeError,
              "y must be a sequence or integer (if model is implicit)");
      }
    }
  if (!PySequence_Check(px))
    {
      PYERR(PyExc_TypeError, "x must be a sequence");
    }
  if (pwe != NULL && !PySequence_Check(pwe) && !PyNumber_Check(pwe))
    {
      PYERR(PyExc_TypeError, "we must be a sequence or a number");
    }
  if (pwd != NULL && !PySequence_Check(pwd) && !PyNumber_Check(pwd))
    {
      PYERR(PyExc_TypeError, "wd must be a sequence or a number");
    }
  if (fjacb != NULL && !PyCallable_Check(fjacb))
    {
      PYERR(PyExc_TypeError, "fjacb must be callable");
    }
  if (fjacd != NULL && !PyCallable_Check(fjacd))
    {
      PYERR(PyExc_TypeError, "fjacd must be callable");
    }
  if (extra_args != NULL && !PySequence_Check(extra_args))
    {
      PYERR(PyExc_TypeError, "extra_args must be a sequence");
    }
  if (pifixx != NULL && !PySequence_Check(pifixx))
    {
      PYERR(PyExc_TypeError, "ifixx must be a sequence");
    }
  if (pifixb != NULL && !PySequence_Check(pifixb))
    {
      PYERR(PyExc_TypeError, "ifixb must be a sequence");
    }
  if (pstpb != NULL && !PySequence_Check(pstpb))
    {
      PYERR(PyExc_TypeError, "stpb must be a sequence");
    }
  if (pstpd != NULL && !PySequence_Check(pstpd))
    {
      PYERR(PyExc_TypeError, "stpd must be a sequence");
    }
  if (psclb != NULL && !PySequence_Check(psclb))
    {
      PYERR(PyExc_TypeError, "sclb must be a sequence");
    }
  if (pscld != NULL && !PySequence_Check(pscld))
    {
      PYERR(PyExc_TypeError, "scld must be a sequence");
    }
  if (pwork != NULL && !PyArray_Check(pwork))
    {
      PYERR(PyExc_TypeError, "work must be an array");
    }
  if (piwork != NULL && !PyArray_Check(piwork))
    {
      PYERR(PyExc_TypeError, "iwork must be an array");
    }

  /* start processing the arguments and check for errors on the way */

  /* check for implicit model */

  implicit = (job % 10 == 1);

  if (!implicit)
    {
      if ((y =
           (PyArrayObject *) PyArray_CopyFromObject(py, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "y could not be made into a suitable array");
        }
      n = y->dimensions[y->nd - 1];     /* pick the last dimension */
      if ((x =
           (PyArrayObject *) PyArray_CopyFromObject(px, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "x could not be made into a suitable array");
        }
      if (n != x->dimensions[x->nd - 1])
        {
          PYERR(PyExc_ValueError,
                "x and y don't have matching numbers of observations");
        }
      if (y->nd == 1)
        {
          nq = 1;
        }
      else
        {
          nq = y->dimensions[0];
        }

      ldx = ldy = n;
    }
  else
    {                           /* we *do* have an implicit model */
      ldy = 1;
      nq = (int)PyInt_AsLong(py);
      dim1[0] = 1;

      /* initialize y to a dummy array; never referenced */
      y = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);

      if ((x =
           (PyArrayObject *) PyArray_CopyFromObject(px, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "x could not be made into a suitable array");
        }

      n = x->dimensions[x->nd - 1];
      ldx = n;
    }

  if (x->nd == 1)
    {
      m = 1;
    }
  else
    {
      m = x->dimensions[0];
    }                           /* x, y */

  if ((beta =
       (PyArrayObject *) PyArray_CopyFromObject(initbeta, NPY_DOUBLE, 1,
                                                1)) == NULL)
    {
      PYERR(PyExc_ValueError,
            "initbeta could not be made into a suitable array");
    }
  np = beta->dimensions[0];

  if (pwe == NULL)
    {
      ldwe = ld2we = 1;
      dim1[0] = n;
      we = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      ((double *)(we->data))[0] = -1.0;
    }
  else if (PyNumber_Check(pwe) && !PyArray_Check(pwe))
    {
      /* we is a single weight, set the first value of we to -pwe */
      PyObject *tmp;
      double val;

      tmp = PyNumber_Float(pwe);
      if (tmp == NULL)
        PYERR(PyExc_ValueError, "could not convert we to a suitable array");
      val = PyFloat_AsDouble(tmp);
      Py_DECREF(tmp);

      dim3[0] = nq;
      dim3[1] = 1;
      dim3[2] = 1;
      we = (PyArrayObject *) PyArray_SimpleNew(3, dim3, NPY_DOUBLE);
      if (implicit)
        {
          ((double *)(we->data))[0] = val;
        }
      else
        {
          ((double *)(we->data))[0] = -val;
        }
      ldwe = ld2we = 1;
    }
  else if (PySequence_Check(pwe))
    {
      /* we needs to be turned into an array */

      if ((we =
           (PyArrayObject *) PyArray_CopyFromObject(pwe, NPY_DOUBLE, 1,
                                                    3)) == NULL)
        {
          PYERR(PyExc_ValueError, "could not convert we to a suitable array");
        }

      if (we->nd == 1 && nq == 1)
        {

          ldwe = n;
          ld2we = 1;
        }
      else if (we->nd == 1 && we->dimensions[0] == nq)
        {
          /* we is a rank-1 array with diagonal weightings to be broadcast 
           * to all observations */
          ldwe = 1;
          ld2we = 1;
        }
      else if (we->nd == 3 && we->dimensions[0] == nq
               && we->dimensions[1] == nq && we->dimensions[2] == 1)
        {
          /* we is a rank-3 array with the covariant weightings 
             to be broadcast to all observations */
          ldwe = 1;
          ld2we = nq;
        }
      else if (we->nd == 2 && we->dimensions[0] == nq
               && we->dimensions[1] == nq)
        {
          /* we is a rank-2 array with the full covariant weightings 
             to be broadcast to all observations */
          ldwe = 1;
          ld2we = nq;
        }

      else if (we->nd == 2 && we->dimensions[0] == nq
               && we->dimensions[1] == n)
        {
          /* we is a rank-2 array with the diagonal elements of the 
             covariant weightings for each observation */
          ldwe = n;
          ld2we = 1;
        }
      else if (we->nd == 3 && we->dimensions[0] == nq
               && we->dimensions[1] == nq && we->dimensions[2] == n)
        {
          /* we is the full specification of the covariant weights
             for each observation */
          ldwe = n;
          ld2we = nq;
        }
      else
        {
          PYERR(PyExc_ValueError, "could not convert we to a suitable array");
        }
    }                           /* we */

  if (pwd == NULL)
    {
      ldwd = ld2wd = 1;

      dim1[0] = m;
      wd = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      ((double *)(wd->data))[0] = -1.0;
    }
  else if (PyNumber_Check(pwd) && !PyArray_Check(pwd))
    {
      /* wd is a single weight, set the first value of wd to -pwd */
      PyObject *tmp;
      double val;

      tmp = PyNumber_Float(pwd);
      if (tmp == NULL)
        PYERR(PyExc_ValueError, "could not convert wd to a suitable array");
      val = PyFloat_AsDouble(tmp);
      Py_DECREF(tmp);

      dim3[0] = 1;
      dim3[1] = 1;
      dim3[2] = m;
      wd = (PyArrayObject *) PyArray_SimpleNew(3, dim3, NPY_DOUBLE);
      ((double *)(wd->data))[0] = -val;
      ldwd = ld2wd = 1;
    }
  else if (PySequence_Check(pwd))
    {
      /* wd needs to be turned into an array */

      if ((wd =
           (PyArrayObject *) PyArray_CopyFromObject(pwd, NPY_DOUBLE, 1,
                                                    3)) == NULL)
        {
          PYERR(PyExc_ValueError, "could not convert wd to a suitable array");
        }

      if (wd->nd == 1 && m == 1)
        {
          ldwd = n;
          ld2wd = 1;
        }
      else if (wd->nd == 1 && wd->dimensions[0] == m)
        {
          /* wd is a rank-1 array with diagonal weightings to be broadcast 
           * to all observations */
          ldwd = 1;
          ld2wd = 1;
        }

      else if (wd->nd == 3 && wd->dimensions[0] == m
               && wd->dimensions[1] == m && wd->dimensions[2] == 1)
        {
          /* wd is a rank-3 array with the covariant wdightings 
             to be broadcast to all observations */
          ldwd = 1;
          ld2wd = m;
        }
      else if (wd->nd == 2 && wd->dimensions[0] == m
               && wd->dimensions[1] == m)
        {
          /* wd is a rank-2 array with the full covariant weightings 
             to be broadcast to all observations */
          ldwd = 1;
          ld2wd = m;
        }

      else if (wd->nd == 2 && wd->dimensions[0] == m
               && wd->dimensions[1] == n)
        {
          /* wd is a rank-2 array with the diagonal elements of the 
             covariant weightings for each observation */
          ldwd = n;
          ld2wd = 1;
        }
      else if (wd->nd == 3 && wd->dimensions[0] == m
               && wd->dimensions[1] == m && wd->dimensions[2] == n)
        {
          /* wd is the full specification of the covariant weights
             for each observation */
          ldwd = n;
          ld2wd = m;
        }
      else
        {
          PYERR(PyExc_ValueError, "could not convert wd to a suitable array");
        }

    }                           /* wd */


  if (pifixb == NULL)
    {
      dim1[0] = np;
      ifixb = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_INT);
      *(int *)(ifixb->data) = -1;      /* set first element negative */
    }
  else
    {
      /* pifixb is a sequence as checked before */

      if ((ifixb =
           (PyArrayObject *) PyArray_CopyFromObject(pifixb, NPY_INT, 1,
                                                    1)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert ifixb to a suitable array");
        }

      if (ifixb->dimensions[0] != np)
        {
          PYERR(PyExc_ValueError,
                "could not convert ifixb to a suitable array");
        }
    }                           /* ifixb */

  if (pifixx == NULL)
    {
      dim2[0] = m;
      dim2[1] = 1;
      ifixx = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_INT);
      *(int *)(ifixx->data) = -1;      /* set first element negative */
      ldifx = 1;
    }
  else
    {
      /* pifixx is a sequence as checked before */

      if ((ifixx =
           (PyArrayObject *) PyArray_CopyFromObject(pifixx, NPY_INT, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert ifixx to a suitable array");
        }

      if (ifixx->nd == 1 && ifixx->dimensions[0] == m)
        {
          ldifx = 1;
        }
      else if (ifixx->nd == 1 && ifixx->dimensions[0] == n && m == 1)
        {
          ldifx = n;
        }
      else if (ifixx->nd == 2 && ifixx->dimensions[0] == m
               && ifixx->dimensions[1] == n)
        {
          ldifx = n;
        }
      else
        {
          PYERR(PyExc_ValueError,
                "could not convert ifixx to a suitable array");
        }
    }                           /* ifixx */

  if (errfile != NULL)
    {
      /* call FORTRAN's OPEN to open the file with a logical unit of 18 */
      lunerr = 18;
      F_FUNC(dluno,DLUNO)(&lunerr, errfile, lerrfile);
    }

  if (rptfile != NULL)
    {
      /* call FORTRAN's OPEN to open the file with a logical unit of 19 */
      lunrpt = 19;
      F_FUNC(dluno,DLUNO)(&lunrpt, rptfile, lrptfile);
    }

  if (pstpb == NULL)
    {
      dim1[0] = np;
      stpb = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      *(double *)(stpb->data) = 0.0;
    }
  else                          /* pstpb is a sequence */
    {
      if ((stpb =
           (PyArrayObject *) PyArray_CopyFromObject(pstpb, NPY_DOUBLE, 1,
                                                    1)) == NULL
          || stpb->dimensions[0] != np)
        {
          PYERR(PyExc_ValueError,
                "could not convert stpb to a suitable array");
        }
    }                           /* stpb */

  if (pstpd == NULL)
    {
      dim2[0] = 1;
      dim2[1] = m;
      stpd = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
      *(double *)(stpd->data) = 0.0;
      ldstpd = 1;
    }
  else
    {
      if ((stpd =
           (PyArrayObject *) PyArray_CopyFromObject(pstpd, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert stpb to a suitable array");
        }

      if (stpd->nd == 1 && stpd->dimensions[0] == m)
        {
          ldstpd = 1;
        }
      else if (stpd->nd == 1 && stpd->dimensions[0] == n && m == 1)
        {
          ldstpd = n;
        }
      else if (stpd->nd == 2 && stpd->dimensions[0] == n &&
               stpd->dimensions[1] == m)
        {
          ldstpd = n;
        }
    }                           /* stpd */

  if (psclb == NULL)
    {
      dim1[0] = np;
      sclb = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      *(double *)(sclb->data) = 0.0;
    }
  else                          /* psclb is a sequence */
    {
      if ((sclb =
           (PyArrayObject *) PyArray_CopyFromObject(psclb, NPY_DOUBLE, 1,
                                                    1)) == NULL
          || sclb->dimensions[0] != np)
        {
          PYERR(PyExc_ValueError,
                "could not convert sclb to a suitable array");
        }
    }                           /* sclb */

  if (pscld == NULL)
    {
      dim2[0] = 1;
      dim2[1] = n;
      scld = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
      *(double *)(scld->data) = 0.0;
      ldscld = 1;
    }
  else
    {
      if ((scld =
           (PyArrayObject *) PyArray_CopyFromObject(pscld, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert stpb to a suitable array");
        }

      if (scld->nd == 1 && scld->dimensions[0] == m)
        {
          ldscld = 1;
        }
      else if (scld->nd == 1 && scld->dimensions[0] == n && m == 1)
        {
          ldscld = n;
        }
      else if (scld->nd == 2 && scld->dimensions[0] == n &&
               scld->dimensions[1] == m)
        {
          ldscld = n;
        }
    }                           /* scld */

  if (job % 10 < 2)
    {
      /* ODR, not OLS */

      lwork =
        18 + 11 * np + np * np + m + m * m + 4 * n * nq + 6 * n * m +
        2 * n * nq * np + 2 * n * nq * m + nq * nq + 5 * nq + nq * (np + m) +
        ldwe * ld2we * nq;

      isodr = 1;
    }
  else
    {
      /* OLS, not ODR */

      lwork =
        18 + 11 * np + np * np + m + m * m + 4 * n * nq + 2 * n * m +
        2 * n * nq * np + 5 * nq + nq * (np + m) + ldwe * ld2we * nq;

      isodr = 0;
    }

  liwork = 20 + np + nq * (np + m);

  if ((job / 10000) % 10 >= 1)
    {
      /* fit is a restart, make sure work and iwork are input */

      if (pwork == NULL || piwork == NULL)
        {
          PYERR(PyExc_ValueError,
                "need to input work and iwork arrays to restart");
        }
    }

  if ((job / 1000) % 10 >= 1)
    {
      /* delta should be supplied, make sure the user does */

      if (pwork == NULL)
        {
          PYERR(PyExc_ValueError,
                "need to input work array for delta initialization");
        }
    }

  if (pwork != NULL)
    {
      if ((work =
           (PyArrayObject *) PyArray_CopyFromObject(pwork, NPY_DOUBLE, 1,
                                                    1)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert work to a suitable array");
        }
      if (work->dimensions[0] < lwork)
        {
            printf("%d %d\n", work->dimensions[0], lwork);
          PYERR(PyExc_ValueError, "work is too small");
        }
    }
  else
    {
      /* initialize our own work array */
      dim1[0] = lwork;
      work = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
    }                           /* work */

  if (piwork != NULL)
    {
      if ((iwork =
           (PyArrayObject *) PyArray_CopyFromObject(piwork, NPY_INT, 1,
                                                    1)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "could not convert iwork to a suitable array");
        }

      if (iwork->dimensions[0] < liwork)
        {
          PYERR(PyExc_ValueError, "iwork is too small");
        }
    }
  else
    {
      /* initialize our own iwork array */
      dim1[0] = liwork;
      iwork = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_INT);
    }                           /* iwork */

  /* check if what JOB requests can be done with what the user has 
     input into the function */

  if ((job / 10) % 10 >= 2)
    {
      /* derivatives are supposed to be supplied */

      if (fjacb == NULL || fjacd == NULL)
        {
          PYERR(PyExc_ValueError,
                "need fjacb and fjacd to calculate derivatives");
        }
    }

  /* setup the global data for the callback */
  odr_global.fcn = fcn;
  Py_INCREF(fcn);
  odr_global.fjacb = fjacb;
  Py_XINCREF(fjacb);
  odr_global.fjacd = fjacd;
  Py_XINCREF(fjacd);
  odr_global.pyBeta = (PyObject *) beta;
  Py_INCREF(beta);
  odr_global.extra_args = extra_args;
  Py_XINCREF(extra_args);

  /* now call DODRC */
  F_FUNC(dodrc,DODRC)(fcn_callback, &n, &m, &np, &nq, (double *)(beta->data),
        (double *)(y->data), &ldy, (double *)(x->data), &ldx,
        (double *)(we->data), &ldwe, &ld2we,
        (double *)(wd->data), &ldwd, &ld2wd,
        (int *)(ifixb->data), (int *)(ifixx->data), &ldifx,
        &job, &ndigit, &taufac, &sstol, &partol, &maxit,
        &iprint, &lunerr, &lunrpt,
        (double *)(stpb->data), (double *)(stpd->data), &ldstpd,
        (double *)(sclb->data), (double *)(scld->data), &ldscld,
        (double *)(work->data), &lwork, (int *)(iwork->data), &liwork,
        &info);

  result = gen_output(n, m, np, nq, ldwe, ld2we,
                      beta, work, iwork, isodr, info, full_output);

  if (result == NULL)
    PYERR(PyExc_RuntimeError, "could not generate output");

  if (lunerr != -1)
    {
      F_FUNC(dlunc,DLUNC)(&lunerr);
    }
  if (lunrpt != -1)
    {
      F_FUNC(dlunc,DLUNC)(&lunrpt);
    }

  Py_DECREF(odr_global.fcn);
  Py_XDECREF(odr_global.fjacb);
  Py_XDECREF(odr_global.fjacd);
  Py_DECREF(odr_global.pyBeta);
  Py_XDECREF(odr_global.extra_args);

  odr_global.fcn = odr_global.fjacb = odr_global.fjacd = odr_global.pyBeta =
    odr_global.extra_args = NULL;

  Py_DECREF(beta);
  Py_DECREF(y);
  Py_DECREF(x);
  Py_DECREF(we);
  Py_DECREF(wd);
  Py_DECREF(ifixb);
  Py_DECREF(ifixx);
  Py_DECREF(stpb);
  Py_DECREF(stpd);
  Py_DECREF(sclb);
  Py_DECREF(scld);
  Py_DECREF(work);
  Py_DECREF(iwork);

  return result;

fail:


  if (lunerr != -1)
    {
      F_FUNC(dlunc,DLUNC)(&lunerr);
    }
  if (lunrpt != -1)
    {
      F_FUNC(dlunc,DLUNC)(&lunrpt);
    }

  Py_XDECREF(beta);
  Py_XDECREF(y);
  Py_XDECREF(x);
  Py_XDECREF(we);
  Py_XDECREF(wd);
  Py_XDECREF(ifixb);
  Py_XDECREF(ifixx);
  Py_XDECREF(stpb);
  Py_XDECREF(stpd);
  Py_XDECREF(sclb);
  Py_XDECREF(scld);
  Py_XDECREF(work);
  Py_XDECREF(iwork);

  return NULL;
}

static void check_args(int n, int m, int np, int nq,
                       PyArrayObject * beta,
                       PyArrayObject * y, int ldy,
                       PyArrayObject * x, int ldx,
                       PyArrayObject * we, int ldwe, int ld2we,
                       PyArrayObject * wd, int ldwd, int ld2wd,
                       PyArrayObject * ifixb, PyArrayObject * ifixx,
                       int ldifx, int job, int ndigit, double taufac,
                       double sstol, double partol, int maxit,
                       PyArrayObject * stpb, PyArrayObject * stpd,
                       int ldstpd, PyArrayObject * sclb,
                       PyArrayObject * scld, int ldscld,
                       PyArrayObject * work, int lwork,
                       PyArrayObject * iwork, int liwork, int info)
{
  PyObject *printdict;

  printdict =
    Py_BuildValue
    ("{s:i,s:i,s:i,s:i,s:O,s:O,s:i,s:O,s:i,s:O,s:i,s:i,s:O,s:i,s:i,s:O,s:O,s:i,s:i,s:i,s:d,s:d,s:d,s:i,s:O,s:O,s:i,s:O,s:O,s:i,s:O,s:i,s:O,s:i,s:i}",
     "n", n, "m", m, "np", np, "nq", nq, "beta", (PyObject *) beta, "y",
     (PyObject *) y, "ldy", ldy, "x", (PyObject *) x, "ldx", ldx, "we",
     (PyObject *) we, "ldwe", ldwe, "ld2we", ld2we, "wd", (PyObject *) wd,
     "ldwd", ldwd, "ld2wd", ld2wd, "ifixb", (PyObject *) ifixb, "ifixx",
     (PyObject *) ifixx, "ldifx", ldifx, "job", job, "ndigit", ndigit,
     "taufac", taufac, "sstol", sstol, "partol", partol, "maxit", maxit,
     "stpb", (PyObject *) stpb, "stpd", (PyObject *) stpd, "ldstpd", ldstpd,
     "sclb", (PyObject *) sclb, "scld", (PyObject *) scld, "ldscld", ldscld,
     "work", (PyObject *) work, "lwork", lwork, "iwork", (PyObject *) iwork,
     "liwork", liwork, "info", info);
  if (printdict == NULL)
    {
      PyErr_Print();
      return;
    }

  PyObject_Print(printdict, stdout, Py_PRINT_RAW);
  printf("\n");
  Py_XDECREF(printdict);
}

static char odr__doc__[] =
  "odr(fcn, beta0, y, x,\nwe=None, wd=None, fjacb=None, fjacd=None,\nextra_args=None, ifixx=None, ifixb=None, job=0, iprint=0,\nerrfile=None, rptfile=None, ndigit=0,\ntaufac=0.0, sstol=-1.0, partol=-1.0,\nmaxit=-1, stpb=None, stpd=None,\nsclb=None, scld=None, work=None, iwork=None,\nfull_output=0)";

static PyMethodDef methods[] = {
  {"odr", (PyCFunction) odr, METH_VARARGS | METH_KEYWORDS, odr__doc__},
  {NULL, NULL},
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_odrpack",
    NULL,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit___odrpack(void)
{
    PyObject *m, *s, *d;

    m = PyModule_Create(&moduledef);
    import_array();

    d = PyModule_GetDict(m);
    odr_error = PyErr_NewException("odr.odrpack.odr_error", NULL, NULL);
    odr_stop = PyErr_NewException("odr.odrpack.odr_stop", NULL, NULL);
    PyDict_SetItemString(d, "odr_error", odr_error);
    PyDict_SetItemString(d, "odr_stop", odr_stop);

    return m;
}
#else
PyMODINIT_FUNC init__odrpack(void)
{
  PyObject *m, *d;

  import_array();

  m = Py_InitModule("__odrpack", methods);
  d = PyModule_GetDict(m);
  odr_error = PyErr_NewException("odr.odrpack.odr_error", NULL, NULL);
  odr_stop = PyErr_NewException("odr.odrpack.odr_stop", NULL, NULL);
  PyDict_SetItemString(d, "odr_error", odr_error);
  PyDict_SetItemString(d, "odr_stop", odr_stop);
}
#endif
