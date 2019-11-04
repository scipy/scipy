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

#define PY_SSIZE_T_CLEAN
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
  PyObject *result = NULL;
  PyArrayObject *result_array = NULL;
  PyArrayObject *pyXplusD;
  void *beta_dst;

  arg01 = PyTuple_New(2);

  if (*m != 1)
    {
      npy_intp dim2[2];
      dim2[0] = *m;
      dim2[1] = *n;
      pyXplusD = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);
      memcpy(PyArray_DATA(pyXplusD), (void *)xplusd, (*m) * (*n) * sizeof(double));
    }
  else
    {
      npy_intp dim1[1];
      dim1[0] = *n;
      pyXplusD = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      memcpy(PyArray_DATA(pyXplusD), (void *)xplusd, (*n) * sizeof(double));
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

  beta_dst = (PyArray_DATA((PyArrayObject *) odr_global.pyBeta));
  if (beta != beta_dst) {
      memcpy(beta_dst, (void *)beta, (*np) * sizeof(double));
  }

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
          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }
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

      memcpy((void *)f, PyArray_DATA(result_array), (*n) * (*nq) * sizeof(double));
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
          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }
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

          if (PyArray_NDIM(result_array) != 3)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "Beta Jacobian is not rank-3");
            }
        }
      else if (*nq == 1)
        {
          /* result_array should be rank-2 */

          if (PyArray_NDIM(result_array) != 2)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "Beta Jacobian is not rank-2");
            }
        }

      memcpy((void *)fjacb, PyArray_DATA(result_array),
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
          if (PyErr_ExceptionMatches(odr_stop))
            {
              /* stop, don't fail */
              *istop = 1;

              Py_DECREF(arglist);
              return;
            }
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

          if (PyArray_NDIM(result_array) != 3)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-3");
            }
        }
      else if (*nq == 1 && *m != 1)
        {
          /* result_array should be rank-2 */

          if (PyArray_NDIM(result_array) != 2)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-2");
            }
        }
      else if (*nq == 1 && *m == 1)
        {
          /* result_array should be rank-1 */

          if (PyArray_NDIM(result_array) != 1)
            {
              Py_DECREF(result_array);
              PYERR2(odr_error, "xplusd Jacobian is not rank-1");
            }
        }

      memcpy((void *)fjacd, PyArray_DATA(result_array),
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
      /* fatal error in fcn call; return NULL to propagate the exception */

      return NULL;
  }

  lwkmn = PyArray_DIMS(work)[0];

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

  dim1[0] = PyArray_DIMS(beta)[0];
  sd_beta = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
  dim2[0] = PyArray_DIMS(beta)[0];
  dim2[1] = PyArray_DIMS(beta)[0];
  cov_beta = (PyArrayObject *) PyArray_SimpleNew(2, dim2, NPY_DOUBLE);

  memcpy(PyArray_DATA(sd_beta), (void *)((double *)(PyArray_DATA(work)) + sd),
         np * sizeof(double));
  memcpy(PyArray_DATA(cov_beta), (void *)((double *)(PyArray_DATA(work)) + vcv),
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

      memcpy(PyArray_DATA(deltaA), (void *)((double *)(PyArray_DATA(work)) + delta),
             m * n * sizeof(double));
      memcpy(PyArray_DATA(epsA), (void *)((double *)(PyArray_DATA(work)) + eps),
             nq * n * sizeof(double));
      memcpy(PyArray_DATA(xplusA), (void *)((double *)(PyArray_DATA(work)) + xplus),
             m * n * sizeof(double));
      memcpy(PyArray_DATA(fnA), (void *)((double *)(PyArray_DATA(work)) + fn),
             nq * n * sizeof(double));

      res_var = *((double *)(PyArray_DATA(work)) + rvar);
      sum_square = *((double *)(PyArray_DATA(work)) + wss);
      sum_square_delta = *((double *)(PyArray_DATA(work)) + wssde);
      sum_square_eps = *((double *)(PyArray_DATA(work)) + wssep);
      inv_condnum = *((double *)(PyArray_DATA(work)) + rcond);
      rel_error = *((double *)(PyArray_DATA(work)) + eta);

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
  Py_ssize_t lerrfile = 0, lrptfile = 0;
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
      n = PyArray_DIMS(y)[PyArray_NDIM(y) - 1];     /* pick the last dimension */
      if ((x =
           (PyArrayObject *) PyArray_CopyFromObject(px, NPY_DOUBLE, 1,
                                                    2)) == NULL)
        {
          PYERR(PyExc_ValueError,
                "x could not be made into a suitable array");
        }
      if (n != PyArray_DIMS(x)[PyArray_NDIM(x) - 1])
        {
          PYERR(PyExc_ValueError,
                "x and y don't have matching numbers of observations");
        }
      if (PyArray_NDIM(y) == 1)
        {
          nq = 1;
        }
      else
        {
          nq = PyArray_DIMS(y)[0];
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

      n = PyArray_DIMS(x)[PyArray_NDIM(x) - 1];
      ldx = n;
    }

  if (PyArray_NDIM(x) == 1)
    {
      m = 1;
    }
  else
    {
      m = PyArray_DIMS(x)[0];
    }                           /* x, y */

  if ((beta =
       (PyArrayObject *) PyArray_CopyFromObject(initbeta, NPY_DOUBLE, 1,
                                                1)) == NULL)
    {
      PYERR(PyExc_ValueError,
            "initbeta could not be made into a suitable array");
    }
  np = PyArray_DIMS(beta)[0];

  if (pwe == NULL)
    {
      ldwe = ld2we = 1;
      dim1[0] = n;
      we = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      ((double *)(PyArray_DATA(we)))[0] = -1.0;
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
          ((double *)(PyArray_DATA(we)))[0] = val;
        }
      else
        {
          ((double *)(PyArray_DATA(we)))[0] = -val;
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

      if (PyArray_NDIM(we) == 1 && nq == 1)
        {

          ldwe = n;
          ld2we = 1;
        }
      else if (PyArray_NDIM(we) == 1 && PyArray_DIMS(we)[0] == nq)
        {
          /* we is a rank-1 array with diagonal weightings to be broadcast 
           * to all observations */
          ldwe = 1;
          ld2we = 1;
        }
      else if (PyArray_NDIM(we) == 3 && PyArray_DIMS(we)[0] == nq
               && PyArray_DIMS(we)[1] == nq && PyArray_DIMS(we)[2] == 1)
        {
          /* we is a rank-3 array with the covariant weightings 
             to be broadcast to all observations */
          ldwe = 1;
          ld2we = nq;
        }
      else if (PyArray_NDIM(we) == 2 && PyArray_DIMS(we)[0] == nq
               && PyArray_DIMS(we)[1] == nq)
        {
          /* we is a rank-2 array with the full covariant weightings 
             to be broadcast to all observations */
          ldwe = 1;
          ld2we = nq;
        }

      else if (PyArray_NDIM(we) == 2 && PyArray_DIMS(we)[0] == nq
               && PyArray_DIMS(we)[1] == n)
        {
          /* we is a rank-2 array with the diagonal elements of the 
             covariant weightings for each observation */
          ldwe = n;
          ld2we = 1;
        }
      else if (PyArray_NDIM(we) == 3 && PyArray_DIMS(we)[0] == nq
               && PyArray_DIMS(we)[1] == nq && PyArray_DIMS(we)[2] == n)
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
      ((double *)(PyArray_DATA(wd)))[0] = -1.0;
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
      ((double *)(PyArray_DATA(wd)))[0] = -val;
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

      if (PyArray_NDIM(wd) == 1 && m == 1)
        {
          ldwd = n;
          ld2wd = 1;
        }
      else if (PyArray_NDIM(wd) == 1 && PyArray_DIMS(wd)[0] == m)
        {
          /* wd is a rank-1 array with diagonal weightings to be broadcast 
           * to all observations */
          ldwd = 1;
          ld2wd = 1;
        }

      else if (PyArray_NDIM(wd) == 3 && PyArray_DIMS(wd)[0] == m
               && PyArray_DIMS(wd)[1] == m && PyArray_DIMS(wd)[2] == 1)
        {
          /* wd is a rank-3 array with the covariant wdightings 
             to be broadcast to all observations */
          ldwd = 1;
          ld2wd = m;
        }
      else if (PyArray_NDIM(wd) == 2 && PyArray_DIMS(wd)[0] == m
               && PyArray_DIMS(wd)[1] == m)
        {
          /* wd is a rank-2 array with the full covariant weightings 
             to be broadcast to all observations */
          ldwd = 1;
          ld2wd = m;
        }

      else if (PyArray_NDIM(wd) == 2 && PyArray_DIMS(wd)[0] == m
               && PyArray_DIMS(wd)[1] == n)
        {
          /* wd is a rank-2 array with the diagonal elements of the 
             covariant weightings for each observation */
          ldwd = n;
          ld2wd = 1;
        }
      else if (PyArray_NDIM(wd) == 3 && PyArray_DIMS(wd)[0] == m
               && PyArray_DIMS(wd)[1] == m && PyArray_DIMS(wd)[2] == n)
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
      *(int *)(PyArray_DATA(ifixb)) = -1;      /* set first element negative */
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

      if (PyArray_DIMS(ifixb)[0] != np)
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
      *(int *)(PyArray_DATA(ifixx)) = -1;      /* set first element negative */
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

      if (PyArray_NDIM(ifixx) == 1 && PyArray_DIMS(ifixx)[0] == m)
        {
          ldifx = 1;
        }
      else if (PyArray_NDIM(ifixx) == 1 && PyArray_DIMS(ifixx)[0] == n && m == 1)
        {
          ldifx = n;
        }
      else if (PyArray_NDIM(ifixx) == 2 && PyArray_DIMS(ifixx)[0] == m
               && PyArray_DIMS(ifixx)[1] == n)
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
      *(double *)(PyArray_DATA(stpb)) = 0.0;
    }
  else                          /* pstpb is a sequence */
    {
      if ((stpb =
           (PyArrayObject *) PyArray_CopyFromObject(pstpb, NPY_DOUBLE, 1,
                                                    1)) == NULL
          || PyArray_DIMS(stpb)[0] != np)
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
      *(double *)(PyArray_DATA(stpd)) = 0.0;
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

      if (PyArray_NDIM(stpd) == 1 && PyArray_DIMS(stpd)[0] == m)
        {
          ldstpd = 1;
        }
      else if (PyArray_NDIM(stpd) == 1 && PyArray_DIMS(stpd)[0] == n && m == 1)
        {
          ldstpd = n;
        }
      else if (PyArray_NDIM(stpd) == 2 && PyArray_DIMS(stpd)[0] == n &&
               PyArray_DIMS(stpd)[1] == m)
        {
          ldstpd = n;
        }
    }                           /* stpd */

  if (psclb == NULL)
    {
      dim1[0] = np;
      sclb = (PyArrayObject *) PyArray_SimpleNew(1, dim1, NPY_DOUBLE);
      *(double *)(PyArray_DATA(sclb)) = 0.0;
    }
  else                          /* psclb is a sequence */
    {
      if ((sclb =
           (PyArrayObject *) PyArray_CopyFromObject(psclb, NPY_DOUBLE, 1,
                                                    1)) == NULL
          || PyArray_DIMS(sclb)[0] != np)
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
      *(double *)(PyArray_DATA(scld)) = 0.0;
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

      if (PyArray_NDIM(scld) == 1 && PyArray_DIMS(scld)[0] == m)
        {
          ldscld = 1;
        }
      else if (PyArray_NDIM(scld) == 1 && PyArray_DIMS(scld)[0] == n && m == 1)
        {
          ldscld = n;
        }
      else if (PyArray_NDIM(scld) == 2 && PyArray_DIMS(scld)[0] == n &&
               PyArray_DIMS(scld)[1] == m)
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
      if (PyArray_DIMS(work)[0] < lwork)
        {
            printf("%ld %d\n", PyArray_DIMS(work)[0], lwork);
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

      if (PyArray_DIMS(iwork)[0] < liwork)
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
   F_FUNC(dodrc,DODRC)(fcn_callback, &n, &m, &np, &nq, (double *)(PyArray_DATA(beta)),
         (double *)(PyArray_DATA(y)), &ldy, (double *)(PyArray_DATA(x)), &ldx,
         (double *)(PyArray_DATA(we)), &ldwe, &ld2we,
         (double *)(PyArray_DATA(wd)), &ldwd, &ld2wd,
         (int *)(PyArray_DATA(ifixb)), (int *)(PyArray_DATA(ifixx)), &ldifx,
         &job, &ndigit, &taufac, &sstol, &partol, &maxit,
         &iprint, &lunerr, &lunrpt,
         (double *)(PyArray_DATA(stpb)), (double *)(PyArray_DATA(stpd)), &ldstpd,
         (double *)(PyArray_DATA(sclb)), (double *)(PyArray_DATA(scld)), &ldscld,
         (double *)(PyArray_DATA(work)), &lwork, (int *)(PyArray_DATA(iwork)), &liwork,
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


PyObject *set_exceptions(PyObject * self, PyObject * args, PyObject * kwds)
{
    PyObject *exc_error, *exc_stop;

    if (!PyArg_ParseTuple(args, "OO", &exc_error, &exc_stop))
	return NULL;

    Py_INCREF(exc_stop);
    Py_INCREF(exc_error);
    odr_stop = exc_stop;
    odr_error = exc_error;

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef methods[] = {
  {"_set_exceptions", (PyCFunction) set_exceptions, METH_VARARGS, NULL},
  {"odr", (PyCFunction) odr, METH_VARARGS | METH_KEYWORDS, NULL},
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
    PyObject *m;
    import_array();
    m = PyModule_Create(&moduledef);
    return m;
}
#else
PyMODINIT_FUNC init__odrpack(void)
{
    PyObject *m;
    import_array();
    m = Py_InitModule("__odrpack", methods);
}
#endif
