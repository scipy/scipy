
/* MULTIPACK module by Travis Oliphant

Copyright (c) 2002 Travis Oliphant all rights reserved
oliphant.travis@ieee.org
Permission to use, modify, and distribute this software is given under the
terms of the SciPy (BSD style) license.  See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
*/


/* This extension module is a collection of wrapper functions around
common FORTRAN code in the package ODEPACK.

The wrappers are meant to be nearly direct translations between the
FORTAN code and Python.  Some parameters like sizes do not need to be
passed since they are available from the objects.

It is anticipated that a pure Python module be written to call these lower
level routines and make a simpler user interface.  All of the routines define
default values for little-used parameters so that even the raw routines are
quite useful without a separate wrapper.

FORTRAN Outputs that are not either an error indicator or the sought-after
results are placed in a dictionary and returned as an optional member of
the result tuple when the full_output argument is non-zero.
*/

#include "Python.h"
#include "numpy/npy_3kcompat.h"

#include "numpy/arrayobject.h"

#define PYERR(errobj,message) {\
    PyErr_SetString(errobj,message); \
    goto fail; \
}
#define PYERR2(errobj,message) { \
    PyErr_Print(); \
    PyErr_SetString(errobj, message); \
    goto fail; \
}

static PyObject *multipack_python_function = NULL;
static PyObject *multipack_python_jacobian = NULL;
static PyObject *multipack_extra_arguments = NULL;    /* a tuple */
static int multipack_jac_transpose = 1;
static int multipack_jac_type;


static
PyObject *call_python_function(PyObject *func, npy_intp n, double *x,
                               PyObject *args, PyObject *error_obj)
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

    /* Build sequence argument from inputs */
    sequence = (PyArrayObject *) PyArray_SimpleNewFromData(1, &n, NPY_DOUBLE,
                                                           (char *) x);
    if (sequence == NULL) {
        goto fail;
    }

    /* Build argument list */
    if ((arg1 = PyTuple_New(1)) == NULL) {
        Py_DECREF(sequence);
        return NULL;
    }
    PyTuple_SET_ITEM(arg1, 0, (PyObject *)sequence);
    /* arg1 now owns sequence reference */
    if ((arglist = PySequence_Concat(arg1, args)) == NULL) {
        goto fail;
    }

    Py_DECREF(arg1);    /* arglist has a reference to sequence, now. */
    arg1 = NULL;

    /*
     * Call function object --- variable passed to routine.  Extra
     * arguments are in another passed variable.
     */
    if ((result = PyEval_CallObject(func, arglist))==NULL) {
        goto fail;
    }

    result_array = (PyArrayObject *) PyArray_ContiguousFromObject(result,
                                                        NPY_DOUBLE, 0, 0);
    if (result_array == NULL) {
        goto fail;
    }

    Py_DECREF(result);
    Py_DECREF(arglist);
    return (PyObject *) result_array;

fail:
    Py_XDECREF(arglist);
    Py_XDECREF(result);
    Py_XDECREF(arg1);
    return NULL;
}


static PyObject *odepack_error;

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

void
ode_function(int *n, double *t, double *y, double *ydot)
{
    /*
    This is the function called from the Fortran code it should
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
            "The size of the array returned by func (%ld) does not match "
            "the size of y0 (%d).",
            PyArray_Size((PyObject *)result_array), *n);
        *n = -1;
        Py_DECREF(arglist);
        Py_DECREF(result_array);
        return;
    }

    memcpy(ydot, PyArray_DATA(result_array), (*n)*sizeof(double));
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


int
ode_jacobian_function(int *n, double *t, double *y, int *ml, int *mu,
                      double *pd, int *nrowpd)
{
    /*
        This is the function called from the Fortran code it should
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
    if ((arglist = PySequence_Concat(arg1, multipack_extra_arguments)) == NULL) {
        *n = -1;
        Py_DECREF(arg1);
        return -1;
    }
    Py_DECREF(arg1);    /* arglist has reference */

    result_array = (PyArrayObject *)call_python_function(multipack_python_jacobian,
                                                         *n, y, arglist, odepack_error);
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
     *  multipack_jac_transpose is !col_deriv, so if multipack_jac_transpose
     *  is 0, the array created by the user is already in Fortran order, and
     *  a transpose is not needed when it is copied to pd.
     */

    if ((multipack_jac_type == 1) && !multipack_jac_transpose) {
        /* Full Jacobian, no transpose needed, so we can use memcpy. */
        memcpy(pd, PyArray_DATA(result_array), (*n)*(*nrowpd)*sizeof(double));
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
        copy_array_to_fortran(pd, *nrowpd, m, *n,
            (double *) PyArray_DATA(result_array),
            !multipack_jac_transpose);
    }

    Py_DECREF(arglist);
    Py_DECREF(result_array);
    return 0;
}


int
setup_extra_inputs(PyArrayObject **ap_rtol, PyObject *o_rtol,
                   PyArrayObject **ap_atol, PyObject *o_atol,
                   PyArrayObject **ap_tcrit, PyObject *o_tcrit,
                   int *numcrit, int neq)
{
    int itol = 0;
    double tol = 1.49012e-8;
    npy_intp one = 1;

    /* Setup tolerances */
    if (o_rtol == NULL) {
        *ap_rtol = (PyArrayObject *) PyArray_SimpleNew(1, &one, NPY_DOUBLE);
        if (*ap_rtol == NULL) {
            PYERR2(odepack_error,"Error constructing relative tolerance.");
        }
        *(double *) PyArray_DATA(*ap_rtol) = tol;                /* Default */
    }
    else {
        *ap_rtol = (PyArrayObject *) PyArray_ContiguousFromObject(o_rtol,
                                                            NPY_DOUBLE, 0, 1);
        if (*ap_rtol == NULL) {
            PYERR2(odepack_error,"Error converting relative tolerance.");
        }
        /* XXX Fix the following. */
        if (PyArray_NDIM(*ap_rtol) == 0); /* rtol is scalar */
        else if (PyArray_DIMS(*ap_rtol)[0] == neq) {
            itol |= 2;      /* Set rtol array flag */
        }
        else {
            PYERR(odepack_error, "Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
        }
    }

    if (o_atol == NULL) {
        *ap_atol = (PyArrayObject *) PyArray_SimpleNew(1, &one, NPY_DOUBLE);
        if (*ap_atol == NULL) {
            PYERR2(odepack_error,"Error constructing absolute tolerance");
        }
        *(double *)PyArray_DATA(*ap_atol) = tol;
    }
    else {
        *ap_atol = (PyArrayObject *) PyArray_ContiguousFromObject(o_atol, NPY_DOUBLE, 0, 1);
        if (*ap_atol == NULL) {
            PYERR2(odepack_error,"Error converting absolute tolerance.");
        }
        /* XXX Fix the following. */
        if (PyArray_NDIM(*ap_atol) == 0); /* atol is scalar */
        else if (PyArray_DIMS(*ap_atol)[0] == neq) {
            itol |= 1;        /* Set atol array flag */
        }
        else {
            PYERR(odepack_error,"Tolerances must be an array of the same length as the\n     number of equations or a scalar.");
        }
    }
    itol++;             /* increment to get correct value */

    /* Setup t-critical */
    if (o_tcrit != NULL) {
        *ap_tcrit = (PyArrayObject *) PyArray_ContiguousFromObject(o_tcrit, NPY_DOUBLE, 0, 1);
        if (*ap_tcrit == NULL) {
            PYERR2(odepack_error,"Error constructing critical times.");
        }
        *numcrit = PyArray_Size((PyObject *) (*ap_tcrit));
    }
    return itol;

fail:       /* Needed for use of PYERR */
    return -1;
}


int
compute_lrw_liw(int *lrw, int *liw, int neq, int jt, int ml, int mu,
                int mxordn, int mxords)
{
    int lrn, lrs, nyh, lmat;

    if (jt == 1 || jt == 2) {
        lmat = neq*neq + 2;
    }
    else if (jt == 4 || jt == 5) {
        lmat = (2*ml + mu + 1)*neq + 2;
    }
    else {
        PYERR(odepack_error,"Incorrect value for jt");
    }

    if (mxordn < 0) {
        PYERR(odepack_error,"Incorrect value for mxordn");
    }
    if (mxords < 0) {
        PYERR(odepack_error,"Incorrect value for mxords");
    }
    nyh = neq;

    lrn = 20 + nyh*(mxordn+1) + 3*neq;
    lrs = 20 + nyh*(mxords+1) + 3*neq + lmat;

    *lrw = PyArray_MAX(lrn,lrs);
    *liw = 20 + neq;
    return 0;

fail:
    return -1;
}

static char doc_odeint[] =
    "[y,{infodict,}istate] = odeint(fun, y0, t, args=(), Dfun=None, "
    "col_deriv=0, ml=, mu=, full_output=0, rtol=, atol=, tcrit=, h0=0.0, "
    "hmax=0.0, hmin=0.0, ixpr=0.0, mxstep=0.0, mxhnil=0, mxordn=0, "
    "mxords=0)\n  yprime = fun(y,t,...)";

/*
 *  The macros INIT_JAC_FUNC and RESTORE_JAC_FUNC are only
 *  used in odepack_odeint().
 */

#define INIT_JAC_FUNC(fun,Dfun,arg,col_deriv,errobj,jac_type) { \
    store_multipack_globals[0] = multipack_python_function; \
    store_multipack_globals[1] = multipack_extra_arguments; \
    store_multipack_globals[2] = multipack_python_jacobian; \
    store_multipack_globals3 = multipack_jac_transpose; \
    if (arg == NULL) { \
        if ((arg = PyTuple_New(0)) == NULL) goto fail; \
    } \
    else {\
        Py_INCREF(arg);   /* We decrement on exit. */ \
    } \
    if (!PyTuple_Check(arg)) { \
        PYERR(errobj,"Extra Arguments must be in a tuple"); \
    } \
    /* Set up callback functions */ \
    if (!PyCallable_Check(fun) || (Dfun != Py_None && !PyCallable_Check(Dfun))) { \
        PYERR(errobj,"The function and its Jacobian must be callable functions."); \
    } \
    multipack_python_function = fun; \
    multipack_extra_arguments = arg; \
    multipack_python_jacobian = Dfun; \
    multipack_jac_transpose = !(col_deriv); \
    multipack_jac_type = jac_type; \
}

#define RESTORE_JAC_FUNC() \
    multipack_python_function = store_multipack_globals[0]; \
    multipack_extra_arguments = store_multipack_globals[1]; \
    multipack_python_jacobian = store_multipack_globals[2]; \
    multipack_jac_transpose = store_multipack_globals3;


static PyObject *
odepack_odeint(PyObject *dummy, PyObject *args, PyObject *kwdict)
{
    PyObject *fcn, *y0, *p_tout, *o_rtol = NULL, *o_atol = NULL;
    PyArrayObject *ap_y = NULL, *ap_yout = NULL;
    PyArrayObject *ap_rtol = NULL, *ap_atol = NULL;
    PyArrayObject *ap_tout = NULL;
    PyObject *extra_args = NULL;
    PyObject *Dfun = Py_None;
    int neq, itol = 1, itask = 1, istate = 1, iopt = 0, lrw, *iwork, liw, jt = 4;
    double *y, t, *tout, *rtol, *atol, *rwork;
    double h0 = 0.0, hmax = 0.0, hmin = 0.0;
    int ixpr = 0, mxstep = 0, mxhnil = 0, mxordn = 12, mxords = 5, ml = -1, mu = -1;
    PyObject *o_tcrit = NULL;
    PyArrayObject *ap_tcrit = NULL;
    PyArrayObject *ap_hu = NULL, *ap_tcur = NULL, *ap_tolsf = NULL, *ap_tsw = NULL;
    PyArrayObject *ap_nst = NULL, *ap_nfe = NULL, *ap_nje = NULL, *ap_nqu = NULL;
    PyArrayObject *ap_mused = NULL;
    int imxer = 0, lenrw = 0, leniw = 0, col_deriv = 0;
    npy_intp out_sz = 0, dims[2];
    int k, ntimes, crit_ind = 0;
    int allocated = 0, full_output = 0, numcrit = 0;
    double *yout, *yout_ptr, *tout_ptr, *tcrit;
    double *wa;
    static char *kwlist[] = {"fun", "y0", "t", "args", "Dfun", "col_deriv",
                             "ml", "mu", "full_output", "rtol", "atol", "tcrit",
                             "h0", "hmax", "hmin", "ixpr", "mxstep", "mxhnil",
                             "mxordn", "mxords", NULL};

    PyObject *store_multipack_globals[4];
    int store_multipack_globals3;

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
    if (Dfun == Py_None) {
        /* set jt for internally generated */
        jt++;
    }
    if (ml < 0 && mu < 0) {
        /* neither ml nor mu given, mark jt for full jacobian */
        jt -= 3;
    }
    if (ml < 0) {
        /* if one but not both are given */
        ml = 0;
    }
    if (mu < 0) {
        mu = 0;
    }

    INIT_JAC_FUNC(fcn, Dfun, extra_args, col_deriv, odepack_error, jt);

    /* Initial input vector */
    ap_y = (PyArrayObject *) PyArray_ContiguousFromObject(y0, NPY_DOUBLE, 0, 0);
    if (ap_y == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(ap_y) > 1) {
        PyErr_SetString(PyExc_ValueError, "Initial condition y0 must be one-dimensional.");
        goto fail;
    }
    y = (double *) PyArray_DATA(ap_y);
    neq = PyArray_Size((PyObject *) ap_y);
    dims[1] = neq;

    /* Set of output times for integration */
    ap_tout = (PyArrayObject *) PyArray_ContiguousFromObject(p_tout, NPY_DOUBLE, 0, 0);
    if (ap_tout == NULL) {
        goto fail;
    }
    if (PyArray_NDIM(ap_tout) > 1) {
        PyErr_SetString(PyExc_ValueError, "Output times t must be one-dimensional.");
        goto fail;
    }
    tout = (double *) PyArray_DATA(ap_tout);
    ntimes = PyArray_Size((PyObject *)ap_tout);
    dims[0] = ntimes;
    t = tout[0];

    /* Setup array to hold the output evaluations*/
    ap_yout= (PyArrayObject *) PyArray_SimpleNew(2,dims,NPY_DOUBLE);
    if (ap_yout== NULL) {
        goto fail;
    }
    yout = (double *) PyArray_DATA(ap_yout);
    /* Copy initial vector into first row of output */
    memcpy(yout, y, neq*sizeof(double));  /* copy initial value to output */
    yout_ptr = yout + neq;    /* set output pointer to next position */

    itol = setup_extra_inputs(&ap_rtol, o_rtol, &ap_atol, o_atol, &ap_tcrit,
                              o_tcrit, &numcrit, neq);
    if (itol < 0 ) {
        goto fail;  /* Something didn't work */
    }
    rtol = (double *) PyArray_DATA(ap_rtol);
    atol = (double *) PyArray_DATA(ap_atol);
    if (o_tcrit != NULL) {
        tcrit = (double *)(PyArray_DATA(ap_tcrit));
    }

    /* Find size of working arrays*/
    if (compute_lrw_liw(&lrw, &liw, neq, jt, ml, mu, mxordn, mxords) < 0) {
        goto fail;
    }

    if ((wa = (double *)malloc(lrw*sizeof(double) + liw*sizeof(int)))==NULL) {
        PyErr_NoMemory();
        goto fail;
    }
    allocated = 1;
    rwork = wa;
    iwork = (int *)(wa + lrw);

    iwork[0] = ml;
    iwork[1] = mu;      /* ignored if not needed */

    if (h0 != 0.0 || hmax != 0.0 || hmin != 0.0 || ixpr != 0 || mxstep != 0 ||
            mxhnil != 0 || mxordn != 0 || mxords != 0) {
        rwork[4] = h0;
        rwork[5] = hmax;
        rwork[6] = hmin;
        iwork[4] = ixpr;
        iwork[5] = mxstep;
        iwork[6] = mxhnil;
        iwork[7] = mxordn;
        iwork[8] = mxords;
        iopt = 1;
    }
    istate = 1;
    k = 1;

    /* If full output make some useful output arrays */
    if (full_output) {
        out_sz = ntimes-1;
        ap_hu = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tcur = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tolsf = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_tsw = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_DOUBLE);
        ap_nst = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nfe = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nje = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_nqu = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        ap_mused = (PyArrayObject *) PyArray_SimpleNew(1, &out_sz, NPY_INT);
        if (ap_hu == NULL || ap_tcur == NULL || ap_tolsf == NULL ||
                ap_tsw == NULL || ap_nst == NULL || ap_nfe == NULL ||
                ap_nje == NULL || ap_nqu == NULL || ap_mused == NULL) {
            goto fail;
        }
    }

    if (o_tcrit != NULL) {
        /* There are critical points */
        itask = 4;
        rwork[0] = *tcrit;
    }
    while (k < ntimes && istate > 0) {    /* loop over desired times */

        tout_ptr = tout + k;
        /* Use tcrit if relevant */
        if (itask == 4 && *tout_ptr > *(tcrit + crit_ind)) {
            crit_ind++;
            rwork[0] = *(tcrit+crit_ind);
        }
        if (crit_ind >= numcrit) {
            itask = 1;  /* No more critical values */
        }

        LSODA(ode_function, &neq, y, &t, tout_ptr, &itol, rtol, atol, &itask,
              &istate, &iopt, rwork, &lrw, iwork, &liw,
              ode_jacobian_function, &jt);
        if (full_output) {
            *((double *)PyArray_DATA(ap_hu) + (k-1)) = rwork[10];
            *((double *)PyArray_DATA(ap_tcur) + (k-1)) = rwork[12];
            *((double *)PyArray_DATA(ap_tolsf) + (k-1)) = rwork[13];
            *((double *)PyArray_DATA(ap_tsw) + (k-1)) = rwork[14];
            *((int *)PyArray_DATA(ap_nst) + (k-1)) = iwork[10];
            *((int *)PyArray_DATA(ap_nfe) + (k-1)) = iwork[11];
            *((int *)PyArray_DATA(ap_nje) + (k-1)) = iwork[12];
            *((int *)PyArray_DATA(ap_nqu) + (k-1)) = iwork[13];
            if (istate == -5 || istate == -4) {
                imxer = iwork[15];
            }
            else {
                imxer = -1;
            }
            lenrw = iwork[16];
            leniw = iwork[17];
            *((int *)PyArray_DATA(ap_mused) + (k-1)) = iwork[18];
        }
        if (PyErr_Occurred()) {
            goto fail;
        }
        memcpy(yout_ptr, y, neq*sizeof(double));  /* copy integration result to output*/
        yout_ptr += neq;
        k++;
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
        return Py_BuildValue("N{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:i,s:i,s:i,s:N}i",
                    PyArray_Return(ap_yout),
                    "hu", PyArray_Return(ap_hu),
                    "tcur", PyArray_Return(ap_tcur),
                    "tolsf", PyArray_Return(ap_tolsf),
                    "tsw", PyArray_Return(ap_tsw),
                    "nst", PyArray_Return(ap_nst),
                    "nfe", PyArray_Return(ap_nfe),
                    "nje", PyArray_Return(ap_nje),
                    "nqu", PyArray_Return(ap_nqu),
                    "imxer", imxer,
                    "lenrw", lenrw,
                    "leniw", leniw,
                    "mused", PyArray_Return(ap_mused),
                    istate);
    }
    else {
        return Py_BuildValue("Ni", PyArray_Return(ap_yout), istate);
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
    if (allocated) {
        free(wa);
    }
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


static struct PyMethodDef odepack_module_methods[] = {
    {"odeint", (PyCFunction) odepack_odeint, METH_VARARGS|METH_KEYWORDS, doc_odeint},
    {NULL, NULL, 0, NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_odepack",
    NULL,
    -1,
    odepack_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *
PyInit__odepack(void)
{
    PyObject *m, *d, *s;

    m = PyModule_Create(&moduledef);
    import_array();
    d = PyModule_GetDict(m);

    s = PyUString_FromString(" 1.9 ");
    PyDict_SetItemString(d, "__version__", s);
    odepack_error = PyErr_NewException ("odepack.error", NULL, NULL);
    Py_DECREF(s);
    PyDict_SetItemString(d, "error", odepack_error);
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module odepack");
    }
    return m;
}
#else
PyMODINIT_FUNC init_odepack(void)
{
    PyObject *m, *d, *s;
    m = Py_InitModule("_odepack", odepack_module_methods);
    import_array();
    d = PyModule_GetDict(m);

    s = PyUString_FromString(" 1.9 ");
    PyDict_SetItemString(d, "__version__", s);
    odepack_error = PyErr_NewException ("odepack.error", NULL, NULL);
    Py_DECREF(s);
    if (PyErr_Occurred()) {
        Py_FatalError("can't initialize module odepack");
    }
}
#endif    
