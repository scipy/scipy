# cython: language_level=3, boundscheck=False

from libc.string cimport memcpy
import numpy as np
cimport numpy as np
cimport cython


np.import_array()

ctypedef np.float64_t float64_t


cdef extern from "tnc.h":
    ctypedef void tnc_callback(double[], void*) except *
    ctypedef int tnc_function(double[], double*, double[], void*) except 1

    cdef enum tnc_rc:
        TNC_MINRC
        TNC_ENOMEM
        TNC_EINVAL
        TNC_INFEASIBLE
        TNC_LOCALMINIMUM
        TNC_FCONVERGED
        TNC_XCONVERGED
        TNC_MAXFUN
        TNC_LSFAIL
        TNC_CONSTANT
        TNC_NOPROGRESS
        TNC_USERABORT

    cdef int tnc(int n, double x[], double *f, double g[],
            tnc_function *function, void *state,
            double low[], double up[], double scale[], double offset[],
            int messages, int maxCGit, int maxnfeval, double eta, double stepmx,
            double accuracy, double fmin, double ftol, double xtol, double pgtol,
            double rescale, int *nfeval, int *niter, tnc_callback *tnc_callback) except *


cdef struct s_pytnc_state:
  void *py_function
  void *py_callback
  int n
  int failed
ctypedef s_pytnc_state pytnc_state


cdef int function(double x[], double *f, double g[], void *state) except 1:
    # return 0 if no error, 1 if error
    # However, if the user function raises Exceptions Cython should deal
    # with it.

    cdef:
        pytnc_state *py_state
        int n
        double *x_data
        double *g_data

    py_state = <pytnc_state *>state
    n = py_state.n

    if py_state.failed:
        # perhaps the callback code had an exception?
        # this will cause the tnc code to have an LS_USERABORT
        return 1

    # ensures we're working on a copy of the data, just in case user function
    # mutates it
    xcopy = np.empty(n, dtype=np.float64)
    x_data = <float64_t *>np.PyArray_DATA(xcopy)
    memcpy(x_data, x, sizeof(double) * n)

    fx, gx = (<object>py_state.py_function)(xcopy)

    if not np.isscalar(fx):
        try:
            fx = np.asarray(fx).item()
        except (TypeError, ValueError) as e:
            raise ValueError(
                "The user-provided objective function "
                "must return a scalar value."
            ) from e

    f[0] = <double> fx

    gx = np.asfarray(gx)
    if gx.size != n:
        raise ValueError("tnc: gradient must have shape (len(x0),)")

    if not gx.flags['C_CONTIGUOUS']:
        gx = np.ascontiguousarray(gx, dtype=np.float64)

    g_data = <float64_t *>np.PyArray_DATA(gx)
    memcpy(g, g_data, n * sizeof(double))

    return 0


cdef void callback_function(double x[], void *state) except *:
    cdef:
        pytnc_state *py_state
        int n
        double *x_data

    py_state = <pytnc_state *>state
    n = py_state.n

    # ensures we're working on a copy of the data, just in case user function
    # mutates it
    try:
        xcopy = np.empty(n, dtype=np.float64)
        x_data = <float64_t *>np.PyArray_DATA(xcopy)
        memcpy(x_data, x, sizeof(double) * n)

        # TODO examine val to see if we should halt?
        val = (<object>py_state.py_callback)(xcopy)
    except BaseException as exc:
        py_state.failed = 1
        raise exc


def tnc_minimize(func_and_grad,
                   np.ndarray[np.float64_t] x0,
                   np.ndarray[np.float64_t] low,
                   np.ndarray[np.float64_t] up,
                   np.ndarray[np.float64_t] scale,
                   np.ndarray[np.float64_t] offset,
                   int messages,
                   int maxCGit,
                   int maxfun,
                   double eta,
                   double stepmx,
                   double accuracy,
                   double fmin,
                   double ftol,
                   double xtol,
                   double pgtol,
                   double rescale,
                   callback=None):

    cdef:
        pytnc_state py_state
        int n
        int rc
        double f = np.inf
        int nfeval = 0
        int niter = 0
        double *x_data
        double *g_data
        double *low_data = NULL
        double *up_data = NULL
        double *scale_data = NULL
        double *offset_data = NULL
        tnc_callback *callback_functionP = NULL

    # initialise structure, use memset?
    py_state.py_callback = NULL
    py_state.failed = 0

    if not callable(func_and_grad):
        raise TypeError("tnc: function must be callable")

    if callback is not None:
        if not callable(callback):
            raise TypeError("tnc: callback must be callable or None.")
        py_state.py_callback = <void*> callback
        callback_functionP = callback_function

    n = x0.size
    py_state.n = n

    n3 = np.size(scale)
    if n3:
        scale_data = <float64_t *>np.PyArray_DATA(scale)

    n1 = np.size(low)
    if n1:
        low_data = <float64_t *>np.PyArray_DATA(low)

    n2 = np.size(up)
    if n2:
        up_data = <float64_t *>np.PyArray_DATA(up)

    n4 = np.size(offset)
    if n4:
        offset_data = <float64_t *>np.PyArray_DATA(offset)

    if (n1 != n2 or
        n != n1 or
        (scale_data != NULL and n != n3) or
        (offset_data != NULL and n != n4)):
        raise ValueError("tnc: vector sizes must be equal")

    x = np.copy(x0, order="C")
    g = np.zeros_like(x, dtype= np.float64)
    x_data = <float64_t *>np.PyArray_DATA(x)
    g_data = <float64_t *>np.PyArray_DATA(g)

    py_state.py_function = <void*> func_and_grad

    rc = tnc(n, x_data, &f, g_data, function, <void*> &py_state, low_data, up_data,
             scale_data, offset_data, messages, maxCGit, maxfun, eta, stepmx,
             accuracy, fmin, ftol, xtol, pgtol, rescale,
             &nfeval, &niter, callback_functionP);

    if py_state.failed or rc == TNC_ENOMEM:
        return None

    return rc, nfeval, niter, x, f, g
