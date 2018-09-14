cimport cpython
from . cimport c_zeros


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef scipy_zeros_parameters *myparams = <scipy_zeros_parameters *> params
    cdef tuple args = <tuple> myparams.args if myparams.args is not NULL else ()
    cdef callback_type_tuple f = myparams.function

    return f(x, args)  # callback_type_tuple takes a double and a tuple


# cythonized way to call scalar bisect
cdef double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef double root
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    root = c_zeros.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return root


# cythonized way to call scalar ridder
cdef double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef double root
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    root = c_zeros.ridder(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return root


# cythonized way to call scalar brenth
cdef double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef double root
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    root = c_zeros.brenth(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return root


# cythonized way to call scalar brentq
cdef double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output):
    cdef double root
    cdef scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    root = c_zeros.brentq(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
    if full_output is not NULL:
        full_output.funcalls = myparams.funcalls
        full_output.iterations = myparams.iterations
        full_output.error_num = myparams.error_num
    return root
