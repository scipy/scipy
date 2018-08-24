from __future__ import division, print_function, absolute_import
import warnings
cimport cpython
from . cimport c_zeros


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_func(double x, void *params):
    cdef c_zeros.scipy_newton_parameters *myparams = <c_zeros.scipy_newton_parameters *> params
    cdef tuple args = <tuple> myparams.args if myparams.args is not NULL else ()
    cdef callback_type_tuple func = myparams.function

    return func(x, args)  # callback_type_tuple takes a double and a tuple


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime(double x, void *params):
    cdef c_zeros.scipy_newton_parameters *myparams = <c_zeros.scipy_newton_parameters *> params
    cdef tuple args = <tuple> myparams.args if myparams.args is not NULL else ()
    cdef callback_type_tuple fprime = myparams.function_derivative

    return fprime(x, args)  # callback_type_tuple takes a double and a tuple


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime2(double x, void *params):
    cdef c_zeros.scipy_newton_parameters *myparams = <c_zeros.scipy_newton_parameters *> params
    cdef tuple args = <tuple> myparams.args if myparams.args is not NULL else ()
    cdef callback_type_tuple fprime2 = myparams.function_second_derivative

    return fprime2(x, args)  # callback_type_tuple takes a double and a tuple


# Newton-Raphson method
cdef double newton(callback_type_tuple func, double p0, callback_type_tuple fprime, tuple args, double tol, int maxiter):
    cdef c_zeros.scipy_newton_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = func
    myparams.function_derivative = fprime
    return c_zeros.newton(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros.default_parameters *> &myparams, tol, maxiter)


# Secant method
cdef double secant(callback_type_tuple func, double p0, tuple args, double tol, int maxiter):
    cdef c_zeros.scipy_newton_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = func
    return c_zeros.secant(scipy_newton_functions_func, p0, <c_zeros.default_parameters *> &myparams, tol, maxiter)


# Halley's method
cdef double halley(callback_type_tuple func, double p0, callback_type_tuple fprime, tuple args, double tol, int maxiter, callback_type_tuple fprime2):
    cdef c_zeros.scipy_newton_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = func
    myparams.function_derivative = fprime
    myparams.function_second_derivative = fprime2
    return c_zeros.halley(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros.default_parameters *> &myparams, tol, maxiter, scipy_newton_functions_fprime2)


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef c_zeros.scipy_zeros_parameters *myparams = <c_zeros.scipy_zeros_parameters *> params
    cdef tuple args = <tuple> myparams.args if myparams.args is not NULL else ()
    cdef callback_type_tuple f = myparams.function

    return f(x, args)  # callback_type_tuple takes a double and a tuple


# cythonized way to call scalar bisect
cdef double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter):
    cdef c_zeros.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    return c_zeros.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)


# cythonized way to call scalar ridder
cdef double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter):
    cdef c_zeros.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    return c_zeros.ridder(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)


# cythonized way to call scalar brenth
cdef double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter):
    cdef c_zeros.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    return c_zeros.brenth(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)


# cythonized way to call scalar brentq
cdef double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter):
    cdef c_zeros.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = <cpython.PyObject *> args
    myparams.function = f
    return c_zeros.brentq(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros.default_parameters *> &myparams)
