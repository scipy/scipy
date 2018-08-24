from __future__ import division, print_function, absolute_import
import warnings
from . cimport c_zeros_struct


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_func(double x, void *params):
    cdef c_zeros_struct.scipy_newton_parameters *myparams = <c_zeros_struct.scipy_newton_parameters *> params
    cdef void* args = myparams.args
    cdef callback_type func = myparams.function

    return func(x, args)  # callback_type takes a double and a struct


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime(double x, void *params):
    cdef c_zeros_struct.scipy_newton_parameters *myparams = <c_zeros_struct.scipy_newton_parameters *> params
    cdef void* args = myparams.args
    cdef callback_type fprime = myparams.function_derivative

    return fprime(x, args)  # callback_type takes a double and a struct


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime2(double x, void *params):
    cdef c_zeros_struct.scipy_newton_parameters *myparams = <c_zeros_struct.scipy_newton_parameters *> params
    cdef void* args = myparams.args
    cdef callback_type fprime2 = myparams.function_second_derivative

    return fprime2(x, args)  # callback_type takes a double and a struct


# Newton-Raphson method
cdef double newton(callback_type func, double p0, callback_type fprime, void *args, double tol, int maxiter):
    cdef c_zeros_struct.scipy_newton_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = func
    myparams.function_derivative = fprime
    return c_zeros_struct.newton(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros_struct.default_parameters *> &myparams, tol, maxiter)


# Secant method
cdef double secant(callback_type func, double p0, void *args, double tol, int maxiter):
    cdef c_zeros_struct.scipy_newton_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = func
    return c_zeros_struct.secant(scipy_newton_functions_func, p0, <c_zeros_struct.default_parameters *> &myparams, tol, maxiter)


# Halley's method
cdef double halley(callback_type func, double p0, callback_type fprime, void *args, double tol, int maxiter, callback_type fprime2):
    cdef c_zeros_struct.scipy_newton_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = func
    myparams.function_derivative = fprime
    myparams.function_second_derivative = fprime2
    return c_zeros_struct.halley(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros_struct.default_parameters *> &myparams, tol, maxiter, scipy_newton_functions_fprime2)


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef c_zeros_struct.scipy_zeros_parameters *myparams = <c_zeros_struct.scipy_zeros_parameters *> params
    cdef void* args = myparams.args
    cdef callback_type f = myparams.function

    return f(x, args)  # callback_type takes a double and a struct


# cythonized way to call scalar bisect
cdef double bisect(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter):
    cdef c_zeros_struct.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    return c_zeros_struct.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_struct.default_parameters *> &myparams)


# cythonized way to call scalar ridder
cdef double ridder(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter):
    cdef c_zeros_struct.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    return c_zeros_struct.ridder(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_struct.default_parameters *> &myparams)


# cythonized way to call scalar brenth
cdef double brenth(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter):
    cdef c_zeros_struct.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    return c_zeros_struct.brenth(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_struct.default_parameters *> &myparams)


# cythonized way to call scalar brentq
cdef double brentq(callback_type f, double xa, double xb, void *args, double xtol, double rtol, int iter):
    cdef c_zeros_struct.scipy_zeros_parameters myparams
    # create params struct
    myparams.args = args
    myparams.function = f
    return c_zeros_struct.brentq(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_struct.default_parameters *> &myparams)
