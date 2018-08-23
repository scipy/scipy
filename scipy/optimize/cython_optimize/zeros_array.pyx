from __future__ import division, print_function, absolute_import
import warnings
from . cimport c_zeros_array

DEF MAXARGS = 10


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_func(double x, void *params):
    cdef c_zeros_array.scipy_newton_parameters *myparams = <c_zeros_array.scipy_newton_parameters *> params
    cdef int n = myparams.n, i
    cdef double[MAXARGS] myargs
    cdef callback_type_array func = myparams.function

    myargs[0] = x
    for i in range(n):
        myargs[i+1] = myparams.args[i]

    return func(n, myargs)  # callback_type_array takes a double and a tuple


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime(double x, void *params):
    cdef c_zeros_array.scipy_newton_parameters *myparams = <c_zeros_array.scipy_newton_parameters *> params
    cdef int n = myparams.n, i
    cdef double[MAXARGS] myargs
    cdef callback_type_array fprime = myparams.function_derivative

    myargs[0] = x
    for i in range(n):
        myargs[i+1] = myparams.args[i]

    return fprime(n, myargs)  # callback_type_array takes a double and a tuple


# callback function wrapper that extracts function, args from params struct
cdef double scipy_newton_functions_fprime2(double x, void *params):
    cdef c_zeros_array.scipy_newton_parameters *myparams = <c_zeros_array.scipy_newton_parameters *> params
    cdef int n = myparams.n, i
    cdef double[MAXARGS] myargs
    cdef callback_type_array fprime2 = myparams.function_second_derivative

    myargs[0] = x
    for i in range(n):
        myargs[i+1] = myparams.args[i]

    return fprime2(n, myargs)  # callback_type_array takes a double and a tuple


# Newton-Raphson method
cdef double newton(callback_type_array func, double p0, callback_type_array fprime, int n, double* args, double tol, int maxiter):
    cdef c_zeros_array.scipy_newton_parameters myparams
    # create params struct
    myparams.n = n
    myparams.args = args
    myparams.function = func
    myparams.function_derivative = fprime
    return c_zeros_array.newton(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros_array.default_parameters *> &myparams, tol, maxiter)


# Secant method
cdef double secant(callback_type_array func, double p0, int n, double* args, double tol, int maxiter):
    cdef c_zeros_array.scipy_newton_parameters myparams
    # create params struct
    myparams.n = n
    myparams.args = args
    myparams.function = func
    return c_zeros_array.secant(scipy_newton_functions_func, p0, <c_zeros_array.default_parameters *> &myparams, tol, maxiter)


# Halley's method
cdef double halley(callback_type_array func, double p0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2):
    cdef c_zeros_array.scipy_newton_parameters myparams
    # create params struct
    myparams.n = n
    myparams.args = args
    myparams.function = func
    myparams.function_derivative = fprime
    myparams.function_second_derivative = fprime2
    return c_zeros_array.halley(scipy_newton_functions_func, p0, scipy_newton_functions_fprime, <c_zeros_array.default_parameters *> &myparams, tol, maxiter, scipy_newton_functions_fprime2)


# callback function wrapper that extracts function, args from params struct
cdef double scipy_zeros_functions_func(double x, void *params):
    cdef c_zeros_array.scipy_zeros_parameters *myparams = <c_zeros_array.scipy_zeros_parameters *> params
    cdef int n = myparams.n, i
    cdef double[MAXARGS] myargs
    cdef callback_type_array f = myparams.function

    myargs[0] = x
    for i in range(n):
        myargs[i+1] = myparams.args[i]

    return f(n, myargs)


# cythonized way to call scalar bisect
cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter):
    cdef c_zeros_array.scipy_zeros_parameters myparams
    # create params struct
    myparams.n = n
    myparams.args = args
    myparams.function = f
    return c_zeros_array.bisect(scipy_zeros_functions_func, xa, xb, xtol, rtol, iter, <c_zeros_array.default_parameters *> &myparams)
