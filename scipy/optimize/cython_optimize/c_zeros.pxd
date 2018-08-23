cimport cpython

ctypedef double (*callback_type_tuple)(double, tuple)

cdef extern from "../Zeros/zeros.h":
    ctypedef double (*callback_type)(double, void*)
    ctypedef struct default_parameters:
        int funcalls
        int iterations
        int error_num

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_tuple function
    cpython.PyObject *args

ctypedef struct scipy_newton_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_tuple function
    callback_type_tuple function_derivative
    callback_type_tuple function_second_derivative
    cpython.PyObject *args


cdef extern from "Newton/newton.c":
    double newton(callback_type func, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter)


cdef extern from "Newton/secant.c":
    double secant(callback_type func, double p0, default_parameters *params, double tol, int maxiter)


cdef extern from "Newton/halley.c":
    double halley(callback_type func, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter, callback_type fprime2)


cdef extern from "../Zeros/bisect.c":
    double bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/ridder.c":
    double ridder(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brenth.c":
    double brenth(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brentq.c":
    double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)
