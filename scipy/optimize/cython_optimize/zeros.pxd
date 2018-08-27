cimport cpython

ctypedef double (*callback_type_tuple)(double, tuple)

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

cdef double newton(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter)

cdef double secant(callback_type_tuple func, double x0, tuple args, double tol, int maxiter)

cdef double halley(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args, double tol, int maxiter, callback_type_tuple fprime2)

cdef double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

cdef double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

cdef double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)

cdef double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)
