ctypedef double (*callback_type_array)(int, double*)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_array function
    int n
    double* args

ctypedef struct scipy_newton_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_array function
    callback_type_array function_derivative
    callback_type_array function_second_derivative
    int n
    double* args

cdef double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter)

cdef double secant(callback_type_array func, double x0, int n, double* args, double tol, int maxiter)

cdef double halley(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2)

cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)
