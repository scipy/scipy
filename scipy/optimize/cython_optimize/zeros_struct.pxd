ctypedef double (*callback_type)(double, void*)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type function
    void* args

ctypedef struct scipy_newton_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type function
    callback_type function_derivative
    callback_type function_second_derivative
    void* args

cdef double newton(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

cdef double secant(callback_type func, double x0, void* args, double tol, int maxiter, scipy_newton_parameters *full_output)

cdef double halley(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, callback_type fprime2, scipy_newton_parameters *full_output)

cdef double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)
