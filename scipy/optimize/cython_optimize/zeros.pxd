ctypedef double (*callback_type)(double, void*)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type function
    void* args
    double root

cdef double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)
