ctypedef double (*callback_type_array)(int, double*)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_array function
    int n
    double* args

cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)
