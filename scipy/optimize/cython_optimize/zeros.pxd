ctypedef double (*callback_type)(double, void*)

ctypedef struct zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type function
    void* args
    double root

cdef double bisect(callback_type f, double xa, double xb, void* args,
                   double xtol, double rtol, int iter,
                   zeros_parameters *full_output) nogil

cdef double ridder(callback_type f, double xa, double xb, void* args,
                   double xtol, double rtol, int iter,
                   zeros_parameters *full_output) nogil

cdef double brenth(callback_type f, double xa, double xb, void* args,
                   double xtol, double rtol, int iter,
                   zeros_parameters *full_output) nogil

cdef double brentq(callback_type f, double xa, double xb, void* args,
                   double xtol, double rtol, int iter,
                   zeros_parameters *full_output) nogil
