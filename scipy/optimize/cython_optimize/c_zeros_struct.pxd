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
    callback_type function
    void* args

cdef extern from "../Zeros/bisect.c":
    double bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/ridder.c":
    double ridder(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brenth.c":
    double brenth(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brentq.c":
    double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)
