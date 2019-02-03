cdef extern from "../Zeros/zeros.h":
    ctypedef double (*callback_type)(double, void*)
    ctypedef struct default_parameters:
        int funcalls
        int iterations
        int error_num

cdef extern from "../Zeros/bisect.c" nogil:
    double bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)

cdef extern from "../Zeros/ridder.c" nogil:
    double ridder(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)

cdef extern from "../Zeros/brenth.c" nogil:
    double brenth(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)

cdef extern from "../Zeros/brentq.c" nogil:
    double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)
