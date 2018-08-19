ctypedef double (*callback_type_array)(int, double)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_array function
    int n
    double[] args

ctypedef double (*callback_type)(double, scipy_zeros_parameters*)

cdef extern from "../Zeros/zeros.h":
    ctypedef struct default_parameters:
        int funcalls
        int iterations
        int error_num

cdef extern from "../Zeros/bisect.c":
    double bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/ridder.c":
    double ridder(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brenth.c":
    double brenth(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)


cdef extern from "../Zeros/brentq.c":
    double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params)
