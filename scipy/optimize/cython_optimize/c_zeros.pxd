ctypedef double (*callback_type_tup)(double, tuple)

cdef extern from "../zeros.c":
    # a parameter structure that has the function and the args tuple
    # this is a contract with the zero solver functions that can't be broken
    ctypedef struct scipy_zeros_parameters:
        int funcalls
        int iterations
        int error_num
        callback_type_tup function
        void *args

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
