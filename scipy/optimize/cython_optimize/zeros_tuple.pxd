cimport cpython

ctypedef double (*callback_type_tuple)(double, tuple)

ctypedef struct scipy_zeros_parameters:
    int funcalls
    int iterations
    int error_num
    callback_type_tuple function
    cpython.PyObject *args

cdef double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double ridder(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brenth(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)

cdef double brentq(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter, scipy_zeros_parameters *full_output)
