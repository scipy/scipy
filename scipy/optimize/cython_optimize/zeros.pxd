
ctypedef double (*callback_type)(double, tuple);

cdef double newton(callback_type func, double x0, callback_type fprime, tuple args)
