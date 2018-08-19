ctypedef double (*callback_type_tuple)(double, tuple);

cdef double newton(callback_type_tuple func, double x0, callback_type_tuple fprime, tuple args)

cdef double bisect(callback_type_tuple f, double xa, double xb, tuple args, double xtol, double rtol, int iter)
