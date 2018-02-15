
ctypedef double (*callback_type_tup)(double, tuple);

cdef double newton(callback_type_tup func, double x0, callback_type_tup fprime, tuple args)

cdef double bisect(callback_type_tup f, double xa, double xb, tuple args, double xtol, double rtol, int iter)
