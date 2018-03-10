ctypedef double (*callback_type)(double, void*);

cdef double newton(callback_type func, double x0, callback_type_tup fprime, tuple args)

cdef double bisect(callback_type f, double xa, double xb, tuple args, double xtol, double rtol, int iter)
