cdef double newton(callback_type func, double x0, callback_type fprime, void* args)

cdef double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)
