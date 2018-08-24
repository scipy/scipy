ctypedef double (*callback_type)(double, void*)

cdef double newton(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter)

cdef double secant(callback_type func, double x0, void* args, double tol, int maxiter)

cdef double halley(callback_type func, double x0, callback_type fprime, void* args, double tol, int maxiter, callback_type fprime2)

cdef double bisect(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

cdef double ridder(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

cdef double brenth(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)

cdef double brentq(callback_type f, double xa, double xb, void* args, double xtol, double rtol, int iter)
