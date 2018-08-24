ctypedef double (*callback_type_array)(int, double*)

cdef double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter)

cdef double secant(callback_type_array func, double x0, int n, double* args, double tol, int maxiter)

cdef double halley(callback_type_array func, double x0, callback_type_array fprime, int n, double* args, double tol, int maxiter, callback_type_array fprime2)

cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double ridder(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double brenth(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)

cdef double brentq(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)
