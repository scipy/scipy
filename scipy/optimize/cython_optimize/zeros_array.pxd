ctypedef double (*callback_type_array)(int, double*)

cdef double newton(callback_type_array func, double x0, callback_type_array fprime, int n, double* args)

cdef double bisect(callback_type_array f, double xa, double xb, int n, double* args, double xtol, double rtol, int iter)
