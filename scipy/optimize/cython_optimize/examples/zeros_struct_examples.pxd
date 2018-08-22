cdef double f_solarcell(double i, void *args)

cdef double fprime(double i, void *args)

cdef double fprime2(double i, void *args)

cdef double solarcell_newton(tuple args)

cdef double solarcell_secant(tuple args)

cdef double solarcell_halley(tuple args)

cdef double solarcell_bisect(tuple args)
