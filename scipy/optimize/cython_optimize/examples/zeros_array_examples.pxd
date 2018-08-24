cdef double f_solarcell(int n, double* args)

cdef double fprime(int n, double* args)

cdef double fprime2(int n, double* args)

cdef double solarcell_newton(tuple args)

cdef double solarcell_secant(tuple args)

cdef double solarcell_halley(tuple args)

cdef double solarcell_bisect(tuple args)

cdef double solarcell_ridder(tuple args)

cdef double solarcell_brenth(tuple args)

cdef double solarcell_brentq(tuple args)
