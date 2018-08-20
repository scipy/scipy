ctypedef struct test_params:
    double voltage
    double light_current
    double dark_current
    double series_resistance
    double shunt_resistance
    double thermal_voltage

cdef double f_solarcell(double i, void *args)

cdef double fprime(double i, void *args)

cdef double solarcell_newton(tuple args)

cdef double solarcell_bisect(tuple args)
