from .base cimport OdeSolver

cdef class RungeKutta(OdeSolver):
    cdef public:
        object C
        object A
        object B
        object E
        object P
        object order
        object n_stages

cdef class RK23(RungeKutta):
    pass

cdef class RK45(RungeKutta):
    pass
