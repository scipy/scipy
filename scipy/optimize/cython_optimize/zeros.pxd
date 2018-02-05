
ctypedef float (*callback_type)(float, tuple);

cdef float newton(callback_type func, float x0, callback_type fprime, tuple args)
