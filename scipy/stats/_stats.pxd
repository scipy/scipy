# destined to be used in a LowLevelCallable
cdef double _geninvgauss_pdf(double x, void *user_data) nogil except *
cdef double _genhyperbolic_pdf(double x, void *user_data) nogil except *
cdef double _genhyperbolic_logpdf(double x, void *user_data) nogil except *
