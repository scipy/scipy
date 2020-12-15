# destined to be used in a LowLevelCallable
cdef double _geninvgauss_pdf(double x, void *user_data) nogil except *
cdef double _genstudentized_range_cdf(int n, double[2] x, void *user_data) nogil
cdef double _genstudentized_range_cdf_asymptomatic(double z, void *user_data) nogil
cdef double _genstudentized_range_pdf(int n, double[2] x, void *user_data) nogil