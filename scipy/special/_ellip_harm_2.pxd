from . cimport sf_error

cdef double _F_integrand(double t, void *data) nogil
cdef double _F_integrand1(double t, void *data) nogil
cdef double _F_integrand2(double t, void *data) nogil
cdef double _F_integrand3(double t, void *data) nogil
cdef double _F_integrand4(double t, void *data) nogil
cdef void _set_action(sf_error.sf_error_t code, sf_error.sf_action_t action) noexcept nogil
