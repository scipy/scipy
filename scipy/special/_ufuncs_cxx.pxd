from . cimport sf_error
cdef void _set_action(sf_error.sf_error_t, sf_error.sf_action_t) nogil
cdef void *_export_faddeeva_dawsn
cdef void *_export_faddeeva_dawsn_complex
cdef void *_export_faddeeva_erf
cdef void *_export_faddeeva_erfc
cdef void *_export_faddeeva_erfcx
cdef void *_export_faddeeva_erfcx_complex
cdef void *_export_faddeeva_erfi
cdef void *_export_faddeeva_erfi_complex
cdef void *_export_faddeeva_log_ndtr
cdef void *_export_faddeeva_ndtr
cdef void *_export_faddeeva_w
cdef void *_export_wrightomega