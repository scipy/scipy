# Defines for mio_utils -*- python -*- like

cpdef size_t cproduct(tup)

cdef class FileReadOpts:
    cdef:
        int chars_as_strings
        int mat_dtype
        int squeeze_me

