# Defines for mio_utils -*- python -*- like

cdef size_t cproduct(tup)

cdef class FileReadOpts:
    cdef:
        int chars_as_strings
        int mat_dtype
        int squeeze_me

