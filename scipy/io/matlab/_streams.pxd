# -*- python -*- or rather like

cdef class GenericStream:
    cdef object fobj

    cpdef int seek(self, long int offset, int whence=*) except -1
    cpdef long int tell(self) except -1
    cdef int read_into(self, void *buf, size_t n) except -1
    cdef object read_string(self, size_t n, void **pp, int copy=*)
    cpdef int all_data_read(self) except *

cpdef GenericStream make_stream(object fobj)
