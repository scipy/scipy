# -*- python -*- or rather like

cdef class Memholder:
    cdef void *ptr

cdef class GenericStream:
    cdef object fobj

    cpdef int seek(self, long int offset, int whence=*) except -1
    cpdef long int tell(self) except -1
    cdef int read_into(self, void *buf, size_t n) except -1
    cdef Memholder read_alloc(self, size_t n)

cpdef GenericStream make_stream(object fobj)
