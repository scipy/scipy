# -*- python -*- or rather like

cdef class GenericStream:
    cdef object fobj

    cpdef int seek(self, long int offset, int whence=*) except -1
    cpdef long int tell(self) except -1
    cpdef int readinto(self, char[:] buf) except -1
    cdef object read_string(self, size_t n, void **pp, int copy=*)


cpdef GenericStream make_stream(object fobj)
