cimport cython
from libc cimport stdio, stdlib

@cython.final
cdef class MessageStream:
    cdef stdio.FILE *handle
    cdef bytes _filename
    cdef bint _removed
    cdef size_t _memstream_size
    cdef char *_memstream_ptr
    cpdef close(self)
