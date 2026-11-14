cimport cython
from libc cimport stdio, stdlib

@cython.final
cdef class MessageStream:
    cdef stdio.FILE *handle
    cdef bint _is_memstream
    cdef size_t _memstream_size
    cdef char *_memstream_ptr
    cpdef close(self)
