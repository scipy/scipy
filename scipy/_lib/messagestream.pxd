cimport cython
from libc cimport stdio, stdlib

@cython.final
cdef class MessageStream:
    # Borrowed: valid only until `close()`, after which the stream may be
    # reused by another MessageStream. Write to it and read the messages
    # with `get()`; never close or reopen it, or change its mode or buffering.
    cdef stdio.FILE *handle
    cdef bint _is_memstream
    cdef bint _too_large
    cdef size_t _memstream_size
    cdef char *_memstream_ptr
    cpdef close(self)
