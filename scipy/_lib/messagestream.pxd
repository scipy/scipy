cimport cython
from libc cimport stdio, stdlib

@cython.final
cdef class MessageStream:
    # The stream that messages are written to. Owned by this MessageStream
    # until `close()`, which sets it to NULL; the stream may then be reused
    # by another MessageStream. Code given a copy of it only borrows it: it
    # must not use the copy after `close()`, nor close or reopen the stream
    # or change its mode or buffering.
    cdef stdio.FILE *handle
    cdef bint _is_memstream
    cdef bint _too_large
    cdef size_t _memstream_size
    cdef char *_memstream_ptr
    cpdef close(self)
