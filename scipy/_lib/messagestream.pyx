cimport cython
from libc cimport stdio, stdlib
from cpython cimport PyBytes_FromStringAndSize
from cpython.pythread cimport (
    PyThread_type_lock, PyThread_allocate_lock, PyThread_acquire_lock,
    PyThread_release_lock, WAIT_LOCK)

import os
import tempfile

cdef extern from "messagestream.h":
    stdio.FILE *messagestream_open_memstream(char **, size_t *)


# Without ``open_memstream`` (e.g., on Windows), every stream is a temporary
# file, and creating one is slow: it costs more than a small Qhull call and
# serializes concurrent processes (gh-26693). So instead of closing them,
# `close` rewinds such streams and keeps up to `_MAX_SPARE` of them here for
# reuse by later streams. The limit is kept small because the Windows C
# runtime allows only ~512 open ``FILE*`` streams per process (gh-21958).
cdef enum:
    _MAX_SPARE = 16

cdef stdio.FILE *_spare[_MAX_SPARE]
cdef int _n_spare = 0
cdef PyThread_type_lock _spare_lock = PyThread_allocate_lock()
if _spare_lock == NULL:
    raise MemoryError("thread lock allocation failed")


cdef stdio.FILE *_take_spare() noexcept nogil:
    global _n_spare
    cdef stdio.FILE *handle = NULL
    PyThread_acquire_lock(_spare_lock, WAIT_LOCK)
    if _n_spare > 0:
        _n_spare -= 1
        handle = _spare[_n_spare]
    PyThread_release_lock(_spare_lock)
    return handle


cdef bint _put_spare(stdio.FILE *handle) noexcept nogil:
    global _n_spare
    cdef bint kept = False
    PyThread_acquire_lock(_spare_lock, WAIT_LOCK)
    if _n_spare < _MAX_SPARE:
        _spare[_n_spare] = handle
        _n_spare += 1
        kept = True
    PyThread_release_lock(_spare_lock)
    return kept


@cython.final
cdef class MessageStream:
    """
    Capture messages emitted to FILE* streams. Do this by directing them
    to a temporary file, residing in memory (if possible) or on disk.
    """

    def __cinit__(self):
        # Try first in-memory files, if available
        self._memstream_ptr = NULL
        self.handle = messagestream_open_memstream(&self._memstream_ptr,
                                                   &self._memstream_size)
        if self.handle != NULL:
            self._is_memstream = 1
            return

        # Fall back to temporary files, reusing one if possible
        self.handle = _take_spare()
        if self.handle != NULL:
            return

        # The OS deletes the file when the last handle to it is closed, also
        # if the process is killed, so there is nothing to clean up.
        with tempfile.TemporaryFile(prefix='scipy-') as f:
            fd = os.dup(f.fileno())

        self.handle = stdio.fdopen(fd, 'wb+')
        if self.handle == NULL:
            os.close(fd)
            raise OSError("Failed to open a temporary file")

    def __dealloc__(self):
        self.close()

    def get(self):
        cdef long pos
        cdef size_t nread
        cdef char *buf = NULL
        cdef bytes obj

        pos = stdio.ftell(self.handle)
        if pos <= 0:
            return ""

        if self._is_memstream:
            stdio.fflush(self.handle)
            obj = PyBytes_FromStringAndSize(self._memstream_ptr, pos)
        else:
            buf = <char*>stdlib.malloc(pos)
            if buf == NULL:
                raise MemoryError()

            try:
                stdio.rewind(self.handle)
                nread = stdio.fread(buf, 1, pos, self.handle)
                if nread != <size_t>pos:
                    raise OSError("failed to read messages from buffer")

                obj = PyBytes_FromStringAndSize(buf, nread)
            finally:
                stdlib.free(buf)

        return obj.decode('latin1')

    def clear(self):
        stdio.rewind(self.handle)

    cpdef close(self):
        if self.handle == NULL:
            return

        if self._is_memstream:
            stdio.fclose(self.handle)
            stdlib.free(self._memstream_ptr)
            self._memstream_ptr = NULL
        else:
            stdio.rewind(self.handle)
            if not _put_spare(self.handle):
                stdio.fclose(self.handle)

        self.handle = NULL
