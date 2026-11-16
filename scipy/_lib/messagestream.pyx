cimport cython
from libc cimport stdio, stdlib
from cpython cimport PyBytes_FromStringAndSize

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
# Reuse only rewinds a file, so it keeps its largest size until it is closed;
# files larger than `_MAX_SPARE_SIZE` bytes (e.g., from Qhull's trace
# options) are therefore closed instead of kept.
cdef enum:
    _MAX_SPARE = 16
    _MAX_SPARE_SIZE = 65536

cdef stdio.FILE *_spare[_MAX_SPARE]
cdef int _n_spare = 0
# Guards `_spare` and `_n_spare`. Hold it only around plain C statements (no
# calls into Python, no other locks), so that it cannot deadlock.
cdef cython.pymutex _spare_lock


cdef stdio.FILE *_take_spare() noexcept nogil:
    """Return a spare stream, now owned by the caller, or NULL if none."""
    global _n_spare
    cdef stdio.FILE *handle = NULL
    with _spare_lock:
        if _n_spare > 0:
            _n_spare -= 1
            handle = _spare[_n_spare]
    return handle


cdef inline bint _past_max_size(stdio.FILE *handle) noexcept nogil:
    # Checks the position, which is cheap, rather than the file size, which
    # costs a slow filesystem query on Windows. Callers only write at the
    # position, and spare files are at most `_MAX_SPARE_SIZE` bytes, so a
    # file grows larger than that only if the position does too.
    cdef long pos = stdio.ftell(handle)
    return pos < 0 or pos > _MAX_SPARE_SIZE


cdef void _recycle(stdio.FILE *handle, bint too_large) noexcept nogil:
    """Take ownership of a file stream: keep it for reuse, or close it."""
    global _n_spare
    if too_large or _past_max_size(handle):
        stdio.fclose(handle)
        return
    stdio.rewind(handle)
    with _spare_lock:
        if _n_spare < _MAX_SPARE:
            _spare[_n_spare] = handle
            _n_spare += 1
            handle = NULL
    if handle != NULL:
        stdio.fclose(handle)


# Test hooks below, not public API: they let the code above also be tested
# where ``open_memstream`` exists.
cdef bint _use_memstream = True


def _force_tempfile(bint enabled):
    """For tests only: make new streams use temporary files."""
    global _use_memstream
    _use_memstream = not enabled


def _spare_count():
    """For tests only: return the number of streams kept for reuse."""
    cdef int n
    with _spare_lock:
        n = _n_spare
    return n


def _close_spares():
    """For tests only: close all streams kept for reuse."""
    cdef stdio.FILE *handle = _take_spare()
    while handle != NULL:
        stdio.fclose(handle)
        handle = _take_spare()


@cython.final
cdef class MessageStream:
    """
    Capture messages emitted to FILE* streams. Do this by directing them
    to a temporary file, residing in memory (if possible) or on disk.
    """

    def __cinit__(self):
        # Try first in-memory files, if available
        self._memstream_ptr = NULL
        if _use_memstream:
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
        with tempfile.TemporaryFile(prefix='scipy-errmsg-') as f:
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
        if _past_max_size(self.handle):
            self._too_large = 1
        stdio.rewind(self.handle)

    cpdef close(self):
        if self.handle == NULL:
            return

        if self._is_memstream:
            stdio.fclose(self.handle)
            stdlib.free(self._memstream_ptr)
            self._memstream_ptr = NULL
        else:
            _recycle(self.handle, self._too_large)

        self.handle = NULL
