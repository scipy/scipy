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
# Statistics, see `_fhcache_info`
cdef Py_ssize_t _hits = 0
cdef Py_ssize_t _misses = 0
cdef Py_ssize_t _n_messagestream = 0
cdef Py_ssize_t _n_tempfile = 0
cdef Py_ssize_t _close_errors = 0
# Guards `_spare`, `_n_spare`, `_hits`, `_misses`, `_n_messagestream`,
# `_n_tempfile` and `_close_errors`: take it for every read or write of
# them. Hold it only around plain C statements (no calls into Python, no
# other locks), so that it cannot deadlock.
cdef cython.pymutex _spare_lock


cdef stdio.FILE *_take_spare() noexcept nogil:
    """Return a spare stream, now owned by the caller, or NULL if none."""
    global _n_spare
    cdef stdio.FILE *handle = NULL
    with _spare_lock:
        if _n_spare > 0:
            _n_spare -= 1
            handle = _spare[_n_spare]
            _spare[_n_spare] = NULL
    return handle


cdef void _close_tempfile(stdio.FILE *handle) noexcept nogil:
    """Close a temporary file stream."""
    global _n_tempfile, _close_errors
    # Even if fclose fails (e.g., when flushing fails), it releases the
    # stream and its file, so there is nothing left to clean up. The failure
    # is only counted, for `_fhcache_info`.
    cdef bint failed = stdio.fclose(handle) != 0
    with _spare_lock:
        _n_tempfile -= 1
        if failed:
            _close_errors += 1


cdef enum:  # how a MessageStream got its stream, for `_count_open`
    _GOT_MEMSTREAM, _GOT_SPARE, _GOT_NEW_TEMPFILE


cdef void _count_open(int how) noexcept nogil:
    global _n_messagestream, _hits, _misses, _n_tempfile
    with _spare_lock:
        _n_messagestream += 1
        if how == _GOT_SPARE:
            _hits += 1
        elif how == _GOT_NEW_TEMPFILE:
            _misses += 1
            _n_tempfile += 1


cdef void _count_close() noexcept nogil:
    global _n_messagestream
    with _spare_lock:
        _n_messagestream -= 1


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
        _close_tempfile(handle)
        return
    stdio.rewind(handle)
    with _spare_lock:
        if _n_spare < _MAX_SPARE:
            _spare[_n_spare] = handle
            _n_spare += 1
            handle = NULL
    if handle != NULL:
        _close_tempfile(handle)


def _fhcache_info():
    """Return statistics about the temporary files of MessageStream.

    The returned dict has these keys:

    - ``hits``: streams taken from the spare files kept for reuse
    - ``misses``: temporary files created because no spare was available
    - ``max_spare``: the maximum number of spare files kept
    - ``n_spare``: the number of spare files currently kept
    - ``n_messagestream``: MessageStream objects that currently have a stream,
      in memory or a temporary file
    - ``n_tempfile``: temporary files currently open, in use or spare
    - ``close_errors``: failed attempts to close a temporary file

    ``n_tempfile - n_spare`` can never exceed ``n_messagestream``; if it
    does, a temporary file was lost without being closed. If
    ``n_messagestream`` keeps growing instead, MessageStream objects are kept
    alive elsewhere, e.g., by a reference cycle (gh-24493).
    """
    cdef Py_ssize_t hits, misses, n_messagestream, n_tempfile, close_errors
    cdef int n_spare
    with _spare_lock:
        hits = _hits
        misses = _misses
        n_spare = _n_spare
        n_messagestream = _n_messagestream
        n_tempfile = _n_tempfile
        close_errors = _close_errors
    return {"hits": hits, "misses": misses, "max_spare": _MAX_SPARE,
            "n_spare": n_spare, "n_messagestream": n_messagestream,
            "n_tempfile": n_tempfile, "close_errors": close_errors}


def _fhcache_size():
    """Return the number of spare temporary files kept for reuse."""
    cdef int n
    with _spare_lock:
        n = _n_spare
    return n


def _fhcache_clear():
    """Close the spare temporary files kept for reuse.

    For example before forking, or to release resources. Streams in use by
    a MessageStream are not affected.
    """
    cdef stdio.FILE *handle = _take_spare()
    while handle != NULL:
        _close_tempfile(handle)
        handle = _take_spare()


# Whether new MessageStreams try ``open_memstream`` before falling back to
# temporary files. `_force_tempfile` turns this off, so that the temporary
# files and their reuse can also be tested where ``open_memstream`` exists.
cdef bint _try_memstream = True


def _force_tempfile(bint enabled):
    """Make new MessageStreams use temporary files even where
    ``open_memstream`` exists (used by the tests).
    """
    global _try_memstream
    _try_memstream = not enabled


@cython.final
cdef class MessageStream:
    """
    Capture messages emitted to FILE* streams. Do this by directing them
    to a temporary file, residing in memory (if possible) or on disk.
    """

    def __cinit__(self):
        # Try first in-memory files, if available
        self._memstream_ptr = NULL
        if _try_memstream:
            self.handle = messagestream_open_memstream(&self._memstream_ptr,
                                                       &self._memstream_size)
        if self.handle != NULL:
            self._is_memstream = 1
            _count_open(_GOT_MEMSTREAM)
            return

        # Fall back to temporary files, reusing one if possible
        self.handle = _take_spare()
        if self.handle != NULL:
            _count_open(_GOT_SPARE)
            return

        # The OS deletes the file when the last handle to it is closed, also
        # if the process is killed, so there is nothing to clean up.
        with tempfile.TemporaryFile(prefix='scipy-errmsg-') as f:
            fd = os.dup(f.fileno())

        self.handle = stdio.fdopen(fd, 'wb+')
        if self.handle == NULL:
            os.close(fd)
            raise OSError("Failed to open a temporary file")
        _count_open(_GOT_NEW_TEMPFILE)

    def __dealloc__(self):
        self.close()

    def get(self):
        cdef long pos
        cdef size_t nread
        cdef char *buf = NULL
        cdef bytes obj

        if self.handle == NULL:
            raise ValueError("I/O operation on closed MessageStream")
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
        if self.handle == NULL:
            raise ValueError("I/O operation on closed MessageStream")
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
        _count_close()
