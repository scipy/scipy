cimport cython
from libc cimport stdio, stdlib
from cpython cimport PyBytes_FromStringAndSize

import os
import sys
import tempfile

cdef extern from "messagestream.h":
    stdio.FILE *messagestream_open_memstream(char **, size_t *)


@cython.final
cdef class MessageStream:
    """
    Capture messages emitted to FILE* streams. Do this by directing them
    to a temporary file, residing in memory (if possible) or on disk.
    """

    def __init__(self):
        # Try first in-memory files, if available
        self._memstream_ptr = NULL
        self.handle = messagestream_open_memstream(&self._memstream_ptr,
                                                   &self._memstream_size)
        if self.handle != NULL:
            self._removed = 1
            return

        # Fall back to temporary files
        fd, filename = tempfile.mkstemp(prefix='scipy-')
        os.close(fd)
        self._filename = filename.encode(sys.getfilesystemencoding())
        self.handle = stdio.fopen(self._filename, "wb+")
        if self.handle == NULL:
            stdio.remove(self._filename)
            raise IOError("Failed to open file {0}".format(self._filename))
        self._removed = 0

        # Use a posix-style deleted file, if possible
        if stdio.remove(self._filename) == 0:
            self._removed = 1

    def __del__(self):
        self.close()

    def get(self):
        cdef long pos
        cdef size_t nread
        cdef char *buf = NULL
        cdef bytes obj

        pos = stdio.ftell(self.handle)
        if pos <= 0:
            return ""

        if self._memstream_ptr != NULL:
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
                    raise IOError("failed to read messages from buffer")

                obj = PyBytes_FromStringAndSize(buf, nread)
            finally:
                stdlib.free(buf)

        return obj.decode('latin1')

    def clear(self):
        stdio.rewind(self.handle)

    def close(self):
        if self.handle != NULL:
            stdio.fclose(self.handle)
            self.handle = NULL

        if self._memstream_ptr != NULL:
            stdlib.free(self._memstream_ptr)
            self._memstream_ptr = NULL

        if not self._removed:
            stdio.remove(self._filename)
            self._removed = 1
