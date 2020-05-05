# -*- python -*- or near enough

import sys
import zlib

from cpython cimport PyBytes_FromStringAndSize, \
    PyBytes_AS_STRING, PyBytes_Size

from .pyalloc cimport pyalloc_v

from libc.stdio cimport fread, fseek, ftell
from libc.string cimport memcpy

cdef extern from "Python.h":
    void *PyCObject_Import(char *, char *) except NULL
    ctypedef struct PyTypeObject:
        pass
    ctypedef struct PyObject:
        pass
    ctypedef struct FILE


DEF _BLOCK_SIZE = 131072

BLOCK_SIZE = _BLOCK_SIZE  # public

cdef class GenericStream:

    def __init__(self, fobj):
        self.fobj = fobj

    cpdef int seek(self, long int offset, int whence=0) except -1:
        self.fobj.seek(offset, whence)
        return 0

    cpdef long int tell(self) except -1:
        return self.fobj.tell()

    def read(self, n_bytes):
        return self.fobj.read(n_bytes)

    cpdef int all_data_read(self) except *:
        return 1

    cdef int read_into(self, void *buf, size_t n) except -1:
        """ Read n bytes from stream into pre-allocated buffer `buf`
        """
        cdef char *p
        cdef size_t read_size, count

        # Read data to buf in _BLOCK_SIZE blocks
        count = 0
        p = <char*>buf
        while count < n:
            read_size = min(n - count, _BLOCK_SIZE)
            data = self.fobj.read(read_size)
            read_size = len(data)
            if read_size == 0:
                break
            memcpy(p, <const char*>data, read_size)
            p += read_size
            count += read_size

        if count != n:
            raise IOError('could not read bytes')
        return 0

    cdef object read_string(self, size_t n, void **pp, int copy=True):
        """Make new memory, wrap with object"""
        if copy != True:
            data = self.fobj.read(n)
            if PyBytes_Size(data) != n:
                raise IOError('could not read bytes')
            pp[0] = <void*>PyBytes_AS_STRING(data)
            return data

        cdef object d_copy = pyalloc_v(n, pp)
        self.read_into(pp[0], n)
        return d_copy


cdef class ZlibInputStream(GenericStream):
    """
    File-like object uncompressing bytes from a zlib compressed stream.

    Parameters
    ----------
    stream : file-like
        Stream to read compressed data from.
    max_length : int
        Maximum number of bytes to read from the stream.

    Notes
    -----
    Some matlab files contain zlib streams without valid Z_STREAM_END
    termination.  To get round this, we use the decompressobj object, that
    allows you to decode an incomplete stream.  See discussion at
    https://bugs.python.org/issue8672

    """

    cdef ssize_t _max_length
    cdef object _decompressor
    cdef bytes _buffer
    cdef size_t _buffer_size
    cdef size_t _buffer_position
    cdef size_t _total_position
    cdef size_t _read_bytes

    def __init__(self, fobj, ssize_t max_length):
        self.fobj = fobj

        self._max_length = max_length
        self._decompressor = zlib.decompressobj()
        self._buffer = b''
        self._buffer_size = 0
        self._buffer_position = 0
        self._total_position = 0
        self._read_bytes = 0

    cdef inline void _fill_buffer(self) except *:
        cdef size_t read_size
        cdef bytes block

        if self._buffer_position < self._buffer_size:
            return

        read_size = min(_BLOCK_SIZE, self._max_length - self._read_bytes)

        block = self.fobj.read(read_size)
        self._read_bytes += len(block)

        self._buffer_position = 0
        if not block:
            self._buffer = self._decompressor.flush()
        else:
            self._buffer = self._decompressor.decompress(block)
        self._buffer_size = len(self._buffer)

    cdef int read_into(self, void *buf, size_t n) except -1:
        """Read n bytes from stream into pre-allocated buffer `buf`
        """
        cdef char *dstp
        cdef char *srcp
        cdef size_t read_size, count, size

        dstp = <char*>buf
        count = 0
        while count < n:
            self._fill_buffer()
            if self._buffer_size == 0:
                break

            srcp = <char*>self._buffer
            srcp += self._buffer_position

            size = min(n - count, self._buffer_size - self._buffer_position)
            memcpy(dstp, srcp, size)

            count += size
            dstp += size
            self._buffer_position += size

        self._total_position += count

        if count != n:
            raise IOError('could not read bytes')

        return 0

    cdef object read_string(self, size_t n, void **pp, int copy=True):
        """Make new memory, wrap with object"""
        cdef object d_copy = pyalloc_v(n, pp)
        self.read_into(pp[0], n)
        return d_copy

    def read(self, n_bytes):
        cdef void *p
        return self.read_string(n_bytes, &p)

    cpdef int all_data_read(self) except *:
        if self._read_bytes < self._max_length:
            # we might still have checksum bytes to read
            self._fill_buffer()
        return (self._max_length == self._read_bytes) and \
               (self._buffer_size == self._buffer_position)

    cpdef long int tell(self) except -1:
        if self._total_position == -1:
            raise IOError("Invalid file position.")
        return self._total_position

    cpdef int seek(self, long int offset, int whence=0) except -1:
        cdef ssize_t new_pos, size
        if whence == 1:
            new_pos = <ssize_t>self._total_position + offset
        elif whence == 0:
            new_pos = offset
        elif whence == 2:
            raise IOError("Zlib stream cannot seek from file end")
        else:
            raise ValueError("Invalid value for whence")

        if new_pos < self._total_position:
            raise IOError("Zlib stream cannot seek backwards")

        while self._total_position < new_pos:
            self._fill_buffer()
            if self._buffer_size == 0:
                break

            size = min(new_pos - self._total_position,
                       self._buffer_size - self._buffer_position)

            self._total_position += size
            self._buffer_position += size

        return 0


def _read_into(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef char * d_ptr
    # use bytearray because bytes() is immutable
    my_str = bytearray(b' ' * n)
    d_ptr = my_str
    st.read_into(d_ptr, n)
    return bytes(my_str)


def _read_string(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef void *d_ptr
    cdef object obj = st.read_string(n, &d_ptr, True)
    # use bytearray because bytes() is immutable
    my_str = bytearray(b'A' * n)
    cdef char *mys_ptr = my_str
    memcpy(mys_ptr, d_ptr, n)
    return bytes(my_str)


cpdef GenericStream make_stream(object fobj):
    """ Make stream of correct type for file-like `fobj`
    """
    if isinstance(fobj, GenericStream):
        return fobj
    return GenericStream(fobj)
