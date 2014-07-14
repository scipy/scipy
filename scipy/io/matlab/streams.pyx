# -*- python -*- or near enough

import sys
import zlib

from cpython cimport PyBytes_FromStringAndSize, \
    PyBytes_AS_STRING, PyBytes_Size

from pyalloc cimport pyalloc_v

cdef extern from "stdlib.h" nogil:
    void *malloc(size_t size)
    void *memcpy(void *str1, void *str2, size_t n)
    void free(void *ptr)

cdef extern from "Python.h":
    void *PyCObject_Import(char *, char *) except NULL
    ctypedef struct PyTypeObject:
        pass
    ctypedef struct PyObject:
        pass
    ctypedef struct FILE
    size_t fread (void *ptr, size_t size, size_t n, FILE* fptr)
    int fseek (FILE * fptr, long int offset, int whence)
    long int ftell (FILE *stream)

cdef extern from "py3k.h":
    # From:
    # http://svn.pyamf.org/pyamf/tags/release-0.4rc2/cpyamf/util.pyx
    # (MIT license) - with thanks
    void PycString_IMPORT()
    int StringIO_cread "PycStringIO->cread" (object, char **, Py_ssize_t)
    int StringIO_creadline "PycStringIO->creadline" (object, char **)
    int StringIO_cwrite "PycStringIO->cwrite" (object, char *, Py_ssize_t)
    object StringIO_cgetvalue "PycStringIO->cgetvalue" (obj)
    bint PycStringIO_InputCheck(object O)
    bint PycStringIO_OutputCheck(object O)

    FILE* npy_PyFile_Dup(object file, char *mode) except NULL
    int npy_PyFile_DupClose(object file, FILE *handle) except -1
    int npy_PyFile_Check(object file)

       
# initialize cStringIO
PycString_IMPORT

DEF BLOCK_SIZE = 131072

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

    cdef int read_into(self, void *buf, size_t n) except -1:
        """ Read n bytes from stream into pre-allocated buffer `buf`
        """
        cdef char *p
        cdef size_t read_size, count

        # Read data to buf in BLOCK_SIZE blocks
        count = 0
        p = <char*>buf
        while count < n:
            read_size = min(n - count, BLOCK_SIZE)
            data = self.fobj.read(read_size)
            read_size = len(data)
            if read_size == 0:
                break
            memcpy(p, <char*>data, read_size)
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
    http://bugs.python.org/issue8672

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

    cdef _fill_buffer(self):
        cdef size_t read_size
        cdef bytes block

        if self._buffer_position < self._buffer_size:
            return

        read_size = min(BLOCK_SIZE, self._max_length - self._read_bytes)

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

    cpdef int all_data_read(self):
        return (self._max_length == self._read_bytes) and \
               (self._buffer_size == self._buffer_position)

    cpdef long int tell(self):
        return self._total_position

    cpdef int seek(self, long int offset, int whence=0) except -1:
        if whence == 1:
            new_pos = self._total_position + offset
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


cdef class cStringStream(GenericStream):
    
    cpdef int seek(self, long int offset, int whence=0) except -1:
        cdef char *ptr
        if whence == 1 and offset >=0: # forward, from here
            StringIO_cread(self.fobj, &ptr, offset)
            return 0
        else: # use python interface
            return GenericStream.seek(self, offset, whence)

    cdef int read_into(self, void *buf, size_t n) except -1:
        """ Read n bytes from stream into pre-allocated buffer `buf`
        """
        cdef:
            size_t n_red
            char* d_ptr
        n_red = StringIO_cread(self.fobj, &d_ptr, n)
        if n_red != n:
            raise IOError('could not read bytes')
        memcpy(buf, <void *>d_ptr, n)
        return 0

    cdef object read_string(self, size_t n, void **pp, int copy=True):
        """ Make new memory, wrap with object

        It's not obvious to me how to avoid a copy
        """
        cdef:
            char *d_ptr
            object obj
        cdef size_t n_red = StringIO_cread(self.fobj, &d_ptr, n)
        if n_red != n:
            raise IOError('could not read bytes')
        obj = pyalloc_v(n, pp)
        memcpy(pp[0], d_ptr, n)
        return obj

   
cdef class FileStream(GenericStream):
    cdef FILE* file

    def __init__(self, fobj):
        self.fobj = fobj
        self.file = npy_PyFile_Dup(fobj, "rb")

    def __del__(self):
        npy_PyFile_DupClose(self.fobj, self.file)

    cpdef int seek(self, long int offset, int whence=0) except -1:
        cdef int ret
        """ move `offset` bytes in stream

        Parameters
        ----------
        offset : long int
           number of bytes to move.  Positive for forward in file,
           negative for backward
        whence : int
           `whence` can be:
           
           * 0 - from beginning of file (`offset` should be >=0)
           * 1 - from current file position
           * 2 - from end of file (`offset` nearly always <=0)

        Returns
        -------
        ret : int
        """
        ret = fseek(self.file, offset, whence)
        if ret:
            raise IOError('Failed seek')
        return ret

    cpdef long int tell(self):
        return ftell(self.file)

    cdef int read_into(self, void *buf, size_t n) except -1:
        """ Read n bytes from stream into pre-allocated buffer `buf`
        """
        cdef:
            size_t n_red
            char* d_ptr
        n_red = fread(buf, 1, n, self.file)
        if n_red != n:
            raise IOError('Could not read bytes')
        return 0

    cdef object read_string(self, size_t n, void **pp, int copy=True):
        """ Make new memory, wrap with object """
        cdef object obj = pyalloc_v(n, pp)
        cdef size_t n_red = fread(pp[0], 1, n, self.file)
        if n_red != n:
            raise IOError('could not read bytes')
        return obj

def _read_into(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef char * d_ptr
    my_str = b' ' * n
    d_ptr = my_str
    st.read_into(d_ptr, n)
    return my_str


def _read_string(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef char *d_ptr
    cdef object obj = st.read_string(n, <void **>&d_ptr, True)
    my_str = b'A' * n
    cdef char *mys_ptr = my_str
    memcpy(mys_ptr, d_ptr, n)
    return my_str


cpdef GenericStream make_stream(object fobj):
    """ Make stream of correct type for file-like `fobj`
    """
    if npy_PyFile_Check(fobj):
        if sys.version_info[0] >= 3:
            return GenericStream(fobj)
        else:
            return FileStream(fobj)
    elif PycStringIO_InputCheck(fobj) or PycStringIO_OutputCheck(fobj):
        return cStringStream(fobj)
    elif isinstance(fobj, GenericStream):
        return fobj
    return GenericStream(fobj)


