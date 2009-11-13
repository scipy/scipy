''' Cython mio5 routines (-*- python -*- like)

'''

import sys

import numpy as np
cimport numpy as cnp

cdef extern from "Python.h":
    void *PyCObject_Import(char *, char *) except NULL
    ctypedef struct PyTypeObject:
        pass
    ctypedef struct PyObject:
        pass
    ctypedef struct FILE
    FILE *PyFile_AsFile(object)
    size_t fread (void *ptr, size_t size, size_t n, FILE* fptr)
    int fseek (FILE * fptr, long int offset, int whence)


cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass


cdef extern from "cStringIO.h":
    struct PycStringIO_CAPI:
        int (*cread)(object Input, char ** bufptr, Py_ssize_t n)
        int (*creadline)(object Input, char ** bufptr)
        int (*cwrite)(object Output, char * out_str, Py_ssize_t n)
        object (*cgetvalue)(object Input)
        object (*NewOutput)(int size)
        PyObject *(*NewInput)(object in_str)
        PyTypeObject *InputType, *OutputType
    cdef PycStringIO_CAPI* PycStringIO
    bint PycStringIO_InputCheck(object O)
    bint PycStringIO_OutputCheck(object O)


# initialize cStringIO
PycStringIO = <PycStringIO_CAPI*>PyCObject_Import("cStringIO", "cStringIO_CAPI")


cdef enum FileType:
    generic
    real_file
    cstringio


cdef int miDOUBLE = 9
cdef int miMATRIX = 14

sys_is_le = sys.byteorder == 'little'
native_code = sys_is_le and '<' or '>'
swapped_code = sys_is_le and '>' or '<'


def bsu4(u4): # for testing
    return byteswap_u4(u4)
    
            
cdef cnp.uint32_t byteswap_u4(cnp.uint32_t u4):
    return ((u4 << 24) |
           ((u4 << 8) & 0xff0000U) |
           ((u4 >> 8 & 0xff00u)) |
           (u4 >> 24))


cdef class CReader:
    cdef public int is_swapped, little_endian
    cdef int stream_type
    cdef object mat_stream, dtypes, codecs, current_getter
    cdef FILE* _file
    cdef object sreader
    
    def __new__(self, preader):
        self.mat_stream = preader.mat_stream
        self.dtypes = preader.dtypes
        self.codecs = preader.codecs
        self.current_getter = preader.current_getter
        # infer endianness from double dtype
        f8_code = self.dtypes[miDOUBLE].byteorder
        self.is_swapped = f8_code == swapped_code
        if self.is_swapped:
            self.little_endian = not sys_is_le
        else:
            self.little_endian = sys_is_le
        self.stream_type = generic
        if isinstance(self.mat_stream, file):
            self.stream_type = real_file
            self._file = PyFile_AsFile(self.mat_stream)
        
            
    def read_element(self, int copy=True):
        cdef cnp.uint16_t mdtype_sde, byte_count_sde
        cdef cnp.uint32_t mdtype, byte_count
        cdef int mod8
        cdef char* tag_ptr
        cdef char tag_bytes[8]
        cdef char* data_ptr
        cdef cnp.uint32_t* u4_ptr
        # First read 8 bytes.  The 8 bytes can be in one of two formats.
        # For the first - standard format - the 8 bytes are two uint32
        # values, of which the first is the integer code for the matlab
        # data type (*mdtype*), and the second is the number of bytes of
        # that data type that follow (*byte_count*).  Thus, if the
        # ``mdtype`` is 4 (miDOUBLE), and the ``byte_count`` is 12, then
        # there will follow 3 double values.  The alternative format is
        # "small data element". The first four bytes contain the
        # ``byte_count`` and the ``mdtype``, but as uint16.  The
        # arrangement of the ``byte_count`` and ``mdtype`` is a little
        # complex, see below. The following 4 bytes of the 8 bytes
        # contain the data.  For example, the ``mdtype`` might be 2
        # (miUINT8), and the byte count is 3, and the data is in a
        # string ``tag``, then the contained matrix is length 3, type
        # uint8, where values are ``tag[4], tag[5], tag[6]``.
        #
        # The following paragraph describes the extraction of ``mdtype``
        # and ``byte_count`` for the small data element format.  The
        # following is somewhat contrary to the matlab documentation,
        # but seems to be true of actual .mat files.
        #
        # If the *file* is big endian, then the first four bytes of the
        # tag are two big-endian uint16 values, first ``byte_count`` and
        # second ``mdtype``.  If the *file* is little-endian then the
        # first four bytes are two little-endian uint16 values, first
        # ``mdtype`` and second ``byte_count``.   
        if self.stream_type == real_file: # really a file object
            fread(tag_bytes, 8, 1, self._file)
            tag_ptr = <char *>tag_bytes
        else: # use python interface to get data
            # python variable needed to hold bytes
            tag = self.mat_stream.read(8)
            tag_ptr = tag
        u4_ptr = <cnp.uint32_t *>tag_ptr
        if self.is_swapped:
            mdtype = byteswap_u4(u4_ptr[0])
        else:
            mdtype = u4_ptr[0]
        # The most significant two bytes of a U4 *mdtype* will always be
        # 0, if they are not, this must be SDE format
        byte_count_sde = mdtype >> 16
        if byte_count_sde: # small data element format
            mdtype_sde = mdtype & 0xffff
            if byte_count_sde > 4:
                raise ValueError('Too many bytes (%d) for sde format'
                                 % byte_count_sde)
            if mdtype_sde == miMATRIX:
                raise TypeError('Cannot have matrix in SDE format')
            data_ptr = tag_ptr + 4
            mdtype = mdtype_sde
            byte_count = byte_count_sde
        else: # regular element
            byte_count = u4_ptr[1]
            if self.is_swapped:
                byte_count = byteswap_u4(byte_count)
            # Deal with miMATRIX type (cannot pass byte string)
            if mdtype == miMATRIX:
                return self.current_getter(byte_count).get_array()
            # All other types can be read from string
            data = self.mat_stream.read(byte_count)
            data_ptr = data
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.mat_stream.seek(8 - mod8, 1)

        if mdtype in self.codecs: # encoded char data
            codec = self.codecs[mdtype]
            if not codec:
                raise TypeError('Do not support encoding %d' % mdtype)
            tmp = data_ptr[:byte_count]
            el = tmp.decode(codec)
        else: # numeric data
            if mdtype not in self.dtypes:
                raise ValueError(self.is_swapped)
            dt = self.dtypes[mdtype]
            el_count = byte_count // dt.itemsize
            el = np.ndarray(shape=(el_count,),
                            dtype=dt,
                            buffer=data_ptr[:byte_count])
            if copy:
                el = el.copy()
        return el
