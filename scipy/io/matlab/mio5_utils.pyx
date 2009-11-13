''' Cython mio5 routines (-*- python -*- like)

'''

import sys

import cStringIO
import zlib

import numpy as np
cimport numpy as cnp

import scipy.io.matlab.mio5 as mio5

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
    long int ftell (FILE *stream)


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


cdef int miDOUBLE = mio5.miDOUBLE
cdef int miMATRIX = mio5.miMATRIX
cdef int miCOMPRESSED = mio5.miCOMPRESSED

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

# Relationship of matlab types to matrix getter classes
getter_dispatch_table = {
    mio5.mxSPARSE_CLASS: mio5.Mat5SparseMatrixGetter,
    mio5.mxCHAR_CLASS: mio5.Mat5CharMatrixGetter,
    mio5.mxCELL_CLASS: mio5.Mat5CellMatrixGetter,
    mio5.mxSTRUCT_CLASS: mio5.Mat5StructMatrixGetter,
    mio5.mxOBJECT_CLASS: mio5.Mat5ObjectMatrixGetter,
    mio5.mxFUNCTION_CLASS: mio5.Mat5FunctionGetter}
for _mc in mio5.mx_numbers:
    getter_dispatch_table[_mc] = mio5.Mat5NumericMatrixGetter


cdef class CReader:
    cdef public int is_swapped, little_endian
    cdef public object mat_stream,
    cdef object stream_type, dtypes, codecs
    cdef FILE* _file
    cdef object preader
    
    def __new__(self, preader):
        self.preader = preader
        self.mat_stream = preader.mat_stream
        self.dtypes = preader.dtypes
        self.codecs = preader.codecs
        # infer endianness from double dtype
        f8_code = self.dtypes[miDOUBLE].byteorder
        self.is_swapped = f8_code == swapped_code
        if self.is_swapped:
            self.little_endian = not sys_is_le
        else:
            self.little_endian = sys_is_le
        fobj = self.mat_stream
        if isinstance(fobj, file):
            self.stream_type = real_file
            self._file = PyFile_AsFile(fobj)
        elif PycStringIO_InputCheck(fobj) or PycStringIO_OutputCheck(fobj):
            self.stream_type = cstringio
        else:
            self.stream_type = generic
            
    def read_element(self, int copy=True):
        cdef cnp.uint16_t mdtype_sde, byte_count_sde
        cdef cnp.uint32_t mdtype, byte_count, bc_in
        cdef int mod8
        cdef char* data_ptr
        cdef int read_ret
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
        read_ret = self.get_2_u4s(&mdtype, &bc_in)
        if self.is_swapped:
            mdtype = byteswap_u4(mdtype)
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
            data_ptr = <char *>&bc_in
            mdtype = mdtype_sde
            byte_count = byte_count_sde
        else: # regular element
            if self.is_swapped:
                byte_count = byteswap_u4(bc_in)
            else:
                byte_count = bc_in
            # Deal with miMATRIX type (cannot pass byte string)
            if mdtype == miMATRIX:
                return self.current_getter(byte_count).get_array()
            # All other types can be read from string
            data = self.mat_stream.read(byte_count)
            data_ptr = data
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.move_forward(8 - mod8)

        if mdtype in self.codecs: # encoded char data
            codec = self.codecs[mdtype]
            if not codec:
                raise TypeError('Do not support encoding %d' % mdtype)
            tmp = data_ptr[:byte_count]
            el = tmp.decode(codec)
        else: # numeric data
            dt = self.dtypes[mdtype]
            el_count = byte_count // dt.itemsize
            el = np.ndarray(shape=(el_count,),
                            dtype=dt,
                            buffer=data_ptr[:byte_count])
            if copy:
                el = el.copy()
        return el

    cdef int get_2_u4s(self, cnp.uint32_t* v1, cnp.uint32_t* v2):
        ''' Get 2 as-yet unswapped uint32 values from stream '''
        cdef char* tag_ptr
        cdef char tag_bytes[8]
        cdef cnp.uint32_t* u4_ptr
        cdef size_t n_red
        if self.stream_type == real_file: # really a file object
            n_red = fread(tag_bytes, 1, 8, self._file)
            if n_red < 8:
                return -1
            tag_ptr = <char *>tag_bytes
        elif self.stream_type == cstringio:
            n_red = PycStringIO.cread(self.mat_stream, &tag_ptr, 8)
            if n_red < 8:
                return -1
        else:# use python interface to get data
            # python variable needed to hold bytes
            tag = self.mat_stream.read(8)
            tag_ptr = tag
            if len(tag) < 8:
                return -1
        u4_ptr = <cnp.uint32_t *>tag_ptr
        v1[0] = u4_ptr[0]
        v2[0] = u4_ptr[1]
        return 1

    cdef void move_forward(self, size_t n):
        ''' Move forward in stream n bytes '''
        cdef char* ptr
        if self.stream_type == real_file: # really a file object
            fseek(self._file, n, 1)
        elif self.stream_type == cstringio:
            PycStringIO.cread(self.mat_stream, &ptr, n)
        else:# use python interface
            self.mat_stream.seek(n, 1)

    cdef size_t tell(self):
        if self.stream_type == real_file: # really a file object
            return ftell(self._file)
        else:# use python interface
            self.mat_stream.tell()

    def matrix_getter_factory(self):
        ''' Returns reader for next matrix at top level '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef int read_ret
        cdef size_t next_pos
        read_ret = self.get_2_u4s(&mdtype, &byte_count)
        if self.is_swapped:
            mdtype = byteswap_u4(mdtype)
            byte_count = byteswap_u4(byte_count)
        next_pos = self.mat_stream.tell() + byte_count
        if mdtype == miCOMPRESSED:
            array_reader = self.preader
            stream = cStringIO.StringIO(
                zlib.decompress(self.mat_stream.read(byte_count)))
            reader = mio5.Mat5ArrayReader(
                stream,
                array_reader.dtypes,
                array_reader.processor_func,
                array_reader.codecs,
                array_reader.class_dtypes,
                array_reader.struct_as_record)
            getter = reader.matrix_getter_factory()
        elif not mdtype == miMATRIX:
            raise TypeError, \
                  'Expecting miMATRIX type here, got %d' %  mdtype
        else:
            getter = self.current_getter(byte_count)
        getter.next_position = next_pos
        return getter

    def current_getter(self, size_t byte_count):
        ''' Return matrix getter for current stream position

        Returns matrix getters at top level and sub levels
        '''
        cdef cnp.uint32_t af_mdtype, af_byte_count
        cdef cnp.uint32_t flags_class, nzmax
        cdef cnp.uint16_t mc
        if not byte_count: # an empty miMATRIX can contain no bytes
            return mio5.Mat5EmptyMatrixGetter(self)
        # Read and discard mdtype and byte_count
        self.get_2_u4s(&af_mdtype, &af_byte_count)
        # get array flags and nzmax
        self.get_2_u4s(&flags_class, &nzmax)
        if self.is_swapped:
            flags_class = byteswap_u4(flags_class)
            nzmax = byteswap_u4(nzmax)
        header = {}
        mc = flags_class & 0xFF
        header['mclass'] = mc
        header['is_logical'] = flags_class >> 9 & 1
        header['is_global'] = flags_class >> 10 & 1
        header['is_complex'] = flags_class >> 11 & 1
        header['nzmax'] = nzmax
        ''' Here I am playing with a binary block read of
        untranslatable data. I am not using this at the moment because
        reading it has the side effect of making opposite ending mat
        files unwritable on the round trip.
        
        if mc == mxFUNCTION_CLASS:
            # we can't read these, and want to keep track of the byte
            # count - so we need to avoid the following unpredictable
            # length element reads
            return Mat5BinaryBlockGetter(self,
                                         header,
                                         af,
                                         byte_count)
        '''
        header['dims'] = self.read_element()
        header['name'] = self.read_element().tostring()
        return getter_dispatch_table[mc](self.preader, header)
