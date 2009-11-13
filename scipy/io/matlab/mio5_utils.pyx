''' Cython mio5 routines (-*- python -*- like)

'''

import sys

import cStringIO
import zlib
from copy import copy as pycopy

import numpy as np
cimport numpy as cnp

import scipy.io.matlab.miobase as miob
import scipy.io.matlab.mio5 as mio5
import scipy.sparse

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
cdef int miINT32 = mio5.miINT32
cdef int miINT8 = mio5.miINT8
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


cdef size_t cproduct(tup):
    cdef size_t res = 1
    cdef int i
    for i in range(len(tup)):
        res *= tup[i]
    return res
           

cdef class ElementHeader:
    cdef int mclass
    cdef int is_logical
    cdef public int is_global
    cdef int is_complex
    cdef size_t nzmax
    cdef object dims
    cdef public object name


cdef class CReader:
    cdef public int is_swapped, little_endian
    cdef readonly object mat_stream, 
    cdef object stream_type, dtypes, class_dtypes, codecs
    cdef int struct_as_record
    cdef FILE* _file
    cdef object preader
    
    def __new__(self, preader):
        self.preader = preader
        self.dtypes = preader.dtypes
        self.class_dtypes = preader.class_dtypes
        self.codecs = preader.codecs
        self.struct_as_record = preader.struct_as_record
        # infer endianness from double dtype
        f8_code = self.dtypes[miDOUBLE].byteorder
        self.is_swapped = f8_code == swapped_code
        if self.is_swapped:
            self.little_endian = not sys_is_le
        else:
            self.little_endian = sys_is_le
        self.set_stream(preader.mat_stream)
        
    def set_stream(self, fobj):
        self.mat_stream = fobj
        fobj = self.mat_stream
        if isinstance(fobj, file):
            self.stream_type = real_file
            self._file = PyFile_AsFile(fobj)
        elif PycStringIO_InputCheck(fobj) or PycStringIO_OutputCheck(fobj):
            self.stream_type = cstringio
        else:
            self.stream_type = generic

    def read_tag(self):
        cdef cnp.uint32_t mdtype, byte_count
        cdef char tag_data[4]
        cdef int tag_res
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_data)
        if tag_res < 1:
            if tag_res == -2:
                raise IOError('Error in SDE format data')
            raise IOError('Error reading data from stream')
        if tag_res == 2: # sde format
            pdata = tag_data
        else:
            pdata = None
        return (mdtype, byte_count, tag_data)

    cdef int cread_tag(self,
                     cnp.uint32_t *mdtype_ptr,
                     cnp.uint32_t *byte_count_ptr,
                     char *data_ptr):
        ''' Read tag mdtype and byte_count

        Does necessary swapping and takes account of SDE formats

        Data may be returned in data_ptr, if this was an SDE

        Returns 1 for success, full format; 2 for success, SDE format; -1
        for failed data read; -2 for error detected in SDE processing
        '''
        cdef cnp.uint16_t mdtype_sde, byte_count_sde
        cdef cnp.uint32_t mdtype, byte_count
        cdef cnp.uint32_t* u4_ptr = <cnp.uint32_t*>data_ptr
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
        read_ret = self.get_2_u4s(&mdtype, &byte_count)
        if read_ret == -1:
            return -1 # 
        if self.is_swapped:
            mdtype = byteswap_u4(mdtype)
        # The most significant two bytes of a U4 *mdtype* will always be
        # 0, if they are not, this must be SDE format
        byte_count_sde = mdtype >> 16
        if byte_count_sde: # small data element format
            mdtype_sde = mdtype & 0xffff
            if byte_count_sde > 4 or mdtype_sde == miMATRIX:
                return -2
            u4_ptr[0] = byte_count
            mdtype_ptr[0] = mdtype_sde
            byte_count_ptr[0] = byte_count_sde
            return 2
        # regular element
        if self.is_swapped:
            byte_count_ptr[0] = byteswap_u4(byte_count)
        else:
            byte_count_ptr[0] = byte_count
        mdtype_ptr[0] = mdtype
        u4_ptr[0] = 0
        return 1
            
    def read_element(self, int copy=True):
        cdef cnp.uint32_t mdtype, byte_count
        cdef int mod8
        cdef char tag_data[4]
        cdef char* data_ptr
        cdef int tag_res
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_data)
        if tag_res < 1:
            if tag_res == -2:
                raise IOError('Error in SDE format data')
            raise IOError('Error reading data from stream')
        if tag_res == 1: # full format
            if mdtype == miMATRIX:
                raise TypeError('Not expecting matrix here')
            # All other types can be read from string
            data = self.mat_stream.read(byte_count)
            data_ptr = data
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.move_forward(8 - mod8)
        else: # SDE format
            data_ptr = tag_data
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

    def read_dims(self):
        ''' Read array dimensions as list '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef cnp.int32_t* dims_ptr
        cdef cnp.int32_t dim
        cdef int n_dims, i
        cdef int mod8
        cdef char tag_data[4], *data_ptr
        cdef int tag_res
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_data)
        if tag_res < 1:
            if tag_res == -2:
                raise IOError('Error in SDE format data')
            raise IOError('Error reading data from stream')
        if mdtype != miINT32:
            raise TypeError('Expecting miINT32 as data type')
        if tag_res == 1: # full format
            # All other types can be read from string
            data = self.mat_stream.read(byte_count)
            data_ptr = data
            dims_ptr = <cnp.int32_t *>data_ptr
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.move_forward(8 - mod8)
        else: # SDE format
            dims_ptr = <cnp.int32_t *>tag_data
        n_dims = byte_count / 4
        dims = []
        for i in range(n_dims):
            dim = dims_ptr[i]
            if self.is_swapped:
                dim = byteswap_u4(dim)
            dims.append(dim)
        return dims
    
    def read_name(self):
        ''' Read matrix name as python string '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef int mod8
        cdef char tag_data[4], *data_ptr
        cdef int tag_res
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_data)
        if tag_res < 1:
            if tag_res == -2:
                raise IOError('Error in SDE format data')
            raise IOError('Error reading data from stream')
        if mdtype != miINT8:
            raise TypeError('Expecting miINT8 as data type')
        if tag_res == 1: # full format
            # All other types can be read from string
            data = self.mat_stream.read(byte_count)
            data_ptr = data
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.move_forward(8 - mod8)
        else: # SDE format
            data_ptr = tag_data
        return data_ptr[:byte_count]

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

    def read_full_tag(self):
        cdef cnp.uint32_t mdtype, byte_count
        cdef size_t next_pos
        read_ret = self.get_2_u4s(&mdtype, &byte_count)
        if self.is_swapped:
            mdtype = byteswap_u4(mdtype)
            byte_count = byteswap_u4(byte_count)
        return mdtype, byte_count

    def get_mi_matrix(self):
        mdtype, byte_count = self.read_full_tag()
        if mdtype != mio5.miMATRIX:
            raise TypeError('Expecting matrix here')
        if byte_count == 0: # empty matrix
            mat_dtype = 'f8'
            arr = np.array([[]])
        else:
            header = self.read_header()
            arr, mdtype = self.get_array(header)
        return miob.process_element(arr, self.preader, mat_dtype)

    def read_header(self):
        ''' Return matrix header for current stream position

        Returns matrix headers at top level and sub levels
        '''
        cdef:
            cnp.uint32_t af_mdtype, af_byte_count
            cnp.uint32_t flags_class, nzmax
            cnp.uint16_t mc
            int ret
            ElementHeader header
        # Read and discard mdtype and byte_count
        ret = self.get_2_u4s(&af_mdtype, &af_byte_count)
        if ret < 0:
            raise IOError('Cannot read from stream')
        # get array flags and nzmax
        ret = self.get_2_u4s(&flags_class, &nzmax)
        if ret < 0:
            raise IOError('Cannot read from stream')
        if self.is_swapped:
            flags_class = byteswap_u4(flags_class)
            nzmax = byteswap_u4(nzmax)
        header = ElementHeader()
        mc = flags_class & 0xFF
        header.mclass = mc
        header.is_logical = flags_class >> 9 & 1
        header.is_global = flags_class >> 10 & 1
        header.is_complex = flags_class >> 11 & 1
        header.nzmax = nzmax
        header.dims = self.read_dims()
        header.name = self.read_name()
        return header

    def get_array(self, ElementHeader header):
        cdef:
            int mc
            size_t length
        mc = header.mclass
        mat_dtype = None
        if mc in mio5.mx_numbers: # numeric
            if header.is_logical:
                mat_dtype = np.dtype('bool')
            else:
                mat_dtype = self.class_dtypes[mc]
            if header.is_complex:
                # avoid array copy to save memory
                res = self.read_element(copy=False)
                res_j = self.read_element(copy=False)
                res = res + (res_j * 1j)
            else:
                res = self.read_element()
            arr = np.ndarray(shape=header.dims,
                              dtype=res.dtype,
                              buffer=res,
                              order='F')
        elif mc == mio5.mxSPARSE_CLASS:
            rowind = self.read_element()
            indptr = self.read_element()
            if header.is_complex:
                # avoid array copy to save memory
                data   = self.read_element(copy=False)
                data_j = self.read_element(copy=False)
                data = data + (data_j * 1j)
            else:
                data = self.read_element()
            ''' From the matlab (TM) API documentation, last found here:
            http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/
            rowind are simply the row indices for all the (nnz) non-zero
            entries in the sparse array.  rowind has nzmax entries, so
            may well have more entries than nnz, the actual number of
            non-zero entries, but rowind[nnz:] can be discarded and
            should be 0. indptr has length (number of columns + 1), and
            is such that, if D = diff(colind), D[j] gives the number of
            non-zero entries in column j. Because rowind values are
            stored in column order, this gives the column corresponding
            to each rowind
            '''
            M,N = header.dims
            indptr = indptr[:N+1]
            nnz = indptr[-1]
            rowind = rowind[:nnz]
            data   = data[:nnz]
            arr = scipy.sparse.csc_matrix(
                (data,rowind,indptr),
                shape=(M,N))
        elif mc == mio5.mxCHAR_CLASS:
            res = self.read_element()
            # Convert non-string types to unicode
            if isinstance(res, np.ndarray):
                if res.dtype.type == np.uint16:
                    codec = mio5.miUINT16_codec
                    if self.codecs['uint16_len'] == 1:
                        res = res.astype(np.uint8)
                elif res.dtype.type in (np.uint8, np.int8):
                    codec = 'ascii'
                else:
                    raise TypeError, 'Did not expect type %s' % res.dtype
                res = res.tostring().decode(codec)
            arr = np.ndarray(shape=header.dims,
                             dtype=np.dtype('U1'),
                             buffer=np.array(res),
                             order='F').copy()
        elif mc == mio5.mxCELL_CLASS:
            arr = self._read_cells(header)
        elif mc == mio5.mxSTRUCT_CLASS:
            arr = self._read_struct(header)
        elif mc == mio5.mxOBJECT_CLASS: # like structs, but with classname
            classname = self.read_element().tostring()
            arr = self._read_struct(header)
            arr = mio5.MatlabObject(result, classname)
        elif mc == mio5.mxFUNCTION_CLASS:
            raise mio5.MatReadError('Cannot read matlab functions')
        return arr, mat_dtype

    def _read_cells(self, ElementHeader header):
        cdef cnp.ndarray[object, ndim=1] result
        # Account for fortran indexing of cells
        tupdims = tuple(header.dims[::-1])
        length = cproduct(tupdims)
        result = np.empty(length, dtype=object)
        for i in range(length):
            result[i] = self.get_mi_matrix()
        return result.reshape(tupdims).T
        
    def _read_struct(self, ElementHeader header):
        cdef:
            int namelength, length, i, j, n_names
            cnp.ndarray[object, ndim=1] result
        namelength = self.read_element()[0]
        names = self.read_element()
        field_names = []
        n_names = len(names)
        for i from 0 <= i < n_names by namelength:
            name = names[i:i+namelength].tostring().strip('\x00')
            field_names.append(name)
        tupdims = tuple(header.dims[::-1])
        length = cproduct(tupdims)
        if self.struct_as_record: # to record arrays
            if not len(field_names):
                # If there are no field names, there is no dtype
                # representation we can use, falling back to empty
                # object
                return np.empty(tupdims, dtype=object).T
            dtype = [(field_name, object) for field_name in field_names]
            rec_res = np.empty(length, dtype=dtype)
            for i in range(length):
                for field_name in field_names:
                    rec_res[i][field_name] = self.get_mi_matrix()
            return rec_res.reshape(tupdims).T
        # Backward compatibility with previous format
        obj_template = mio5.mat_struct()
        obj_template._fieldnames = field_names
        result = np.empty(length, dtype=object)
        for i in range(length):
            item = pycopy(obj_template)
            for name in field_names:
                item.__dict__[name] = self.get_mi_matrix()
            result[i] = item
        return result.reshape(tupdims).T


