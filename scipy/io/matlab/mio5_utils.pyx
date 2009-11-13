''' Cython mio5 utility routines (-*- python -*- like)

'''

import sys

from copy import copy as pycopy

import numpy as np
cimport numpy as cnp

# Constant from numpy - max number of array dimensions
DEF _MAT_MAXDIMS = 32

# Max length of matlab variable names - in fact its 63 at the moment
DEF _MAT_MAX_NAME_LEN = 128

cimport streams
import scipy.io.matlab.miobase as miob
from scipy.io.matlab.mio_utils import FileReadOpts, process_element
cimport mio_utils as cmio_utils
import scipy.io.matlab.mio5_params as mio5p
import scipy.sparse


cdef extern from "stdlib.h" nogil:
    void *memcpy(void *str1, void *str2, size_t n)


cdef:
    int miDOUBLE = mio5p.miDOUBLE
    int miINT32 = mio5p.miINT32
    int miINT8 = mio5p.miINT8
    int miMATRIX = mio5p.miMATRIX
    int miCOMPRESSED = mio5p.miCOMPRESSED
    # Numeric matrix classes
    int mxDOUBLE_CLASS = mio5p.mxDOUBLE_CLASS
    int mxSINGLE_CLASS = mio5p.mxSINGLE_CLASS
    int mxINT8_CLASS = mio5p.mxINT8_CLASS
    int mxUINT8_CLASS = mio5p.mxUINT8_CLASS
    int mxINT16_CLASS = mio5p.mxINT16_CLASS
    int mxUINT16_CLASS = mio5p.mxUINT16_CLASS
    int mxINT32_CLASS = mio5p.mxINT32_CLASS
    int mxUINT32_CLASS = mio5p.mxUINT32_CLASS
    int mxINT64_CLASS = mio5p.mxINT64_CLASS
    int mxUINT64_CLASS = mio5p.mxUINT64_CLASS
    # other matrix classes
    int mxSPARSE_CLASS = mio5p.mxSPARSE_CLASS
    int mxCHAR_CLASS = mio5p.mxCHAR_CLASS
    int mxCELL_CLASS = mio5p.mxCELL_CLASS
    int mxSTRUCT_CLASS = mio5p.mxSTRUCT_CLASS
    int mxOBJECT_CLASS = mio5p.mxOBJECT_CLASS
    int mxFUNCTION_CLASS = mio5p.mxFUNCTION_CLASS


sys_is_le = sys.byteorder == 'little'
native_code = sys_is_le and '<' or '>'
swapped_code = sys_is_le and '>' or '<'


cpdef cnp.uint32_t byteswap_u4(cnp.uint32_t u4):
    return ((u4 << 24) |
           ((u4 << 8) & 0xff0000U) |
           ((u4 >> 8 & 0xff00u)) |
           (u4 >> 24))


cdef class VarHeader5:
    cdef char name_ptr[_MAT_MAX_NAME_LEN]
    cdef int name_len
    cdef int mclass
    cdef object dims
    cdef cnp.int32_t dims_ptr[_MAT_MAXDIMS]
    cdef int n_dims
    cdef int is_complex
    cdef int is_logical
    cdef public int is_global
    cdef size_t nzmax

    property name:
        def __get__(self):
            return self.name_ptr[:self.name_len]



cdef class VarReader5:
    cdef public int is_swapped, little_endian
    cdef readonly object mat_stream
    cdef object dtypes, class_dtypes, codecs, uint16_codec
    cdef int struct_as_record
    cdef object preader
    cdef cmio_utils.FileReadOpts read_opts
    cdef streams.GenericStream cstream
    
    def __new__(self, preader):
        self.preader = preader
        self.dtypes = preader.dtypes
        self.class_dtypes = preader.class_dtypes
        self.codecs = preader.codecs
        self.struct_as_record = preader.struct_as_record
        self.uint16_codec = preader.uint16_codec
        self.is_swapped = preader.byte_order == swapped_code
        if self.is_swapped:
            self.little_endian = not sys_is_le
        else:
            self.little_endian = sys_is_le
        self.set_stream(preader.mat_stream)
        self.read_opts = FileReadOpts(
            preader.chars_as_strings,
            preader.mat_dtype,
            preader.squeeze_me)
        
    def set_stream(self, fobj):
        self.mat_stream = fobj
        self.cstream = streams.make_stream(fobj)
        
    def read_tag(self):
        cdef cnp.uint32_t mdtype, byte_count
        cdef char tag_data[4]
        cdef int tag_res
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_data)
        if tag_res == 2: # sde format
            pdata = tag_data
        else:
            pdata = None
        return (mdtype, byte_count, pdata)

    cdef int cread_tag(self,
                     cnp.uint32_t *mdtype_ptr,
                     cnp.uint32_t *byte_count_ptr,
                     char *data_ptr) except -1:
        ''' Read tag mdtype and byte_count

        Does necessary swapping and takes account of SDE formats

        Data may be returned in data_ptr, if this was an SDE

        Returns 1 for success, full format; 2 for success, SDE format; -1
        if error arises
        '''
        cdef cnp.uint16_t mdtype_sde, byte_count_sde
        cdef cnp.uint32_t mdtype
        cdef cnp.uint32_t* u4_ptr = <cnp.uint32_t*>data_ptr
        cdef cnp.uint32_t u4s[2]
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
        self.cstream.read_into(<void *>u4s, 8)
        if self.is_swapped:
            mdtype = byteswap_u4(u4s[0])
        else:
            mdtype = u4s[0]
        # The most significant two bytes of a U4 *mdtype* will always be
        # 0, if they are not, this must be SDE format
        byte_count_sde = mdtype >> 16
        if byte_count_sde: # small data element format
            mdtype_sde = mdtype & 0xffff
            if byte_count_sde > 4:
                raise ValueError('Error in SDE format data')
                return -1
            u4_ptr[0] = u4s[1]
            mdtype_ptr[0] = mdtype_sde
            byte_count_ptr[0] = byte_count_sde
            return 2
        # regular element
        if self.is_swapped:
            byte_count_ptr[0] = byteswap_u4(u4s[1])
        else:
            byte_count_ptr[0] = u4s[1]
        mdtype_ptr[0] = mdtype
        u4_ptr[0] = 0
        return 1

    def read_element(self, int copy=True):
        cdef cnp.uint32_t mdtype, byte_count
        cdef char *data_ptr
        cdef size_t el_count
        cdef cnp.dtype dt
        cdef object data
        data = self.read_element_string(
            &mdtype, &byte_count, <void **>&data_ptr)
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
                            buffer=data)
        return el
            
    cdef object read_element_string(self,
                                    cnp.uint32_t *mdtype_ptr,
                                    cnp.uint32_t *byte_count_ptr,
                                    void **pp,
                                    int copy=True):
        cdef cnp.uint32_t mdtype, byte_count
        cdef char tag_data[4]
        cdef object data
        cdef int mod8
        cdef int tag_res = self.cread_tag(mdtype_ptr,
                                          byte_count_ptr,
                                          tag_data)
        mdtype = mdtype_ptr[0]
        byte_count = byte_count_ptr[0]
        if tag_res == 1: # full format
            if mdtype == miMATRIX:
                raise TypeError('Not expecting matrix here')
            data = self.cstream.read_string(
                byte_count,
                pp,
                copy)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
        else: # SDE format, make safer home for data
            data = streams.pyalloc_v(byte_count, <void **>pp)
            memcpy(pp[0], tag_data, byte_count)
        return data

    cdef int read_dims_into(self, cnp.int32_t *dims) except -1:
        ''' Read array dimensions into memory '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef int n_dims
        cdef int mod8
        cdef int i, res
        res = self.cread_tag(&mdtype, &byte_count, <char *>dims)
        if mdtype != miINT32:
            raise TypeError('Expecting miINT32 as data type')
            return -1
        if res == 1: # full format
            res = self.cstream.read_into(<void *>dims, byte_count)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
        n_dims = byte_count / 4
        if self.is_swapped:
            for i in range(n_dims):
                dims[i] = byteswap_u4(dims[i])
        return n_dims
    
    cdef cnp.uint32_t read_name_into(self, char *name_ptr) except -1:
        ''' Read matrix name into provided memory '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef int name_len
        cdef int mod8
        cdef int res
        res = self.cread_tag(&mdtype, &byte_count, name_ptr)
        if mdtype != miINT8:
            raise TypeError('Expecting miINT8 as data type')
            return -1
        if res == 1: # full format
            res = self.cstream.read_into(<void *>name_ptr, byte_count)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
        return byte_count

    def read_full_tag(self):
        cdef cnp.uint32_t mdtype, byte_count
        self.cread_full_tag(&mdtype, &byte_count)
        return mdtype, byte_count

    cdef void cread_full_tag(self,
                        cnp.uint32_t* mdtype,
                        cnp.uint32_t* byte_count):
        cdef cnp.uint32_t u4s[2]
        self.cstream.read_into(<void *>u4s, 8)
        if self.is_swapped:
            mdtype[0] = byteswap_u4(u4s[0])
            byte_count[0] = byteswap_u4(u4s[1])
        else:
            mdtype[0] = u4s[0]
            byte_count[0] = u4s[1]

    cpdef VarHeader5 read_header(self):
        ''' Return matrix header for current stream position

        Returns matrix headers at top level and sub levels
        '''
        cdef:
            cdef cnp.uint32_t u4s[2]
            cnp.uint32_t af_mdtype, af_byte_count
            cnp.uint32_t flags_class, nzmax
            cnp.uint16_t mc
            int ret, i
            VarHeader5 header
        # Read and discard mdtype and byte_count
        self.cstream.read_into(<void *>u4s, 8)
        # get array flags and nzmax
        self.cstream.read_into(<void *>u4s, 8)
        if self.is_swapped:
            flags_class = byteswap_u4(u4s[0])
            nzmax = byteswap_u4(u4s[1])
        else:
            flags_class = u4s[0]
            nzmax = u4s[1]
        header = VarHeader5()
        mc = flags_class & 0xFF
        # Might want to detect unreadable mclasses (such as functions)
        # here and divert to something that preserves the binary data.
        # For now we bail later in read_array_mdtype, via
        # read_mi_matrix, giving a MatReadError
        header.mclass = mc
        header.is_logical = flags_class >> 9 & 1
        header.is_global = flags_class >> 10 & 1
        header.is_complex = flags_class >> 11 & 1
        header.nzmax = nzmax
        header.n_dims = self.read_dims_into(header.dims_ptr)
        # convert dims to list
        header.dims = []
        for i in range(header.n_dims):
            header.dims.append(header.dims_ptr[i])
        header.name_len = self.read_name_into(header.name_ptr)
        return header

    def array_from_header(self, VarHeader5 header):
        ''' Only called at the top level '''
        arr, mat_dtype = self.read_array_mdtype(header)
        if header.mclass == mxSPARSE_CLASS:
            # no current processing makes sense for sparse
            return arr
        return process_element(arr, self.read_opts, mat_dtype)

    cdef cnp.ndarray read_mi_matrix(self):
        ''' Read header with matrix at sub-levels '''
        cdef:
            VarHeader5 header
            cnp.uint32_t mdtype, byte_count
            cnp.ndarray arr
        # read full tag
        self.cread_full_tag(&mdtype, &byte_count)
        if mdtype != miMATRIX:
            raise TypeError('Expecting matrix here')
        if byte_count == 0: # empty matrix
            mat_dtype = 'f8'
            arr = np.array([[]])
        else:
            header = self.read_header()
            arr, mat_dtype = self.read_array_mdtype(header)
        if header.mclass == mxSPARSE_CLASS:
            # no current processing makes sense for sparse
            return arr
        return process_element(arr, self.read_opts, mat_dtype)

    cpdef read_array_mdtype(self, VarHeader5 header):
        cdef:
            int mc
        mc = header.mclass
        mat_dtype = None
        if (mc == mxDOUBLE_CLASS
            or mc == mxSINGLE_CLASS
            or mc == mxINT8_CLASS
            or mc == mxUINT8_CLASS
            or mc == mxINT16_CLASS
            or mc == mxUINT16_CLASS
            or mc == mxINT32_CLASS
            or mc == mxUINT32_CLASS
            or mc == mxINT64_CLASS
            or mc == mxUINT64_CLASS): # numeric matrix
            if header.is_logical:
                mat_dtype = np.dtype('bool')
            else:
                mat_dtype = self.class_dtypes[mc]
            arr = self._read_numeric(header)
        elif mc == mxSPARSE_CLASS:
            arr = self._read_sparse(header)
        elif mc == mxCHAR_CLASS:
            arr = self._read_char(header)
        elif mc == mxCELL_CLASS:
            arr = self._read_cells(header)
        elif mc == mxSTRUCT_CLASS:
            arr = self._read_struct(header)
        elif mc == mxOBJECT_CLASS: # like structs, but with classname
            classname = self.read_element().tostring()
            arr = self._read_struct(header)
            arr = mio5p.MatlabObject(arr, classname)
        elif mc == mxFUNCTION_CLASS:
            raise miob.MatReadError('Cannot read matlab functions')
        return arr, mat_dtype

    cdef cnp.ndarray _read_numeric(self, VarHeader5 header):
        cdef:
            cnp.ndarray res, res_j
        if header.is_complex:
            # avoid array copy to save memory
            res = self.read_element(copy=False)
            res_j = self.read_element(copy=False)
            res = res + (res_j * 1j)
        else:
            res = self.read_element()
        return np.ndarray(shape=header.dims,
                          dtype=res.dtype,
                          buffer=res,
                          order='F')

    cdef object _read_sparse(self, VarHeader5 header):
        cdef cnp.ndarray rowind, indptr, data, data_j
        cdef size_t M, N, nnz
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
        return scipy.sparse.csc_matrix(
            (data,rowind,indptr),
            shape=(M,N))
                
    cdef cnp.ndarray _read_char(self, VarHeader5 header):
        res = self.read_element()
        # Convert non-string types to unicode
        if isinstance(res, np.ndarray):
            if res.dtype.type == np.uint16:
                codec = self.uint16_codec
                if self.codecs['uint16_len'] == 1:
                    res = res.astype(np.uint8)
            elif res.dtype.type in (np.uint8, np.int8):
                codec = 'ascii'
            else:
                raise TypeError, 'Did not expect type %s' % res.dtype
            res = res.tostring().decode(codec)
        return np.ndarray(shape=header.dims,
                          dtype=np.dtype('U1'),
                          buffer=np.array(res),
                          order='F').copy()
                             
    cdef cnp.ndarray _read_cells(self, VarHeader5 header):
        cdef:
            size_t length, i
            cnp.ndarray[object, ndim=1] result
        # Account for fortran indexing of cells
        tupdims = tuple(header.dims[::-1])
        length = cmio_utils.cproduct(tupdims)
        result = np.empty(length, dtype=object)
        for i in range(length):
            result[i] = self.read_mi_matrix()
        return result.reshape(tupdims).T
        
    cdef cnp.ndarray _read_struct(self, VarHeader5 header):
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
        length = cmio_utils.cproduct(tupdims)
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
                    rec_res[i][field_name] = self.read_mi_matrix()
            return rec_res.reshape(tupdims).T
        # Backward compatibility with previous format
        obj_template = mio5p.mat_struct()
        obj_template._fieldnames = field_names
        result = np.empty(length, dtype=object)
        for i in range(length):
            item = pycopy(obj_template)
            for name in field_names:
                item.__dict__[name] = self.read_mi_matrix()
            result[i] = item
        return result.reshape(tupdims).T


