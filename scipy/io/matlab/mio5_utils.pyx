''' Cython mio5 utility routines (-*- python -*- like)

'''

import sys

from copy import copy as pycopy

import numpy as np
cimport numpy as cnp

# Constant from numpy - max number of array dimensions
DEF _MAT_MAXDIMS = 32

cimport streams
import scipy.io.matlab.miobase as miob
from scipy.io.matlab.mio_utils import FileReadOpts, process_element
cimport mio_utils as cmio_utils
import scipy.io.matlab.mio5_params as mio5p
import scipy.sparse

from python_string cimport PyString_Size, PyString_FromString, \
    PyString_FromStringAndSize


cdef:
    int miINT8 = mio5p.miINT8
    int miUINT8 = mio5p.miUINT8
    int miUINT16 = mio5p.miUINT16
    int miDOUBLE = mio5p.miDOUBLE
    int miINT32 = mio5p.miINT32
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
    cdef readonly object name
    cdef int mclass
    cdef object dims
    cdef cnp.int32_t dims_ptr[_MAT_MAXDIMS]
    cdef int n_dims
    cdef int is_complex
    cdef int is_logical
    cdef public int is_global
    cdef size_t nzmax


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

    cdef object read_element(self,
                             cnp.uint32_t *mdtype_ptr,
                             cnp.uint32_t *byte_count_ptr,
                             void **pp,
                             int copy=True):
        ''' Read data element into string buffer, return buffer

        The element is the atom of the matlab file format.

        We use a python string as a container for the read memory.  If
        `copy` is True, this string is writable (from python),
        otherwise, it may not be. 

        See ``read_element_into`` for routine to read element into a
        pre-allocated block of memory.  
        '''
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
            data = self.cstream.read_string(
                byte_count,
                pp,
                copy)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
        else: # SDE format, make safer home for data
            data = PyString_FromStringAndSize(tag_data, byte_count)
        return data

    cdef void read_element_into(self,
                                cnp.uint32_t *mdtype_ptr,
                                cnp.uint32_t *byte_count_ptr,
                                void *ptr):
        ''' Read element into pre-allocated memory in `ptr`

        Compare ``read_element``.
        '''
        cdef:
           int mod8
        cdef int res = self.cread_tag(
            mdtype_ptr,
            byte_count_ptr,
            <char *>ptr)
        cdef cnp.uint32_t byte_count = byte_count_ptr[0]
        if res == 1: # full format
            res = self.cstream.read_into(ptr, byte_count)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
    
    cpdef inline cnp.ndarray read_numeric(self, int copy=True):
        ''' Read numeric data element into ndarray

        Reads element, then casts to ndarray. 

        The type of the array is given by the ``mdtype`` returned via
        ``read_element``. 
        '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef char *data_ptr
        cdef size_t el_count
        cdef cnp.dtype dt
        cdef object data = self.read_element(
            &mdtype, &byte_count, <void **>&data_ptr, copy)
        dt = self.dtypes[mdtype]
        el_count = byte_count // dt.itemsize
        return np.ndarray(shape=(el_count,),
                        dtype=dt,
                        buffer=data)
            
    cdef inline object read_int8_string(self):
        ''' Read, return int8 type string

        int8 type strings used for variable names, class names of
        objects, and field names of structs and objects.

        Specializes ``read_element``
        '''
        cdef:
            cnp.uint32_t mdtype, byte_count
            void *ptr
            object data
        data = self.read_element(&mdtype, &byte_count, &ptr)
        if mdtype != miINT8:
            raise TypeError('Expecting miINT8 as data type')
        return data

    cdef int read_into_int32s(self, cnp.int32_t *int32p) except -1:
        ''' Read int32 values into pre-allocated memory

        Byteswap as necessary.  Specializes ``read_element_into``

        Parameters
        ----------
        int32p : int32 pointer

        Returns
        -------
        n_ints : int
           Number of integers read
        '''
        cdef:
            cnp.uint32_t mdtype, byte_count
            int i
        self.read_element_into(&mdtype, &byte_count, <void *>int32p)
        if mdtype != miINT32:
            raise TypeError('Expecting miINT32 as data type')
            return -1
        cdef int n_ints = byte_count // 4
        if self.is_swapped:
            for i in range(n_ints):
                int32p[i] = byteswap_u4(int32p[i])
        return n_ints

    def read_full_tag(self):
        ''' Python method for reading full u4, u4 tag from stream

        Returns
        -------
        mdtype : int32
           matlab data type code
        byte_count : int32
           number of data bytes following

        Notes
        -----
        Assumes tag is in fact full, that is, is not a small data
        element.  This means it can skip some checks and makes it
        slightly faster than ``read_tag``
        '''
        cdef cnp.uint32_t mdtype, byte_count
        self.cread_full_tag(&mdtype, &byte_count)
        return mdtype, byte_count

    cdef void cread_full_tag(self,
                        cnp.uint32_t* mdtype,
                        cnp.uint32_t* byte_count):
        ''' C method for reading full u4, u4 tag from stream'''
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
            cnp.uint32_t mdtype, byte_count
            cnp.uint32_t flags_class, nzmax
            cnp.uint16_t mc
            int ret, i
            void *ptr
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
        header.n_dims = self.read_into_int32s(header.dims_ptr)
        # convert dims to list
        header.dims = []
        for i in range(header.n_dims):
            header.dims.append(header.dims_ptr[i])
        header.name = self.read_int8_string()
        return header

    cdef inline size_t size_from_header(self, VarHeader5 header):
        ''' Supporting routine for calculating array sizes from header

        Probably unnecessary optimization that uses integers stored in
        header rather than ``header.dims`` that is a python list.

        Parameters
        ----------
        header : VarHeader5
           array header

        Returns
        -------
        size : size_t
           size of array (product of dims)
        '''
        # calculate number of items in array from dims product
        cdef size_t size = 1
        for i in range(header.n_dims):
            size *= header.dims_ptr[i]
        return size

    def array_from_header(self, VarHeader5 header):
        '''Read array from stream, given array header

        Apply standard processing of array given ``read_opts`` options
        set in ``self``.

        Only called at the top level - that is, at when we start
        reading of each variable in the mat file. 
        '''
        arr, mat_dtype = self.read_array_mdtype(header)
        if header.mclass == mxSPARSE_CLASS:
            # no current processing makes sense for sparse
            return arr
        return process_element(arr, self.read_opts, mat_dtype)

    cdef cnp.ndarray read_mi_matrix(self):
        ''' Read header with matrix at sub-levels

        Combines ``read_header`` and functionality of
        ``array_from_header``.  Applies standard processing of array
        given ``read_opts`` options set in self. 
        '''
        cdef:
            VarHeader5 header
            cnp.uint32_t mdtype, byte_count
            object arr, mat_dtype
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
        ''' Read array of any type, return array and mat_dtype

        Parameters
        ----------
        header : VarHeader5
           array header object

        Returns
        -------
        arr : array or sparse array
           read array
        mat_dtype : None or dtype
           dtype of array as it would be read into matlab workspace, or
           None if it does not make sense that this would differ from
           the dtype of `arr`
        '''
        cdef:
            object arr
        cdef int mc = header.mclass
        cdef object mat_dtype = None
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
            arr = self.read_real_complex(header)
        elif mc == mxSPARSE_CLASS:
            arr = self.read_sparse(header)
        elif mc == mxCHAR_CLASS:
            arr = self.read_char(header)
        elif mc == mxCELL_CLASS:
            arr = self.read_cells(header)
        elif mc == mxSTRUCT_CLASS:
            arr = self.read_struct(header)
        elif mc == mxOBJECT_CLASS: # like structs, but with classname
            classname = self.read_int8_string()
            arr = self.read_struct(header)
            arr = mio5p.MatlabObject(arr, classname)
        elif mc == mxFUNCTION_CLASS:
            raise miob.MatReadError('Cannot read matlab functions')
        return arr, mat_dtype

    cpdef cnp.ndarray read_real_complex(self, VarHeader5 header):
        ''' Read real / complex matrices from stream '''
        cdef:
            cnp.ndarray res, res_j
        if header.is_complex:
            # avoid array copy to save memory
            res = self.read_numeric(False)
            res_j = self.read_numeric(False)
            res = res + (res_j * 1j)
        else:
            res = self.read_numeric()
        return res.reshape(header.dims[::-1]).T

    cdef object read_sparse(self, VarHeader5 header):
        ''' Read sparse matrices from stream '''
        cdef cnp.ndarray rowind, indptr, data, data_j
        cdef size_t M, N, nnz
        rowind = self.read_numeric()
        indptr = self.read_numeric()
        if header.is_complex:
            # avoid array copy to save memory
            data   = self.read_numeric(False)
            data_j = self.read_numeric(False)
            data = data + (data_j * 1j)
        else:
            data = self.read_numeric()
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
                
    cpdef cnp.ndarray read_char(self, VarHeader5 header):
        ''' Read char matrices from stream as arrays

        Matrices of char are likely to be converted to matrices of
        string by later processing in ``process_element``
        '''
        cdef:
            cnp.uint32_t mdtype, byte_count
            char *data_ptr
            size_t el_count
            cnp.dtype dt
            object data, res, codec
            cnp.ndarray arr
        cdef size_t length = self.size_from_header(header)
        data = self.read_element(
            &mdtype, &byte_count, <void **>&data_ptr, True)
        # Character data can be of apparently numerical types,
        # specifically np.uint8, np.int8, np.uint16.  np.unit16 can have
        # a length 1 type encoding, like ascii, or length 2 type encoding
        if mdtype == miUINT16:
            codec = self.uint16_codec
            if self.codecs['uint16_len'] == 1: # need LSBs only
                arr = np.ndarray(shape=(length,),
                                  dtype=self.dtypes[mdtype],
                                  buffer=data)
                data = arr.astype(np.uint8).tostring()
        elif mdtype == miINT8 or mdtype == miUINT8:
            codec = 'ascii'
        elif mdtype in self.codecs: # encoded char data
            codec = self.codecs[mdtype]
            if not codec:
                raise TypeError('Do not support encoding %d' % mdtype)
        else:
            raise ValueError('Type %d does not appear to be char type'
                             % mdtype)
        uc_str = data.decode(codec)
        # cast to array to deal with 2, 4 byte width characters
        arr = np.array(uc_str)
        if self.little_endian:
            dtc = '<U1'
        else:
            dtc = '>U1'
        return np.ndarray(shape=header.dims,
                          dtype=dtc,
                          buffer=arr,
                          order='F')
                             
    cpdef cnp.ndarray read_cells(self, VarHeader5 header):
        ''' Read cell array from stream '''
        cdef:
            size_t i
            cnp.ndarray[object, ndim=1] result
        # Account for fortran indexing of cells
        tupdims = tuple(header.dims[::-1])
        cdef size_t length = self.size_from_header(header)
        result = np.empty(length, dtype=object)
        for i in range(length):
            result[i] = self.read_mi_matrix()
        return result.reshape(tupdims).T
        
    cpdef cnp.ndarray read_struct(self, VarHeader5 header):
        ''' Read struct or object array from stream

        Objects are just structs with an extra field *classname*,
        defined before (this here) struct format structure
        '''
        cdef:
            cnp.int32_t namelength
            int i, n_names
            cnp.ndarray rec_res
            cnp.ndarray[object, ndim=1] result
            object name, field_names, dt, tupdims
        # Read field names into list
        cdef int res = self.read_into_int32s(&namelength)
        if res != 1:
            raise ValueError('Only one value for namelength')
        cdef object names = self.read_int8_string()
        field_names = []
        n_names = PyString_Size(names) // namelength
        cdef char *n_ptr = names
        for i in range(n_names):
            name = PyString_FromString(n_ptr)
            field_names.append(name)
            n_ptr += namelength
        # Prepare struct array
        tupdims = tuple(header.dims[::-1])
        cdef size_t length = self.size_from_header(header)
        if self.struct_as_record: # to record arrays
            if not n_names:
                # If there are no field names, there is no dtype
                # representation we can use, falling back to empty
                # object
                return np.empty(tupdims, dtype=object).T
            dt = [(field_name, object) for field_name in field_names]
            rec_res = np.empty(length, dtype=dt)
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


