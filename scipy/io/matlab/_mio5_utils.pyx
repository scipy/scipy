''' Cython mio5 utility routines (-*- python -*- like)

'''

# Programmer's notes
# ------------------
# Routines here have been reasonably optimized.

# The char matrix reading is not very fast, but it's not usually a
# bottleneck. See comments in ``read_char`` for possible ways to go if you
# want to optimize.

import sys

from copy import copy as pycopy

cimport cython

from libc.stdlib cimport calloc, free
from libc.string cimport strcmp

from cpython cimport Py_INCREF
from cpython cimport PyObject

cdef extern from "Python.h":
    unicode PyUnicode_FromString(const char *u)
    ctypedef struct PyTypeObject:
        pass

from cpython cimport PyBytes_Size

import numpy as np
cimport numpy as cnp

cdef extern from "numpy/arrayobject.h":
    PyTypeObject PyArray_Type
    cnp.ndarray PyArray_NewFromDescr(PyTypeObject *subtype,
                                     cnp.dtype newdtype,
                                     int nd,
                                     cnp.npy_intp* dims,
                                     cnp.npy_intp* strides,
                                     void* data,
                                     int flags,
                                     object parent)

cdef extern from "numpy_rephrasing.h":
    void PyArray_Set_BASE(cnp.ndarray arr, object obj)

# Numpy must be initialized before any code using the numpy C-API
# directly
cnp.import_array()

# Constant from numpy - max number of array dimensions
DEF _MAT_MAXDIMS = 32
# max number of integer indices of matlab data types (miINT8 etc)
DEF _N_MIS = 20
# max number of integer indices of matlab class types (mxINT8_CLASS etc)
DEF _N_MXS = 20

from . cimport _streams
from scipy.io.matlab._mio_utils import squeeze_element, chars_to_strings
import scipy.io.matlab._mio5_params as mio5p
from scipy.sparse import csc_array


cdef enum:
    miINT8 = 1
    miUINT8 = 2
    miINT16 = 3
    miUINT16 = 4
    miINT32 = 5
    miUINT32 = 6
    miSINGLE = 7
    miDOUBLE = 9
    miINT64 = 12
    miUINT64 = 13
    miMATRIX = 14
    miCOMPRESSED = 15
    miUTF8 = 16
    miUTF16 = 17
    miUTF32 = 18

cdef enum: # see comments in mio5_params
    mxCELL_CLASS = 1
    mxSTRUCT_CLASS = 2
    mxOBJECT_CLASS = 3
    mxCHAR_CLASS = 4
    mxSPARSE_CLASS = 5
    mxDOUBLE_CLASS = 6
    mxSINGLE_CLASS = 7
    mxINT8_CLASS = 8
    mxUINT8_CLASS = 9
    mxINT16_CLASS = 10
    mxUINT16_CLASS = 11
    mxINT32_CLASS = 12
    mxUINT32_CLASS = 13
    mxINT64_CLASS = 14
    mxUINT64_CLASS = 15
    mxFUNCTION_CLASS = 16
    mxOPAQUE_CLASS = 17 # This appears to be a function workspace
    mxOBJECT_CLASS_FROM_MATRIX_H = 18

cdef bint sys_is_le = sys.byteorder == 'little'
swapped_code = '>' if sys_is_le else '<'

cdef cnp.dtype OPAQUE_DTYPE = mio5p.OPAQUE_DTYPE
cdef cnp.dtype BOOL_DTYPE = np.dtype(np.bool_)


cpdef cnp.uint32_t byteswap_u4(cnp.uint32_t u4) noexcept:
    return ((u4 << 24) |
           ((u4 << 8) & 0xff0000U) |
           ((u4 >> 8 & 0xff00u)) |
           (u4 >> 24))


cdef class VarHeader5:
    cdef readonly object name 
    cdef readonly int mclass
    cdef readonly object dims
    cdef cnp.int32_t dims_ptr[_MAT_MAXDIMS]
    cdef int n_dims
    cdef int check_stream_limit
    cdef int is_complex
    cdef readonly int is_logical
    cdef public int is_global
    cdef readonly size_t nzmax

    def set_dims(self, dims):
        """ Allow setting of dimensions from python

        This is for constructing headers for tests
        """
        self.dims = dims
        self.n_dims = len(dims)
        for i, dim in enumerate(dims):
            self.dims_ptr[i] = <cnp.int32_t>int(dim)


cdef class VarReader5:
    """Initialize from file reader object

    preader needs the following fields defined:

    * mat_stream (file-like)
    * byte_order (str)
    * uint16_codec (str)
    * struct_as_record (bool)
    * chars_as_strings (bool)
    * mat_dtype (bool)
    * squeeze_me (bool)
    """

    cdef public int is_swapped, little_endian
    cdef int struct_as_record
    cdef object codecs, uint16_codec
    # c-optimized version of reading stream
    cdef _streams.GenericStream cstream
    # pointers to stuff in preader.dtypes
    cdef PyObject* dtypes[_N_MIS]
    # pointers to stuff in preader.class_dtypes
    cdef PyObject* class_dtypes[_N_MXS]
    # element processing options
    cdef:
        int mat_dtype
        int squeeze_me
        int chars_as_strings

    def __cinit__(self, preader):
        byte_order = preader.byte_order
        self.is_swapped = byte_order == swapped_code
        if self.is_swapped:
            self.little_endian = not sys_is_le
        else:
            self.little_endian = sys_is_le
        # option affecting reading of matlab struct arrays
        self.struct_as_record = preader.struct_as_record
        # store codecs for text matrix reading
        self.codecs = mio5p.MDTYPES[byte_order]['codecs'].copy()
        self.uint16_codec = preader.uint16_codec
        uint16_codec = self.uint16_codec
        # Set length of miUINT16 char encoding
        self.codecs['uint16_len'] = len("  ".encode(uint16_codec)) \
                - len(" ".encode(uint16_codec))
        self.codecs['uint16_codec'] = uint16_codec
        # set c-optimized stream object from python file-like object
        self.cstream = _streams.make_stream(preader.mat_stream)
        # options for element processing
        self.mat_dtype = preader.mat_dtype
        self.chars_as_strings = preader.chars_as_strings
        self.squeeze_me = preader.squeeze_me
        # copy refs to dtypes into object pointer array. We only need the
        # integer-keyed dtypes
        for key, dt in mio5p.MDTYPES[byte_order]['dtypes'].items():
            if isinstance(key, str):
                continue
            self.dtypes[key] = <PyObject*>dt
        # copy refs to class_dtypes into object pointer array
        for key, dt in mio5p.MDTYPES[byte_order]['classes'].items():
            if isinstance(key, str):
                continue
            self.class_dtypes[key] = <PyObject*>dt

    def set_stream(self, fobj):
        ''' Set stream of best type from file-like `fobj`

        Called from Python when initiating a variable read
        '''
        self.cstream = _streams.make_stream(fobj)

    def read_tag(self):
        ''' Read tag mdtype and byte_count

        Does necessary swapping and takes account of SDE formats.

        See also ``read_full_tag`` method.

        Returns
        -------
        mdtype : int
           matlab data type code
        byte_count : int
           number of bytes following that comprise the data
        tag_data : None or str
           Any data from the tag itself.  This is None for a full tag,
           and string length `byte_count` if this is a small data
           element.
        '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef char tag_ptr[4]
        cdef int tag_res
        cdef object tag_data = None
        tag_res = self.cread_tag(&mdtype, &byte_count, tag_ptr)
        if tag_res == 2: # sde format
            tag_data = tag_ptr[:byte_count]
        return (mdtype, byte_count, tag_data)

    cdef int cread_tag(self,
                     cnp.uint32_t *mdtype_ptr,
                     cnp.uint32_t *byte_count_ptr,
                     char *data_ptr) except -1:
        ''' Read tag mdtype and byte_count

        Does necessary swapping and takes account of SDE formats

        Data may be returned in data_ptr, if this was an SDE.
        data_ptr must point to a buffer of size >= 4 bytes.

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

        Parameters
        ----------
        mdtype_ptr : uint32_t*
           pointer to uint32_t value to which we write the mdtype value
        byte_count_ptr : uint32_t*
           pointer to uint32_t value to which we write the byte count
        pp : void**
           pointer to void*. pp[0] will be set to point to the start of
           the returned string memory
        copy : int
           If not 0, do any copies required to allow memory to be freely
           altered without interfering with other objects.  Otherwise
           return string that should not be written to, therefore saving
           unnecessary copies

        Return
        ------
        data : str
           Python string object containing read data

        Notes
        -----
        See ``read_element_into`` for routine to read element into a
        pre-allocated block of memory.
        '''
        cdef cnp.uint32_t byte_count
        cdef char tag_data[4]
        cdef object data
        cdef int mod8
        cdef int tag_res = self.cread_tag(mdtype_ptr,
                                          byte_count_ptr,
                                          tag_data)
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
            data = tag_data[:byte_count]
            pp[0] = <char *>data
        return data

    cdef int read_element_into(self,
                               cnp.uint32_t *mdtype_ptr,
                               cnp.uint32_t *byte_count_ptr,
                               void *ptr,
                               cnp.uint32_t max_byte_count) except -1:
        ''' Read element into pre-allocated memory in `ptr`

        Parameters
        ----------
        mdtype_ptr : uint32_t*
           pointer to uint32_t value to which we write the mdtype value
        byte_count_ptr : uint32_t*
           pointer to uint32_t value to which we write the byte count
        ptr : void*
           memory buffer into which to read
        max_byte_count : uint32_t
           size of the buffer pointed to by ptr

        Returns
        -------
        void

        Notes
        -----
        Compare ``read_element``.
        '''
        cdef:
           int mod8
        if max_byte_count < 4:
            raise ValueError('Unexpected amount of data to read (malformed input file?)')
        cdef int res = self.cread_tag(
            mdtype_ptr,
            byte_count_ptr,
            <char *>ptr)
        cdef cnp.uint32_t byte_count = byte_count_ptr[0]
        if res == 1: # full format
            if byte_count > max_byte_count:
                raise ValueError('Unexpected amount of data to read (malformed input file?)')
            res = self.cstream.read_into(ptr, byte_count)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.cstream.seek(8 - mod8, 1)
        return 0

    cpdef cnp.ndarray read_numeric(self, int copy=True, size_t nnz=-1):
        ''' Read numeric data element into ndarray

        Reads element, then casts to ndarray.

        The type of the array is usually given by the ``mdtype`` returned via
        ``read_element``.  Sparse logical arrays are an exception, where the
        type of the array may be ``np.bool_`` even if the ``mdtype`` claims the
        data is of float64 type.

        Parameters
        ----------
        copy : bool, optional
            Whether to copy the array before returning.  If False, return array
            backed by bytes read from file.
        nnz : int, optional
            Number of non-zero values when reading numeric data from sparse
            matrices.  -1 if not reading sparse matrices, or to disable check
            for bytes data instead of declared data type (see Notes).

        Returns
        -------
        arr : array
            Numeric array

        Notes
        -----
        MATLAB apparently likes to store sparse logical matrix data as bytes
        instead of miDOUBLE (float64) data type, even though the data element
        still declares its type as miDOUBLE.  We can guess this has happened by
        looking for the length of the data compared to the expected number of
        elements, using the `nnz` input parameter.
        '''
        cdef cnp.uint32_t mdtype, byte_count
        cdef void *data_ptr
        cdef cnp.npy_intp el_count
        cdef cnp.ndarray el
        cdef object data = self.read_element(
            &mdtype, &byte_count, <void **>&data_ptr, copy)
        cdef cnp.dtype dt = <cnp.dtype>self.dtypes[mdtype]
        if dt.itemsize != 1 and nnz != -1 and byte_count == nnz:
            el_count = <cnp.npy_intp> nnz
            dt = BOOL_DTYPE
        else:
            el_count = byte_count // dt.itemsize
        cdef int flags = 0
        if copy:
            flags = cnp.NPY_ARRAY_WRITEABLE
        Py_INCREF(<object> dt)
        el = PyArray_NewFromDescr(&PyArray_Type,
                                   dt,
                                   1,
                                   &el_count,
                                   NULL,
                                   <void*>data_ptr,
                                   flags,
                                   <object>NULL)
        Py_INCREF(<object> data)
        PyArray_Set_BASE(el, data)
        return el

    cdef inline object read_int8_string(self):
        ''' Read, return int8 type string

        int8 type strings used for variable names, class names of
        objects, and field names of structs and objects.

        Specializes ``read_element``
        '''
        cdef:
            cnp.uint32_t mdtype, byte_count, i
            void* ptr
            unsigned char* byte_ptr
            object data
        data = self.read_element(&mdtype, &byte_count, &ptr)
        if mdtype == miUTF8:  # Some badly-formed .mat files have utf8 here
            byte_ptr = <unsigned char*> ptr
            for i in range(byte_count):
                if byte_ptr[i] > 127:
                    raise ValueError('Non ascii int8 string')
        elif mdtype != miINT8:
            raise TypeError('Expecting miINT8 as data type')
        return data

    cdef int read_into_int32s(self, cnp.int32_t *int32p, cnp.uint32_t max_byte_count) except -1:
        ''' Read int32 values into pre-allocated memory

        Byteswap as necessary.  Specializes ``read_element_into``

        Parameters
        ----------
        int32p : int32 pointer
        max_count : uint32_t

        Returns
        -------
        n_ints : int
           Number of integers read
        '''
        cdef:
            cnp.uint32_t mdtype, byte_count, n_ints
            int i, check_ints=0
        self.read_element_into(&mdtype, &byte_count, <void *>int32p, max_byte_count)
        if mdtype == miUINT32:
            check_ints = 1
        elif mdtype != miINT32:
            raise TypeError('Expecting miINT32 as data type')
        n_ints = byte_count // 4
        if self.is_swapped:
            for i in range(n_ints):
                int32p[i] = byteswap_u4(int32p[i])
        if check_ints:
            for i in range(n_ints):
                if int32p[i] < 0:
                    raise ValueError('Expecting miINT32, got miUINT32 with '
                                     'negative values')
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

    cdef int cread_full_tag(self,
                            cnp.uint32_t* mdtype,
                            cnp.uint32_t* byte_count) except -1:
        ''' C method for reading full u4, u4 tag from stream'''
        cdef cnp.uint32_t u4s[2]
        self.cstream.read_into(<void *>u4s, 8)
        if self.is_swapped:
            mdtype[0] = byteswap_u4(u4s[0])
            byte_count[0] = byteswap_u4(u4s[1])
        else:
            mdtype[0] = u4s[0]
            byte_count[0] = u4s[1]
        return 0

    cpdef VarHeader5 read_header(self, int check_stream_limit):
        ''' Return matrix header for current stream position

        Returns matrix headers at top level and sub levels

        Parameters
        ----------
        check_stream_limit : if True, then if the returned header
        is passed to array_from_header, it will be verified that
        the length of the uncompressed data is not overlong (which
        can indicate .mat file corruption)
        '''
        cdef:
            cdef cnp.uint32_t u4s[2]
            cnp.uint32_t flags_class, nzmax
            cnp.uint16_t mc
            int i
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
        header.mclass = mc
        header.check_stream_limit = check_stream_limit
        header.is_logical = flags_class >> 9 & 1
        header.is_global = flags_class >> 10 & 1
        header.is_complex = flags_class >> 11 & 1
        header.nzmax = nzmax
        # all miMATRIX types except the mxOPAQUE_CLASS have dims and a
        # name.
        if mc == mxOPAQUE_CLASS:
            # Only extract variable name for mxOPAQUE_CLASS
            # Remaining metadata is used in read_opaque()
            header.name = self.read_int8_string()
            return header
        header.n_dims = self.read_into_int32s(header.dims_ptr, sizeof(header.dims_ptr))
        if header.n_dims > _MAT_MAXDIMS:
            raise ValueError('Too many dimensions (%d) for numpy arrays'
                             % header.n_dims)
        # convert dims to list
        header.dims = [header.dims_ptr[i] for i in range(header.n_dims)]
        header.name = self.read_int8_string()
        return header

    cdef inline size_t size_from_header(self, VarHeader5 header) noexcept:
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
           size of array referenced by header (product of dims)
        '''
        # calculate number of items in array from dims product
        cdef size_t size = 1
        cdef int i
        for i in range(header.n_dims):
            size *= header.dims_ptr[i]
        return size

    cdef read_mi_matrix(self, int process=1, SubsystemReader5 subsystem=None, is_opaque=0):
        ''' Read header with matrix at sub-levels

        Combines ``read_header`` and functionality of
        ``array_from_header``.  Applies standard processing of array
        given options set in self.

        Parameters
        ----------
        process : int, optional
           If not zero, apply post-processing on returned array

        Returns
        -------
        arr : ndarray or sparse csc_array
        '''
        cdef:
            VarHeader5 header
            cnp.uint32_t mdtype, byte_count
        # read full tag
        self.cread_full_tag(&mdtype, &byte_count)
        if mdtype != miMATRIX:
            raise TypeError('Expecting matrix here')
        if byte_count == 0: # empty matrix
            if process and self.squeeze_me:
                return np.array([])
            else:
                return np.array([[]])
        header = self.read_header(False)
        return self.array_from_header(header, process, subsystem, is_opaque)

    cpdef array_from_header(self, VarHeader5 header, int process=1, SubsystemReader5 subsystem=None, is_opaque=0):
        ''' Read array of any class, given matrix `header`

        Parameters
        ----------
        header : VarHeader5
           array header object
        process : int, optional
           If not zero, apply post-processing on returned array

        Returns
        -------
        arr : array or sparse csc_array
           read array
        '''
        cdef:
            object arr
            cnp.dtype mat_dtype
        cdef int mc = header.mclass
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
            arr = self.read_real_complex(header, subsystem, is_opaque)
            if process and self.mat_dtype: # might need to recast
                if header.is_logical:
                    mat_dtype = BOOL_DTYPE
                else:
                    mat_dtype = <object>self.class_dtypes[mc]
                arr = arr.astype(mat_dtype)
        elif mc == mxSPARSE_CLASS:
            arr = self.read_sparse(header)
            # no current processing makes sense for sparse
            process = False
        elif mc == mxCHAR_CLASS:
            arr = self.read_char(header)
            if process and self.chars_as_strings:
                arr = chars_to_strings(arr)
        elif mc == mxCELL_CLASS:
            arr = self.read_cells(header, subsystem, is_opaque)
        elif mc == mxSTRUCT_CLASS:
            arr = self.read_struct(header, subsystem, is_opaque)
        elif mc == mxOBJECT_CLASS: # like structs, but with classname
            classname = self.read_int8_string().decode('latin1')
            arr = self.read_struct(header, subsystem, is_opaque)
            arr = mio5p.MatlabObject(arr, classname)
        elif mc == mxFUNCTION_CLASS: # just a matrix of struct type
            arr = self.read_mi_matrix(process, subsystem)
            # arr = mio5p.MatlabFunction(arr)
            # to make them more re-writeable - don't squeeze
            process = 0
        elif mc == mxOPAQUE_CLASS:
            arr = self.read_opaque(header, subsystem)
            # to make them more re-writeable - don't squeeze
            process = 0
        # ensure we have read checksum.
        read_ok = self.cstream.all_data_read()
        if header.check_stream_limit and not read_ok:
            raise ValueError('Did not fully consume compressed contents' +
                             ' of an miCOMPRESSED element. This can' +
                             ' indicate that the .mat file is corrupted.')
        if process and self.squeeze_me:
            return squeeze_element(arr)
        return arr

    def shape_from_header(self, VarHeader5 header):
        cdef int mc = header.mclass
        cdef tuple shape
        if mc == mxSPARSE_CLASS:
            shape = tuple(header.dims)
        elif mc == mxCHAR_CLASS:
            shape = tuple(header.dims)
            if self.chars_as_strings:
                shape = shape[:-1]
        else:
            shape = tuple(header.dims)
        if self.squeeze_me:
            shape = tuple([x for x in shape if x != 1])
        return shape

    cpdef int is_mcos_reference(self, arr):
        ''' Check if array is a reference to an object
        Object references are uint32 arrays with a specific structure.
        1. The first value is the magic value 0xDD000000
        2. The second value is ndims, followed by size of each dimension
        3. The next prod(dims) values are object IDs 
        4. The last value is always the class ID
        '''
        cdef:
            cnp.uint32_t n_dims
            cnp.uint32_t total_objs
        
        if arr.dtype != np.uint32:
            return 0
        if arr.size < 6:
            return 0
        if arr[0] != 0xdd000000:
            # First value of object references are zeros
            return 0
        n_dims = arr[1]
        if n_dims < 2:
            # Object is min 2D
            return 0
        total_objs = np.prod(arr[2 : 2 + n_dims])
        if total_objs <= 0:
            # At least 1 object should be present
            return 0
        if np.any(arr[2 + n_dims : 2 + n_dims + total_objs] <= 0):
            # All objects should have object_id > 0
            return 0
        if arr.size != 3 + n_dims + total_objs:
            return 0
        
        # print('Reference Found!')
        return 1

    cpdef cnp.ndarray read_real_complex(self, VarHeader5 header, SubsystemReader5 subsystem, int is_opaque=0):
        ''' Read real / complex matrices from stream '''
        cdef:
            cnp.ndarray res, res_j
            cdef cnp.uint32_t class_id, ndims
            cnp.ndarray object_ids, dims
        if header.is_complex:
            # avoid array copy to save memory
            res = self.read_numeric(False)
            res_j = self.read_numeric(False)
            # Use c8 for f4s and c16 for f8 input. Just ``res = res + res_j *
            # 1j`` upcasts to c16 regardless of input type.
            if res.itemsize == 4:
                res = res.astype('c8')
            else:
                res = res.astype('c16')
            res.imag = res_j
        else:
            res = self.read_numeric()
        
        if is_opaque and self.is_mcos_reference(res):
            # Read object array if reference
            ndims = res[1]
            dims = res[2 : 2 + ndims]
            total_objs = np.prod(dims)
            object_ids = res[2 + ndims : 2 + ndims + total_objs]
            class_id = res[-1]
            return subsystem.read_mcos_arrays(object_ids, class_id, dims)
        return res.reshape(header.dims[::-1]).T

    cdef object read_sparse(self, VarHeader5 header):
        ''' Read sparse matrices from stream '''
        cdef cnp.ndarray rowind, indptr, data, data_j
        cdef size_t M, N, nnz
        rowind = self.read_numeric()
        indptr = self.read_numeric()
        M, N = header.dims[0], header.dims[1]
        indptr = indptr[:N+1]
        nnz = indptr[-1]
        if header.is_complex:
            # avoid array copy to save memory
            data   = self.read_numeric(False)
            data_j = self.read_numeric(False)
            data = data + (data_j * 1j)
        elif header.is_logical:
            data = self.read_numeric(True, nnz)
        else:
            data = self.read_numeric()

        # From the matlab (TM) API documentation, last found here:
        # https://www.mathworks.com/help/pdf_doc/matlab/apiext.pdf
        # rowind are simply the row indices for all the (nnz) non-zero
        # entries in the sparse array.  rowind has nzmax entries, so
        # may well have more entries than nnz, the actual number of
        # non-zero entries, but rowind[nnz:] can be discarded and
        # should be 0. indptr has length (number of columns + 1), and
        # is such that, if D = diff(colind), D[j] gives the number of
        # non-zero entries in column j. Because rowind values are
        # stored in column order, this gives the column corresponding
        # to each rowind

        return csc_array((data[:nnz], rowind[:nnz], indptr), shape=(M, N))

    cpdef cnp.ndarray read_char(self, VarHeader5 header):
        ''' Read char matrices from stream as arrays

        Matrices of char are likely to be converted to matrices of
        string by later processing in ``array_from_header``
        '''

        # Notes to friendly fellow-optimizer
        #
        # This routine is not much optimized.  If I was going to do it,
        # I'd store the codecs as an object pointer array, as for the
        # .dtypes, I might use python_string.PyBytes_Decode for decoding,
        # I'd do something with pointers to pull the LSB out of the uint16
        # dtype, without using an intermediate array, I guess I'd consider
        # using the numpy C-API for array creation. I'd try and work out
        # how to deal with UCS-2 and UCS-4 builds of python, and how numpy
        # deals with unicode strings passed as memory,
        #
        # My own unicode introduction here:
        # http://matthew-brett.github.com/pydagogue/python_unicode.html

        cdef:
            cnp.uint32_t mdtype, byte_count
            char *data_ptr
            object data, codec
            cnp.ndarray arr
            cnp.dtype dt
        cdef size_t length = self.size_from_header(header)
        data = self.read_element(
            &mdtype, &byte_count, <void **>&data_ptr, True)
        # There are mat files in the wild that have 0 byte count strings, but
        # maybe with non-zero length.
        if byte_count == 0:
            arr = np.array(' ' * length, dtype='U')
            return np.ndarray(shape=header.dims,
                              dtype='U1',
                              buffer=arr,
                              order='F')
        # Character data can be of apparently numerical types,
        # specifically np.uint8, np.int8, np.uint16.  np.unit16 can have
        # a length 1 type encoding, like ascii, or length 2 type
        # encoding
        dt = <cnp.dtype>self.dtypes[mdtype]
        if mdtype == miUINT16:
            codec = self.uint16_codec
            if self.codecs['uint16_len'] == 1: # need LSBs only
                arr = np.ndarray(shape=(length,),
                                  dtype=dt,
                                  buffer=data)
                data = arr.astype(np.uint8).tobytes()
        elif mdtype == miINT8 or mdtype == miUINT8:
            codec = 'ascii'
        elif mdtype in self.codecs: # encoded char data
            codec = self.codecs[mdtype]
            if not codec:
                raise TypeError('Do not support encoding %d' % mdtype)
        else:
            raise ValueError('Type %d does not appear to be char type'
                             % mdtype)
        uc_str = data.decode(codec, 'replace')
        # cast to array to deal with 2, 4 byte width characters
        arr = np.array(uc_str, dtype='U')
        # could take this to numpy C-API level, but probably not worth
        # it
        return np.ndarray(shape=header.dims,
                          dtype='U1',
                          buffer=arr,
                          order='F')

    cpdef cnp.ndarray read_cells(self, VarHeader5 header, SubsystemReader5 subsystem, int is_opaque):
        ''' Read cell array from stream '''
        cdef:
            size_t i
            cnp.ndarray[object, ndim=1] result
        # Account for fortran indexing of cells
        tupdims = tuple(header.dims[::-1])
        cdef size_t length = self.size_from_header(header)
        result = np.empty(length, dtype=object)
        for i in range(length):
            result[i] = self.read_mi_matrix(process=1, subsystem=subsystem, is_opaque=is_opaque)
        return result.reshape(tupdims).T

    def read_fieldnames(self):
        '''Read fieldnames for struct-like matrix.'''
        cdef int n_names
        return self.cread_fieldnames(&n_names)

    cdef inline object cread_fieldnames(self, int *n_names_ptr):
        cdef:
            cnp.int32_t namelength
            int i, n_names
            list field_names
            object name
        # Read field names into list
        cdef int res = self.read_into_int32s(&namelength, 4)
        if res != 1:
            raise ValueError('Only one value for namelength')
        cdef object names = self.read_int8_string()
        field_names = []
        n_names = PyBytes_Size(names) // namelength
        # Make n_duplicates and pointer arrays
        cdef:
            int *n_duplicates
        n_duplicates = <int *>calloc(n_names, sizeof(int))
        cdef:
            char *names_ptr = names
            char *n_ptr = names
            int j, dup_no
        for i in range(n_names):
            name = PyUnicode_FromString(n_ptr)
            # Check if this is a duplicate field, rename if so
            dup_no = 0
            for j in range(i):
                if strcmp(n_ptr, names_ptr + j * namelength) == 0: # the same
                    n_duplicates[j] += 1
                    dup_no = n_duplicates[j]
                    break
            if dup_no != 0:
                name = '_%d_%s' % (dup_no, name)
            field_names.append(name)
            n_ptr += namelength
        free(n_duplicates)
        n_names_ptr[0] = n_names
        return field_names

    cpdef cnp.ndarray read_struct(self, VarHeader5 header, SubsystemReader5 subsystem, int is_opaque):
        ''' Read struct or object array from stream

        Objects are just structs with an extra field *classname*,
        defined before (this here) struct format structure
        '''
        cdef:
            int i, n_names
            cnp.ndarray[object, ndim=1] result
            object dt, tupdims
        # Read field names into list
        cdef object field_names = self.cread_fieldnames(&n_names)
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
                    rec_res[i][field_name] = self.read_mi_matrix(process=1, subsystem=subsystem, is_opaque=is_opaque)
            return rec_res.reshape(tupdims).T
        # Backward compatibility with previous format
        obj_template = mio5p.mat_struct()
        obj_template._fieldnames = field_names
        result = np.empty(length, dtype=object)
        for i in range(length):
            item = pycopy(obj_template)
            for name in field_names:
                item.__dict__[name] = self.read_mi_matrix(process=1, subsystem=subsystem, is_opaque=is_opaque)
            with cython.boundscheck(False):
                result[i] = item
        return result.reshape(tupdims).T

    cpdef object read_opaque(self, VarHeader5 hdr, SubsystemReader5 subsystem):
        ''' Read mxOPAQUE_CLASS arrays
        mxOPAQUE Class has the format:
        1. Array Flag
        2. Variable Name
        3. Type System (e.g. MCOS or java)
        4. Class Name (e.g. string)
        5. Object Reference: A metadata array that points to object contents in subsystem
        '''
        cdef:
            object type_system
            object class_name
            int is_opaque = 1 # Flag used to check for object references in mxOPAQUE_CLASS arrays
        
        type_system = PyUnicode_FromString(self.read_int8_string())
        class_name = PyUnicode_FromString(self.read_int8_string())
        # print(type_system, class_name)
        
        if type_system == "MCOS":
            if class_name == "FileWrapper__":
                is_opaque = 0 
                #! Needs to be fixed as subsystem can also contain reference within struct/cell fields
                #! Not checking object reference for FileWrapper__ class currently as it is not built yet
                
        res = self.read_mi_matrix(process=0, subsystem=subsystem, is_opaque=is_opaque)
        # If an object reference is detected during read_mi_matrix(), then the corresponding object extracted from subsystem is returned
        # The check is implemented within read_real_complex()
        return res

cdef class SubsystemReader5:
    ''' Read subsystem data'''
    cdef:
        object byte_order
        int raw_data
        cnp.ndarray mcos_metadata
        cnp.ndarray mcos_field_contents
        list mcos_names
        
    def __cinit__(self, byte_order, int raw_data=0):
        self.byte_order = '<u4' if byte_order == '<' else '>u4' # Subsystem MCOS metadata uses uint32 integers
        self.raw_data = raw_data
        self.mcos_metadata = None
        self.mcos_field_contents = None
        self.mcos_names = None
    
    def set_fields(self,  arr):
        '''Caches the field contents within subsystem
        For MCOS arrays, the field contents is a cell array
        Ignore first two cells as it is metadata/padding
        '''
        if "MCOS" in arr.dtype.names:
            # Assumes struct dimensions are not squeezed
            self.mcos_field_contents = arr["MCOS"][0,0][2:]
            self.mcos_metadata = arr["MCOS"][0,0][0,0][:,0]
            version = self.mcos_metadata[0]
            if version > 4 or version < 1:
                # As far as I am aware, the only difference between version 2 to 4 is the introduction of additional offsets for parsing other type of data
                # Can change logic here
                # Use version number to determine num_offsets
                raise TypeError("Only compatible with Version 4 MCOS objects")
            
            self.set_mcos_names()

    cdef set_mcos_names(self):
        '''Extracts list of all field names and class names from MCOS FileWrapper metadata
        1. Field names and class names are stored as a list of null terminated strings
        2, Second uint32 number gives the number of strings
        3. Names start immediately after offsets
        '''
        cdef:
            cnp.uint32_t i, start, end, count
            list names
    
        end = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=8)
        num_strings = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=4)
        names = []
        count = 0
        i = 40
        
        while i < end and count < num_strings:
            start = i
            while i < end and self.mcos_metadata[i] != 0:
                i += 1
            names.append(PyUnicode_FromString(self.mcos_metadata[start:i].tobytes()))
            i += 1
            count += 1
        
        self.mcos_names = names
    
    cdef cnp.ndarray get_object_dependencies(self, cnp.uint32_t object_id):
        ''' Get object dependencies for the current subsystem
        Format: (class_id, 0, 0, type1_id, type2_id, dep_id)
        Start of object deps region is determined by offset 3
        '''
        cdef:
            cnp.uint32_t byte_offset
            cnp.ndarray[cnp.uint32_t, ndim=1] object_ids
        
        byte_offset = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=16)
        byte_offset = byte_offset + object_id * 24

        object_ids = np.frombuffer(self.mcos_metadata, 
        dtype=self.byte_order, 
        count=6, 
        offset=byte_offset)
        
        return object_ids

    cdef object get_class_name(self, cnp.uint32_t class_id):
        ''' Get class name for the current subsystem
        Format: (handle_class_idx, class_idx, 0, 0)
        1. idx points to index of class name in mcos_names (indexed from 1)
        2. Start of class name region is determined by offset 1
        '''
        cdef:
            cnp.uint32_t byte_offset, handle_idx, class_idx
            object class_name
            object handle_name
        
        byte_offset = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=8)
        byte_offset = byte_offset + class_id * 16

        handle_idx, class_idx = np.frombuffer(
            self.mcos_metadata,
            dtype=self.byte_order,
            count=2,
            offset=byte_offset
        )

        class_name = self.mcos_names[class_idx - 1]
        if handle_idx > 0:
            handle_name = self.names[class_idx - 1]
            class_name = f"{handle_name}.{class_name}"

        return class_name

    cdef cnp.ndarray get_ids(self, cnp.uint32_t type_id, cnp.uint32_t byte_offset, cnp.uint32_t nbytes):
        ''' Get object ids for the current subsystem
        Helper method to get IDs for metadata regions in the format (nblocks, subblock1, subblock2, subblock3 ...)
        Input Parameters:
            1. type_id: ID indicating object number
            2. byte_offset: Offset to start of metadata region
            3. nbytes: Number of bytes for each subblock 
        
        Output:
            nump.ndarray containing all subblock data for the given type_id
        '''
        cdef cnp.uint32_t nblocks
        
        while type_id > 0:
            nblocks = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=byte_offset)
            byte_offset = byte_offset + 4 + nblocks * nbytes
            if ((nblocks * nbytes) + 4) % 8 != 0:
                byte_offset += 4 # Padding bytes
            type_id -= 1
        
        nblocks = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=byte_offset)
        byte_offset += 4
        content_ids = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=nblocks * (nbytes//4), offset=byte_offset)
        return content_ids.reshape(nblocks, nbytes//4)

    cdef get_handle_class_instance(self, cnp.uint32_t type2_id):
        ''' Searches for the object ID of the handle class instance
        Input:
            1. Type 2 ID: Type 2 ID of the handle class instance
        Output:
            1. Class ID: Class ID of the handle class instance
            2. Object ID: Object ID of the handle class instance
        
        '''
        cdef:
            cnp.uint32_t start,end, class_id, object_id
            cnp.ndarray object_deps
            int i = 0
            int nblocks

        start, end = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=2, offset=16)
        object_deps = np.frombuffer(self.mcos_metadata[start:end]).reshape(-1, 6)

        # Find block with type2_id
        nblocks = object_deps.shape[0]
        while i < nblocks:
            if object_deps[i, 4] == type2_id:
                class_id = object_deps[i, 0]
                object_id = <cnp.uint32_t> i
                return class_id, object_id
            i += 1

        raise ValueError("Could not find handle class instance for object")
        
    cdef cnp.ndarray extract_properties(self, cnp.ndarray[cnp.uint32_t, ndim=1] object_deps):
        '''Extract properties for a given object ID and its dependency IDs
        '''
        cdef:
            cnp.uint32_t byte_offset, obj_type_id
            object field_names
            object field_dtype
            cnp.ndarray field_ids
            cnp.ndarray handle_ids
            cnp.uint32_t class_id, object_id
            int i, j

        # Determine object type (1 or 2)
        if object_deps[3] == 0 and object_deps[4] != 0:
            obj_type_id = object_deps[4]
            byte_offset = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=20)
        elif object_deps[3] != 0 and object_deps[4] == 0:
            obj_type_id = object_deps[3]
            byte_offset = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=12)
        else:
            # Possible type 3 exists but has not been spotted in the wild yet
            raise ValueError("Could not determine object type. Got dependencies: {object_deps}")
        
        # Get Field Content metadata for type 1 or type 2 objects
        # Field content metadata format: (field_name_idx, field_type, field_value)
        # If field_type == 1, then field_value is the index of the field value in mcos_field_contents
        # If field_type == 2, then field_value is the boolean value
        field_ids = self.get_ids(obj_type_id, byte_offset, 12)
        field_names = [self.mcos_names[i-1] for i in field_ids[:, 0]]
        field_dtype = [(name, object) for name in field_names]
        # print("Field Names: ", field_names)

        # Add Handle Class instances (if any)
        byte_offset = np.frombuffer(self.mcos_metadata, dtype=self.byte_order, count=1, offset=24)
        handle_ids = self.get_ids(object_deps[5], byte_offset, 4)[:,0] # Gets object IDs of handle instances
        if handle_ids.size != 0:
            i = 0
            while i < handle_ids.size:
                field_dtype.append((f"handle_{i+1}", object)) # Appending handle instances as properties to class
                i += 1
        
        # Extract fields
        # TODO: Check object representation
        obj_props = np.empty((1, ), dtype=field_dtype)
        i = 0
        
        while i < field_ids.shape[0]:
            if field_ids[i, 1] == 1:
                obj_props[0][field_names[i-1]] = self.mcos_field_contents[field_ids[i, 2], 0]
            elif field_ids[i, 1] == 2:
                obj_props[0, i] = np.array(field_ids[2], dtype=np.bool_)
            else:
                raise ValueError("Could not determine property type. Got type: {field_ids[i, 1]}")
            i += 1

        # Extract Handles
        j = 0
        while j < handle_ids.size:
            class_id, object_id = self.get_handle_class_instance(handle_ids[j])
            obj_props[0, i] = self.read_mcos_arrays(object_id, class_id, dims = np.array([1, 1]))

        return obj_props

    cdef cnp.ndarray read_mcos_arrays(self, cnp.ndarray object_ids, int class_id, cnp.ndarray dims):
        ''' Read MCOS arrays from subsystem
        '''
        cdef:
            int i = 0
            cnp.ndarray[cnp.uint32_t, ndim=1] object_deps
            object class_name 

        # Extract props for all objects in object array
        props_list = []
        while i < object_ids.size:
            object_deps = self.get_object_dependencies(object_ids[i])
            obj_props = self.extract_properties(object_deps)
            props_list.append(obj_props)
            i += 1
        
        obj_props = np.array(props_list).reshape(dims)
        obj_default_props = self.mcos_field_contents[-1, 0][class_id, 0]
        class_name = self.get_class_name(class_id)
        obj_props = convert_object(obj_props, obj_default_props, class_name, self.byte_order)

        # Purpose of these fields are unknown
        _u3 = self.mcos_field_contents[-3, 0][class_id, 0]
        _u2 = self.mcos_field_contents[-2, 0][class_id, 0]

        # TODO: Build object representation
        # print("Class Name: ", class_name)
        # print("Default Properties: ", obj_default_props)
        # print("U3: ", _u3)
        # print("U2: ", _u2)
        return obj_props


#* Below functions just placeholders for now, just for readability
def mat_datetime(props):
    
    data = props["data"][0,0]
    if data.size == 0:
        return props
    
    # if "tz" in props.dtype.names:
    #     tz = props["tz"].item()
    #     py_tz = pytz.timezone(tz)
    #     now_utc = datetime.now(timezone.utc)
    #     now_local = now_utc.astimezone(py_tz)
    #     offset = int(now_local.utcoffset().total_seconds()) * 1000
    # else:
    #     offset = 0

    millis = data.real + data.imag * 1e3
    millis = np.array(millis, dtype="float64")
    props["data"][0,0] = millis.astype("datetime64[ms]")
    return props

def mat_duration(props, defaults):
    """Initialize the MatDuration object"""

    millis = props["millis"][0]
    if millis.size == 0:
        return props

    if "fmt" in props.dtype.names:
        fmt = props["fmt"]
    else:
        fmt = defaults["fmt"]

    if fmt == "s":
        count = millis / 1000  # Seconds
        dur = count.astype("timedelta64[s]")
    elif fmt == "m":
        count = millis / 60000  # Minutes
        dur = count.astype("timedelta64[m]")
    elif fmt == "h":
        count = millis / 3600000  # Hours
        dur = count.astype("timedelta64[h]")
    elif fmt == "d":
        count = millis / 86400000  # Days
        dur = count.astype("timedelta64[D]")
    else:
        count = millis
        dur = count.astype("timedelta64[ms]")
        # Default case

    props["millis"] = dur
    return props

def parse_string(data, byte_order):
    """Parse string data from MATLAB file"""

    data = data["any"][0]
    version = data[0, 0]
    if version != 1:
        print("Unsupported version for string data")

    ndims = data[0, 1]
    shape = data[0, 2 : 2 + ndims]
    num_strings = np.prod(shape)
    char_counts = data[0, 2 + ndims : 2 + ndims + num_strings]
    offset = 2 + ndims + num_strings  # start of string data
    byte_data = data[0, offset:].tobytes()

    strings = []
    pos = 0
    encoding = "utf-16-le" if byte_order[0] == "<" else "utf-16-be"
    for char_count in char_counts:
        byte_length = char_count * 2  # UTF-16 encoding
        extracted_string = byte_data[pos : pos + byte_length].decode(encoding)
        strings.append(extracted_string)
        pos += byte_length

    return np.reshape(strings, shape, order="F")

def convert_object(props, defaults, class_name, byte_order):
    if class_name == "datetime":
        return mat_datetime(props)
    elif class_name == "duration":
        return mat_duration(props, defaults)
    elif class_name == "string":
        if "any" in props.dtype.names:
            props = parse_string(props, byte_order)
    
    return props
