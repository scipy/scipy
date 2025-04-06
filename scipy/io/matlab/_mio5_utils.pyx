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
    cdef readonly object class_name
    cdef int object_id
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
        if mc == mxOPAQUE_CLASS:
            header.name = self.read_int8_string()
            type_system_name = self.read_int8_string()
            if type_system_name != b'MCOS':
                raise ValueError('Expecting Type as MATLAB Class Object System')
            header.class_name = PyUnicode_FromString(self.read_int8_string())
            # object metadata stored as a Nx1 uint32 array
            obj_metadata = self.read_mi_matrix()
            header.n_dims = obj_metadata[1,0]
            if header.n_dims > _MAT_MAXDIMS:
                raise ValueError('Too many dimensions (%d) for numpy arrays'
                                 % header.n_dims)
            header.dims = [obj_metadata[2 + i, 0].item() for i in range(header.n_dims)]
            header.object_id = obj_metadata[-2, 0]
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

    cdef read_mi_matrix(self, int process=1, SubsystemReader5 subsystem=None, int is_subsystem=0):
        ''' Read header with matrix at sub-levels

        Combines ``read_header`` and functionality of
        ``array_from_header``.  Applies standard processing of array
        given options set in self.

        Parameters
        ----------
        process : int, optional
           If not zero, apply post-processing on returned array
        subsystem : SubsystemReader5, optional
           If not None, read miMATRIX arrays from subsystem
        is_subsystem : int, optional
           If not zero, check if miMATRIX array is an object reference

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
        return self.array_from_header(header, process, subsystem, is_subsystem)

    cpdef array_from_header(self, VarHeader5 header, int process=1, SubsystemReader5 subsystem=None, int is_subsystem=0):
        ''' Read array of any class, given matrix `header`

        Parameters
        ----------
        header : VarHeader5
           array header object
        process : int, optional
           If not zero, apply post-processing on returned array
        subsystem : SubsystemReader5, optional
            If not None, array from subsystem
        is_subsystem : int, optional
            If not zero, check if array is an object reference

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
            arr = self.read_real_complex(header, subsystem, is_subsystem)
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
            arr = self.read_cells(header)
        elif mc == mxSTRUCT_CLASS:
            arr = self.read_struct(header)
        elif mc == mxOBJECT_CLASS: # like structs, but with classname
            classname = self.read_int8_string().decode('latin1')
            arr = self.read_struct(header)
            arr = mio5p.MatlabObject(arr, classname)
        elif mc == mxFUNCTION_CLASS: # just a matrix of struct type
            arr = self.read_mi_matrix()
            arr = mio5p.MatlabFunction(arr)
            # to make them more re-writeable - don't squeeze
            process = 0
        elif mc == mxOPAQUE_CLASS:
            arr = self.read_opaque(header, subsystem)
            # arr = mio5p.MatlabOpaque(arr)
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

    cpdef cnp.ndarray read_real_complex(self, VarHeader5 header, SubsystemReader5 subsystem=None, int is_subsystem=0):
        ''' Read real / complex matrices from stream '''
        cdef:
            cnp.ndarray res, res_j
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
        if is_subsystem:
            if self.check_array_object(res):
                nested_obj_hdr = VarHeader5()
                nested_obj_hdr.class_name = subsystem.class_names[res[-1]-1]
                nested_obj_hdr.n_dims = res[1]
                nested_obj_hdr.dims = [res[2+i] for i in range(nested_obj_hdr.n_dims)]
                nested_obj_hdr.object_id = res[-2]
                res = self.read_opaque(nested_obj_hdr, subsystem)
                return res
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

    cpdef cnp.ndarray read_struct(self, VarHeader5 header):
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
            with cython.boundscheck(False):
                result[i] = item
        return result.reshape(tupdims).T

    cdef int check_array_object(self, cnp.ndarray arr):
        ''' Checks if numeric array is a MATLAB object reference'''
        #! Need to add this check within all functions calling miMATRIX, e.g. struct or cell
        cdef:
            cnp.uint32_t ref_value, ndims

        if not isinstance(arr, cnp.ndarray):
            return 0
        if arr.dtype != np.uint32:
            return 0

        if arr.shape[0] < 6:
            return 0

        ref_value = arr[0]
        if ref_value != 0xDD000000:
            return 0

        ndims = arr[1]
        shapes = arr[2:2+ndims]
        if ndims <= 1 or len(shapes) != ndims:
            return 0

        return 1

    cdef extract_opaque_from_field(self, int nfields, SubsystemReader5 subsystem, int is_subsystem):
        ''' Extract opaque class type from field '''
        cdef:
            object dt, field_names
            cnp.uint32_t* field_ids
            cnp.uint32_t* field_content_ids
            int process=1

        try:
            field_ids = <cnp.uint32_t*>calloc(nfields, sizeof(cnp.uint32_t))
            field_content_ids = <cnp.uint32_t*>calloc(nfields, sizeof(cnp.uint32_t))
            subsystem.get_field_ids(nfields, field_ids, field_content_ids)
            field_names = [subsystem.names[field_ids[i]-1] for i in range(nfields)]

            dt = [(field_name, object) for field_name in field_names]
            rec_res = np.empty(1, dtype=dt)

            # Save current stream (useful for miCOMPRESSED arrays)
            stream = self.cstream
            for i in range(nfields):
                self.set_stream(subsystem.get_field_content_stream(field_content_ids[i]))
                # MATLAB Datatype object arrays are stored as array of field contents
                # e.g. MATLAB datetime array is stored as a 1x1 datetime object array
                # However, the field contents are stored as a MxN array instead
                rec_res[0][field_names[i]] = self.read_mi_matrix(process, subsystem, is_subsystem)
                #! Need to fix this np.array structure

        finally:
            if field_ids is not NULL:
                free(field_ids)
            if field_content_ids is not NULL:
                free(field_content_ids)

        # Restore the original stream
        self.set_stream(stream)

        return rec_res

    cpdef cnp.ndarray read_opaque(self, VarHeader5 hdr, SubsystemReader5 subsystem):
        ''' Read opaque class type'''
        # Neither res nor the return value of this function are cdef'd as
        # cnp.ndarray, because that only adds useless checks with current
        # Cython (0.23.4).

        cdef:
            int length
            int class_id, type1_id, type2_id, dependency_id, object_ids
            int ndeps, nfields
            object class_name
            int is_subsystem = 1
            int index=0

        length = np.prod(hdr.dims)
        tupdims = tuple(hdr.dims[::-1])
        object_id = hdr.object_id
        obj_res = np.empty(length, dtype=object)
        class_name = hdr.class_name

        for i in range (object_id - length + 1, object_id + 1):
            dependency_id = subsystem.object_list[object_id-1][3]
            # For nested objects or object arrays
            # Dependency ID gives the number of objects used to construct
            ndeps = dependency_id - hdr.object_id

            class_id, type1_id, type2_id = subsystem.object_list[i-1][:3]
            nfields = subsystem.get_num_fields(type1_id, type2_id)
            obj_res[index] = self.extract_opaque_from_field(nfields, subsystem, is_subsystem)
            #* Add Option to parse class raw data
            i += ndeps
            index += 1

        obj_res = obj_res.reshape(tupdims).T
        res = np.empty(1, dtype=[('class', object), ('properties', object)])
        res[0]['class'] = class_name
        res[0]['properties'] = obj_res
        return res


cdef class SubsystemReader5:
    cdef int is_swapped # endian indicator derived from subsystem
    cdef list names # list of all field names and class names
    cdef list class_names # list of class names ordered by class_id from names
    cdef list object_list # list of all objects with their attributes like type_id, dependency_id, and class_id
    cdef list field_content_pos # list of byte markers for field contents
    cdef _streams.GenericStream cstream
    cdef cnp.uint64_t type1_pos, type2_pos

    def __cinit__(self, mat_stream, is_compressed=True):
        ''' Initialize SubsystemReader5'''
        self.cstream = _streams.make_stream(mat_stream)
        if is_compressed:
            self.cstream.seek(8,1) # Read miMATRIX Tag
        self.cstream.seek(48, 1) # Skip Array Headers
        self.is_swapped = self.read_byte_order()
        self.names = []
        self.class_names = []
        self.object_list = []
        self.field_content_pos = []
        self.type1_pos = 0
        self.type2_pos = 0

        self.initialize_subsystem()

    cdef int read_byte_order(self):
        '''Reads byte order from subsystem'''
        cdef char byte_order[8]

        self.cstream.read_into(<void *>byte_order, 8)
        # Check the 3rd and 4th bytes for byte order
        if byte_order[2] == ord('I') and byte_order[3] == ord('M'):
            return '<' == swapped_code
        else:
            return '>' == swapped_code

    cdef inline int cread_full_tag(self,
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

    cdef read_int32_blocks(self, cnp.uint32_t* arr, int size):
        '''Read int32 integers into a pre-allocated array.'''
        self.cstream.read_into(<void *>arr, size * 4)
        if self.is_swapped:
            for i in range(size):
                arr[i] = byteswap_u4(arr[i])

    cdef cread_names(self, cnp.uint32_t *fc_len_ptr, cnp.uint32_t *byte_end_ptr):
        '''Read field and class names from the stream.'''
        cdef:
            int i
            cnp.uint32_t byte_count
            object names
            char *names_ptr
            object name

        # Calculate the number of bytes to read
        byte_count = byte_end_ptr[0] - 40 # 40 is the offset to the start of cell 1

        names = self.cstream.read_string(
                        byte_count,
                        <void **>&names_ptr,
                        copy=True)

        # Names are stored continuously and null terminated
        for i in range(fc_len_ptr[0]):
            name = PyUnicode_FromString(names_ptr)
            self.names.append(name)
            while names_ptr[0] != 0:
                names_ptr += 1
            names_ptr += 1

    cdef cread_class_names(self, cnp.uint32_t *region_ptr):
        '''Read class names from the stream.'''
        cdef:
            int i
            cnp.uint32_t byte_count
            int n_blocks
            cnp.uint32_t u4s[4]
            int class_index

        # Calculate the number of bytes to read
        byte_count = region_ptr[1] - region_ptr[0]
        n_blocks = byte_count // 16
        # Read and discard first block
        self.cstream.read_into(<void *>u4s, 16)

        # Class names stored as (0, class_index, 0, 0)
        # class_index points to class name in list obtained in cread_names()
        # indexing starts from 1
        for i in range(0, n_blocks-1):
            self.cstream.read_into(<void *>u4s, 16)
            if self.is_swapped:
                class_index = byteswap_u4(u4s[1])
            else:
                class_index = u4s[1]
            self.class_names.append(self.names[class_index - 1])

        # Saving this region offset as it contains field IDs for type 1 objects
        self.type1_pos = self.cstream.tell()

    cdef cread_object_ids(self, cnp.uint32_t *region_ptr):
        '''Read object IDs from the stream.'''
        cdef:
            int i
            cnp.uint32_t byte_count
            int n_blocks
            cnp.uint32_t u4s[6]
            int class_id, type1_id, type2_id, dependency_id

        byte_count = region_ptr[2] - region_ptr[1]
        self.cstream.seek(byte_count, 1) # Skip to Region 3

        # Calculate the number of bytes to read
        byte_count = region_ptr[3] - region_ptr[2]
        n_blocks = byte_count // 24
        # Read and discard first block
        self.cstream.read_into(<void *>u4s, 24)

        # Format is (class_id, 0, 0, type1_id, type2_id, dependency_id)
        for i in range(0, n_blocks-1):
            self.cstream.read_into(<void *>u4s, 24)
            if self.is_swapped:
                class_id = byteswap_u4(u4s[0])
                type1_id = byteswap_u4(u4s[3])
                type2_id = byteswap_u4(u4s[4])
                dependency_id = byteswap_u4(u4s[5])
            else:
                class_id = u4s[0]
                type1_id = u4s[3]
                type2_id = u4s[4]
                dependency_id = u4s[5]
            self.object_list.append(
                (class_id, type1_id, type2_id, dependency_id)
            )

        # Saving this region offset as it contains field IDs for type 2 objects
        self.type2_pos = self.cstream.tell()

    cdef extract_field_content_pos(self, cnp.uint32_t *region_ptr, cnp.uint64_t *byte_end_ptr):
        '''Extracting field content positions from the stream.'''
        #! Note: Need to check with compressed MAT-files
        cdef:
            cnp.uint32_t mdtype, byte_count
            int nbytes

        nbytes = region_ptr[7] - region_ptr[3] + 8
        self.cstream.seek(nbytes, 1) # Skip to start of Cell 3 (field contents)

        # Field IDs are 0-indexed
        while self.cstream.tell() < byte_end_ptr[0]:
            # Read the next 8 bytes
            self.cread_full_tag(&mdtype, &byte_count)
            if mdtype != miMATRIX:
                raise TypeError('Expecting matrix here')
            self.field_content_pos.append(self.cstream.tell()-8)
            self.cstream.seek(byte_count, 1)

    cdef initialize_subsystem(self):
        '''Initialize the subsystem reader and extract metadata'''
        cdef:
            cnp.uint32_t unknown_flag, fc_len
            cnp.uint32_t offsets[8]
            cnp.uint64_t total_byte_count

        self.cstream.seek(136, 1) # Skip to object metadata
        self.cread_full_tag(&unknown_flag, &fc_len)
        if unknown_flag != miMATRIX:
            raise TypeError('Expecting matrix here')
        total_byte_count = self.cstream.tell() + fc_len

        self.cstream.seek(96, 1) # Skip to cell 1 contents in array
        self.cread_full_tag(&unknown_flag, &fc_len)
        if unknown_flag != 0x00000004:
            raise ValueError('Unknown flag in Subsystem Table of Contents'
             'Expecting 4, got %d' % unknown_flag)

        self.read_int32_blocks(offsets, 8)
        self.cread_names(&fc_len, offsets)
        self.cread_class_names(offsets)
        self.cread_object_ids(offsets)
        self.extract_field_content_pos(offsets, &total_byte_count)
    
    cdef int get_num_fields(self, int type1_id, int type2_id):
        '''Get number of fields for object ID'''
        cdef:
            int i, mod8
            int obj_type_id
            int nfields = 0
            cnp.uint32_t u4s[3]

        if type1_id != 0 and type2_id == 0:
            self.cstream.seek(self.type1_pos)
            obj_type_id = type1_id
        elif type1_id == 0 and type2_id != 0:
            self.cstream.seek(self.type2_pos)
            obj_type_id = type2_id
        else:
            raise ValueError('Unknown object dependency type. Feature not supported yet'
            'Received type1_id = %d, type2_id = %d' % (type1_id, type2_id))

        self.cstream.seek(8, 1) # Skip first 8 bytes
        # Seek to block corresponding to object_ID
        while obj_type_id - 1 > 0:
            self.cstream.read_into(<void *>&nfields, 4)
            if self.is_swapped:
                nfields = byteswap_u4(nfields)
            mod8 = (nfields*12 - 4) % 8
            self.cstream.seek(nfields*12 + mod8, 1)
            obj_type_id -= 1

        self.cstream.read_into(<void *>&nfields, 4)
        if self.is_swapped:
                nfields = byteswap_u4(nfields)
        return nfields

    cdef get_field_ids(self, int nfields, cnp.uint32_t *field_ids, cnp.uint32_t *field_content_ids):
        '''Get field IDs and field content IDs for a given object ID'''
        cdef:
            cnp.uint32_t u4s[3]
            int i = 0

        while nfields > 0:
            self.cstream.read_into(<void *>u4s, 12)
            if self.is_swapped:
                field_ids[i] = byteswap_u4(u4s[0])
                field_content_ids[i] = byteswap_u4(u4s[2])
            else:
                field_ids[i] = u4s[0]
                field_content_ids[i] = u4s[2]
            i += 1
            nfields -= 1

    cdef get_field_content_stream(self, int field_content_id):
        '''Get field content stream for a given field content ID
        Used for setting stream when MAT-file is compressed
        '''

        cdef:
            cnp.uint32_t byte_count
            cnp.uint64_t pos
            cnp.uint32_t u4s[2]
            object stream

        self.cstream.seek(self.field_content_pos[field_content_id])
        stream = _streams.make_stream(self.cstream)
        return stream
