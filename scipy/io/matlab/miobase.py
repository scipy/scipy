# Authors: Travis Oliphant, Matthew Brett

"""
Base classes for matlab (TM) file stream reading
"""
import warnings

import numpy as np

from scipy.ndimage import doccer

import byteordercodes as boc

class MatReadError(Exception): pass

doc_dict = \
    {'file_arg':
         '''file_name : string
   Name of the mat file (do not need .mat extension if
   appendmat==True) If name not a full path name, search for the
   file on the sys.path list and use the first one found (the
   current directory is searched first).  Can also pass open
   file-like object''',
     'append_arg':
         '''appendmat : {True, False} optional
   True to append the .mat extension to the end of the given
   filename, if not already present''',
     'basename_arg':
         '''base_name : string, optional, unused
   base name for unnamed variables.  The code no longer uses
   this.  We deprecate for this version of scipy, and will remove
   it in future versions''',
     'load_args':
         '''byte_order : {None, string}, optional
   None by default, implying byte order guessed from mat
   file. Otherwise can be one of ('native', '=', 'little', '<',
   'BIG', '>')
mat_dtype : {False, True} optional
   If True, return arrays in same dtype as would be loaded into
   matlab (instead of the dtype with which they are saved)
squeeze_me : {False, True} optional
   whether to squeeze unit matrix dimensions or not
chars_as_strings : {True, False} optional
   whether to convert char arrays to string arrays
matlab_compatible : {False, True}
   returns matrices as would be loaded by matlab (implies
   squeeze_me=False, chars_as_strings=False, mat_dtype=True,
   struct_as_record=True)''',
     'struct_arg':
         '''struct_as_record : {False, True} optional
   Whether to load matlab structs as numpy record arrays, or as
   old-style numpy arrays with dtype=object.  Setting this flag to
   False replicates the behaviour of scipy version 0.6 (returning
   numpy object arrays).  The preferred setting is True, because it
   allows easier round-trip load and save of matlab files.  In a
   future version of scipy, we will change the default setting to
   True, and following versions may remove this flag entirely.  For
   now, we set the default to False, for backwards compatibility, but
   issue a warning.  Note that non-record arrays cannot be exported
   via savemat.''',
     'matstream_arg':
         '''mat_stream : file-like
   object with file API, open for reading''',
     'long_fields':
         '''long_field_names : boolean, optional, default=False
   * False - maximum field name length in a structure is 31 characters
     which is the documented maximum length
   * True - maximum field name length in a structure is 63 characters
     which works for Matlab 7.6''',
     'do_compression':
         '''do_compression : {False, True} bool, optional
   Whether to compress matrices on write. Default is False''',
     'oned_as':
         '''oned_as : {'column', 'row'} string, optional
   If 'column', write 1D numpy arrays as column vectors
   If 'row', write 1D numpy arrays as row vectors''',
     'unicode_strings':
         '''unicode_strings : {True, False} boolean, optional
   If True, write strings as Unicode, else matlab usual encoding'''}

docfiller = doccer.filldoc(doc_dict)


def small_product(arr):
    ''' Faster than product for small arrays '''
    res = 1
    for e in arr:
        res *= e
    return res


def get_matfile_version(fileobj):
    ''' Return major, minor tuple depending on apparent mat file type

    Where:

     #. 0,x -> version 4 format mat files
     #. 1,x -> version 5 format mat files
     #. 2,x -> version 7.3 format mat files (HDF format)

    Parameters
    ----------
    fileobj : {file-like}
              object implementing seek() and read()

    Returns
    -------
    major_version : {0, 1, 2}
                    major matlab file format version
    minor_version : int
                    major matlab file format version

    Notes
    -----
    Has the side effect of setting the file read pointer to 0
    '''
    # Mat4 files have a zero somewhere in first 4 bytes
    fileobj.seek(0)
    mopt_bytes = np.ndarray(shape=(4,),
                           dtype=np.uint8,
                           buffer = fileobj.read(4))
    if 0 in mopt_bytes:
        fileobj.seek(0)
        return (0,0)

    # For 5 format or 7.3 format we need to read an integer in the
    # header. Bytes 124 through 128 contain a version integer and an
    # endian test string
    fileobj.seek(124)
    tst_str = fileobj.read(4)
    fileobj.seek(0)
    maj_ind = int(tst_str[2] == 'I')
    maj_val = ord(tst_str[maj_ind])
    min_val = ord(tst_str[1-maj_ind])
    ret = (maj_val, min_val)
    if maj_val in (1, 2):
        return ret
    else:
        raise ValueError('Unknown mat file type, version %s, %s'
                         % ret)


class MatReadError(Exception): pass


def matdims(arr, oned_as='column'):
    ''' Determine equivalent matlab dimensions for given array 
    
    Parameters
    ----------
    arr : ndarray
    oned_as : {'column', 'row'} string, optional

    Returns
    -------
    dims : shape as matlab expects

    Examples
    --------
    >>> matdims(np.array(1)) # numpy scalar
    (1, 1)
    >>> matdims(np.array([1])) # 1d array, 1 element
    (1, 1)
    >>> matdims(np.array([1,2])) # 1d array, 2 elements
    (2, 1)
    >>> matdims(np.array([[2],[3]])) # 2d array, column vector
    (2, 1)
    >>> matdims(np.array([[2,3]])) # 2d array, row vector
    (1, 2)
    >>> matdims(np.array([[[2,3]]])) # 3d array, rowish vector
    (1, 1, 2)
    >>> matdims(np.array([])) # empty 1d array
    (0, 0)
    >>> matdims(np.array([[]])) # empty 2d
    (0, 0)
    >>> matdims(np.array([[[]]])) # empty 3d
    (0, 0, 0)

    Optional argument flips 1d shape behavior

    >>> matdims(np.array([1,2]), 'row') # 1d array, 2 elements
    (1, 2)

    The argument has to make sense though

    >>> matdims(np.array([1,2]), 'bizarre')
    Traceback (most recent call last):
       ...
    ValueError: 1D option "bizarre" is strange

    Notes
    -----
    We had to decide what shape a 1 dimensional array would be by
    default.  ``np.atleast_2d`` thinks it is a row vector.  The
    default for a vector in matlab (e.g. ``>> 1:12``) is a row vector.

    Versions of scipy up to and including 0.7 resulted (accidentally)
    in 1d arrays being read as column vectors.  For the moment, we
    maintain the same tradition here.
    '''
    if arr.size == 0: # empty
        return (0,) * np.max([arr.ndim, 2])
    shape = arr.shape
    if shape == (): # scalar
        return (1,1)
    if len(shape) == 1: # 1D
        if oned_as == 'column':
            return shape + (1,)
        elif oned_as == 'row':
            return (1,) + shape
        else:
            raise ValueError('1D option "%s" is strange'
                             % oned_as)
    return shape


class ByteOrder(object):
    ''' Namespace for byte ordering '''
    little_endian = boc.sys_is_le
    native_code = boc.native_code
    swapped_code = boc.swapped_code
    to_numpy_code = boc.to_numpy_code

ByteOrder = np.deprecate_with_doc("""
We no longer use the ByteOrder class, and deprecate it; we will remove
it in future versions of scipy.  Please use the
scipy.io.matlab.byteordercodes module instead.
""")(ByteOrder)


class MatStreamAgent(object):
    ''' Base object for readers / getters from mat file streams

    Attaches to initialized stream

    Base class for "getters" - which do store state of what they are
    reading on initialization, and therefore need to be initialized
    before each read, and "readers" which do not store state, and only
    need to be initialized once on object creation

    Implements common array reading functions

    Inputs mat_steam - MatFileReader object
    '''

    def __init__(self, mat_stream):
        self.mat_stream = mat_stream

    def read_dtype(self, a_dtype):
        ''' Generic get of byte stream data of known type

        Inputs
        a_dtype     - dtype of array

        a_dtype is assumed to be correct endianness
        '''
        num_bytes = a_dtype.itemsize
        arr = np.ndarray(shape=(),
                         dtype=a_dtype,
                         buffer=self.mat_stream.read(num_bytes),
                         order='F')
        return arr

    def read_ztstring(self, num_bytes):
        return self.mat_stream.read(num_bytes).strip('\x00')


class MatFileReader(MatStreamAgent):
    """ Base object for reading mat files

    To make this class functional, you will need to override the
    following methods:

    set_dtypes              - sets data types defs from byte order
    matrix_getter_factory   - gives object to fetch next matrix from stream
    guess_byte_order        - guesses file byte order from file
    """

    @docfiller
    def __init__(self, mat_stream,
                 byte_order=None,
                 mat_dtype=False,
                 squeeze_me=False,
                 chars_as_strings=True,
                 matlab_compatible=False,
                 struct_as_record=None
                 ):
        '''
        Initializer for mat file reader

        mat_stream : file-like
            object with file API, open for reading
    %(load_args)s
        '''
        # Initialize stream
        self.mat_stream = mat_stream
        self.dtypes = {}
        if not byte_order:
            byte_order = self.guess_byte_order()
        self.order_code = byte_order # sets dtypes and other things too
        if matlab_compatible:
            self.set_matlab_compatible()
        else:
            self._squeeze_me = squeeze_me
            self._chars_as_strings = chars_as_strings
            self._mat_dtype = mat_dtype
            self.processor_func = self.get_processor_func()

    def set_matlab_compatible(self):
        ''' Sets options to return arrays as matlab (tm) loads them '''
        self._mat_dtype = True
        self._squeeze_me = False
        self._chars_as_strings = False
        self.processor_func = self.get_processor_func()

    def get_mat_dtype(self):
        return self._mat_dtype
    def set_mat_dtype(self, mat_dtype):
        self._mat_dtype = mat_dtype
        self.processor_func = self.get_processor_func()
    mat_dtype = property(get_mat_dtype,
                         set_mat_dtype,
                         None,
                         'get/set mat_dtype property')

    def get_squeeze_me(self):
        return self._squeeze_me
    def set_squeeze_me(self, squeeze_me):
        self._squeeze_me = squeeze_me
        self.processor_func = self.get_processor_func()
    squeeze_me = property(get_squeeze_me,
                          set_squeeze_me,
                          None,
                          'get/set squeeze me property')

    def get_chars_as_strings(self):
        return self._chars_as_strings
    def set_chars_as_strings(self, chars_as_strings):
        self._chars_as_strings = chars_as_strings
        self.processor_func = self.get_processor_func()
    chars_as_strings = property(get_chars_as_strings,
                                set_chars_as_strings,
                                None,
                                'get/set chars_as_strings property')

    def get_order_code(self):
        return self._order_code
    def set_order_code(self, order_code):
        order_code = boc.to_numpy_code(order_code)
        self._order_code = order_code
        self.set_dtypes()
    order_code = property(get_order_code,
                          set_order_code,
                          None,
                          'get/set order code')

    def set_dtypes(self):
        ''' Set dtype endianness. In this case we have no dtypes '''
        pass

    def convert_dtypes(self, dtype_template):
        dtypes = dtype_template.copy()
        for k in dtypes:
            dtypes[k] = np.dtype(dtypes[k]).newbyteorder(self.order_code)
        return dtypes

    def matrix_getter_factory(self):
        assert False, 'Not implemented'

    def file_header(self):
        return {}

    def guess_byte_order(self):
        ''' As we do not know what file type we have, assume native '''
        return ByteOrder.native_code

    def get_processor_func(self):
        ''' Processing to apply to read matrices

        Function applies options to matrices. We have to pass this
        function into the reader routines because Mat5 matrices
        occur as submatrices - in cell arrays, structs and objects -
        so we will not see these in the main variable getting routine
        here.

        The read array is the first argument.
        The getter, passed as second argument to the function, must
        define properties, iff mat_dtype option is True:

        mat_dtype    - data type when loaded into matlab (tm)
                       (None for no conversion)

        func returns the processed array
        '''

        def func(arr, getter):
            if arr.dtype.kind == 'U' and self.chars_as_strings:
                # Convert char array to string or array of strings
                dims = arr.shape
                if len(dims) >= 2: # return array of strings
                    n_dims = dims[:-1]
                    last_dim = dims[-1]
                    str_arr = arr.reshape(
                        (small_product(n_dims),
                         last_dim))
                    dtstr = 'U%d' % (last_dim and last_dim or 1)
                    arr = np.empty(n_dims, dtype=dtstr)
                    for i in range(0, n_dims[-1]):
                        arr[...,i] = self.chars_to_str(str_arr[i])
                else: # return string
                    arr = self.chars_to_str(arr)
            if self.mat_dtype:
                # Apply options to replicate matlab's (TM)
                # load into workspace
                if getter.mat_dtype is not None:
                    arr = arr.astype(getter.mat_dtype)
            if self.squeeze_me:
                arr = np.squeeze(arr)
                if not arr.size:
                    arr = np.array([])
                elif not arr.shape and arr.dtype.isbuiltin: # 0d coverted to scalar
                    arr = arr.item()
            return arr
        return func

    def chars_to_str(self, str_arr):
        ''' Convert string array to string '''
        dt = np.dtype('U' + str(small_product(str_arr.shape)))
        return np.ndarray(shape=(),
                          dtype = dt,
                          buffer = str_arr.copy()).item()

    def get_variables(self, variable_names=None):
        ''' get variables from stream as dictionary

        variable_names   - optional list of variable names to get

        If variable_names is None, then get all variables in file
        '''
        if isinstance(variable_names, basestring):
            variable_names = [variable_names]
        self.mat_stream.seek(0)
        mdict = self.file_header()
        mdict['__globals__'] = []
        while not self.end_of_stream():
            getter = self.matrix_getter_factory()
            name = getter.name
            if variable_names and name not in variable_names:
                getter.to_next()
                continue
            try:
                res = getter.get_array()
            except MatReadError, err:
                warnings.warn(
                    'Unreadable variable "%s", because "%s"' % \
                    (name, err),
                    Warning, stacklevel=2)
                res = "Read error: %s" % err
                getter.to_next()
            mdict[name] = res
            if getter.is_global:
                mdict['__globals__'].append(name)
            if variable_names:
                variable_names.remove(name)
                if not variable_names:
                    break
        return mdict

    def end_of_stream(self):
        b = self.mat_stream.read(1)
        curpos = self.mat_stream.tell()
        self.mat_stream.seek(curpos-1)
        return len(b) == 0


class MatMatrixGetter(MatStreamAgent):
    """ Base class for matrix getters

    Getters are stateful versions of agents, and record state of
    current read on initialization, so need to be created for each
    read - one-shot objects.

    MatrixGetters are initialized with the content of the matrix
    header

    Accepts
    array_reader - array reading object (see below)
    header       - header dictionary for matrix being read
    """

    def __init__(self, array_reader, header):
        super(MatMatrixGetter, self).__init__(array_reader.mat_stream)
        self.array_reader = array_reader
        self.dtypes = array_reader.dtypes
        self.header = header
        self.name = header['name']

    def get_array(self):
        ''' Gets an array from matrix, and applies any necessary processing '''
        arr = self.get_raw_array()
        return self.array_reader.processor_func(arr, self)

    def get_raw_array(self):
        assert False, 'Not implemented'

    def to_next(self):
        self.mat_stream.seek(self.next_position)


class MatArrayReader(MatStreamAgent):
    ''' Base class for array readers

    The array_reader contains information about the current reading
    process, such as byte ordered dtypes and the processing function
    to apply to matrices as they are read, as well as routines for
    reading matrix compenents.
    '''

    def __init__(self, mat_stream, dtypes, processor_func):
        self.mat_stream = mat_stream
        self.dtypes = dtypes
        self.processor_func = processor_func

    def matrix_getter_factory(self):
        assert False, 'Not implemented'


class MatStreamWriter(object):
    ''' Base object for writing to mat files '''
    def __init__(self, file_stream, arr, name, oned_as):
        self.file_stream = file_stream
        self.arr = arr
        dt = self.arr.dtype
        if not dt.isnative:
            self.arr = self.arr.astype(dt.newbyteorder('='))
        self.name = name
        self.oned_as = oned_as

    def rewind(self):
        self.file_stream.seek(0)

    def arr_dtype_number(self, num):
        ''' Return dtype for given number of items per element'''
        return np.dtype(self.arr.dtype.str[:2] + str(num))

    def arr_to_chars(self):
        ''' Convert string array to char array '''
        dims = list(self.arr.shape)
        if not dims:
            dims = [1]
        dims.append(int(self.arr.dtype.str[2:]))
        self.arr = np.ndarray(shape=dims,
                              dtype=self.arr_dtype_number(1),
                              buffer=self.arr)

    def write_bytes(self, arr):
        self.file_stream.write(arr.tostring(order='F'))

    def write_string(self, s):
        self.file_stream.write(s)


class MatFileWriter(object):
    ''' Base class for Mat file writers '''
    def __init__(self, file_stream):
        self.file_stream = file_stream
