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
   issue a warning.''',
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

'''

 Note on architecture
======================

There are two set of parameters relevant for reading files.  The first
are *file read parameters* - containing options that are common for
reading the whole file, and therefore every variable within that
file. At the moment these are:

* dtypes (derived from order code)
* processor_func
* oned_as
* byte_order
* main_stream
* struct_as_record (matlab 5 files)
* class_dtypes (derived from order code, matlab 5 files)
* codecs (matlab 5 files)
* uint16_codec (matlab 5 files)

The other set of parameters are those that apply only the the current
variable being read - *variable read parameters*:

* header
* var_stream
* next_position

Then, there can be, for each element in a matrix, *element read
parameters*.  An element is, for example, one element in a Matlab cell
array.  At the moment these are:

* mat_dtype

The file-reading object contains the *file read parameters*.  There is
also a variable-reading object that contains the *file read parameters*
and the *variable read parameters*.  

'''


def convert_dtypes(dtype_template, order_code):
    ''' Convert dtypes in mapping to given order

    Parameters
    ----------
    dtype_template : mapping
       mapping with values returning numpy dtype from ``np.dtype(val)``
    order_code : str
       an order code suitable for using in ``dtype.newbyteorder()``

    Returns
    -------
    dtypes : mapping
       mapping where values have been replaced by
       ``np.dtype(val).newbyteorder(order_code)``
       
    '''
    dtypes = dtype_template.copy()
    for k in dtypes:
        dtypes[k] = np.dtype(dtypes[k]).newbyteorder(order_code)
    return dtypes


def read_dtype(mat_stream, a_dtype):
    ''' Generic get of byte stream data of known type

    Parameters
    ----------
    mat_stream : file-like object
    a_dtype : dtype
       dtype of array to read.  `a_dtype` is assumed to be correct
       endianness

    Returns
    -------
    arr : array
    '''
    num_bytes = a_dtype.itemsize
    arr = np.ndarray(shape=(),
                     dtype=a_dtype,
                     buffer=mat_stream.read(num_bytes),
                     order='F')
    return arr


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

ByteOrder = np.deprecate(ByteOrder, message="""
We no longer use the ByteOrder class, and deprecate it; we will remove
it in future versions of scipy.  Please use the
scipy.io.matlab.byteordercodes module instead.
""")


class FileReadParameters(object):
    ''' Container for file read parameters '''
    def __init__(self, byte_order, mat_dtype, squeeze_me,
                 chars_as_strings, struct_as_record):
        self.byte_order = byte_order
        self.mat_dtype = mat_dtype
        self.squeeze_me = squeeze_me
        self.chars_as_strings = chars_as_strings
        self.struct_as_record = struct_as_record
        
    @classmethod
    def from_object(klass, obj):
        return klass(
            obj.byte_order,
            obj.mat_dtype,
            obj.squeeze_me,
            obj.chars_as_strings,
            obj.struct_as_record)


def process_element(arr, file_read_opts, mat_dtype):
    if file_read_opts.chars_as_strings and arr.dtype.kind == 'U':
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
                arr[...,i] = chars_to_str(str_arr[i])
        else: # return string
            arr = chars_to_str(arr)
    if file_read_opts.mat_dtype:
        # Apply options to replicate matlab's (TM)
        # load into workspace
        if mat_dtype is not None:
            arr = arr.astype(mat_dtype)
    if file_read_opts.squeeze_me:
        arr = np.squeeze(arr)
        if not arr.size:
            arr = np.array([])
        elif not arr.shape and arr.dtype.isbuiltin: # 0d coverted to scalar
            arr = arr.item()
    return arr


def chars_to_str(str_arr):
    ''' Convert string array to string '''
    dt = np.dtype('U' + str(small_product(str_arr.shape)))
    return np.ndarray(shape=(),
                      dtype = dt,
                      buffer = str_arr.copy()).item()


class MatFileReader(object):
    """ Base object for reading mat files

    To make this class functional, you will need to override the
    following methods:

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
        else:
            byte_order = boc.to_numpy_code(byte_order)
        self._byte_order = byte_order
        self._struct_as_record = struct_as_record
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

    def get_byte_order(self):
        'Read only byte_order property'
        return self._byte_order
    byte_order = property(get_byte_order,
                          None,
                          None,
                          'get byte_order code')

    def get_mat_dtype(self):
        return self._mat_dtype
    mat_dtype = property(get_mat_dtype,
                         None,
                         None,
                         'get mat_dtype property')

    def get_squeeze_me(self):
        return self._squeeze_me
    squeeze_me = property(get_squeeze_me,
                          None,
                          None,
                          'get squeeze me property')

    def get_chars_as_strings(self):
        return self._chars_as_strings
    chars_as_strings = property(get_chars_as_strings,
                                None,
                                None,
                                'get chars_as_strings property')

    def get_struct_as_record(self):
        return self._struct_as_record
    struct_as_record = property(get_struct_as_record,
                                None,
                                None,
                                'get struct_as_record property')

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
        params = FileReadParameters.from_object(self)
        def func(arr, getter):
            return process_element(arr, params, getter.mat_dtype)
        return func
    
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
            reader = self.get_reader()
            name = reader.name
            if variable_names and name not in variable_names:
                self.mat_stream.seek(var_params.next_position)
                continue
            try:
                res = self.get_variable(reader)
            except MatReadError, err:
                warnings.warn(
                    'Unreadable variable "%s", because "%s"' % \
                    (name, err),
                    Warning, stacklevel=2)
                res = "Read error: %s" % err
            self.mat_stream.seek(reader.next_position)
            mdict[name] = res
            if reader.is_global:
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

    def get_reader(self):
        raise NotImplementedError

    def get_variable(self, reader):
        raise NotImplementedError


class MatMatrixGetter(object):
    """ Base class for matrix getters

    Getters do store state of what they are reading on initialization,
    and therefore need to be initialized before each read.  They are
    one-shot objects.

    The getter stores the:
       mat_stream (file-like)
          stream from which to read data
       array_reader (object for reading arrays from the stream)
          used to (base)
             contain reference to 'processor_func'
             pass in mat_stream
             pass in dtypes
          used to (mio4)
             mat_stream used for read_array
             dtypes not used
          mio5 also:
             pass reference to class_dtypes, codecs
             provide read_element implementation
             provide struct_as_record read flag to struct getter
       dtypes (little- or big-endian versions of matlab dtypes)
       header (information about the current matrix being read)
       name (from the header - name of matrix being read)

    Also has:
       next_position (position for next matrix at base level)

    Does:
       get_array (just a shell to call get_raw_array, and run processor
       func)
       get_raw_array (actual work for array reading)
       to_next (move to ``next_position`` above)
       
    Parameters
    ----------
    array_reader : array reading object
       Containing stream, dtypes.
    header : mapping
       header dictionary for matrix being read
    """

    def __init__(self, array_reader, header):
        self.mat_stream = array_reader.mat_stream
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
