# Authors: Travis Oliphant, Matthew Brett

"""
Base classes for matlab (TM) file stream reading
"""

import sys

from numpy import *

try:
    import scipy.sparse
    have_sparse = 1
except ImportError:
    have_sparse = 0


def small_product(arr):
    ''' Faster than product for small arrays '''
    res = 1
    for e in arr:
        res *= e
    return res

class ByteOrder(object):
    ''' Namespace for byte ordering '''
    little_endian = sys.byteorder == 'little'
    native_code = little_endian and '<' or '>'
    swapped_code = little_endian and '>' or '<'
    
    def to_numpy_code(code):
        if code is None:
            return ByteOrder.native_code
        if code in ('little', '<', 'l', 'L'):
            return '<'
        elif code in ('BIG', '>', 'B', 'b'):
            return '>'
        elif code in ('native', '='):
            return ByteOrder.native_code
        elif code in ('swapped'):
            return ByteOrder.swapped_code
        else:
            raise ValueError, 'We cannot handle byte order %s' % byte_order
    to_numpy_code = staticmethod(to_numpy_code)


class MatStreamAgent(object):
    ''' Base object for readers / getters from mat file streams

    Attaches to initialized stream

    Base class for "getters" - which do store state of what they are
    reading on itialization, and therefore need to be initialized
    before each read, and "readers" which do not store state, and only
    need to be initialized once on object creation

    Implements common array reading functions

    Inputs @mat_steam - MatFileReader object
    '''

    def __init__(self, mat_stream):
        self.mat_stream = mat_stream

    def read_dtype(self, a_dtype):
        ''' Generic get of byte stream data of known type

        Inputs
        @a_dtype     - dtype of array
        a_dtype is assumed to be correct endianness
        '''
        num_bytes = a_dtype.itemsize
        arr = ndarray(shape=(),
                      dtype=a_dtype,
                      buffer=self.mat_stream.read(num_bytes),
                      order='F')
        return arr

    def read_ztstring(self, num_bytes):
        return self.mat_stream.read(num_bytes).strip('\x00')

        
class MatFileReader(MatStreamAgent):
    """ Base object for reading mat files

    @initialized byte stream object  - file io interface object
    @byte_order         - byte order ('native', 'little', 'BIG')
                          in ('native', '=')
                          or in ('little', '<')
                          or in ('BIG', '>')
    @base_name          - base name for unnamed variables (unused in code)
    @mat_dtype          - return arrays in same dtype as loaded into matlab
                          (instead of the dtype with which they are saved)
    @squeeze_me         - whether to squeeze unit dimensions or not
    @chars_as_strings   - whether to convert char arrays to string arrays

    To make this class functional, you will need to override the
    following methods:

    set_dtypes              - sets data types defs from byte order
    matrix_getter_factory   - gives object to fetch next matrix from stream
    format_looks_right      - returns True if format looks correct for
                              this file type (Mat4, Mat5)
    guess_byte_order        - guesses file byte order from file
    """

    def __init__(self, mat_stream,
                 byte_order=None,
                 base_name='raw',
                 mat_dtype=False,
                 squeeze_me=True,
                 chars_as_strings=True,
                 ):
        # Initialize stream
        self.mat_stream = mat_stream
        self.dtypes = {}
        if not byte_order:
            byte_order = self.guess_byte_order()
        self.order_code = byte_order # sets dtypes and other things too
        self.base_name = base_name
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
                                'get/set squeeze me property')
    
    def get_order_code(self):
        return self._order_code
    def set_order_code(self, order_code):
        order_code = ByteOrder.to_numpy_code(order_code)
        self._order_code = order_code
        self.set_dtypes()
    order_code = property(get_order_code,
                          set_order_code,
                          None,
                          'get/set order code')

    def set_dtypes(self):
        assert False, 'Not implemented'

    def convert_dtypes(self, dtype_template):
        dtypes = dtype_template.copy()
        for k in dtypes:
            dtypes[k] = dtype(dtypes[k]).newbyteorder(
                self.order_code)
        return dtypes
    
    def matrix_getter_factory(self):
        assert False, 'Not implemented'
    
    def format_looks_right(self):
        "Return True if the format looks right for this object"
        assert False, 'Not implemented'

    def file_header(self):
        return {}
    
    def guess_byte_order(self):
        assert 0, 'Not implemented'

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
                    dtt = self.order_code + 'U'
                    n_dims = dims[:-1]
                    str_arr = reshape(arr,
                                    (small_product(n_dims),
                                     dims[-1]))
                    arr = empty(n_dims, dtype=object)
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
                arr = squeeze(arr)
                if not arr.size:
                    arr = array([])
                elif not arr.shape: # 0d coverted to scalar
                    arr = arr.item()
            return arr
        return func

    def chars_to_str(self, str_arr):
        ''' Convert string array to string '''
        dt = dtype('U' + str(small_product(str_arr.shape)))
        return ndarray(shape=(),
                       dtype = dt,
                       buffer = str_arr.copy()).item()

    def get_variables(self, variable_names=None):
        ''' get variables from stream as dictionary

        @variable_names   - optional list of variable names to get

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
            res = getter.get_array()
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
        self.mat_stream.seek(-1,1)
        return len(b) == 0


class MatMatrixGetter(MatStreamAgent):
    """ Base class for matrix getters

    Getters are stateful versions of agents, and record state of
    current read on initialization, so need to be created for each
    read - one-shot objects.

    MatrixGetters are initialized with the content of the matrix
    header

    Accepts
    @array_reader - array reading object (see below)
    @header       - header dictionary for matrix being read
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
    def __init__(self, file_stream, arr, name):
        self.file_stream = file_stream
        self.arr = arr
        dt = self.arr.dtype
        if not dt.isnative:
            self.arr = self.arr.astype(dt.newbyteorder('='))
        self.name = name
        
    def arr_dtype_number(self, num):
        ''' Return dtype for given number of items per element'''
        return dtype(self.arr.dtype.str[:2] + str(num))

    def arr_to_chars(self):
        ''' Convert string array to char array '''
        dims = list(self.arr.shape)
        if not dims:
            dims = [1]
        dims.append(int(self.arr.dtype.str[2:]))
        self.arr = ndarray(shape=dims,
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


class DownCaster(object):
    ''' Downcasts arrays '''

    def __init__(self,
                 type_list=None,
                 rtol=1.0000000000000001e-05,
                 atol=1e-08):
        ''' Set types for which we are attempting to downcast '''
        def_dict = self.default_dt_dict()
        if type_list is None:
            self.dt_dict = def_dict
        else:
            dt_dict = {}
            for T in type_list:
                T = dtype(T).type
                dt_dict[T] = def_dict[T]
            self.dt_dict = dt_dict
        self.rtol = rtol
        self.atol = atol

    def eps(self, dt):
        ''' Calculate machine precision for datatype

        Machine precision defined as difference between X and smallest
        encodable number greater than X, where X is usually 1.

        Input can be datatype, in which case X=1, or X.
        '''
        try:
            dt = dtype(dt)
            start = array(1, dt)
        except TypeError:
            start = array(dt)
            dt = start.dtype
        two = array(2, dt)
        e = start.copy()
        while (e / two + start) > start:
            e = e / two
        return e
    
    def default_dt_dict(self):
        d_dict = {}
        for sc_type in ('complex','float'):
            t_list = sctypes[sc_type]
            for T in t_list:
                dt = dtype(T)
                d_dict[T] = {
                    'kind': dt.kind,
                    'size': dt.itemsize}
        for T in sctypes['int']:
            dt = dtype(T)
            sz = dt.itemsize
            bits = sz*8-1
            end = 2**bits
            d_dict[T] = {
                'kind': dt.kind,
                'size': sz,
                'min': -end,
                'max': end-1
                }
        for T in sctypes['uint']:
            dt = dtype(T)
            sz = dt.itemsize
            bits = sz*8
            end = 2**bits
            d_dict[T] = {
                'kind': dt.kind,
                'size': sz,
                'min': 0,
                'max': end
            }
        return d_dict
    
    def storage_criterion(self, maxstorage, kinds, cmp_func=lambda x, y: x <= y):
        D = {}
        for k, v in self.dt_dict.items():
            if v['kind'] in kinds:
                sz = v['size']
                if cmp_func(sz, maxstorage):
                    D[k] = sz
        I = D.items()
        I.sort(lambda x, y: cmp(x[1], y[1]))
        return I

    def smaller_same_kind(self, arr):
        dts = self.storage_criterion(arr.dtype.itemsize,
                                      (arr.dtype.kind,),
                                      lambda x, y: x < y)
        ret_arr = arr
        for T in dts:
            test_arr = arr.astype(T)
            if allclose(test_arr, arr, self.rtol, self.atol):
                ret_arr = test_arr
            else:
                break
        return ret_arr

        
    def smallest_int_type(self, mx, mn):
        dt = None
        for k, v in self.dt_dict.items():
            if v['kind'] in ('i', 'u'):
                if v['max'] >= mx and v['min'] <= mn:
                    c_sz = v['size']
                    if dt is None or c_sz < sz:
                        dt = k
                        sz = c_sz
        return dt

    def downcast(self, arr):
        dtk = arr.dtype.kind
        if dtk == 'c':
            return self.downcast_complex(arr)
        elif dtk == 'f':
            return self.downcast_float(arr)
        elif dtk in ('u', 'i'):
            return self.downcast_integer(arr)
        else:
            raise TypeError, 'Do not recognize array kind %s' % dtk
            
    def downcast_complex(self, arr):
        # can we downcast to float?
        fts = self.dt_arrs['float']
        flts = flts[flts['storage'] <= arr.dtype.itemsize]
        test_arr = arr.astype(flt[0]['type'])
        if allclose(arr, test_arr, self.rtol, self.atol):
            return self.downcast_float(test_arr)
        # try downcasting to another complex type
        return self.smaller_same_kind(arr)
    
    def downcast_float(self, arr):
        # Try integer
        test_arr = self.downcast_integer(arr)
        if allclose(arr, test_arr, self.rtol, self.atol):
            return test_arr
        # Otherwise descend the float types
        return self.smaller_same_kind(arr)

    def downcast_integer(self, arr):
        mx = amax(arr)
        mn = amin(arr)
        idt = self.smallest_int_type(mx, mn)
        return arr.astype(idt)
