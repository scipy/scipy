''' Classes for read / write of matlab (TM) 4 files
'''
import sys
import warnings

import numpy as np

import scipy.sparse

from miobase import MatFileReader, docfiller, matdims, \
     read_dtype, convert_dtypes, arr_to_chars, arr_dtype_number, \
     MatWriteError

from mio_utils import squeeze_element, chars_to_strings


SYS_LITTLE_ENDIAN = sys.byteorder == 'little'

miDOUBLE = 0
miSINGLE = 1
miINT32 = 2
miINT16 = 3
miUINT16 = 4
miUINT8 = 5

mdtypes_template = {
    miDOUBLE: 'f8',
    miSINGLE: 'f4',
    miINT32: 'i4',
    miINT16: 'i2',
    miUINT16: 'u2',
    miUINT8: 'u1',
    'header': [('mopt', 'i4'),
               ('mrows', 'i4'),
               ('ncols', 'i4'),
               ('imagf', 'i4'),
               ('namlen', 'i4')],
    'U1': 'U1',
    }

np_to_mtypes = {
    'f8': miDOUBLE,
    'c32': miDOUBLE,
    'c24': miDOUBLE,
    'c16': miDOUBLE,
    'f4': miSINGLE,
    'c8': miSINGLE,
    'i4': miINT32,
    'i2': miINT16,
    'u2': miUINT16,
    'u1': miUINT8,
    'S1': miUINT8,
    }

# matrix classes
mxFULL_CLASS = 0
mxCHAR_CLASS = 1
mxSPARSE_CLASS = 2

order_codes = {
    0: '<',
    1: '>',
    2: 'VAX D-float', #!
    3: 'VAX G-float',
    4: 'Cray', #!!
    }

class VarHeader4(object):
    # Mat4 variables never logical or global
    is_logical = False
    is_global = False

    def __init__(self,
                 name,
                 dtype,
                 mclass,
                 dims,
                 is_complex):
        self.name = name
        self.dtype = dtype
        self.mclass = mclass
        self.dims = dims
        self.is_complex = is_complex
    

class VarReader4(object):
    ''' Class to read matlab 4 variables '''
    
    def __init__(self, file_reader):
        self.file_reader = file_reader
        self.mat_stream = file_reader.mat_stream
        self.dtypes = file_reader.dtypes
        self.chars_as_strings = file_reader.chars_as_strings
        self.squeeze_me = file_reader.squeeze_me
        
    def read_header(self):
        ''' Reads and return header for variable '''
        data = read_dtype(self.mat_stream, self.dtypes['header'])
        name = self.mat_stream.read(int(data['namlen'])).strip('\x00')
        if data['mopt'] < 0 or  data['mopt'] > 5000:
            ValueError, 'Mat 4 mopt wrong format, byteswapping problem?'
        M,rest = divmod(data['mopt'], 1000)
        O,rest = divmod(rest,100)
        P,rest = divmod(rest,10)
        T = rest
        if O != 0:
            raise ValueError, 'O in MOPT integer should be 0, wrong format?'
        dims = (data['mrows'], data['ncols'])
        is_complex = data['imagf'] == 1
        dtype = self.dtypes[P]
        return VarHeader4(
            name,
            dtype,
            T,
            dims,
            is_complex)
        
    def array_from_header(self, hdr, process=True):
        mclass = hdr.mclass
        if mclass == mxFULL_CLASS:
            arr = self.read_full_array(hdr)
        elif mclass == mxCHAR_CLASS:
            arr = self.read_char_array(hdr)
            if process and self.chars_as_strings:
                arr = chars_to_strings(arr)
        elif mclass == mxSPARSE_CLASS:
            # no current processing (below) makes sense for sparse
            return self.read_sparse_array(hdr)
        else:
            raise TypeError, 'No reader for class code %s' % mclass
        if process and self.squeeze_me:
            return squeeze_element(arr)
        return arr
    
    def read_sub_array(self, hdr, copy=True):
        ''' Mat4 read always uses header dtype and dims
        hdr : object
           object with attributes 'dtype', 'dims'
        copy : bool
           copies array if True (default True)
           (buffer is usually read only)

        self.dtype is assumed to be correct endianness
        '''
        dt = hdr.dtype
        dims = hdr.dims
        num_bytes = dt.itemsize
        for d in dims:
            num_bytes *= d
        arr = np.ndarray(shape=dims,
                         dtype=dt,
                         buffer=self.mat_stream.read(num_bytes),
                         order='F')
        if copy:
            arr = arr.copy()
        return arr

    def read_full_array(self, hdr):
        ''' Full (rather than sparse matrix) getter
        '''
        if hdr.is_complex:
            # avoid array copy to save memory
            res = self.read_sub_array(hdr, copy=False)
            res_j = self.read_sub_array(hdr, copy=False)
            return res + (res_j * 1j)
        return self.read_sub_array(hdr)

    def read_char_array(self, hdr):
        ''' Ascii text matrix (char matrix) reader

        '''
        arr = self.read_sub_array(hdr).astype(np.uint8)
        # ascii to unicode
        S = arr.tostring().decode('ascii')
        return np.ndarray(shape=hdr.dims,
                          dtype=np.dtype('U1'),
                          buffer = np.array(S)).copy()

    def read_sparse_array(self, hdr):
        ''' Read sparse matrix type

        Matlab (TM) 4 real sparse arrays are saved in a N+1 by 3 array
        format, where N is the number of non-zero values.  Column 1 values
        [0:N] are the (1-based) row indices of the each non-zero value,
        column 2 [0:N] are the column indices, column 3 [0:N] are the
        (real) values.  The last values [-1,0:2] of the rows, column
        indices are shape[0] and shape[1] respectively of the output
        matrix. The last value for the values column is a padding 0. mrows
        and ncols values from the header give the shape of the stored
        matrix, here [N+1, 3].  Complex data is saved as a 4 column
        matrix, where the fourth column contains the imaginary component;
        the last value is again 0.  Complex sparse data do _not_ have the
        header imagf field set to True; the fact that the data are complex
        is only detectable because there are 4 storage columns
        '''
        res = self.read_sub_array(hdr)
        tmp = res[:-1,:]
        dims = res[-1,0:2]
        I = np.ascontiguousarray(tmp[:,0],dtype='intc') #fixes byte order also
        J = np.ascontiguousarray(tmp[:,1],dtype='intc')
        I -= 1  # for 1-based indexing
        J -= 1
        if res.shape[1] == 3:
            V = np.ascontiguousarray(tmp[:,2],dtype='float')
        else:
            V = np.ascontiguousarray(tmp[:,2],dtype='complex')
            V.imag = tmp[:,3]
        return scipy.sparse.coo_matrix((V,(I,J)), dims)


class MatFile4Reader(MatFileReader):
    ''' Reader for Mat4 files '''
    @docfiller
    def __init__(self, mat_stream, *args, **kwargs):
        ''' Initialize matlab 4 file reader

    %(matstream_arg)s
    %(load_args)s
        '''
        super(MatFile4Reader, self).__init__(mat_stream, *args, **kwargs)
        self._matrix_reader = None
        
    def guess_byte_order(self):
        self.mat_stream.seek(0)
        mopt = read_dtype(self.mat_stream, np.dtype('i4'))
        self.mat_stream.seek(0)
        if mopt < 0 or mopt > 5000:
            return SYS_LITTLE_ENDIAN and '>' or '<'
        return SYS_LITTLE_ENDIAN and '<' or '>'

    def initialize_read(self):
        ''' Run when beginning read of variables

        Sets up readers from parameters in `self`
        '''
        self.dtypes = convert_dtypes(mdtypes_template, self.byte_order)
        self._matrix_reader = VarReader4(self)
    
    def read_var_header(self):
        ''' Read header, return header, next position

        Header has to define at least .name and .is_global

        Parameters
        ----------
        None

        Returns
        -------
        header : object
           object that can be passed to self.read_var_array, and that
           has attributes .name and .is_global
        next_position : int
           position in stream of next variable
        '''
        hdr = self._matrix_reader.read_header()
        n = reduce(lambda x, y: x*y, hdr.dims, 1) # fast product
        remaining_bytes = hdr.dtype.itemsize * n
        if hdr.is_complex and not hdr.mclass == mxSPARSE_CLASS:
            remaining_bytes *= 2
        next_position = self.mat_stream.tell() + remaining_bytes
        return hdr, next_position
    
    def read_var_array(self, header, process=True):
        ''' Read array, given `header`

        Parameters
        ----------
        header : header object
           object with fields defining variable header
        process : {True, False} bool, optional
           If True, apply recursive post-processing during loading of
           array. 
        
        Returns
        -------
        arr : array
           array with post-processing applied or not according to
           `process`. 
        '''
        return self._matrix_reader.array_from_header(header, process)

    def get_variables(self, variable_names=None):
        ''' get variables from stream as dictionary

        variable_names   - optional list of variable names to get

        If variable_names is None, then get all variables in file
        '''
        if isinstance(variable_names, basestring):
            variable_names = [variable_names]
        self.mat_stream.seek(0)
        # set up variable reader
        self.initialize_read()
        mdict = {}
        while not self.end_of_stream():
            hdr, next_position = self.read_var_header()
            name = hdr.name
            if variable_names and name not in variable_names:
                self.mat_stream.seek(next_position)
                continue
            mdict[name] = self.read_var_array(hdr)
            self.mat_stream.seek(next_position)
            if variable_names:
                variable_names.remove(name)
                if len(variable_names) == 0:
                    break
        return mdict


def arr_to_2d(arr, oned_as='row'):
    ''' Make ``arr`` exactly two dimensional

    If `arr` has more than 2 dimensions, then, for the sake of
    compatibility with previous versions of scipy, we reshape to 2D
    preserving the last dimension and increasing the first dimension.
    In future versions we will raise an error, as this is at best a very
    counterinituitive thing to do.

    Parameters
    ----------
    arr : array
    oned_as : {'row', 'column'}
       Whether to reshape 1D vectors as row vectors or column vectors.
       See documentation for ``matdims`` for more detail

    Returns
    -------
    arr2d : array
       2D version of the array
    '''
    dims = matdims(arr, oned_as)
    if len(dims) > 2:
        warnings.warn('Matlab 4 files only support <=2 '
                      'dimensions; the next version of scipy will '
                      'raise an error when trying to write >2D arrays '
                      'to matlab 4 format files',
                      DeprecationWarning,
                      )
        return arr.reshape((-1,dims[-1]))
    return arr.reshape(dims)


class VarWriter4(object):
    def __init__(self, file_writer):
        self.file_stream = file_writer.file_stream
        self.oned_as = file_writer.oned_as

    def write_bytes(self, arr):
        self.file_stream.write(arr.tostring(order='F'))

    def write_string(self, s):
        self.file_stream.write(s)

    def write_header(self, name, shape, P=0,  T=0, imagf=0):
        ''' Write header for given data options

        Parameters
        ----------
        name : str
        shape : sequence
           Shape of array as it will be read in matlab
        P      - mat4 data type
        T      - mat4 matrix class
        imagf  - complex flag
        '''
        header = np.empty((), mdtypes_template['header'])
        M = not SYS_LITTLE_ENDIAN
        O = 0
        header['mopt'] = (M * 1000 +
                          O * 100 +
                          P * 10 +
                          T)
        header['mrows'] = shape[0]
        header['ncols'] = shape[1]
        header['imagf'] = imagf
        header['namlen'] = len(name) + 1
        self.write_bytes(header)
        self.write_string(name + '\0')

    def write(self, arr, name):
        ''' Write matrix `arr`, with name `name`

        Parameters
        ----------
        arr : array-like
           array to write
        name : str
           name in matlab workspace
        '''
        # we need to catch sparse first, because np.asarray returns an
        # an object array for scipy.sparse
        if scipy.sparse.issparse(arr):
            self.write_sparse(arr, name)
            return
        arr = np.asarray(arr)
        dt = arr.dtype
        if not dt.isnative:
            arr = arr.astype(dt.newbyteorder('='))
        dtt = dt.type
        if dtt is np.object_:
            raise TypeError, 'Cannot save object arrays in Mat4'
        elif dtt is np.void:
            raise TypeError, 'Cannot save void type arrays'
        elif dtt in (np.unicode_, np.string_):
            self.write_char(arr, name)
            return
        self.write_numeric(arr, name)

    def write_numeric(self, arr, name):
        arr = arr_to_2d(arr, self.oned_as)
        imagf = arr.dtype.kind == 'c'
        try:
            P = np_to_mtypes[arr.dtype.str[1:]]
        except KeyError:
            if imagf:
                arr = arr.astype('c128')
            else:
                arr = arr.astype('f8')
            P = miDOUBLE
        self.write_header(name,
                          arr.shape,
                          P=P,
                          T=mxFULL_CLASS,
                          imagf=imagf)
        if imagf:
            self.write_bytes(arr.real)
            self.write_bytes(arr.imag)
        else:
            self.write_bytes(arr)

    def write_char(self, arr, name):
        arr = arr_to_chars(arr)
        arr = arr_to_2d(arr, self.oned_as)
        dims = arr.shape
        self.write_header(
            name,
            dims,
            P=miUINT8,
            T=mxCHAR_CLASS)
        if arr.dtype.kind == 'U':
            # Recode unicode to ascii
            n_chars = np.product(dims)
            st_arr = np.ndarray(shape=(),
                                dtype=arr_dtype_number(arr, n_chars),
                                buffer=arr)
            st = st_arr.item().encode('ascii')
            arr = np.ndarray(shape=dims, dtype='S1', buffer=st)
        self.write_bytes(arr)

    def write_sparse(self, arr, name):
        ''' Sparse matrices are 2D
        
        See docstring for VarReader4.read_sparse_array
        '''
        A = arr.tocoo() #convert to sparse COO format (ijv)
        imagf = A.dtype.kind == 'c'
        ijv = np.zeros((A.nnz + 1, 3+imagf), dtype='f8')
        ijv[:-1,0] = A.row
        ijv[:-1,1] = A.col
        ijv[:-1,0:2] += 1 # 1 based indexing
        if imagf:
            ijv[:-1,2] = A.data.real
            ijv[:-1,3] = A.data.imag
        else:
            ijv[:-1,2] = A.data
        ijv[-1,0:2] = A.shape
        self.write_header(
            name,
            ijv.shape,
            P=miDOUBLE,
            T=mxSPARSE_CLASS)
        self.write_bytes(ijv)


class MatFile4Writer(object):
    ''' Class for writing matlab 4 format files '''
    def __init__(self, file_stream, oned_as=None):
        self.file_stream = file_stream
        if oned_as is None:
            oned_as = 'row'
        self.oned_as = oned_as
        self._matrix_writer = None

    def put_variables(self, mdict, write_header=None):
        ''' Write variables in `mdict` to stream

        Parameters
        ----------
        mdict : mapping
           mapping with method ``items`` return name, contents pairs
           where ``name`` which will appeak in the matlab workspace in
           file load, and ``contents`` is something writeable to a
           matlab file, such as a numpy array.
        write_header : {None, True, False}
           If True, then write the matlab file header before writing the
           variables.  If None (the default) then write the file header
           if we are at position 0 in the stream.  By setting False
           here, and setting the stream position to the end of the file,
           you can append variables to a matlab file
        '''
        # there is no header for a matlab 4 mat file, so we ignore the
        # ``write_header`` input argument.  It's there for compatibility
        # with the matlab 5 version of this method
        self._matrix_writer = VarWriter4(self)
        for name, var in mdict.items():
            self._matrix_writer.write(var, name)
