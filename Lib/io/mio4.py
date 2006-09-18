''' Classes for read / write of matlab 4 files
'''

from numpy import *

from miobase import *

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

class Mat4Header(object):
    ''' Place holder for Mat4 header

        Defines:
        next_position - start position of next matrix
        name
        dims - shape of matrix as stored (see sparse reader)
        dtype - numpy dtype of matrix
        mclass - matlab code for class of matrix
        is_char    - True if these are char data
        is_numeric - True if these are numeric data
        is_complex - True if data are complex
        original_dtype - data type in matlab workspace
    '''
    def __init__(self):
        self.next_position = None
        self.name = ''
        self.dims = ()
        self.dtype = None
        self.mclass = None
        self.is_char = None
        self.is_numeric = None
        self.is_complex = None
        self.original_dtype = None
        

class Mat4ArrayReader(MatArrayReader):
    ''' Class for reading Mat4 arrays
    '''
    
    def __init__(self, *args, **kwargs):
        super(Mat4ArrayReader,self).__init__(*args, **kwargs)
        self._getter_classes = {
            mxFULL_CLASS: Mat4FullGetter,
            mxCHAR_CLASS: Mat4CharGetter,
            mxSPARSE_CLASS: Mat4SparseGetter,
            }
        
    def read_header(self):
        ''' Read and return Mat4 matrix header
        '''
        header = Mat4Header()
        data = self.read_array(self.dtypes['header'])
        header.name = self.read_ztstring(data['namlen'])
        if data['mopt'] < 0 or  data['mopt'] > 5000:
            ValueError, 'Mat 4 mopt wrong format, byteswapping problem?'
        M,rest = divmod(data['mopt'], 1000)
        O,rest = divmod(rest,100)
        P,rest = divmod(rest,10)
        T = rest
        if O != 0:
            raise ValueError, 'O in MOPT integer should be 0, wrong format?'
        header.dtype = self.dtypes[P]
        header.mclass = T
        header.dims = (data['mrows'], data['ncols'])
        header.is_complex = data['imagf'] == 1
        remaining_bytes = header.dtype.itemsize * product(header.dims)
        if header.is_complex and not header.mclass == mxSPARSE_CLASS:
            remaining_bytes *= 2
        header.next_position = self.mat_stream.tell() + remaining_bytes
        return header

    def matrix_getter_factory(self):
        header = self.read_header()
        return self._getter_classes[header.mclass](self, header)


class Mat4MatrixGetter(MatMatrixGetter):

    # Mat4 variables never global or logical
    is_global = False
    is_logical = False
    
    def read_hdr_array(self, *args, **kwargs):
        ''' Mat4 read array always uses header dtype and dims '''
        return self.read_array(
            self.header.dtype, self.dims, *args, **kwargs)


class Mat4FullGetter(Mat4MatrixGetter):
    def get_raw_array(self):
        self.header.is_numeric = True
        self.header.original_dtype = dtype(float64)
        if self.header.is_complex:
            # avoid array copy to save memory
            res = self.read_hdr_array(copy=False)
            res_j = self.read_hdr_array(copy=False)
            return res + (res_j * 1j)
        else:
            return self.read_hdr_array()


class Mat4CharGetter(Mat4MatrixGetter):
    def get_raw_array(self):
        self.header.is_char = True
        arr = self.read_hdr_array().astype(uint8)
        # ascii to unicode
        S = arr.tostring().decode('ascii')
        return ndarray(shape=self.dims,
                       dtype=dtype('U1'),
                       buffer = array(S)).copy()


class Mat4SparseGetter(Mat4MatrixGetter):
    ''' Read sparse matrix type 

    Matlab 4 real sparse arrays are saved in a N+1 by 3 array format,
    where N is the number of non-zero values.  Column 1 values [0:N]
    are the (1-based) row indices of the each non-zero value, column 2
    [0:N] are the column indices, column 3 [0:N] are the (real)
    values.  The last values [-1:0:2] of the rows, column indices are
    shape[0] and shape[1] respectively of the output matrix. The last
    value for the values column is a padding 0. mrows and ncols values
    from the header give the shape of the stored matrix, here [N+1,
    3].  Complex data is saved as a 4 column matrix, where the fourth
    column contains the imaginary component of the data; the last
    value is again 0
    '''
    def get_raw_array(self):
        self.header.original_dtype = dtype(float64)
        res = self.read_hdr_array()
        tmp = res[:-1,:]
        dims = res[-1,0:2]
        ij = transpose(tmp[:,0:2]) - 1 # for matlab 1-based indexing
        vals = tmp[:,2]
        if res.shape[1] == 4:
            self.header.is_complex = True
            vals = vals + res[:-1,3] * 1j
        if have_sparse:
            return scipy.sparse.csc_matrix((vals,ij), dims)
        return (dims, ij, vals)

    
class MatFile4Reader(MatFileReader):
    ''' Reader for Mat4 files '''
    def __init__(self, mat_stream, *args, **kwargs):
        self._array_reader = Mat4ArrayReader(
            mat_stream,
            None,
            None,
            )
        super(MatFile4Reader, self).__init__(mat_stream, *args, **kwargs)
        self._array_reader.processor_func = self.processor_func

    def set_dtypes(self):
        self.dtypes = self.convert_dtypes(mdtypes_template)
        self._array_reader.dtypes = self.dtypes

    def matrix_getter_factory(self):
        return self._array_reader.matrix_getter_factory()

    def format_looks_right(self):
        # Matlab 4 files have a zero somewhere in first 4 bytes
        self.mat_stream.seek(0)
        mopt_bytes = self.read_bytes(4)
        self.mat_stream.seek(0)
        return 0 in mopt_bytes
    
    def guess_byte_order(self):
        self.mat_stream.seek(0)
        mopt = self.read_array(dtype('i4'))
        self.mat_stream.seek(0)
        if mopt < 0 or mopt > 5000:
            return ByteOrder.swapped_code
        return ByteOrder.native_code


class MatFile4Writer(MatFileWriter):
    codec = 'ascii'
    
    def arr_to_matrix(self, arr):
        ''' Convert numeric array to matlab format '''
        dts = arr.dtype.str[1:]
        if not dts in np_to_mtypes:
            arr = arr.astype('f8')
        return atleast_2d(arr)
        
    def matrix_header(self, var, name):
        ''' Return header for matrix array '''
        header = empty((), mdtypes_template['header'])
        dt = var.dtype.str[1:]
        M = not ByteOrder.little_endian
        O = 0
        P = np_to_mtypes[dt]
        T = dt == 'S1' # could also be sparse -> 2
        header['mopt'] = M*1000+O*100+P*10+T
        dims = var.shape
        header['mrows'] = dims[0]
        header['ncols'] = dims[1]
        header['imagf'] = var.dtype.kind == 'c'
        header['namlen'] = len(name) + 1
        return header
    
    def put_variable(self, var, name):
        arr = array(var)
        if arr.dtype.hasobject:
            raise TypeError, 'Cannot save object arrays in Mat4'
        if have_sparse:
            if scipy.sparse.issparse(arr):
                raise TypeError, 'Cannot save sparse arrays yet'
        if arr.dtype.kind in ('U', 'S'):
            arr = self.str_to_chars(arr)
        else:
            arr = self.arr_to_matrix(arr)
        dims = arr.shape
        if len(dims) > 2:
            dims = [product(dims[:-1]), dims[-1]]
            arr = reshape(arr, dims)
        if arr.dtype.kind == 'U':
            # Recode unicode to ascii
            dt = 'U' + str(product(dims))
            st_arr = ndarray(shape=(), dtype=dt, buffer=arr)
            st = st_arr.item().encode('ascii')
            arr = ndarray(shape=dims, dtype='S1', buffer=st)
        header = self.matrix_header(arr, name)
        self.write_bytes(header)
        self.write_string(name + '\0')
        if header['imagf']:
            self.write_bytes(arr.real)
            self.write_bytes(arr.imag)
        else:
            self.write_bytes(arr)
            
    def put_variables(self, mdict):
        for name, var in mdict.items():
            self.put_variable(var, name)
