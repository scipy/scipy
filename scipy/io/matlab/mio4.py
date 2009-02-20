''' Classes for read / write of matlab (TM) 4 files
'''
import sys
import warnings

import numpy as np

import scipy.sparse

from miobase import MatFileReader, MatArrayReader, MatMatrixGetter, \
     MatFileWriter, MatStreamWriter, docfiller, matdims


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


class Mat4ArrayReader(MatArrayReader):
    ''' Class for reading Mat4 arrays
    '''

    def matrix_getter_factory(self):
        ''' Read header, return matrix getter '''
        data = self.read_dtype(self.dtypes['header'])
        header = {}
        header['name'] = self.read_ztstring(int(data['namlen']))
        if data['mopt'] < 0 or  data['mopt'] > 5000:
            ValueError, 'Mat 4 mopt wrong format, byteswapping problem?'
        M,rest = divmod(data['mopt'], 1000)
        O,rest = divmod(rest,100)
        P,rest = divmod(rest,10)
        T = rest
        if O != 0:
            raise ValueError, 'O in MOPT integer should be 0, wrong format?'
        header['dtype'] = self.dtypes[P]
        header['mclass'] = T
        header['dims'] = (data['mrows'], data['ncols'])
        header['is_complex'] = data['imagf'] == 1
        remaining_bytes = header['dtype'].itemsize * np.product(header['dims'])
        if header['is_complex'] and not header['mclass'] == mxSPARSE_CLASS:
            remaining_bytes *= 2
        next_pos = self.mat_stream.tell() + remaining_bytes
        if T == mxFULL_CLASS:
            getter = Mat4FullGetter(self, header)
        elif T == mxCHAR_CLASS:
            getter = Mat4CharGetter(self, header)
        elif T == mxSPARSE_CLASS:
            getter = Mat4SparseGetter(self, header)
        else:
            raise TypeError, 'No reader for class code %s' % T
        getter.next_position = next_pos
        return getter


class Mat4MatrixGetter(MatMatrixGetter):

    # Mat4 variables never global or logical
    is_global = False
    is_logical = False

    def read_array(self, copy=True):
        ''' Mat4 read array always uses header dtype and dims
        copy        - copies array if True
        (buffer is usually read only)
        a_dtype is assumed to be correct endianness
        '''
        dt = self.header['dtype']
        dims = self.header['dims']
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


class Mat4FullGetter(Mat4MatrixGetter):
    def __init__(self, array_reader, header):
        super(Mat4FullGetter, self).__init__(array_reader, header)
        if header['is_complex']:
            self.mat_dtype = np.dtype(np.complex128)
        else:
            self.mat_dtype = np.dtype(np.float64)

    def get_raw_array(self):
        if self.header['is_complex']:
            # avoid array copy to save memory
            res = self.read_array(copy=False)
            res_j = self.read_array(copy=False)
            return res + (res_j * 1j)
        return self.read_array()


class Mat4CharGetter(Mat4MatrixGetter):
    def get_raw_array(self):
        arr = self.read_array().astype(np.uint8)
        # ascii to unicode
        S = arr.tostring().decode('ascii')
        return np.ndarray(shape=self.header['dims'],
                          dtype=np.dtype('U1'),
                          buffer = np.array(S)).copy()


class Mat4SparseGetter(Mat4MatrixGetter):
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
    def get_raw_array(self):
        res = self.read_array()
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

    def guess_byte_order(self):
        self.mat_stream.seek(0)
        mopt = self.read_dtype(np.dtype('i4'))
        self.mat_stream.seek(0)
        if mopt < 0 or mopt > 5000:
            return SYS_LITTLE_ENDIAN and '>' or '<'
        return SYS_LITTLE_ENDIAN and '<' or '>'


class Mat4MatrixWriter(MatStreamWriter):

    def write_header(self, P=0,  T=0, imagf=0, dims=None):
        ''' Write header for given data options
        P      - mat4 data type
        T      - mat4 matrix class
        imagf  - complex flag
        dims   - matrix dimensions
        '''
        if dims is None:
            dims = self.arr.shape
        header = np.empty((), mdtypes_template['header'])
        M = not SYS_LITTLE_ENDIAN
        O = 0
        header['mopt'] = (M * 1000 +
                          O * 100 +
                          P * 10 +
                          T)
        header['mrows'] = dims[0]
        header['ncols'] = dims[1]
        header['imagf'] = imagf
        header['namlen'] = len(self.name) + 1
        self.write_bytes(header)
        self.write_string(self.name + '\0')

    def arr_to_2d(self):
        dims = matdims(self.arr, self.oned_as)
        self.arr.shape = dims
        if len(dims) > 2:
            self.arr = self.arr.reshape(-1,dims[-1])

    def write(self):
        assert False, 'Not implemented'


class Mat4NumericWriter(Mat4MatrixWriter):

    def write(self):
        self.arr_to_2d()
        imagf = self.arr.dtype.kind == 'c'
        try:
            P = np_to_mtypes[self.arr.dtype.str[1:]]
        except KeyError:
            if imagf:
                self.arr = self.arr.astype('c128')
            else:
                self.arr = self.arr.astype('f8')
            P = miDOUBLE
        self.write_header(P=P,
                          T=mxFULL_CLASS,
                          imagf=imagf)
        if imagf:
            self.write_bytes(self.arr.real)
            self.write_bytes(self.arr.imag)
        else:
            self.write_bytes(self.arr)


class Mat4CharWriter(Mat4MatrixWriter):

    def write(self):
        self.arr_to_chars()
        self.arr_to_2d()
        dims = self.arr.shape
        self.write_header(P=miUINT8,
                          T=mxCHAR_CLASS)
        if self.arr.dtype.kind == 'U':
            # Recode unicode to ascii
            n_chars = np.product(dims)
            st_arr = np.ndarray(shape=(),
                                dtype=self.arr_dtype_number(n_chars),
                                buffer=self.arr)
            st = st_arr.item().encode('ascii')
            self.arr = np.ndarray(shape=dims, dtype='S1', buffer=st)
        self.write_bytes(self.arr)


class Mat4SparseWriter(Mat4MatrixWriter):

    def write(self):
        ''' Sparse matrices are 2D
        See docstring for Mat4SparseGetter
        '''
        A = self.arr.tocoo() #convert to sparse COO format (ijv)
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
        self.write_header(P=miDOUBLE,
                          T=mxSPARSE_CLASS,
                          dims=ijv.shape)
        self.write_bytes(ijv)


def matrix_writer_factory(stream, arr, name, oned_as='row'):
    ''' Factory function to return matrix writer given variable to write
    stream      - file or file-like stream to write to
    arr         - array to write
    name        - name in matlab (TM) workspace
    '''
    if scipy.sparse.issparse(arr):
        return Mat4SparseWriter(stream, arr, name, oned_as)
    arr = np.array(arr)
    dtt = arr.dtype.type
    if dtt is np.object_:
        raise TypeError, 'Cannot save object arrays in Mat4'
    elif dtt is np.void:
        raise TypeError, 'Cannot save void type arrays'
    elif dtt in (np.unicode_, np.string_):
        return Mat4CharWriter(stream, arr, name, oned_as)
    else:
        return Mat4NumericWriter(stream, arr, name, oned_as)


class MatFile4Writer(MatFileWriter):
    ''' Class for writing matlab 4 format files '''
    def __init__(self, file_stream, oned_as=None):
        self.file_stream = file_stream
        if oned_as is None:
            oned_as = 'row'
        self.oned_as = oned_as

    def put_variables(self, mdict):
        for name, var in mdict.items():
            matrix_writer_factory(self.file_stream, 
                                  var, 
                                  name, 
                                  self.oned_as).write()
