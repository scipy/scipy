''' Classes for read / write of matlab (TM) 5 files
'''

# Small fragments of current code adapted from matfile.py by Heiko
# Henkelmann

import sys
import zlib
from cStringIO import StringIO
from copy import copy as pycopy

import numpy as np

from miobase import MatFileReader, MatArrayReader, MatMatrixGetter, \
     MatFileWriter, MatStreamWriter, spsparse

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

mxCELL_CLASS = 1
mxSTRUCT_CLASS = 2
# The March 2008 edition of "Matlab 7 MAT-File Format" says that
# mxOBJECT_CLASS = 3, whereas matrix.h says that mxLOGICAL = 3.
# Matlab 2008a appears to save logicals as type 9, so we assume that
# the document is correct.  See type 18, below.
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
# The following are not in the March 2008 edition of "Matlab 7
# MAT-File Format," but were guessed from matrix.h.
mxINT64_CLASS = 14
mxUINT64_CLASS = 15
mxFUNCTION_CLASS = 16
# Not doing anything with these at the moment.
mxOPAQUE_CLASS = 17
mxOBJECT_CLASS_FROM_MATRIX_H = 18

mdtypes_template = {
    miINT8: 'i1',
    miUINT8: 'u1',
    miINT16: 'i2',
    miUINT16: 'u2',
    miINT32: 'i4',
    miUINT32: 'u4',
    miSINGLE: 'f4',
    miDOUBLE: 'f8',
    miINT64: 'i8',
    miUINT64: 'u8',
    miUTF8: 'u1',
    miUTF16: 'u2',
    miUTF32: 'u4',
    'file_header': [('description', 'S116'),
                    ('subsystem_offset', 'i8'),
                    ('version', 'u2'),
                    ('endian_test', 'S2')],
    'tag_full': [('mdtype', 'u4'), ('byte_count', 'u4')],
    'tag_smalldata':[('byte_count_mdtype', 'u4'), ('data', 'S4')],
    'array_flags': [('data_type', 'u4'),
                    ('byte_count', 'u4'),
                    ('flags_class','u4'),
                    ('nzmax', 'u4')],
    'U1': 'U1',
    }

mclass_dtypes_template = {
    mxINT8_CLASS: 'i1',
    mxUINT8_CLASS: 'u1',
    mxINT16_CLASS: 'i2',
    mxUINT16_CLASS: 'u2',
    mxINT32_CLASS: 'i4',
    mxUINT32_CLASS: 'u4',
    mxINT64_CLASS: 'i8',
    mxUINT64_CLASS: 'u8',
    mxSINGLE_CLASS: 'f4',
    mxDOUBLE_CLASS: 'f8',
    }


np_to_mtypes = {
    'f8': miDOUBLE,
    'c32': miDOUBLE,
    'c24': miDOUBLE,
    'c16': miDOUBLE,
    'f4': miSINGLE,
    'c8': miSINGLE,
    'i1': miINT8,
    'i2': miINT16,
    'i4': miINT32,
    'i8': miINT64,
    'u1': miUINT8,
    'u2': miUINT16,
    'u4': miUINT32,
    'u8': miUINT64,
    'S1': miUINT8,
    'U1': miUTF16,
    }


np_to_mxtypes = {
    'f8': mxDOUBLE_CLASS,
    'c32': mxDOUBLE_CLASS,
    'c24': mxDOUBLE_CLASS,
    'c16': mxDOUBLE_CLASS,
    'f4': mxSINGLE_CLASS,
    'c8': mxSINGLE_CLASS,
    'i8': mxINT64_CLASS,
    'i4': mxINT32_CLASS,
    'i2': mxINT16_CLASS,
    'u8': mxUINT64_CLASS,
    'u2': mxUINT16_CLASS,
    'u1': mxUINT8_CLASS,
    'S1': mxUINT8_CLASS,
    }



''' Before release v7.1 (release 14) matlab (TM) used the system
default character encoding scheme padded out to 16-bits. Release 14
and later use Unicode. When saving character data, R14 checks if it
can be encoded in 7-bit ascii, and saves in that format if so.'''

codecs_template = {
    miUTF8: {'codec': 'utf_8', 'width': 1},
    miUTF16: {'codec': 'utf_16', 'width': 2},
    miUTF32: {'codec': 'utf_32','width': 4},
    }

miUINT16_codec = sys.getdefaultencoding()

mx_numbers = (
    mxDOUBLE_CLASS,
    mxSINGLE_CLASS,
    mxINT8_CLASS,
    mxUINT8_CLASS,
    mxINT16_CLASS,
    mxUINT16_CLASS,
    mxINT32_CLASS,
    mxUINT32_CLASS,
    mxINT64_CLASS,
    mxUINT64_CLASS,
    )

class Mat5ArrayReader(MatArrayReader):
    ''' Class to get Mat5 arrays

    Provides element reader functions, header reader, matrix reader
    factory function
    '''

    def __init__(self, mat_stream, dtypes, processor_func, codecs, class_dtypes, struct_as_record):
        super(Mat5ArrayReader, self).__init__(mat_stream,
                                              dtypes,
                                              processor_func)
        self.codecs = codecs
        self.class_dtypes = class_dtypes
        self.struct_as_record = struct_as_record

    def read_element(self, copy=True):
        raw_tag = self.mat_stream.read(8)
        tag = np.ndarray(shape=(),
                         dtype=self.dtypes['tag_full'],
                         buffer=raw_tag)
        mdtype = tag['mdtype'].item()
        # Byte count if this is small data element
        byte_count = mdtype >> 16
        if byte_count: # small data element format
            if byte_count > 4:
                raise ValueError, 'Too many bytes for sde format'
            mdtype = mdtype & 0xFFFF
            if mdtype == miMATRIX:
                raise TypeError('Cannot have matrix in SDE format')
            raw_str = raw_tag[4:byte_count+4]
        else: # regular element
            byte_count = tag['byte_count'].item()
            # Deal with miMATRIX type (cannot pass byte string)
            if mdtype == miMATRIX:
                return self.current_getter(byte_count).get_array()
            # All other types can be read from string
            raw_str = self.mat_stream.read(byte_count)
            # Seek to next 64-bit boundary
            mod8 = byte_count % 8
            if mod8:
                self.mat_stream.seek(8 - mod8, 1)

        if mdtype in self.codecs: # encoded char data
            codec = self.codecs[mdtype]
            if not codec:
                raise TypeError, 'Do not support encoding %d' % mdtype
            el = raw_str.decode(codec)
        else: # numeric data
            dt = self.dtypes[mdtype]
            el_count = byte_count // dt.itemsize
            el = np.ndarray(shape=(el_count,),
                            dtype=dt,
                            buffer=raw_str)
            if copy:
                el = el.copy()

        return el

    def matrix_getter_factory(self):
        ''' Returns reader for next matrix at top level '''
        tag = self.read_dtype(self.dtypes['tag_full'])
        mdtype = tag['mdtype'].item()
        byte_count = tag['byte_count'].item()
        next_pos = self.mat_stream.tell() + byte_count
        if mdtype == miCOMPRESSED:
            getter = Mat5ZArrayReader(self, byte_count).matrix_getter_factory()
        elif not mdtype == miMATRIX:
            raise TypeError, \
                  'Expecting miMATRIX type here, got %d' %  mdtype
        else:
            getter = self.current_getter(byte_count)
        getter.next_position = next_pos
        return getter

    def current_getter(self, byte_count):
        ''' Return matrix getter for current stream position

        Returns matrix getters at top level and sub levels
        '''
        if not byte_count: # an empty miMATRIX can contain no bytes
            return Mat5EmptyMatrixGetter(self)
        af = self.read_dtype(self.dtypes['array_flags'])
        header = {}
        flags_class = af['flags_class']
        mc = flags_class & 0xFF
        header['mclass'] = mc
        header['is_logical'] = flags_class >> 9 & 1
        header['is_global'] = flags_class >> 10 & 1
        header['is_complex'] = flags_class >> 11 & 1
        header['nzmax'] = af['nzmax']
        header['dims'] = self.read_element()
        header['name'] = self.read_element().tostring()
        if mc in mx_numbers:
            return Mat5NumericMatrixGetter(self, header)
        if mc == mxSPARSE_CLASS:
            return Mat5SparseMatrixGetter(self, header)
        if mc == mxCHAR_CLASS:
            return Mat5CharMatrixGetter(self, header)
        if mc == mxCELL_CLASS:
            return Mat5CellMatrixGetter(self, header)
        if mc == mxSTRUCT_CLASS:
            return Mat5StructMatrixGetter(self, header)
        if mc == mxOBJECT_CLASS:
            return Mat5ObjectMatrixGetter(self, header)
        if mc == mxFUNCTION_CLASS:
            return Mat5FunctionMatrixGetter(self, header)
        raise TypeError, 'No reader for class code %s' % mc


class Mat5ZArrayReader(Mat5ArrayReader):
    ''' Getter for compressed arrays

    Reads and uncompresses gzipped stream on init, providing wrapper
    for this new sub-stream.
    '''
    def __init__(self, array_reader, byte_count):
        '''Reads and uncompresses gzipped stream'''
        data = array_reader.mat_stream.read(byte_count)
        super(Mat5ZArrayReader, self).__init__(
            StringIO(zlib.decompress(data)),
            array_reader.dtypes,
            array_reader.processor_func,
            array_reader.codecs,
            array_reader.class_dtypes,
            array_reader.struct_as_record)


class Mat5MatrixGetter(MatMatrixGetter):
    ''' Base class for getting Mat5 matrices

    Gets current read information from passed array_reader
    '''

    def __init__(self, array_reader, header):
        super(Mat5MatrixGetter, self).__init__(array_reader, header)
        self.class_dtypes = array_reader.class_dtypes
        self.codecs = array_reader.codecs
        self.is_global = header['is_global']
        self.mat_dtype = None

    def read_element(self, *args, **kwargs):
        return self.array_reader.read_element(*args, **kwargs)


class Mat5EmptyMatrixGetter(Mat5MatrixGetter):
    ''' Dummy class to return empty array for empty matrix
    '''
    def __init__(self, array_reader):
        self.array_reader = array_reader
        self.mat_stream = array_reader.mat_stream
        self.data_position = self.mat_stream.tell()
        self.header = {}
        self.name = ''
        self.is_global = False
        self.mat_dtype = 'f8'

    def get_raw_array(self):
        return np.array([[]])


class Mat5NumericMatrixGetter(Mat5MatrixGetter):

    def __init__(self, array_reader, header):
        super(Mat5NumericMatrixGetter, self).__init__(array_reader, header)
        if header['is_logical']:
            self.mat_dtype = np.dtype('bool')
        else:
            self.mat_dtype = self.class_dtypes[header['mclass']]

    def get_raw_array(self):
        if self.header['is_complex']:
            # avoid array copy to save memory
            res = self.read_element(copy=False)
            res_j = self.read_element(copy=False)
            res = res + (res_j * 1j)
        else:
            res = self.read_element()
        return np.ndarray(shape=self.header['dims'],
                          dtype=res.dtype,
                          buffer=res,
                          order='F')


class Mat5SparseMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        rowind = self.read_element()
        indptr = self.read_element()
        if self.header['is_complex']:
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
        may well have more entries than nnz, the actual number
        of non-zero entries, but rowind[nnz:] can be discarded
        and should be 0. indptr has length (number of columns + 1),
        and is such that, if D = diff(colind), D[j] gives the number
        of non-zero entries in column j. Because rowind values are
        stored in column order, this gives the column corresponding to
        each rowind
        '''
        M,N = self.header['dims']
        indptr = indptr[:N+1]
        nnz = indptr[-1]
        rowind = rowind[:nnz]
        data   = data[:nnz]
        if spsparse:
            return spsparse.csc_matrix((data,rowind,indptr), shape=(M,N))
        else:
            return ((M,N), data, rowind, indptr)


class Mat5CharMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        res = self.read_element()
        # Convert non-string types to unicode
        if isinstance(res, np.ndarray):
            if res.dtype.type == np.uint16:
                codec = miUINT16_codec
                if self.codecs['uint16_len'] == 1:
                    res = res.astype(np.uint8)
            elif res.dtype.type in (np.uint8, np.int8):
                codec = 'ascii'
            else:
                raise TypeError, 'Did not expect type %s' % res.dtype
            res = res.tostring().decode(codec)
        return np.ndarray(shape=self.header['dims'],
                          dtype=np.dtype('U1'),
                          buffer=np.array(res),
                          order='F').copy()


class Mat5CellMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        # Account for fortran indexing of cells
        tupdims = tuple(self.header['dims'][::-1])
        length = np.product(tupdims)
        result = np.empty(length, dtype=object)
        for i in range(length):
            result[i] = self.get_item()
        return result.reshape(tupdims).T

    def get_item(self):
        return self.read_element()

class mat_struct(object):
    ''' Placeholder for holding read data from structs '''
    pass

class Mat5StructMatrixGetter(Mat5MatrixGetter):
    def __init__(self, array_reader, header):
        super(Mat5StructMatrixGetter, self).__init__(array_reader, header)
        self.struct_as_record = array_reader.struct_as_record

    def get_raw_array(self):
        namelength = self.read_element()[0]
        names = self.read_element()
        field_names = [names[i:i+namelength].tostring().strip('\x00')
                       for i in xrange(0,len(names),namelength)]
        tupdims = tuple(self.header['dims'][::-1])
        length = np.product(tupdims)
        if self.struct_as_record:
            result = np.empty(length, dtype=[(field_name, object)
                                             for field_name in field_names])
            for i in range(length):
                for field_name in field_names:
                    result[i][field_name] = self.read_element()
        else: # Backward compatibility with previous format
            self.obj_template = mat_struct()
            self.obj_template._fieldnames = field_names
            result = np.empty(length, dtype=object)
            for i in range(length):
                item = pycopy(self.obj_template)
                for name in field_names:
                    item.__dict__[name] = self.read_element()
                result[i] = item

        return result.reshape(tupdims).T

class MatlabObject(object):
    ''' Class to contain read data from matlab objects '''
    def __init__(self, classname, field_names):
        self.__dict__['classname'] = classname
        self.__dict__['mobj_recarray'] = np.empty((1,1), dtype=[(field_name, object)
                                            for field_name in field_names])

    def __getattr__(self, name):
        mobj_recarray = self.__dict__['mobj_recarray']
        if name in mobj_recarray.dtype.fields:
            return mobj_recarray[0,0][name]
        else:
            raise AttributeError, "no field named %s in MatlabObject"%(name)

    def __setattr__(self, name, value):
        if name in self.__dict__['mobj_recarray'].dtype.fields:
            self.__dict__['mobj_recarray'][0,0][name] = value
        else:
            self.__dict__[name] = value


class Mat5ObjectMatrixGetter(Mat5MatrixGetter):
    def get_array(self):
        '''Matlab ojects are essentially structs, with an extra field, the classname.'''
        classname = self.read_element().tostring()
        namelength = self.read_element()[0]
        names = self.read_element()
        field_names = [names[i:i+namelength].tostring().strip('\x00')
                       for i in xrange(0,len(names),namelength)]
        result = MatlabObject(classname, field_names)

        for field_name in field_names:
            result.__setattr__(field_name, self.read_element())

        return result


class MatlabFunctionMatrix:
    ''' Opaque object representing an array of function handles. '''
    def __init__(self, arr):
        self.arr = arr

class Mat5FunctionMatrixGetter(Mat5CellMatrixGetter):
    def get_array(self):
        return MatlabFunctionMatrix(self.get_raw_array())


class MatFile5Reader(MatFileReader):
    ''' Reader for Mat 5 mat files
    Adds the following attribute to base class

    uint16_codec       - char codec to use for uint16 char arrays
                          (defaults to system default codec)
   '''

    def __init__(self,
                 mat_stream,
                 byte_order=None,
                 mat_dtype=False,
                 squeeze_me=False,
                 chars_as_strings=True,
                 matlab_compatible=False,
                 struct_as_record=False,
                 uint16_codec=None
                 ):
        '''
        mat_stream : file-like
                     object with file API, open for reading
        byte_order : {None, string}
                      specification of byte order, one of:
                      ('native', '=', 'little', '<', 'BIG', '>')
        mat_dtype : {True, False} boolean
                     If True, return arrays in same dtype as loaded into matlab
                     otherwise return with dtype with which they were saved
        squeeze_me : {False, True} boolean
                     If True, squeezes dimensions of size 1 from arrays
        chars_as_strings : {True, False} boolean
                     If True, convert char arrays to string arrays
        matlab_compatible : {False, True} boolean
                     If True, returns matrices as would be loaded by matlab
                     (implies squeeze_me=False, chars_as_strings=False
                     mat_dtype=True, struct_as_record=True)
        struct_as_record : {False, True} boolean
                     If True, return strutures as numpy records,
                     otherwise, return as custom object (for
                     compatibility with scipy 0.6)
        uint16_codec : {None, string}
                     Set codec to use for uint16 char arrays
                     (e.g. 'utf-8').  Use system default codec if None
        '''
        self.codecs = {}
        # Missing inputs to array reader set later (processor func
        # below, dtypes, codecs via our own set_dtype function, called
        # from parent __init__)
        self._array_reader = Mat5ArrayReader(
            mat_stream,
            None,
            None,
            None,
            None,
            struct_as_record
            )
        super(MatFile5Reader, self).__init__(
            mat_stream,
            byte_order,
            mat_dtype,
            squeeze_me,
            chars_as_strings,
            matlab_compatible,
            )
        self._array_reader.processor_func = self.processor_func
        self.uint16_codec = uint16_codec

    def get_uint16_codec(self):
        return self._uint16_codec
    def set_uint16_codec(self, uint16_codec):
        if not uint16_codec:
            uint16_codec = sys.getdefaultencoding()
        # Set length of miUINT16 char encoding
        self.codecs['uint16_len'] = len("  ".encode(uint16_codec)) \
                               - len(" ".encode(uint16_codec))
        self.codecs['uint16_codec'] = uint16_codec
        self._array_reader.codecs = self.codecs
        self._uint16_codec = uint16_codec
    uint16_codec = property(get_uint16_codec,
                            set_uint16_codec,
                            None,
                            'get/set uint16_codec')

    def set_dtypes(self):
        ''' Set dtypes and codecs '''
        self.dtypes = self.convert_dtypes(mdtypes_template)
        self.class_dtypes = self.convert_dtypes(mclass_dtypes_template)
        codecs = {}
        postfix = self.order_code == '<' and '_le' or '_be'
        for k, v in codecs_template.items():
            codec = v['codec']
            try:
                " ".encode(codec)
            except LookupError:
                codecs[k] = None
                continue
            if v['width'] > 1:
                codec += postfix
            codecs[k] = codec
        self.codecs.update(codecs)
        self.update_array_reader()

    def update_array_reader(self):
        self._array_reader.codecs = self.codecs
        self._array_reader.dtypes = self.dtypes
        self._array_reader.class_dtypes = self.class_dtypes

    def matrix_getter_factory(self):
        return self._array_reader.matrix_getter_factory()

    def guess_byte_order(self):
        ''' Guess byte order.
        Sets stream pointer to 0 '''
        self.mat_stream.seek(126)
        mi = self.mat_stream.read(2)
        self.mat_stream.seek(0)
        return mi == 'IM' and '<' or '>'

    def file_header(self):
        ''' Read in mat 5 file header '''
        hdict = {}
        hdr = self.read_dtype(self.dtypes['file_header'])
        hdict['__header__'] = hdr['description'].item().strip(' \t\n\000')
        v_major = hdr['version'] >> 8
        v_minor = hdr['version'] & 0xFF
        hdict['__version__'] = '%d.%d' % (v_major, v_minor)
        return hdict


class Mat5MatrixWriter(MatStreamWriter):

    mat_tag = np.zeros((), mdtypes_template['tag_full'])
    mat_tag['mdtype'] = miMATRIX

    def __init__(self, file_stream, arr, name, is_global=False):
        super(Mat5MatrixWriter, self).__init__(file_stream, arr, name)
        self.is_global = is_global

    def write_dtype(self, arr):
        self.file_stream.write(arr.tostring())

    def write_element(self, arr, mdtype=None):
        # write tag, data
        if mdtype is None:
            mdtype = np_to_mtypes[arr.dtype.str[1:]]
        byte_count = arr.size*arr.itemsize
        if byte_count <= 4:
            self.write_smalldata_element(arr, mdtype, byte_count)
        else:
            self.write_regular_element(arr, mdtype, byte_count)

    def write_smalldata_element(self, arr, mdtype, byte_count):
        # write tag with embedded data
        tag = np.zeros((), mdtypes_template['tag_smalldata'])
        tag['byte_count_mdtype'] = (byte_count << 16) + mdtype
        # if arr.tostring is < 4, the element will be zero-padded as needed.
        tag['data'] = arr.tostring(order='F')
        self.write_dtype(tag)

    def write_regular_element(self, arr, mdtype, byte_count):
        # write tag, data
        tag = np.zeros((), mdtypes_template['tag_full'])
        tag['mdtype'] = mdtype
        tag['byte_count'] = byte_count
        padding = (8 - tag['byte_count']) % 8
        self.write_dtype(tag)
        self.write_bytes(arr)
        # pad to next 64-bit boundary
        self.write_bytes(np.zeros((padding,),'u1'))

    def write_header(self, mclass,
                     is_global=False,
                     is_complex=False,
                     is_logical=False,
                     nzmax=0):
        ''' Write header for given data options
        mclass      - mat5 matrix class
        is_global   - True if matrix is global
        is_complex  - True if matrix is complex
        is_logical  - True if matrix is logical
        nzmax        - max non zero elements for sparse arrays
        '''
        self._mat_tag_pos = self.file_stream.tell()
        self.write_dtype(self.mat_tag)
        # write array flags (complex, global, logical, class, nzmax)
        af = np.zeros((), mdtypes_template['array_flags'])
        af['data_type'] = miUINT32
        af['byte_count'] = 8
        flags = is_complex << 3 | is_global << 2 | is_logical << 1
        af['flags_class'] = mclass | flags << 8
        af['nzmax'] = nzmax
        self.write_dtype(af)
        # write array shape
        if self.arr.ndim < 2:
            new_arr = np.atleast_2d(self.arr)
            if type(new_arr) != type(self.arr):
                raise ValueError("Array should be 2-dimensional.")
            self.arr = new_arr
        self.write_element(np.array(self.arr.shape, dtype='i4'))
        # write name
        self.write_element(np.array([ord(c) for c in self.name], 'i1'))

    def update_matrix_tag(self):
        curr_pos = self.file_stream.tell()
        self.file_stream.seek(self._mat_tag_pos)
        self.mat_tag['byte_count'] = curr_pos - self._mat_tag_pos - 8
        self.write_dtype(self.mat_tag)
        self.file_stream.seek(curr_pos)

    def write(self):
        assert False, 'Not implemented'


class Mat5NumericWriter(Mat5MatrixWriter):
    def write(self):
        imagf = self.arr.dtype.kind == 'c'
        try:
            mclass = np_to_mxtypes[self.arr.dtype.str[1:]]
        except KeyError:
            if imagf:
                self.arr = self.arr.astype('c128')
            else:
                self.arr = self.arr.astype('f8')
            mclass = mxDOUBLE_CLASS
        self.write_header(mclass=mclass,is_complex=imagf)
        if imagf:
            self.write_element(self.arr.real)
            self.write_element(self.arr.imag)
        else:
            self.write_element(self.arr)
        self.update_matrix_tag()

class Mat5CharWriter(Mat5MatrixWriter):
    codec='ascii'
    def write(self):
        self.arr_to_chars()
        self.write_header(mclass=mxCHAR_CLASS)
        if self.arr.dtype.kind == 'U':
            # Recode unicode using self.codec
            n_chars = np.product(self.arr.shape)
            st_arr = np.ndarray(shape=(),
                                dtype=self.arr_dtype_number(n_chars),
                                buffer=self.arr)
            st = st_arr.item().encode(self.codec)
            self.arr = np.ndarray(shape=(len(st)), dtype='u1', buffer=st)
        self.write_element(self.arr,mdtype=miUTF8)
        self.update_matrix_tag()

class Mat5UniCharWriter(Mat5CharWriter):
    codec='UTF8'


class Mat5SparseWriter(Mat5MatrixWriter):
    def write(self):
        ''' Sparse matrices are 2D
        '''
        A = self.arr.tocsc() # convert to sparse CSC format
        A.sort_indices()     # MATLAB expects sorted row indices
        is_complex = (A.dtype.kind == 'c')
        nz = A.nnz
        self.write_header(mclass=mxSPARSE_CLASS,
                          is_complex=is_complex,
                          nzmax=nz)
        self.write_element(A.indices.astype('i4'))
        self.write_element(A.indptr.astype('i4'))
        self.write_element(A.data.real)
        if is_complex:
            self.write_element(A.data.imag)
        self.update_matrix_tag()


class Mat5CompositeWriter(Mat5MatrixWriter):
    def __init__(self, file_stream, arr, name, is_global=False, unicode_strings=False):
        super(Mat5CompositeWriter, self).__init__(file_stream, arr, name, is_global)
        self.unicode_strings = unicode_strings


class Mat5CellWriter(Mat5CompositeWriter):
    def write(self):
        self.write_header(mclass=mxCELL_CLASS)
        # loop over data, column major
        A = np.atleast_2d(self.arr).flatten('F')
        MWG = Mat5WriterGetter(self.file_stream, self.unicode_strings)
        for el in A:
            MW = MWG.matrix_writer_factory(el, '')
            MW.write()
        self.update_matrix_tag()

class Mat5FunctionWriter(Mat5CompositeWriter):
    def __init__(self, file_stream, arr, name, is_global=False, unicode_strings=False):
        super(Mat5FunctionWriter, self).__init__(file_stream, arr.arr, name, is_global)

    def write(self):
        self.write_header(mclass=mxFUNCTION_CLASS)
        # loop over data, column major
        A = np.atleast_2d(self.arr).flatten('F')
        MWG = Mat5WriterGetter(self.file_stream, self.unicode_strings)
        for el in A:
            MW = MWG.matrix_writer_factory(el, '')
            MW.write()
        self.update_matrix_tag()


class Mat5StructWriter(Mat5CompositeWriter):
    def write(self):
        self.write_header(mclass=mxSTRUCT_CLASS)

        # write fieldnames
        fieldnames = [f[0] for f in self.arr.dtype.descr]
        self.write_element(np.array([32], dtype='i4'))
        self.write_element(np.array(fieldnames, dtype='S32'), mdtype=miINT8)

        A = np.atleast_2d(self.arr).flatten('F')
        MWG = Mat5WriterGetter(self.file_stream, self.unicode_strings)
        for el in A:
            for f in fieldnames:
                MW = MWG.matrix_writer_factory(el[f], '')
                MW.write()
        self.update_matrix_tag()

class Mat5ObjectWriter(Mat5CompositeWriter):
    def __init__(self, file_stream, arr, name, is_global=False, unicode_strings=False):
        super(Mat5ObjectWriter, self).__init__(file_stream, arr.__dict__['mobj_recarray'], name, is_global)
        self.classname = arr.classname

    def write(self):
        self.write_header(mclass=mxOBJECT_CLASS)

        # write classnames
        self.write_element(np.array(self.classname, dtype='S'), mdtype=miINT8)

        # write fieldnames
        fieldnames = [f[0] for f in self.arr.dtype.descr]
        self.write_element(np.array([32], dtype='i4'))
        self.write_element(np.array(fieldnames, dtype='S32'), mdtype=miINT8)

        A = np.atleast_2d(self.arr).flatten('F')
        MWG = Mat5WriterGetter(self.file_stream, self.unicode_strings)
        for el in A:
            for f in fieldnames:
                MW = MWG.matrix_writer_factory(el[f], '')
                MW.write()
        self.update_matrix_tag()


class Mat5WriterGetter(object):
    ''' Wraps stream and options, provides methods for getting Writer objects '''
    def __init__(self, stream, unicode_strings):
        self.stream = stream
        self.unicode_strings = unicode_strings

    def rewind(self):
        self.stream.seek(0)

    def matrix_writer_factory(self, arr, name, is_global=False):
        ''' Factory function to return matrix writer given variable to write
        stream      - file or file-like stream to write to
        arr         - array to write
        name        - name in matlab (TM) workspace
        '''
        if spsparse:
            if spsparse.issparse(arr):
                return Mat5SparseWriter(self.stream, arr, name, is_global)

        if isinstance(arr, MatlabFunctionMatrix):
            return Mat5FunctionWriter(self.stream, arr, name, is_global, self.unicode_strings)
        if isinstance(arr, MatlabObject):
            return Mat5ObjectWriter(self.stream, arr, name, is_global, self.unicode_strings)

        arr = np.array(arr)
        if arr.dtype.hasobject:
            if arr.dtype.fields == None:
                return Mat5CellWriter(self.stream, arr, name, is_global, self.unicode_strings)
            else:
                return Mat5StructWriter(self.stream, arr, name, is_global, self.unicode_strings)
        if arr.dtype.kind in ('U', 'S'):
            if self.unicode_strings:
                return Mat5UniCharWriter(self.stream, arr, name, is_global)
            else:
                return Mat5CharWriter(self.stream, arr, name, is_global)
        else:
            return Mat5NumericWriter(self.stream, arr, name, is_global)

class MatFile5Writer(MatFileWriter):
    ''' Class for writing mat5 files '''
    def __init__(self, file_stream,
                 do_compression=False,
                 unicode_strings=False,
                 global_vars=None):
        super(MatFile5Writer, self).__init__(file_stream)
        self.do_compression = do_compression
        if global_vars:
            self.global_vars = global_vars
        else:
            self.global_vars = []
        self.writer_getter = Mat5WriterGetter(
            StringIO(),
            unicode_strings)
        # write header
        import os, time
        hdr =  np.zeros((), mdtypes_template['file_header'])
        hdr['description']='MATLAB 5.0 MAT-file Platform: %s, Created on: %s' % (
                            os.name,time.asctime())
        hdr['version']= 0x0100
        hdr['endian_test']=np.ndarray(shape=(),dtype='S2',buffer=np.uint16(0x4d49))
        file_stream.write(hdr.tostring())

    def get_unicode_strings(self):
        return self.write_getter.unicode_strings
    def set_unicode_strings(self, unicode_strings):
        self.writer_getter.unicode_strings = unicode_strings
    unicode_strings = property(get_unicode_strings,
                               set_unicode_strings,
                               None,
                               'get/set unicode strings property')

    def put_variables(self, mdict):
        for name, var in mdict.items():
            if name[0] == '_':
                continue
            is_global = name in self.global_vars
            self.writer_getter.rewind()
            self.writer_getter.matrix_writer_factory(
                var,
                name,
                is_global,
                ).write()
            stream = self.writer_getter.stream
            if self.do_compression:
                str = zlib.compress(stream.getvalue(stream.tell()))
                tag = np.empty((), mdtypes_template['tag_full'])
                tag['mdtype'] = miCOMPRESSED
                tag['byte_count'] = len(str)
                self.file_stream.write(tag.tostring() + str)
            else:
                self.file_stream.write(stream.getvalue(stream.tell()))
