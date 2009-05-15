''' Classes for read / write of matlab (TM) 5 files

The matfile specification last found here:

http://www.mathworks.com/access/helpdesk/help/pdf_doc/matlab/matfile_format.pdf

(as of December 5 2008)
'''

# Small fragments of current code adapted from matfile.py by Heiko
# Henkelmann

import os
import time
import sys
import zlib
from StringIO import StringIO
from cStringIO import StringIO as cStringIO
from copy import copy as pycopy
import warnings

import numpy as np

import scipy.sparse

import byteordercodes
from miobase import MatFileReader, MatArrayReader, MatMatrixGetter, \
     MatFileWriter, MatStreamWriter, docfiller, matdims, \
     MatReadError

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
mxOPAQUE_CLASS = 17 # This appears to be a function workspace
# https://www-old.cae.wisc.edu/pipermail/octave-maintainers/2007-May/002824.html
mxOBJECT_CLASS_FROM_MATRIX_H = 18

mxmap = { # Sometimes good for debug prints
    mxCELL_CLASS: 'mxCELL_CLASS',
    mxSTRUCT_CLASS: 'mxSTRUCT_CLASS',
    mxOBJECT_CLASS: 'mxOBJECT_CLASS',
    mxCHAR_CLASS: 'mxCHAR_CLASS',
    mxSPARSE_CLASS: 'mxSPARSE_CLASS',
    mxDOUBLE_CLASS: 'mxDOUBLE_CLASS',
    mxSINGLE_CLASS: 'mxSINGLE_CLASS',
    mxINT8_CLASS: 'mxINT8_CLASS',
    mxUINT8_CLASS: 'mxUINT8_CLASS',
    mxINT16_CLASS: 'mxINT16_CLASS',
    mxUINT16_CLASS: 'mxUINT16_CLASS',
    mxINT32_CLASS: 'mxINT32_CLASS',
    mxUINT32_CLASS: 'mxUINT32_CLASS',
    mxINT64_CLASS: 'mxINT64_CLASS',
    mxUINT64_CLASS: 'mxUINT64_CLASS',
    mxFUNCTION_CLASS: 'mxFUNCTION_CLASS',
    mxOPAQUE_CLASS: 'mxOPAQUE_CLASS',
    mxOBJECT_CLASS_FROM_MATRIX_H: 'mxOBJECT_CLASS_FROM_MATRIX_H',
}

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


class mat_struct(object):
    ''' Placeholder for holding read data from structs

    We will deprecate this method of holding struct information in a
    future version of scipy, in favor of the recarray method (see
    loadmat docstring)
    '''
    pass


class MatlabObject(np.ndarray):
    ''' ndarray Subclass to contain matlab object '''
    def __new__(cls, input_array, classname=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.classname = classname
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self,obj):
        # reset the attribute from passed original object
        self.classname = getattr(obj, 'classname', None)
        # We do not need to return anything


class MatlabFunction(np.ndarray):
    ''' Subclass to signal this is a matlab function '''
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)


class MatlabBinaryBlock(object):
    ''' Class to contain matlab unreadable blocks '''
    def __init__(self, binaryblock, endian):
        self.binaryblock = binaryblock
        self.endian = endian


class Mat5ArrayReader(MatArrayReader):
    ''' Class to get Mat5 arrays

    Provides element reader functions, header reader, matrix reader
    factory function
    '''

    def __init__(self,
                 mat_stream,
                 dtypes,
                 processor_func,
                 codecs,
                 class_dtypes,
                 struct_as_record):
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
        ''' Here I am playing with a binary block read of
        untranslatable data. I am not using this at the moment because
        reading it has the side effect of making opposite ending mat
        files unwritable on the round trip.
        
        if mc == mxFUNCTION_CLASS:
            # we can't read these, and want to keep track of the byte
            # count - so we need to avoid the following unpredictable
            # length element reads
            return Mat5BinaryBlockGetter(self,
                                         header,
                                         af,
                                         byte_count)
        '''
        header['dims'] = self.read_element()
        header['name'] = self.read_element().tostring()
        # maybe a dictionary mapping here as a dispatch table
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
            return Mat5FunctionGetter(self, header)
        raise TypeError, 'No reader for class code %s' % mc


class Mat5ZArrayReader(Mat5ArrayReader):
    ''' Getter for compressed arrays

    Sets up reader for gzipped stream on init, providing wrapper
    for this new sub-stream.

    '''
    def __init__(self, array_reader, byte_count):
        super(Mat5ZArrayReader, self).__init__(
            cStringIO(zlib.decompress(
                        array_reader.mat_stream.read(byte_count))),
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
        return scipy.sparse.csc_matrix(
            (data,rowind,indptr),
            shape=(M,N))


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
            if not len(field_names):
                # If there are no field names, there is no dtype
                # representation we can use, falling back to empty
                # object
                return np.empty(tupdims, dtype=object).T
            dtype = [(field_name, object) for field_name in field_names]
            result = np.empty(length, dtype=dtype)
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


class Mat5ObjectMatrixGetter(Mat5StructMatrixGetter):
    def get_raw_array(self):
        '''Matlab objects are like structs, with an extra classname field'''
        classname = self.read_element().tostring()
        result = super(Mat5ObjectMatrixGetter, self).get_raw_array()
        return MatlabObject(result, classname)


class Mat5FunctionGetter(Mat5ObjectMatrixGetter):
    ''' Class to provide warning and message string for unreadable
    matlab function data
    '''
    def get_raw_array(self):
        raise MatReadError('Cannot read matlab functions')


class Mat5BinaryBlockGetter(object):
    ''' Class to read in unreadable binary blocks

    This class could be used to read in matlab functions
    '''

    def __init__(self,
                 array_reader,
                 header,
                 array_flags,
                 byte_count):
        self.array_reader = array_reader
        self.header = header
        self.array_flags = array_flags
        arr_str = array_flags.tostring()
        self.binaryblock = array_reader.mat_stream.read(
            byte_count-len(array_flags.tostring()))
        stream = StringIO(self.binaryblock)
        reader = Mat5ArrayReader(
            stream,
            array_reader.dtypes,
            lambda x : None,
            array_reader.codecs,
            array_reader.class_dtypes,
            False)
        self.header['dims'] = reader.read_element()
        self.header['name'] = reader.read_element().tostring()
        self.name = self.header['name']
        self.is_global = header['is_global']

    def get_array(self):
        dt = self.array_reader.dtypes[miINT32]
        endian = byteordercodes.to_numpy_code(dt.byteorder)
        data = self.array_flags.tostring() + self.binaryblock
        return MatlabBinaryBlock(data, endian)

               
class MatFile5Reader(MatFileReader):
    ''' Reader for Mat 5 mat files
    Adds the following attribute to base class

    uint16_codec       - char codec to use for uint16 char arrays
                          (defaults to system default codec)
   '''
    @docfiller
    def __init__(self,
                 mat_stream,
                 byte_order=None,
                 mat_dtype=False,
                 squeeze_me=False,
                 chars_as_strings=True,
                 matlab_compatible=False,
                 struct_as_record=None, # default False, for now
                 uint16_codec=None
                 ):
        '''Initializer for matlab 5 file format reader

    %(matstream_arg)s
    %(load_args)s
    %(struct_arg)s
    uint16_codec : {None, string}
        Set codec to use for uint16 char arrays (e.g. 'utf-8').
        Use system default codec if None
        '''
        # Deal with deprecations
        if struct_as_record is None:
            warnings.warn("Using struct_as_record default value (False)" +
                          " This will change to True in future versions",
                          FutureWarning, stacklevel=2)
            struct_as_record = False
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
    ''' Generic matlab matrix writing class '''
    mat_tag = np.zeros((), mdtypes_template['tag_full'])
    mat_tag['mdtype'] = miMATRIX
    default_mclass = None # default class for header writing
    def __init__(self,
                 file_stream,
                 arr,
                 name,
                 is_global=False,
                 unicode_strings=False,
                 long_field_names=False,
                 oned_as='column'):
        super(Mat5MatrixWriter, self).__init__(file_stream, 
                                               arr, 
                                               name,
                                               oned_as)
        self.is_global = is_global
        self.unicode_strings = unicode_strings
        self.long_field_names = long_field_names
        self.oned_as = oned_as

    def write_dtype(self, arr):
        self.file_stream.write(arr.tostring())

    def write_element(self, arr, mdtype=None):
        ''' write tag and data '''
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

    def write_header(self, mclass=None,
                     is_global=False,
                     is_complex=False,
                     is_logical=False,
                     nzmax=0,
                     shape=None):
        ''' Write header for given data options
        mclass      - mat5 matrix class
        is_global   - True if matrix is global
        is_complex  - True if matrix is complex
        is_logical  - True if matrix is logical
        nzmax        - max non zero elements for sparse arrays
        shape : {None, tuple} optional
            directly specify shape if this is not the same as for
            self.arr
        '''
        if mclass is None:
            mclass = self.default_mclass
        if shape is None:
            shape = matdims(self.arr, self.oned_as)
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
        self.write_element(np.array(shape, dtype='i4'))
        # write name
        self.write_element(np.array([ord(c) for c in self.name], 'i1'))

    def update_matrix_tag(self):
        curr_pos = self.file_stream.tell()
        self.file_stream.seek(self._mat_tag_pos)
        self.mat_tag['byte_count'] = curr_pos - self._mat_tag_pos - 8
        self.write_dtype(self.mat_tag)
        self.file_stream.seek(curr_pos)

    def write(self):
        raise NotImplementedError

    def make_writer_getter(self):
        ''' Make writer getter for this stream '''
        return Mat5WriterGetter(self.unicode_strings,
                                self.long_field_names,
                                self.oned_as)


class Mat5NumericWriter(Mat5MatrixWriter):
    default_mclass = None # can be any numeric type
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
    default_mclass = mxCHAR_CLASS
    def write(self):
        self.arr_to_chars()
        # We have to write the shape directly, because we are going
        # recode the characters, and the resulting stream of chars
        # may have a different length
        shape = self.arr.shape
        self.write_header(shape=shape)
        # We need to do our own transpose (not using the normal
        # write routines that do this for us)
        arr = self.arr.T.copy()
        if self.arr.dtype.kind == 'U' and arr.size:
            # Recode unicode using self.codec
            n_chars = np.product(shape)
            st_arr = np.ndarray(shape=(),
                                dtype=self.arr_dtype_number(n_chars),
                                buffer=arr)
            st = st_arr.item().encode(self.codec)
            arr = np.ndarray(shape=(len(st),),
                             dtype='u1',
                             buffer=st)
        self.write_element(arr, mdtype=miUTF8)
        self.update_matrix_tag()


class Mat5UniCharWriter(Mat5CharWriter):
    codec='UTF8'


class Mat5SparseWriter(Mat5MatrixWriter):
    default_mclass = mxSPARSE_CLASS
    def write(self):
        ''' Sparse matrices are 2D
        '''
        A = self.arr.tocsc() # convert to sparse CSC format
        A.sort_indices()     # MATLAB expects sorted row indices
        is_complex = (A.dtype.kind == 'c')
        nz = A.nnz
        self.write_header(is_complex=is_complex,
                          nzmax=nz)
        self.write_element(A.indices.astype('i4'))
        self.write_element(A.indptr.astype('i4'))
        self.write_element(A.data.real)
        if is_complex:
            self.write_element(A.data.imag)
        self.update_matrix_tag()


class Mat5CellWriter(Mat5MatrixWriter):
    default_mclass = mxCELL_CLASS
    def write(self):
        self.write_header()
        self._write_items()

    def _write_items(self):
        # loop over data, column major
        A = np.atleast_2d(self.arr).flatten('F')
        MWG = self.make_writer_getter()
        for el in A:
            MW = MWG.matrix_writer_factory(self.file_stream, el)
            MW.write()
        self.update_matrix_tag()


class Mat5BinaryBlockWriter(Mat5MatrixWriter):
    ''' class to write untranslatable binary blocks '''
    def write(self):
        # check endian
        # write binary block as is
        pass

class Mat5StructWriter(Mat5CellWriter):
    ''' class to write matlab structs

    Differs from cell writing class in writing field names,
    and in mx class
    '''
    default_mclass = mxSTRUCT_CLASS

    def _write_items(self):
        # write fieldnames
        fieldnames = [f[0] for f in self.arr.dtype.descr]
        length = max([len(fieldname) for fieldname in fieldnames])+1
        max_length = (self.long_field_names and 64) or 32
        if length > max_length:
            raise ValueError(
                "Field names are restricted to %d characters"
                 % (max_length-1))
        self.write_element(np.array([length], dtype='i4'))
        self.write_element(
            np.array(fieldnames, dtype='S%d'%(length)),
            mdtype=miINT8)
        A = np.atleast_2d(self.arr).flatten('F')
        MWG = self.make_writer_getter()
        for el in A:
            for f in fieldnames:
                MW = MWG.matrix_writer_factory(self.file_stream, el[f])
                MW.write()
        self.update_matrix_tag()


class Mat5ObjectWriter(Mat5StructWriter):
    ''' class to write matlab objects

    Same as writing structs, except different mx class, and extra
    classname element after header
    '''
    default_mclass = mxOBJECT_CLASS
    def write(self):
        self.write_header()
        self.write_element(np.array(self.arr.classname, dtype='S'),
                           mdtype=miINT8)
        self._write_items()


class Mat5WriterGetter(object):
    ''' Wraps options, provides methods for getting Writer objects '''
    @docfiller
    def __init__(self, 
                 unicode_strings=True, 
                 long_field_names=False,
                 oned_as='column'):
        ''' Initialize writer getter

        Parameters
        ----------
        unicode_strings : bool
           If True, write unicode strings
        %(long_fields)s
        %(oned_as)s
        '''
        self.unicode_strings = unicode_strings
        self.long_field_names = long_field_names
        self.oned_as = oned_as

    def to_writeable(self, source):
        ''' Convert input object ``source`` to something we can write

        Parameters
        ----------
        source : object

        Returns
        -------
        arr : ndarray

        Examples
        --------
        >>> mwg = Mat5WriterGetter()
        >>> mwg.to_writeable(np.array([1])) # pass through ndarrays
        array([1])
        >>> expected = np.array([(1, 2)], dtype=[('a', '|O8'), ('b', '|O8')])
        >>> np.all(mwg.to_writeable({'a':1,'b':2}) == expected)
        True
        >>> np.all(mwg.to_writeable({'a':1,'b':2, '_c':3}) == expected)
        True
        >>> np.all(mwg.to_writeable({'a':1,'b':2, 100:3}) == expected)
        True
        >>> np.all(mwg.to_writeable({'a':1,'b':2, '99':3}) == expected)
        True
        >>> class klass(object): pass
        >>> c = klass
        >>> c.a = 1
        >>> c.b = 2
        >>> np.all(mwg.to_writeable({'a':1,'b':2}) == expected)
        True
        >>> mwg.to_writeable([])
        array([], dtype=float64)
        >>> mwg.to_writeable(())
        array([], dtype=float64)
        >>> mwg.to_writeable(None)

        >>> mwg.to_writeable('a string').dtype
        dtype('|S8')
        >>> mwg.to_writeable(1)
        array(1)
        >>> mwg.to_writeable([1])
        array([1])
        >>> mwg.to_writeable([1])
        array([1])
        >>> mwg.to_writeable(object()) # not convertable

        dict keys with legal characters are convertible

        >>> mwg.to_writeable({'a':1})['a']
        array([1], dtype=object)

        but not with illegal characters

        >>> mwg.to_writeable({'1':1}) is None
        True
        >>> mwg.to_writeable({'_a':1}) is None
        True
        '''
        if isinstance(source, np.ndarray):
            return source
        if source is None:
            return None
        # Objects that have dicts
        if hasattr(source, '__dict__'):
            source = dict((key, value) for key, value in source.__dict__.items()
                          if not key.startswith('_'))
        # Mappings or object dicts
        if hasattr(source, 'keys'):
            dtype = []
            values = []
            for field, value in source.items():
                if (isinstance(field, basestring) and 
                    not field[0] in '_0123456789'):
                    dtype.append((field,object))
                    values.append(value)
            if dtype:
                return np.array( [tuple(values)] ,dtype)
            else:
                return None
        # Next try and convert to an array
        narr = np.asanyarray(source)
        if narr.dtype.type in (np.object, np.object_) and \
           narr.shape == () and narr == source:
            # No interesting conversion possible
            return None
        return narr

    def matrix_writer_factory(self, stream, arr, name='', is_global=False):
        ''' Factory function to return matrix writer given variable to write

        Parameters
        ----------
        stream : fileobj
            stream to write to
        arr : array-like
            array-like object to create writer for
        name : string
            name as it will appear in matlab workspace
            default is empty string
        is_global : {False, True} optional
            whether variable will be global on load into matlab

        Returns
        -------
        writer : matrix writer object
        '''
        # First check if these are sparse
        if scipy.sparse.issparse(arr):
            return Mat5SparseWriter(stream, arr, name, is_global)
        # Try to convert things that aren't arrays
        narr = self.to_writeable(arr)
        if narr is None:
            raise TypeError('Could not convert %s (type %s) to array'
                            % (arr, type(arr)))
        args = (stream,
                narr,
                name,
                is_global,
                self.unicode_strings,
                self.long_field_names,
                self.oned_as)
        if isinstance(narr, MatlabBinaryBlock):
            return Mat5BinaryBlockWriter(*args)
        if isinstance(narr, MatlabObject):
            return Mat5ObjectWriter(*args)
        if narr.dtype.fields: # struct array
            return Mat5StructWriter(*args)
        if narr.dtype.hasobject: # cell array
            return Mat5CellWriter(*args)
        if narr.dtype.kind in ('U', 'S'):
            if self.unicode_strings:
                return Mat5UniCharWriter(*args)
            else:
                return Mat5CharWriter(*args)
        else:
            return Mat5NumericWriter(*args)


class MatFile5Writer(MatFileWriter):
    ''' Class for writing mat5 files '''
    @docfiller
    def __init__(self, file_stream,
                 do_compression=False,
                 unicode_strings=False,
                 global_vars=None,
                 long_field_names=False,
                 oned_as=None):
        ''' Initialize writer for matlab 5 format files 

        Parameters
        ----------
        %(do_compression)s
        %(unicode_strings)s
        global_vars : None or sequence of strings, optional
            Names of variables to be marked as global for matlab
        %(long_fields)s
        %(oned_as)s
        '''
        super(MatFile5Writer, self).__init__(file_stream)
        self.do_compression = do_compression
        if global_vars:
            self.global_vars = global_vars
        else:
            self.global_vars = []
        # deal with deprecations
        if oned_as is None:
            warnings.warn("Using oned_as default value ('column')" +
                          " This will change to 'row' in future versions",
                          FutureWarning, stacklevel=2)
            oned_as = 'column'
        self.writer_getter = Mat5WriterGetter(
            unicode_strings,
            long_field_names,
            oned_as)
        # write header
        hdr =  np.zeros((), mdtypes_template['file_header'])
        hdr['description']='MATLAB 5.0 MAT-file Platform: %s, Created on: %s' \
            % (os.name,time.asctime())
        hdr['version']= 0x0100
        hdr['endian_test']=np.ndarray(shape=(),
                                      dtype='S2',
                                      buffer=np.uint16(0x4d49))
        file_stream.write(hdr.tostring())

    def get_unicode_strings(self):
        return self.writer_getter.unicode_strings
    def set_unicode_strings(self, unicode_strings):
        self.writer_getter.unicode_strings = unicode_strings
    unicode_strings = property(get_unicode_strings,
                               set_unicode_strings,
                               None,
                               'get/set unicode strings property')

    def get_long_field_names(self):
        return self.writer_getter.long_field_names
    def set_long_field_names(self, long_field_names):
        self.writer_getter.long_field_names = long_field_names
    long_field_names = property(get_long_field_names,
                                set_long_field_names,
                                None,
                                'enable writing 32-63 character field '
                                'names for Matlab 7.6+')

    def get_oned_as(self):
        return self.writer_getter.oned_as
    def set_oned_as(self, oned_as):
        self.writer_getter.oned_as = oned_as
    oned_as = property(get_oned_as,
                       set_oned_as,
                       None,
                       'get/set oned_as property')

    def put_variables(self, mdict):
        for name, var in mdict.items():
            if name[0] == '_':
                continue
            is_global = name in self.global_vars
            if self.do_compression:
                stream = StringIO()
                mat_writer = self.writer_getter.matrix_writer_factory(
                    stream,
                    var,
                    name,
                    is_global)
                mat_writer.write()
                out_str = zlib.compress(stream.getvalue())
                tag = np.empty((), mdtypes_template['tag_full'])
                tag['mdtype'] = miCOMPRESSED
                tag['byte_count'] = len(out_str)
                self.file_stream.write(tag.tostring() + out_str)
            else: # not compressing
                mat_writer = self.writer_getter.matrix_writer_factory(
                    self.file_stream,
                    var,
                    name,
                    is_global)
                mat_writer.write()
