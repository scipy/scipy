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
from cStringIO import StringIO
import warnings

import numpy as np

import scipy.sparse

from miobase import MatFileReader, \
     MatFileWriter, MatStreamWriter, docfiller, matdims, \
     read_dtype, convert_dtypes

# Reader object for matlab 5 format variables
from mio5_utils import VarReader5

# Constants and helper objects
from mio5_params import MatlabObject, \
    miINT8, miUINT8, miINT16, miUINT16, miINT32, miUINT32, \
    miSINGLE, miDOUBLE, miINT64, miUINT64, miMATRIX, \
    miCOMPRESSED, miUTF8, miUTF16, miUTF32, \
    mxCELL_CLASS, mxSTRUCT_CLASS, mxOBJECT_CLASS, mxCHAR_CLASS, \
    mxSPARSE_CLASS, mxDOUBLE_CLASS, mxSINGLE_CLASS, mxINT8_CLASS, \
    mxUINT8_CLASS, mxINT16_CLASS, mxUINT16_CLASS, mxINT32_CLASS, \
    mxUINT32_CLASS, mxINT64_CLASS, mxUINT64_CLASS


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


def convert_codecs(template, byte_order):
    ''' Convert codec template mapping to byte order

    Set codecs not on this system to None

    Parameters
    ----------
    template : mapping
       key, value are respectively codec name, and root name for codec
       (without byte order suffix)
    byte_order : {'<', '>'}
       code for little or big endian

    Returns
    -------
    codecs : dict
       key, value are name, codec (as in .encode(codec))
    '''
    codecs = {}
    postfix = byte_order == '<' and '_le' or '_be'
    for k, v in template.items():
        codec = v['codec']
        try:
            " ".encode(codec)
        except LookupError:
            codecs[k] = None
            continue
        if v['width'] > 1:
            codec += postfix
        codecs[k] = codec
    return codecs.copy()


class MatFile5Reader(MatFileReader):
    ''' Reader for Mat 5 mat files
    Adds the following attribute to base class
    
    uint16_codec - char codec to use for uint16 char arrays
        (defaults to system default codec)

    Uses variable reader that has the following stardard interface (see
    abstract class in ``miobase``::
    
       __init__(self, file_reader)
       read_header(self)
       array_from_header(self)

    and added interface::

       set_stream(self, stream)
       read_full_tag(self)
       
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
        super(MatFile5Reader, self).__init__(
            mat_stream,
            byte_order,
            mat_dtype,
            squeeze_me,
            chars_as_strings,
            matlab_compatible,
            struct_as_record
            )
        # Set uint16 codec
        if not uint16_codec:
            uint16_codec = sys.getdefaultencoding()
        self.uint16_codec = uint16_codec
        # placeholders for dtypes, codecs - see initialize_read
        self.dtypes = None
        self.class_dtypes = None
        self.codecs = None
        # placeholders for readers - see initialize_read method
        self._file_reader = None
        self._matrix_reader = None

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
        hdr = read_dtype(self.mat_stream, self.dtypes['file_header'])
        hdict['__header__'] = hdr['description'].item().strip(' \t\n\000')
        v_major = hdr['version'] >> 8
        v_minor = hdr['version'] & 0xFF
        hdict['__version__'] = '%d.%d' % (v_major, v_minor)
        return hdict

    def initialize_read(self):
        ''' Run when beginning read of variables

        Sets up readers from parameters in `self`
        '''
        self.dtypes = convert_dtypes(mdtypes_template, self.byte_order)
        self.class_dtypes = convert_dtypes(mclass_dtypes_template,
                                           self.byte_order)
        self.codecs = convert_codecs(codecs_template, self.byte_order)
        uint16_codec = self.uint16_codec
        # Set length of miUINT16 char encoding
        self.codecs['uint16_len'] = len("  ".encode(uint16_codec)) \
                               - len(" ".encode(uint16_codec))
        self.codecs['uint16_codec'] = uint16_codec
        # reader for top level stream.  We need this extra top-level
        # reader because we use the matrix_reader object to contain
        # compressed matrices (so they have their own stream)
        self._file_reader = VarReader5(self)
        # reader for matrix streams 
        self._matrix_reader = VarReader5(self)

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
        mdtype, byte_count = self._file_reader.read_full_tag()
        assert byte_count > 0
        next_pos = self.mat_stream.tell() + byte_count
        if mdtype == miCOMPRESSED: # make new stream from compressed data
            stream = StringIO(zlib.decompress(self.mat_stream.read(byte_count)))
            self._matrix_reader.set_stream(stream)
            mdtype, byte_count = self._matrix_reader.read_full_tag()
        else:
            self._matrix_reader.set_stream(self.mat_stream)
        if not mdtype == miMATRIX:
            raise TypeError, \
                'Expecting miMATRIX type here, got %d' %  mdtype
        header = self._matrix_reader.read_header()
        return header, next_pos
            
    
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
