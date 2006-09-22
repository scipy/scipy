''' Classes for read / write of matlab (TM) 5 files
'''

# Small fragments of current code adapted from matfile.py by Heiko
# Henkelmann

## Notice in matfile.py file

# Copyright (c) 2003 Heiko Henkelmann

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import zlib
from copy import copy as pycopy
from cStringIO import StringIO
from numpy import *

from miobase import *

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
    mxSINGLE_CLASS: 'f4',
    mxDOUBLE_CLASS: 'f8',
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
    )

class mat_struct(object):
    ''' Placeholder for holding read data from structs '''
    pass

class mat_obj(object):
    ''' Placeholder for holding read data from objects '''
    pass

class Mat5ArrayReader(MatArrayReader):
    ''' Class to get Mat5 arrays

    Provides element reader functions, header reader, matrix reader
    factory function
    '''

    def __init__(self, mat_stream, dtypes, processor_func, codecs, class_dtypes):
        super(Mat5ArrayReader, self).__init__(mat_stream,
                                              dtypes,
                                              processor_func,
                                              )
        self.codecs = codecs
        self.class_dtypes = class_dtypes

    def read_element(self, copy=True):
        raw_tag = self.mat_stream.read(8)
        tag = ndarray(shape=(),
                      dtype=self.dtypes['tag_full'],
                      buffer = raw_tag)
        mdtype = tag['mdtype']
        byte_count = mdtype >> 16
        if byte_count: # small data element format
            if byte_count > 4:
                raise ValueError, 'Too many bytes for sde format'
            mdtype = mdtype & 0xFFFF
            dt = self.dtypes[mdtype]
            el_count = byte_count / dt.itemsize
            return ndarray(shape=(el_count,),
                           dtype=dt,
                           buffer=raw_tag[4:])
        byte_count = tag['byte_count']
        if mdtype == miMATRIX:
            return self.getter_from_bytes(byte_count).get_array()
        if mdtype in self.codecs: # encoded char data
           raw_str = self.mat_stream.read(byte_count)
           codec = self.codecs[mdtype]
           if not codec:
               raise TypeError, 'Do not support encoding %d' % mdtype
           el = raw_str.decode(codec)
        else: # numeric data
            dt = self.dtypes[mdtype]
            el_count = byte_count / dt.itemsize
            el = ndarray(shape=(el_count,),
                         dtype=dt,
                         buffer=self.mat_stream.read(byte_count))
            if copy:
                el = el.copy()
        mod8 = byte_count % 8
        if mod8:
            self.mat_stream.seek(8 - mod8, 1)
        return el

    def matrix_getter_factory(self):
        ''' Returns reader for next matrix '''
        tag = self.read_dtype(self.dtypes['tag_full'])
        mdtype = tag['mdtype']
        byte_count = tag['byte_count']
        if mdtype == miCOMPRESSED:
            return Mat5ZArrayReader(self, byte_count).matrix_getter_factory()
        if not mdtype == miMATRIX:
            raise TypeError, \
                  'Expecting miMATRIX type here, got %d' %  mdtype
        return self.getter_from_bytes(byte_count)

    def getter_from_bytes(self, byte_count):
        ''' Return matrix getter for current stream position '''
        # Apparently an empty miMATRIX can contain no bytes
        if not byte_count:
            return Mat5EmptyMatrixGetter(self)
        af = self.read_dtype(self.dtypes['array_flags'])
        header = {}
        flags_class = af['flags_class']
        header['next_position'] = self.mat_stream.tell() + byte_count
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
        raise TypeError, 'No reader for class code %s' % mc


class Mat5ZArrayReader(Mat5ArrayReader):
    ''' Getter for compressed arrays

    Reads and uncompresses gzipped stream on init, providing wrapper
    for this new sub-stream.  Sets next_position for main stream to
    allow skipping over this variable (although we have to read and
    uncompress the whole thing anyway to get the name)
    '''
    def __init__(self, array_reader, byte_count):
        '''Reads and uncompresses gzipped stream'''
        data = array_reader.mat_stream.read(byte_count)
        super(Mat5ZArrayReader, self).__init__(
            StringIO(zlib.decompress(data)),
            array_reader.dtypes,
            array_reader.processor_func,
            array_reader.codecs,
            array_reader.class_dtypes)
        self._next_position = array_reader.mat_stream.tell()
        
    def getter_from_bytes(self, byte_count):
        ''' Set next_position to current position in parent stream
        
        self.next_position is only used by the get_variables routine
        of the main file reading loop, so must refer to the position
        in the main stream, not the compressed stream.
        '''
        getter = super(Mat5ZArrayReader, self).getter_from_bytes(byte_count)
        getter.header['next_position'] = self._next_position
        return getter
    

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
        self.is_global = False
        self.mat_dtype = 'f8'
    
    def get_raw_array(self):
        return array([[]])


class Mat5NumericMatrixGetter(Mat5MatrixGetter):

    def __init__(self, array_reader, header):
        super(Mat5NumericMatrixGetter, self).__init__(array_reader, header)
        if header['is_logical']:
            self.mat_dtype = dtype('bool')
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
        return ndarray(shape=self.header['dims'],
                       dtype=res.dtype,
                       buffer=res,
                       order='F')
        

class Mat5SparseMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        rowind  = self.read_element()
        colind = self.read_element()
        if self.header['is_complex']:
            # avoid array copy to save memory
            res = self.read_element(copy=False)
            res_j = self.read_element(copy=False)
            res = res + (res_j * 1j)
        else:
            res = self.read_element()
        ''' From the matlab (TM) API documentation, last found here:
        http://www.mathworks.com/access/helpdesk/help/techdoc/matlab_external/
        @rowind are simply the row indices for all the (@res) non-zero
        entries in the sparse array.  @rowind has nzmax entries, so
        may well have more entries than len(@res), the actual number
        of non-zero entries, but @rowind[len(res):] can be discarded
        and should be 0. @colind has length (number of columns + 1),
        and is such that, if D = diff(@colind), D[j] gives the number
        of non-zero entries in column j. Because @rowind values are
        stored in column order, this gives the column corresponding to
        each @rowind
        '''
        cols = empty((len(res)), dtype=rowind.dtype)
        col_counts = diff(colind)
        start_row = 0
        for i in where(col_counts)[0]:
            end_row = start_row + col_counts[i]
            cols[start_row:end_row] = i
            start_row = end_row
        ij = vstack((rowind[:len(res)], cols))
        if have_sparse:
            result = scipy.sparse.csc_matrix((res,ij),
                                             self.header['dims'])
        else:
            result = (dims, ij, res)
        return result


class Mat5CharMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        res = self.read_element()
        # Convert non-string types to unicode
        if isinstance(res, ndarray):
            if res.dtype.type == uint16:
                codec = miUINT16_codec
                if self.codecs['uint16_len'] == 1:
                    res = res.astype(uint8)
            elif res.dtype.type in (uint8, int8):
                codec = 'ascii'
            else:
                raise TypeError, 'Did not expect type %s' % res.dtype
            res = res.tostring().decode(codec)
        return ndarray(shape=self.header['dims'],
                       dtype=dtype('U1'),
                       buffer=array(res),
                       order='F').copy()


class Mat5CellMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        # Account for fortran indexing of cells
        tupdims = tuple(self.header['dims'][::-1])
        length = product(tupdims)
        result = empty(length, dtype=object)
        for i in range(length):
            result[i] = self.get_item()
        result = transpose(reshape(result,tupdims))
        return result

    def get_item(self):
        return self.read_element()


class Mat5StructMatrixGetter(Mat5CellMatrixGetter):
    def __init__(self, *args, **kwargs):
        super(Mat5StructMatrixGetter, self).__init__(*args, **kwargs)
        self.obj_template = mat_struct()
        
    def get_raw_array(self):
        namelength = self.read_element()
        # get field names
        names = self.read_element()
        splitnames = [names[i:i+namelength] for i in \
                      xrange(0,len(names),namelength)]
        self.obj_template._fieldnames = [x.tostring().strip('\x00')
                                        for x in splitnames]
        return super(Mat5StructMatrixGetter, self).get_raw_array()

    def get_item(self):
        item = pycopy(self.obj_template)
        for element in item._fieldnames:
            item.__dict__[element]  = self.read_element()
        return item


class Mat5ObjectMatrixGetter(Mat5StructMatrixGetter):
    def __init__(self, *args, **kwargs):
        super(Mat5StructMatrixGetter, self).__init__(*args, **kwargs)
        self.obj_template = mat_obj()

    def get_raw_array(self):
        self.obj_template._classname = self.read_element().tostring()
        return super(Mat5ObjectMatrixGetter, self).get_raw_array()


class MatFile5Reader(MatFileReader):
    ''' Reader for Mat 5 mat files

    Adds the following attribute to base class
    
    @uint16_codec       - char codec to use for uint16 char arrays
                          (defaults to system default codec)
   '''

    def __init__(self,
                 mat_stream,
                 byte_order=None,
                 base_name='raw',
                 matlab_compatible=False,
                 squeeze_me=True,
                 chars_as_strings=True,
                 uint16_codec=None
                 ):
        self.codecs = {}
        self._array_reader = Mat5ArrayReader(
            mat_stream,
            None,
            None,
            None,
            None,
            )
        super(MatFile5Reader, self).__init__(
            mat_stream,
            byte_order,
            base_name,
            matlab_compatible,
            squeeze_me,
            chars_as_strings)
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
        self.mat_stream.seek(126)
        mi = self.mat_stream.read(2)
        self.mat_stream.seek(0)
        return mi == 'IM' and '<' or '>'

    def file_header(self):
        ''' Read in mat 5 file header '''
        hdict = {}
        hdr = self.read_dtype(self.dtypes['file_header'])
        hdict['__header__'] = hdr['description'].strip(' \t\n\000')
        v_major = hdr['version'] >> 8
        v_minor = hdr['version'] & 0xFF
        hdict['__version__'] = '%d.%d' % (v_major, v_minor)
        return hdict
        
    def format_looks_right(self):
        # Mat4 files have a zero somewhere in first 4 bytes
        self.mat_stream.seek(0)
        mopt_bytes = ndarray(shape=(4,),
                             dtype=uint8,
                             buffer = self.mat_stream.read(4))
        self.mat_stream.seek(0)
        return 0 not in mopt_bytes


class Mat5Writer(MatFileWriter):
    pass
