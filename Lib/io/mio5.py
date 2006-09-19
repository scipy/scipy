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

from numpy import *

from bytestream import ByteStream
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
    'tag_mdtype': 'u4',
    'tag_byte_count': 'u4',
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

class Mat5Tag(object):
    ''' Placeholder for holding tag information '''
    pass

class Mat5Header(object):
    ''' Placeholder for Mat5 header

    Defines:
    next_position - start position of next matrix
    name
    dtype - numpy dtype of matrix
    mclass - matlab (TM) code for class of matrix
    dims - shape of matrix as stored (see sparse reader)
    is_complex - True if data are complex
    is_char    - True if these are char data
    is_global  - is a global variable in matlab (TM) workspace
    is_numeric - is basic numeric matrix
    original_dtype - data type when saved from matlab (TM)
    '''
    def __init__(self):
        self.next_position = None
        self.is_empty = False
        self.flags = None
        self.is_complex = False
        self.is_global = False
        self.is_logical = False
        self.mclass = 0
        self.is_numeric = None
        self.original_dtype = None
        self.is_char = None
        self.dims = ()
        self.name = ''

class Mat5ArrayFlags(object):
    ''' Place holder for array flags '''
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

    def read_tag(self):
        tag = Mat5Tag()
        # Check for small data element first
        tag.mdtype = int(self.read_array(self.dtypes['tag_mdtype']))
        byte_count = tag.mdtype >> 16
        if byte_count: # small data element format
            if byte_count > 4:
                raise ValueError, 'Too many bytes for sde format'
            tag.byte_count = byte_count
            tag.mdtype = tag.mdtype & 0xFFFF
            tag.skip = 4 - byte_count
        else: # standard tag format
            tag.byte_count = self.read_array(
                self.dtypes['tag_byte_count'])
            tag.skip = tag.byte_count % 8 and 8 - tag.byte_count % 8
        return tag
    
    def read_element(self, copy=True):
        tag = self.read_tag()
        if tag.mdtype == miMATRIX:
            header = self.read_header(tag)
            return self.header_to_getter(header).get_array()
        if tag.mdtype in self.codecs: # encoded char data
           raw_str = self.read_bytes(tag.byte_count)
           codec = self.codecs[tag.mdtype]
           if not codec:
               raise TypeError, 'Do not support encoding %d' % tag.mdtype
           el = raw_str.tostring().decode(codec)
        else: # numeric data
            try:
                dt = self.dtypes[tag.mdtype]
            except KeyError:
                raise TypeError, 'Do not know matlab (TM) data code %d' \
                      % tag.mdtype
            el_count = tag.byte_count / dt.itemsize
            el = self.read_array(dt, a_shape=(el_count), copy=copy)
        if tag.skip:
            self.mat_stream.seek(tag.skip, 1)
        return el

    def read_header(self, tag):
        ''' Read header from Mat5 matrix
        '''
        if not tag.mdtype == miMATRIX:
            raise TypeError, \
                  'Expecting miMATRIX type here, got %d' %  tag.mdtype
        header = Mat5Header()
        # Note - there is no skip value for the miMATRIX type; the
        # number of bytes field in the tag points to the next variable
        # in the file. This is always aligned to 8 byte boundaries
        # (except for the miCOMPRESSED type)
        header.next_position = (self.mat_stream.pos +
                                tag.byte_count)
        # Apparently an empty miMATRIX can contain no bytes
        header.is_empty = tag.byte_count == 0
        if header.is_empty:
            return header
        header.flags = self.read_array_flags()
        header.is_complex = header.flags.is_complex
        header.is_global = header.flags.is_global
        header.is_logical = header.flags.is_logical
        header.mclass = header.flags.mclass
        header.dims = self.read_element()
        header.name = self.read_element().tostring()
        return header
    
    def read_array_flags(self):
        flags = Mat5ArrayFlags()
        af = self.read_array(self.dtypes['array_flags'])
        flags_class = af['flags_class']
        flags.mclass = flags_class & 0xFF
        flags.is_logical = flags_class >> 9 & 1
        flags.is_global = flags_class >> 10 & 1
        flags.is_complex = flags_class >> 11 & 1
        flags.nzmax = af['nzmax']
        return flags

    def matrix_getter_factory(self):
        ''' Returns reader for next matrix '''
        tag = self.read_tag()
        if tag.mdtype == miCOMPRESSED:
            return Mat5ZArrayReader(self, tag).matrix_getter_factory()
        header = self.read_header(tag)
        return self.header_to_getter(header)
    
    def header_to_getter(self, header):
        if header.is_empty:
            return Mat5EmptyMatrixGetter(self, header)
        mc = header.mclass
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
    def __init__(self, array_reader, tag):
        '''Reads and uncompresses gzipped stream'''
        data = array_reader.read_bytes(tag.byte_count)
        super(Mat5ZArrayReader, self).__init__(
            ByteStream(zlib.decompress(data.tostring())),
            array_reader.dtypes,
            array_reader.processor_func,
            array_reader.codecs,
            array_reader.class_dtypes)
        self.next_position = array_reader.mat_stream.tell()
        
    def header_to_getter(self, header):
        ''' Set next_position to current position in parent stream '''
        header.next_position = self.next_position
        return super(Mat5ZArrayReader, self).header_to_getter(header)
        

class Mat5MatrixGetter(MatMatrixGetter):
    ''' Base class for getting Mat5 matrices

    Gets current read information from passed array_reader
    '''
    
    def __init__(self, array_reader, header):
        ''' Accepts @array_reader and @header '''
        super(Mat5MatrixGetter, self).__init__(array_reader, header)
        self.class_dtypes = array_reader.class_dtypes
        self.codecs = array_reader.codecs
        self.is_global = header.is_global

    def read_tag(self):
        return self.array_reader.read_tag()
    
    def read_element(self, *args, **kwargs):
        return self.array_reader.read_element(*args, **kwargs)
    

class Mat5EmptyMatrixGetter(Mat5MatrixGetter):
    ''' Dummy class to return empty array for empty matrix '''
    def get_raw_array(self):
        return array([[]])


class Mat5NumericMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        self.header.is_numeric = True
        self.header.original_dtype = self.class_dtypes[self.header.mclass]
        if self.header.is_complex:
            # avoid array copy to save memory
            res = self.read_element(copy=False)
            res_j = self.read_element(copy=False)
            res = res + (res_j * 1j)
        else:
            res = self.read_element()
        return ndarray(shape=self.dims,
                       dtype=res.dtype,
                       buffer=res,
                       order='F')
        

class Mat5SparseMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        rowind  = self.read_element()
        colind = self.read_element()
        if self.header.is_complex:
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
                                             self.dims)
        else:
            result = (dims, ij, res)
        return result


class Mat5CharMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        self.header.is_char = True
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
        return ndarray(shape=self.dims,
                       dtype=dtype('U1'),
                       buffer=array(res),
                       order='F').copy()


class Mat5CellMatrixGetter(Mat5MatrixGetter):
    def get_raw_array(self):
        # Account for fortran indexing of cells
        tupdims = tuple(self.dims[::-1]) 
        length = product(self.dims)
        result = empty(length, dtype=object)
        for i in range(length):
            result[i] = self.get_item()
        result = transpose(reshape(result,tupdims))
        return result

    def get_item(self):
        return self.read_element()


class Mat5StructMatrixGetter(Mat5CellMatrixGetter):
    obj_template = mat_struct()
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
    obj_template = mat_obj()
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
        mi = self.read_bytes(2).tostring()
        self.mat_stream.seek(0)
        return mi == 'IM' and '<' or '>'

    def file_header(self):
        ''' Read in mat 5 file header '''
        hdict = {}
        hdr = self.read_array(self.dtypes['file_header'])
        hdict['__header__'] = hdr['description'].strip(' \t\n\000')
        v_major = hdr['version'] >> 8
        v_minor = hdr['version'] & 0xFF
        hdict['__version__'] = '%d.%d' % (v_major, v_minor)
        return hdict
        
    def format_looks_right(self):
        # Mat4 files have a zero somewhere in first 4 bytes
        self.mat_stream.seek(0)
        mopt_bytes = self.read_bytes(4)
        self.mat_stream.seek(0)
        return 0 not in mopt_bytes


class Mat5Writer(MatFileWriter):
    pass
