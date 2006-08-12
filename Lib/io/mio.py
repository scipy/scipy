## Automatically adapted for scipy Oct 05, 2005 by convertcode.py

# Author: Travis Oliphant

import struct, os, sys
import types
from tempfile import mkstemp
import zlib

from numpy import array, asarray, empty, obj2sctype, product, reshape, \
    squeeze, transpose, zeros, vstack, ndarray, shape, diff, where, uint8
import numpyio

try:
    import scipy.sparse
    have_sparse = 1
except ImportError:
    have_sparse = 0

LittleEndian = (sys.byteorder == 'little')

_unit_imag = {'f': array(1j,'F'), 'd': 1j}

__all__ = ['fopen','loadmat','savemat']

def getsize_type(mtype):
    if mtype in ['B','uchar','byte','unsigned char','integer*1', 'int8']:
        mtype = 'B'
    elif mtype in ['S1', 'char', 'char*1']:
        mtype = 'B'
    elif mtype in ['b', 'schar', 'signed char']:
        mtype = 'b'
    elif mtype in ['h','short','int16','integer*2']:
        mtype = 'h'
    elif mtype in ['H','ushort','uint16','unsigned short']:
        mtype = 'H'
    elif mtype in ['i','int']:
        mtype = 'i'
    elif mtype in ['I','uint','uint32','unsigned int']:
        mtype = 'I'
    elif mtype in ['l','long','int32','integer*4']:
        mtype = 'l'
    elif mtype in ['f','float','float32','real*4', 'real']:
        mtype = 'f'
    elif mtype in ['d','double','float64','real*8', 'double precision']:
        mtype = 'd'
    elif mtype in ['F','complex float','complex*8','complex64']:
        mtype = 'F'
    elif mtype in ['D','complex*16','complex128','complex','complex double']:
        mtype = 'D'
    else:
        mtype = obj2sctype(mtype)

    newarr = empty((1,),mtype)
    return newarr.itemsize, newarr.dtype.char

class fopen(object):
    """Class for reading and writing binary files into numpy arrays.

    Inputs:

      file_name -- The complete path name to the file to open.
      permission -- Open the file with given permissions: ('r', 'H', 'a')
                    for reading, writing, or appending.  This is the same
                    as the mode argument in the builtin open command.
      format -- The byte-ordering of the file:
                (['native', 'n'], ['ieee-le', 'l'], ['ieee-be', 'B']) for
                native, little-endian, or big-endian respectively.

    Attributes (Read only):

      bs -- non-zero if byte-swapping is performed on read and write.
      format -- 'native', 'ieee-le', or 'ieee-be'
      closed -- non-zero if the file is closed.
      mode -- permissions with which this file was opened
      name -- name of the file
    """

#    Methods:
#
#      read -- read data from file and return numpy array
#      write -- write to file from numpy array
#      fort_read -- read Fortran-formatted binary data from the file.
#      fort_write -- write Fortran-formatted binary data to the file.
#      rewind -- rewind to beginning of file
#      size -- get size of file
#      seek -- seek to some position in the file
#      tell -- return current position in file
#      close -- close the file

    def __init__(self,file_name,permission='rb',format='n'):
        if 'b' not in permission: permission += 'b'
        if isinstance(file_name, basestring):
            self.file = file(file_name, permission)
        elif isinstance(file_name, file) and not file_name.closed:
            # first argument is an open file
            self.file = file_name
        else:
            raise TypeError, 'Need filename or open file as input'
        self.setformat(format)
        self.zbuffer = None
        
    def __del__(self):
        try:
            self.file.close()
        except:
            pass

    def close(self):
        self.file.close()

    def seek(self, *args):
        self.file.seek(*args)

    def tell(self):
        self.file.tell()
        
    def raw_read(self, size=-1):
        """Read raw bytes from file as string."""
        return self.file.read(size)

    def raw_write(self, str):
        """Write string to file as raw bytes."""
        return self.file.write(str)

    def setformat(self, format):
        """Set the byte-order of the file."""
        if format in ['native','n','default']:
            self.bs = False
            self.format = 'native'
        elif format in ['ieee-le','l','little-endian','le']:
            self.bs = not LittleEndian
            self.format = 'ieee-le'
        elif format in ['ieee-be','B','big-endian','be']:
            self.bs = LittleEndian
            self.format = 'ieee-be'
        else:
            raise ValueError, "Unrecognized format: " + format
        return

    def write(self,data,mtype=None,bs=None):
        """Write to open file object the flattened numpy array data.

        Inputs:

          data -- the numpy array to write.
          mtype -- a string indicating the binary type to write.
                   The default is the type of data. If necessary a cast is made.
                   unsigned byte  : 'B', 'uchar', 'byte' 'unsigned char', 'int8',
                                    'integer*1'
                   character      : 'S1', 'char', 'char*1'
                   signed char    : 'b', 'schar', 'signed char'
                   short          : 'h', 'short', 'int16', 'integer*2'
                   unsigned short : 'H', 'ushort','uint16','unsigned short'
                   int            : 'i', 'int'
                   unsigned int   : 'I', 'uint32','uint','unsigned int'
                   long           : 'l', 'long', 'int32', 'integer*4'
                   float          : 'f', 'float', 'float32', 'real*4'
                   double         : 'd', 'double', 'float64', 'real*8'
                   complex float  : 'F', 'complex float', 'complex*8', 'complex64'
                   complex double : 'D', 'complex', 'complex double', 'complex*16',
                                    'complex128'
        """
        if bs is None:
            bs = self.bs
        else:
            bs = (bs == 1)
        if isinstance(data, str):
            N, buf = len(data), buffer(data)
            data = ndarray(shape=(N,),dtype='B',buffer=buf)
        else:
            data = asarray(data)
        if mtype is None:
            mtype = data.dtype.char
        howmany,mtype = getsize_type(mtype)
        count = product(data.shape)
        numpyio.fwrite(self.file,count,data,mtype,bs)
        return

    fwrite = write

    def read(self,count,stype,rtype=None,bs=None,c_is_b=0):
        """Read data from file and return it in a numpy array.

        Inputs:

          count -- an integer specifying the number of elements of type
                   stype to read or a tuple indicating the shape of
                   the output array.
          stype -- The data type of the stored data (see fwrite method).
          rtype -- The type of the output array.  Same as stype if None.
          bs -- Whether or not to byteswap (or use self.bs if None)
          c_is_b --- If non-zero then the count is an integer
                   specifying the total number of bytes to read
                   (must be a multiple of the size of stype).

        Outputs: (output,)

          output -- a numpy array of type rtype.
        """
        if bs is None:
            bs = self.bs
        else:
            bs = (bs == 1)
        howmany,stype = getsize_type(stype)
        shape = None
        if c_is_b:
            if count % howmany != 0:
                raise ValueError, "When c_is_b is non-zero then " \
                      "count is bytes\nand must be multiple of basic size."
            count = count / howmany
        elif type(count) in [types.TupleType, types.ListType]:
            shape = list(count)
            # allow -1 to specify unknown dimension size as in reshape
            minus_ones = shape.count(-1)
            if minus_ones == 0:
                count = product(shape)
            elif minus_ones == 1:
                now = self.tell()
                self.seek(0,2)
                end = self.tell()
                self.seek(now)
                remaining_bytes = end - now
                know_dimensions_size = -product(count) * getsize_type(stype)[0]
                unknown_dimension_size, illegal = divmod(remaining_bytes,
                                                         know_dimensions_size)
                if illegal:
                    raise ValueError("unknown dimension doesn't match filesize")
                shape[shape.index(-1)] = unknown_dimension_size
                count = product(shape)
            else:
                raise ValueError(
                    "illegal count; can only specify one unknown dimension")
            shape = tuple(shape)
        if rtype is None:
            rtype = stype
        else:
            howmany,rtype = getsize_type(rtype)
        if count == 0:
            return zeros(0,rtype)
        retval = numpyio.fread(self.file, count, stype, rtype, bs)
        if shape is not None:
            retval = resize(retval, shape)
        return retval

    fread = read

    def rewind(self,howmany=None):
        """Rewind a file to its beginning or by a specified amount.
        """
        if howmany is None:
            self.seek(0)
        else:
            self.seek(-howmany,1)

    def size(self):
        """Return the size of the file.
        """
        try:
            sz = self.thesize
        except AttributeError:
            curpos = self.tell()
            self.seek(0,2)
            sz = self.tell()
            self.seek(curpos)
            self.thesize = sz
        return sz

    def fort_write(self,fmt,*args):
        """Write a Fortran binary record.

        Inputs:

          fmt -- If a string then it represents the same format string as
                 used by struct.pack.  The remaining arguments are passed
                 to struct.pack.

                 If fmt is an array, then this array will be written as
                 a Fortran record using the output type args[0].

          *args -- Arguments representing data to write.
        """
        if self.format == 'ieee-le':
            nfmt = "<i"
        elif self.format == 'ieee-be':
            nfmt = ">i"
        else:
            nfmt = "i"
        if isinstance(fmt, basestring):
            if self.format == 'ieee-le':
                fmt = "<"+fmt
            elif self.format == 'ieee-be':
                fmt = ">"+fmt
            str = apply(struct.pack,(fmt,)+args)
            strlen = struct.pack(nfmt,len(str))
            self.write(strlen)
            self.write(str)
            self.write(strlen)
        elif type(fmt) == type(array([0])):
            if len(args) > 0:
                sz,mtype = getsize_type(args[0])
            else:
                sz,mtype = getsize_type(fmt.dtype.char)
            count = product(fmt.shape)
            strlen = struct.pack(nfmt,count*sz)
            self.write(strlen)
            numpyio.fwrite(self.file,count,fmt,mtype,self.bs)
            self.write(strlen)
        else:
            raise TypeError, "Unknown type in first argument"

    def fort_read(self,fmt,dtype=None):
        """Read a Fortran binary record.

        Inputs:

          fmt -- If dtype is not given this represents a struct.pack
                 format string to interpret the next record.  Otherwise this
                 argument is ignored.
          dtype -- If dtype is not None, then read in the next record as
                   an array of type dtype.

        Outputs: (data,)

          data -- If dtype is None, then data is a tuple containing the output
                  of struct.unpack on the next Fortan record.
                  If dtype is a datatype string, then the next record is
                  read in as a 1-D array of type datatype.
        """
        lookup_dict = {'ieee-le':"<",'ieee-be':">",'native':''}
        if dtype is None:
            fmt = lookup_dict[self.format] + fmt
            numbytes = struct.calcsize(fmt)
            nn = struct.calcsize("i");
            if (self.raw_read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            strdata = self.raw_read(numbytes)
            if strdata == '':
                raise ValueError, "Unexpected end of file..."
            data = struct.unpack(fmt,strdata)
            if (self.raw_read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            return data
        else:  # Ignore format string and read in next record as an array.
            fmt = lookup_dict[self.format] + "i"
            nn = struct.calcsize(fmt)
            nbytestr = self.raw_read(nn)
            if nbytestr == '':
                raise ValueError, "Unexpected end of file..."
            nbytes = struct.unpack(fmt,nbytestr)[0]
            howmany, dtype = getsize_type(dtype)
            ncount = nbytes / howmany
            if ncount*howmany != nbytes:
                self.rewind(4)
                raise ValueError, "A mismatch between the type requested and the data stored."
            if ncount < 0:
                raise ValueError, "Negative number of bytes to read:\n    file is probably not opened with correct endian-ness."
            if ncount == 0:
                raise ValueError, "End of file?  Zero-bytes to read."
            retval = numpyio.fread(self.file, ncount, dtype, dtype, self.bs)
            if len(retval) == 1:
                retval = retval[0]
            if (self.raw_read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            return retval
        

class CompressedFopen(fopen):
    """ File container for temporary buffer to decompress data """
    def __init__(self, *args, **kwargs):
        fd, fname = mkstemp()
        super(CompressedFopen, self).__init__(
            os.fdopen(fd, 'w+b'), *args, **kwargs)
        self.file_name = fname
        
    def fill(self, bytes):
        """ Uncompress buffer in @bytes and write to file """
        self.rewind()
        self.raw_write(zlib.decompress(bytes))
        self.rewind()

    def __del__(self):
        try:
            self.file.truncate(0)
        except:
            pass
        try:
            self.close()
        except:
            pass
        try:
            os.remove(self.file_name)
        except:
            pass
        
#### MATLAB Version 5 Support ###########

# Portions of code borrowed and (heavily) adapted
#    from matfile.py by Heiko Henkelmann

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


class mat_struct:    # dummy structure holder
    pass

class mat_obj:    # dummy object holder
    pass

miINT8 = 1
miUINT8 = 2
miINT16 = 3
miUINT16 = 4
miINT32 = 5
miUINT32 = 6
miSINGLE = 7
miDOUBLE = 9
miINT64 =12
miUINT64 = 13
miMATRIX = 14
miCOMPRESSED = 15
miUTF8 = 16
miUTF16 = 17
miUTF32 = 18

miNumbers = (
    miINT8,
    miUINT8,
    miINT16,
    miUINT16,
    miINT32,
    miUINT32,
    miSINGLE,
    miDOUBLE,
    miINT64,
    miUINT64,
    )

miDataTypes = {
    miINT8 : ('miINT8', 1,'b'),
    miUINT8 : ('miUINT8', 1,'B'),
    miINT16 : ('miINT16', 2,'h'),
    miUINT16 :('miUINT16',2,'H'),
    miINT32 : ('miINT32',4,'l'),
    miUINT32 : ('miUINT32',4,'I'),
    miSINGLE : ('miSINGLE',4,'f'),
    miDOUBLE : ('miDOUBLE',8,'d'),
    miINT64 : ('miINT64',8,'q'),
    miUINT64 : ('miUINT64',8,'Q'),
    miMATRIX : ('miMATRIX',0,None),
    miUTF8 : ('miUTF8',1,'b'),
    miUTF16 : ('miUTF16',2,'h'),
    miUTF32 : ('miUTF32',4,'l'),
    }

''' Before release v7.1 (release 14) matlab used the system default
character encoding scheme padded out to 16-bits. Release 14 and later
use Unicode. When saving character data, matlab R14 checks if it can
be encoded in 7-bit ascii, and saves in that format if so.'''
miCodecs = {
    miUINT8: 'ascii',
    miUINT16: sys.getdefaultencoding(),
    miUTF8: 'utf8',
    miUTF16: 'utf16',
    miUTF32: 'utf32',
    } 

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

mxArrays = (
    mxCHAR_CLASS,
    mxDOUBLE_CLASS,
    mxSINGLE_CLASS,
    mxINT8_CLASS,
    mxUINT8_CLASS,
    mxINT16_CLASS,
    mxUINT16_CLASS,
    mxINT32_CLASS,
    mxUINT32_CLASS,
    )
    
def _parse_header(fid, hdict):
    correct_endian = (ord('M')<<8) + ord('I')
                 # if this number is read no BS
    fid.seek(126)  # skip to endian detector
    endian_test = fid.read(1,'int16')
    if (endian_test == correct_endian): openstr = 'n'
    else:  # must byteswap
        if LittleEndian:
            openstr = 'B'
        else: openstr = 'l'
    fid.setformat(openstr)  # change byte-order if necessary
    fid.rewind()
    hdict['__header__'] = fid.raw_read(124).strip(' \t\n\000')
    vers = fid.read(1,'int16')
    hdict['__version__'] = '%d.%d' % (vers >> 8, vers & 0xFF)
    fid.seek(2,1)  # move to start of data
    return

def _skip_padding(fid, numbytes, rowsize):
    """ Skip to next row or @rowsize after previous read of @numbytes """
    mod = numbytes % rowsize
    if mod:
        skip = rowsize-mod
        fid.seek(skip,1)

def _parse_array_flags(fid):
    # first 8 bytes are always miUINT32 and 8 --- just a check
    dtype, nbytes = fid.read(2,'I')
    if (dtype != miUINT32) or (nbytes != 8):
        raise IOError, "Invalid MAT file. Perhaps a byte-order problem."

    # read array flags.
    rawflags = fid.read(2,'I')
    class_ = rawflags[0] & 0xFF
    flags = (rawflags[0] & 0xFFFF) >> 8
    # Global and logical fields are currently ignored
    if (flags & 8): cmplx = 1
    else: cmplx = 0
    if class_ == mxSPARSE_CLASS:
        nzmax = rawflags[1]
    else:
        nzmax = None
    return class_, cmplx, nzmax

def _parse_mimatrix(fid,bytes):
    dclass, cmplx, nzmax =_parse_array_flags(fid)
    dims = _get_element(fid)
    name = _get_element(fid).tostring()
    tupdims = tuple(dims[::-1])
    if dclass in mxArrays:
        result, unused, dtype =_get_element(fid, return_name_dtype=True)
        if dclass == mxCHAR_CLASS:
            en = miCodecs[dtype]
            try:
                " ".encode(en)
            except LookupError:
                raise ValueError, 'Character encoding %s not supported' % en
            if dtype == miUINT16:
                char_len = len("  ".encode(en)) - len(" ".encode(en))
                if char_len == 1: # Need to downsample from 16 bit
                    result = result.astype(uint8)
            result = squeeze(transpose(reshape(result,tupdims)))
            dims = result.shape
            if len(dims) >= 2: # return array of strings
                n_dims = dims[:-1]
                string_arr = reshape(result, (product(n_dims), dims[-1]))
                result = empty(n_dims, dtype=object)
                for i in range(0, n_dims[-1]):
                    result[...,i] = string_arr[i].tostring().decode(en)
            else: # return string
                result = result.tostring().decode(en)
        else:
            if cmplx:
                imag  =_get_element(fid)
                try:
                    result = result + _unit_imag[imag.dtype.char] * imag
                except KeyError:
                    result = result + 1j*imag
            result = squeeze(transpose(reshape(result,tupdims)))

    elif dclass == mxCELL_CLASS:
        length = product(dims)
        result = empty(length, dtype=object)
        for i in range(length):
            result[i] = _get_element(fid)
        result = squeeze(transpose(reshape(result,tupdims)))
        if not result.shape:
            result = result.item()

    elif dclass == mxSTRUCT_CLASS:
        length = product(dims)
        result = zeros(length, object)
        namelength = _get_element(fid)
        # get field names
        names = _get_element(fid)
        splitnames = [names[i:i+namelength] for i in \
                      xrange(0,len(names),namelength)]
        fieldnames = [x.tostring().strip('\x00')
                              for x in splitnames]
        for i in range(length):
            result[i] = mat_struct()
            for element in fieldnames:
                result[i].__dict__[element]  = _get_element(fid)
        result = squeeze(transpose(reshape(result,tupdims)))
        if not result.shape:
            result = result.item()

        # object is like a structure with but with a class name
    elif dclass == mxOBJECT_CLASS:
        class_name = _get_element(fid).tostring()
        length = product(dims)
        result = zeros(length, object)
        namelength = _get_element(fid)
        # get field names
        names = _get_element(fid)
        splitnames = [names[i:i+namelength] for i in \
                      xrange(0,len(names),namelength)]
        fieldnames = [x.tostring().strip('\x00')
                              for x in splitnames]
        for i in range(length):
            result[i] = mat_obj()
            result[i]._classname = class_name
            for element in fieldnames:
                result[i].__dict__[element] = _get_element(fid)
        result = squeeze(transpose(reshape(result,tupdims)))
        if not result.shape:
            result = result.item()

    elif dclass == mxSPARSE_CLASS:
        rowind  = _get_element(fid)
        colind = _get_element(fid)
        res = _get_element(fid)
        if cmplx:
            imag = _get_element(fid)
            try:
                res = res + _unit_imag[imag.dtype.char] * imag
            except (KeyError,AttributeError):
                res = res + 1j*imag
        ''' From the matlab API documentation, last found here:
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
            result = scipy.sparse.csc_matrix((res,ij), [dims[0],dims[1]])
        else:
            result = (dims, ij, res)

    return result, name

# Return a Python object for the element
def _get_element(fid, return_name_dtype=False):
    """ Return a python object from next element in @fid

    @fid    - fopen object for matfile
    @return_name_dtype - if True, return tuple of (element, name, dtype)
                         if False, return element only
    """
    name = None
    test = fid.raw_read(1)
    if len(test) == 0:  # nothing left
        raise EOFError
    else:
        fid.rewind(1)
    # get the data tag
    raw_tag = int(fid.read(1,'I'))
    
    # check for small data element format
    numbytes = raw_tag >> 16
    if numbytes > 0:  # small data element format
        if numbytes > 4:
            raise IOError, "Problem with MAT file: " \
                  "too many bytes in small data element format."
        dtype = int(raw_tag & 0xFFFF)
        el = fid.read(numbytes,miDataTypes[dtype][2],c_is_b=1)
        fid.seek(4-numbytes,1)  # skip padding
    else:
        # otherwise parse tag
        dtype = raw_tag
        numbytes = fid.read(1,'I')
        
        if dtype == miCOMPRESSED: # compressed data type
            if not fid.zbuffer:
                fid.zbuffer = CompressedFopen(format=fid.format)
            fid.zbuffer.fill(fid.raw_read(numbytes))
            _skip_padding(fid, numbytes, 8)
            return _get_element(fid.zbuffer, return_name_dtype)
        if dtype != miMATRIX:  # basic data type
            try:
                el = fid.read(numbytes,miDataTypes[dtype][2],c_is_b=1)
            except KeyError:
                raise ValueError, "Unknown data type"
            _skip_padding(fid, numbytes, 8)
        else:
            # handle miMatrix type
            el, name = _parse_mimatrix(fid,numbytes)

    if return_name_dtype:
        return el, name, dtype
    return el

def _loadv5(fid,basename):
    # return a dictionary from a Matlab version 5-7.1 file
    # always contains the variable __header__
    mdict = {}
    _parse_header(fid,mdict)
    var = 0
    while 1:  # file pointer to start of next data
        try:
            var = var + 1
            el, varname, unused = _get_element(fid, return_name_dtype=True)
            if varname is None:
                varname = '%s_%04d' % (basename,var)
            mdict[varname] = el
        except EOFError:
            break
    return mdict

### END MATLAB v5 support #############

def loadmat(name, mdict=None, appendmat=1, basename='raw'):
    """Load the MATLAB(tm) mat file.

    If name is a full path name load it in.  Otherwise search for the file
    on the sys.path list and load the first one found (the current directory
    is searched first).

    v4 (Level 1.0), v6 and v7.1 matfiles are supported.  

    Inputs:

      name -- name of the mat file (don't need .mat extension if appendmat=1)
      dict -- the dictionary to insert into.  If none the variables will be
              returned in a dictionary.
      appendmat -- non-zero to append the .mat extension to the end of the
                   given filename.
      basename -- for MATLAB(tm) v5 matfiles raw data will have this basename.

    Outputs:

      If dict is None, then a dictionary of names and objects representing the
      stored arrays is returned.
    """

    if appendmat and name[-4:] == ".mat":
        name = name[:-4]
    if os.sep in name:
        full_name = name
        if appendmat:
            full_name = name + ".mat"
    else:
        full_name = None
        junk,name = os.path.split(name)
        for path in sys.path:
            test_name = os.path.join(path,name)
            if appendmat:
                test_name += ".mat"
            try:
                fid = open(test_name,'rb')
                fid.close()
                full_name = test_name
                break
            except IOError:
                pass
        if full_name is None:
            raise IOError, "%s not found on the path." % name

    fid = fopen(full_name,'rb')
    test_vals = fid.fread(4,'byte')

    if not (0 in test_vals):       # MATLAB version 5 format
        fid.rewind()
        thisdict = _loadv5(fid,basename)
        if mdict is not None:
            mdict.update(thisdict)
            return
        else:
            return thisdict
        
    # The remainder of this function is the v4 codepath
    testtype = struct.unpack('i',test_vals.tostring())
    # Check to see if the number is positive and less than 5000.
    if testtype[0] < 0 or testtype[0] > 4999:
        # wrong byte-order
        if LittleEndian:
            format = 'ieee-be'
        else:
            format = 'ieee-le'
    else:  # otherwise we are O.K.
        if LittleEndian:
            format = 'ieee-le'
        else:
            format = 'ieee-be'

    fid.setformat(format)

    length = fid.size()
    fid.rewind()  # back to the begining

    defnames = []
    thisdict = {}
    while 1:
        if (fid.tell() == length):
            break
        header = fid.fread(5,'int')
        if len(header) != 5:
            fid.close()
            print "Warning: Read error in file."
            break
        M,rest = divmod(int(header[0]),1000) # int is for workaround numpy 0.9.9 bug
        O,rest = divmod(rest,100)
        P,rest = divmod(rest,10)
        T = rest

        if (M > 1):
            fid.close()
            raise ValueError, "Unsupported binary format."
        if (O != 0):
            fid.close()
            raise ValuError, "Hundreds digit of first integer should be zero."

        if (T not in [0,1]):
            fid.close()
            raise ValueError, "Cannot handle sparse matrices, yet."

        storage = {0:'d',1:'f',2:'i',3:'h',4:'H',5:'B'}[P]

        varname = fid.fread(header[-1],'char')[:-1]
        varname = varname.tostring()
        defnames.append(varname)
        numels = header[1]*header[2]
        if T == 0:             # Text data
            data = atleast_1d(fid.fread(numels,storage))
            if header[3]:  # imaginary data
                data2 = fid.fread(numels,storage)
                if data.dtype.char == 'f' and data2.dtype.char == 'f':
                    new = empty(data.shape,'F')
                    new.real = data
                    new.imag = data2
                    data = new
                    del(new)
                    del(data2)
            if len(data) > 1:
                data=data.reshape((header[2], header[1])                )
                thisdict[varname] = transpose(squeeze(data))
            else:
                thisdict[varname] = data
        else:
            data = atleast_1d(fid.fread(numels,storage,'char'))
            if len(data) > 1:
                data=data.reshape((header[2], header[1]))
                thisdict[varname] = transpose(squeeze(data))
            else:
                thisdict[varname] = data

    fid.close()
    if mdict is not None:
        print "Names defined = ", defnames
        mdict.update(thisdict)
    else:
        return thisdict


def savemat(filename, mdict):
    """Save a dictionary of names and arrays into the MATLAB-style .mat file.

    This saves the arrayobjects in the given dictionary to a matlab Version 4
    style .mat file.
    """
    storage = {'D':0,'d':0,'F':1,'f':1,'l':2,'i':2,'h':3,'B':5}
    if filename[-4:] != ".mat":
        filename = filename + ".mat"
    fid = fopen(filename,'wb')
    M = not LittleEndian
    O = 0
    for variable in mdict.keys():
        var = mdict[variable]
        if not isinstance(var, ndarray):
            continue
        if var.dtype.char == 'S1':
            T = 1
        else:
            T = 0
        if var.dtype.char == 'b':
            var = var.astype('h')
        P = storage[var.dtype.char]
        fid.fwrite([M*1000+O*100+P*10+T],'int')

        if len(var.shape) == 1:
            var=var.reshape((len(var), 1))
        var = transpose(var)

        if len(var.shape) > 2:
            var=var.reshape((product(var.shape[:-1]), var.shape[-1]))

        imagf = var.dtype.char in ['F', 'D']
        fid.fwrite([var.shape[1], var.shape[0], imagf, len(variable)+1],'int')
        fid.fwrite(variable+'\x00','char')
        if imagf:
            fid.fwrite(var.real)
            fid.fwrite(var.imag)
        else:
            fid.fwrite(var)
    fid.close()
    return
