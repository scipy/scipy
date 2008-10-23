## Automatically adapted for scipy Oct 05, 2005 by convertcode.py

# Author: Travis Oliphant

import struct
import sys
import types

from numpy import *
import numpyio

import warnings
warnings.warn('fopen module is deprecated, please use npfile instead',
              DeprecationWarning, stacklevel=2)

LittleEndian = (sys.byteorder == 'little')

__all__ = ['fopen']

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
    elif mtype in ['u4','int32','integer*4']:
        mtype = 'u4'
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
        return self.file.tell()

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
                   int32           : 'u4', 'int32', 'integer*4'
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
        count = product(data.shape,axis=0)
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
                count = product(shape,axis=0)
            elif minus_ones == 1:
                now = self.tell()
                self.seek(0,2)
                end = self.tell()
                self.seek(now)
                remaining_bytes = end - now
                know_dimensions_size = -product(count,axis=0) * getsize_type(stype)[0]
                unknown_dimension_size, illegal = divmod(remaining_bytes,
                                                         know_dimensions_size)
                if illegal:
                    raise ValueError("unknown dimension doesn't match filesize")
                shape[shape.index(-1)] = unknown_dimension_size
                count = product(shape,axis=0)
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
            count = product(fmt.shape,axis=0)
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
