import numpyio
import struct
from scipy.numeric import *

def getsize_type(mtype):
    if mtype in ['b','uchar','byte','unsigned char','integer*1', 'int8']:
        mtype = 'b'
    elif mtype in ['c', 'char','char*1']:
        mtype = 'c'
    elif mtype in ['1','schar', 'signed char']:
        mtype = '1'
    elif mtype in ['s','short','int16','integer*2']:
        mtype = 's'
    elif mtype in ['i','int']:
        mtype = 'i'
    elif mtype in ['l','long','int32','integer*4']:
        mtype = 'l'
    elif mtype in ['f','float','float32','real*4']:
        mtype = 'f'
    elif mtype in ['d','double','float64','real*8']:
        mtype = 'd'
    elif mtype in ['F','complex float','complex*8','complex64']:
        mtype = 'F'
    elif mtype in ['D','complex*16','complex128','complex','complex double']:
        mtype = 'D'
    else:
        raise TypeError, 'Bad datatype -- ' + mtype

    argout = (array(0,mtype).itemsize(),mtype)
    return argout

class fopen:
    """Class for reading and writing binary files into Numeric arrays.

    Inputs:

      file_name -- The complete path name to the file to open.
      permission -- Open the file with given permissions: ('r', 'w', 'a')
                    for reading, writing, or appending.  This is the same
                    as the mode argument in the builtin open command.
      format -- The byte-ordering of the file:
                (['native', 'n'], ['ieee-le', 'l'], ['ieee-be', 'b']) for
                native, little-endian, or big-endian respectively.

    Methods:

      fread -- read data from file and return Numeric array
      fwrite -- write to file from scipy.numeric array
      fort_read -- read Fortran-formatted binary data from the file.
      fort_write -- write Fortran-formatted binary data to the file.
      rewind -- rewind to beginning of file
      size -- get size of file
      seek -- seek to some position in the file
      tell -- return current position in file
      close -- close the file

    Attributes (Read only):

      bs -- non-zero if byte-swapping is performed on read and write.
      format -- 'native', 'ieee-le', or 'ieee-be'
      fid -- the file object
      closed -- non-zero if the file is closed.
      mode -- permissions with which this file was opened
      name -- name of the file
      
    """
    
    def __init__(self,file_name,permission='r',format='n'):
        if type(file_name) == type(''):
            self.__dict__['fid'] = open(file_name,permission)
        elif 'fileno' in file_name.__methods__:  # first argument is an open file
            self.__dict__['fid'] = file_name 
        if format in ['native','n']:
            self.__dict__['bs'] = 0
            self.__dict__['format'] = 'native'
        elif format in ['ieee-le','l']:
            self.__dict__['bs'] = not LittleEndian
            self.__dict__['format'] = 'ieee-le'
        elif format in ['ieee-be','b']:
            self.__dict__['bs'] = LittleEndian
            self.__dict__['format'] = 'ieee-be'

        self.__dict__['seek'] = self.fid.seek
        self.__dict__['tell']= self.fid.tell
        self.__dict__['close'] = self.fid.close
        self.__dict__['fileno'] = self.fid.fileno
        self.__dict__['mode'] = self.fid.mode
        self.__dict__['closed'] = self.fid.closed
        self.__dict__['name'] = self.fid.name

    def __setattr__(self, attribute):
        raise SyntaxError, "There are no user-settable attributes."            

    def __del__(self):
        try:
            self.fid.close()
        except:
            pass
    
    def fwrite(self,data,mtype=None):
        """Write to open file object the flattened Numeric array data.

        Inputs:

          data -- the Numeric array to write.
          mtype -- a string indicating the binary type to write.
                   The default is the type of data. If necessary a cast is made.
                   unsigned byte  : 'b', 'uchar', 'byte' 'unsigned char', 'int8',
                                    'integer*1'
                   character      : 'c', 'char', 'char*1'
                   signed char    : '1', 'schar', 'signed char'
                   short          : 's', 'short', 'int16', 'integer*2'
                   int            : 'i', 'int'
                   long           : 'l', 'long', 'int32', 'integer*4'
                   float          : 'f', 'float', 'float32', 'real*4'
                   double         : 'd', 'double', 'float64', 'real*8'
                   complex float  : 'F', 'complex float', 'complex*8', 'complex64'
                   complex double : 'D', 'complex', 'complex double', 'complex*16',
                                    'complex128'

        Outputs: (val,)

          val -- The number of bytes written.
        """
        data = asarray(data)
        if mtype is None:
            mtype = data.typecode()
        howmany,mtype = getsize_type(mtype)
        count = product(data.shape)
        val = numpyio.fwrite(self.fid,count,data,mtype,self.bs)
        return val

    def fread(self,count,stype,rtype=None):
        """Read data from file and return it in a Numeric array.

        Inputs:

          count -- an integer specifying the number of bytes to read or
                   a tuple indicating the shape of the output array.
          stype -- The data type of the stored data (see fwrite method).
          rtype -- The type of the output array.  Same as stype if None.

        Outputs: (output,)

          output -- a Numeric array of type rtype.
        """
        shape = None
        if type(count) in [types.Tupletype, types.ListType]:
            shape = tuple(count)
            count = product(shape)
        howmany,stype = getsize_type(stype)
        if rtype is None:
            rtype = stype
        else:
            howmany,rtype = getsize_type(rtype)
        retval = numpyio.fread(self.fid, count, stype, rtype, self.bs)
        if len(retval) == 1:
            retval = retval[0]
        if shape is not None:
            retval = resize(retval, shape)
        return retval

    def rewind(self,howmany=None):
        """Rewind a file to it's beginning or by a specified amount.
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
            self.fid.seek(0,2)
            sz = self.fid.tell()
            self.fid.seek(curpos)
            self.__dict__['thesize'] = sz
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
        if type(fmt) == type(''):
            if self.format == 'ieee-le':
                fmt = "<"+fmt
            elif self.format == 'ieee-be':
                fmt = ">"+fmt
            str = apply(struct.pack,(fmt,)+args)
            strlen = struct.pack(nfmt,len(str))
            self.fid.write(strlen)
            self.fid.write(str)
            self.fid.write(strlen)
        elif type(fmt) == type(array([0])):
            sz,mtype = getsize_type(args[0])
            count = product(fmt.shape)
            strlen = struct.pack(nfmt,count*sz)
            self.fid.write(strlen)
            numpyio.fwrite(self.fid,count,fmt,mtype,self.bs)
            self.fid.write(strlen)
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
            self.fid.read(nn)
            data = struct.unpack(fmt,self.fid.read(numbytes))
            self.fid.read(nn)
            return data
        else:  # Ignore format string and read in next record as an array.
            fmt = lookup_dict[self.format] + "i"
            nn = struct.calcsize(fmt)
            nbytes = struct.unpack(fmt,self.fid.read(nn))[0]
            howmany, dtype = getsize_type(dtype)
            ncount = nbytes / howmany
            if ncount*howmany != nbytes:
                self.rewind(4)
                raise ValueError, "A mismatch between the type requested and the data stored."
            retval = numpyio.fread(self.fid, ncount, dtype, dtype, self.bs)
            if len(retval) == 1:
                retval = retval[0]
            self.fid.read(nn)
            return retval
                                      












