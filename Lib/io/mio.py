# Author: Travis Oliphant

from Numeric import *
from MLab import squeeze
from scipy import r1array
from fastumath import *
import numpyio
import struct, os, sys

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
      fwrite -- write to file from Numeric array
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
        """
        data = asarray(data)
        if mtype is None:
            mtype = data.typecode()
        howmany,mtype = getsize_type(mtype)
        count = product(data.shape)
        numpyio.fwrite(self.fid,count,data,mtype,self.bs)
        return 

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
        if type(count) in [types.TupleType, types.ListType]:
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
                                      
def loadmat(name, dict=None, appendmat=1):
    """Load the MATLAB mat file saved in level 1.0 format.

    If name is a full path name load it in.  Otherwise search for the file
    on the sys.path list and load the first one found (the current directory
    is searched first).

    Only Level 1.0 MAT files are supported so far.

    Inputs:

      name -- name of the mat file (don't need .mat extension)
      dict -- the dictionary to insert into.  If none the variables will be
              returned in a dictionary.
      appendmat -- non-zero to append the .mat extension to the end of the
                   given filename.

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
                fid = open(test_name,'r')
                fid.close()
                full_name = test_name
            except IOError:
                pass
        if full_name is None:
            raise IOError, "%s not found on the path." % name

    fid = fopen(full_name,'r')
    test_vals = fid.fread(4,'byte')
    if not (0 in test_vals):
        fid.close()
        raise ValueError, "Version 5.0 file format not supported."

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

    fid.close()
    fid = fopen(full_name, 'r', format)

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

        M,rest = divmod(header[0],1000)
        O,rest = divmod(rest,100)
        P,rest = divmod(rest,10)
        T = rest

        if (M > 1):
            fid.close()
            raise ValueError, "Unsupported binary format."
        if (O != 0):
            fid.close()
            raise ValuError, "Hundreds digit of first integer should be zero."

        if (P == 4):
            fid.close()
            raise ValueError, "No support for 16-bit unsigned integers."

        if (T not in [0,1]):
            fid.close()
            raise ValueError, "Cannot handle sparse matrices, yet."

        storage = {0:'d',1:'f',2:'i',3:'s',5:'b'}[P]

        varname = fid.fread(header[-1],'char')[:-1]
        varname = varname.tostring()
        defnames.append(varname)
        numels = header[1]*header[2]
        if T == 0:             # Text data
            data = r1array(fid.fread(numels,storage))
            if header[3]:  # imaginary data
                data2 = fid.fread(numels,storage)
                if data.typecode() == 'f' and data2.typecode() == 'f':
                    new = zeros(data.shape,'F')
                    new.real = data
                    new.imag = data2
                    data = new
                    del(new)
                    del(data2)
            if len(data) > 1:
                data.shape = (header[2], header[1])                
                thisdict[varname] = transpose(squeeze(data))
            else:
                thisdict[varname] = data
        else:
            data = r1array(fid.fread(numels,storage,'char'))
            if len(data) > 1:
                data.shape = (header[2], header[1])
                thisdict[varname] = transpose(squeeze(data))
            else:
                thisdict[varname] = data

    fid.close()
    if dict is not None:
        print "Names defined = ", defnames
        dict.update(thisdict)
    else:
        return thisdict


def savemat(filename, dict):
    """Save a dictionary of names and arrays into the MATLAB-style .mat file.

    This saves the arrayobjects in the given dictionary to a Matlab Version 4
    style .mat file.
    """
    storage = {'D':0,'d':0,'F':1,'f':1,'l':2,'i':2,'s':3,'b':5}
    if filename[-4:] != ".mat":
        filename = filename + ".mat"
    fid = fopen(filename,'w')
    M = not LittleEndian
    O = 0
    for variable in dict.keys():
        var = dict[variable]
        if type(var) is not ArrayType:
            continue
        if var.typecode() == 'c':
            T = 1
        else:
            T = 0
        if var.typecode() == '1':
            var = var.astype('s')
        P = storage[var.typecode()]
        fid.fwrite([M*1000+O*100+P*10+T],'int')

        if len(var.shape) == 1:
            var.shape = (len(var), 1)
        var = transpose(var)

        if len(var.shape) > 2:
            var.shape = (product(var.shape[:-1]), var.shape[-1])

        imagf = var.typecode() in ['F', 'D']
        fid.fwrite([var.shape[1], var.shape[0], imagf, len(variable)+1],'int')
        fid.fwrite(variable+'\x00','char')
        if imagf:
            fid.fwrite(var.real)
            fid.fwrite(var.imag)
        else:
            fid.fwrite(var)
    fid.close()
    return
        
        
            




