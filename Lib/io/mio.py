# Author: Travis Oliphant

from MLab import squeeze
from scipy_base import *
from scipy_base.fastumath import *
import numpyio
import struct, os, sys
import types
try:
    import scipy.sparse
    have_sparse = 1
except ImportError:
    have_sparse = 0

if sys.version_info[0] < 2 or sys.version_info[1] < 3:
    False = 0
    True = 1

__all__ = ['fopen','loadmat','savemat']

def getsize_type(mtype):
    if mtype in ['b','uchar','byte','unsigned char','integer*1', 'int8']:
        mtype = 'b'
    elif mtype in ['c', 'char','char*1']:
        mtype = 'c'
    elif mtype in ['1','schar', 'signed char']:
        mtype = '1'
    elif mtype in ['s','short','int16','integer*2']:
        mtype = 's'
    elif mtype in ['w','ushort','uint16','unsigned short']:
        mtype = 'w'
    elif mtype in ['i','int']:
        mtype = 'i'
    elif mtype in ['u','uint','uint32','unsigned int']:
        mtype = 'u'
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

      read -- read data from file and return Numeric array
      write -- write to file from Numeric array
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
    
    def __init__(self,file_name,permission='rb',format='n'):
        if 'b' not in permission: permission += 'b'
        if type(file_name) in (types.StringType, types.UnicodeType):
            self.__dict__['fid'] = open(file_name,permission)
        elif 'fileno' in file_name.__methods__:  # first argument is an open file
            self.__dict__['fid'] = file_name 
        if format in ['native','n','default']:
            self.__dict__['bs'] = 0
            self.__dict__['format'] = 'native'
        elif format in ['ieee-le','l','little-endian','le']:
            self.__dict__['bs'] = not LittleEndian
            self.__dict__['format'] = 'ieee-le'
        elif format in ['ieee-be','b','big-endian','be']:
            self.__dict__['bs'] = LittleEndian
            self.__dict__['format'] = 'ieee-be'
        else:
            raise ValueError, "Unrecognized format: " + format

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

    def setformat(self, format):
        if format in ['native','n','default']:
            self.__dict__['bs'] = False
            self.__dict__['format'] = 'native'
        elif format in ['ieee-le','l','little-endian','le']:
            self.__dict__['bs'] = not LittleEndian
            self.__dict__['format'] = 'ieee-le'
        elif format in ['ieee-be','b','big-endian','be']:
            self.__dict__['bs'] = LittleEndian
            self.__dict__['format'] = 'ieee-be'
        else:
            raise ValueError, "Unrecognized format: " + format
        return
        
    def write(self,data,mtype=None,bs=None):
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
                   unsigned short : 'w', 'ushort','uint16','unsigned short'
                   int            : 'i', 'int'
                   unsigned int   : 'u', 'uint32','uint','unsigned int'
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
        data = asarray(data)
        if mtype is None:
            mtype = data.typecode()
        howmany,mtype = getsize_type(mtype)
        count = product(data.shape)
        numpyio.fwrite(self.fid,count,data,mtype,bs)
        return 

    fwrite = write

    def read(self,count,stype,rtype=None,bs=None,c_is_b=0):
        """Read data from file and return it in a Numeric array.

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

          output -- a Numeric array of type rtype.
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
                now = self.fid.tell()
                self.fid.seek(0,2)
                end = self.fid.tell()
                self.fid.seek(now)
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
        retval = numpyio.fread(self.fid, count, stype, rtype, bs)
        if len(retval) == 1:
            retval = retval[0]
        if shape is not None:
            retval = resize(retval, shape)
        return retval

    fread = read

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
        if type(fmt) in (types.StringType, types.UnicodeType):
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
            if len(args) > 0:
                sz,mtype = getsize_type(args[0])
            else:
                sz,mtype = getsize_type(fmt.typecode())
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
            if (self.fid.read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            strdata = self.fid.read(numbytes)
            if strdata == '':
                raise ValueError, "Unexpected end of file..."
            data = struct.unpack(fmt,strdata)
            if (self.fid.read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            return data
        else:  # Ignore format string and read in next record as an array.
            fmt = lookup_dict[self.format] + "i"
            nn = struct.calcsize(fmt)
            nbytestr = self.fid.read(nn)
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
            retval = numpyio.fread(self.fid, ncount, dtype, dtype, self.bs)
            if len(retval) == 1:
                retval = retval[0]
            if (self.fid.read(nn) == ''):
                raise ValueError, "Unexpected end of file..."
            return retval


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
    miINT8 : ('miINT8', 1,'1'),
    miUINT8 : ('miUINT8', 1,'b'),
    miINT16 : ('miINT16', 2,'s'),
    miUINT16 :('miUINT16',2,'w'),
    miINT32 : ('miINT32',4,'l'),
    miUINT32 : ('miUINT32',4,'u'),
    miSINGLE : ('miSINGLE',4,'f'),
    miDOUBLE : ('miDOUBLE',8,'d'),
    miINT64 : ('miINT64',8,'c'),   # Handle 64 bit integers as string data.
    miUINT64 : ('miUINT64',8,'c'),
    miMATRIX : ('miMATRIX',0,None),
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

def _parse_header(fid, dict):
    correct_endian = (ord('M')<<8) + ord('I')
                 # if this number is read no BS    
    fid.seek(126)  # skip to endian detector
    endian_test = fid.read(1,'int16')
    if (endian_test == correct_endian): openstr = 'n'
    else:  # must byteswap
        if LittleEndian:
            openstr = 'b'
        else: openstr = 'l'
    fid.setformat(openstr)  # change byte-order if necessary
    fid.rewind()    
    dict['__header__'] = fid.fid.read(124).strip(' \t\n\000')
    vers = fid.read(1,'int16')
    dict['__version__'] = '%d.%d' % (vers >> 8, vers & 255)
    fid.seek(2,1)  # move to start of data
    return

def _parse_array_flags(fid):
    # first 8 bytes are always miUINT32 and 8 --- just a check
    dtype, nbytes = fid.read(2,'u')
    if (dtype != miUINT32) or (nbytes != 8):
        raise IOError, "Invalid MAT file. Perhaps a byte-order problem."

    # read array flags.
    rawflags = fid.read(2,'u')
    class_ = rawflags[0] & 255
    flags = (rawflags[0] & 65535) >> 8
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
    dims = _get_element(fid)[0]
    name = ''.join(asarray(_get_element(fid)[0]).astype('c'))
    if dclass in mxArrays:
        result, unused =_get_element(fid)
        if type == mxCHAR_CLASS:
            result = ''.join(asarray(result).astype('c'))
        else:
            if cmplx:
                imag, unused =_get_element(fid)
                result = result + cast[imag.typecode()](1j) * imag
            result = squeeze(transpose(reshape(result,dims[::-1])))

    elif dclass == mxCELL_CLASS:
        length = product(dims)
        result = zeros(length, PyObject)
        for i in range(length):
            sa, unused = _get_element(fid)
            result[i]= sa
        result = squeeze(transpose(reshape(result,dims[::-1])))
        if rank(result)==0: result = result.toscalar()
        
    elif dclass == mxSTRUCT_CLASS:
        length = product(dims)
        result = zeros(length, PyObject)
        namelength = _get_element(fid)[0]
        # get field names
        names = _get_element(fid)[0]
        splitnames = [names[i:i+namelength] for i in \
                      xrange(0,len(names),namelength)]
        fieldnames = [''.join(asarray(x).astype('c')).strip('\x00')
                              for x in splitnames]
        for i in range(length):
            result[i] = mat_struct()
            for element in fieldnames:
                val,unused = _get_element(fid)
                result[i].__dict__[element] = val
        result = squeeze(transpose(reshape(result,dims[::-1])))
        if rank(result)==0: result = result.toscalar()        

        # object is like a structure with but with a class name
    elif dclass == mxOBJECT_CLASS:
        class_name = ''.join(asarray(_get_element(fid)[0]).astype('c'))
        length = product(dims)
        result = zeros(length, PyObject)
        namelength = _get_element(fid)[0]
        # get field names
        names = _get_element(fid)[0]
        splitnames = [names[i:i+namelength] for i in \
                      xrange(0,len(names),namelength)]
        fieldnames = [''.join(asarray(x).astype('c')).strip('\x00')
                              for x in splitnames]
        for i in range(length):
            result[i] = mat_obj()
            result[i]._classname = class_name
            for element in fieldnames:
                val,unused = _get_element(fid)
                result[i].__dict__[element] = val
        result = squeeze(transpose(reshape(result,dims[::-1])))
        if rank(result)==0: result = result.toscalar()        
         
    elif dclass == mxSPARSE_CLASS:
        rowind, unused = _get_element(fid)
        colind, unused = _get_element(fid)
        res, unused = _get_element(fid)
        if cmplx:
            imag, unused = _get_element(fid)
            res = res + cast[imag.typecode()](1j)*imag
        if have_sparse:
            spmat = scipy.sparse.spmatrix(dims[1],dims[0],typecode=res.typecode())
            spmat.data = res
            # stored with zero-based indices --- spmatrix needs 1 based
            #  indices internally
            spmat.index = [rowind+1, colind+1]
            spmat.nzmax = len(res)
            spmat.lastel = len(res)-1
            result = spmat.transp(inplace=1)
            result = spmat
        else:
            result = (rowind, colind, res)

    return result, name
    
# Return a Python object for the element
def _get_element(fid):

    test = fid.fid.read(1)
    if len(test) == 0:  # nothing left
        raise EOFError
    else:
        fid.rewind(1)
    # get the data tag
    raw_tag = fid.read(1,'u')
    
    # check for compressed
    numbytes = raw_tag >> 16
    if numbytes > 0:  # compressed format
        if numbytes > 4:
            raise IOError, "Problem with MAT file: " \
                  "too many bytes in compressed format."
        dtype = raw_tag & 65535
        el = fid.read(numbytes,miDataTypes[dtype][2],c_is_b=1)
        fid.seek(4-numbytes,1)  # skip padding
        return el, None

    # otherwise parse tag
    dtype = raw_tag
    numbytes = fid.read(1,'u')
    if dtype != miMATRIX:  # basic data type 
        try:
            outarr = fid.read(numbytes,miDataTypes[dtype][2],c_is_b=1)
        except KeyError:
            raise ValueError, "Unknown data type"
        mod8 = numbytes%8
        if mod8:       # skip past padding
            skip = 8-mod8
            fid.seek(skip,1)
        return outarr, None

    # handle miMatrix type
    el, name = _parse_mimatrix(fid,numbytes)
    return el, name

def _loadv5(fid,basename):
    # return a dictionary from a Matlab version 5 file
    # always contains the variable __header__
    dict = {}
    _parse_header(fid,dict)
    var = 0
    while 1:  # file pointer to start of next data
        try:
            var = var + 1
            el, varname = _get_element(fid)
            if varname is None:
                varname = '%s_%04d' % (basename,var)
            dict[varname] = el
        except EOFError:
            break
    return dict

### END MATLAB v5 support #############

def loadmat(name, dict=None, appendmat=1, basename='raw'):
    """Load the MATLAB(tm) mat file.

    If name is a full path name load it in.  Otherwise search for the file
    on the sys.path list and load the first one found (the current directory
    is searched first).

    Both v4 (Level 1.0) and v5 matfiles are supported.

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
            except IOError:
                pass
        if full_name is None:
            raise IOError, "%s not found on the path." % name

    fid = fopen(full_name,'rb')
    test_vals = fid.fread(4,'byte')

    if not (0 in test_vals):       # MATLAB version 5 format
        fid.rewind()
        thisdict = _loadv5(fid,basename)
        if dict is not None:
            dict.update(thisdict)
            return
        else:
            return thisdict
        

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

        if (T not in [0,1]):
            fid.close()
            raise ValueError, "Cannot handle sparse matrices, yet."

        storage = {0:'d',1:'f',2:'i',3:'s',4:'w',5:'b'}[P]

        varname = fid.fread(header[-1],'char')[:-1]
        varname = varname.tostring()
        defnames.append(varname)
        numels = header[1]*header[2]
        if T == 0:             # Text data
            data = atleast_1d(fid.fread(numels,storage))
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
            data = atleast_1d(fid.fread(numels,storage,'char'))
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

    This saves the arrayobjects in the given dictionary to a matlab Version 4
    style .mat file.
    """
    storage = {'D':0,'d':0,'F':1,'f':1,'l':2,'i':2,'s':3,'b':5}
    if filename[-4:] != ".mat":
        filename = filename + ".mat"
    fid = fopen(filename,'wb')
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
        
        
            




