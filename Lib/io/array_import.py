"""This module allows for the loading of an array from an ASCII
Text File

"""

import Numeric
from Numeric import array, take, concatenate, Float
import types, re
import numpyio
default = None
_READ_BUFFER_SIZE = 1024*1024
#_READ_BUFFER_SIZE = 1000
#_READ_BUFFER_SIZE = 160

# ASCII Text object stream with automatic (un)compression and URL access.
#
# Adapted from
#    TextFile class Written by: Konrad Hinsen <hinsen@cnrs-orleans.fr>
#
# Written by Trent Oliphant and Travis Oliphant
#   with support from Agilent, Inc.
# 

import os, sys

def convert_separator(sep):
    newsep = ''
    for k in sep:
        if k in '.^$*+?{[\\|()':
            newsep = newsep + '\\' + k
        else:
            newsep = newsep + k
    return newsep

def build_numberlist(lines):
    if lines is default:
        linelist = [-1]
    else:
        linelist = []
        errstr = "Argument lines must be a sequence of integers and/or range tuples."
        try:
            for num in lines[:-1]:  # handle all but last element 
                if type(num) not in [types.IntType, types.TupleType]:
                    raise ValueError, errstr
                if isinstance(num, types.IntType):
                    linelist.append(num)
                else:
                    if not 1 < len(num) < 4:
                        raise ValueError, "Tuples must be valid range tuples."                        
                    linelist.extend(range(*num))
        except TypeError:
            raise ValueError, errstr
        num = lines[-1]
        if type(num) is types.IntType:
            linelist.append(num)
        elif type(num) is types.TupleType:
            if [types.IntType]*len(num) != map(type, num):
                if len(num) > 1 and num[1] is not None:                        
                    raise ValueError, errstr
            if len(num) == 1:
                linelist.extend([num[0],-1])
            elif len(num) == 2:
                if num[1] is None:
                    linelist.extend([num[0], -1])
                else:
                    linelist.extend(range(*num))
            elif len(num) == 3:
                if num[1] is None:
                    linelist.extend([num[0], -num[2]])
                else:
                    linelist.extend(range(*num))
            else:
                raise ValueError, errstr
    return linelist

def get_open_file(fileobject, mode='rb'):
    if isinstance(fileobject, types.StringType):
        fileobject = os.path.expanduser(fileobject)
        if mode[0]=='r' and not os.path.exists(fileobject):
            raise IOError, (2, 'No such file or directory: '
                            + fileobject)
        else:
            try:
                file = open(fileobject, mode)
            except IOError, details:
                file = None
                if type(details) == type(()):
                    details = details + (fileobject,)
                raise IOError, details
    else:
        file = fileobject

    return file
            

class ascii_stream:

    """Text files with line iteration

    Ascii_stream instances can be used like normal read-only file objects
    (i.e. by calling readline() and readlines()), but can
    also be used as sequences of lines in for-loops.

    Finally, ascii_stream objects accept file names that start with '~' or
    '~user' to indicate a home directory(for reading only).

    Constructor: ascii_stream(|fileobject|, |lines|,|comment|),
    where |fileobject| is either an open python file object or 
    the name of the file, |lines| is a sequence of integers
    or tuples(indicating ranges) of lines to be read, |comment| is the
    comment line identifier """

    def __init__(self, fileobject, lines=default, comment="#",
                 linesep='\n'):
        if not isinstance(comment, types.StringType):
            raise ValueError, "Comment must be a string."
        self.linelist = build_numberlist(lines)
        self.comment = comment
        self.lencomment = len(comment)
        self.file = get_open_file(fileobject, mode='r')
        self._pos = self.file.tell()
        self._lineindex = 0
        if self.linelist[-1] < 0:
            self._linetoget = self.linelist[-1]
        else:
            self._linetoget = 0
        self._oldbuflines = 0
        self._linesplitter = linesep
        self._buffer = self.readlines(_READ_BUFFER_SIZE)
        self._totbuflines = len(self._buffer)

    def readlines(self, sizehint):
        buffer = self.file.read(sizehint)
        lines = buffer.split(self._linesplitter)
        if len(buffer) < sizehint:  # EOF
            return lines
        else:
            if len(lines) < 2:
                raise ValueError, "Buffer size too small."
            backup = len(lines[-1])
            self.file.seek(-backup, 1)
            return lines[:-1]

    def __del__(self):
        if hasattr(self.file,'close'):
            self.file.close()

    def __getitem__(self, item):
        while 1:
            line = self.readnextline()
            if not line:
                raise IndexError
            if len(line) < self.lencomment or line[:self.lencomment] != self.comment:
                break
        return line
        
    def readnextline(self):
        if self.linelist[self._lineindex] >= 0:
            self._linetoget = self.linelist[self._lineindex]
            self._lineindex += 1
        else:
            self._linetoget = self._linetoget - self.linelist[self._lineindex]
        while self._linetoget >= self._totbuflines:
            self._buffer = self.readlines(_READ_BUFFER_SIZE)
            self._oldbuflines = self._totbuflines
            self._totbuflines += len(self._buffer)
            if (self._totbuflines == self._oldbuflines):
                return None
        line = self._buffer[self._linetoget - self._oldbuflines]
        return line    

    def close(self):
        self.file.close()

    def flush(self):
        self.file.flush()


def move_past_spaces(firstline):
    ind = 0
    firstline = firstline.lstrip()
    while firstline[ind] not in [' ','\n','\t','\v','\f','\r']:
        ind += 1
    return firstline[ind:], ind

def process_line(line, separator, collist, atype, missing):
    strlist = []
    for mysep in separator[:-1]:
        if mysep is None:
            newline, ind = move_past_spaces(line)
            strlist.append(line[:ind])
            line = newline
        else:
            ind = line.find(mysep)
            strlist.append(line[:ind])
            line = line[ind+len(mysep):]
    strlist.extend(line.split(separator[-1]))
    arlist = array(strlist,'O')
    if collist[-1] < 0:
        if len(collist) == 1:
            toconvlist = arlist[::-collist[-1]]
        else:
            toconvlist = take(arlist,collist[:-1])
            toconvlist = concatenate((toconvlist,
                                      arlist[(collist[-2]+1)::(-collist[-1])]))
    else:
        toconvlist = take(arlist, collist)

    return numpyio.convert_objectarray(toconvlist, atype, missing)

def getcolumns(stream, columns, separator):
    firstline = stream._buffer[0]
    collist = build_numberlist(columns)
    colsize = len(collist)
    val = process_line(firstline, separator, collist, Float, 0)
    return len(val), collist

def read_array(fileobject, separator=default, columns=default, comment="#",
               lines=default, atype=Numeric.Float, linesep='\n',
               rowsize=10000, missing=0):
    """Return an array represented ascii_formatted data in |fileobject|.

    Inputs:

      fileobject -- An open file object or a string for a valid filename.
      separator -- a string or a tuple of strings to indicate the column
                   separators.  If the length of the string tuple is less
                   than the total number of columns, then the last separator
                   is assumed to be the separator for the rest of the columns.
      columns -- a tuple of integers and range-tuples which describe the
                 columns to read from the file.  A negative entry in the
                 last column specifies the negative skip value to the end.
                 Example:  columns=(1, 4, (5, 9), (11, 15, 3), 17, -2)
                         will read [1,4,5,6,7,8,11,14,17,18,20,22,...]
      lines   -- a tuple with the same structure as columns which indicates
                 the lines to read. 
      comment -- the comment character (line will be ignored even if it is
                 specified by the lines tuple)
      linesep -- separator between rows.
      missing -- value to insert in array when conversion to number fails.
      atype -- the typecode of the output array.
      rowsize -- the allocation row size (array grows by this amount as
               data is read in).

    Outputs:

      file -- The open file.
    """

    if type(separator) in [types.StringType, type(default)]:
        sep = (separator,)
    else:
        sep = tuple(separator)
    ascii_object = ascii_stream(fileobject, lines=lines, comment=comment, linesep=linesep)
    colsize, collist = getcolumns(ascii_object, columns, sep)
    a = Numeric.zeros((rowsize, colsize),atype)
    row = 0
    block_row = 0
    for line in ascii_object:        
        a[row] = process_line(line, sep, collist, atype, missing)
        row += 1
        block_row += 1
        if block_row >= rowsize:
            print a.shape[0], rowsize
            a.resize((a.shape[0] + rowsize,colsize))
            block_row = 0
    if a.shape[0] != row:
        a.resize((row,colsize))
    if a.shape[0] == 1 or a.shape[1] == 1:
        a = Numeric.ravel(a)
    return a

def write_array(fileobject, arr, separator=" ", linesep='\n',
                precision=5, suppress_small=0):
    """Write a rank-2 or less array to file represented by fileobject.

    Inputs:

      fileobject -- An open file object or a string to a valid filename.
      arr -- The array to write.
      separator -- separator to write between elements of the array.
      linesep -- separator to write between rows of array
      precision -- number of digits after the decimal place to write.
      suppress_small -- non-zero to suppress small digits and not write
                        them in scientific notation.

    Outputs:

      file -- The open file.
    """
      
    file = get_open_file(fileobject, mode='wa')
    rank = len(arr.shape)
    if rank > 2:
        raise ValueError, "Can-only write up to 2-D arrays."

    if rank == 0:
        file.write(str(a[0])+linesep)
    for k in range(arr.shape[0]):
        astr = Numeric.array2string(arr[k], max_line_width=sys.maxint,
                                    precision=precision,
                                    suppress_small=suppress_small,
                                    separator=' '+separator,
                                    array_output=0)
        astr = astr[1:-1]
        file.write(astr)
        file.write(linesep)
    return file
                 
