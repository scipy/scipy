# Authors:  Travis Oliphant, Trent Oliphant
# with support from Lee Barford's group at Agilent, Inc.
#

"""This module allows for the loading of an array from an ASCII
Text File

"""

__all__ = ['read_array', 'write_array']

# Standard library imports.
import os
import re
import sys
import types

# Numpy imports.
import numpy

from numpy import array, take, concatenate, asarray, real, imag, \
  deprecate_with_doc
# Sadly, this module is still written with typecodes in mind.
from numpy.oldnumeric import Float

# Local imports.
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
# Written by Travis Oliphant and Trent Oliphant
#   with support from Agilent, Inc.
#


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
    try:
        # this is the duck typing check: if fileobject
        # can be used is os.path.expanduser, it is a string
        # otherwise it is a fileobject
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
    except (AttributeError, TypeError):
        # it is assumed that the fileobject is a python
        # file object if it can not be used in os.path.expanduser
        file = fileobject

    return file


class ascii_stream(object):
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
        self.should_close_file = not (self.file is fileobject)
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
            if buffer == '':
                return []
            else:
                return lines
        else:
            if len(lines) < 2:
                raise ValueError, "Buffer size too small."
            backup = len(lines[-1])
            self.file.seek(-backup, 1)
            return lines[:-1]

    def __del__(self):
        if hasattr(getattr(self, 'file', None),'close') and self.should_close_file:
            self.file.close()

    def __getitem__(self, item):
        while 1:
            line = self.readnextline()
            if line is None:
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


def extract_columns(arlist, collist, atype, missing):
    if collist[-1] < 0:
        if len(collist) == 1:
            toconvlist = arlist[::-collist[-1]]
        else:
            toconvlist = take(arlist,collist[:-1],0)
            toconvlist = concatenate((toconvlist,
                                      arlist[(collist[-2]-collist[-1])::(-collist[-1])]))
    else:
        toconvlist = take(arlist, collist,0)

    return numpyio.convert_objectarray(toconvlist, atype, missing)


# Given a string representing one line, a separator tuple, a list of
#  columns to read for each element of the atype list and a missing
#  value to insert when conversion fails.

# Regular expressions for detecting complex numbers and for dealing
#  with spaces between the real and imaginary parts

_obj = re.compile(r"""
      ([0-9.eE]+)            # Real part
      ([\t ]*)               # Space between real and imaginary part
      ([+-])                 # +/- sign
      ([\t ]*)               # 0 or more spaces
      (([0-9.eE]+[iIjJ])
      |([iIjJ][0-9.eE]+))    # Imaginary part
      """, re.VERBOSE)

_not_warned = 1
def process_line(line, separator, collist, atype, missing):
    global _not_warned
    strlist = []
    line = _obj.sub(r"\1\3\5",line)  # remove spaces between real
                                     # and imaginary parts of complex numbers

    if _not_warned:
        warn = 0
        if (_obj.search(line) is not None):
            warn = 1
            for k in range(len(atype)):
                if atype[k] in numpy.typecodes['Complex']:
                    warn = 0
        if warn:
            numpy.disp("Warning: Complex data detected, but no requested typecode was complex.")
            _not_warned = 0
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
    N = len(atype)
    vals = [None]*N
    for k in range(len(atype)):
        vals[k] = extract_columns(arlist, collist[k], atype[k], missing)
    return vals

def getcolumns(stream, columns, separator):
    global _not_warned
    comment = stream.comment
    lenc = stream.lencomment
    k, K = stream.linelist[0], len(stream._buffer)
    while k < K:
        firstline = stream._buffer[k]
        if firstline != '' and firstline[:lenc] != comment:
            break
        k = k + 1
    if k == K:
        raise ValueError, "First line to read not within %d lines of top." % K
    firstline = stream._buffer[k]
    N = len(columns)
    collist = [None]*N
    colsize = [None]*N
    for k in range(N):
        collist[k] = build_numberlist(columns[k])
    _not_warned = 0
    val = process_line(firstline, separator, collist, [Float]*N, 0)
    for k in range(N):
        colsize[k] = len(val[k])
    return colsize, collist

def convert_to_equal_lists(cols, atype):
    if not isinstance(cols, types.ListType):
        cols = [cols]
    if not isinstance(atype, types.ListType):
        atype = [atype]
    N = len(cols) - len(atype)
    if N > 0:
        atype.extend([atype[-1]]*N)
    elif N < 0:
        cols.extend([cols[-1]]*(-N))
    return cols, atype


@deprecate_with_doc("""
The functionality of read_array is in numpy.loadtxt which allows the same
functionality using different syntax.
""")
def read_array(fileobject, separator=default, columns=default, comment="#",
               lines=default, atype=Float, linesep='\n',
               rowsize=10000, missing=0):
    """Return an array or arrays from ascii_formatted data in |fileobject|.

    Inputs:

      fileobject -- An open file object or a string for a valid filename.
                    The string can be prepended by "~/" or "~<name>/" to
                    read a file from the home directory.
      separator -- a string or a tuple of strings to indicate the column
                   separators.  If the length of the string tuple is less
                   than the total number of columns, then the last separator
                   is assumed to be the separator for the rest of the columns.
      columns -- a tuple of integers and range-tuples which describe the
                 columns to read from the file.  A negative entry in the
                 last column specifies the negative skip value to the end.
                 Example:  columns=(1, 4, (5, 9), (11, 15, 3), 17, -2)
                         will read [1,4,5,6,7,8,11,14,17,19,21,23,...]
                 If multiple arrays are to be returned, then this argument
                 should be an ordered list of such tuples.  There should be
                 one entry in the list for each arraytype in the atype list.
      lines   -- a tuple with the same structure as columns which indicates
                 the lines to read.
      comment -- the comment character (line will be ignored even if it is
                 specified by the lines tuple)
      linesep -- separator between rows.
      missing -- value to insert in array when conversion to number fails.
      atype -- the typecode of the output array.  If multiple outputs are
               desired, then this should be a list of typecodes.  The columns
               to fill the array represented by the given typecode is
               determined from the columns argument.  If the length of atype
               does not match the length of the columns list, then, the
               smallest one is expanded to match the largest by repeatedly
               copying the last entry.
      rowsize -- the allocation row size (array grows by this amount as
               data is read in).

    Output -- the 1 or 2d array, or a tuple of output arrays of different
              types, sorted in order of the first column to be placed
              in the output array.

    """

    global _not_warned
    # Make separator into a tuple of separators.
    if type(separator) in [types.StringType, type(default)]:
        sep = (separator,)
    else:
        sep = tuple(separator)
    # Create ascii_object from |fileobject| argument.
    ascii_object = ascii_stream(fileobject, lines=lines, comment=comment, linesep=linesep)
    columns, atype = convert_to_equal_lists(columns, atype)
    numout = len(atype)
    # Get the number of columns to read and expand the columns argument
    colsize, collist = getcolumns(ascii_object, columns, sep)
    # Intialize the output arrays
    outrange = range(numout)
    outarr = []
    typecodes = "".join(numpy.typecodes.values())
    for k in outrange:
        if not atype[k] in typecodes:
            raise ValueError, "One of the array types is invalid, k=%d" % k
        outarr.append(numpy.zeros((rowsize, colsize[k]),atype[k]))
    row = 0
    block_row = 0
    _not_warned = 1
    for line in ascii_object:
        if line.strip() == '':
            continue
        vals = process_line(line, sep, collist, atype, missing)
        for k in outrange:
            outarr[k][row] = vals[k]
        row += 1
        block_row += 1
        if block_row >= rowsize:
            for k in outrange:
                outarr[k].resize((outarr[k].shape[0] + rowsize,colsize[k]))
            block_row = 0
    for k in outrange:
        if outarr[k].shape[0] != row:
            outarr[k].resize((row,colsize[k]))
        a = outarr[k]
        if a.shape[0] == 1 or a.shape[1] == 1:
            outarr[k] = numpy.ravel(a)
    if len(outarr) == 1:
        return outarr[0]
    else:
        return tuple(outarr)


# takes 1-d array and returns a string
def str_array(arr, precision=5,col_sep=' ',row_sep="\n",ss=0):
    thestr = []
    arr = asarray(arr)
    N,M = arr.shape
    thistype = arr.dtype.char
    nofloat = (thistype in '1silbwu') or (thistype in 'Oc')
    cmplx = thistype in 'FD'
    fmtstr = "%%.%de" % precision
    cmpnum = pow(10.0,-precision)
    for n in xrange(N):
        theline = []
        for m in xrange(M):
            val = arr[n,m]
            if ss and abs(val) < cmpnum:
                val = 0*val
            if nofloat or val==0:
                thisval = str(val)
            elif cmplx:
                rval = real(val)
                ival = imag(val)
                thisval = eval('fmtstr % rval')
                if (ival >= 0):
                    istr = eval('fmtstr % ival')
                    thisval = '%s+j%s' % (thisval, istr)
                else:
                    istr = eval('fmtstr % abs(ival)')
                    thisval = '%s-j%s' % (thisval, istr)
            else:
                thisval = eval('fmtstr % val')
            theline.append(thisval)
        strline = col_sep.join(theline)
        thestr.append(strline)
    return row_sep.join(thestr)


@deprecate_with_doc("""

This function is replaced by numpy.savetxt which allows the same functionality
through a different syntax.
""")
def write_array(fileobject, arr, separator=" ", linesep='\n',
                precision=5, suppress_small=0, keep_open=0):
    """Write a rank-2 or less array to file represented by fileobject.

    Inputs:

      fileobject -- An open file object or a string to a valid filename.
      arr -- The array to write.
      separator -- separator to write between elements of the array.
      linesep -- separator to write between rows of array
      precision -- number of digits after the decimal place to write.
      suppress_small -- non-zero to round small numbers down to 0.0
      keep_open = non-zero to return the open file, otherwise, the file is closed.
    Outputs:

      file -- The open file (if keep_open is non-zero)
    """
    # XXX: What to when appending to files ? 'wa' does not do what one might
    # expect, and opening a file twice to create it first is not easily doable
    # with get_open_file ?
    file = get_open_file(fileobject, mode='w')
    rank = numpy.rank(arr)
    if rank > 2:
        raise ValueError, "Can-only write up to 2-D arrays."

    if rank == 0:
        h = 1
        arr = numpy.reshape(arr, (1,1))
    elif rank == 1:
        h = numpy.shape(arr)[0]
        arr = numpy.reshape(arr, (h,1))
    else:
        h = numpy.shape(arr)[0]
        arr = asarray(arr)

    for ch in separator:
        if ch in '0123456789-+FfeEgGjJIi.':
            raise ValueError, "Bad string for separator"

    astr = str_array(arr, precision=precision,
                     col_sep=separator, row_sep=linesep,
                     ss = suppress_small)
    file.write(astr)
    file.write('\n')
    if keep_open:
        return file
    else:
        if file is sys.stdout or file is sys.stderr:
            return
        file.close()
    return
