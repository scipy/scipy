"""pstat.py module

############################################################
#######  Written by:  Gary Strangman             ###########
#######  Last modified:  July 9, 2001            ###########
#######  (deleted unnecessary regex dependence)  ###########
############################################################

### February  2002,  all list-only manipulation routines are eliminated,
###                   this version requires Numeric.

This module provides some useful list and array manipulation routines
modeled after those found in the |Stat package by Gary Perlman, plus a
number of other useful list/file manipulation functions.  The list-based
functions include:

      abut (source,*args):  
      simpleabut (source, addon):
  **  colex (arr,cnums):
  **  collapse (arr,keepcols,collapsecols,sterr=None,ns=None):
      dm (arr,criterion):
      get (namepattern,verbose=1):
      getstrings (namepattern):
      linexand (arr,columnlist,valuelist):
      linexor (arr,columnlist,valuelist):
      linedelimited (inlist,delimiter):
      lineincols (inlist,colsize):
      lineincustcols (inlist,colsizes):
      list2string (inlist):
      printcc (lst,extra=2):
      printincols (arr,colsize):
      pl (arr):
      printl(arr):
      put (outlist,fname,writetype='w'):
      replace (lst,oldval,newval):
  **  recode (inlist,listmap,cols='all'):
      remap (arr,criterion):
      rowcompare (row1, row2)
      rowsame (row1, row2)
      sortby(arr,sortcols):
      unique (inlist):
      writecc (arr,file,extra=2,writetype='w'):
      writedelimited (arr, delimiter, file, writetype='w'):  

Currently, the code is all but completely un-optimized.  In many cases, the
array versions of functions amount simply to aliases to built-in array
functions/methods.  Their inclusion here is for function name consistency.
Additions, suggestions, or comments are welcome (strang@nmr.mgh.harvard.edu).

(Note: this file has been modified, and therefore comments should be sent
to scipy-users@scipy.org or scipy-dev@scipy.org).
"""

from Numeric import asarray
import scipy_base.fastumath
from scipy_base import unique
import string, sys, os, copy, math, stats
from types import *

def simpleabut (source, addon):
    """\nConcatenates two lists as columns and returns the result.  2D lists
    are also accomodated for either argument (source or addon).  This DOES NOT
    repeat either list to make the 2 lists of equal length.  Beware of list pairs
    with different lengths ... the resulting list will be the length of the
    SHORTEST list passed.
    
    Returns: a list of lists as long as the SHORTER of source and addon, with
         source on the 'left' and addon on the 'right'.\n
         """

    if type(source) not in [ListType,TupleType]:
        source = [source]
    if type(addon) not in [ListType,TupleType]:
        addon = [addon]
    minlen = min(len(source),len(addon))
    list = copy.deepcopy(source)                # start abut process
    if type(source[0]) not in [ListType,TupleType]:
        if type(addon[0]) not in [ListType,TupleType]:
            for i in range(minlen):
                list[i] = [source[i]] + [addon[i]]      # source/addon = column
        else:
            for i in range(minlen):
                list[i] = [source[i]] + addon[i]        # addon=list-of-lists
    else:
        if type(addon[0]) not in [ListType,TupleType]:
            for i in range(minlen):
                list[i] = source[i] + [addon[i]]        # source=list-of-lists
        else:
            for i in range(minlen):
                list[i] = source[i] + addon[i]  # source/addon = list-of-lists
    source = list
    return source

def getstrings (namepattern,verbose=1):
    """\nLoads a set of text files, with all elements left as string type.
Uses UNIX-style wildcards (i.e., function uses glob).

Format:  getstrings (namepattern, verbose=1)
Returns: a list of strings, one per line in each text file specified by
         namepattern\n"""

    import glob

    fnames = glob.glob(namepattern)
    if len(fnames) == 0:
        if verbose:
            print 'NO FILENAMES MATCH PATTERN !!'
        return None
    print fnames
    elements = []
    for filename in fnames:
        file = open(filename)
        newelements = map(string.split,file.readlines())
        elements = elements + newelements
    return elements

def linedelimited (inlist,delimiter):
    """\nReturns a string composed of elements in inlist, with each element
separated by 'delimiter.'  Used by function writedelimited.  Use '\t' for
tab-delimiting.

Format:  linedelimited (inlist,delimiter)\n"""

    outstr = ''
    for item in inlist:
        if type(item) <> StringType:
            item = str(item)
        outstr = outstr + item + delimiter
    outstr = outstr[0:-1]
    return outstr


def lineincols (inlist,colsize):
    """\nReturns a string composed of elements in inlist, with each element
right-aligned in columns of (fixed) colsize.

Format:  lineincols (inlist,colsize)   where colsize is an integer\n"""

    outstr = ''
    for item in inlist:
        if type(item) <> StringType:
            item = str(item)
        size = len(item)
        if size <= colsize:
            for i in range(colsize-size):
                outstr = outstr + ' '
            outstr = outstr + item
        else:
            outstr = outstr + item[0:colsize+1]
    return outstr


def lineincustcols (inlist,colsizes):
    """\nReturns a string composed of elements in inlist, with each element
right-aligned in a column of width specified by a sequence colsizes.  The
length of colsizes must be greater than or equal to the number of columns in
inlist.

Format:  lineincustcols (inlist,colsizes)
Returns: formatted string created from inlist\n"""

    outstr = ''
    for i in range(len(inlist)):
        if type(inlist[i]) <> StringType:
            item = str(inlist[i])
        else:
            item = inlist[i]
        size = len(item)
        if size <= colsizes[i]:
            for j in range(colsizes[i]-size):
                outstr = outstr + ' '
            outstr = outstr + item
        else:
            outstr = outstr + item[0:colsizes[i]+1]
    return outstr


def list2string (inlist):
    """\nConverts a 1D list to a single long string for file output, using
the string.join function.

Format:  list2string (inlist)
Returns: the string created from inlist\n"""

    def makestr (item):
        if type(item) <> StringType:
            item = str(item)
        return item
    stringlist = map(makestr,inlist)
    return string.join(stringlist)


def printcc (lst,extra=2):
    """\nPrints a list of lists in columns, customized by the max size of items
within the columns (max size of items in col, plus 'extra' number of spaces).
Use 'dashes' or '\n' in the list(oflists) to print dashes or blank lines,
respectively.

Format:  printcc (lst,extra=2)
Returns: None\n"""

    def makestr (x):
        if type(x) <> StringType:
            x = str(x)
        return x
 
    if type(lst[0]) not in [ListType,TupleType]:
        lst = [lst]
    rowstokill = []
    list2print = copy.deepcopy(lst)
    for i in range(len(lst)):
        if lst[i] == ['\n'] or lst[i]=='\n' or lst[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()   # delete blank rows from the end
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = colex(list2print,col)
        items = map(makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in lst:
        if row == ['\n'] or row == '\n':
            print
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            print lineincustcols(dashes,maxsize)
        else:
            print lineincustcols(row,maxsize)
    return None


def printincols (listoflists,colsize):
    """\nPrints a list of lists in columns of (fixed) colsize width, where
colsize is an integer.

Format:  printincols (listoflists,colsize)
Returns: None\n"""

    for row in listoflists:
        print lineincols(row,colsize)
    return None


def pl (listoflists):
    """\nPrints a list of lists, 1 list (row) at a time.

Format:  pl(listoflists)
Returns: None\n"""

    for row in listoflists:
        print row
    return None


def printl(listoflists):
    """Alias for pl."""

    pl(listoflists)
    return


def replace (lst,oldval,newval):
    """\nReplaces all occurrences of 'oldval' with 'newval'.

Format:  replace (lst,oldval,newval)"""

    for i in range(len(lst)):
        if type(lst[i]) not in [ListType,TupleType]:
            if lst[i]==oldval: lst[i]==newval
        else:
            for j in range(len(lst[i])):
                if lst[i][j]==oldval:  lst[i][j]=newval
    return lst


def remap (listoflists,criterion):

    function = 'lines = map(lambda x: '+criterion+',listoflists)'
    exec function in globals(), locals()
    return lines


def sortby(listoflists,sortcols):
    """/nSorts a list of lists on the column(s) specified in the sequence
sortcols.

Format:  sortby(listoflists,sortcols)
Returns: sorted list, unchanged column ordering\n"""

    newlist = abut(colex(listoflists,sortcols),listoflists)
    newlist.sort()
    try:
        numcols = len(sortcols)
    except TypeError:
        numcols = 1
    crit = '[' + str(numcols) + ':]'
    newlist = colex(newlist,crit)
    return newlist

def writecc (listoflists,file,writetype='w',extra=2):
    """\nWrites a list of lists to a file in columns, customized by the max
size of items within the columns (max size of items in col, +2 characters)
to specified file.  File-overwrite is the default.

Format:  writecc (listoflists,file,writetype='w',extra=2)
Returns: None\n"""

    def makestr (x):
        if type(x) <> StringType:
            x = str(x)
        return x

    if type(listoflists[0]) not in [ListType,TupleType]:
        listoflists = [listoflists]
    outfile = open(file,writetype)
    rowstokill = []
    list2print = copy.deepcopy(listoflists)
    for i in range(len(listoflists)):
        if listoflists[i] == ['\n'] or listoflists[i]=='\n' or listoflists[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = colex(list2print,col)
        items = map(makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in listoflists:
        if row == ['\n'] or row == '\n':
            outfile.write('\n')
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            outfile.write(lineincustcols(dashes,maxsize))
        else:
            outfile.write(lineincustcols(row,maxsize))
        outfile.write('\n')
    outfile.close()
    return None


def writedelimited (listoflists, delimiter, file, writetype='w'):
    """\nWrites a list of lists in columns, separated by character(s) delimiter
to specified file.  File-overwrite is the default.

Format:  writedelimited (listoflists,delimiter,filename,writetype='w')
Returns: None\n"""

    def makestr (x):
        if type(x) <> StringType:
            x = str(x)
        return x

    if type(listoflists[0]) not in [ListType,TupleType]:
        listoflists = [listoflists]
    outfile = open(file,writetype)
    rowstokill = []
    list2print = copy.deepcopy(listoflists)
    for i in range(len(listoflists)):
        if listoflists[i] == ['\n'] or listoflists[i]=='\n' or listoflists[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for row in listoflists:
        if row == ['\n'] or row == '\n':
            outfile.write('\n')
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '------'
            outfile.write(linedelimited(dashes,delimiter))
        else:
            outfile.write(linedelimited(row,delimiter))
        outfile.write('\n')
    outfile.close()
    return None

#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================
#===================   PipeSTAT-like ARRAY FUNCTIONS  =====================

import Numeric
N = Numeric

def abut (source, *args):
    """\nLike the |Stat abut command.  It concatenates two arrays column-wise
    and returns the result.  CAUTION:  If one array is shorter, it will be
    repeated until it is as long as the other.

    Format:  aabut (source, args)    where args=any # of arrays
    Returns: an array as long as the LONGEST array past, source appearing on the
    'left', arrays in <args> attached on the 'right'.\n"""

    source = asarray(source)
    if len(source.shape)==1:
        width = 1
        source = N.resize(source,[source.shape[0],width])
    else:
        width = source.shape[1]
    for addon in args:
        if len(addon.shape)==1:
            width = 1
            addon = N.resize(addon,[source.shape[0],width])
        else:
            width = source.shape[1]
        if len(addon) < len(source):
            addon = N.resize(addon,[source.shape[0],addon.shape[1]])
        elif len(source) < len(addon):
            source = N.resize(source,[addon.shape[0],source.shape[1]])
        source = N.concatenate((source,addon),1)
    return source


def colex (a,indices,axis=1):
    """\nExtracts specified indices (a list) from passed array, along passed
    axis (column extraction is default).  BEWARE: A 1D array is presumed to be a
    column-array (and that the whole array will be returned as a column).
 
    Returns: the columns of a specified by indices\n"""

    if type(indices) not in [ListType,TupleType,N.ArrayType]:
        indices = [indices]
    if len(N.shape(a)) == 1:
        cols = N.resize(a,[a.shape[0],1])
    else:
        cols = N.take(a,indices,axis)
    return cols


def collapse (a,keepcols,collapsecols,stderr=0,ns=0,cfcn=None):
    """Averages data in collapsecol, keeping all unique items in keepcols
    (using unique, which keeps unique LISTS of column numbers), retaining
    the unique sets of values in keepcols, the mean for each.  If the sterr or
    N of the mean are desired, set either or both parameters to 1.

    Returns: unique 'conditions' specified by the contents of columns specified
    by keepcols, abutted with the mean(s) of column(s) specified by
    collapsecols
    """
    if cfcn is None:
        cfcn = stats.mean
    a = asarray(a)
    if keepcols == []:
        avgcol = colex(a,collapsecols)
        means = cfcn(avgcol)
        return means
    else:
        if type(keepcols) not in [ListType,TupleType,N.ArrayType]:
            keepcols = [keepcols]
        values = colex(a,keepcols)   # so that "item" can be appended (below)
        uniques = unique(values)  # get a LIST, so .sort keeps rows intact
        uniques.sort()
        newlist = []
        for item in uniques:
            if type(item) not in [ListType,TupleType,N.ArrayType]:
                item =[item]
            tmprows = linexand(a,keepcols,item)
            for col in collapsecols:
                avgcol = colex(tmprows,col)
                item.append(cfcn(avgcol))
                if sterr:
                    if len(avgcol)>1:
                        item.append(stats.sterr(avgcol))
                    else:
                        item.append('N/A')
                if ns:
                    item.append(len(avgcol))
                newlist.append(item)
        try:
            new_a = N.array(newlist)
        except TypeError:
            new_a = N.array(newlist,'O')
        return new_a


def dm (a,criterion):
    """\nReturns rows from the passed list of lists that meet the criteria in
    the passed criterion expression (a string).

    Format:  dm (a,criterion)   where criterion is like 'x[2]==37'\n"""

    function = 'lines = filter(lambda x: '+criterion+',a)'
    exec function in globals(), locals()
    
    try:
        lines = N.array(lines)
    except:
        lines = N.array(lines,'O')
    return lines


def isstring(x):
    return isinstance(x,StringType)

import glob
def get (namepattern):
    """Loads an array from 2D text files (specified by a UNIX-style
    wildcard filename pattern).  ONLY 'GET' FILES WITH EQUAL NUMBERS OF COLUMNS
    ON EVERY ROW (otherwise returned array will be zero-dimensional).

    Returns: an array of integers, floats or objects (type='O'), depending on the
    contents of the files specified by namepattern
    """

    fnames = glob.glob(namepattern)
    if len(fnames) == 0:
        print 'NO FILENAMES MATCH PATTERN !!'
        return None
    print fnames
    elements = []
    for filename in fnames:
        file = open(filename)
        newelements = map(string.split,file.readlines())
        for i in range(len(newelements)):
            for j in range(len(newelements[i])):
                try:
                    newelements[i][j] = string.atof(newelements[i][j])
                except:
                    pass
        elements = elements + newelements
    for row in elements:
        if N.add.reduce(N.array(map(isstring,row)))==len(row):
            print "A row of strings was found.  Returning a LIST."
            return elements
    try:
       elements = N.array(elements)
    except TypeError:
       elements = N.array(elements,'O')
    return elements

def linexand (a,columnlist,valuelist):
    """Returns the rows of an array where col (from columnlist) = val
    (from valuelist).  One value is required for each column in columnlist.

    Returns: the rows of a where columnlist[i]=valuelist[i] for ALL i\n"""

    a = asarray(a)
    if type(columnlist) not in [ListType,TupleType,N.ArrayType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType,N.ArrayType]:
        valuelist = [valuelist]
    criterion = ''
    for i in range(len(columnlist)):
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' and'
    criterion = criterion[0:-3]         # remove the "and" after the last crit
    return adm(a,criterion)


def linexor (a,columnlist,valuelist):
    """Returns the rows of an array where col (from columnlist) = val
    (from valuelist).  One value is required for each column in columnlist.
    The exception is if either columnlist or valuelist has only 1 value, in
    which case that item will be expanded to match the length of the other list.
        
    Returns: the rows of a where columnlist[i]=valuelist[i] for ANY i
    """
    a = asarray(a)
    if type(columnlist) not in [ListType,TupleType,N.ArrayType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType,N.ArrayType]:
        valuelist = [valuelist]
    criterion = ''
    if len(columnlist) == 1 and len(valuelist) > 1:
        columnlist = columnlist*len(valuelist)
    elif len(valuelist) == 1 and len(columnlist) > 1:
        valuelist = valuelist*len(columnlist)
    for i in range(len(columnlist)):
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' or'
    criterion = criterion[0:-2]         # remove the "or" after the last crit
    return adm(a,criterion)

def put (outarray,fname,writetype='w'):
    """Sends passed 1D or 2D array to an output file and closes the file.
    
    Returns: None\n"""
    
    outfile = open(fname,writetype)
    try:
        tmp = len(outarray[0])
    except:
        outarray = outarray[Newaxis,:]
    for row in outarray:
        outfile.write(N.array2string(row)[1:-1])
        outfile.write('\n')
    outfile.close()

def recode (a,listmap,col='all'):
    """Remaps the values in an array to a new set of values (useful when
    you need to recode data from (e.g.) strings to numbers as most stats packages
    require.  Can work on SINGLE columns, or 'all' columns at once.
    
    Format:  arecode (a,listmap,col='all')
    Returns: a version of array a where listmap[i][0] = (instead) listmap[i][1]
    """
    a = asarray(a)
    ashape = a.shape
    if col == 'all':
        work = a.flat
    else:
        work = acolex(a,col)
        work = work.flat
    for pair in listmap:
        if type(pair[1]) == StringType or work.typecode()=='O' or a.typecode()=='O':
            work = N.array(work,'O')
            a = N.array(a,'O')
            for i in range(len(work)):
                if work[i]==pair[0]:
                    work[i] = pair[1]
            if col == 'all':
                return N.reshape(work,ashape)
            else:
                return N.concatenate([a[:,0:col],work[:,N.NewAxis],a[:,col+1:]],1)
        else:   # must be a non-Object type array and replacement
            work = N.where(N.equal(work,pair[0]),pair[1],work)
            return N.concatenate([a[:,0:col],work[:,N.NewAxis],a[:,col+1:]],1)

def rowcompare(row1, row2):
    """Compares two rows from an array, regardless of whether it is an
    array of numbers or of python objects (which requires the cmp function).
    
    Format:  rowcompare(row1,row2)
    Returns: an array of equal length containing 1s where the two rows had
    identical elements and 0 otherwise.
    """
    
    if row1.typecode()=='O' or row2.typecode=='O':
        cmpvect = N.logical_not(abs(N.array(map(cmp,row1,row2)))) # cmp fcn gives -1,0,1
    else:
        cmpvect = N.equal(row1,row2)
    return cmpvect    

def rowsame(row1, row2):
    """Compares two rows from an array, regardless of whether it is an
    array of numbers or of python objects (which requires the cmp function).
    
    Format:  rowsame(row1,row2)
    Returns: 1 if the two rows are identical, 0 otherwise.
    """
    
    cmpval = N.alltrue(rowcompare(row1,row2))
    return cmpval

def sortrows(a,axis=0):
    """Sorts a 2D array, keeping rows together (i.e., assumes data are in
    rows.
    
    Returns: sorted version of a
    """
    
    #    def _sort(s1, s2):
    #        for i in range(s1.shape[0]):
    #            if cmp(s1[i], s2[i]) != 0:
    #                return cmp(s1[i], s2[i])
    #        return 0
    x = asarray(a)
    if axis != 0:
        x = swapaxes(x, axis, 0)
    l = list(x)
    l.sort()           # or l.sort(_sort)
    y = array(l)
    if axis != 0:
        y = swapaxes(y, axis, 0)
    return y
    
