
"""\npstat.py module

#################################################
#######  Written by:  Gary Strangman  ###########
#######  Last modified:  July 25, 1998  #########
#################################################

This module provides some useful list and array manipulation routines
modeled after those found in the |Stat package by Gary Perlman, plus a
number of other useful list/file manipulation functions.  The list-based
functions include:

      abut (source,*args):  
      simpleabut (source, addon):
      colex (listoflists,cnums):
      collapse (listoflists,keepcols,collapsecols,sterr=None,ns=None):
      dm (listoflists,criterion):
      get (namepattern,verbose=1):
      getstrings (namepattern):
      linexand (listoflists,columnlist,valuelist):
      linexor (listoflists,columnlist,valuelist):
      linedelimited (inlist,delimiter):
      lineincols (inlist,colsize):
      lineincustcols (inlist,colsizes):
      list2string (inlist):
      printcc (lst,extra=2):
      printincols (listoflists,colsize):
      pl (listoflists):
      printl(listoflists):
      put (outlist,fname,writetype='w'):
      replace (lst,oldval,newval):
      recode (inlist,listmap,cols='all'):
      remap (listoflists,criterion):
      sortby(listoflists,sortcols):
      unique (inlist):
      writecc (listoflists,file,extra=2,writetype='w'):
      writedelimited (listoflists, delimiter, file, writetype='w'):

Some of these functions have alternate versions which are defined only if
Numeric (NumPy) can be imported.  These functions are generally named as
above, with an 'a' prefix.

      aabut (source, *args):
      acolex (a,indices,axis=1):
      acollapse (a,keepcols,collapsecols,sterr=0,ns=0):
      adm (a,criterion):
      aget (namepattern,verbose=1):
      alinexand (a,columnlist,valuelist):
      alinexor (a,columnlist,valuelist):
      aput (outarray,fname,writetype='w'):
      arecode (a,listmap,col='all'):
      arowcompare (row1, row2):
      arowsame (row1, row2):
      aunique(inarray):

Currently, the code is all but completely un-optimized.  In many cases, the
array versions of functions amount simply to aliases to built-in array
functions/methods.  Their inclusion here is for function name consistency.
Additions, suggestions, or comments are welcome (strang@nmr.mgh.harvard.edu).
"""

import string, sys, os, copy, math, stats, regex
from types import *



###===========================  LIST FUNCTIONS  ==========================
###
### Here are the list functions, defined for all systems.
### Array functions (for NumPy-enabled computers) appear below.
###

def abut (source,*args):
    """\nLike the |Stat abut command.  It concatenates two lists column-wise
and returns the result.  2D lists are also accomodated for either argument
(source or addon).  CAUTION:  If one list is shorter, it will be repeated
until it is as long as the longest list.

Format:  abut(source, args)   where args=any # of lists
Returns: a list of lists as long as the LONGEST list past, source on the
         'left', lists in <args> attached on the 'right'.\n"""

    if type(source) not in [ListType,TupleType]:
        source = [source]
    for addon in args:
        if type(addon) not in [ListType,TupleType]:
            addon = [addon]
        if len(addon) < len(source):            # is source list longer?
            if len(source) % len(addon) == 0:   # are they integer multiples?
                repeats = len(source)/len(addon)    # repeat addon n times
                origadd = copy.deepcopy(addon)
                for i in range(repeats-1):
                    addon = addon + origadd
            else:
                repeats = len(source)/len(addon)+1  # repeat addon x times,
                origadd = copy.deepcopy(addon)      #    x is NOT an integer
                for i in range(repeats-1):
                    addon = addon + origadd
                    addon = addon[0:len(source)]
        elif len(source) < len(addon):          # is addon list longer?
            if len(addon) % len(source) == 0:   # are they integer multiples?
                repeats = len(addon)/len(source)    # repeat source n times
                origsour = copy.deepcopy(source)
                for i in range(repeats-1):
                    source = source + origsour
            else:
                repeats = len(addon)/len(source)+1  # repeat source x times,
                origsour = copy.deepcopy(source)    #   x is NOT an integer
                for i in range(repeats-1):
                    source = source + origsour
                source = source[0:len(addon)]

        source = simpleabut(source,addon)
    return source


def simpleabut (source, addon):
    """\nConcatenates two lists as columns and returns the result.  2D lists
are also accomodated for either argument (source or addon).  This DOES NOT
repeat either list to make the 2 lists of equal length.  Beware of list pairs
with different lengths ... the resulting list will be the length of the
SHORTEST list passed.

Format:  simpleabut(source,addon)  where source, addon=list (or list-of-lists)
Returns: a list of lists as long as the SHORTER of source and addon, with
         source on the 'left' and addon on the 'right'.\n"""

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


def colex (listoflists,cnums):
    """\nExtracts from listoflists the columns specified in the list 'cnums'
(cnums can be an integer, a sequence of integers, or an expression that
corresponds to a slice operation on the variable x ... e.g., x[3:] will colex
columns 3 onward from the listoflists).

Format:  colex (listoflists,cnums)
Returns: a lists of lists corresponding to the columns from listoflists
         specified by cnums, in the order the column numbers appear in cnums.\n"""

    global index
    column = 0
    if type(cnums) in [ListType,TupleType]:   # if multiple columns to get
        index = cnums[0]
        column = map(lambda x: x[index], listoflists)
        for col in cnums[1:]:
            index = col
            column = abut(column,map(lambda x: x[index], listoflists))
    elif type(cnums) == StringType:           # if an 'x[3:]' type expr.
        execstring = 'column = map(lambda x: x'+cnums+', listoflists)'
        exec(execstring)
    else:                                     # else it's just 1 col to get
        index = cnums
        column = map(lambda x: x[index], listoflists)
    return column


def collapse (listoflists,keepcols,collapsecols,sterr=None,ns=None):
    """\nAverages data in collapsecol, keeping all unique items in keepcols (using
unique, which keeps unique LISTS of column numbers), retaining the unique
sets of values in keepcols, the mean for each.  Setting a value for sterr
and/or ns will retain the sterr and N of the calculated means.

Format:  collapse (listoflists,keepcols,collapsecols,sterr=None,ns=None)
Returns: a list of lists with all unique permutations of entries appearing in
         columns ("conditions") specified by keepcols, abutted with the mean
         of each column specified by collapsecols.\n"""

    if type(keepcols) not in [ListType,TupleType]:
        keepcols = [keepcols]
    if type(collapsecols) not in [ListType,TupleType]:
        collapsecols = [collapsecols]
    if keepcols == []:
        means = [0]*len(collapsecols)
        for i in range(len(collapsecols)):
            avgcol = colex(listoflists,collapsecols[i])
            means[i] = stats.mean(avgcol)
        return means
    else:
        values = colex(listoflists,keepcols)
        uniques = unique(values)
        uniques.sort()
        newlist = []
        if type(keepcols) not in [ListType,TupleType]:  keepcols = [keepcols]
        for item in uniques:
            if type(item) not in [ListType,TupleType]:  item =[item]
            tmprows = linexand(listoflists,keepcols,item)
            for col in collapsecols:
                avgcol = colex(tmprows,col)
                item.append(stats.mean(avgcol))
                if sterr:
                    if len(avgcol)>1:
                        item.append(stats.sterr(avgcol))
                    else:
                        item.append('N/A')
                if ns:
                    item.append(len(avgcol))
                newlist.append(item)
        return newlist


def dm (listoflists,criterion):
    """\nReturns rows from the passed list of lists that meet the criteria in
the passed criterion expression (a string as a function of x; e.g., 'x[3]>=9'
will return all rows where the 4th column>=9 and "x[2]=='N'" will return rows
with column 2 equal to the string 'N').

Format:  dm (listoflists, criterion)
Returns: rows from listoflists that meet the specified criterion.\n"""

    function = 'lines = filter(lambda x: '+criterion+',listoflists)'
    exec(function)
    return lines


def get (namepatterns,verbose=1):
    """Loads a list of lists from text files (specified by a UNIX-style
wildcard filename pattern) and converts all numeric values to floats.  Uses
the glob module for filename pattern conversion.

Format:  get (namepatterns,verbose=1)
Returns: a 1D or 2D list of lists from whitespace delimited text files
         specified by namepatterns; numbers that can be converted to floats
         are so converted\n"""

    import glob

    fnames = glob.glob(namepatterns)
    if len(fnames) == 0:
        if verbose:
            print 'NO FILENAMES MATCH PATTERN !!'
        return None

    print fnames                        # so user knows what has been loaded
    elements = []
    for i in range(len(fnames)):
        file = open(fnames[i])
        newelements = map(string.split,file.readlines())
        for i in range(len(newelements)):
            for j in range(len(newelements[i])):
                try:
                    newelements[i][j] = string.atof(newelements[i][j])
                except ValueError:
                    pass
        elements = elements + newelements
    if len(elements)==1:  elements = elements[0]
    return elements


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


def linexand (listoflists,columnlist,valuelist):
    """\nReturns the rows of a list of lists where col (from columnlist) = val
(from valuelist).  One value is required for each column in columnlist.

Format:  linexand (listoflists,columnlist,valuelist)
Returns: the rows of listoflists where columnlist[i]=valuelist[i] for ALL i \n"""

    if type(columnlist) not in [ListType,TupleType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType]:
        valuelist = [valuelist]
    criterion = ''
    for i in range(len(columnlist)):
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' and'
    criterion = criterion[0:-3]         # remove the "and" after the last crit
    function = 'lines = filter(lambda x: '+criterion+',listoflists)'
    exec(function)
    return lines


def linexor (listoflists,columnlist,valuelist):
    """\nReturns the rows of a list of lists where col (from columnlist) = val
(from valuelist).  One value is required for each column in columnlist.  If
only one value exists for columnlist but multiple values appear in valuelist,
the valuelist values are all assumed to pertain to the same column.

Format:  linexor (listoflists,columnlist,valuelist)
Returns: the rows of listoflists where columnlist[i]=valuelist[i] for ANY i \n"""

    if type(columnlist) not in [ListType,TupleType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType]:
        valuelist = [valuelist]
    criterion = ''
    if len(columnlist) == 1 and len(valuelist) > 1:
        columnlist = columnlist*len(valuelist)
    for i in range(len(columnlist)):          # build an exec string
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' or'
    criterion = criterion[0:-2]         # remove the "or" after the last crit
    function = 'lines = filter(lambda x: '+criterion+',listoflists)'
    exec(function)
    return lines


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


def put (outlist,fname,writetype='w'):
    """\nWrites a passed mixed-type (str and/or numbers) list to an output
file, and then closes the file.  Default is overwrite the destination file.

Format:  put (outlist,fname,writetype='w')
Returns: None\n"""

    if type(outlist) in [N.ArrayType]:
        aput(outlist,fname,writetype)
        return
    if type(outlist[0]) not in [ListType,TupleType]:  # 1D list
        outfile = open(fname,writetype)
        outlist = list2string(outlist)
        outfile.write(outlist)
        outfile.write('\n')
        outfile.close()
    else:                                             # 2D list (list-of-lists)
        outfile = open(fname,writetype)
        for row in outlist:
            outfile.write(list2string(row))
            outfile.write('\n')
        outfile.close()
    return None


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


def recode (inlist,listmap,cols='all'):
    """\nChanges the values in a list to a new set of values (useful when
you need to recode data from (e.g.) strings to numbers.\n

Format:  recode (inlist,listmap,cols='all')  cols=recode cols, listmap=2D list
Returns: inlist with the appropriate values replaced with new ones"""

    lst = copy.deepcopy(inlist)
    if cols != 'all':
        if type(cols) not in [ListType,TupleType]:
            cols = [cols]
        for col in cols:
            for row in range(len(lst[col])):
                try:
                    idx = colex(listmap,0).index(lst[row][col])
                    lst[row][col] = listmap[idx][1]
                except ValueError:
                    pass
    else:
        for row in range(len(lst)):
            for col in range(len(lst[row])):
                try:
                    idx = colex(listmap,0).index(lst[row][col])
                    lst[row][col] = listmap[idx][1]
                except ValueError:
                    pass
    return lst


def remap (listoflists,criterion):

    function = 'lines = map(lambda x: '+criterion+',listoflists)'
    exec(function)
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


def unique (inlist):
    """\nReturns all unique items in the passed list.  If the a list-of-lists
is passed, unique LISTS are found (i.e., items in the first dimension are
compared).

Format:  unique (inlist)
Returns: the unique elements (or rows) in inlist\n"""

    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques


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


try:                         # DEFINE THESE *ONLY* IF NUMERIC IS AVAILABLE
 import Numeric
 N = Numeric

 def aabut (source, *args):
    """\nLike the |Stat abut command.  It concatenates two arrays column-wise
and returns the result.  CAUTION:  If one array is shorter, it will be
repeated until it is as long as the other.

Format:  aabut (source, args)    where args=any # of arrays
Returns: an array as long as the LONGEST array past, source appearing on the
         'left', arrays in <args> attached on the 'right'.\n"""

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


 def acolex (a,indices,axis=1):
    """\nExtracts specified indices (a list) from passed array, along passed
axis (column extraction is default).  BEWARE: A 1D array is presumed to be a
column-array (and that the whole array will be returned as a column).

Format:  acolex (a,indices,axis=1)
Returns: the columns of a specified by indices\n"""

    if type(indices) not in [ListType,TupleType,N.ArrayType]:
        indices = [indices]
    if len(N.shape(a)) == 1:
        cols = N.resize(a,[a.shape[0],1])
    else:
        cols = N.take(a,indices,axis)
    return cols


 def acollapse (a,keepcols,collapsecols,sterr=0,ns=0):
    """\nAverages data in collapsecol, keeping all unique items in keepcols
(using unique, which keeps unique LISTS of column numbers), retaining
the unique sets of values in keepcols, the mean for each.  If the sterr or
N of the mean are desired, set either or both parameters to 1.

Format:  acollapse (a,keepcols,collapsecols,sterr=0,ns=0)
Returns: unique 'conditions' specified by the contents of columns specified
         by keepcols, abutted with the mean(s) of column(s) specified by
         collapsecols\n"""

    if keepcols == []:
        avgcol = acolex(a,collapsecols)
        means = N.sum(avgcol)/float(len(avgcol))
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
            tmprows = alinexand(a,keepcols,item)
            for col in collapsecols:
                avgcol = acolex(tmprows,col)
                item.append(stats.mean(avgcol))
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


 def adm (a,criterion):
    """\nReturns rows from the passed list of lists that meet the criteria in
the passed criterion expression (a string).

Format:  adm (a,criterion)   where criterion is like 'x[2]==37'\n"""

    function = 'lines = filter(lambda x: '+criterion+',a)'
    exec(function)
    try:
        lines = N.array(lines)
    except:
        lines = N.array(lines,'O')
    return lines


 def isstring(x):
    if type(x)==StringType:
        return 1
    else:
        return 0


 def aget (namepattern):
    """\nLoads an array from 2D text files (specified by a UNIX-style
wildcard filename pattern).  ONLY 'GET' FILES WITH EQUAL NUMBERS OF COLUMNS
ON EVERY ROW (otherwise returned array will be zero-dimensional).

Format:  aget (namepattern)
Returns: an array of integers, floats or objects (type='O'), depending on the
         contents of the files specified by namepattern\n"""

    import glob

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


 def alinexand (a,columnlist,valuelist):
    """\nReturns the rows of an array where col (from columnlist) = val
(from valuelist).  One value is required for each column in columnlist.

Format:  alinexand (a,columnlist,valuelist)
Returns: the rows of a where columnlist[i]=valuelist[i] for ALL i\n"""

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


 def alinexor (a,columnlist,valuelist):
    """\nReturns the rows of an array where col (from columnlist) = val
(from valuelist).  One value is required for each column in columnlist.
The exception is if either columnlist or valuelist has only 1 value, in
which case that item will be expanded to match the length of the other list.

Format:  alinexor (a,columnlist,valuelist)
Returns: the rows of a where columnlist[i]=valuelist[i] for ANY i\n"""

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


 def aput (outarray,fname,writetype='w'):
    """\nSends passed 1D or 2D array to an output file and closes the file.

Format:  aput (outarray,fname,writetype='w')
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
    return None


 def arecode (a,listmap,col='all'):
    """\nRemaps the values in an array to a new set of values (useful when
you need to recode data from (e.g.) strings to numbers as most stats packages
require.  Can work on SINGLE columns, or 'all' columns at once.

Format:  arecode (a,listmap,col='all')
Returns: a version of array a where listmap[i][0] = (instead) listmap[i][1]\n"""

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


 def arowcompare(row1, row2):
    """\nCompares two rows from an array, regardless of whether it is an
array of numbers or of python objects (which requires the cmp function).

Format:  arowcompare(row1,row2)
Returns: an array of equal length containing 1s where the two rows had
         identical elements and 0 otherwise.\n"""

    if row1.typecode()=='O' or row2.typecode=='O':
        cmpvect = N.logical_not(abs(N.array(map(cmp,row1,row2)))) # cmp fcn gives -1,0,1
    else:
        cmpvect = N.equal(row1,row2)
    return cmpvect


 def arowsame(row1, row2):
    """\nCompares two rows from an array, regardless of whether it is an
array of numbers or of python objects (which requires the cmp function).

Format:  arowsame(row1,row2)
Returns: 1 if the two rows are identical, 0 otherwise.\n"""

    cmpval = N.alltrue(arowcompare(row1,row2))
    return cmpval


 def asortrows(a,axis=0):
    """\nSorts a 2D array, keeping rows together (i.e., assumes data are in
rows.

Format:  asortrows(a,axis=0)
Returns: sorted version of a\n"""

#    def _sort(s1, s2):
#        for i in range(s1.shape[0]):
#            if cmp(s1[i], s2[i]) != 0:
#                return cmp(s1[i], s2[i])
#        return 0

    if axis != 0:
        x = swapaxes(x, axis, 0)
    l = list(x)
    l.sort()           # or l.sort(_sort)
    y = array(l)
    if axis != 0:
        y = swapaxes(y, axis, 0)
    return y


 def aunique(inarray):
    """\nReturns unique items in the FIRST dimension of the passed array. Only
works on arrays NOT including string items (e.g., type 'O' or 'c').

Format:  aunique (inarray)\n"""

    uniques = N.array([inarray[0]])
    if len(uniques.shape) == 1:            # IF IT'S A 1D ARRAY
        for item in inarray[1:]:
            if N.add.reduce(N.equal(uniques,item).flat) == 0:
                try:
                    uniques = N.concatenate([uniques,N.array[N.NewAxis,:]])
                except TypeError:
                    uniques = N.concatenate([uniques,N.array([item])])
    else:                                  # IT MUST BE A 2+D ARRAY
        if inarray.typecode() != 'O':  # not an Object array
            for item in inarray[1:]:
                if not N.sum(N.alltrue(N.equal(uniques,item),1)):
                    try:
                        uniques = N.concatenate( [uniques,item[N.NewAxis,:]] )
                    except TypeError:    # the item to add isn't a list
                        uniques = N.concatenate([uniques,N.array([item])])
                else:
                    pass  # this item is already in the uniques array
        else:   # must be an Object array, alltrue/equal functions don't work
            for item in inarray[1:]:
                newflag = 1
                for unq in uniques:  # NOTE: cmp --> 0=same, -1=<, 1=>
                    test = N.sum(abs(N.array(map(cmp,item,unq))))
                    if test == 0:   # if item identical to any 1 row in uniques
                        newflag = 0 # then not a novel item to add
                        break
                if newflag == 1:
                    try:
                        uniques = N.concatenate( [uniques,item[N.NewAxis,:]] )
                    except TypeError:    # the item to add isn't a list
                        uniques = N.concatenate([uniques,N.array([item])])
    return uniques



except ImportError:    # IF NUMERIC ISN'T AVAILABLE, SKIP ALL arrayfuncs
 pass
