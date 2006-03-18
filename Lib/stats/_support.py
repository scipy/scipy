from numpy import asarray
import stats
import numpy as N
from types import ListType, TupleType, StringType
import copy

def abut (source, *args):
    """\nLike the |Stat abut command.  It concatenates two arrays column-wise
    and returns the result.  CAUTION:  If one array is shorter, it will be
    repeated until it is as long as the other.

    Format:  abut (source, args)    where args=any # of arrays
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


def unique(inarray):
    """Returns unique items in the FIRST dimension of the passed array. Only
    works on arrays NOT including string items (e.g., type 'O' or 'c').
    """
    inarray = asarray(inarray)
    uniques = N.array([inarray[0]])
    if len(uniques.shape) == 1:            # IF IT'S A 1D ARRAY
        for item in inarray[1:]:
            if N.add.reduce(N.equal(uniques,item).flat) == 0:
                try:
                    uniques = N.concatenate([uniques,N.array[N.newaxis,:]])
                except TypeError:
                    uniques = N.concatenate([uniques,N.array([item])])
    else:                                  # IT MUST BE A 2+D ARRAY
        if inarray.typecode() != 'O':  # not an Object array
            for item in inarray[1:]:
                if not N.sum(N.alltrue(N.equal(uniques,item),1)):
                    try:
                        uniques = N.concatenate( [uniques,item[N.newaxis,:]] )
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
                        uniques = N.concatenate( [uniques,item[N.newaxis,:]] )
                    except TypeError:    # the item to add isn't a list
                        uniques = N.concatenate([uniques,N.array([item])])
    return uniques

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

def printcc (lst,extra=2):
    """\nPrints a list of lists in columns, customized by the max size of items
within the columns (max size of items in col, plus 'extra' number of spaces).
Use 'dashes' or '\n' in the list(oflists) to print dashes or blank lines,
respectively.

Format:  printcc (lst,extra=2)
Returns: None\n"""

    def makestr (x):
        if type(x) != StringType:
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
                if stderr:
                    if len(avgcol)>1:
                        item.append(stats.stderr(avgcol))
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


def makestr (item):
    if type(item) != StringType:
        item = str(item)
    return item

def lineincustcols (inlist,colsizes):
    """\nReturns a string composed of elements in inlist, with each element
right-aligned in a column of width specified by a sequence colsizes.  The
length of colsizes must be greater than or equal to the number of columns in
inlist.

Format:  lineincustcols (inlist,colsizes)
Returns: formatted string created from inlist\n"""

    outstr = ''
    for i in range(len(inlist)):
        if type(inlist[i]) != StringType:
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

    stringlist = map(makestr,inlist)
    return "".join(stringlist)
