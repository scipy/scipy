import cStringIO,operator, types
import timeseries as ts

class fmtfunc_wrapper:
    """wraps a formatting function such that it handles masked values

:IVariables:
    - `fmtfunc` : formatting function.
    - `mask_rep` : string to use for masked values
    """
    def __init__ (self, fmtfunc, mask_rep):
        self.f = fmtfunc
        self.mr = mask_rep

    def __call__ (self, item):
        "Execute the call behavior."

        if hasattr(item, "_mask") and item._mask:
            return self.mr
        else:
            return self.f(item)



def report(*tseries, **kwargs):
    """generate a table report of *tseries with dates in the left column.

:Parameters:
        - `*tseries` : time series objects. Must all be at the same frequency, but
          do not need to be aligned.
        - `dates` (DateArray, *[None]*) : dates at which values of all the series
          will be output. If not specified, data will be output from the minimum
          start_date to the maximum end_date of all the time series objects
        - `header_row` (list, *[None]*) : optional list of column headers. Length
          must be equal to len(tseries) (no date column header specified) or
          len(tseries)+1 (first header is assumed to be date column header)
        - `header_char` : Character to be used for the row separator line
        - `justify` (List of strings or single string, *[None]*) : Determines how are
          data justified in their column. If not specified, the date column and string
          columns are left justified, and everything else is right justified. If a
          string is specified, it must be one of 'left', 'right', or 'center' and all
          columns will be justified the same way. If a list is specified, each column
          will be justified according to the specification for that column in the list
          (specifying the justification for the date column is optional).
        - `prefix` (string, *['']*) : A string prepended to each printed row.
        - `postfix` (string, *['']*) : A string appended to each printed row.
        - `mask_rep` (string, *['--']*): String used to represent masked values in
          output
        - `datefmt` (string, *[None]*) : Formatting string used for displaying the
          dates in the date column. If None, str() is simply called on the dates
        - `fmtfunc` (List of functions or single function, *[None]*) : A function or
          list of functions for formatting each data column in the report. If not
          specified, str() is simply called on each item. If a list of functions is
          provided, there must be exactly one function for each column
        - `wrapfunc` (function, *[lambda x:x]*): A function f(text) for wrapping text;
          each element in the table is first wrapped by this function. Useful functions
          for this are wrap_onspace, wrap_onspace_strict, and wrap_always (which are
          part of this module). Eg wrapfunc=lambda x: wrap_onspace(x, 10)
          
:Examples:

    import numpy as np
    import timeseries as ts
    from timeseries.io import report as r

    series1 = ts.time_series(np.random.uniform(-100,100,15), start_date=ts.thisday('b')-15)
    series2 = ts.time_series(np.random.uniform(-100,100,13), start_date=ts.thisday('b')-10)
    series3 = ts.time_series(['string1', 'another string', 'yet another string']*3, start_date=ts.thisday('b')-10)
    
    darray = ts.date_array(start_date=ts.thisday('b')-8, end_date=ts.thisday('b')-3)
    
    # print all values of series1 and series2 and show 2 decimal places.
    # show masked values as "N/A"
    print r.report(series1, series2, fmtfunc=lambda x:'%.2f' % x, mask_rep='N/A')

    # same thing, but format one column one with 2 decimal places, and column two with 4
    print r.report(series1, series2, fmtfunc=[(lambda x:'%.2f' % x), (lambda x:'%.4f' % x)], mask_rep='N/A')

    # print an html table of the data over a specified range
    print "<table>" + \
          r.report(series1, series2, series3, dates=darray,
                   delim="</td><td>", prefix="<tr><td>", postfix="</td></tr>") + \
          "</table>"
    
    # print a table with columns 10 characters wide when possible, but don't break up a word
    print r.report(series1, series3, series2, wrapfunc=lambda x: r.wrap_onspace(x, 10))"""
    
    dates = kwargs.pop('dates', None)
    header_row = kwargs.pop('header_row', None)
    header_char = kwargs.pop('header_char', '-')
    delim = kwargs.pop('delim', ' | ')
    justify = kwargs.pop('justify', None)
    prefix = kwargs.pop('prefix', '')
    postfix = kwargs.pop('postfix', '')
    mask_rep = kwargs.pop('mask_rep', '--')
    datefmt = kwargs.pop('datefmt', None)
    fmtfunc = kwargs.pop('fmtfunc', str)
    
    if type(fmtfunc) != types.ListType:
        fmtfunc = [fmtfunc_wrapper(fmtfunc, mask_rep)]*len(tseries)
    else:
        fmtfunc = [fmtfunc_wrapper(f, mask_rep) for f in fmtfunc]
    
    wrapfunc = kwargs.pop('wrapfunc', lambda x:x)

    if len(kwargs) > 0:
        raise KeyError("Unrecognized keyword(s): %s" % (", ".join(kwargs.keys())))

    if header_row is not None:
        hasHeader=True
        if len(header_row) == len(tseries)+1:
            # label for date column included
            rows = [header_row]
        elif len(header_row) == len(tseries):
            # label for date column not included
            rows = [['']+header_row]
    else:
        hasHeader=False
        rows=[]
        
    if justify is not None:
        if type(justify) == types.StringType:
            # justify all columns the the same way
            justify = [justify for x in range(len(tseries)+1)]
        else: #assume it is a list or tuple, etc
            if len(justify) == len(tseries):
                # justification for date column not included, so set that
                # to left by default
                justify = ['left'] + justify
    else:
        # default column justification
        justify = ['left']
        for ser in tseries:
            if str(ser.dtype)[:2] == '|S': justify.append('left')
            else: justify.append('right')
        
        
    if datefmt is None:
        def datefmt_func(date): return str(date)
    else:
        def datefmt_func(date): return date.strfmt(datefmt)

    tseries = ts.align_series(*tseries)
    
    if dates is None:
        dates = ts.date_array(start_date=tseries[0].start_date,
                              end_date=tseries[0].end_date)
    
    for d in dates:
        rows.append([datefmt_func(d)]+[fmtfunc[i](ser[d]) for i, ser in enumerate(tseries)])

    return indent(rows, hasHeader=hasHeader, headerChar=header_char,
                  delim=delim, justify=justify, separateRows=False,
                  prefix=prefix, postfix=postfix, wrapfunc=wrapfunc)
   



# written by George Sakkis
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/267662
def indent(rows, hasHeader=False, headerChar='-', delim=' | ', justify=None,
           separateRows=False, prefix='', postfix='', wrapfunc=lambda x:x):
    """Indents a table by column.
       - rows: A sequence of sequences of items, one sequence per row.
       - hasHeader: True if the first row consists of the columns' names.
       - headerChar: Character to be used for the row separator line
         (if hasHeader==True or separateRows==True).
       - delim: The column delimiter.
       - justify: Determines how are data justified in their column. 
         Valid values are 'left','right' and 'center'.
       - separateRows: True if rows are to be separated by a line
         of 'headerChar's.
       - prefix: A string prepended to each printed row.
       - postfix: A string appended to each printed row.
       - wrapfunc: A function f(text) for wrapping text; each element in
         the table is first wrapped by this function."""

    
    def rowWrapper(row):
        newRows = [wrapfunc(item).split('\n') for item in row]
        return [[substr or '' for substr in item] for item in map(None,*newRows)]
    # break each logical row into one or more physical ones
    logicalRows = [rowWrapper(row) for row in rows]
    # columns of physical rows
    columns = map(None,*reduce(operator.add,logicalRows))
    numCols = len(columns)
    colNums = list(range(numCols))
    
    if justify is None:
        justify = ['left' for x in range(numCols)]
    
    # get the maximum of each column by the string length of its items
    maxWidths = [max([len(str(item)) for item in column]) for column in columns]
    rowSeparator = headerChar * (len(prefix) + len(postfix) + sum(maxWidths) + \
                                 len(delim)*(len(maxWidths)-1))
    # select the appropriate justify method
    justify_funcs = {'center':str.center, 'right':str.rjust, 'left':str.ljust}
   
    
    output=cStringIO.StringIO()
    if separateRows: print >> output, rowSeparator
    for physicalRows in logicalRows:
        for row in physicalRows:
            print >> output, \
                prefix \
                + delim.join([justify_funcs[justify[colNum].lower()](str(item),width) for (colNum,item,width) in zip(colNums,row,maxWidths)]) \
                + postfix
        if separateRows or hasHeader: print >> output, rowSeparator; hasHeader=False
    return output.getvalue()

# written by Mike Brown
# http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/148061
def wrap_onspace(text, width):
    """
    A word-wrap function that preserves existing line breaks
    and most spaces in the text. Expects that existing line
    breaks are posix newlines (\n).
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line[line.rfind('\n')+1:])
                         + len(word.split('\n',1)[0]
                              ) >= width)],
                   word),
                  text.split(' ')
                 )

import re
def wrap_onspace_strict(text, width):
    """Similar to wrap_onspace, but enforces the width constraint:
       words longer than width are split."""
    wordRegex = re.compile(r'\S{'+str(width)+r',}')
    return wrap_onspace(wordRegex.sub(lambda m: wrap_always(m.group(),width),text),width)

import math
def wrap_always(text, width):
    """A simple word-wrap function that wraps text on exactly width characters.
       It doesn't split the text in words."""
    return '\n'.join([ text[width*i:width*(i+1)] \
                       for i in xrange(int(math.ceil(1.*len(text)/width))) ])
