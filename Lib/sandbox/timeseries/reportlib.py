"""
Reporting functions

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com
:version: $Id: tdates.py 2641 2007-01-30 18:40:17Z mattknox_ca $

Ideas borrowed from:

- George Sakkis
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/267662

- Mike Brown
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/148061
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: mattknox_ca $)"
__version__ = '1.0'
__revision__ = "$Revision: 2641 $"
__date__     = '$Date: 2007-01-30 13:40:17 -0500 (Tue, 30 Jan 2007) $'

import sys
import cStringIO, operator, types, copy
import tseries as ts
import tdates as td

__all__ = [
    'Report', 'wrap_onspace', 'wrap_onspace_strict',
    'wrap_always']

class fmtfunc_wrapper:
    """wraps a formatting function such that it handles masked values

:IVariables:
    - `fmtfunc` : formatting function.
    - `mask_rep` : string to use for masked values
    """
    def __init__ (self, fmtfunc, mask_rep):
        if fmtfunc is None:
            self.f = str
        else:
            self.f = fmtfunc
        self.mr = mask_rep

    def __call__ (self, item):
        "Execute the call behavior."

        if hasattr(item, "_mask") and item._mask:
            return self.mr
        else:
            return self.f(item)


_default_options = {
    'dates':None,
    'header_row':None,
    'header_char':'-',
    'row_char':None,
    'footer_label':None,
    'footer_char':'-',
    'footer_func':None,
    'delim':' | ',
    'justify':None,
    'prefix':'',
    'postfix':'',
    'mask_rep':'--',
    'datefmt':None,
    'fmtfunc':str,
    'wrapfunc':lambda x:x,
    'col_width':None,
    'nls':'\n',
    'output':sys.stdout,
    'fixed_width':True
}

class Report(object):
    """Create a tabular TimeSeries report with dates in the left column.
All instance variables are optional and simply serve as the defaults when calling
the report. Parameters for calling the report are the exact same as for
initialization. When calling the report, new options specified will not be saved
to the instance.

:IVariables:
        - `*tseries` : time series objects. Must all be at the same frequency, but
          do not need to be aligned.
          
        - `dates` (DateArray, *[None]*) : dates at which values of all the series
          will be output. If not specified, data will be output from the minimum
          start_date to the maximum end_date of all the time series objects
          
        - `header_row` (list, *[None]*) : List of column headers. Specifying
          the header for the date column is optional.
          
        - `header_char` (string, *['-']*): Character to be used for the row separator
          line between the header and first row of data. None for no separator. This
          is ignored if `header_row` is None.
          
        - `row_char` (string, *[None]*): Character to be used for the row separator
          line between each row of data. None for no separator
          
        - `footer_func` (List of functions or single function, *[None]*) : A function or
          list of functions for summarizing each data column in the report. For example,
          ma.sum to get the sum of the column. If a list of functions is provided
          there must be exactly one function for each column. Do not specify a function
          for the Date column.
          
        - `footer_char` (string, *['-']*): Character to be used for the row separator
          line between the last row of data and the footer. None for no separator. This
          is ignored if `footer_func` is None.
          
        - `footer_label` (string, *[None]*) : label for the footer row. This goes at the
          end of the date column. This is ignored if footer_func is None.
          
        - `justify` (List of strings or single string, *[None]*) : Determines how are
          data justified in their column. If not specified, the date column and string
          columns are left justified, and everything else is right justified. If a
          string is specified, it must be one of 'left', 'right', or 'center' and all
          columns will be justified the same way. If a list is specified, each column
          will be justified according to the specification for that column in the list
          Specifying the justification for the date column is optional.
          
        - `prefix` (string, *['']*) : A string prepended to each printed row.
        
        - `postfix` (string, *['']*) : A string appended to each printed row.
        
        - `mask_rep` (string, *['--']*): String used to represent masked values in
          output
          
        - `datefmt` (string, *[None]*) : Formatting string used for displaying the
          dates in the date column. If None, str() is simply called on the dates
          
        - `fmtfunc` (List of functions or single function, *[None]*) : A function or
          list of functions for formatting each data column in the report. If not
          specified, str() is simply called on each item. If a list of functions is
          provided, there must be exactly one function for each column. Do not specify
          a function for the Date column, that is handled by the datefmt argument
          
        - `wrapfunc` (List of functions or single function, *[lambda x:x]*): A function
          f(text) for wrapping text; each element in the column is first wrapped by this
          function. Instances of  wrap_onspace, wrap_onspace_strict, and wrap_always
          (which are part of this module) work well for this. Eg. wrapfunc=wrap_onspace(10)
          If a list is specified, each column will be wrapped according to the
          specification for that column in the list. Specifying a function for the Date
          column is optional
          
        - `col_width` (list of integers or single integer, *[None]*): use this to specify
          a width for all columns (single integer), or each column individually (list
          of integers). The column will be at least as wide as col_width, but may be
          larger if cell contents exceed col_width. If specifying a list, you may
          optionally specify the width for the Date column as the first entry
          
        - `output` (buffer, *[sys.stdout]*): `output` must have a write method.
          
        - `fixed_width` (boolean, *[True]*): If True, columns are fixed width (ie.
          cells will be padded with spaces to ensure all cells in a given column are
          the same width). If False, `col_width` will be ignored and cells will not
          be padded.
          
:Examples:

    import numpy as np
    import timeseries as ts
    import maskedarray as ma
    from timeseries import Report, wrap_onspace

    series1 = ts.time_series(np.random.uniform(-100,100,15), start_date=ts.thisday('b')-15)
    series2 = ts.time_series(np.random.uniform(-100,100,13), start_date=ts.thisday('b')-10)
    series3 = ts.time_series(['string1', 'another string', 'yet another string']*3, start_date=ts.thisday('b')-10)
    
    darray = ts.date_array(start_date=ts.thisday('b')-8, end_date=ts.thisday('b')-3)
    
    txt_o = open('myfile.txt', 'w')
    html_o = open('myfile.html', 'w')
    
    # report containing only numerical series, showing 2 decimal places
    num_report = Report(series1, series2, fmtfunc=lambda x:'%.2f' % x)
    
    # report containing some string and numerical data
    mixed_report = Report(series1, series2, series3)

    # output a csv report suitable for excel to sys.stdout, show masked values as "N/A"
    num_report(delim=', ', mask_rep='N/A')

    # format one column one with 2 decimal places, and column two with 4.
    # Add a sum footer. Write the output to txt_o
    num_report(fmtfunc=[(lambda x:'%.2f' % x), (lambda x:'%.4f' % x)],
                 footer_func=ma.sum, footer_label='sum', output=txt_o)

    # create an html table of the data over a specified range.
    # Wrap text in cells to width 10. Output to html_o
    html_o.write("<table>")
    mixed_report(series1, series2, series3, dates=darray,
               delim="</td><td>", prefix="<tr><td>", postfix="</td></tr>",
               wrapfunc=wrap_onspace(10, nls='<BR>'), output=html_o)
    html_o.write("</table>")"""

    def __init__(self, *tseries, **kwargs):
        
        self.options = {}
        self.tseries = None
        if len(tseries) > 0:
            self.tseries = tseries
        self.options = self.__make_dict(**kwargs)
        
    def __make_dict(self, **kwargs):
    
        option_dict = copy.copy(self.options)
        
        option_list = list(_default_options)

        for x in [kw for kw in option_list if kw in kwargs]:
            option_dict[x] = kwargs.pop(x)
            
        if len(kwargs) > 0:
            raise KeyError("Unrecognized keyword(s): %s" % (", ".join(kwargs.keys())))
            
        return option_dict
        
    def set_series(self, *tseries):
        """set new time series for the report

:Paramaters:
    - `*tseries` : the TimeSeries objects to be used in the report"""
        self.tseries = tseries

    def set_options(self, **kwargs):
        """set new options or modify options in the report

:Paramaters:
    - `**kwargs` : the options to be used in the report. See the __doc__
      string for the Report class for valid options"""
        self.options = self.__make_dict(**kwargs)

    
    def __call__(self, *tseries, **kwargs):
        """generate a report

:Paramaters:
    - `*tseries` : the TimeSeries objects to be used in the report. If
      omitted, the previously set TimeSeries objects will be used
    - `**kwargs` : the options to be used in the report. See the __doc__
      string for the Report class for valid options. If omitted, the
      previously set options will be used"""

        option_dict = self.__make_dict(**kwargs)
        if len(tseries) == 0:
            tseries = self.tseries
        
        def option(kw):
            return option_dict.get(kw, _default_options[kw])
            
        dates = option('dates')
        header_row = option('header_row')
        header_char = option('header_char')
        row_char = option('row_char')
        footer_label = option('footer_label')
        footer_char = option('footer_char')
        footer_func = option('footer_func')
        delim = option('delim')
        justify = option('justify')
        prefix = option('prefix')
        postfix = option('postfix')
        mask_rep = option('mask_rep')
        datefmt = option('datefmt')
        fmtfunc = option('fmtfunc')
        wrapfunc = option('wrapfunc')
        col_width = option('col_width')
        nls=option('nls')
        output=option('output')
        fixed_width=option('fixed_width')
        
        if header_row is not None:
            has_header=True
            if len(header_row) == len(tseries)+1:
                # label for date column included
                rows = [header_row]
            elif len(header_row) == len(tseries):
                # label for date column not included
                rows = [['']+header_row]
            else:
                raise ValueError("mismatch with number of headers and series")
        else:
            has_header=False
            rows=[]

        if fixed_width:
            if justify is not None:
                _justify = kwargs.pop('justify')
                if isinstance(justify, str):
                    # justify all columns the the same way
                    justify = [justify for x in range(len(tseries)+1)]
                elif isinstance(justify, list): #assume it is a list or tuple, etc
                    if len(justify) == len(tseries):
                        # justification for date column not included, so set that
                        # to left by default
                        justify = ['left'] + justify
                else:
                    raise ValueError("invalid `justify` specification")
            else:
                # default column justification
                justify = ['left']
                for ser in tseries:
                    if str(ser.dtype)[:2] == '|S': justify.append('left')
                    else: justify.append('right')
        else:
            justify = [None for x in range(len(tseries)+1)]

        if datefmt is None:
            def datefmt_func(date): return str(date)
        else:
            def datefmt_func(date): return date.strfmt(datefmt)

        if dates is None:
            tseries = ts.align_series(*tseries)
            dates = td.date_array(start_date=tseries[0].start_date,
                                  end_date=tseries[0].end_date)
        else:
            tseries = ts.align_series(start_date=dates[0], end_date=dates[-1], *tseries)

        if isinstance(fmtfunc, list):
            fmtfunc = [fmtfunc_wrapper(f, mask_rep) for f in fmtfunc]
        else:
            fmtfunc = [fmtfunc_wrapper(fmtfunc, mask_rep)]*len(tseries)

        def wrapfunc_default(func):
            if func is None: return lambda x:x
            else: return func
            
        if isinstance(wrapfunc, list):
            if len(wrapfunc) == len(tseries):
                wrapfunc = [lambda x: x] + wrapfunc
            wrapfunc = [wrapfunc_default(func) for func in wrapfunc]
        else:
            wrapfunc = [wrapfunc_default(wrapfunc) for x in range(len(tseries)+1)]
    
            
        if isinstance(col_width, list):
            if len(col_width) == len(tseries):
                col_width = [None] + col_width
        else:
            col_width = [col_width for x in range(len(tseries)+1)]

        ############################################################
        # temporary hack to handle singletons for different types of
        # behaviour until we finalize how they will be handled
        ############################################################
        if ts.time_series([1,2], start_date=ts.thisday('b'))[0].ndim == 0:
            def getval(series, date): return series[date]
        else:
            def getval(series, date):
                temp = series[date]
                if temp is ts.tsmasked:
                    return temp
                else:
                    return temp.series[0]
        ############################################################

        for d in dates:
            rows.append([datefmt_func(d)]+[fmtfunc[i](getval(ser, d)) for i, ser in enumerate(tseries)])

        if footer_func is not None:
            has_footer=True
            if not isinstance(footer_func, list):
                footer_func = [footer_func]*len(tseries)

            if footer_label is None: footer_label = ['']
            else: footer_label = [footer_label]

            footer_data = []
            for i, ser in enumerate(tseries):
                if footer_func[i] is None:
                    footer_data.append('')
                else:
                    footer_data.append(fmtfunc[i](footer_func[i](ser[dates])))

            rows.append(footer_label + footer_data)
        else:
            has_footer=False
            
            
        def rowWrapper(row):
            newRows = [wrapfunc[i](item).split('\n') for i, item in enumerate(row)]
            return [[(substr or '') for substr in item] for item in map(None,*newRows)]
        # break each logical row into one or more physical ones
        logicalRows = [rowWrapper(row) for row in rows]
        numLogicalRows = len(logicalRows)
        # columns of physical rows
        columns = map(None,*reduce(operator.add,logicalRows))
        numCols = len(columns)
        colNums = list(range(numCols))

        # get the maximum of each column by the string length of its items
        maxWidths = [max(col_width[i], *[len(str(item)) for item in column])
                        for i, column in enumerate(columns)]

        def getSeparator(char, separate):
            if char is not None and separate:
                return char * (len(prefix) + len(postfix) + sum(maxWidths) + \
                                             len(delim)*(len(maxWidths)-1))
            else:
                return None

        header_separator = getSeparator(header_char, has_header)
        footer_separator = getSeparator(footer_char, has_footer)
        row_separator = getSeparator(row_char, True)

        # select the appropriate justify method
        justify_funcs = {'center':str.center, 'right':str.rjust, 'left':str.ljust,
                          'none':(lambda text, width: text)}

        if has_header and has_footer:
            data_start = 1
            data_end = numLogicalRows-3
        elif has_header:
            data_start = 1
            data_end = numLogicalRows-2
        elif has_footer:
            data_start = 0
            data_end = numLogicalRows-3
        else:
            data_start = 0
            data_end = numLogicalRows-2

        for rowNum, physicalRows in enumerate(logicalRows):
            for row in physicalRows:
                output.write(prefix \
                           + delim.join([justify_funcs[str(justify[colNum]).lower()](str(item),width) for (colNum,item,width) in zip(colNums,row,maxWidths)]) \
                           + postfix + nls)

            if row_separator and (data_start <= rowNum <= data_end):
                output.write(row_separator + nls)
            elif header_separator and rowNum < data_start:
                output.write(header_separator + nls)
            elif footer_separator and rowNum == data_end + 1:
                output.write(footer_separator + nls)


class wrap_onspace(object):
    """A callable word-wrap class that preserves existing line breaks
and most spaces in the text.
    
:IVariables:
    - `width` (int): width to wrap at. Won't split up words wider than `width`
    - `nls` (str, *['\n']*): New line separator. Assumes existing line
      breaks use this new line separator as well.
      
:Parameters (for __call__ method):
    - `text` (str): text to wrap"""

    def __init__(self, width, nls='\n'):
        self.width = width
        self.nls = nls
        
    def __call__(self, text):
    
        width = self.width
        nls = self.nls

        def break_or_space(line, word, width):
            temp_idx = (len(line[line.rfind(nls)+1:]) + len(word.split(nls,1)[0]) >= width)
            if temp_idx:
                return nls
            else:
                return ' '

        return reduce(lambda line, word, width=width: '%s%s%s' %
                      (line,
                       break_or_space(line, word, width),
                       word),
                      text.split(' ')
                     )    


import re
class wrap_onspace_strict(object):
    """A callable word-wrap class similar to wrap_onspace, but
enforces the width constraint: words longer than width are split.
    
:IVariables:
    - `width` (int): width to wrap at. Will split up words wider than `width`
    - `nls` (str, *['\n']*): New line separator. Assumes existing line
      breaks use this new line separator as well.
      
:Parameters (for __call__ method):
    - `text` (str): text to wrap"""

    def __init__(self, width, nls='\n'):
        self.width = width
        self.nls = nls
        
    def __call__(self, text):
    
        width = self.width
        nls = self.nls

        wordRegex = re.compile(r'\S{'+str(width)+r',}')
        return wrap_onspace(wordRegex.sub(lambda m: wrap_always(m.group(),width, nls=nls),text),width, nls=nls)


import math
class wrap_always(object):
    """A callable word-wrap class that wraps text on exactly width
characters. It doesn't split the text into words.
    
:IVariables:
    - `width` (int): width to wrap at.
    - `nls` (str, *['\n']*): New line separator.
      
:Parameters (for __call__ method):
    - `text` (str): text to wrap"""

    def __init__(self, width, nls='\n'):
        self.width = width
        self.nls = nls
        
    def __call__(self, text):
    
        width = self.width
        nls = self.nls
        return nls.join([ text[width*i:width*(i+1)] \
                           for i in xrange(int(math.ceil(1.*len(text)/width))) ])
