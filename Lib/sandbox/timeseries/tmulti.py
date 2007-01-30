# pylint: disable-msg=W0201, W0212
"""
Support for multi-variable time series, through masked recarrays.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'


import sys

import numpy
from numpy import bool_, complex_, float_, int_, str_, object_
import numpy.core.fromnumeric as fromnumeric
import numpy.core.numeric as numeric
from numpy.core.numeric import ndarray
import numpy.core.numerictypes as ntypes
import numpy.core.umath as umath
from numpy.core.defchararray import chararray
from numpy.core.records import find_duplicate
from numpy.core.records import format_parser, recarray, record 
from numpy.core.records import fromarrays as recfromarrays

import maskedarray as MA
#import numpy.core.ma as MA
#reload(MA)
#MaskedArray = MA.MaskedArray
from maskedarray.core import MaskedArray, MAError, default_fill_value, \
    masked_print_option
from maskedarray.core import masked, nomask, getmask, getmaskarray, make_mask,\
    make_mask_none, mask_or, masked_array, filled

import maskedarray.mrecords as MR
#reload(MR)
from maskedarray.mrecords import _checknames, _guessvartypes, openfile,\
    MaskedRecords
from maskedarray.mrecords import fromrecords as mrecfromrecords

from tseries import TimeSeries, time_series, _getdatalength
from tdates import Date, DateArray, date_array

#ndarray = numeric.ndarray
_byteorderconv = numpy.core.records._byteorderconv
_typestr = ntypes._typestr

reserved_fields = MR.reserved_fields + ['_dates']

import warnings
                    format='%(name)-15s %(levelname)s %(message)s',)

__all__ = [
'MultiTimeSeries','fromarrays','fromrecords','fromtextfile',           
]

def _getformats(data):
    """Returns the formats of each array of arraylist as a comma-separated 
    string."""
    if isinstance(data, record):
        return ",".join([desc[1] for desc in data.dtype.descr])
    
    formats = ''
    for obj in data:
        obj = numeric.asarray(obj)
#        if not isinstance(obj, ndarray):
##        if not isinstance(obj, ndarray):
#            raise ValueError, "item in the array list must be an ndarray."
        formats += _typestr[obj.dtype.type]
        if issubclass(obj.dtype.type, ntypes.flexible):
            formats += `obj.itemsize`
        formats += ','
    return formats[:-1]    


    


class MultiTimeSeries(TimeSeries, MaskedRecords, object):
    """
    
:IVariables:
    - `__localfdict` : Dictionary
        Dictionary of local fields (`f0_data`, `f0_mask`...)
    - `__globalfdict` : Dictionary
        Dictionary of global fields, as the combination of a `_data` and a `_mask`.
        (`f0`)
    """
    _defaultfieldmask = nomask
    _defaulthardmask = False
    def __new__(cls, data, dates=None, mask=nomask, dtype=None, 
                freq=None, observed=None, start_date=None, 
                hard_mask=False, fill_value=None,
#                offset=0, strides=None,
                formats=None, names=None, titles=None,
                byteorder=None, aligned=False):
        tsoptions = dict(fill_value=fill_value, hard_mask=hard_mask,)
        mroptions = dict(fill_value=fill_value, hard_mask=hard_mask, 
                         formats=formats, names=names, titles=titles, 
                         byteorder=byteorder, aligned=aligned)
        #
        if isinstance(data, MultiTimeSeries):
            cls._defaultfieldmask = data._series._fieldmask
            cls._defaulthardmask = data._series._hardmask | hard_mask
            cls._fill_value = data._series._fill_value
            return data._data.view(cls)
#        elif isinstance(data, TimeSeries):
#            cls._defaultfieldmask = data._series._fieldmask
#            cls._defaulthardmask = data._series._hardmask | hard_mask
#            cls._fill_value = data._series._fill_value
            
        # .......................................
        _data = MaskedRecords(data, mask=mask, dtype=dtype, **mroptions)
        if dates is None:
            length = _getdatalength(data)
            newdates = date_array(start_date=start_date, length=length,
                                  freq=freq)                 
        elif not hasattr(dates, 'freq'):
            newdates = date_array(dlist=dates, freq=freq)
        else:
            newdates = dates
#            _data = data
#            if hasattr(data, '_mask') :
#                mask = mask_or(data._mask, mask)
        cls._defaultdates = newdates    
        cls._defaultobserved = observed  
        cls._defaultfieldmask = _data._fieldmask
#        assert(_datadatescompat(data,newdates))
        #
        return _data.view(cls)
    
        #..................................
    def __array_wrap__(self, obj, context=None):
        """Special hook for ufuncs.
Wraps the numpy array and sets the mask according to context.
        """
#        mclass = self.__class__
        #..........
        if context is None:
#            return mclass(obj, mask=self._mask, copy=False)
            return MaskedArray(obj, mask=self._mask, copy=False,
                               dtype=obj.dtype,
                               fill_value=self.fill_value, )
        #..........
        (func, args) = context[:2]
 
#        return mclass(obj, copy=False, mask=m)
        return MultiTimeSeries(obj, copy=False, mask=m,)
#                           dtype=obj.dtype, fill_value=self._fill_value)
    
        
    def __array_finalize__(self,obj):
        if isinstance(obj, MultiTimeSeries):
            self.__dict__.update(_dates=obj._dates,
                                 _series=obj._series,
                                 _data=obj._series._data,
                                 _fieldmask=obj._series._fieldmask,
                                 _hardmask=obj._series._hardmask,
                                 _fill_value=obj._fill_value                                 
                                 )
        else:     
            self.__dict__.update(_data = obj.view(recarray),
                                 _dates = self._defaultdates,
                                 _series = MaskedRecords(obj, dtype=obj.dtype),
                                 _fieldmask = self._defaultfieldmask,
                                 _hardmask = self._defaulthardmask,
                                 fill_value = self._fill_value
                                )
            MultiTimeSeries._defaultfieldmask = nomask
            MultiTimeSeries._defaulthardmask = False
        return
    #......................................................
    def __getattribute__(self, attr):
        try:
            # Returns a generic attribute
            return object.__getattribute__(self,attr)
        except AttributeError: 
            # OK, so attr must be a field name
            pass
        # Get the list of fields ......
        _names = self.dtype.names
        _local = self.__dict__
        _mask = _local['_fieldmask']
        if attr in _names:
            _data = _local['_data']
            obj = numeric.asarray(_data.__getattribute__(attr)).view(MaskedArray)
            obj._mask = make_mask(_mask.__getattribute__(attr))
            return obj
        elif attr == '_mask':
            if self.size > 1:
                return _mask.view((bool_, len(self.dtype))).all(1)
            return _mask.view((bool_, len(self.dtype)))
        raise AttributeError,"No attribute '%s' !" % attr
            
    def __setattr__(self, attr, val):
        newattr = attr not in self.__dict__
        try:
            # Is attr a generic attribute ?
            ret = object.__setattr__(self, attr, val)
        except:
            # Not a generic attribute: exit if it's not a valid field
            fielddict = self.dtype.names or {}
            if attr not in fielddict:
                exctype, value = sys.exc_info()[:2]
                raise exctype, value
        else:
            if attr not in list(self.dtype.names) + ['_mask']:
                return ret
            if newattr:         # We just added this one
                try:            #  or this setattr worked on an internal
                                #  attribute. 
                    object.__delattr__(self, attr)
                except:
                    return ret
        # Case #1.: Basic field ............
        base_fmask = self._fieldmask
        _names = self.dtype.names
        if attr in _names:
            fval = filled(val)
            mval = getmaskarray(val)
            if self._hardmask:
                mval = mask_or(mval, base_fmask.__getattr__(attr))
            self._data.__setattr__(attr, fval)
            base_fmask.__setattr__(attr, mval)
            return
        elif attr == '_mask':
            if self._hardmask:
                val = make_mask(val)
                if val is not nomask:
#                    mval = getmaskarray(val)
                    for k in _names:
                        m = mask_or(val, base_fmask.__getattr__(k))
                        base_fmask.__setattr__(k, m)
                else:
                return
            else:
                mval = getmaskarray(val)
                for k in _names:
                    base_fmask.__setattr__(k, mval)  
                return
    #............................................
    def __getitem__(self, indx):
        """Returns all the fields sharing the same fieldname base.
    The fieldname base is either `_data` or `_mask`."""
        _localdict = self.__dict__
        # We want a field ........
        if indx in self.dtype.names:
            obj = _localdict['_series'][indx].view(TimeSeries)
            obj._mask = make_mask(_localdict['_fieldmask'][indx])
            return obj
        # We want some elements ..
        indx = super(MultiTimeSeries, self)._TimeSeries__checkindex(indx)
        return MultiTimeSeries(_localdict['_series'][indx], 
                               dates=_localdict['_dates'][indx],
#                              mask=_localdict['_fieldmask'][indx],
                               dtype=self.dtype)
        
    def __getslice__(self, i, j):
        """Returns the slice described by [i,j]."""
        _localdict = self.__dict__
        return MultiTimeSeries(_localdict['_data'][i:j], 
                               mask=_localdict['_fieldmask'][i:j],
                               dates=_localdict['_dates'][i:j],
                               dtype=self.dtype)      
        
    def __setslice__(self, i, j, value):
        """Sets the slice described by [i,j] to `value`."""
        _localdict = self.__dict__
        d = _localdict['_data']
        t = _localdict['_dates']
        m = _localdict['_fieldmask']
        names = self.dtype.names
        if value is masked:
            for n in names:
                m[i:j][n] = masked
        elif not self._hardmask:
            fval = filled(value)
            mval = getmaskarray(value)
            for n in names:
                d[n][i:j] = fval
                m[n][i:j] = mval
        else:
            mindx = getmaskarray(self)[i:j]
            val = masked_array(value, mask=mindx, keep_mask=True)
            valmask = getmask(value)
            if valmask is nomask:
                for n in names:
                    mval = mask_or(m[n][i:j], valmask)
                    d[n][i:j][~mval] = filled(value)
            elif valmask.size > 1:
                for n in names:
                    mval = mask_or(m[n][i:j], valmask)
                    d[n][i:j][~mval] = fval[~mval]
                    m[n][i:j] = mask_or(m[n][i:j], mval) 
            
        return MultiTimeSeries(d, mask=m, dates=t[i:j], dtype=self.dtype)      
            
        
    #......................................................
    def __str__(self):
        """x.__str__() <==> str(x)
Calculates the string representation, using masked for fill if it is enabled. 
Otherwise, fills with fill value.
        """
        if self.size > 1:
            mstr = ["(%s)" % ",".join([str(i) for i in s])  
                    for s in zip(*[getattr(self,f) for f in self.dtype.names])]
            return "[%s]" % ", ".join(mstr)
        else:
            mstr = numeric.asarray(self._data.item(), dtype=object_)
            mstr[list(self._fieldmask)] = masked_print_option
            return str(mstr)
    
    def __repr__(self):
        """x.__repr__() <==> repr(x)
Calculates the repr representation, using masked for fill if it is enabled. 
Otherwise fill with fill value.
        """
        _names = self.dtype.names
        _dates = self._dates
        if numeric.size(_dates) > 2 and self._dates.isvalid():
            timestr = "[%s ... %s]" % (str(_dates[0]),str(_dates[-1]))
        else:
            timestr = str(_dates)
        fmt = "%%%is : %%s" % (max([len(n) for n in _names])+4,)
        reprstr = [fmt % (f,getattr(self,f)) for f in self.dtype.names]
        reprstr.insert(0,'multitimeseries(')
        reprstr.extend([fmt % ('dates', timestr),
                        fmt % ('    fill_value', self._fill_value), 
                         '               )'])
        return str("\n".join(reprstr))
    #......................................................
    def view(self, obj):
        """Returns a view of the mrecarray."""
        try:
            if issubclass(obj, ndarray):
                return ndarray.view(self, obj)
        except TypeError:
            pass
        dtype = numeric.dtype(obj)
        if dtype.fields is None:
            return self.__array__().view(dtype)
        return ndarray.view(self, obj)            
    #............................................
#    def harden_mask(self):
#        "Forces the mask to hard"
#        self._hardmask = True
#    def soften_mask(self):
#        "Forces the mask to soft"
#        self._hardmask = False
    #.............................................
    def copy(self):
        "Returns a copy of the argument."
        _localdict = self.__dict__
        return MultiTimeSeries(_localdict['_data'].copy(),
                       dates=_localdict['_dates'].copy(),
                        mask=_localdict['_fieldmask'].copy(),
                       dtype=self.dtype)
    #.............................................
    def addfield(self, newfield, newfieldname=None):
        MaskedRecords.addfield(self, newfield, newfieldname)

#####---------------------------------------------------------------------------
#---- --- Constructors ---
#####---------------------------------------------------------------------------

def fromarrays(arraylist, dates=None,
               dtype=None, shape=None, formats=None,
               names=None, titles=None, aligned=False, byteorder=None):
    """Creates a mrecarray from a (flat) list of masked arrays.

:Parameters:
    - `arraylist` : Sequence
      A list of (masked) arrays. Each element of the sequence is first converted
      to a masked array if needed. If a 2D array is passed as argument, it is
      processed line by line
    - `dtype` : numeric.dtype
      Data type descriptor.
    - `shape` : Integer *[None]*
      Number of records. If None, `shape` is defined from the shape of the first
      array in the list.
    - `formats` :
      (Description to write)
    - `names` : 
      (description to write)
    - `titles`:
      (Description to write)
    - `aligned`: Boolen *[False]*
      (Description to write, not used anyway)   
    - `byteorder`: Boolen *[None]*
      (Description to write, not used anyway)
       

    """
    arraylist = [MA.asarray(x) for x in arraylist]
    # Define/check the shape.....................
    if shape is None or shape == 0:
        shape = arraylist[0].shape
    if isinstance(shape, int):
        shape = (shape,)
    # Define formats from scratch ...............
    if formats is None and dtype is None:
        formats = _getformats(arraylist)
    # Define the dtype ..........................
    if dtype is not None:
        descr = numeric.dtype(dtype)
        _names = descr.names
    else:
        parsed = format_parser(formats, names, titles, aligned, byteorder)
        _names = parsed._names
        descr = parsed._descr
    # Determine shape from data-type.............
    if len(descr) != len(arraylist):
        msg = "Mismatch between the number of fields (%i) and the number of "\
              "arrays (%i)"
        raise ValueError, msg % (len(descr), len(arraylist))
    d0 = descr[0].shape
    nn = len(d0)
    if nn > 0:
        shape = shape[:-nn]
    # Make sure the shape is the correct one ....
    for k, obj in enumerate(arraylist):
        nn = len(descr[k].shape)
        testshape = obj.shape[:len(obj.shape)-nn]
        if testshape != shape:
            raise ValueError, "Array-shape mismatch in array %d" % k
    # Reconstruct the descriptor, by creating a _data and _mask version
    return MultiTimeSeries(arraylist, dtype=descr)

def __getdates(dates=None, newdates=None, length=None, freq=None, 
               start_date=None):
    """Determines new dates (private function not meant to be used)."""
    if dates is None:
        if newdates is not None:
            if not hasattr(newdates, 'freq'):
                newdates = date_array(dlist=newdates, freq=freq)
        else:
            newdates = date_array(start_date=start_date, length=length,
                                  freq=freq)               
    elif not hasattr(dates, 'freq'):
        newdates = date_array(dlist=dates, freq=freq)
    else:
        newdates = dates
    return newdates

#..............................................................................
def fromrecords(reclist, dates=None, freq=None, start_date=None,
                dtype=None, shape=None, formats=None, names=None,
                titles=None, aligned=False, byteorder=None):
    """Creates a MaskedRecords from a list of records.

    The data in the same field can be heterogeneous, they will be promoted
    to the highest data type.  This method is intended for creating
    smaller record arrays.  If used to create large array without formats
    defined, it can be slow.

    If formats is None, then this will auto-detect formats. Use a list of
    tuples rather than a list of lists for faster processing.
    """    
    # reclist is in fact a mrecarray .................
    if isinstance(reclist, MultiTimeSeries):
        mdescr = reclist.dtype
        shape = reclist.shape
        return MultiTimeSeries(reclist, dtype=mdescr)
    # No format, no dtype: create from to arrays .....
    _data = mrecfromrecords(reclist, dtype=dtype, shape=shape, formats=formats, 
                            names=names, titles=titles, aligned=aligned, 
                            byteorder=byteorder)
    _dtype = _data.dtype
    # Check the names for a '_dates' .................
    newdates = None
    _names = list(_dtype.names)
    reserved = [n for n in _names if n.lower() in ['dates', '_dates']]
    if len(reserved) > 0:
        newdates = _data[reserved[-1]]
        [_names.remove(n) for n in reserved]
        _dtype = numeric.dtype([t for t in _dtype.descr \
                                    if t[0] not in reserved ])
        _data = [_data[n] for n in _names]
    #
    newdates = __getdates(dates=dates, newdates=newdates, length=len(_data), 
                          freq=freq, start_date=start_date)
    #
    return MultiTimeSeries(_data, dates=newdates, dtype=_dtype,
                           names=_names)


def fromtextfile(fname, delimitor=None, commentchar='#', missingchar='',
                 dates_column=None, varnames=None, vartypes=None,
                 dates=None):
    """Creates a multitimeseries from data stored in the file `filename`.

:Parameters:
    - `filename` : file name/handle
      Handle of an opened file.  
    - `delimitor` : Character *None*
      Alphanumeric character used to separate columns in the file.
      If None, any (group of) white spacestring(s) will be used.
    - `commentchar` : String *['#']*
      Alphanumeric character used to mark the start of a comment.
    - `missingchar` : String *['']*
      String indicating missing data, and used to create the masks.
    - `datescol` : Integer *[None]*
      Position of the columns storing dates. If None, a position will be 
      estimated from the variable names.
    - `varnames` : Sequence *[None]*
      Sequence of the variable names. If None, a list will be created from
      the first non empty line of the file.
    - `vartypes` : Sequence *[None]*
      Sequence of the variables dtypes. If None, the sequence will be estimated
      from the first non-commented line.  
    
    
    Ultra simple: the varnames are in the header, one line"""
    # Try to open the file ......................
    f = openfile(fname)
    # Get the first non-empty line as the varnames
    while True:
        line = f.readline()
        firstline = line[:line.find(commentchar)].strip()
        _varnames = firstline.split(delimitor)
        if len(_varnames) > 1:
            break
    if varnames is None:
        varnames = _varnames
    # Get the data ..............................
    _variables = MA.asarray([line.strip().split(delimitor) for line in f
                                  if line[0] != commentchar and len(line) > 1])
    (nvars, nfields) = _variables.shape
    # Check if we need to get the dates..........
    if dates_column is None:
        dates_column = [i for (i,n) in enumerate(list(varnames))
                            if n.lower() in ['_dates','dates']]
    elif isinstance(dates_column,(int,float)):
        if dates_column > nfields:
            raise ValueError,\
                  "Invalid column number: %i > %i" % (dates_column, nfields)
        dates_column = [dates_column,]
    if len(dates_column) > 0:
        cols = range(nfields)
        [cols.remove(i) for i in dates_column]
        newdates = date_array(_variables[:,dates_column[-1]])
        _variables = _variables[:,cols]    
        varnames = [varnames[i] for i in cols]
        if vartypes is not None:
            vartypes = [vartypes[i] for i in cols]
        nfields -= len(dates_column)
    else:
        newdates = None
    # Try to guess the dtype ....................
    if vartypes is None:
        vartypes = _guessvartypes(_variables[0])
    else:
        vartypes = [numeric.dtype(v) for v in vartypes]
        if len(vartypes) != nfields:
            msg = "Attempting to %i dtypes for %i fields!"
            msg += " Reverting to default."
            warnings.warn(msg % (len(vartypes), nfields))
            vartypes = _guessvartypes(_variables[0])
    # Construct the descriptor ..................
    mdescr = [(n,f) for (n,f) in zip(varnames, vartypes)]
    # Get the data and the mask .................
    # We just need a list of masked_arrays. It's easier to create it like that:
    _mask = (_variables.T == missingchar)
    _datalist = [masked_array(a,mask=m,dtype=t)
                     for (a,m,t) in zip(_variables.T, _mask, vartypes)]
    #
    newdates = __getdates(dates=dates, newdates=newdates, length=nvars, 
                          freq=None, start_date=None)
    return MultiTimeSeries(_datalist, dates=newdates, dtype=mdescr)
    


################################################################################

    
################################################################################
#if 1:
#    from tseries import aligned
#    import maskedarray.testutils
#    from maskedarray.testutils import *
#    
#    if 1:
#        mlist = ['2005-%02i' % i for i in range(1,13)]
#        mlist += ['2006-%02i' % i for i in range(1,13)]
#        mdata = numpy.arange(24)
#        mser1 = time_series(mdata, mlist, observed='SUMMED')
#        #
#    if 1:
#        mlist2 = ['2004-%02i' % i for i in range(1,13)]
#        mlist2 += ['2005-%02i' % i for i in range(1,13)]
#        mser2 = time_series(mdata, mlist2, observed='SUMMED')
#        #
#        (malg1,malg2) = aligned(mser1, mser2)
#    if 1:
#        mrec = MR.fromarrays([mser1])
#        descr = [('A',N.float_),('B',N.float_)]
#
#if 1:
#    import numpy as N
#    if 1:        
##    def setup(self):       
##        "Generic setup" 
#        d = N.arange(5)
#        m = MA.make_mask([1,0,0,1,1])
#        base_d = N.r_[d,d[::-1]].reshape(2,-1).T
#        base_m = N.r_[[m, m[::-1]]].T
#        base = MA.array(base_d, mask=base_m)    
#        mrec = MR.fromarrays(base.T,)
#        dlist = ['2007-%02i' % (i+1) for i in d]
#        dates = date_array(dlist)
#        ts = time_series(mrec,dates)
#        mts = MultiTimeSeries(mrec,dates)
#        
#        logmts = N.log(mts)
#        self_data = [d, m, mrec, dlist, dates, ts, mts]
#        #
#        mts.addfield(masked_array(d+10, mask=m[::-1]))
#        assert('f2' in mts.dtype.names)
#        assert_equal(mts.f2, d+10)
#        assert_equal(mts.f2._mask, m[::-1])

