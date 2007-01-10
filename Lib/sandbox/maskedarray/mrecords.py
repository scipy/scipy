"""mrecords
Defines a class of record arrays supporting masked arrays.

:author: Pierre Gerard-Marchant
:contact: pierregm_at_uga_dot_edu
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import sys
import types

import numpy as N
from numpy import bool_, complex_, float_, int_, str_, object_
import numpy.core.numeric as numeric
import numpy.core.numerictypes as ntypes
from numpy.core.defchararray import chararray
from numpy.core.records import find_duplicate
from numpy import bool_

from numpy.core.records import format_parser, record, recarray

ndarray = numeric.ndarray
_byteorderconv = N.core.records._byteorderconv
_typestr = ntypes._typestr

import maskedarray as MA
reload(MA)
filled = MA.filled
getmaskarray = MA.getmaskarray
masked = MA.masked
nomask = MA.nomask
mask_or = MA.mask_or
masked_array = MA.masked_array

import warnings
import logging
logging.basicConfig(level=logging.DEBUG,
                    format='%(name)-15s %(levelname)s %(message)s',)


class mrecarray(ndarray):
    """
    
:IVariables:
    - `__localfdict` : Dictionary
        Dictionary of local fields (`f0_data`, `f0_mask`...)
    - `__globalfdict` : Dictionary
        Dictionary of global fields, as the combination of a `_data` and a `_mask`.
        (`f0`)
    """
    __localfdict = {}
    __globalfdict = {}
    def __new__(subtype, shape, dtype=None, buf=None, offset=0, strides=None,
                formats=None, names=None, titles=None,
                byteorder=None, aligned=False, hard_mask=False):
        
        if dtype is not None:
            descr = numeric.dtype(dtype)
        else:
            descr = format_parser(formats, names, titles, aligned, byteorder)._descr

        if buf is None:
            mrec = ndarray.__new__(subtype, shape, (record, descr))
        else:
            mrec = ndarray.__new__(subtype, shape, (record, descr),
                                   buffer=buf, offset=offset,
                                   strides=strides)           
        # Stores the field names in directories..
        mrec.__localfdict = dict(descr.fields)
        mrec.__globalfdict = {}
        keys = sorted(mrec.__localfdict.keys())
        for i in range(len(keys)//2):
            ikey = keys[2*i]
            nkey = "_".join(ikey.split('_')[:-1])
            (dfield, mfield) = ("%s_data" % nkey, "%s_mask" % nkey)
            mrec.__globalfdict[nkey] = dict(_data=mrec.__localfdict[dfield],
                                            _mask=mrec.__localfdict[mfield])
        mrec._hardmask = hard_mask
        return mrec
    
    def __getallfields(self,fieldname):
        """Returns all the fields sharing the same fieldname base.
    The fieldname base is either `_data` or `_mask`."""
        logging.debug('__getallfields(%s)' % fieldname)
        (names, formats, offsets, objs) = ([], [], [], [])
        fkeyname = '%%s%s' % fieldname
        for s in self._mrecarray__globalfdict.keys():
            fkey = fkeyname % s
            fattr =  self._mrecarray__localfdict[fkey]
            names.append(s)
            obj = self.__getobj(ndarray.__getitem__(self,fkey))
            objs.append(obj)
            formats.append(fattr[0])
            offsets.append(fattr[1])
        descr = [(n,f) for (n,f) in zip(names, formats)]
        return N.core.records.fromarrays(objs, dtype=descr)
    
    def _getdata(self):
        """Returns all the `_data` fields."""
        return self.__getallfields('_data')
    
    def _getfieldmask(self):
        """Returns a recarray of the mask of each field."""
        return self.__getallfields('_mask')
    
    def _getmask(self):
        """Returns the mask of the mrecarray.
    An element of the mrecarray is considered masked when all the corresponding
    fields are masked."""
        nbfields = len(self.dtype )//2
        recmask = self._getfieldmask().view(bool_).reshape(-1,nbfields)
        return recmask.all(1)        
    
    def __getobj(self, obj, viewtype=ndarray):
        "Returns an object as a view of a ndarray, or as itself."
        if (isinstance(obj, ndarray) and obj.dtype.isbuiltin):
            return obj.view(viewtype)
        return obj
    #......................................................
    def __getattribute__(self, attr):
        try:
            # Returns a generic attribute
            return object.__getattribute__(self,attr)
        except AttributeError: 
            # OK, so attr must be a field name
            pass
        # Get the list of fields ......
        fdict = ndarray.__getattribute__(self,'_mrecarray__localfdict') or {}
        # Case #1: attr is a basic field
        if attr in fdict.keys():
            fattr = fdict[attr]
            obj = self.getfield(*fattr)
            if obj.dtype.fields:
                return obj
            if obj.dtype.char in 'SU':
                return obj.view(chararray)
            return obj.view(ndarray)
        # Case #2: attr is acompund field
        elif ("%s_data" % attr) in fdict.keys():
            data = self.getfield(*fdict["%s_data" % attr ][:2])
            mask = self.getfield(*fdict["%s_mask" % attr ][:2])    
            return MA.masked_array(data.view(ndarray), 
                              mask=mask.view(ndarray), 
                              copy=False)
        # Case #3/4/5: attr is a generic field
        elif attr == '_data':    
            func = ndarray.__getattribute__(self,'_mrecarray__getallfields')
            return func.__call__('_data')        
        elif attr == '_fieldmask':    
            func = ndarray.__getattribute__(self,'_mrecarray__getallfields')
            return func.__call__('_mask')            
        elif attr == '_mask':
            logging.debug('__getattribute__: all fields %s' % attr)
            func = ndarray.__getattribute__(self,'_getmask')
            return func.__call__()            
        # Case #6: attr is not a field at all !
        else:
            raise AttributeError, "record array has no attribute %s" % attr

# Save the dictionary
#  If the attr is a field name and not in the saved dictionary
#  Undo any "setting" of the attribute and do a setfield
# Thus, you can't create attributes on-the-fly that are field names. 

    def __setattr__(self, attr, val):
        # gets some status on attr: an existing field ? a new attribute ?
        fdict = ndarray.__getattribute__(self,'_mrecarray__localfdict') or {}
        gdict = ndarray.__getattribute__(self,'_mrecarray__globalfdict') or {}
        attrlist = fdict.keys() + ['_data', '_fieldmask', '_mask']
        isvalidattr = (attr in attrlist ) or ('%s_data' % attr in attrlist)
        newattr = attr not in self.__dict__
        
        try:
            # Is attr a generic attribute ?
            ret = object.__setattr__(self, attr, val)
        except:
            # Not a generic attribute: exit if it's not a valid field
            if not isvalidattr:
                exctype, value = sys.exc_info()[:2]
                raise exctype, value
        else:
            if not isvalidattr:
                return ret
            if newattr:         # We just added this one
                try:            #  or this setattr worked on an internal
                                #  attribute. 
                    object.__delattr__(self, attr)
                except:
                    return ret
        
        # Case #1.: Basic field ............
        if attr in fdict.keys():
            return self.setfield(val, *fdict[attr][:2])
        # Case #2 Compund field ............
        elif ("%s_data" % attr) in fdict.keys():
            data = self.setfield(filled(val), *fdict["%s_data" % attr ][:2])
            mask = self.setfield(getmaskarray(val), *fdict["%s_mask" % attr ][:2])             
            return 
        elif attr == '_data':    
            fval = filled(val)
            for k in gdict.keys():
                self.setfield(fval, *gdict["%s_data" % k ][:2])
            return
#            func = ndarray.__getattribute__(self,'_mrecarray__getallfields')
#            return func.__call__('_data')        
        elif attr == '_fieldmask':    
            mval = getmaskarray(val)
            for k in gdict.keys():
                self.setfield(mval, *gdict["%s_mask" % k ][:2])
            return
#            func = ndarray.__getattribute__(self,'_mrecarray__getallfields')
#            return func.__call__('_mask')            
        elif attr == '_mask':
            logging.debug(" setattr _mask to %s [%s]" % (val,(val is nomask)))
            if self._hardmask:
                logging.debug("setattr: object has hardmask")
                if val is not nomask:
                    mval = getmaskarray(val)
                    for k in gdict.keys():
                        fkey = fdict["%s_mask" % k ][:2]
                        m = mask_or(mval, self.getfield(*fkey))
                        logging.debug("setattr: set %s to %s" % (k,m))
                        self.setfield(m, *fkey)
            else:
                mval = getmaskarray(val)
                for k in gdict.keys():
                    self.setfield(mval, *fdict["%s_mask" % k ][:2])            
            logging.debug('__getattribute__: all fields %s' % attr)
            return
#            func = ndarray.__getattribute__(self,'_getmask')
#            return func.__call__()            
            

    #......................................................        
    def __getitem__(self, indx):
        logging.debug('__getitem__ got %s' % indx)
        try:
            obj = ndarray.__getitem__(self, indx)
        except ValueError:
            if indx in self.__globalfdict.keys():
                objd = ndarray.__getitem__(self, "%s_data" % indx)
                objm = ndarray.__getitem__(self, "%s_mask" % indx)
                return MA.masked_array(objd.view(ndarray), 
                                  mask=objm.view(ndarray))
            elif indx in ['_data', '_fieldmask']:
                return self.__getallfields(indx)
            elif indx == '_mask':
                return self._getmask()
            else:
                msg = "Cannot do anything w/ indx '%s'!" % indx
                raise ValueError, msg
            
#        logging.debug('__getitem__ send %s' % type(self.__getobj(obj)))
        return self.__getobj(obj)
    #.......................................................
    def field(self,attr, val=None):
        """Sets the field `attr` to the new value `val`.
    If `val` is None, returns the corresponding field.
    """
        if isinstance(attr,int):
            names = ndarray.__getattribute__(self,'dtype').names
            attr = names[attr]
        
        fdict = ndarray.__getattribute__(self,'_mrecarray__localfdict') or {}
        f = fdict[attr]
        # Case #1: just retrieve the data .......
        if val is None:
            try:
                return self.__getattribute__(attr)
            except:
                raise ValueError, "Unable to retrieve field '%s'" % attr 
        # Case #2: set the field to a new value ..
        else:
            try:
                return self.__setattribute__(attr)
            except:
                raise ValueError, "Unable to set field '%s'" % attr 
    #......................................................
    def view(self, obj):
        """Returns a view of the mrecarray."""
        try:
            if issubclass(obj, ndarray):
                logging.debug('direct view as %s' % obj)
                return ndarray.view(self, obj)
        except TypeError:
            pass
        dtype = numeric.dtype(obj)
        if dtype.fields is None:
            return self.__array__().view(dtype)
        return ndarray.view(self, obj)            
    #............................................
    def harden_mask(self):
        "Forces the mask to hard"
        self._hardmask = True
    def soften_mask(self):
        "Forces the mask to soft"
        self._hardmask = False

#####---------------------------------------------------------------------------
#---- --- Constructors ---
#####---------------------------------------------------------------------------
def _splitfields(descr):
    """Creates a new descriptor from the descriptor `descr`.
    The initial fields are renamed by adding a `_data` suffix to the name.
    Their dtype is kept.
    New fields are also created from the initial ones by adding a `_mask` suffix
    to the name.
    The dtype of these latter is set to `bool_`
    """
    mdescr = []
    for (n,d) in descr.descr:
        mdescr.append( ("%s_data" % n, d) )
        mdescr.append( ("%s_mask" % n, bool_) )
    return numeric.dtype(mdescr)

def fromarrays(arraylist, dtype=None, shape=None, formats=None,
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
        # go through each object in the list to see if it is an ndarray
        # and determine the formats.
        formats = ''
        for obj in arraylist:
            if not isinstance(obj, ndarray):
                raise ValueError, "item in the array list must be an ndarray."
            formats += _typestr[obj.dtype.type]
            if issubclass(obj.dtype.type, ntypes.flexible):
                formats += `obj.itemsize`
            formats += ','
        formats = formats[:-1]
        
    logging.debug("fromarrays: formats",formats)
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
    mdescr = _splitfields(descr)
    _array = mrecarray(shape, mdescr)
    _names = mdescr.names
    # Populate the record array (makes a copy)
    for i in range(len(arraylist)):
        logging.debug("fromarrays: i:%i-%s/%s" % \
                      (i, arraylist[i]._data, MA.getmaskarray(arraylist[i])))
        logging.debug("fromarrays: i:%i-%s/%s" % \
                      (i,_names[2*i], _names[2*i+1]))
        _array[_names[2*i]] = arraylist[i]._data
        _array[_names[2*i+1]] = getmaskarray(arraylist[i])
    return _array
#..............................................................................
def fromrecords(reclist, dtype=None, shape=None, formats=None, names=None,
                titles=None, aligned=False, byteorder=None):
    """Creates a mrecarray from a list of records.

    The data in the same field can be heterogeneous, they will be promoted
    to the highest data type.  This method is intended for creating
    smaller record arrays.  If used to create large array without formats
    defined

        r=fromrecords([(2,3.,'abc')]*100000)

        it can be slow.

        If formats is None, then this will auto-detect formats. Use list of
        tuples rather than list of lists for faster processing.

    >>> r=fromrecords([(456,'dbe',1.2),(2,'de',1.3)],names='col1,col2,col3')
    >>> print r[0]
    (456, 'dbe', 1.2)
    >>> r.col1
    array([456,   2])
    >>> r.col2
    chararray(['dbe', 'de'])
    >>> import cPickle
    >>> print cPickle.loads(cPickle.dumps(r))
    recarray[
    (456, 'dbe', 1.2),
    (2, 'de', 1.3)
    ]
    """    
    # Case #1: reclist is in fact a mrecarray ........
    if isinstance(reclist, mrecarray):
        mdescr = reclist.dtype
        shape = reclist.shape
        _array = mrecarray(shape, mdescr)
        for (i,r) in enumerate(reclist):
            _array[i] = r
        return _array
    
    # No format, no dtype: create from to arrays .....
    nfields = len(reclist[0])
    if formats is None and dtype is None:  # slower
        if isinstance(reclist, recarray):
            arrlist = [reclist.field(i) for i in range(len(reclist.dtype))]
            if names is None:
                names = reclist.dtype.names
        else:
            obj = numeric.array(reclist,dtype=object)
            arrlist = [numeric.array(obj[...,i].tolist()) 
                               for i in xrange(nfields)]
        return fromarrays(arrlist, formats=formats, shape=shape, names=names,
                          titles=titles, aligned=aligned, byteorder=byteorder)
    # Construct the descriptor .......................
    if dtype is not None:
        descr = numeric.dtype(dtype)
        _names = descr.names
    else:
        parsed = format_parser(formats, names, titles, aligned, byteorder)
        _names = parsed._names
        descr = parsed._descr
    mdescr = _splitfields(descr)

    try:
        retval = numeric.array(reclist, dtype = descr)
    except TypeError:  # list of lists instead of list of tuples
        if (shape is None or shape == 0):
            shape = len(reclist)*2
        if isinstance(shape, (int, long)):
            shape = (shape*2,)
        if len(shape) > 1:
            raise ValueError, "Can only deal with 1-d array."
        _array = recarray(shape, mdescr)
        raise NotImplementedError,"I should really test that..."
        for k in xrange(_array.size):
            _array[k] = tuple(reclist[k])
        return _array
    else:
        if shape is not None and retval.shape != shape:
            retval.shape = shape
    #
    tmp = retval.view(recarray)
    _array = mrecarray(shape, mdescr)
    for n in tmp.dtype.names:
        _array['%s_data' % n] = tmp[n]
        _array['%s_mask' % n] = nomask
    return _array

def _guessvartypes(arr):        
    """Tries to guess the dtypes of the str_ ndarray `arr`, by testing element-wise
    conversion. Returns a list of dtypes.
    The array is first converted to ndarray. If the array is 2D, the test is 
    performed on the first line. An exception is raised if the file is 3D or more.
    """
    vartypes = []
    arr = numeric.asarray(arr)
    if len(arr.shape) == 2 :
        arr = arr[0]
    elif len(arr.shape) > 2:
        raise ValueError, "The array should be 2D at most!"
    # Start the conversion loop .......
    for f in arr:
        try:
            val = int(f)
        except ValueError:
            try:
                val = float(f)
            except ValueError:
                try: 
                    val = complex(f)
                except ValueError:
                    print "str_!"
                    vartypes.append(arr.dtype)
                else:
                    vartypes.append(complex_)
            else:
                vartypes.append(float_)
        else:
            vartypes.append(int_)
    return vartypes

def openfile(fname):
    "Opens the file handle of file `fname`"
    # A file handle ...................
    if hasattr(fname, 'readline'):
        return fname
    # Try to open the file and guess its type
    try:
        f = open(fname)
    except IOError:
        raise IOError, "No such file: '%s'" % fname
    if f.readline()[:2] != "\\x":
        f.seek(0,0)
        return f
    raise NotImplementedError, "Wow, binary file" 
    

def fromtextfile(fname, delimitor=None, commentchar='#', missingchar='',
                 varnames=None, vartypes=None):
    """Creates a mrecarray from data stored in the file `filename`.

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
        logging.debug("_VARNAMES:%s-%s"% (_varnames,len(_varnames)))
        if len(_varnames) > 1:
            break
    if varnames is None:
        varnames = _varnames
    # Get the data ..............................
    _variables = [line.strip().split(delimitor) for line in f
                   if line[0] != commentchar and len(line) > 1]
    _variables = N.array(_variables)
    #_variables = MA.masked_equal(_variables,'')
    (nvars, nfields) = _variables.shape
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
    return fromarrays(_datalist, dtype=mdescr)
    


################################################################################
from maskedarray.testutils import assert_equal, assert_array_equal
if 1:
    if 0:
        fcontent = """#
'One (S)','Two (I)','Three (F)','Four (M)','Five (-)','Six (C)'
'strings',1,1.0,'mixed column',,1
'with embedded "double quotes"',2,2.0,1.0,,1
'strings',3,3.0E5,3,,1
'strings',4,-1e-10,,,1
"""    
        import os
        from datetime import datetime
        fname = 'tmp%s' % datetime.now().strftime("%y%m%d%H%M%S%s")
        f = open(fname, 'w')
        f.write(fcontent)
        f.close()
        mrectxt = fromtextfile(fname,delimitor=',',varnames='ABCDEFG')        
        os.unlink(fname)
        #
        assert(isinstance(mrectxt, mrecarray))
        assert_equal(mrectxt.F, [1,1,1,1])
        assert_equal(mrectxt.E._mask, [1,1,1,1])
        assert_equal(mrectxt.C, [1,2,3.e+5,-1e-10])
#...............................................................................


        
    
        

