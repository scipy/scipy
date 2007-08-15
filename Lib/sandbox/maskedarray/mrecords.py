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

import numpy
from numpy import bool_, complex_, float_, int_, str_, object_
import numpy.core.numeric as numeric
import numpy.core.numerictypes as ntypes
from numpy.core.defchararray import chararray
from numpy.core.records import find_duplicate

from numpy.core.records import format_parser, record, recarray
from numpy.core.records import fromarrays as recfromarrays

ndarray = numeric.ndarray
_byteorderconv = numpy.core.records._byteorderconv
_typestr = ntypes._typestr

import maskedarray as MA
from maskedarray import masked, nomask, mask_or, filled, getmask, getmaskarray, \
    masked_array, make_mask
from maskedarray import MaskedArray
from maskedarray.core import default_fill_value, masked_print_option

import warnings

reserved_fields = ['_data','_mask','_fieldmask', 'dtype']

def _getformats(data):
    """Returns the formats of each array of arraylist as a comma-separated 
    string."""
    if hasattr(data,'dtype'):
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

def _checknames(descr, names=None):
    """Checks that the field names of the descriptor `descr` are not some 
    reserved keywords. If this is the case, a default 'f%i' is substituted.
    If the argument `names` is not None, updates the field names to valid names.    
    """    
    ndescr = len(descr)
    default_names = ['f%i' % i for i in range(ndescr)]
    if names is None:
        new_names = default_names
    else:
        if isinstance(names, (tuple, list)):
            new_names = names
        elif isinstance(names, str):
            new_names = names.split(',')
        else:
            raise NameError, "illegal input names %s" % `names`
        nnames = len(new_names)
        if nnames < ndescr:
            new_names += default_names[nnames:]
    ndescr = []
    for (n, d, t) in zip(new_names, default_names, descr.descr):
        if n in reserved_fields:
            if t[0] in reserved_fields: 
                ndescr.append((d,t[1]))
            else:
                ndescr.append(t)
        else:
            ndescr.append((n,t[1]))
    return numeric.dtype(ndescr)
    


class MaskedRecords(MaskedArray, object):
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
    def __new__(cls, data, mask=nomask, dtype=None, 
                hard_mask=False, fill_value=None,
#                offset=0, strides=None,
                formats=None, names=None, titles=None,
                byteorder=None, aligned=False):
        # Get the new descriptor ................
        if dtype is not None:
            descr = numeric.dtype(dtype)
        else:
            if formats is None:
                formats = _getformats(data)
            parsed = format_parser(formats, names, titles, aligned, byteorder)
            descr = parsed._descr
        if names is not None:
            descr = _checknames(descr,names)
        _names = descr.names    
        mdescr = [(n,'|b1') for n in _names]
        # get the shape .........................
   