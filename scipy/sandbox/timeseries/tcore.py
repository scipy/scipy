"""
A collection of tools for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy
import numpy.core.numeric as numeric

import maskedarray as MA

#####---------------------------------------------------------------------------
#---- --- Generic functions ---
#####---------------------------------------------------------------------------
def first_unmasked_val(a):
    "Returns the first unmasked value in a 1d maskedarray."
    (i,j) = MA.extras.flatnotmasked_edges(a)
    return a[i]

def last_unmasked_val(a):
    "Returns the last unmasked value in a 1d maskedarray."
    (i,j) = MA.extras.flatnotmasked_edges(a)
    return a[j]

def reverse_dict(d):
    "Reverses the keys and values of a dictionary."
    alt = []
    tmp = [alt.extend([(w,k) for w in v]) for (k,v) in d.iteritems()]
    return dict(alt)



#####---------------------------------------------------------------------------
#---- --- Misc functions ---
#####---------------------------------------------------------------------------
#http://aspn.activestate.com/ASPN/Mail/Message/python-tutor/2302348
def flatten_sequence(iterable):
    """Flattens a compound of nested iterables."""
    itm = iter(iterable)
    for elm in itm:
        if hasattr(elm,'__iter__') and not isinstance(elm, basestring):
            for f in flatten_sequence(elm):
                yield f
        else:
            yield elm

def flatargs(*args):
    "Flattens the arguments."
    if not hasattr(args, '__iter__'):
        return args
    else:
        return flatten_sequence(args)
