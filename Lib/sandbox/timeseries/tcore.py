"""
A collection of tools for timeseries

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'

import numpy
import numpy.core.numeric as numeric

from scipy.interpolate import fitpack

import maskedarray as MA
from maskedarray import masked, nomask, getmask

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
#---- --- Option conversion ---
#####---------------------------------------------------------------------------
fmtobs_dict = {'UNDEFINED': ['UNDEF','UNDEFINED',None],
               'BEGINNING': ['BEGIN','BEGINNING'],
               'ENDING': ['END','ENDING'],
               'AVERAGED': ['AVERAGE','AVERAGE','MEAN'],
               'SUMMED': ['SUM','SUMMED'],
               'MAXIMUM': ['MAX','MAXIMUM','HIGH'],
               'MINIMUM': ['MIN','MINIMUM','LOW']}

obs_dict = {None:None,
            "UNDEFINED":None,
            "UNDEF":None,
            "BEGIN": first_unmasked_val,
            "BEGINNING": first_unmasked_val,
            "END": last_unmasked_val,
            "ENDING": last_unmasked_val,
            "AVERAGED": MA.average,
            "AVERAGE": MA.average,
            "MEAN": MA.average,
            "SUMMED": MA.sum,
            "SUM": MA.sum,
            "MAXIMUM": MA.maximum,
            "MAX": MA.maximum,
            "MINIMUM": MA.minimum,
            "MIN": MA.minimum,
            }
obsDict = obs_dict
fmtobs_revdict = reverse_dict(fmtobs_dict)

#
def fmtObserv(obStr):
    "Converts a possible 'Observed' string into acceptable values."
    if obStr is None:
        return fmtobs_revdict[None]
    elif obStr.upper() in fmtobs_revdict:
        return fmtobs_revdict[obStr.upper()]
    else:
        raise ValueError("Invalid value for observed attribute: %s " % str(obStr))



freq_dict = { 1000: ['A','Y','ANNUAL','ANNUALLY','YEAR','YEARLY'],
              2000: ['Q','QUARTER','QUARTERLY',],
              3000: ['M','MONTH','MONTHLY',],
              4000: ['W','WEEK','WEEKLY',],
              5000: ['B','BUSINESS','BUSINESSLY'],
              6000: ['D','DAY','DAILY',],
              7000: ['H','HOUR','HOURLY',],
              8000: ['T','MINUTE','MINUTELY',],
              9000: ['S','SECOND','SECONDLY',],
             -9999: ['U','UNDEF','UNDEFINED'],
                }
freq_revdict = reverse_dict(freq_dict)

def freq_fromstr(freq_asstr):
    "Converts a frequency given as string to the corresponding integer."
    freq_asstr = freq_asstr.upper()
    if freq_asstr not in freq_revdict.keys():
        raise ValueError, "Invalid frequency string %s" % freq_asstr
    return freq_revdict[freq_asstr]
    
def freq_tostr(freq_asint):
    "Converts a frequency given as integer to the corresponding symbol."
    if freq_asint not in freq_dict.keys():
        raise ValueError, "Invalid frequency representation %s" % freq_asint
    return freq_dict[freq_asint][0]

def check_freq(freq):
    "Converts a possible 'frequency' string to acceptable values."
    if freq is None:
        return None
    elif isinstance(freq, int):
        if freq not in freq_dict.keys():
            raise ValueError("Invalid frequency: %s " % str(freq))
        return freq
    elif freq.upper() in freq_revdict.keys():
        return freq_revdict[freq.upper()]
    else:
        raise ValueError("Invalid frequency: %s " % str(freq))
    
def check_freqstr(freq):
    if freq is None:
        return None
    elif isinstance(freq, int):
        if freq not in freq_dict.keys():
            raise ValueError("Invalid frequency: %s " % str(freq))
        return freq_dict[freq][0]
    elif freq.upper() in freq_revdict.keys():
        return freq_dict[freq_revdict[freq.upper()]][0]
    else:
        raise ValueError("Invalid frequency: %s " % str(freq))    
fmtFreq = check_freqstr
        
#class DateSpec:
#    "Fake data type for date variables."
#    def __init__(self, freq):
#        self.freq = fmtFreq(freq)
#        
#    def __hash__(self): 
#        return hash(self.freq)
#    
#    def __eq__(self, other):
#        if hasattr(other, "freq"): 
#            return self.freq == other.freq
#        else: 
#            return False
#    def __str__(self): 
#        return "Date(%s)" % str(self.freq) 
#
## define custom numpy types.
## Note: A more robust approach would register these as actual valid numpy types
## this is just a hack for now
#numpy.dateA = DateSpec("Annual")
#numpy.dateB = DateSpec("Business")
#numpy.dateD = DateSpec("Daily")
#numpy.dateH = DateSpec("Hourly")
#numpy.dateM = DateSpec("Monthly")
#numpy.dateQ = DateSpec("Quarterly")
#numpy.dateS = DateSpec("Secondly")
#numpy.dateT = DateSpec("Minutely")
#numpy.dateW = DateSpec("Weekly")
#numpy.dateU = DateSpec("Undefined")
#
#
#freq_type_mapping = {'A': numpy.dateA,
#                     'B': numpy.dateB,
#                     'D': numpy.dateD,
#                     'H': numpy.dateH,
#                     'M': numpy.dateM,
#                     'Q': numpy.dateQ,
#                     'S': numpy.dateS,
#                     'T': numpy.dateT,
#                     'W': numpy.dateW,
#                     'U': numpy.dateU,
#                     }
#        
#def freqToType(freq):
#    "Returns the Date dtype corresponding to the given frequency."
#    return freq_type_mapping[fmtFreq(freq)]
#
#def isDateType(dtype):
#    "Returns True whether the argument is the dtype of a Date."
#    #TODO: That looks messy. We should simplify that
#    if len([x for x in freq_type_mapping.values() if x == dtype]) > 0: 
#        return True
#    else: 
#        return False

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
        

#####---------------------------------------------------------------------------
#---- --- Functions for filling in masked values in a masked array ---
#####---------------------------------------------------------------------------

def forward_fill(marr, maxgap=None):
    """forward_fill(marr, maxgap=None)

Forward fill masked values in a 1-d array when there are <= maxgap
consecutive masked values. If maxgap is None, then forward fill all
masked values."""

    if numeric.ndim(marr) > 1:
        raise ValueError,"The input array should be 1D only!"

    marr = MA.array(marr, copy=True)
    if getmask(marr) is nomask or marr.size == 0:
        return marr

    currGap = 0

    if maxgap is not None:
        for i in range(1, marr.size):
            if marr._mask[i]:
                currGap += 1
                if currGap <= maxgap and not marr._mask[i-1]:
                    marr._data[i] = marr._data[i-1]
                    marr._mask[i] = False
                elif currGap == maxgap + 1:
                    marr._mask[i-maxgap:i] = True
            else:
                currGap = 0               
    else:
        for i in range(1, marr.size):
            if marr._mask[i] and not marr._mask[i-1]:
                marr._data[i] = marr._data[i-1]
                marr._mask[i] = False
    return marr


def backward_fill(marr, maxgap=None):
    """backward_fill(marr, maxgap=None)

backward fill masked values in a 1-d array when there are <= maxgap
consecutive masked values. If maxgap is None, then backward fill all
masked values."""
    return forward_fill(marr[::-1], maxgap=maxgap)[::-1]
    

def interp_masked1d(marr, kind='linear'):
    """interp_masked1d(marr, king='linear')

Interpolate masked values in marr according to method kind.
kind must be one of 'constant', 'linear', 'cubic', quintic'
"""
    if numeric.ndim(marr) > 1: 
        raise ValueError("array must be 1 dimensional!")
    #
    marr = MA.array(marr, copy=True)
    if getmask(marr) is nomask: 
        return marr
    #
    unmaskedIndices = (~marr._mask).nonzero()[0]
    if unmaskedIndices.size < 2: 
        return marr
    #    
    kind = kind.lower()
    if kind == 'constant': 
        return forward_fill(marr)
    try:
        k = {'linear' : 1,
             'cubic' : 3,
             'quintic' : 5}[kind.lower()]
    except KeyError:
        raise ValueError("Unsupported interpolation type.")
    
    first_unmasked, last_unmasked = MA.extras.flatnotmasked_edges(marr)
    
    vals = marr.data[unmaskedIndices]
    
    tck = fitpack.splrep(unmaskedIndices, vals, k=k)
    
    maskedIndices = marr._mask.nonzero()[0]
    interpIndices = maskedIndices[(maskedIndices > first_unmasked) & \
                                  (maskedIndices < last_unmasked)]
    marr[interpIndices] = fitpack.splev(interpIndices, tck).astype(marr.dtype)
    return marr
