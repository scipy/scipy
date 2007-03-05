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

from scipy.interpolate import fitpack

import maskedarray as MA

from cseries import TSER_CONSTANTS

"""add constants in cfame.FAME_CONSTANTS dictionary to global namespace
for this module"""

_g = globals()
_g.update(TSER_CONSTANTS)


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

_weekly_prefixes = ['W','WEEK','WEEKLY']
_week_end_map = {
    FR_WKSUN:'SUNDAY',
    FR_WKSAT:'SATURDAY',
    FR_WKFRI:'FRIDAY',
    FR_WKTHU:'THURSDAY',
    FR_WKWED:'WEDNESDAY',
    FR_WKTUE:'TUESDAY',
    FR_WKMON:'MONDAY'}

def _gen_weekly_strs(day):
    result = []
    for pr in _weekly_prefixes:
        result += [pr+'-'+day_str for day_str in (day[:3], day)]
    return result

freq_dict = { FR_ANN: ['A','Y','ANNUAL','ANNUALLY','YEAR','YEARLY'],
              FR_QTR: ['Q','QUARTER','QUARTERLY',],
              FR_MTH: ['M','MONTH','MONTHLY',],
              FR_BUS: ['B','BUSINESS','BUSINESSLY'],
              FR_DAY: ['D','DAY','DAILY',],
              FR_HR: ['H','HOUR','HOURLY',],
              FR_MIN: ['T','MINUTE','MINUTELY',],
              FR_SEC: ['S','SECOND','SECONDLY',],
              FR_UND: ['U','UNDEF','UNDEFINED'],
                }
                
for _freq, day_str in _week_end_map.iteritems():
    freq_dict[_freq] = _gen_weekly_strs(day_str)
freq_dict[FR_WK] += _weekly_prefixes
    
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
        
