import numpy
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
#---- --- Option conversion ---
#####---------------------------------------------------------------------------
obs_dict = {"UNDEFINED":None,
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
#
def fmtObserv(obStr):
    "Converts a possible 'Observed' string into acceptable values."
    if obStr is None:
        return None
    elif obStr.upper() in obs_dict.keys():
        return obStr.upper()    
    else:
        raise ValueError("Invalid value for observed attribute: %s " % str(obStr))


fmtfreq_dict = {'A': ['ANNUAL','ANNUALLY','YEAR','YEARLY'],
                'B': ['BUSINESS','BUSINESSLYT'],
                'D': ['DAY','DAILY',],
                'H': ['HOUR','HOURLY',],
                'M': ['MONTH','MONTHLY',],
                'Q': ['QUARTER','QUARTERLY',],
                'S': ['SECOND','SECONDLY',],
                'T': ['MINUTE','MINUTELY',],
                'W': ['WEEK','WEEKLY',],
                'U': ['UNDEF','UNDEFINED'],
                }
fmtfreq_revdict = reverse_dict(fmtfreq_dict)

def fmtFreq (freqStr):
    "Converts a possible 'frequency' string to acceptable values."
    if freqStr is None:
        return None    
    elif freqStr.upper() in fmtfreq_dict.keys():
        return freqStr[0].upper()
    elif freqStr.upper() in fmtfreq_revdict.keys():
        return fmtfreq_revdict[freqStr.upper()]
    else:
        raise ValueError("Invalid frequency: %s " % str(freqStr))
        
class DateSpec:
    "Fake data type for date variables."
    def __init__(self, freq):
        self.freq = fmtFreq(freq)
        
    def __hash__(self): 
        return hash(self.freq)
    
    def __eq__(self, other):
        if hasattr(other, "freq"): 
            return self.freq == other.freq
        else: 
            return False
    def __str__(self): 
        return "Date(%s)" % str(self.freq)
    
    

# define custom numpy types.
# Note: A more robust approach would register these as actual valid numpy types
# this is just a hack for now
numpy.dateA = DateSpec("Annual")
numpy.dateB = DateSpec("Business")
numpy.dateD = DateSpec("Daily")
numpy.dateH = DateSpec("Hourly")
numpy.dateM = DateSpec("Monthly")
numpy.dateQ = DateSpec("Quarterly")
numpy.dateS = DateSpec("Secondly")
numpy.dateT = DateSpec("Minutely")
numpy.dateW = DateSpec("Weekly")
numpy.dateU = DateSpec("Undefined")


freq_type_mapping = {'A': numpy.dateA,
                     'B': numpy.dateB,
                     'D': numpy.dateD,
                     'H': numpy.dateH,
                     'M': numpy.dateM,
                     'Q': numpy.dateQ,
                     'S': numpy.dateS,
                     'T': numpy.dateT,
                     'W': numpy.dateW,
                     'U': numpy.dateU,
                     }
        
def freqToType(freq):
    return freq_type_mapping[fmtFreq(freq)]

def isDateType(dtype):
    #TODO: That looks messy. We should simplify that
    if len([x for x in freq_type_mapping.values() if x == dtype]) > 0: 
        return True
    else: 
        return False

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
        

