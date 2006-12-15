import numpy
from numpy import ma


#############################################################
############## generally applicable functions ###############
#############################################################
def apply_along_axis(func1d, axis, arr, *args):
    """ Execute func1d(arr[i],*args) where func1d takes 1-D arrays
        and arr is an N-d array.  i varies so as to apply the function
        along the given axis for each 1-d subarray in arr.
        
        Slightly modified version of the standard numpy version to work with masked arrays.
    """

    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis,nd))
    ind = [0]*(nd-1)
    i = numpy.zeros(nd,'O')
    indlist = range(nd)
    indlist.remove(axis)
    i[axis] = slice(None,None)
    outshape = numpy.asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    res = func1d(arr[tuple(i.tolist())],*args)
    #  if res is a number, then we have a smaller output array
    if not hasattr(res,'shape') or len(res.shape) == 0:
        outarr = ma.zeros(outshape,ma.asarray(res).dtype)
        outarr[ind] = res
        Ntot = numpy.product(outshape)
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= outshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist,ind)
            res = func1d(arr[tuple(i.tolist())],*args)
            outarr[ind] = res
            k += 1
        return outarr
    else:
        Ntot = numpy.product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = len(res)
        outarr = ma.zeros(outshape,ma.asarray(res).dtype)
        outarr[tuple(i.tolist())] = res
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= holdshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())],*args)
            outarr[tuple(i.tolist())] = res
            k += 1
        return outarr


def first_unmasked(m):
    return __unmasked(m, False, 0)
    
def last_unmasked(m):
    return __unmasked(m, False, -1)

def first_unmasked_val(m):
    return __unmasked(m, True, 0)
    
def last_unmasked_val(m):
    return __unmasked(m, True, -1)


def __unmasked(m, get_val, relpos):
    idx = numpy.where(m.mask == False)
    if len(idx) != 0 and len(idx[0]) != 0:
        idx = idx[0][relpos]
    else:
        idx = None
        
    if get_val:
        if idx is None: return ma.masked
        else: return m[idx]
    else:
        return idx
#############################################################

#converts possible strings for frequency into acceptable values             
def fmtFreq (freqStr):
    if freqStr is None:
        return None    
    elif freqStr.upper() in ("A","ANNUAL","B","BUSINESS","D","DAILY","M","MONTHLY","Q","QUARTERLY","S","SECONDLY"):
        return freqStr[0].upper()
    else:
        raise ValueError("Invalid frequency: "+str(freqStr))
        

obsDict = {
            "UNDEFINED":None,
            "BEGINNING":first_unmasked_val,
            "END":last_unmasked_val,
            "AVERAGED":ma.average,
            "SUMMED":ma.sum,
            "MAXIMUM":ma.maximum,
            "MINIMUM":ma.minimum
          }

#converts possible strings for observed into acceptable values
def fmtObserv(obStr):

    obsVals = list(obsDict)

    if obStr is None:
        return None
    elif obStr.upper() in obsVals:
        return obStr.upper()    
    elif obStr.upper() in ("UNDEFINED", "BEGIN", "END", "AVERAGE", "SUM", "MAX", "MIN"):
        obStr = obStr.upper()
        for x in obsVals:
            if obStr[:2] == x[:2]:
                return x       
    else:
        raise ValueError("Invalid value for observed attribute: "+str(obStr))
        
def freqToType(freq):
    return freqTypeMapping[fmtFreq(freq)]


# fake data type for date variables
class DateSpec:
    def __init__(self, freq):
        self.freq = fmtFreq(freq)
        
    def __hash__(self): return hash(self.freq)
    
    def __eq__(self, other):
        if hasattr(other, "freq"): return self.freq == other.freq
        else: return False
        
    def __str__(self): return "Date(" + str(self.freq) + ")"
    
    

# define custom numpy types.
# Note: A more robust approach would register these as actual valid numpy types
# this is just a hack for now
numpy.dateS = DateSpec("Secondly")
numpy.dateD = DateSpec("Daily")
numpy.dateB = DateSpec("Business")
numpy.dateM = DateSpec("Monthly")
numpy.dateQ = DateSpec("Quarterly")
numpy.dateA = DateSpec("Annual")

freqTypeMapping = {
    'S':numpy.dateS,
    'D':numpy.dateD,
    'B':numpy.dateB,
    'M':numpy.dateM,
    'Q':numpy.dateQ,
    'A':numpy.dateA
}

def isDateType(dtype):
    if len([x for x in (numpy.dateS,numpy.dateD,numpy.dateB,numpy.dateM,numpy.dateQ,numpy.dateA) if x == dtype]) > 0: return True
    else: return False

