import numpy

#converts possible strings for frequency into acceptable values             
def fmtFreq (freqStr):
    if freqStr is None:
        return None    
    elif freqStr.upper() in ("A","ANNUAL","B","BUSINESS","D","DAILY","M","MONTHLY","Q","QUARTERLY","S","SECONDLY"):
        return freqStr[0].upper()
    else:
        raise ValueError("Invalid frequency: "+str(freqStr))
        

#converts possible strings for observed into acceptable values
def fmtObserv(obStr):

    obsVals = (   "UNDEFINED",
                  "BEGINNING",
                  "END",
                  "AVERAGED",
                  "SUMMED",
                  "ANNUALIZED",
                  "FORMULA",
                  "HIGH",
                  "LOW")

    if obStr is None:
        return None
    elif obStr.upper() in obsVals:
        return obStr.upper()    
    elif obStr.upper() in ("UNDEFINED", "BEGIN", "END", "AVERAGE", "SUM", "ANNUAL" , "FORMULA", "HIGH", "LOW"):
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

