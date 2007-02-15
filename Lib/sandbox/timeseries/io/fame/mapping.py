

# ---------------------------
# For fametype mapping
import types
import numpy
from timeseries import TimeSeries, Date, DateArray, freq_fromstr


# ---------------------------
# Fame specific constants

HRMODE = 1 # READ        
HCMODE = 2 # CREATE      
HOMODE = 3 # OVERWRITE       
HUMODE = 4 # UPDATE      
HSMODE = 5 # SHARED      
HWMODE = 6 # WRITE       
HDMODE = 7 # DIRECT WRITE    

#** FAME Data Object Classes **

HSERIE = 1 # SERIES   
HSCALA = 2 # SCALAR   
HFRMLA = 3 # FORMULA  
HITEM  = 4 # ITEM     
HGLNAM = 5 # GLNAME   
HGLFOR = 6 # GLFORMULA    

#** FAME Data Object Types **

HUNDFT = 0 # Undefined    
HNUMRC = 1 # NUMERIC  
HNAMEL = 2 # NAMELIST 
HBOOLN = 3 # BOOLEAN  
HSTRNG = 4 # STRING   
HPRECN = 5 # PRECISION    
HDATE  = 6 # General DATE 
HRECRD = 7 # RECORD   

#** FAME Frequencies **

HUNDFX = 0 # Undefined            
HDAILY = 8 # DAILY            
HBUSNS = 9 # BUSINESS              
HWKSUN = 16 #WEEKLY (SUNDAY)
HMONTH = 129 # MONTHLY             
HCASEX = 232 # CASE
HSEC   = 226 # SECONDLY
HMIN   = 227 # MINUTELY
HHOUR  = 228 # HOURLY
HQTOCT = 160 # QUARTERLY (OCTOBER)
HQTNOV = 161 # QUARTERLY (NOVEMBER)
HQTDEC = 162 # QUARTERLY (DECEMBER)
HANJAN = 192 # ANNUAL (JANUARY)
HANFEB = 193 # ANNUAL (FEBRUARY)  
HANMAR = 194 # ANNUAL (MARCH) 
HANAPR = 195 # ANNUAL (APRIL) 
HANMAY = 196 # ANNUAL (MAY)   
HANJUN = 197 # ANNUAL (JUNE)  
HANJUL = 198 # ANNUAL (JULY)  
HANAUG = 199 # ANNUAL (AUGUST)
HANSEP = 200 # ANNUAL (SEPTEMBER) 
HANOCT = 201 # ANNUAL (OCTOBER)
HANNOV = 202 # ANNUAL (NOVEMBER)  
HANDEC = 203 # ANNUAL (DECEMBER)  

#** FAME BASIS Attribute Settings **

HBSUND = 0 # Undefined
HBSDAY = 1 # DAILY
HBSBUS = 2 # BUSINESS

#** FAME OBSERVED Attribute Settings **

HOBUND = 0 # Undefined
HOBBEG = 1 # BEGINNING
HOBEND = 2 # ENDING
HOBAVG = 3 # AVERAGED
HOBSUM = 4 # SUMMED
HOBANN = 5 # ANNUALIZED
HOBFRM = 6 # FORMULA
HOBHI  = 7 # HIGH
HOBLO  = 8 # LOW

def reverse_dict(d):
    return dict([(y, x) for x, y in d.iteritems()])

basisMapping = { HBSUND:"UNDEFINED",
                 HBSDAY:"D",
                 HBSBUS:"B"}
basisReverseMapping = reverse_dict(basisMapping)

observedMapping = { HOBUND:"UNDEFINED",
                  HOBBEG: "BEGINNING",
                  HOBEND: "ENDING",
                  HOBAVG: "AVERAGED",
                  HOBSUM: "SUMMED",
                  HOBANN: "ANNUALIZED",
                  HOBFRM: "FORMULA",
                  HOBHI: "MAXIMUM",
                  HOBLO: "MINIMUM"  }
                  
observedReverseMapping = reverse_dict(observedMapping)

freqMapping = { HDAILY:"D",
                HBUSNS:"B",
                HMONTH:"M",
                HWKSUN:"W",
                HSEC  :"S",
                HMIN  :"T",
                HHOUR :"H",
                HQTOCT:"Q",
                HQTNOV:"Q",
                HQTDEC:"Q",
                HANJAN:"A",
                HANFEB:"A",
                HANMAR:"A",
                HANAPR:"A",
                HANMAY:"A",
                HANJUN:"A",
                HANJUL:"A",
                HANAUG:"A",
                HANSEP:"A",
                HANOCT:"A",
                HANNOV:"A",
                HANDEC:"A" }
                
freqMapping = dict([(x, freq_fromstr(val)) for x, val in freqMapping.iteritems()])

freqReverseMapping = {  "D" : HDAILY,
                        "B" : HBUSNS,
                        "M" : HMONTH,
                        "W" : HWKSUN,
                        "S" : HSEC,
                        "T" : HMIN,
                        "H" : HHOUR,
                        "Q" : HQTDEC,
                        "A" : HANDEC}
                        
freqReverseMapping = dict([(freq_fromstr(x), val) for x, val in freqReverseMapping.iteritems()])

value_adjust = {
    'A':1849,
    'Q':7396,
    'M':22188,
    'W':96477,
    'B':482381,
    'D':675333,
    'H':87648,
    'T':5258880,
    'S':315532800}

value_adjust = dict([(freq_fromstr(x), val) for x, val in value_adjust.iteritems()])


def fametype_fromdata(data):
    """determine fame type code from a data object"""
    
    if isinstance(data, DateArray) or isinstance(data, Date):
        return freqReverseMapping[data.freq]
    elif hasattr(data, 'dtype'):
        dtypeStr = str(data.dtype)
        
        if dtypeStr[:5] == "float":
            if int(dtypeStr[5:]) > 32: return HPRECN
            else: return HNUMRC
        elif dtypeStr[:3] == "int":
            if int(dtypeStr[3:]) > 32: return HPRECN
            else: return HNUMRC
        elif dtypeStr[:4] == "uint":
            if int(dtypeStr[4:]) >= 32: return HPRECN
            else: return HNUMRC
        elif dtypeStr[:2] == "|S" or dtypeStr == 'object':
            return HSTRNG
        elif dtypeStr == "bool":
            return HBOOLN
        else:
            raise ValueError("Unsupported dtype for fame database: %s", dtypeStr)
    
    elif type(data) == types.StringType:
        return HSTRNG
    elif type(data) in (types.IntType, types.FloatType):
        return HPRECN
    elif type(data) == types.BooleanType:
        return HBOOLN
    elif type(data) == types.ListType:
        return HNAMEL
    else:
        raise ValueError("Unrecognized data type")

def fametype_tonumpy(fametype):
    if fametype >= 8:
        # date types
        return numpy.int32
    elif fametype == HNAMEL:
        return None
    else:
        typeMap = {
            HNUMRC:numpy.float32,
            HBOOLN:numpy.int32,
            HSTRNG:numpy.object_,
            HPRECN:numpy.float64}
        return typeMap[fametype]

