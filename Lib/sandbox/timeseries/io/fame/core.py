import re, thread

import numpy
import maskedarray as ma
import timeseries as ts

import cfame
from cfame import FAME_CONSTANTS

__all__ = [
    'FameDb', 'set_option', 'license_expires', 'DBError'
           ]

_g = globals()
for var, val in FAME_CONSTANTS.iteritems():
    _g[var] = val
    
__all__.extend(list(FAME_CONSTANTS))

def reverse_dict(d):
    return dict([(y, x) for x, y in d.iteritems()])

def convert_dict(d, key_func=lambda x:x, val_func=lambda x:x):
    return dict([(key_func(key), val_func(val)) for key, val in d.iteritems()])

basis_map = { HBSUND:"U",
              HBSDAY:"D",
              HBSBUS:"B"}
basis_revmap = reverse_dict(basis_map)

observed_map = { HOBUND:"UNDEFINED",
                 HOBBEG:"BEGINNING",
                 HOBEND:"ENDING",
                 HOBAVG:"AVERAGED",
                 HOBSUM:"SUMMED",
                 HOBANN:"ANNUALIZED",
                 HOBFRM:"FORMULA",
                 HOBHI:"MAXIMUM",
                 HOBLO:"MINIMUM"}
observed_revmap = reverse_dict(observed_map)
observed_revmap['HIGH'] = HOBHI
observed_revmap['LOW'] = HOBLO

def translate_basis(basis):
    "translate user specified basis to FAME constant"

    if isinstance(basis, str):
        freq = ts.freq_tostr(ts.freq_fromstr(basis))
        try:
            return basis_revmap[freq]
        except KeyError:
            raise ValueError("Basis must be " + \
                             "'DAILY', 'BUSINESS', or 'UNDEFINED'")
    else:
        if basis in basis_map: return basis
        elif basis == ts.freq_fromstr('D'): return HBSDAY
        elif basis == ts.freq_fromstr('B'): return HBSBUS
        elif basis == ts.freq_fromstr('U'): return HBSUND
        else:
            raise ValueError("Invalid Basis value")

def translate_observed(observed):
    "translate user specified observed to FAME constant"
    if isinstance(observed, str):
        return observed_revmap[ts.fmtObserv(observed)]
    elif observed in (observed_map):
        return observed
    else:
        raise ValueError("Invalid Observed value")

freq_map = { HDAILY:"D",
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
freq_map = convert_dict(freq_map, val_func=ts.freq_fromstr)

freq_revmap = { "D" : HDAILY,
                "B" : HBUSNS,
                "M" : HMONTH,
                "W" : HWKSUN,
                "S" : HSEC,
                "T" : HMIN,
                "H" : HHOUR,
                "Q" : HQTDEC,
                "A" : HANDEC}
freq_revmap = convert_dict(freq_revmap, key_func=ts.freq_fromstr)

date_value_adjust = { 'A':1849,
                      'Q':7396,
                      'M':22188,
                      'W':96477,
                      'B':482381,
                      'D':675333,
                      'H':87648,
                      'T':5258880,
                      'S':315532800}
date_value_adjust = convert_dict(date_value_adjust, key_func=ts.freq_fromstr)

def fametype_fromdata(data):
    """determine fame type code from a data object"""
    
    if isinstance(data, ts.DateArray) or isinstance(data, ts.Date):
        return freq_revmap[data.freq]
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
    
    elif isinstance(data, str):
        return HSTRNG
    elif isinstance(data, (int, float)):
        return HPRECN
    elif isinstance(data, bool):
        return HBOOLN
    elif isinstance(data, list):
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

class CaseInsensitiveDict(dict):
    def __init__(self, data={}):
        for i, v in data.iteritems():
            self[i.upper()] = v
            
    def __getitem__(self, key):
        if hasattr(key, 'upper'): key = key.upper()
        return super(CaseInsensitiveDict, self).__getitem__(key)
        
    def __setitem__(self, key, item):
        if hasattr(key, 'upper'): key = key.upper()
        super(CaseInsensitiveDict, self).__setitem__(key, item)


def _single_or_multi_func(dbkey, name, func, *args, **kwargs):
    if isinstance(name, str):
        single_obj = True
        name = [name]
    else:
        single_obj = False

    result = {}
    for n in name:
        result[n] = func(dbkey, n, *args, **kwargs)

    if single_obj:
        return result.values()[0]

    return result
    
def _famedate_to_tsdate(fame_date, freqstr):
    "convert integer fame date to a timeseries Date"
    value = fame_date + date_value_adjust[ts.freq_fromstr(freqstr)]
    return ts.Date(freq=freqstr, value=value)


def _fame_params_from_pyobj_scalar(pyobj):
    return {
        'cls':HSCALA,
        'freq':HUNDFX,
        'type':fametype_fromdata(pyobj),
        'basis':HBSUND,
        'observed':HOBUND}

def _fame_params_from_pyobj_tser(pyobj):

    if hasattr(pyobj, "observed"):
        fame_observed = observed_revmap[pyobj.observed]
        if fame_observed == 0: fame_observed = HOBEND
    else:
        fame_observed = HOBEND

    return {
        'cls':HSERIE,
        'freq':freq_revmap[pyobj.freq],
        'type':fametype_fromdata(pyobj._data),
        'basis':HBSDAY,
        'observed':fame_observed}

def _fame_params_from_pyobj_cser(pyobj):
    if hasattr(pyobj, "_data"):
        fame_data = pyobj._data
    else:
        fame_data = pyobj

    return {
        'cls':HSERIE,
        'freq':HCASEX,
        'type':fametype_fromdata(fame_data),
        'basis':HBSUND,
        'observed':HOBUND}


class DBError(Exception): pass

class FameDb(object):
    """Fame database object.

:Construction:
    x = FameDb(conn_str, mode='r')

:Paramaters:
    - `conn_str` (str) : valid connection string. Can be a physical path,
    or channel specification, etc. See FAME documentation on cfmopdb for
    valid connection strings.
    - `mode` (str, *['r']*) : method of access to the database. Can be one
    of the following:
        'r' => read only
        's' => shared
        'o' => overwrite
        'c' => create
        'u' => update
        'w' => write
        'd' => direct

Notes
    - For changes to be posted, you must explictly use the "post" or
      "close" methods (changes are posted on close)."""

    def __init__(self, conn_str, mode='r'):
        mode = mode.lower()
        if mode == 'r':
            intmode = HRMODE
        elif mode == 's':
            intmode = HSMODE
        elif mode == 'u':
            intmode = HUMODE
        elif mode == 'w':
            intmode = HWMODE
        elif mode == 'd':
            intmode = HDMODE
        elif mode == 'c':
            intmode = HCMODE
        elif mode == 'o':
            intmode = HOMODE
        else:
            raise ValueError, "Database access mode not supported."
        self.mode = mode
        
        self.dbkey = cf_open(conn_str, intmode)

        
    def read(self, name,
             start_date=None, end_date=None,
             start_case=None, end_case=None, max_string_len=65):
    
        """read specified object(s) from database

:Parameters:
        - `name` (string or list of strings) : names of objects that will be
          read from the database

        - `start_date` (int, *[None]*) : Applies only when reading time series.
          If specified, only data points on or after `start_date` will be read.
          If None, data will be read from the first value of the series.
        - `end_date` (int, *[None]*) : Applies only when reading time series.
          If specified, only data points on or before `end_date` will be read.
          If None, data will be read to the last value of the series.
        - `start_case` (int, *[None]*) : Applies only when reading case series.
          If specified, only data points on or after `start_case` will be read.
          If None, data will be read starting from case index 1
        - `end_case` (int, *[None]*) : Applies only when reading case series.
          If specified, only data points on or before `end_case` will be read.
          If None, data will be read to the last value of the series.
        - `max_string_len` (int, *[65]*) : Applies only when readings strings
           or series of strings. This is the maximum length of string that can
           be read. Lower values result in less memory usage, so you should
           specify this as low as is reasonable for your data.
           
:Return:
        if `name` is a list of strings:
            case insensitive dictionary of the objects
        if `name` is a single string:
            object from database that is stored as `name`"""

        isSingle = False
        if isinstance(name, str):
            names = [name]
            isSingle = True
        else:
            names = name

        items = CaseInsensitiveDict()
        
        #default to -1. This will get the entire range
        _start_case = _end_case = -1
        _start_date = _end_date = -1

        range_freq = None
        if start_date is not None:
            _start_date = start_date.value - date_value_adjust[start_date.freq]
            range_freq = freq_revmap[start_date.freq]

        if end_date is not None:
            if start_date is not None and start_date.freq != end_date.freq:
                raise ValueError("start_date and end_date must be same frequency")
            _end_date = end_date.value - date_value_adjust[end_date.freq]
            if range_freq is None:
                range_freq = freq_revmap[end_date.freq]

        if start_case is not None: _start_case = start_case
        if end_case is not None: _end_case = end_case
       
        if len(set([_start_case, _end_case, _start_date, _end_date, -1])) != 1:
            checkFreq = True
        else:
            checkFreq = False

        for objName in names:
            objName = objName.upper()

            if checkFreq:
                objFreq = self.obj_size(objName)['freq']

                if objFreq == range_freq:
                    start_index, end_index = _start_date, _end_date
                elif objFreq == HCASEX:
                    start_index, end_index = _start_case, _end_case
                else:
                    start_index, end_index = -1, -1
            else:
                start_index, end_index = -1, -1

            result = cf_read(self.dbkey, objName, start_index,
                             end_index, max_string_len)

            if result['type'] == HBOOLN:
                numpyType = numpy.bool_
            else:
                numpyType = fametype_tonumpy(result['type'])

            if result['type'] == HNAMEL:
                pyObj = [x for x in result['data'][1:-1].split(", ") \
                         if x != '']

            elif result['class'] == HSCALA:
                if isinstance(result['data'], str):
                    if result['mask']:
                        pyObj = None
                    else:
                        pyObj = result['data']
                else:
                    if result['mask'][0]:
                        pyObj = None
                    else:
                        pyObj = result['data'][0]
                        if result['type'] >= 8: # date type
                            value = pyObj+ \
                               date_value_adjust[freq_map[result['type']]]
                            pyObj = ts.Date(
                                        freq=freq_map[result['type']],
                                        value=value)
                        else:
                            pyObj = numpyType(pyObj)

            elif result['class'] == HSERIE:
                
                if 'data' in result:
                    vals = result['data']
                    mask = result['mask']
                    if not mask.any(): mask = ma.nomask
                else:
                    vals = []
                    mask = ma.nomask
                    
                if result['type'] >= 8: # date type
                    valadj = date_value_adjust[freq_map[result['type']]]
                    if len(vals) > 0: vals += valadj
                    data = ts.DateArray(vals,
                                        freq=freq_map[result['type']])
                else:
                    data = numpy.array(vals, dtype=numpyType)
                    
                if result['freq'] == HCASEX:
                    pyObj = ma.array(data, mask=mask)
                else:
                    observed = observed_map[result['observed']]
                    freq = freq_map[result['freq']]

                    if 'data' in result:
                        start_date = ts.Date(
                              freq=freq,
                              value=result['startindex']+date_value_adjust[freq])
                    else:
                        start_date = None
                    
                    pyObj = ts.time_series(data, freq=freq,
                                           start_date=start_date,
                                           observed=observed, mask=mask)

            items[objName] = pyObj

        if isSingle:
            return items.values()[0]
            
        return items


    def write_tser_dict(self, objdict,
                        overwrite=False, assume_exists=False,
                        start_date=None, end_date=None):
        """for each key, value pair in the dictionary `objdict` write value to
the database as key, as a time series (calls FameDb.write_tser on each key,
value pair)

:Parameters:
        - `objdict` (dict) : dictionary of TimeSeries objects to be written. Object
          names for keys and TimeSeries objects for values
        - `overwrite (boolean, *[False]*) : If True, if the key exists in the database it
           will be overwritten. If False, data will be added to series that already exist
           (data in objects in `objdict` will be given priority over pre-existing data in
           the db where there is overlap)
        - `assume_exists` (boolean, *[False]*) : If True, an error will be
           raised if the series does not exist. If False, the series will be
           created if it does not exist already.
        - `start_date` (Date, *[None]*) : If None, data will be written from the start of
           the series. If specified, only data points on or after start_date will be written.
        - `end_date` (Date, *[None]*) : If None, data will be written until the end of
           the series. If specified, only data points on or before end_date will be written.
"""
        for key, obj in objdict.iteritems():
            self.write_tser(key, obj, overwrite=overwrite,
                            assume_exists=assume_exists,
                            start_date=start_date, end_date=end_date)


    def write_cser_dict(self, objdict,
                        overwrite=False, assume_exists=False,
                        zero_represents=1, start_case=None, end_case=None):
        """for each key, value pair in the dictionary `objdict` write value to
the database as key, as a case series (calls FameDb.write_tser on each key,
value pair)

:Parameters:
        - `objdict` (dict) : dictionary of arrays to be written as Case Series.
           Object names for keys and arrays for values
        - `overwrite (boolean, *[False]*) : If True, if the key exists in the database it
           will be overwritten. If False, data will be added to series that already exist
           (data in objects in `objdict` will be given priority over pre-existing data in
           the db where there is overlap)
        - `assume_exists` (boolean, *[False]*) : If True, an error will be
           raised if the series does not exist. If False, the series will be
           created if it does not exist already.
        - `zero_represents` (int, *[1]*) : the case index for FAME that index zero in
           the array represents
        - `start_case` (int, *[None]*) : If None, data will be written from the start of
           the array. If specified, only data points on or after start_case will be written.
        - `end_case` (int, *[None]*) : If None, data will be written until the end of
           the array. If specified, only data points on or before end_case will be written.
"""
        for key, obj in objdict.iteritems():
            self.write_cser(key, obj, overwrite=overwrite,
                            assume_exists=assume_exists,
                            zero_represents=zero_represents,
                            start_case=start_case, end_case=end_case)

    def write_scalar_dict(self, objdict):
        """for each key, value pair in the dictionary `objdict` write value to
the database as key, as a scalar (calls FameDb.write_scalar on each key,
value pair)

:Parameters:
        - `objdict` (dict) : dictionary of items to be written as scalars.
           Object names for keys and scalar items for values
"""
        for key, obj in objdict.iteritems():
            self.write_scalar(key, obj)


    def write_tser(self, name, tser,
                   overwrite=False, assume_exists=False,
                   start_date=None, end_date=None):
        """write `tser` to the database as `name` as a time series.

:Parameters:
        - `name` (string) : database key that the object will be written to
        - `tser` (TimeSeries) : TimeSeries object to be written. Cannot have missing dates.
           Use fill_missing_dates first on your series if you suspect this is the situation.
           TimeSeries must be 1-dimensional
        - `overwrite (boolean, *[False]*) : If True, if `name` exists in the database it
           will be overwritten. If False, data will be added to series that already exist
           (data in `tser` will be given priority over pre-existing data in the db where
           there is overlap)
        - `assume_exists` (boolean, *[False]*) : If True, an error will be
           raised if the series does not exist. If False, the series will be
           created if it does not exist already.
        - `start_date` (Date, *[None]*) : If None, data will be written from the start of
           `tser`. If specified, only data points on or after start_date will be written.
        - `end_date` (Date, *[None]*) : If None, data will be written until the end of
           `tser`. If specified, only data points on or before end_date will be written.
"""
            
        if not isinstance(tser, ts.TimeSeries):
            raise ValueError("tser is not a valid time series")
        elif tser.has_missing_dates():
            raise ValueError("tser must not have any missing dates")
        elif tser.ndim != 1:
            raise ValueError("FAME db only supports 1-dimensional time series")

        exists = self.obj_exists(name)

        if assume_exists and not exists:
            raise DBError("%s does not exist" % name)

        if overwrite or not exists: create = True
        else: create = False
        
        fame_params = _fame_params_from_pyobj_tser(tser)
        
        fame_cls = fame_params['cls']
        fame_type = fame_params['type']
        fame_freq = fame_params['freq']
        fame_basis = fame_params['basis']
        fame_observed = fame_params['observed']
        
        if create:
            if exists: self.delete_obj(name)
            cf_create(self.dbkey, name, fame_cls, fame_freq, fame_type, fame_basis, fame_observed)

        def get_boundary_date(bdate, attr):
            if bdate is not None:
                if bdate.freq != tser.freq:
                    raise ValueError(attr+" frequency must be same as tser frequency")
                if tser.start_date > bdate or tser.end_date < bdate:
                    raise ValueError(attr+" outside range of series")
                return bdate
            else:
                return getattr(tser, attr)
            
        start_date = get_boundary_date(start_date, "start_date")
        end_date = get_boundary_date(end_date, "end_date")
        
        if start_date is not None:

            towrite = tser[start_date:end_date+1]

            start_index = start_date.value
            end_index = end_date.value

            # convert integer types to floats since FAME does not have an integer type
            newType = fametype_tonumpy(fame_type)
            if fame_type >= 8:
                # date type
                fame_data = towrite._data - date_value_adjust[towrite._data.freq]
            elif newType != tser._data.dtype:
                fame_data = towrite._data.astype(newType)
            else:
                fame_data = towrite._data

            if towrite._mask is ma.nomask:
                fame_mask = numpy.zeros(towrite._data.shape, dtype=numpy.bool_)
            else:
                fame_mask = towrite._mask

            start_index -= date_value_adjust[towrite.freq]
            end_index   -= date_value_adjust[towrite.freq]

            cfame.write_series(self.dbkey, name, fame_data, fame_mask, start_index, end_index, fame_type, fame_freq)

    def write_cser(self, name, cser, overwrite=False, assume_exists=False, zero_represents=1, start_case=None, end_case=None):
        """write `cser` to the database as `name` as a case series.

:Parameters:
        - `name` (string) : database key that the object will be written to
        - `cser` (ndarray) : 1-dimensional ndarray (or subclass of ndarray) object to be
           written. If `cser` is a MaskedArray, then masked values will be written as ND.
        - `overwrite (boolean, *[False]*) : If True, if `name` exists in the database it
           will be overwritten. If False, data will be added to series that already exist
           (data in `cser` will be given priority over pre-existing data in the db where
           there is overlap)
        - `assume_exists` (boolean, *[False]*) : If True, an error will be
           raised if the series does not exist. If False, the series will be
           created if it does not exist already.
        - `zero_represents` (int, *[1]*) : the case index for FAME that index zero in
           the array represents
        - `start_case` (int, *[None]*) : If None, data will be written from the start of
           `cser`. If specified, only data points on or after start_case will be written.
        - `end_case` (int, *[None]*) : If None, data will be written until the end of
           `cser`. If specified, only data points on or before end_case will be written.
"""
            
        if not isinstance(cser, numpy.ndarray):
            raise ValueError("cser is not a valid ndarray")
        elif cser.ndim != 1:
            raise ValueError("FAME db only supports 1-dimensional arrays")

        exists = self.obj_exists(name)
        if assume_exists and not exists:
            raise DBError("%s does not exist" % name)

        if overwrite or not exists: create = True
        else: create = False
        
        fame_params = _fame_params_from_pyobj_cser(cser)
        
        fame_cls = fame_params['cls']
        fame_type = fame_params['type']
        fame_freq = fame_params['freq']
        fame_basis = fame_params['basis']
        fame_observed = fame_params['observed']

        if hasattr(cser, "_data"):
            fame_data = cser._data
            if cser._mask is ma.nomask:
                fame_mask = numpy.zeros(fame_data.shape, dtype=numpy.bool_)
            else:
                fame_mask = cser._mask
        else:
            fame_data = cser
            fame_mask = numpy.zeros(fame_data.shape, dtype=numpy.bool_)

        if create:
            if exists: self.delete_obj(name)
            cf_create(self.dbkey, name, fame_cls, fame_freq, fame_type, fame_basis, fame_observed)

        def get_boundary_case(bcase, attr):
            if bcase is not None:
                idx = bcase - zero_represents
                if idx < 0 or idx > cser.size:
                    raise ValueError("%s outside range of series" % attr)
                return bcase
            else:
                if cser.size == 0:
                    return None
                else:
                    if attr == 'start_case':
                        return zero_represents
                    elif attr == 'end_case':
                        return zero_represents + cser.size - 1
                    else:
                        raise ValueError("unexpected argument: %s " % attr)
            
        start_case = get_boundary_case(start_case, "start_case")
        end_case = get_boundary_case(end_case, "end_case")

        if start_case is not None:        
            # convert integer types to floats since FAME does not have an integer type
            s = start_case - zero_represents
            e = end_case - zero_represents
            
            fame_data = fame_data[s:e+1]
            fame_mask = fame_mask[s:e+1]
            newType = fametype_tonumpy(fame_type)
            if fame_type >= 8:
                # date type
                fame_data = fame_data - date_value_adjust[fame_data.freq]
            elif newType != fame_data.dtype:
                fame_data = fame_data.astype(newType)

            cfame.write_series(self.dbkey, name, fame_data, fame_mask, start_case, end_case, fame_type, fame_freq)


    def write_scalar(self, name, scalar):
        """write `scalar` to the database as `name` as a scalar object. If an
object already exists in the database named as `name` then it is
over-written, otherwise it is created.

:Parameters:
        - `name` (string) : database key that the object will be written to
        - `scalar` : one of the following: string, numpy scalar, int, float,
           list of strings (for name lists), Date, boolean"""
        
        fame_params = _fame_params_from_pyobj_scalar(scalar)
        fame_type = fame_params['type']

        if isinstance(scalar, ts.Date):
            fame_data = numpy.int32(scalar.value - date_value_adjust[scalar.freq])
        elif hasattr(scalar, "dtype"):
            if scalar.ndim != 0: raise ValueError("received non-scalar data")
            newType = fametype_tonumpy(fame_type)
            if newType != scalar.dtype: fame_data = scalar.astype(newType)
            else: fame_data = scalar
        elif fame_type == HSTRNG:
            fame_data = scalar
        elif fame_type == HPRECN:
            fame_data = numpy.float64(scalar)
        elif fame_type == HBOOLN:
            fame_data = numpy.int32(scalar)
        elif fame_type == HNAMEL:
            fame_data = "{" + ", ".join(scalar) + "}"
        else:
            raise ValueError("Unrecognized data type")
            
        if self.obj_exists(name): self.delete_obj(name)
        cf_create(self.dbkey, name,
                    fame_params['cls'],
                    fame_params['freq'],
                    fame_params['type'],
                    fame_params['basis'],
                    fame_params['observed'])

        # convert integer types to floats since FAME does not have an integer type
        newType = fametype_tonumpy(fame_type)
        if hasattr(fame_data, 'dtype') and newType != fame_data.dtype:
            fame_data = fame_data.astype(newType)
        
        if fame_type == HNAMEL:
            cf_write_namelist(self.dbkey, name, fame_data)
        else:
            cf_write_scalar(self.dbkey, name, fame_data, fame_type)
            
    def _create_obj(self, name, cls, type, freq=None, basis=None, observed=None):
        """create object in database with specified attributes as `name`.

You must use the fame constants defined in mapping.py for each of the
parameters. Generally speaking, it is easier to use initialize_obj
with a prototype object.
"""
        
        if cls not in (HSERIE, HSCALA):
            raise ValueError("unrecognized object class: "+str(cls))
        
        if freq is None:
            if cls == HSCALA:
                freq = HUNDFX
            else:
                raise ValueError("freq must be specified for series")

        if freq in (HUNDFX, HCASEX):
            basis = HBSUND
            observed = HOBUND
        else:
            if basis is None: basis = HBSDAY
            if observed is None: observed = HOBEND

        cf_create(self.dbkey, name, cls, freq, type, basis, observed)
        
    def initialize_obj(self, name, pyobj):
        """initialize object of appropriate type in database based on the
python object `pyobj` as `name`. Does not write any data to the
database, simply initializes the object in the database."""
        if isinstance(pyobj, ts.TimeSeries):
            param_func = _fame_params_from_pyobj_tser
        elif isinstance(pyobj, numpy.ndarray):
            param_func = _fame_params_from_pyobj_cser
        else:
            param_func = _fame_params_from_pyobj_scalar
            
        fame_params = param_func(pyobj)
        cf_create(self.dbkey, name,
                    fame_params['cls'],
                    fame_params['freq'],
                    fame_params['type'],
                    fame_params['basis'],
                    fame_params['observed'])


    def db_desc(self):
        "get 'description' attribute of database"
        return cf_get_db_attr(self.dbkey, "DESC")
        
    def db_doc(self):
        "get 'doc' attribute of database"
        return cf_get_db_attr(self.dbkey, "DOC")
        
    def db_created(self):
        "get 'created' attribute of database"
        fame_date = cf_get_db_attr(self.dbkey, "CREATED")
        return _famedate_to_tsdate(fame_date, 's')
        
    def db_modified(self):
        "get 'modified' attribute of database"
        fame_date = cf_get_db_attr(self.dbkey, "MODIFIED")
        return _famedate_to_tsdate(fame_date, 's')

    def set_db_desc(self, desc):
        "set description attribute of database"
        cf_set_db_desc(self.dbkey, desc)

    def set_db_doc(self, doc):
        "set doc attribute of database"
        cf_set_db_doc(self.dbkey, doc)

    def db_is_open(self):
        "returns True if database is open. False otherwise"
        return cf_get_db_attr(self.dbkey, "ISOPEN")

    def wildlist(self, exp, wildonly=False):
        """performs a wildlist lookup on the database, using Fame syntax
("?" and "^"), returns a normal python list of strings"""
        res = cf_wildlist(self.dbkey, exp)
            
        if wildonly:
            exp = exp.replace("?", "(.*)")
            exp = exp.replace("^", "(.)")
            exp = exp.replace("$","\$")
            regex = re.compile(exp)
            res = ["".join(regex.match(res[i]).groups()) \
                   for i in range(len(res))]
        return res

    def close(self):
        """Closes the database. Changes will be posted."""
        if self.db_is_open():
            cf_close(self.dbkey)
            
    def post(self):
        cf_post(self.dbkey)

    def restore(self):
        """Discard any changes made to the database since it was last opened or posted."""
        cf_restore(self.dbkey)

    def obj_exists(self, name):
        return cf_exists(self.dbkey, name)
        
    def delete_obj(self, name, must_exist=True):
        """Deletes the specified object(s) from the database"""
        if isinstance(name, str): name = [name]
        [cf_delete_obj(self.dbkey, n) for n in name if must_exist or self.obj_exists(n)]

    def obj_size(self, name):
        """basic information about the size of an object(s) in a database"""
        return _single_or_multi_func(self.dbkey, name, cf_obj_size)
        
    def obj_freq(self, name):
        """get frequency of a FAME time series object in the database"""
        obj_sz = self.obj_size(name)
        fame_freq = obj_sz['freq']
        if fame_freq == 0: return None
        return freq_map[obj_sz['freq']]

    def __ser_date(self, name, date_type):
        """helper method for start_date and end_date"""
        obj_sz = self.obj_size(name)
        fame_freq = obj_sz['freq']
        if fame_freq == 0: return None

        try:
            ts_freq = freq_map[obj_sz['freq']]
        except KeyError:
            raise DBError("unsupported FAME frequency: %i", fame_freq)
        
        if obj_sz[date_type+'_year'] == -1: return None
            
        annDate = ts.Date(freq='A', year=obj_sz[date_type+'_year'])
        return annDate.asfreq(ts_freq, relation='BEFORE') + (obj_sz[date_type+'_period'] - 1)

    def obj_start_date(self, name):
        """get start_date of a FAME time series object"""
        return self.__ser_date(name, 'start')

    def obj_end_date(self, name):
        """get end_date of a FAME time series object"""
        return self.__ser_date(name, 'end')
        
    def obj_created(self, name):
        "get 'created' attribute of object in database"
        fame_date = cf_get_obj_attr(self.dbkey, name, "CREATED")
        return _famedate_to_tsdate(fame_date, 's')

    def obj_modified(self, name):
        "get 'modified' attribute of object in database"
        fame_date = cf_get_obj_attr(self.dbkey, name, "MODIFIED")
        return _famedate_to_tsdate(fame_date, 's')

    def set_obj_desc(self, name, desc):
        "set 'description' attribute of object in database"
        cf_set_obj_desc(self.dbkey, name, desc)

    def set_obj_doc(self, name, doc):
        "set 'documentation' attribute of object in database"
        cf_set_obj_doc(self.dbkey, name, doc)

    def set_obj_basis(self, name, basis):
        "set 'basis' attribute of object in database"
        basis = translate_basis(basis)
        cf_set_obj_basis(self.dbkey, name, basis)

    def set_obj_observed(self, name, observed):
        "set 'observed' attribute of object in database"
        observed = translate_observed(observed)
        cf_set_obj_observed(self.dbkey, name, observed)

    def _whats(self, name):
        """Preforms a fame "whats" command on the provided name(s)

Note: Returns FAME constants which are not directly interpretable
in the context of the timeseries module. For this reason, it is
recommended that you use the obj_* methods to retrieve the desired
information about an object.
"""
        return _single_or_multi_func(self.dbkey, name, cf_whats)
        
    def obj_desc(self, name):
        """get desc attribute for an object"""
        return self._whats(name)['desc']

    def obj_doc(self, name):
        """get doc attribute for an object"""
        return self._whats(name)['doc']

    def rename_obj(self, name, new_name):
        """rename fame object(s) in database"""
        if isinstance(name, str): name = [name]
        if isinstance(new_name, str): new_name = [new_name]

        if len(name) != len(new_name):
            raise ValueError("number of original names does not match " + \
                             "number of new names")
                             
        for n, o in zip(name, new_name):
            cf_rename_obj(self.dbkey, n, o)
        
    def copy_obj(self, target_db, source_name, target_name=None):
        """copy fame object(s) to another destination"""
        if target_name is None: target_name = source_name
        if isinstance(source_name, str): source_name = [source_name]
        if isinstance(target_name, str): target_name = [target_name]

        if len(source_name) != len(target_name):
            raise ValueError("number of source names does not match " + \
                             "number of target names")

        for s, t in zip(source_name, target_name):
            cf_copy(self.dbkey, target_db.dbkey, s, t)


class cFameCall:
    """wrapper for cfame functions that acquires and releases a resource lock.
This is needed because the Fame C api is not thread safe."""

    fameLock = thread.allocate_lock()

    def __init__ (self, func):
        self.f = func
        self.__doc__ = getattr(func, "__doc__", str(func))
        self.__name__ = getattr(func, "__name__", str(func))

    def __call__ (self, *args, **kwargs):
        "Execute the call behavior."
        tmp = self.fameLock.acquire()
        try:
            result = self.f(*args, **kwargs)
            self.fameLock.release()
        except:
            self.fameLock.release()
            raise
            
        return result

cf_open = cFameCall(cfame.open)
cf_set_option = cFameCall(cfame.set_option)
cf_close = cFameCall(cfame.close)
cf_post = cFameCall(cfame.post)
cf_restore = cFameCall(cfame.restore)
cf_obj_size = cFameCall(cfame.obj_size)
cf_whats = cFameCall(cfame.whats)
cf_delete_obj = cFameCall(cfame.delete_obj)
cf_create = cFameCall(cfame.create)
cf_read = cFameCall(cfame.read)
cf_write_scalar = cFameCall(cfame.write_scalar)
cf_write_series = cFameCall(cfame.write_series)
cf_write_namelist = cFameCall(cfame.write_namelist)
cf_wildlist = cFameCall(cfame.wildlist)
cf_exists = cFameCall(cfame.exists)
cf_get_db_attr = cFameCall(cfame.get_db_attr)
cf_get_obj_attr = cFameCall(cfame.get_obj_attr)
cf_copy = cFameCall(cfame.copy)
cf_rename_obj = cFameCall(cfame.rename_obj)
cf_license_expires = cFameCall(cfame.license_expires)
cf_set_db_desc = cFameCall(cfame.set_db_desc)
cf_set_db_doc = cFameCall(cfame.set_db_doc)
cf_set_obj_desc = cFameCall(cfame.set_obj_desc)
cf_set_obj_doc = cFameCall(cfame.set_obj_doc)
cf_set_obj_basis = cFameCall(cfame.set_obj_basis)
cf_set_obj_observed = cFameCall(cfame.set_obj_observed)


set_option = cf_set_option
set_option.__doc__ = \
"""Set an option in the C HLI. See the FAME documentation for cfmsopt for a
listing of allowable option settings.

:Parameters:
    - option (str) : name of the option to set
    - setting (str) : value of the option to set
    
:Example:
    set_option("DBSIZE", "LARGE")
"""

def license_expires():
    """get date that license expires on"""
    fame_date = cf_license_expires()
    adj_val = date_value_adjust[ts.freq_fromstr('D')]
    return ts.Date(freq='D', value=fame_date+adj_val)
