import sys, types, re, os

import timeseries as ts
import cfame
import mapping as mp

import numpy
import maskedarray as ma
import thread

fameLock = thread.allocate_lock()

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

class DBError(Exception): pass

class FameDb(object):
    """Fame database object

:Construction:
    x = FameDb(conn_str, mode='r')

:Paramaters:
    - `conn_str` (str) : valid connection string. Can be a physical path,
    channel specification, etc.
    - `mode` (str, *['r']*) : method of access to the database. Can be one
    of the following:
        'r' => read only
        's' => shared
        'o' => overwrite
        'c' => create
        'u' => update
        'w' => write
        'd' => direct"""
    def __init__(self, conn_str, mode='r'):
        mode = mode.lower()
        if mode == 'r':
            intmode = mp.HRMODE
        elif mode == 's':
            intmode = mp.HSMODE
        elif mode == 'u':
            intmode = mp.HUMODE
        elif mode == 'w':
            intmode = mp.HWMODE
        elif mode == 'd':
            intmode = mp.HDMODE
        elif mode == 'c':
            intmode = mp.HCMODE
        elif mode == 'o':
            intmode = mp.HOMODE
        else:
            raise ValueError, "Database access mode not supported."
        self.mode = mode
        
        try:
            self.dbkey = cf_open(conn_str, intmode)
            self.dbIsOpen = True
        except:
            self.dbIsOpen = False
            raise

        
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


        if not self.dbIsOpen:
            raise DBError("Database is not open")

        isSingle = False
        if isinstance(name, types.StringType):
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
            _start_date = start_date.value - mp.value_adjust[start_date.freq]
            range_freq = mp.freqReverseMapping[start_date.freq]

        if end_date is not None:
            if start_date is not None and start_date.freq != end_date.freq:
                raise ValueError("start_date and end_date must be same frequency")
            _end_date = end_date.value - mp.value_adjust[end_date.freq]
            if range_freq is None:
                range_freq = mp.freqReverseMapping[end_date.freq]

        if start_case is not None: _start_case = start_case
        if end_case is not None: _end_case = end_case
       
        if len(set([_start_case, _end_case, _start_date, _end_date, -1])) != 1:
            checkFreq = True
        else:
            checkFreq = False

        for objName in names:
            objName = objName.upper()

            if checkFreq:
                objFreq = self.get_freq(objName)

                if objFreq == range_freq:
                    start_index, end_index = _start_date, _end_date
                elif objFreq == mp.HCASEX:
                    start_index, end_index = _start_case, _end_case
                else:
                    start_index, end_index = -1, -1
            else:
                start_index, end_index = -1, -1

            result = cf_read(self.dbkey, objName, start_index,
                             end_index, max_string_len)

            if result['type'] == mp.HBOOLN:
                numpyType = numpy.bool_
            else:
                numpyType = mp.fametype_tonumpy(result['type'])

            if result['type'] == mp.HNAMEL:
                pyObj = [x for x in result['data'][1:-1].split(", ") \
                         if x != '']

            elif result['class'] == mp.HSCALA:
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
                               mp.value_adjust[mp.freqMapping[result['type']]]
                            pyObj = ts.Date(
                                        freq=mp.freqMapping[result['type']],
                                        value=value)
                        else:
                            pyObj = numpyType(pyObj)

            elif result['class'] == mp.HSERIE:
                
                if 'data' in result:
                    vals = result['data']
                    mask = result['mask']
                else:
                    vals = []
                    mask = ma.nomask
                    
                if result['type'] >= 8: # date type
                    valadj = mp.value_adjust[mp.freqMapping[result['type']]]
                    if len(vals) > 0: vals += valadj
                    data = ts.DateArray(vals,
                                        freq=mp.freqMapping[result['type']])
                else:
                    data = numpy.array(vals, dtype=numpyType)
                    
                if result['freq'] == mp.HCASEX:
                    pyObj = ma.array(data, mask=mask)
                else:
                    observed = mp.observedMapping[result['observed']]
                    basis = mp.basisMapping[result['basis']]
                    freq = mp.freqMapping[result['freq']]

                    if 'data' in result:
                        start_date = ts.Date(
                              freq=freq,
                              value=result['startindex']+mp.value_adjust[freq])
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
        
        self.__check_writeable()
            
        if not isinstance(tser, ts.TimeSeries):
            raise ValueError("tser is not a valid time series")
        elif tser.has_missing_dates():
            raise ValueError("tser must not have any missing dates")
        elif tser.ndim != 1:
            raise ValueError("FAME db only supports 1-dimensional time series")

        if assume_exists and not self.exists(name):
            raise DBError("%s does not exist" % name)

        if overwrite or not self.exists(name): create = True
        else: create = False

        fame_type = mp.fametype_fromdata(tser._data)
        fame_freq = mp.freqReverseMapping[tser.freq]

        if create:
            
            if hasattr(tser, "basis"):
                fame_basis = mp.basisReverseMapping[tser.basis]
            else:
                fame_basis = mp.HBSDAY

            if hasattr(tser, "observed"):
                fame_observed = mp.observedReverseMapping[tser.observed]
                if fame_observed == 0: fame_observed = mp.HOBEND
            else:
                fame_observed = mp.HOBEND

            if self.exists(name): self.remove(name)
            cf_create(self.dbkey, name, mp.HSERIE, fame_freq, fame_type, fame_basis, fame_observed)

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
            newType = mp.fametype_tonumpy(fame_type)
            if fame_type >= 8:
                # date type
                fame_data = towrite._data - mp.value_adjust[towrite._data.freq]
            elif newType != tser._data.dtype:
                fame_data = towrite._data.astype(newType)
            else:
                fame_data = towrite._data

            if towrite._mask is ma.nomask:
                fame_mask = numpy.zeros(towrite._data.shape, dtype=numpy.bool_)
            else:
                fame_mask = towrite._mask

            start_index -= mp.value_adjust[towrite.freq]
            end_index   -= mp.value_adjust[towrite.freq]

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
        
        self.__check_writeable()
            
        if not isinstance(cser, numpy.ndarray):
            raise ValueError("cser is not a valid ndarray")
        elif cser.ndim != 1:
            raise ValueError("FAME db only supports 1-dimensional arrays")

        if assume_exists and not self.exists(name):
            raise DBError("%s does not exist" % name)

        if overwrite or not self.exists(name): create = True
        else: create = False

        if hasattr(cser, "_data"):
            fame_data = cser._data
            if cser._mask is ma.nomask:
                fame_mask = numpy.zeros(fame_data.shape, dtype=numpy.bool_)
            else:
                fame_mask = cser._mask
        else:
            fame_data = cser
            fame_mask = numpy.zeros(fame_data.shape, dtype=numpy.bool_)
            
        fame_type = mp.fametype_fromdata(fame_data)

        if create:
            if self.exists(name): self.remove(name)
            cf_create(self.dbkey, name, mp.HSERIE, mp.HCASEX, fame_type, mp.HBSUND, mp.HOBUND)

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
            newType = mp.fametype_tonumpy(fame_type)
            if fame_type >= 8:
                # date type
                fame_data = fame_data - mp.value_adjust[fame_data.freq]
            elif newType != fame_data.dtype:
                fame_data = fame_data.astype(newType)

            cfame.write_series(self.dbkey, name, fame_data, fame_mask, start_case, end_case, fame_type, mp.HCASEX)


    def write_scalar(self, name, scalar):
        """write `scalar` to the database as `name` as a scalar object. If an
object already exists in the database named as `name` then it is
over-written, otherwise it is created.

:Parameters:
        - `name` (string) : database key that the object will be written to
        - `scalar` : one of the following: string, numpy scalar, int, float,
           list of strings (for name lists), Date, boolean"""
        
        self.__check_writeable()
        
        fame_type = mp.fametype_fromdata(scalar)

        if isinstance(scalar, ts.Date):
            fame_data = numpy.int32(scalar.value - mp.value_adjust[scalar.freq])
        elif hasattr(scalar, "dtype"):
            if scalar.ndim != 0: raise ValueError("received non-scalar data")
            newType = mp.fametype_tonumpy(fame_type)
            if newType != scalar.dtype: fame_data = scalar.astype(newType)
            else: fame_data = scalar
        elif fame_type == mp.HSTRNG:
            fame_data = scalar
        elif fame_type == mp.HPRECN:
            fame_data = numpy.float64(scalar)
        elif fame_type == mp.HBOOLN:
            fame_data = numpy.int32(scalar)
        elif fame_type == mp.HNAMEL:
            fame_data = "{" + ", ".join(scalar) + "}"
        else:
            raise ValueError("Unrecognized data type")
            
        if self.exists(name): self.remove(name)
        cf_create(self.dbkey, name, mp.HSCALA, mp.HUNDFX, fame_type, mp.HBSUND, mp.HOBUND)

        # convert integer types to floats since FAME does not have an integer type
        newType = mp.fametype_tonumpy(fame_type)
        if hasattr(fame_data, 'dtype') and newType != fame_data.dtype:
            fame_data = fame_data.astype(newType)
        
        if fame_type == mp.HNAMEL:
            cf_write_namelist(self.dbkey, name, fame_data)
        else:
            cf_write_scalar(self.dbkey, name, fame_data, fame_type)



    def wildlist(self, exp, wildonly=False):
        """performs a wildlist lookup on the database, using Fame syntax
("?" and "^"), returns a normal python list of strings"""
        self.__check_readable()
        res = cf_wildlist(self.dbkey, exp)
            
        if wildonly:
            exp = exp.replace("?", "(.*)")
            exp = exp.replace("^", "(.)")
            exp = exp.replace("$","\$")
            regex = re.compile(exp)
            for i in range(len(res)):
                res[i] = "".join(regex.match(res[i]).groups())
        return res

    def exists(self, objName):
        return cf_exists(self.dbkey, objName)

    def close(self):
        if self.dbIsOpen:
            cf_close(self.dbkey)
        self.dbIsOpen = False

    def __del__(self):
        if self.dbIsOpen:
            self.close()


    def __check_writeable(self):
        """Raises error if data base is not writeable"""
        if not self.dbIsOpen:
            raise DBError("Database is not open")
        if self.mode == 'r':
            raise DBError("Cannot write to a read-only database")

    def __check_readable(self):
        """Raises error if data base is not readable"""
        if not self.dbIsOpen:
            raise DBError("Database is not open")


    def remove(self, name, ignoreError=True):
        """Deletes the given series from the database"""
        if type(name) == type(""): name = [name]

        for x in name:
            try:
                cf_remove(self.dbkey, x)
            except:
                if not ignoreError: raise            

    def get_freq(self, name):
        """Finds the frequency of the object stored in the db as `name`"""
        if not self.dbIsOpen:
            raise DBError("Database is not open")

        result = cf_size(self.dbkey, name.upper())
        return result['freq']


    def whats(self, name):
        """Preforms a fame "whats" command on the provided series"""
        if type(name) == type(""): name = [name]

        result = {}
        for dbname in name:
            if not self.dbIsOpen:
                raise DBError("Database is not open")

            result[dbname] = cf_whats(self.dbkey, dbname.upper())

        if len(result) == 1:
            return result.values()[0]
        return result



    def restore(self):
        """Discard any changes made to the database since it was last opened or posted."""
        return cf_restore(self.dbkey)


class cFameCall:
    """wrapper for cfame functions that acquires and releases a resource log.
This is needed because the Fame C api is not thread safe."""

    def __init__ (self, func):
        self.f = func
        self.__doc__ = getattr(func, "__doc__", str(func))
        self.__name__ = getattr(func, "__name__", str(func))
    #
    def __call__ (self, *args, **kwargs):
        "Execute the call behavior."
        tmp = fameLock.acquire()
        try:
            result = self.f(*args, **kwargs)
            fameLock.release()
        except:
            fameLock.release()
            raise
            
        return result

cf_open = cFameCall(cfame.open)
cf_set_option = cFameCall(cfame.set_option)
cf_close = cFameCall(cfame.close)
cf_restore = cFameCall(cfame.restore)
cf_size = cFameCall(cfame.size)
cf_whats = cFameCall(cfame.whats)
cf_remove = cFameCall(cfame.remove)
cf_create = cFameCall(cfame.create)
cf_read = cFameCall(cfame.read)
cf_write_scalar = cFameCall(cfame.write_scalar)
cf_write_series = cFameCall(cfame.write_series)
cf_write_namelist = cFameCall(cfame.write_namelist)
cf_wildlist = cFameCall(cfame.wildlist)
cf_exists = cFameCall(cfame.exists)

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
