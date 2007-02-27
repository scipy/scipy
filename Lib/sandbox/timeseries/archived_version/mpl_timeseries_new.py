"""
Classes to plot TimeSeries w/ matplotlib.

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknow_ca_at_hotmail_dot_com
:date: $Date: 2007-02-02 23:19:06 -0500 (Fri, 02 Feb 2007) $
:version: $Id: tdates.py 2726 2007-02-19 07:37:26Z pierregm $
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author: pierregm $)"
__version__ = '1.0'
__revision__ = "$Revision: 2676 $"
__date__     = '$Date: 2007-02-02 23:19:06 -0500 (Fri, 02 Feb 2007) $'


import matplotlib
from matplotlib import pylab, rcParams
from matplotlib.artist import setp
from matplotlib.axes import Subplot, PolarSubplot
from matplotlib.cbook import flatten
from matplotlib.collections import LineCollection
from matplotlib.contour import ContourSet
from matplotlib.dates import DayLocator, MonthLocator, YearLocator, \
                             DateFormatter
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.mlab import meshgrid
from matplotlib.ticker import Formatter, ScalarFormatter, FuncFormatter, \
                              Locator, FixedLocator

#from matplotlib.transforms import nonsingular

import numpy as N
import maskedarray as MA

import timeseries
from timeseries import date_array, Date, DateArray, TimeSeries
#from tdates import date_array, Date
#import tseries
#from tseries import TimeSeries

import warnings

#####---------------------------------------------------------------------------
#---- --- Matplotlib extensions ---
#####---------------------------------------------------------------------------

def add_generic_subplot(figure_instance, *args, **kwargs):
    """Generalizes the `add_subplot` figure method to generic subplots.
The specific Subplot object class to add is given through the keywords
`SubplotClass` or `class`.

:Parameters:
    `figure_instance` : Figure object
        Figure to which the generic subplot should be attached.
    `args` : Misc
        Miscellaneous arguments to the subplot.
    `kwargs` : Dictionary
        Keywords. Same keywords as `Subplot`, with the addition of
        - `SubplotClass` : Type of subplot
        - `subclass` : Shortcut to `SubplotClass`.
        - any keyword required by the `SubplotClass` subclass.
    """

    key = figure_instance._make_key(*args, **kwargs)
    #TODO: Find why, sometimes, key is not hashable (even if tuple)
    # else, there's a fix below
    try:
        key.__hash__()
    except TypeError:
        key = str(key)
    #        
    if figure_instance._seen.has_key(key):
        ax = figure_instance._seen[key]
        figure_instance.sca(ax)
        return ax
    #
    if not len(args): 
        return
#    if hasattr(args[0], '__array__'):
#        fixedargs = args[1:]
#    else:
#        fixedargs = args
    #
    SubplotClass = kwargs.pop("SubplotClass", Subplot)
    SubplotClass = kwargs.pop("subclass",SubplotClass)
    if isinstance(args[0], Subplot) or isinstance(args[0], PolarSubplot):
        a = args[0]
        assert(a.get_figure() is figure_instance)
#        a.set_figure(figure_instance)
    else:
        ispolar = kwargs.pop('polar', False)
        if ispolar:
            a = PolarSubplot(figure_instance, *args, **kwargs)
        else:
            a = SubplotClass(figure_instance, *args, **kwargs)
            
    figure_instance.axes.append(a)
    figure_instance._axstack.push(a)
    figure_instance.sca(a)
    figure_instance._seen[key] = a
    return a


def nonsingular(vmin, vmax, expander=0.001, tiny=1e-15, increasing=True):
    '''
    Ensure the endpoints of a range are not too close together.

    "too close" means the interval is smaller than 'tiny' times
            the maximum absolute value.

    If they are too close, each will be moved by the 'expander'.
    If 'increasing' is True and vmin > vmax, they will be swapped.
    '''
    #TODO: Remove that when matplotlib incorporate it by default
    swapped = False
    if vmax < vmin:
        vmin, vmax = vmax, vmin
        swapped = True
    if vmax - vmin <= max(abs(vmin), abs(vmax)) * tiny:
        if vmin == 0.0:
            vmin = -expander
            vmax = expander
        else:
            vmin -= expander*abs(vmin)
            vmax += expander*abs(vmax)
    if swapped and not increasing:
        vmin, vmax = vmax, vmin
    return vmin, vmax

##### -------------------------------------------------------------------------
#---- --- Locators ---
##### -------------------------------------------------------------------------

def _year_start(_dates):
    return (_dates.year != (_dates-1).year)

def _quarter_start(_dates):
    return (_dates.quarter != (_dates-1).quarter)
    
def _month_start(_dates):
    return (_dates.month != (_dates-1).month)

def _week_start(_dates):
    return (_dates.day_of_week == 1)

def _get_default_annual_spacing(nyears):
    """Returns a default spacing between consecutive ticks for annual data."""
    if nyears < 11: 
        (min_spacing, maj_spacing) = (1, 1)
    elif nyears < 20: 
        (min_spacing, maj_spacing) = (1, 2)
    elif nyears < 50: 
        (min_spacing, maj_spacing) = (1, 5)
    elif nyears < 100: 
        (min_spacing, maj_spacing) = (5, 10)
    elif nyears < 200: 
        (min_spacing, maj_spacing) = (5, 25)
    elif nyears < 600: 
        (min_spacing, maj_spacing) = (10, 50)
    else:
        factor = nyears // 1000 + 1
        (min_spacing, maj_spacing) = (factor*20, factor*100)
    return (min_spacing, maj_spacing)


def _BreakDown_ParamCheck(locator, formatter):
    if not locator and not formatter:
        raise ValueError("Must specify either locator or formatter")
    
    if locator and formatter:
        raise ValueError("Must specify only one of locator or formatter")

def _Daily_BreakDown(dates, locator=False, formatter=False):
    
    _BreakDown_ParamCheck(locator, formatter)
    
    if dates.freqstr == 'B': periodsperyear = 261
    elif dates.freqstr == 'D': periodsperyear = 365
    else: raise ValueError("unexpected frequency")

    vmin = dates[0].value
    vmax = dates[-1].value
    span = vmax - vmin + 1

    if locator:
        default = N.arange(vmin, vmax+1) 
    else: #formatter
        format = N.empty(dates.size, dtype="|S8")
        format.flat = ''

    if span <= (periodsperyear//12 - 2):
    
        month_start = _month_start(dates)
    
        if locator:
            major = default[month_start]
            minor = default
        else:
            year_start = _year_start(dates)
            year_start[0] = False
            month_start[0] = False
            
            format[:] = '%d'
            format[month_start] = '%d\n%b'
            format[year_start] = '%d\n%b\n%Y'

            if not year_start.any():
                if not month_start.any():
                    if dates.size > 1: idx = 1
                    else: idx = 0
                    format[idx] = '%d\n%b\n%Y'
                else:
                    format[N.where(month_start)[0][0]] = '%d\n%b\n%Y'
            
    elif span <= periodsperyear//4:
    
        month_start = _month_start(dates)

        if locator:
            major = default[month_start]
            minor = default
        else:
            week_start = _week_start(dates)
            year_start = _year_start(dates)
            
            week_start[0] = False
            month_start[0] = False
            year_start[0] = False

            format[week_start] = '%d'
            format[month_start] = '\n\n%b'
            format[year_start] = '\n\n%b\n%Y'

            if not year_start.any():
                if not month_start.any():
                    format[N.where(week_start)[0][0]] = '\n\n%b\n%Y'
                else:
                    format[N.where(month_start)[0][0]] = '\n\n%b\n%Y'
        
    elif span <= 1.15 * periodsperyear:
        month_start = _month_start(dates)

        if locator:
            week_start = _week_start(dates)
            minor_idx = week_start | month_start
            minor_idx[0] = True
            major = default[month_start]
            minor = default[minor_idx]
        else:
            year_start = _year_start(dates)
            month_start[0] = False
            year_start[0] = False
            format[month_start] = '%b'
            format[year_start] = '%b\n%Y'

            if not year_start.any():
                format[N.where(month_start)[0][0]] = '%b\n%Y'
        
    elif span <= 2.5 * periodsperyear:

        year_start = _year_start(dates)
        month_start = _month_start(dates)

        if locator:
            major = default[year_start]
            minor = default[month_start]
        else:
            quarter_start = _quarter_start(dates)
            format[quarter_start] = '%b'
            format[year_start] = '%b\n%Y'

    elif span <= 4 * periodsperyear:

        year_start = _year_start(dates)
        month_start = _month_start(dates)

        if locator:
            major = default[year_start]
            minor = default[month_start]
        else:
            jan = (dates.month == 1)
            jul = (dates.month == 7)
            jan_or_jul = month_start & (jan | jul)
            format[jan_or_jul] = '%b'
            format[year_start] = '%b\n%Y'

    elif span <= 11 * periodsperyear:
    
        year_start = _year_start(dates)

        if locator:
            quarter_start = _quarter_start(dates)
            major = default[year_start]
            minor = default[quarter_start]
        else:
            format[year_start] = '%Y'
        
    else:
    
        (min_anndef, maj_anndef) = _get_default_annual_spacing(span/periodsperyear)
        year_start = _year_start(dates)

        major_idx = year_start & (dates.years % maj_anndef == 0)

        if locator:
            major = default[major_idx]
            minor = default[year_start & (dates.years % min_anndef == 0)]
        else:
            format[major_idx] = '%Y'

    if locator:
        return minor, major
    else:
        return format


def _Monthly_BreakDown(dates, locator=False, formatter=False):

    _BreakDown_ParamCheck(locator, formatter)
    
    if dates.freqstr != 'M': raise ValueError("unexpected frequency")

    periodsperyear = 12

    vmin = dates[0].value
    vmax = dates[-1].value
    span = vmax - vmin + 1

    if locator:
        default = N.arange(vmin, vmax+1) 
    else: #formatter
        format = N.empty(dates.size, dtype="|S8")
        format.flat = ''

    if span <= 1.15 * periodsperyear:
        year_start = _year_start(dates)

        if locator:
            major = default[year_start]
            minor = default
        else:
            year_start[0] = False

            format[:] = '%b'
            format[year_start] = '%b\n%Y'

            if not year_start.any():
                if dates.size > 1: idx = 1
                else: idx = 0
                format[idx] = '%b\n%Y'
        
    elif span <= 2.5 * periodsperyear:

        year_start = _year_start(dates)

        if locator:
            major = default[year_start]
            minor = default
        else:
            quarter_start = _quarter_start(dates)
            format[quarter_start] = '%b'
            format[year_start] = '%b\n%Y'

    elif span <= 4 * periodsperyear:

        year_start = _year_start(dates)

        if locator:
            major = default[year_start]
            minor = default
        else:
            months = dates.month
            format[(months == 1) | (months == 7)] = '%b'
            format[year_start] = '%b\n%Y'

    elif span <= 11 * periodsperyear:
    
        year_start = _year_start(dates)

        if locator:
            quarter_start = _quarter_start(dates)
            major = default[year_start]
            minor = default[quarter_start]
        else:
            format[year_start] = '%Y'
        
    else:
    
        (min_anndef, maj_anndef) = _get_default_annual_spacing(span/periodsperyear)
        year_start = _year_start(dates)

        major_idx = year_start & (dates.years % maj_anndef == 0)

        if locator:
            major = default[major_idx]
            minor = default[year_start & (dates.years % min_anndef == 0)]
        else:
            format[major_idx] = '%Y'

    if locator:
        return minor, major
    else:
        return format


def _Quarterly_BreakDown(dates, locator=False, formatter=False):

    _BreakDown_ParamCheck(locator, formatter)
    
    if dates.freqstr != 'Q': raise ValueError("unexpected frequency")

    periodsperyear = 4

    vmin = dates[0].value
    vmax = dates[-1].value
    span = vmax - vmin + 1

    if locator:
        default = N.arange(vmin, vmax+1) 
    else: #formatter
        format = N.empty(dates.size, dtype="|S8")
        format.flat = ''

    if span <= 3.5 * periodsperyear:

        year_start = _year_start(dates)

        if locator:
            major = default[year_start]
            minor = default
        else:
            year_start[0] = False
            
            format[:] = 'Q%q'
            format[year_start] = 'Q%q\n%Y'

            if not year_start.any():
                if dates.size > 1: idx = 1
                else: idx = 0
                format[idx] = 'Q%q\n%Y'

    elif span <= 11 * periodsperyear:
    
        year_start = _year_start(dates)

        if locator:
            major = default[year_start]
            minor = default
        else:
            format[year_start] = '%Y'
        
    else:
    
        (min_anndef, maj_anndef) = _get_default_annual_spacing(span/periodsperyear)
        year_start = _year_start(dates)

        major_idx = year_start & (dates.years % maj_anndef == 0)

        if locator:
            major = default[major_idx]
            minor = default[year_start & (dates.years % min_anndef == 0)]
        else:
            format[major_idx] = '%Y'

    if locator:
        return minor, major
    else:
        return format


def _Annual_BreakDown(dates, locator=False, formatter=False):

    _BreakDown_ParamCheck(locator, formatter)
    
    if dates.freqstr != 'A': raise ValueError("unexpected frequency")

    vmin = dates[0].value
    vmax = dates[-1].value
    span = vmax - vmin + 1

    if locator:
        default = N.arange(vmin, vmax+1) 
    else: #formatter
        format = N.empty(dates.size, dtype="|S8")
        format.flat = ''

    (min_anndef, maj_anndef) = _get_default_annual_spacing(span)
    year_start = _year_start(dates)

    major_idx = year_start & (dates.years % maj_anndef == 0)

    if locator:
        major = default[major_idx]
        minor = default[year_start & (dates.years % min_anndef == 0)]
    else:
        format[major_idx] = '%Y'

    if locator:
        return minor, major
    else:
        return format

#...............................................................................
class TimeSeries_DateLocator(Locator):
    "Locates the ticks along an axis controlled by a DateArray."

    def __init__(self, freq, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        self.freqstr = freq
        self.base = base
        (self.quarter, self.month, self.day) = (quarter, month, day)
        self.isminor = minor_locator
        self.isdynamic = dynamic_mode
        self.offset = 0
            
    def _initialize_dates(self, start_val, end_val):
        "Returns a DateArray for the current frequency."
        freq = self.freqstr
        dates = date_array(start_date=Date(freq, value=int(start_val)),
                           end_date=Date(freq, value=int(end_val)), 
                           freq=freq)
        return dates

    def _get_default_spacing(self, span):
        "Returns the default ticks spacing."
        raise NotImplementedError('Derived must override')
    
    def _get_default_locs(self, vmin, vmax):
        "Returns the default ticks spacing."
        raise NotImplementedError('Derived must override')
    
    def __call__(self):
        'Return the locations of the ticks.'
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if self.isdynamic:
            locs = self._get_default_locs(vmin, vmax)
        else:
            base = self.base
            (d, m) = divmod(vmin, base)
            vmin = (d+1) * base
            locs = range(vmin, vmax+1, base)
        return locs
    
    def autoscale(self):
        """Sets the view limits to the nearest multiples of base that contain 
    the data.
        """
        self.verify_intervals()
        dmin, dmax = self.dataInterval.get_bounds()
        locs = self._get_default_locs(dmin, dmax)
        (vmin, vmax) = locs[[0, -1]]
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        return nonsingular(vmin, vmax)      

def _generic_get_default_locs(self, vmin, vmax, BreakDownFunc):
        dates = self._initialize_dates(vmin, vmax)
        minor, major = BreakDownFunc(dates, locator=True)
        if self.isminor: return minor
        return major

#...............................................................................
class TimeSeries_AnnualLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by an annual DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'A', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
    
    def _get_default_locs(self, vmin, vmax):
        "Returns the default tick spacing for annual data."
        return _generic_get_default_locs(self, vmin, vmax, _Annual_BreakDown)
    
#...............................................................................
class TimeSeries_QuarterlyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a quarterly DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'Q', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        self.offset=1
    
    def _get_default_locs(self, vmin, vmax):
        "Returns the default ticks spacing."
        return _generic_get_default_locs(self, vmin, vmax, _Quarterly_BreakDown)
    
#...............................................................................
class TimeSeries_MonthlyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a monthly DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'M', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        self.offset = 1
    
    def _get_default_locs(self, vmin, vmax):
        "Returns the default ticks spacing."
        return _generic_get_default_locs(self, vmin, vmax, _Monthly_BreakDown)

#...............................................................................
class TimeSeries_DailyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a daily DateArray."

    def __init__(self, freq, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self, freq, minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        self._cacheddates = None
        
    def _get_default_locs(self, vmin, vmax):
        "Returns the default tick locations for daily data."
        return _generic_get_default_locs(self, vmin, vmax, _Daily_BreakDown)

#...............................................................................
class TimeSeries_YearLocator(TimeSeries_DateLocator):
    """Locates ticks along a Date axis, for each (multiple of) year.
    
:Ivariables:
    - `base` : Integer
      Gives the spacing between two consecutive annual ticks.
    - `quarter` : Integer *[1]*
      Tells on which quarter the ticks should be.
    - `month` : Integer *[1]*
      Tells on which month the ticks should be.
    - `day` : Integer *[1]*
      Tells on which day the ticks should be.    
    """
    def __init__(self, freq, minor_locator=False,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self, freq, minor_locator, False,
                                        base, quarter, month, day)
    
    def __call__(self):
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        freq = self.freqstr
        if freq == 'A':
            return range(vmin, vmax+1, self.base)
        else:
            dates = self._initialize_dates()
            if freq == 'Q':
                locs = (dates.quarters == self.quarter)
            elif freq == 'M':
                locs = (dates.months == self.month)
            elif freq in 'BDU':
                locs = (dates.months == self.month) & (dates.day == self.day)
            if self.base > 1:
                locs &= (locs.cumsum() % self.base == 1)
            return dates.tovalue()[locs]
#...............................................
class TimeSeries_QuarterLocator(TimeSeries_DateLocator):
    """Locates ticks along a Date axis, for each (multiple of) quarter.
    
:Ivariables:
    - `base` : Integer
      Gives the spacing between two consecutive quarter ticks.
    - `month` : Integer *[1]*
      Tells on which month the ticks should be.
    - `day` : Integer *[1]*
      Tells on which day the ticks should be.    
    """
    
    def __init__(self, freq, minor_locator=False,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self, freq, minor_locator, False,
                                        base, quarter, month, day)
    
    def __call__(self):
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        freq = self.freqstr
        if freq == 'A':
            msg = "The current frequency ('%s') is too coarse!" % freq
            raise ValueError, msg
        elif freq == 'Q':
            return range(vmin, vmax+1, self.base)
        else:
            dates = self._initialize_dates()
            values = dates.tovalue()
            if freq == 'M':
                locs = (dates.months % 4 == self.month)
            elif freq in 'BDU':
                locs = (dates.months % 4 == self.month) & (dates.day == self.day)
            if self.base > 1:
                locs &= (locs.cumsum() % self.base == 1)
            return values[locs]
#...............................................................................
class TimeSeries_MonthLocator(TimeSeries_DateLocator):
    """Locates ticks along a Date axis, for each (multiple of) month.
    
:Ivariables:
    - `base` : Integer
      Gives the spacing between two consecutive quarter ticks.
    - `month` : Integer *[1]*
      Tells on which month the ticks should be.
    - `day` : Integer *[1]*
      Tells on which day the ticks should be.    
    """
    
    def __init__(self, freq, minor_locator=False,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self, freq, minor_locator, False,
                                        base, quarter, month, day)
    
    def __call__(self):
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        freq = self.freqstr
        if freq == 'AQ':
            msg = "The current frequency ('%s') is too coarse!" % freq
            raise ValueError, msg
        elif freq == 'M':
            return range(vmin, vmax+1, self.base)
        else:
            dates = self._initialize_dates()
            values = dates.tovalue()
            if freq in 'BDU':
                locs = (dates.months == self.month) & (dates.day == self.day)
            if self.base > 1:
                locs &= (locs.cumsum() % self.base == 1)
            return values[locs]         

#####---------------------------------------------------------------------------
#---- --- Formatter ---
#####---------------------------------------------------------------------------            
class TimeSeries_DateFormatter(Formatter):
    """Formats the ticks along a DateArray axis."""
    
    def __init__(self, freq, minor_locator=False, dynamic_mode=True,):
        self.format = None
        self.freqstr = freq
        self.locs = []
        self.formatdict = {}
        self.isminor = minor_locator
        self.isdynamic = dynamic_mode
        self.offset = 0

            
    def _initialize_dates(self, locs):
        "Returns a DateArray for the current frequency."
        freq = self.freqstr
        dates = date_array(dlist=self.locs, freq=freq)
        return dates
    
    def set_locs(self, locs):
        'Sets the locations of the ticks'
        self.locs = locs
        if len(self.locs) > 0:
            self.verify_intervals()
            d = abs(self.viewInterval.span())
            self._set_format(d)
    #
    def __call__(self, x, pos=0):
        if self.isminor:
            fmt = self.formatdict.pop(x, '')
            if fmt is not '':
                retval = Date(self.freqstr, value=int(x)).strfmt(fmt)
            else:
                retval = ''
        else:
            retval = ''
        return retval
    
def _generic_set_format(self, BreakDownFunc):
    dates = self._initialize_dates(self.locs)
    format = BreakDownFunc(dates, formatter=True)
    return dict([(x,f) for (x,f) in zip(self.locs, format)])

#...............................................................................    
class TimeSeries_AnnualFormatter(TimeSeries_DateFormatter):
    #
    def __init__(self, minor_locator=False, dynamic_mode=True,):
        TimeSeries_DateFormatter.__init__(self, 'A',
                                          minor_locator=minor_locator, 
                                          dynamic_mode=dynamic_mode,)

    def _set_format(self, span):
        self.formatdict = _generic_set_format(self, _Annual_BreakDown)
        
#...............................................................................        
class TimeSeries_QuarterlyFormatter(TimeSeries_DateFormatter):
    #
    def __init__(self, minor_locator=False, dynamic_mode=True,):
        TimeSeries_DateFormatter.__init__(self, 'Q',
                                          minor_locator=minor_locator, 
                                          dynamic_mode=dynamic_mode,)

    def _set_format(self, span):
        self.formatdict = _generic_set_format(self, _Quarterly_BreakDown)
        
#...............................................................................        
class TimeSeries_MonthlyFormatter(TimeSeries_DateFormatter):
    #
    def __init__(self, minor_locator=False, dynamic_mode=True,):
        TimeSeries_DateFormatter.__init__(self, 'M',
                                          minor_locator=minor_locator, 
                                          dynamic_mode=dynamic_mode,)
    #
    def _set_format(self, span):
        self.formatdict = _generic_set_format(self, _Monthly_BreakDown)

#...............................................................................
class TimeSeries_DailyFormatter(TimeSeries_DateFormatter):
    #
    def __init__(self, freq, minor_locator=False, dynamic_mode=True,):
        TimeSeries_DateFormatter.__init__(self, freq, 
                                          minor_locator=minor_locator, 
                                          dynamic_mode=dynamic_mode,)
    #
    def _set_format(self, span):
        self.formatdict = _generic_set_format(self, _Daily_BreakDown)

#####--------------------------------------------------------------------------
#---- --- TimeSeries plots ---
#####--------------------------------------------------------------------------
class TimeSeriesPlot(Subplot, object):
    """Defines a time series based subclass of Subplot."""
    def __init__(self, fig=None, *args, **kwargs):
        """
Accepts the same keywords as a standard subplot, plus a specific `series` keyword.

:Parameters:
    `fig` : Figure
        Base figure.
        
:Keywords:
    `series` : TimeSeries
        Data to plot
        
        """
        # Retrieve the series ...................
        _series = kwargs.pop('series',None)
        Subplot.__init__(self,fig,*args,**kwargs)
#        # Force fig to be defined .....
#        if fig is None:
#            fig = TSFigure(_series)
        # Process options .......................
        if _series is not None:
            assert hasattr(_series, "dates")
            self._series = _series.ravel()
            self.xdata = _series.dates
            self.freqstr = _series.dates.freqstr
            self.xaxis.set_major_locator
            
        else:
            self._series = None
            self.xdata = None
            self.freqstr = None
        self._austoscale = False
        # Get the data to plot 
        self.legendsymbols = []
        self.legendlabels = []
    #............................................
    def set_ydata(self, series=None):
        """Sets the base time series."""
        if self._series is not None:
            print "WARNING ! Base series is being changed."""
        self._series = series.ravel()
        if isinstance(series, TimeSeries):
            self.xdata = self.series.dates
    #....
    def get_ydata(self):
        """Gets the base time series."""
        return self._series
    ydata = property(fget=get_ydata, fset=set_ydata, doc='Time series')
    #............................................    
    def _check_plot_params(self,*args):
        """Defines the plot coordinates (and basic plotting arguments)."""
        remaining = list(args)
        # No args ? Use defaults, if any
        if len(args) == 0:
            if self.xdata is None:
                raise ValueError, "No date information available!"
            return (self.xdata, self.ydata)
        output = []
        while len(remaining) > 0:
            a = remaining.pop(0)
            # The argument is a format: use default dates/
            if isinstance(a,str):
                if self.xdata is None:
                    raise ValueError, "No date information available!"
                else:
                    output.extend([self.xdata, self.ydata, a])
            # The argument is a TimeSeries: use its dates for x
            elif isinstance(a, TimeSeries):
                (x,y) = (a._dates, a._series)
                if len(remaining) > 0 and isinstance(remaining[0], str):
                    b = remaining.pop(0)
                    output.extend([x,y,b])
                else:
                    output.extend([x,y])
            # The argument is a DateArray............
            elif isinstance(a, (Date, DateArray)):
                # Force to current freq
                if self.freqstr is not None:
                    if a.freqstr != self.freqstr:
                        a = a.asfreq(self.freqstr)
                # There's an argument after
                if len(remaining) > 0:
                    #...and it's a format string
                    if isinstance(remaining[0], str):
                        b = remaining.pop(0)
                        if self.ydata is None:
                            raise ValueError, "No data information available!"
                        else:
                            output.extend([a, self.ydata, b])
                    #... and it's another date: use the default
                    elif isinstance(remaining[0], DateArray):
                        if self.ydata is None:
                            raise ValueError, "No data information available!"
                        else:
                            output.extend([a, self.ydata])
                    #... and it must be some data
                    else:   
                        b = remaining.pop(0)
                        if len(remaining) > 0:
                            if isinstance(remaining[0], str):
                                c = remaining.pop(0)
                                output.extend([a,b,c])
                            else:
                                output.extend([a,b])
                     #   continue
                else:
                    if self.ydata is None:
                        raise ValueError, "No data information available!"
                    #else:
                    #    break
            # Otherwise..............................
            elif len(remaining) > 0:
                if isinstance(remaining[0], str):
                    b = remaining.pop(0)
                    if self.xdata is None:
                        raise ValueError, "No date information available!"
                    else:
                        output.extend([self.xdata, a, b])
                    #continue
                elif self.xdata is None:
                    raise ValueError, "No date information available!"
                else:
                    output.extend([self.xdata, a])
                    #continue
        # Reinitialize the plot if needed ...........
        if self.xdata is None:
            self.xdata = output[0]
            self.freqstr = self.xdata.freqstr
        # Force the xdata to the current frequency
        elif output[0].freqstr != self.freqstr:
            output = list(output)
            output[0] = output[0].asfreq(self.freqstr)
        return output
    #............................................
    def tsplot(self,*parms,**kwargs):
        """Plots the data parsed in argument.
This command accepts the same keywords as `matplotlib.plot`."""
        #print "Parameters: %s - %i" % (parms, len(parms))
        parms = self._check_plot_params(*parms)
        self.legendlabels.append(kwargs.get('label',None))
        Subplot.plot(self, *parms,**kwargs)
        pylab.draw_if_interactive()
#    #............................................
#    def ybaseline(self,ybase,**kwargs):
#        """Plots an horizontal baseline on each subplot."""
#        self.axhline(ybase,**kwargs)
    #............................................       
    def format_dateaxis(self,maj_spacing=None, min_spacing=None, 
                        strformat="%Y", rotate=True):
        """Pretty-formats the date axis (x-axis).
        
:Parameters:
    `major` : Integer *[5]* 
        Major tick locator, in years (major tick every `major` years).
    `minor` : Integer *[12]* 
        Minor tick locator, in months (minor ticks every `minor` months).
    `strformat` : String *['%Y']*
        String format for major ticks ("%Y").
        """
        # Get the locator class .................
        if self.freqstr in 'BDU':
            locator = TimeSeries_DailyLocator
            formatter = TimeSeries_DailyFormatter
            self.xaxis.set_major_locator(locator(self.freqstr,
                                                 minor_locator=False,
                                                 dynamic_mode=True))
            self.xaxis.set_minor_locator(locator(self.freqstr,
                                                 minor_locator=True,
                                                 dynamic_mode=True))
            self.xaxis.set_major_formatter(formatter(self.freqstr,
                                                     minor_locator=False,
                                                     dynamic_mode=True))
            self.xaxis.set_minor_formatter(formatter(self.freqstr,
                                                     minor_locator=True,
                                                     dynamic_mode=True))
        else:
            if self.freqstr == 'A':
                locator = TimeSeries_AnnualLocator
                formatter = TimeSeries_AnnualFormatter
            elif self.freqstr == 'Q':
                locator = TimeSeries_QuarterlyLocator
                formatter = TimeSeries_QuarterlyFormatter
            elif self.freqstr == 'M':
                locator = TimeSeries_MonthlyLocator
                formatter = TimeSeries_MonthlyFormatter
            self.xaxis.set_major_locator(locator(minor_locator=False,
                                                 dynamic_mode=True))
            self.xaxis.set_minor_locator(locator(minor_locator=True,
                                                 dynamic_mode=True))
            self.xaxis.set_major_formatter(formatter(minor_locator=False,
                                                     dynamic_mode=True))
            self.xaxis.set_minor_formatter(formatter(minor_locator=True,
                                                     dynamic_mode=True))
        #........................................
#        if rcParams['backend'] == 'PS':
#            rotate = False
#            warnings.warn("dateplot: PS backend detected, rotate disabled")
#        if self.is_last_row():
#            if rotate:
#                setp(self.get_xticklabels(),rotation=45)
#        self.xaxis.set_major_formatter(FuncFormatter(self.dateticks_formatter))
#        self.xaxis.set_minor_formatter(FuncFormatter(self.dateticks_formatter))
#        else:
#            self.set_xticklabels([])
#            self.set_xlabel('')          

TSPlot = TimeSeriesPlot


#####--------------------------------------------------------------------------
#---- --- TimeSeries Figures ---
#####--------------------------------------------------------------------------        
class TimeSeriesFigure(Figure):
    """Time Series Figure: all subplots share the same time series.
    """
    def __init__(self, series=None, **kwargs):
        self._series = series
        Figure.__init__(self,**kwargs)
        fspnum = kwargs.pop('fspnum',None)
        if fspnum is not None:
            self.add_tsplot(fspnum, series=series)
    #.........
    def add_tsplot(self, *args, **kwargs):
        """Adds a `TimeSeriesPlot` subplot to the figure."""
        kwargs.update(SubplotClass=TimeSeriesPlot,
                      series=self._series)        
        return add_generic_subplot(self, *args, **kwargs)
    add_plot = add_tsplot
TSFigure = TimeSeriesFigure
#Figure.add_tsplot = 
#................................................
def tsfigure(series, **figargs):    
    """Creates a new `TimeSeriesFigure` object.
    
:Parameters:
    `series` : TimeSeries object
        Input data.
    `figargs` : Dictionary
        Figure options [`figsize`, `dpi`, `facecolor`, `edgecolor`, `frameon`].
    """
    figargs.update(FigureClass=TSFigure)
    figargs.update(series=series)
#    print "figargs:",figargs
#    num = figargs.pop('num',None)
    fig = pylab.figure(**figargs)
    return fig

def add_tsplot(axes, *args, **kwargs):
    kwargs.update(SubplotClass=TimeSeriesPlot)
    if 'series' not in kwargs.keys():
        kwargs['series'] = None
    return add_generic_subplot(axes, *args, **kwargs)
Figure.add_tsplot = add_tsplot
    

def tsplot(*args, **kwargs):
    # allow callers to override the hold state by passing hold=True|False
    b = pylab.ishold()
    h = kwargs.pop('hold', None)
    if h is not None:
        pylab.hold(h)
    try:
        ret =  pylab.gca().add_tsplot(*args, **kwargs)
        pylab.draw_if_interactive()
    except:
        pylab.hold(b)
        raise

    pylab.hold(b)
    return ret

################################################################################
if __name__ == '__main__':

    da = date_array(start_date=Date(freq='A', year=2003, quarter=3, month=1, day=17), 
                    length=51)
    ser = timeseries.time_series(MA.arange(len(da)), dates=da)
    ser[4] = MA.masked
#    ser_2 = timeseries.time_series(MA.arange(len(da)), dates=da.asfreq('Q'))
    
    pylab.figure()
    pylab.gcf().add_tsplot(111)
    pylab.gca().tsplot(ser, 'ko-')
    pylab.gca().format_dateaxis()
#    pylab.gca().tsplot(ser_2, 'rs')
    pylab.show()
    