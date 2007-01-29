"""
Classes to plot TimeSeries w/ matplotlib.

:author: Pierre GF Gerard-Marchant
:contact: pierregm_at_uga_edu
:date: $Date$
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant ($Author$)"
__version__ = '1.0'
__revision__ = "$Revision$"
__date__     = '$Date$'


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

from matplotlib.transforms import nonsingular

import numpy as N
import maskedarray as MA

import tdates
from tdates import date_array, Date
import tseries
from tseries import TimeSeries

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

##### -------------------------------------------------------------------------
#---- --- Locators ---
##### -------------------------------------------------------------------------

def _get_default_annual_spacing(nyears):
    """Returns a default spacing between consecutive ticks for annual data."""
    if nyears < 20: 
        (min_spacing, maj_spacing) = (1, 2)
    elif nyears < 50: 
        (min_spacing, maj_spacing) = (1, 5)
    elif nyears < 100: 
        (min_spacing, maj_spacing) = (5, 10)
    elif nyears < 200: 
        (min_spacing, maj_spacing) = (5, 20)
    elif nyears < 400: 
        (min_spacing, maj_spacing) = (5, 25)
    elif nyears < 1000: 
        (min_spacing, maj_spacing) = (10, 50)
    else:
        (min_spacing, maj_spacing) = (20, 100)
    return (min_spacing, maj_spacing)

def _get_default_quarterly_spacing(nquarters):
    """Returns a default spacing between consecutive ticks for quarterly data."""
    if nquarters <= 3*4:
        (min_spacing, maj_spacing) = (1,4)
    elif nquarters <= 11*4:
        (min_spacing, maj_spacing) = (1,4)
    else:
        (min_anndef, maj_anndef) = _get_default_annual_spacing(nquarters//4)
        min_spacing = min_anndef * 4
        maj_spacing = maj_anndef * 4
    return (min_spacing, maj_spacing)

def _get_default_monthly_spacing(nmonths):
    """Returns a default spacing between consecutive ticks for monthly data."""
    if nmonths <= 10:
        (min_spacing, maj_spacing) = (1,3)
    elif nmonths <= 2*12:
        (min_spacing, maj_spacing) = (1,6)
    elif nmonths <= 3*12:
        (min_spacing, maj_spacing) = (1,12)
    elif nmonths <= 11*12:
        (min_spacing, maj_spacing) = (3,12)  
    else:
        (min_anndef, maj_anndef) = _get_default_annual_spacing(nmonths//12)
        min_spacing = min_anndef * 12
        maj_spacing = maj_anndef * 12
    return (min_spacing, maj_spacing)

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
        dates = date_array(start_date=Date(freq, value=start_val),
                           end_date=Date(freq, value=end_val), 
                           freq=freq)
        return dates

    def _get_default_spacing(self, span):
        "Returns the default ticks spacing."
        raise NotImplementedError('Derived must override')
    
    def __call__(self):
        'Return the locations of the ticks.'
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if self.isdynamic:
            base = self._get_default_spacing(vmax-vmin+1)
        else:
            base = self.base
        d = vmin // base
        vmin = (d+1) * base + self.offset
        locs = range(vmin, vmax+1, base)
        return locs
    
    def autoscale(self):
        """Sets the view limits to the nearest multiples of base that contain 
    the data.
        """
        self.verify_intervals()
        dmin, dmax = self.dataInterval.get_bounds()
        if self.isdynamic:
            base = self._get_default_spacing(dmax-dmin+1)
        else:
            base = self.base
        (d,m) = divmod(dmin, base)
        if m < base/2:
            vmin = d * base
        else:
            vmin = (d+1) * base
        (d,m) = divmod(dmax, base)
        vmax = (d+1) * base
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        return nonsingular(vmin, vmax)        
    
#...............................................................................
class TimeSeries_AnnualLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by an annual DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'A', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
    
    def _get_default_spacing(self, span):
        "Returns the default tick spacing for annual data."
        (minor, major) = _get_default_annual_spacing(span)
        if self.isminor:
            return minor
        return major
#...............................................................................
class TimeSeries_QuarterlyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a quarterly DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'Q', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        self.offset=1
    
    def _get_default_spacing(self, span):
        "Returns the default tick spacing for quarterly data."
        (minor, major) = _get_default_quarterly_spacing(span)
        if self.isminor:
            return minor
        return major       
#...............................................................................
class TimeSeries_MonthlyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a monthly DateArray."

    def __init__(self, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self,'M', minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        self.offset = 1
    
    def _get_default_spacing(self, span):
        "Returns the default tick spacing for monthly data."
        (minor, major) = _get_default_monthly_spacing(span)
        if self.isminor:
            return minor
        return major
    
#...............................................................................
class TimeSeries_DailyLocator(TimeSeries_DateLocator):
    "Locates the ticks along an axis controlled by a daily DateArray."

    def __init__(self, freq, minor_locator=False, dynamic_mode=True,
                 base=1, quarter=1, month=1, day=1):
        TimeSeries_DateLocator.__init__(self, freq, minor_locator, dynamic_mode,
                                        base, quarter, month, day)
        if self.freqstr == 'B':
            self.daysinyear = 261
        else:
            self.daysinyear = 365
        self._cacheddates = None
        
    def _get_default_locs(self, vmin, vmax):
        "Returns the default tick locations for daily data."
        daysperyear = self.daysinyear
        span = vmax - vmin + 1
        dates = self._initialize_dates(vmin, vmax)
        default = N.arange(vmin, vmax+1) 
        #
        if span <= daysperyear//12:
            minor = default
            major = default[(dates.day_of_week == 1)]
        elif span <= daysperyear//3:
            minor = default[(dates.day_of_week == 1)]
            major = default[(dates.day == 1)]
        elif span <= 1.5 * daysperyear:
            minor = default[(dates.day_of_week == 1)]
            major = default[(dates.day == 1)]
        elif span <= 3 * daysperyear:
            minor = default[(dates.day == 1)]
            major = default[(dates.day_of_year == 1)]
        elif span <= 11 * daysperyear:
            minor = default[(dates.quarter != (dates-1).quarter)]
            major = default[(dates.day_of_year == 1)]
        else:
            (min_anndef, maj_anndef) = _get_default_annual_spacing(span/daysperyear)
            annual = (dates.day_of_year == 1)
            minor = default[annual & (dates.years % min_anndef == 0)]
            major = default[annual & (dates.years % maj_anndef == 0)]
        if self.isminor:
            return minor
        return major

    def __call__(self):
        'Return the locations of the ticks'
        self.verify_intervals()
        vmin, vmax = self.viewInterval.get_bounds()
        if vmax < vmin:
            vmin, vmax = vmax, vmin
        if self.isdynamic:
            locs = self._get_default_locs(vmin, vmax)
        else:
            base = self.base
            (d,m) = divmod(vmin, base)
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
        (vmin, vmax) = locs[[0,-1]]
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        return nonsingular(vmin, vmax)        

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
    
    def __init__(self, freq, fmt=None):
        if fmt is None:
            fmt = Date.default_fmtstr[freq]
        self.fmt = fmt
        self.freqstr = freq
    
    def __call__(self, x, pos=0):
        return Date(self.freqstr, value=int(x)).strfmt(self.fmt)


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
        # At least three arguments ....
        if len(args) >= 3:
            params = args[:3]
        # Two arguments only ..........
        elif len(args) == 2:
            if isinstance(args[1], str):
                # The last argument is a format string
                arg = args[0]
                if isinstance(arg, TimeSeries):
                    params = (arg._dates, arg._series, args[1])
                elif self.xdata is None:
                    raise ValueError, "No date information available!"
                else:
                    params = (self.xdata, args[0], args[1])
            else:
                params = args
        # One argument only ...........
        elif len(args) == 1:
            if isinstance(args[0], str):
                if self.xdata is None:
                    raise ValueError, "No date information available!"
                else:
                    params =  (self.xdata, self.ydata, args[0])
            elif isinstance(args[0], TimeSeries):
                if self.xdata is None:
                    arg = args[0]
                    params = (arg._dates, arg._series)
                else:
                    params = (self.xdata, args[0])
        else:
            params = (self.xdata, self.ydata)
        # Reinitialize the plot if needed 
        if self.xdata is None:
            self.xdata = params[0]
            self.freqstr = self.xdata.freqstr
        # Force the xdata to the current frequency
        elif params[0].freqstr != self.freqstr:
            params = list(params)
            params[0] = params[0].asfreq(self.freqstr)
        return params
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
            self.xaxis.set_major_locator(locator(self.freqstr,
                                                 minor_locator=False,
                                                 dynamic_mode=True))
            self.xaxis.set_minor_locator(locator(self.freqstr,
                                                 minor_locator=True,
                                                 dynamic_mode=True))
        else:
            if self.freqstr == 'A':
                locator = TimeSeries_AnnualLocator
            elif self.freqstr == 'Q':
                locator = TimeSeries_QuarterlyLocator
            elif self.freqstr == 'M':
                locator = TimeSeries_MonthlyLocator
            self.xaxis.set_major_locator(locator(minor_locator=False,
                                                 dynamic_mode=True))
            self.xaxis.set_minor_locator(locator(minor_locator=True,
                                                 dynamic_mode=True))
        #........................................
        self.xaxis.set_major_formatter(TimeSeries_DateFormatter(self.freqstr))
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
#    #............................................
#    def plot_shifts(self,shifts,**kwargs):
#        """Plots regime shifts.
#:param shifts: Shifts/trends to plot.
#:type shifts: `RegimeShift`
#        """
#        self.tsplot(self.xdata,shifts.regimes,**kwargs)
#        for x in shifts.xshifts[0]:
#            self.axvline(self.xdata[x],ls=':',c='#999999',lw=0.5)    
    #............................................
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

    da = date_array(start_date=Date(freq='D', year=2003, quarter=3, month=1, day=17), 
                    length=51)
    ser = tseries.time_series(MA.arange(len(da)), dates=da)
    ser[4] = MA.masked
    ser_2 = tseries.time_series(MA.arange(len(da)), dates=da.asfreq('M'))
    
    pylab.figure()
    pylab.gcf().add_tsplot(111)
    pylab.gca().tsplot(ser, 'ko-')
    pylab.gca().format_dateaxis()
    pylab.gca().tsplot(ser_2, 'rs')
    pylab.show()
    