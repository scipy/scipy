import pylab as pl
import timeseries as ts
from matplotlib.ticker import FixedLocator, FuncFormatter
import numpy as np
import types


class DateAxis:
    def __init__(self, start_date, end_date):

        freq = start_date.freq
        numPers = end_date - start_date + 1
        
        self.xaxis = np.arange(numPers)
        self.xlabel_txt = ''
        
        tser = ts.tser(start_date, end_date)
        
        minor_labels = ['' for x in range(numPers)]
        major_labels = ['' for x in range(numPers)]
        
        minor_ticks = []
        major_ticks = []
        
        years = ts.year(tser)
        year_starts = (years - ts.year(tser-1) != 0)
        year_start_inds = [int(x) for x in np.where(year_starts)[0]]
        
        quarters = ts.quarter(tser)
        quarter_starts = (quarters - ts.quarter(tser-1) != 0)
        quarter_start_inds = [int(x) for x in np.where(quarter_starts)[0]]
        
        months = ts.month(tser)
        month_starts = (months - ts.month(tser-1) != 0)
        month_start_inds = [int(x) for x in np.where(month_starts)[0]]

        
        def yearsOnly(minor_labels, major_labels, minor_ticks, major_ticks):

            numYears = end_date.year() - start_date.year()

            if numYears < 20: minor_spacing, major_spacing = 1, 2
            elif numYears < 50: minor_spacing, major_spacing = 1, 5
            elif numYears < 100: minor_spacing, major_spacing = 5, 10
            elif numYears < 200: minor_spacing, major_spacing = 5, 20
            elif numYears < 400: minor_spacing,  major_spacing = 5, 25
            elif numYears < 1000: minor_spacing,  major_spacing = 10, 50
            else: minor_spacing,  major_spacing = 20, 100

            for x in [y for y in year_start_inds if years[y] % minor_spacing == 0]:
                minor_ticks += [x]
                if years[x] % major_spacing == 0:
                    major_ticks += [x]
                    major_labels[x] += (start_date + x).strfmt('%Y')
                    
            return minor_labels, major_labels, minor_ticks, major_ticks
        
        if freq == 'A':
            minor_labels, major_labels, minor_ticks, major_ticks = \
                yearsOnly(minor_labels, major_labels, minor_ticks, major_ticks)

        elif freq == 'Q':

            if numPers <= 3 * 4:
            
                minor_ticks = list(range(numPers))
                major_ticks = year_start_inds

                for x in range(numPers):
                    minor_labels[x] = (start_date + x).strfmt('Q%q')
                    if year_starts[x]:
                        major_labels[x] += '\n\n'+(start_date + x).strfmt('%Y')
    
            elif numPers <= 11 * 4:

                minor_ticks = list(range(numPers))

                for x in [y for y in quarter_start_inds if quarters[y] == 1]:
                    major_ticks += [x]
                    major_labels[x] += (start_date + x).strfmt('%Y')
                
            else:
                
                minor_labels, major_labels, minor_ticks, major_ticks = \
                    yearsOnly(minor_labels, major_labels, minor_ticks, major_ticks)

        elif freq == 'M':

            if numPers <= 1.5 * 12:

                minor_ticks = list(range(numPers))
                major_ticks = year_start_inds

                for x in range(numPers):
                    minor_labels[x] = (start_date + x).strfmt('%b')
                    if year_starts[x]:
                        major_labels[x] += '\n\n'+(start_date + x).strfmt('%Y')

            elif numPers <= 3 * 12:
                
                minor_ticks = list(range(numPers))

                for x in range(numPers):

                    _date = start_date + x
                    if (_date.month() - 1) % 2 == 0:
                        minor_labels[x] += _date.strfmt('%b')

                        if year_starts[x] == 1:
                            major_ticks += [x]
                            major_labels[x] += '\n\n'+_date.strfmt('%Y')
    
            elif numPers <= 11 * 12:

                minor_ticks = quarter_start_inds

                for x in [y for y in quarter_start_inds if quarters[y] == 1]:
                    major_ticks += [x]
                    major_labels[x] += (start_date + x).strfmt('%Y')
                
            else:
                
                minor_labels, major_labels, minor_ticks, major_ticks = \
                    yearsOnly(minor_labels, major_labels, minor_ticks, major_ticks)

        elif freq in ('B', 'D'):
            
            if freq == 'B': daysInYear = 261
            else: daysInYear = 365

            if numPers <= daysInYear // 3:

                if numPers < daysInYear // 12: spacing = 1
                else: spacing = 5
                
                minor_ticks = list(range(numPers))
                major_ticks = [int(x) for x in month_start_inds]
                
                for x in range(numPers):
                
                    _date = start_date + x
                
                    if x % spacing == 0:  minor_labels[x] = _date.strfmt('%d')
                    if month_starts[x]:
                        major_labels[x] += '\n\n'+_date.strfmt('%b')
                        major_ticks += [x]
                        
                        if year_starts[x]:
                            major_labels[x] += '\n'+_date.strfmt('%Y')

            elif numPers <= 1.5 * daysInYear:
                
                minor_ticks = list(range(numPers))
                major_ticks = month_start_inds

                for x in month_start_inds:
                    _date = start_date + x
                    major_labels[x] += _date.strfmt('%b')
                    if months[x] == 1: major_labels[x] += '\n'+_date.strfmt('%Y')
                        
            elif numPers <= 3 * daysInYear:
                
                minor_ticks = month_start_inds

                for x in [y for y in month_start_inds if (months[y] - 1) % 3 == 0]:

                    _date = start_date + x
                    minor_labels[x] += _date.strfmt('%b')

                    if months[x] == 1:
                        major_ticks += [x]
                        major_labels[x] += '\n\n'+_date.strfmt('%Y')

            elif numPers <= 11 * daysInYear:
                
                minor_ticks = quarter_start_inds

                for x in [y for y in quarter_start_inds if quarters[y] == 1]:
                    major_ticks += [x]
                    major_labels[x] += (start_date + x).strfmt('%Y')
                
            else:

                minor_labels, major_labels, minor_ticks, major_ticks = \
                    yearsOnly(minor_labels, major_labels, minor_ticks, major_ticks)

        else:
            raise ValueError("unsupported frequency: "+freq)
        
        # if date range is such that no year or month is included in the labels, add these to the start
        if years[0] == years[-1] and (start_date - 1).year() == start_date.year():
            
            breaksBeforeYear = 2
        
            if freq in ('B', 'D'):
                if months[0] == months[-1] and (start_date - 1).month() == start_date.month():
                    major_labels[0] += '\n\n'+start_date.strfmt('%b')
                    breaksBeforeYear = 1

            if breaksBeforeYear > 1 and len([x for x in major_labels if x != '']) > 0: breaksBeforeYear += 1
            
            major_labels[0] += ('\n' * breaksBeforeYear)+(start_date).strfmt('%Y')
            
            if 0 not in major_ticks: major_ticks = [0] + major_ticks
            
        self.major_labels = major_labels
        self.major_locator = FixedLocator(major_ticks)
        
        self.minor_labels = minor_labels
        self.minor_locator = FixedLocator(minor_ticks)
       
    def minor_formatter(self, x, pos):
        if type(x) is types.IntType:
            return self.minor_labels[x]
        else:
            return ''

    def major_formatter(self, x, pos):
        if type(x) is types.IntType:
            return self.major_labels[x]
        else:
            return ''


    def set_labels(self, subplot):
        subplot.xaxis.set_minor_locator(self.minor_locator)
        subplot.xaxis.set_minor_formatter(FuncFormatter(self.minor_formatter))

        subplot.xaxis.set_major_locator(self.major_locator)
        subplot.xaxis.set_major_formatter(FuncFormatter(self.major_formatter))
