import numpy as N
import maskedarray as MA
import tseries as TS
import tdates as TD
D = TD.Date(freq='D', year=2007, month=1, day=1)
M = TD.Date(freq='M', year=2007, month=1, day=1)
Y = TD.Date(freq='A', year=2007, month=1, day=1)
TD.Date(freq='Q',year=2004,quarter=3)
TD.Date(freq='D',year=2001,month=1,day=1)
TD.Date('D', '2007-01-01')
TD.Date('D', mxDate=mx.DateTime.now())
TD.Date('D', mxDate=datetime.datetime.now())
data = N.random.uniform(-100,100,600)
today = TD.thisday('B')
series = TS.time_series(data, dtype=np.float_, freq='B', observed='SUMMED',
                        start_date=today-600)
series[0]
series[-30:]
thirtydaysago = today - 30
series[thirtydaysago:]
series[thirtydaysago.tostring():]
series[series<0] = 0
series[series.day_of_week == 4] = 100
weekdays = td.day_of_week(series)
series[weekdays == 4] = 100
series_log = N.log(series)
mlist_1 = ['2005-%02i' % i for i in range(1,10)]
mlist_1 += ['2006-%02i' % i for i in range(2,13)]
mdata_1 = N.arange(len(mlist_1))
mser_1 = TS.time_series(mdata_1, mlist_1, observed='SUMMED')
mser = mser1.asfreq('M')
mser1.has_duplicated_dates()
mser1.has_missing_dates()
mlist_2 = ['2004-%02i' % i for i in range(1,13)]
mlist_2 += ['2005-%02i' % i for i in range(1,13)]
mser_2 = TS.time_series(N.arange(len(mlist_2)), mlist_2, observed='SUMMED')
mser_3 = mser_1 + mser_2
(malg_1,malg_2) = aligned(mser_1, mser_2) 
mser_1_filled = fill_missing_dates(mser_1)
(malg_1,malg_2) = align_series(mser_1_filled, mser_2) 
mser_3 = malg_1 + malg_2
mser_3 = filled(malg_1,0) + filled(malg_2,0)
(malg_1,malg_2) = aligned(mser_1_filled, mser2, 
                          start_date='2004-06', end_date='2006-06')
mseries = series.convert('M',func=ma.average)
mseries_default = series.convert('M')
mseries_2d = series.convert('M',func=None)
data = N.random.uniform(-100,100,100).astype(np.float_)
today = TS.today.asfreq('M') - 110
nseries = TS.TimeSeries(data, freq='m', observed='END',start_date=today)
sixtymonthsago = today-60
nseries[sixtymonthsago:sixtymonthsago+10] = 12
