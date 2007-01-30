import pylab as pl
import timeseries as ts
import timeseries.plotlib as tpl
import numpy as np

numPoints = 100
freq = 'b'

y_ts1 = ts.TimeSeries(100*np.cumprod(1 + np.random.normal(size=numPoints)/100), freq=freq, start_date=ts.thisday(freq)-numPoints)
y_ts2 = ts.TimeSeries(100*np.cumprod(1 + np.random.normal(size=numPoints)/100), freq=freq, start_date=ts.thisday(freq)-numPoints)

y_ts1, y_ts2 = ts.aligned(y_ts1, y_ts2)

sDate, eDate = y_ts1.start_date(), y_ts1.end_date()
dAxis = tpl.DateAxis(sDate, eDate)

pl.clf()

a = pl.subplot(1,1,1)
a.plot(dAxis.xaxis, y_ts1, '-', dAxis.xaxis, y_ts2, '-')
dAxis.set_labels(a)

#a.set_xlim(dAxis.xaxis[0], dAxis.xaxis[-1])

pl.title('Time Series Plot')
pl.show()


