import scipy.interpolate as interpolate
import numpy as np

def scoreatpercentile(data, per):
    """Return the score at the given 'per' percentile of the data.
    
    Example
    >>> scoreatpercentile(randn(100), 50)
    will return the median of the sample.
    """
    cdf = empiricalcdf(data)
    interpolator = interpolate.interp1d(sort(cdf), sort(data))
    return interpolator(per/100.)
    
def percentileofscore(data, score):
    """Percentile-position of score relative to data.
      
    score: Array of scores at which the percentile is computed.
    
    Return percentiles (0-100).
    
    Example
        x = linspace(-2,2,100)
        percentileofscore(randn(50),x)
        
    Return an error if the score is outside the range of data. 
    """
    cdf = empiricalcdf(data)
    interpolator = interpolate.interp1d(sort(data), sort(cdf))
    return interpolator(score)*100.

def empiricalcdf(data, method='Hazen'):
    """Return the empirical cdf.
    
    Methods available:
    	Hazen:       (i-0.5)/N
	    Weibull:     i/(N+1)
    	Chegodayev:  (i-.3)/(N+.4)
    	Cunnane:     (i-.4)/(N+.2)
    	Gringorten:  (i-.44)/(N+.12)
    	California:  (i-1)/N
    
    Where i goes from 1 to N. 
    """
    
    i = np.argsort(np.argsort(data)) + 1.
    N = len(data)
    method = method.lower()
    
    if method == 'weibull':
        cdf = i/(N+1.)
    elif method == 'hazen':
        cdf = (i-0.5)/N
    elif method == 'california':
        cdf = (i-1.)/N
    elif method == 'chegodayev':
        cdf = (i-.3)/(N+.4)
    elif method == 'cunnane':
        cdf = (i-.4)/(N+.2)
    elif method == 'gringorten':
        cdf = (i-.44)/(N+.12)
    else:
        raise 'Unknown method. Choose among Weibull, Hazen, Chegodayev, Cunnane, Gringorten and California.'
    
    return cdf 
