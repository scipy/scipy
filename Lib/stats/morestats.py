
from __future__ import nested_scopes
import distributions
import inspect
from scipy_base import isscalar
from scipy_base import zeros, arange, sort, amin, amax
import stats

__all__ = ['probplot','ppcc_max','ppcc_plot']

def probplot(x, sparams=(), dist='norm', fit=1, plot=None):
    """Return (osm, osr){,(scale,loc,r)} where (osm, osr) are order statistic
    medians and ordered response data respectively so that plot(osm, osr)
    is a probability plot.  If fit==1, then do a regression fit and compute the
    slope (scale), intercept (loc), and correlation coefficient (r), of the
    best straight line through the points.  If fit==0, only (osm, osr) is returned.
    
    sparams is a tuple of shape parameter arguments for the distribution.
    """
    N = len(x)
    Ui = zeros(N)*1.0
    Ui[-1] = 0.5**(1.0/N)
    Ui[0] = 1-Ui[-1]
    i = arange(2,N)
    Ui[1:-1] = (i-0.3175)/(N+0.365)
    try:
        ppf_func = eval('distributions.%sppf'%dist)
    except AttributError:
        raise dist, "is not a valid distribution with a ppf."
    if sparams is None:
        sparams = ()
    if isscalar(sparams):
        sparams = (sparams,)
    if not isinstance(sparams,types.TupleType):
        sparams = tuple(sparams)
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError, "Function has does not have default location", \
              "and scale parameters\n  that are 0.0 and 1.0 respectively."
    if (len(sparams) < len(res[0])-len(res[-1])-1) or \
       (len(sparams) > len(res[0])-3):
        raise ValueError, "Incorrect number of shape parameters."
    osm = ppf_func(Ui,*sparams)
    osr = sort(x)
    if fit or (plot is not None):
        # perform a linear fit.
        slope, intercept, r, prob, sterrest = stats.linregress(osm,osr)
    if plot is not None:
        try: xplt.limits()
        except: pass
        plot.plot(osm, osr, 'o', osm, slope*osm + intercept)
        plot.title('Probability Plot')
        plot.xlabel('Order Statistic Medians')
        plot.ylabel('Ordered Values')
        try: plot.expand_limits(5)
        except: pass
        xmin,xmax= amin(osm),amax(osm)
        ymin,ymax= amin(x),amax(x)
        pos = xmin+0.70*(xmax-xmin), ymin+0.01*(ymax-ymin)
        try: plot.addtext("r^2^=%1.4f" % r, xy=pos,tosys=1)
        except: pass
    if fit:
        return (osm, osr), (slope, intercept, r)
    else:
        return osm, osr

def ppcc_max(x, dist='tukeylambda'):
    """Returns the shape parameter that maximizes the probability plot
    correlation coefficient for the given data to a one-parameter
    family of distributions.

    See also ppcc_plot
    """
    try:
        ppf_func = eval('distributions.%sppf'%dist)
    except AttributError:
        raise dist, "is not a valid distribution with a ppf."
    res = inspect.getargspec(ppf_func)
    if not ('loc' == res[0][-2] and 'scale' == res[0][-1] and \
            0.0==res[-1][-2] and 1.0==res[-1][-1]):
        raise ValueError, "Function has does not have default location", \
              "and scale parameters\n  that are 0.0 and 1.0 respectively."
    if (1 < len(res[0])-len(res[-1])-1) or \
       (1 > len(res[0])-3):
        raise ValueError, "Must be a one-parameter family."
    N = len(x)
    Ui = zeros(N)*1.0
    Ui[-1] = 0.5**(1.0/N)
    Ui[0] = 1-Ui[-1]
    i = arange(2,N)
    Ui[1:-1] = (i-0.3175)/(N+0.365)
    osr = sort(x)
    def tempfunc(shape, mi, yvals, func):
        xvals = func(mi, shape)
        slope, intercept, r, prob, sterrest = stats.linregress(xvals, yvals)
        return 1-r
    return optimize.brent(tempfunc, args=(Ui, osr, ppf_func))

def ppcc_plot(x,a,b,dist='tukeylambda', plot=None, N=80):
    """Returns (shape, ppcc), and optionally plots shape vs. ppcc
    (probability plot correlation coefficient) as a function of shape
    parameter for a one-parameter family of distributions from shape
    value a to b.

    See also ppcc_max
    """
    svals = r_[a:b:complex(N)]
    ppcc = svals*0.0
    k=0
    for sval in svals:
        r1,r2 = probplot(x,sval,dist=dist,fit=1)
        ppcc[k] = r2[-1]
        k += 1
    if plot is not None:
        try: xplt.limits()
        except: pass
        plot.plot(svals, ppcc, 'x')
        plot.title('(%s) PPCC Plot' % dist)
        plot.xlabel('Prob Plot Corr. Coef.',deltay=-0.01)
        plot.ylabel('Shape Values',deltax=-0.01)
        try: plot.expand_limits(5)
        except: pass
    return svals, ppcc
        
    
