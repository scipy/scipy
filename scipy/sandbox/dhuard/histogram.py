import numpy as np
import subprocess

try:
	import flib
except:
	print 'Building the flib fortran library.'
	subprocess.call('f2py -c histogram.f -m flib', shell=True)
	import flib

def histogram(a, bins=10, range=None, normed=False, weights=None, axis=None, strategy=None):
    """histogram(a, bins=10, range=None, normed=False, weights=None, axis=None) 
                                                                   -> H, dict
    
    Return the distribution of sample.
    
    :Parameters:
      - `a` : Array sample.
      - `bins` : Number of bins, or an array of bin edges, in which case the 
                range is not used. If 'Scott' or 'Freeman' is passed, then 
                the named method is used to find the optimal number of bins.
      - `range` : Lower and upper bin edges, default: [min, max].
      - `normed` :Boolean, if False, return the number of samples in each bin,
                if True, return the density.  
      - `weights` : Sample weights. The weights are normed only if normed is 
                True. Should weights.sum() not equal len(a), the total bin count 
                will not be equal to the number of samples.
      - `axis` : Specifies the dimension along which the histogram is computed. 
                Defaults to None, which aggregates the entire sample array.
      - `strategy` : Histogramming method (binsize, searchsorted or digitize).
    
    :Return:
      - `H` : The number of samples in each bin. 
        If normed is True, H is a frequency distribution.
      - dict{ 'edges':      The bin edges, including the rightmost edge.
        'upper':      Upper outliers.
        'lower':      Lower outliers.
        'bincenters': Center of bins. 
        'strategy': the histogramming method employed.}
    
    :Examples:
      >>> x = random.rand(100,10)
      >>> H, D = histogram(x, bins=10, range=[0,1], normed=True)
      >>> H2, D = histogram(x, bins=10, range=[0,1], normed=True, axis=0)
    
    :SeeAlso: histogramnd
    """
    weighted = weights is not None
    
    a = np.asarray(a)
    if axis is None:
        a = np.atleast_1d(a.ravel())
        if weighted:
            weights = np.atleast_1d(weights.ravel())
        axis = 0 
        
    # Define the range    
    if range is None:
        mn, mx = a.min(), a.max()
        if mn == mx:
            mn = mn - .5
            mx = mx + .5
        range = [mn, mx]
    
    # Find the optimal number of bins.
    if bins is None or type(bins) == str:
        bins = _optimize_binning(a, range, bins)
        
    # Compute the bin edges if they are not given explicitely.    
    # For the rightmost bin, we want values equal to the right 
    # edge to be counted in the last bin, and not as an outlier. 
    # Hence, we shift the last bin by a tiny amount.
    if not np.iterable(bins):
        dr = np.diff(range)/bins*1e-10
        edges = np.linspace(range[0], range[1]+dr, bins+1, endpoint=True)
    else:
        edges = np.asarray(bins, float)
    
    dedges = np.diff(edges)
    bincenters = edges[:-1] + dedges/2.
    
    # Number of bins
    nbin = len(edges)-1
    
        # Measure of bin precision.
    decimal = int(-np.log10(dedges.min())+10)
        
    # Choose the fastest histogramming method
    even = (len(set(np.around(dedges, decimal))) == 1)
    if strategy is None:
        if even:
            strategy = 'binsize'
        else:
            if nbin > 30: # approximative threshold
                strategy = 'searchsort'
            else:
                strategy = 'digitize'
    else:
        if strategy not in ['binsize', 'digitize', 'searchsort']:
            raise 'Unknown histogramming strategy.', strategy
        if strategy == 'binsize' and not even:
            raise 'This binsize strategy cannot be used for uneven bins.'
        
    # Parameters for the fixed_binsize functions.
    start = float(edges[0])
    binwidth = float(dedges[0])
       
    # Looping to reduce memory usage
    block = 66600 
    slices = [slice(None)]*a.ndim
    for i in np.arange(0,len(a),block):
        slices[axis] = slice(i,i+block)
        at = a[slices]
        if weighted:
            at = np.concatenate((at, weights[slices]), axis)        
            if strategy == 'binsize':   
                count = np.apply_along_axis(_splitinmiddle,axis,at,
                    flib.weighted_fixed_binsize,start,binwidth,nbin)               
            elif strategy == 'searchsort':
                count = np.apply_along_axis(_splitinmiddle,axis,at, \
                        _histogram_searchsort_weighted, edges)
            elif strategy == 'digitize':
                    count = np.apply_along_axis(_splitinmiddle,axis,at,\
                        _histogram_digitize,edges,normed)
        else:
            if strategy == 'binsize':
                count = np.apply_along_axis(flib.fixed_binsize,axis,at,start,binwidth,nbin)
            elif strategy == 'searchsort':
                count = np.apply_along_axis(_histogram_searchsort,axis,at,edges)
            elif strategy == 'digitize':
                count = np.apply_along_axis(_histogram_digitize,axis,at,None,edges,
                        normed)
                    
        if i == 0:
            total = count
        else:
            total += count
        
    # Outlier count
    upper = total.take(np.array([-1]), axis)
    lower = total.take(np.array([0]), axis)
    
    # Non-outlier count
    core = a.ndim*[slice(None)]
    core[axis] = slice(1, -1)
    hist = total[core]
    
    if normed:
        normalize = lambda x: np.atleast_1d(x/(x*dedges).sum())
        hist = np.apply_along_axis(normalize, axis, hist)

    return hist, {'edges':edges, 'lower':lower, 'upper':upper, \
        'bincenters':bincenters, 'strategy':strategy}
        


def _histogram_fixed_binsize(a, start, width, n):
    """histogram_even(a, start, width, n) -> histogram
    
    Return an histogram where the first bin counts the number of lower
    outliers and the last bin the number of upper outliers. Works only with 
    fixed width bins. 
    
    :Parameters:
      a : array
        Array of samples.
      start : float
        Left-most bin edge.
      width : float
        Width of the bins. All bins are considered to have the same width.
      n : int
        Number of bins. 
    
    :Return:
      H : array
        Array containing the number of elements in each bin. H[0] is the number
        of samples smaller than start and H[-1] the number of samples 
        greater than start + n*width.
    """    
                 
    return flib.fixed_binsize(a, start, width, n)


def _histogram_binsize_weighted(a, w, start, width, n):
    """histogram_even_weighted(a, start, width, n) -> histogram
    
    Return an histogram where the first bin counts the number of lower
    outliers and the last bin the number of upper outliers. Works only with 
    fixed width bins. 
    
    :Parameters:
      a : array
        Array of samples.
      w : array
        Weights of samples.
      start : float
        Left-most bin edge.
      width : float
        Width of the bins. All bins are considered to have the same width.
      n : int
        Number of bins. 
    
    :Return:
      H : array
        Array containing the number of elements in each bin. H[0] is the number
        of samples smaller than start and H[-1] the number of samples 
        greater than start + n*width.
    """    
    return flib.weighted_fixed_binsize(a, w, start, width, n)
       
def _histogram_searchsort(a, bins):
    n = np.sort(a).searchsorted(bins)
    n = np.concatenate([n, [len(a)]])
    count = np.concatenate([[n[0]], n[1:]-n[:-1]])
    return count
    
def _histogram_searchsort_weighted(a, w, bins):
    i = np.sort(a).searchsorted(bins)
    sw = w[np.argsort(a)]
    i = np.concatenate([i, [len(a)]])
    n = np.concatenate([[0],sw.cumsum()])[i]
    count = np.concatenate([[n[0]], n[1:]-n[:-1]])
    return count

def _splitinmiddle(x, function, *args, **kwds):
    x1,x2 = np.hsplit(x, 2)
    return function(x1,x2,*args, **kwds)

def _histogram_digitize(a, w, edges, normed):
    """Internal routine to compute the 1d weighted histogram for uneven bins.
    a: sample
    w: weights
    edges: bin edges
    weighted: Means that the weights are appended to array a. 
    Return the bin count or frequency if normed.
    """
    weighted = w is not None
    nbin = edges.shape[0]+1
    if weighted:
        count = np.zeros(nbin, dtype=w.dtype)
        if normed:    
            count = np.zeros(nbin, dtype=float)
            w = w/w.mean()
    else:
        count = np.zeros(nbin, int)
            
    binindex = np.digitize(a, edges)
        
    # Count the number of identical indices.
    flatcount = np.bincount(binindex, w)
    
    # Place the count in the histogram array.
    count[:len(flatcount)] = flatcount
       
    return count


def _optimize_binning(x, range, method='Freedman'):
    """Find the optimal number of bins.
    Available methods : Freedman, Scott
    """
    N = x.shape[0]
    if method.lower()=='freedman':
        s=np.sort(x) 
        IQR = s[int(N*.75)] - s[int(N*.25)] # Interquantile range (75% -25%)
        width = 2* IQR*N**(-1./3)
        
    elif method.lower()=='scott':
        width = 3.49 * x.std()* N**(-1./3)
    else:
        raise 'Method must be Scott or Freedman', method
    return int(np.diff(range)/width)
