from numpy.testing import *
from histogram import _histogram_fixed_binsize, _histogram_digitize,\
    _histogram_searchsort, histogram,_optimize_binning
import numpy as np
from numpy.random import rand

class test_histogram1d_functions(NumpyTestCase):
    def check_consistency(self):
        n = 100
        r = rand(n)*12-1
        bins = range(11)
        a = _histogram_fixed_binsize(r, bins[0], bins[1]-bins[0], len(bins)-1)
        b = _histogram_digitize(r, None, np.array(bins), False)
        c = _histogram_searchsort(r,bins)
        assert_array_equal(a,b)
        assert_array_equal(c,b)
        
class test_histogram(NumpyTestCase):
    def check_simple(self):
        n=100
        v=rand(n)
        (a,b)=histogram(v)
        #check if the sum of the bins equals the number of samples
        assert_equal(np.sum(a,axis=0),n)
        #check that the bin counts are evenly spaced when the data is from a linear function
        (a,b)=histogram(np.linspace(0,10,100))
        assert_array_equal(a,10)
        #Check the construction of the bin array
        a, b = histogram(v, bins=4, range=[.2,.8])
        assert_array_almost_equal(b['edges'],np.linspace(.2, .8, 5),8)
        #Check the number of outliers
        assert_equal((v<.2).sum(), b['lower'])
        assert_equal((v>.8).sum(),b['upper'])
        #Check the normalization
        bins = [0,.5,.75,1]
        a,b = histogram(v, bins, normed=True)
        assert_almost_equal((a*np.diff(bins)).sum(), 1)
        
    def check_axis(self):
        n,m = 100,20
        v = rand(n,m)
        a,b = histogram(v, bins=5)
        # Check dimension is reduced (axis=None).
        assert_equal(a.ndim, 1)
        #Check total number of count is equal to the number of samples.
        assert_equal(a.sum(), n*m)
        a,b = histogram(v, bins = 7, axis=0)
        # Check shape of new array is ok.
        assert(a.ndim == 2)
        assert_array_equal(a.shape,[7, m])
        # Check normalization is consistent 
        a,b = histogram(v, bins = 7, axis=0, normed=True)
        assert_array_almost_equal((a.T*np.diff(b['edges'])).sum(1), np.ones((m)),5)
        a,b = histogram(v, bins = 7, axis=1, normed=True)
        assert_array_equal(a.shape, [n,7])
        assert_array_almost_equal((a*np.diff(b['edges'])).sum(1), np.ones((n)))
        # Check results are consistent with 1d estimate
        a1, b1 = histogram(v[0,:], bins=b['edges'], normed=True)
        assert_array_almost_equal(a1, a[0,:],7)
            
    def check_weights(self):
        # Check weights = constant gives the same answer as no weights.
        v = rand(100)
        w = np.ones(100)*5
        a,b = histogram(v)
        na,nb = histogram(v, normed=True)
        wa,wb = histogram(v, weights=w)
        nwa,nwb = histogram(v, weights=w, normed=True)
        assert_array_equal(a*5, wa)
        assert_array_almost_equal(na, nwa,8)
        # Check weights are properly applied.
        v = np.linspace(0,10,10)
        w = np.concatenate((np.zeros(5), np.ones(5)))
        wa,wb = histogram(v, bins=np.linspace(0,10.01, 11),weights=w)
        assert_array_almost_equal(wa, w)
        
    def check_strategies(self):
        v = rand(100)
        ae,be = histogram(v, strategy='binsize')
        ab,bb = histogram(v, strategy='digitize')
        as,bs = histogram(v, strategy='searchsort')
        assert_array_equal(ae, ab)
        assert_array_equal(ae, as)
        
        w = rand(100)
        ae,be = histogram(v, weights=w, strategy='binsize')
        ab,bb = histogram(v, weights=w, strategy='digitize')
        as,bs = histogram(v, weights=w, strategy='searchsort')
        assert_array_almost_equal(ae, ab,8)
        assert_array_almost_equal(ae, as,8)
    
    def check_automatic_binning(self):
        v = rand(100)
        h,b = histogram(v, 'Scott')
        h,b = histogram(v, 'Freedman')
                            
        
if __name__ == "__main__":
    NumpyTest().run()
