""" Test functions for stats module

"""

from numpy.testing import *


import numpy
from numpy import typecodes, array
import scipy.stats as stats

def kolmogorov_check(diststr,args=(),N=20,significance=0.01):
    qtest = stats.ksoneisf(significance,N)
    cdf = eval('stats.'+diststr+'.cdf')
    dist = eval('stats.'+diststr)
    # Get random numbers
    kwds = {'size':N}
    vals = numpy.sort(dist.rvs(*args,**kwds))
    cdfvals = cdf(vals,*args)
    q = max(abs(cdfvals-arange(1.0,N+1)/N))
    assert (q < qtest), "Failed q=%f, bound=%f, alpha=%f" % (q, qtest, significance)
    return


# generate test cases to test cdf and distribution consistency
dists = ['uniform','norm','lognorm','expon','beta',
         'powerlaw','bradford','burr','fisk','cauchy','halfcauchy',
         'foldcauchy','gamma','gengamma','loggamma',
         'alpha','anglit','arcsine','betaprime','erlang',
         'dgamma','exponweib','exponpow','frechet_l','frechet_r',
         'gilbrat','f','ncf','chi2','chi','nakagami','genpareto',
         'genextreme','genhalflogistic','pareto','lomax','halfnorm',
         'halflogistic','fatiguelife','foldnorm','ncx2','t','nct',
         'weibull_min','weibull_max','dweibull','maxwell','rayleigh',
         'genlogistic', 'logistic','gumbel_l','gumbel_r','gompertz',
         'hypsecant', 'laplace', 'reciprocal','triang','tukeylambda']

# check function for test generator
def check_distribution(dist, args, alpha):
    D,pval = stats.kstest(dist,'', args=args, N=1000)
    if (pval < alpha):
        D,pval = stats.kstest(dist,'',args=args, N=1000)
        #if (pval < alpha):
        #    D,pval = stats.kstest(dist,'',args=args, N=1000)
        assert (pval > alpha), "D = " + str(D) + "; pval = " + str(pval) + \
               "; alpha = " + str(alpha) + "\nargs = " + str(args)

# nose test generator
def test_all_distributions():
    for dist in dists:
        distfunc = getattr(stats, dist)
        nargs = distfunc.numargs
        alpha = 0.01
        if dist == 'fatiguelife':
            alpha = 0.001
        if dist == 'erlang':
            args = (4,)+tuple(rand(2))
        elif dist == 'frechet':
            args = tuple(2*rand(1))+(0,)+tuple(2*rand(2))
        elif dist == 'triang':
            args = tuple(rand(nargs))
        elif dist == 'reciprocal':
            vals = rand(nargs)
            vals[1] = vals[0] + 1.0
            args = tuple(vals)
        else:
            args = tuple(1.0+rand(nargs))
        yield check_distribution, dist, args, alpha


class TestRandInt(TestCase):
    def test_rvs(self):
        vals = stats.randint.rvs(5,30,size=100)
        assert(numpy.all(vals < 30) & numpy.all(vals >= 5))
        assert(len(vals) == 100)
        vals = stats.randint.rvs(5,30,size=(2,50))
        assert(numpy.shape(vals) == (2,50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.randint.rvs(15,46)
        assert((val >= 15) & (val < 46))
        assert isinstance(val, numpy.ScalarType),`type(val)`
        assert(val.dtype.char in typecodes['AllInteger'])

    def test_pdf(self):
        k = numpy.r_[0:36]
        out = numpy.where((k >= 5) & (k < 30), 1.0/(30-5), 0)
        vals = stats.randint.pmf(k,5,30)
        assert_array_almost_equal(vals,out)

    def test_cdf(self):
        x = numpy.r_[0:36:100j]
        k = numpy.floor(x)
        out = numpy.select([k>=30,k>=5],[1.0,(k-5.0+1)/(30-5.0)],0)
        vals = stats.randint.cdf(x,5,30)
        assert_array_almost_equal(vals, out, decimal=12)

class TestBinom(TestCase):
    def test_rvs(self):
        vals = stats.binom.rvs(10, 0.75, size=(2, 50))
        assert(numpy.all(vals >= 0) & numpy.all(vals <= 10))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.binom.rvs(10, 0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])


class TestBernoulli(TestCase):
    def test_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0) & numpy.all(vals <= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestNBinom(TestCase):
    def test_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestGeom(TestCase):
    def test_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

    def test_pmf(self):
        vals = stats.geom.pmf([1,2,3],0.5)
        assert_array_almost_equal(vals,[0.5,0.25,0.125])

    def test_cdf_sf(self):
        vals = stats.geom.cdf([1,2,3],0.5)
        vals_sf = stats.geom.sf([1,2,3],0.5)
        expected = array([0.5,0.75,0.875])
        assert_array_almost_equal(vals,expected)
        assert_array_almost_equal(vals_sf,1-expected)


class TestHypergeom(TestCase):
    def test_rvs(self):
        vals = stats.hypergeom.rvs(20, 10, 3, size=(2, 50))
        assert(numpy.all(vals >= 0) &
               numpy.all(vals <= 3))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.hypergeom.rvs(20, 3, 10)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestLogser(TestCase):
    def test_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestPoisson(TestCase):
    def test_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestZipf(TestCase):
    def test_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestDLaplace(TestCase):
    def test_rvs(self):
        vals = stats.dlaplace.rvs(1.5 , size=(2, 50))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestRvDiscrete(TestCase):
    def test_rvs(self):
        states = [-1,0,1,2,3,4]
        probability = [0.0,0.3,0.4,0.0,0.3,0.0]
        samples = 1000
        r = stats.rv_discrete(name='sample',values=(states,probability))
        x = r.rvs(size=samples)

        for s,p in zip(states,probability):
            assert abs(sum(x == s)/float(samples) - p) < 0.05

class TestExpon(TestCase):
    def test_zero(self):
        assert_equal(stats.expon.pdf(0),1)

if __name__ == "__main__":
    run_module_suite()
