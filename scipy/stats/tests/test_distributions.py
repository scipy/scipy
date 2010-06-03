""" Test functions for stats module

"""

from numpy.testing import *


import numpy
import numpy as np
from numpy import typecodes, array
import scipy.stats as stats
from scipy.stats.distributions import argsreduce

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
         'hypsecant', 'laplace', 'reciprocal','triang','tukeylambda',
         'vonmises']

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
        elif dist == 'vonmises':
            yield check_distribution, dist, (10,), alpha
            yield check_distribution, dist, (101,), alpha
            args = tuple(1.0+rand(nargs))
        else:
            args = tuple(1.0+rand(nargs))
        yield check_distribution, dist, args, alpha

def check_vonmises_pdf_periodic(k,l,s,x):
    vm = stats.vonmises(k,loc=l,scale=s)
    assert_almost_equal(vm.pdf(x),vm.pdf(x%(2*numpy.pi*s)))
def check_vonmises_cdf_periodic(k,l,s,x):
    vm = stats.vonmises(k,loc=l,scale=s)
    assert_almost_equal(vm.cdf(x)%1,vm.cdf(x%(2*numpy.pi*s))%1)

def test_vonmises_pdf_periodic():
    for k in [0.1, 1, 101]:
        for x in [0,1,numpy.pi,10,100]:
            yield check_vonmises_pdf_periodic, k, 0, 1, x
            yield check_vonmises_pdf_periodic, k, 1, 1, x
            yield check_vonmises_pdf_periodic, k, 0, 10, x

            yield check_vonmises_cdf_periodic, k, 0, 1, x
            yield check_vonmises_cdf_periodic, k, 1, 1, x
            yield check_vonmises_cdf_periodic, k, 0, 10, x

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
        val = stats.randint(15,46).rvs(3)
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
        assert(isinstance(val, int))
        val = stats.binom(10, 0.75).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])


class TestBernoulli(TestCase):
    def test_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0) & numpy.all(vals <= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert(isinstance(val, int))
        val = stats.bernoulli(0.75).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestNBinom(TestCase):
    def test_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert(isinstance(val, int))
        val = stats.nbinom(10, 0.75).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestGeom(TestCase):
    def test_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert(isinstance(val, int))
        val = stats.geom(0.75).rvs(3)
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
        assert(isinstance(val, int))
        val = stats.hypergeom(20, 3, 10).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestLogser(TestCase):
    def test_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert(isinstance(val, int))
        val = stats.logser(0.75).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestPoisson(TestCase):
    def test_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert(numpy.all(vals >= 0))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert(isinstance(val, int))
        val = stats.poisson(0.5).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestZipf(TestCase):
    def test_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert(numpy.all(vals >= 1))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert(isinstance(val, int))
        val = stats.zipf(1.5).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestDLaplace(TestCase):
    def test_rvs(self):
        vals = stats.dlaplace.rvs(1.5 , size=(2, 50))
        assert(numpy.shape(vals) == (2, 50))
        assert(vals.dtype.char in typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert(isinstance(val, int))
        val = stats.dlaplace(1.5).rvs(3)
        assert(isinstance(val, numpy.ndarray))
        assert(val.dtype.char in typecodes['AllInteger'])

class TestRvDiscrete(TestCase):
    def test_rvs(self):
        states = [-1,0,1,2,3,4]
        probability = [0.0,0.3,0.4,0.0,0.3,0.0]
        samples = 1000
        r = stats.rv_discrete(name='sample',values=(states,probability))
        x = r.rvs(size=samples)
        assert(isinstance(x, numpy.ndarray))

        for s,p in zip(states,probability):
            assert abs(sum(x == s)/float(samples) - p) < 0.05

        x = r.rvs()
        assert(isinstance(x, int))

class TestExpon(TestCase):
    def test_zero(self):
        assert_equal(stats.expon.pdf(0),1)
        
    def test_tail(self):  # Regression test for ticket 807
        assert_equal(stats.expon.cdf(1e-18),  1e-18)
        assert_equal(stats.expon.isf(stats.expon.sf(40)),  40)

class TestGenExpon(TestCase):
    def test_pdf_unity_area(self):
        from scipy.integrate import simps
        # PDF should integrate to one
        assert_almost_equal(simps(stats.genexpon.pdf(numpy.arange(0,10,0.01),
                                                     0.5, 0.5, 2.0),
                                  dx=0.01), 1, 1)

    def test_cdf_bounds(self):
        # CDF should always be positive
        cdf = stats.genexpon.cdf(numpy.arange(0, 10, 0.01), 0.5, 0.5, 2.0)
        assert(numpy.all((0 <= cdf) & (cdf <= 1)))
        
class TestExponpow(TestCase):
    def test_tail(self):
        assert_almost_equal(stats.exponpow.cdf(1e-10,  2.),  1e-20)
        assert_almost_equal(stats.exponpow.isf(stats.exponpow.sf(5, .8), .8),  5)


class TestSkellam(TestCase):
    def test_pmf(self):
        #comparison to R
        k = numpy.arange(-10, 15)
        mu1, mu2 = 10, 5
        skpmfR = numpy.array(
                   [4.2254582961926893e-005, 1.1404838449648488e-004,
                    2.8979625801752660e-004, 6.9177078182101231e-004,
                    1.5480716105844708e-003, 3.2412274963433889e-003,
                    6.3373707175123292e-003, 1.1552351566696643e-002,
                    1.9606152375042644e-002, 3.0947164083410337e-002,
                    4.5401737566767360e-002, 6.1894328166820688e-002,
                    7.8424609500170578e-002, 9.2418812533573133e-002,
                    1.0139793148019728e-001, 1.0371927988298846e-001,
                    9.9076583077406091e-002, 8.8546660073089561e-002,
                    7.4187842052486810e-002, 5.8392772862200251e-002,
                    4.3268692953013159e-002, 3.0248159818374226e-002,
                    1.9991434305603021e-002, 1.2516877303301180e-002,
                    7.4389876226229707e-003])
        
        assert_almost_equal(stats.skellam.pmf(k, mu1, mu2), skpmfR, decimal=15)
        
    def test_cdf(self):
        #comparison to R, only 5 decimals
        k = numpy.arange(-10, 15)
        mu1, mu2 = 10, 5
        skcdfR = numpy.array(
                   [6.4061475386192104e-005, 1.7810985988267694e-004,
                    4.6790611790020336e-004, 1.1596768997212152e-003,
                    2.7077485103056847e-003, 5.9489760066490718e-003,
                    1.2286346724161398e-002, 2.3838698290858034e-002,
                    4.3444850665900668e-002, 7.4392014749310995e-002,
                    1.1979375231607835e-001, 1.8168808048289900e-001,
                    2.6011268998306952e-001, 3.5253150251664261e-001,
                    4.5392943399683988e-001, 5.5764871387982828e-001,
                    6.5672529695723436e-001, 7.4527195703032389e-001,
                    8.1945979908281064e-001, 8.7785257194501087e-001,
                    9.2112126489802404e-001, 9.5136942471639818e-001,
                    9.7136085902200120e-001, 9.8387773632530240e-001,
                    9.9131672394792536e-001])

        assert_almost_equal(stats.skellam.cdf(k, mu1, mu2), skcdfR, decimal=5)


class TestHypergeom(TestCase):
    def test_precision(self):
        # comparison number from mpmath

        M,n,N = 2500,50,500
        tot=M;good=n;bad=tot-good
        hgpmf = stats.hypergeom.pmf(2,tot,good,N)

        assert_almost_equal(hgpmf, 0.0010114963068932233, 11)

class TestChi2(TestCase):
    # regression tests after precision improvements, ticket:1041, not verified
    def test_precision(self):
        assert_almost_equal(stats.chi2.pdf(1000, 1000), 8.919133934753128e-003, 14)
        assert_almost_equal(stats.chi2.pdf(100, 100), 0.028162503162596778, 14)

class TestArrayArgument(TestCase): #test for ticket:992
    def test_noexception(self):
        rvs = stats.norm.rvs(loc=(np.arange(5)), scale=np.ones(5), size=(10,5))
        assert_equal(rvs.shape, (10,5))

class TestDocstring(TestCase):
    def test_docstrings(self):
        """See ticket #761"""
        if stats.rayleigh.__doc__ is not None:
            self.failUnless("rayleigh" in stats.rayleigh.__doc__.lower())
        if stats.bernoulli.__doc__ is not None:
            self.failUnless("bernoulli" in stats.bernoulli.__doc__.lower())

class TestEntropy(TestCase):
    def test_entropy_positive(self):
        """See ticket #497"""
        pk = [0.5,0.2,0.3]
        qk = [0.1,0.25,0.65]
        eself = stats.entropy(pk,pk)
        edouble = stats.entropy(pk,qk)
        assert(0.0 == eself)
        assert(edouble >= 0.0)

def TestArgsreduce():
    a = array([1,3,2,1,2,3,3])
    b,c = argsreduce(a > 1, a, 2)

    assert_array_equal(b, [3,2,2,3,3])
    assert_array_equal(c, [2,2,2,2,2])

    b,c = argsreduce(2 > 1, a, 2)
    assert_array_equal(b, a[0])
    assert_array_equal(c, [2])

    b,c = argsreduce(a > 0, a, 2)
    assert_array_equal(b, a)
    assert_array_equal(c, [2] * numpy.size(a))


if __name__ == "__main__":
    run_module_suite()
