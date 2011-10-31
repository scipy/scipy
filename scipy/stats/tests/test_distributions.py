""" Test functions for stats module

"""

from numpy.testing import TestCase, run_module_suite, assert_equal, \
    assert_array_equal, assert_almost_equal, assert_array_almost_equal, \
    assert_allclose, assert_, rand, dec


import numpy
import numpy as np
from numpy import typecodes, array
import scipy.stats as stats
from scipy.stats.distributions import argsreduce

def kolmogorov_check(diststr, args=(), N=20, significance=0.01):
    qtest = stats.ksoneisf(significance, N)
    cdf = eval('stats.'+diststr+'.cdf')
    dist = eval('stats.'+diststr)
    # Get random numbers
    kwds = {'size':N}
    vals = numpy.sort(dist.rvs(*args, **kwds))
    cdfvals = cdf(vals, *args)
    q = max(abs(cdfvals - np.arange(1.0, N+1)/N))
    assert_(q < qtest, msg="Failed q=%f, bound=%f, alpha=%f" % (q, qtest, significance))
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
        assert_(pval > alpha, msg="D = " + str(D) + "; pval = " + str(pval) + \
               "; alpha = " + str(alpha) + "\nargs = " + str(args))

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
        assert_(numpy.all(vals < 30) & numpy.all(vals >= 5))
        assert_(len(vals) == 100)
        vals = stats.randint.rvs(5,30,size=(2,50))
        assert_(numpy.shape(vals) == (2,50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.randint.rvs(15,46)
        assert_((val >= 15) & (val < 46))
        assert_(isinstance(val, numpy.ScalarType), msg=`type(val)`)
        val = stats.randint(15,46).rvs(3)
        assert_(val.dtype.char in typecodes['AllInteger'])

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
        assert_(numpy.all(vals >= 0) & numpy.all(vals <= 10))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.binom.rvs(10, 0.75)
        assert_(isinstance(val, int))
        val = stats.binom(10, 0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])


class TestBernoulli(TestCase):
    def test_rvs(self):
        vals = stats.bernoulli.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0) & numpy.all(vals <= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.bernoulli.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.bernoulli(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

class TestNBinom(TestCase):
    def test_rvs(self):
        vals = stats.nbinom.rvs(10, 0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.nbinom.rvs(10, 0.75)
        assert_(isinstance(val, int))
        val = stats.nbinom(10, 0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

class TestGeom(TestCase):
    def test_rvs(self):
        vals = stats.geom.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.geom.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.geom(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

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
        assert_(numpy.all(vals >= 0) &
               numpy.all(vals <= 3))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.hypergeom.rvs(20, 3, 10)
        assert_(isinstance(val, int))
        val = stats.hypergeom(20, 3, 10).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

    def test_precision(self):
        # comparison number from mpmath
        M = 2500
        n = 50
        N = 500
        tot = M
        good = n
        hgpmf = stats.hypergeom.pmf(2, tot, good, N)
        assert_almost_equal(hgpmf, 0.0010114963068932233, 11)

    def test_precision2(self):
        """Test hypergeom precision for large numbers.  See #1218."""
        # Results compared with those from R.
        oranges = 9.9e4
        pears = 1.1e5
        fruits_eaten = np.array([3, 3.8, 3.9, 4, 4.1, 4.2, 5]) * 1e4
        quantile = 2e4
        res = []
        for eaten in fruits_eaten:
            res.append(stats.hypergeom.sf(quantile, oranges + pears, oranges, eaten))
        expected = np.array([0, 1.904153e-114, 2.752693e-66, 4.931217e-32,
                             8.265601e-11, 0.1237904, 1])
        assert_allclose(res, expected, atol=0, rtol=5e-7)

        # Test with array_like first argument
        quantiles = [1.9e4, 2e4, 2.1e4, 2.15e4]
        res2 = stats.hypergeom.sf(quantiles, oranges + pears, oranges, 4.2e4)
        expected2 = [1, 0.1237904, 6.511452e-34, 3.277667e-69]
        assert_allclose(res2, expected2, atol=0, rtol=5e-7)


class TestLogser(TestCase):
    def test_rvs(self):
        vals = stats.logser.rvs(0.75, size=(2, 50))
        assert_(numpy.all(vals >= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.logser.rvs(0.75)
        assert_(isinstance(val, int))
        val = stats.logser(0.75).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

class TestPoisson(TestCase):
    def test_rvs(self):
        vals = stats.poisson.rvs(0.5, size=(2, 50))
        assert_(numpy.all(vals >= 0))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.poisson.rvs(0.5)
        assert_(isinstance(val, int))
        val = stats.poisson(0.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

class TestZipf(TestCase):
    def test_rvs(self):
        vals = stats.zipf.rvs(1.5, size=(2, 50))
        assert_(numpy.all(vals >= 1))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.zipf.rvs(1.5)
        assert_(isinstance(val, int))
        val = stats.zipf(1.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

class TestDLaplace(TestCase):
    def test_rvs(self):
        vals = stats.dlaplace.rvs(1.5 , size=(2, 50))
        assert_(numpy.shape(vals) == (2, 50))
        assert_(vals.dtype.char in typecodes['AllInteger'])
        val = stats.dlaplace.rvs(1.5)
        assert_(isinstance(val, int))
        val = stats.dlaplace(1.5).rvs(3)
        assert_(isinstance(val, numpy.ndarray))
        assert_(val.dtype.char in typecodes['AllInteger'])

def test_rvgeneric_std():
    """Regression test for #1191"""
    assert_array_almost_equal(stats.t.std([5, 6]), [1.29099445, 1.22474487])

class TestRvDiscrete(TestCase):
    def test_rvs(self):
        states = [-1,0,1,2,3,4]
        probability = [0.0,0.3,0.4,0.0,0.3,0.0]
        samples = 1000
        r = stats.rv_discrete(name='sample',values=(states,probability))
        x = r.rvs(size=samples)
        assert_(isinstance(x, numpy.ndarray))

        for s,p in zip(states,probability):
            assert_(abs(sum(x == s)/float(samples) - p) < 0.05)

        x = r.rvs()
        assert_(isinstance(x, int))

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
        assert_(numpy.all((0 <= cdf) & (cdf <= 1)))

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


class TestGamma(TestCase):

    def test_pdf(self):
        # a few test cases to compare with R
        pdf = stats.gamma.pdf(90, 394, scale=1./5)
        assert_almost_equal(pdf, 0.002312341)

        pdf = stats.gamma.pdf(3, 10, scale=1./5)
        assert_almost_equal(pdf, 0.1620358)


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
            self.assertTrue("rayleigh" in stats.rayleigh.__doc__.lower())
        if stats.bernoulli.__doc__ is not None:
            self.assertTrue("bernoulli" in stats.bernoulli.__doc__.lower())

    def test_no_name_arg(self):
        """If name is not given, construction shouldn't fail.  See #1508."""
        stats.rv_continuous()
        stats.rv_discrete()


class TestEntropy(TestCase):
    def test_entropy_positive(self):
        """See ticket #497"""
        pk = [0.5,0.2,0.3]
        qk = [0.1,0.25,0.65]
        eself = stats.entropy(pk,pk)
        edouble = stats.entropy(pk,qk)
        assert_(0.0 == eself)
        assert_(edouble >= 0.0)

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


class TestFitMethod(TestCase):
    skip = ['ncf']

    @dec.slow
    def test_fit(self):
        for func, dist, args, alpha in test_all_distributions():
            if dist in self.skip:
                continue
            distfunc = getattr(stats, dist)
            res = distfunc.rvs(*args, **{'size':200})
            vals = distfunc.fit(res)
            vals2 = distfunc.fit(res, optimizer='powell')
            # Only check the length of the return
            # FIXME: should check the actual results to see if we are 'close'
            #   to what was created --- but what is 'close' enough
            if dist in ['erlang', 'frechet']:
                assert_(len(vals)==len(args))
                assert_(len(vals2)==len(args))
            else:
                assert_(len(vals) == 2+len(args))
                assert_(len(vals2)==2+len(args))

    @dec.slow
    def test_fix_fit(self):
        for func, dist, args, alpha in test_all_distributions():
            # Not sure why 'ncf', and 'beta' are failing
            # erlang and frechet have different len(args) than distfunc.numargs
            if dist in self.skip + ['erlang', 'frechet', 'beta']:
                continue
            distfunc = getattr(stats, dist)
            res = distfunc.rvs(*args, **{'size':200})
            vals = distfunc.fit(res,floc=0)
            vals2 = distfunc.fit(res,fscale=1)
            assert_(len(vals) == 2+len(args))
            assert_(vals[-2] == 0)
            assert_(vals2[-1] == 1)
            assert_(len(vals2) == 2+len(args))
            if len(args) > 0:
                vals3 = distfunc.fit(res, f0=args[0])
                assert_(len(vals3) == 2+len(args))
                assert_(vals3[0] == args[0])
            if len(args) > 1:
                vals4 = distfunc.fit(res, f1=args[1])
                assert_(len(vals4) == 2+len(args))
                assert_(vals4[1] == args[1])
            if len(args) > 2:
                vals5 = distfunc.fit(res, f2=args[2])
                assert_(len(vals5) == 2+len(args))
                assert_(vals5[2] == args[2])

class TestFrozen(TestCase):
    """Test that a frozen distribution gives the same results as the original object.

    Only tested for the normal distribution (with loc and scale specified) and for the
    gamma distribution (with a shape parameter specified).
    """
    def test_norm(self):
        dist = stats.norm
        frozen = stats.norm(loc=10.0, scale=3.0)

        result_f = frozen.pdf(20.0)
        result = dist.pdf(20.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.cdf(20.0)
        result = dist.cdf(20.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.ppf(0.25)
        result = dist.ppf(0.25, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.isf(0.25)
        result = dist.isf(0.25, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.sf(10.0)
        result = dist.sf(10.0, loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.median()
        result = dist.median(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.mean()
        result = dist.mean(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.var()
        result = dist.var(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.std()
        result = dist.std(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.entropy()
        result = dist.entropy(loc=10.0, scale=3.0)
        assert_equal(result_f, result)

        result_f = frozen.moment(2)
        result = dist.moment(2,loc=10.0, scale=3.0)
        assert_equal(result_f, result)

    def test_gamma(self):
        a = 2.0
        dist = stats.gamma
        frozen = stats.gamma(a)

        result_f = frozen.pdf(20.0)
        result = dist.pdf(20.0, a)
        assert_equal(result_f, result)

        result_f = frozen.cdf(20.0)
        result = dist.cdf(20.0, a)
        assert_equal(result_f, result)

        result_f = frozen.ppf(0.25)
        result = dist.ppf(0.25, a)
        assert_equal(result_f, result)

        result_f = frozen.isf(0.25)
        result = dist.isf(0.25, a)
        assert_equal(result_f, result)

        result_f = frozen.sf(10.0)
        result = dist.sf(10.0, a)
        assert_equal(result_f, result)

        result_f = frozen.median()
        result = dist.median(a)
        assert_equal(result_f, result)

        result_f = frozen.mean()
        result = dist.mean(a)
        assert_equal(result_f, result)

        result_f = frozen.var()
        result = dist.var(a)
        assert_equal(result_f, result)

        result_f = frozen.std()
        result = dist.std(a)
        assert_equal(result_f, result)

        result_f = frozen.entropy()
        result = dist.entropy(a)
        assert_equal(result_f, result)

        result_f = frozen.moment(2)
        result = dist.moment(2, a)
        assert_equal(result_f, result)

    def test_regression_02(self):
        """Regression test for ticket #1293."""
        # Create a frozen distribution.
        frozen = stats.lognorm(1)
        # Call one of its methods that does not take any keyword arguments.
        m1 = frozen.moment(2)
        # Now call a method that takes a keyword argument.
        s = frozen.stats(moments='mvsk')
        # Call moment(2) again.
        # After calling stats(), the following was raising an exception.
        # So this test passes if the following does not raise an exception.
        m2 = frozen.moment(2)
        # The following should also be true, of course.  But it is not
        # the focus of this test.
        assert_equal(m1, m2)

class TestExpect(TestCase):
    """Test for expect method.

    Uses normal distribution and beta distribution for finite bounds, and
    hypergeom for discrete distribution with finite support

    """
    def test_norm(self):
        v = stats.norm.expect(lambda x: (x-5)*(x-5), loc=5, scale=2)
        assert_almost_equal(v, 4, decimal=14)

        m = stats.norm.expect(lambda x: (x), loc=5, scale=2)
        assert_almost_equal(m, 5, decimal=14)

        lb = stats.norm.ppf(0.05, loc=5, scale=2)
        ub = stats.norm.ppf(0.95, loc=5, scale=2)
        prob90 = stats.norm.expect(lambda x: 1, loc=5, scale=2, lb=lb, ub=ub)
        assert_almost_equal(prob90, 0.9, decimal=14)

        prob90c = stats.norm.expect(lambda x: 1, loc=5, scale=2, lb=lb, ub=ub,
                                    conditional=True)
        assert_almost_equal(prob90c, 1., decimal=14)

    def test_beta(self):
        #case with finite support interval
##        >>> mtrue, vtrue = stats.beta.stats(10,5, loc=5., scale=2.)
##        >>> mtrue, vtrue
##        (array(6.333333333333333), array(0.055555555555555552))
        v = stats.beta.expect(lambda x: (x-19/3.)*(x-19/3.), args=(10,5),
                              loc=5, scale=2)
        assert_almost_equal(v, 1./18., decimal=14)

        m = stats.beta.expect(lambda x: x, args=(10,5), loc=5., scale=2.)
        assert_almost_equal(m, 19/3., decimal=14)

        ub = stats.beta.ppf(0.95, 10, 10, loc=5, scale=2)
        lb = stats.beta.ppf(0.05, 10, 10, loc=5, scale=2)
        prob90 = stats.beta.expect(lambda x: 1., args=(10,10), loc=5.,
                                   scale=2.,lb=lb, ub=ub, conditional=False)
        assert_almost_equal(prob90, 0.9, decimal=14)

        prob90c = stats.beta.expect(lambda x: 1, args=(10,10), loc=5,
                                    scale=2, lb=lb, ub=ub, conditional=True)
        assert_almost_equal(prob90c, 1., decimal=14)


    def test_hypergeom(self):
        #test case with finite bounds

        #without specifying bounds
        m_true, v_true = stats.hypergeom.stats(20, 10, 8, loc=5.)
        m = stats.hypergeom.expect(lambda x: x, args=(20, 10, 8), loc=5.)
        assert_almost_equal(m, m_true, decimal=13)

        v = stats.hypergeom.expect(lambda x: (x-9.)**2, args=(20, 10, 8),
                                   loc=5.)
        assert_almost_equal(v, v_true, decimal=14)

        #with bounds, bounds equal to shifted support
        v_bounds = stats.hypergeom.expect(lambda x: (x-9.)**2, args=(20, 10, 8),
                                          loc=5., lb=5, ub=13)
        assert_almost_equal(v_bounds, v_true, decimal=14)

        #drop boundary points
        prob_true = 1-stats.hypergeom.pmf([5, 13], 20, 10, 8, loc=5).sum()
        prob_bounds = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8),
                                          loc=5., lb=6, ub=12)
        assert_almost_equal(prob_bounds, prob_true, decimal=13)

        #conditional
        prob_bc = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8), loc=5.,
                                           lb=6, ub=12, conditional=True)
        assert_almost_equal(prob_bc, 1, decimal=14)

        #check simple integral
        prob_b = stats.hypergeom.expect(lambda x: 1, args=(20, 10, 8),
                                        lb=0, ub=8)
        assert_almost_equal(prob_b, 1, decimal=13)

    def test_poisson(self):
        #poisson, use lower bound only
        prob_bounds = stats.poisson.expect(lambda x: 1, args=(2,), lb=3,
                                      conditional=False)
        prob_b_true = 1-stats.poisson.cdf(2,2)
        assert_almost_equal(prob_bounds, prob_b_true, decimal=14)


        prob_lb = stats.poisson.expect(lambda x: 1, args=(2,), lb=2,
                                       conditional=True)
        assert_almost_equal(prob_lb, 1, decimal=14)





def test_regression_ticket_1316():
    """Regression test for ticket #1316."""
    # The following was raising an exception, because _construct_default_doc()
    # did not handle the default keyword extradoc=None.  See ticket #1316.
    g = stats.distributions.gamma_gen(name='gamma')


def test_regression_ticket_1326():
    """Regression test for ticket #1326."""
    #adjust to avoid nan with 0*log(0)
    assert_almost_equal(stats.chi2.pdf(0.0, 2), 0.5, 14)


def test_regression_tukey_lambda():
    """ Make sure that Tukey-Lambda distribution correctly handles non-positive lambdas.
    """
    x = np.linspace(-5.0, 5.0, 101)

    olderr = np.seterr(divide='ignore')
    try:
        for lam in [0.0, -1.0, -2.0, np.array([[-1.0], [0.0], [-2.0]])]:
            p = stats.tukeylambda.pdf(x, lam)
            assert_((p != 0.0).all())
            assert_(~np.isnan(p).all())

        lam = np.array([[-1.0], [0.0], [2.0]])
        p = stats.tukeylambda.pdf(x, lam)
    finally:
        np.seterr(**olderr)

    assert_(~np.isnan(p).all())
    assert_((p[0] != 0.0).all())
    assert_((p[1] != 0.0).all())
    assert_((p[2] != 0.0).any())
    assert_((p[2] == 0.0).any())


def test_regression_ticket_1421():
    """Regression test for ticket #1421 - correction discrete docs."""
    assert_('pdf(x, mu, loc=0, scale=1)' not in stats.poisson.__doc__)
    assert_('pmf(x,' in stats.poisson.__doc__)

def test_nan_arguments_ticket_835():
    assert_(np.isnan(stats.t.logcdf(np.nan)))
    assert_(np.isnan(stats.t.cdf(np.nan)))
    assert_(np.isnan(stats.t.logsf(np.nan)))
    assert_(np.isnan(stats.t.sf(np.nan)))
    assert_(np.isnan(stats.t.pdf(np.nan)))
    assert_(np.isnan(stats.t.logpdf(np.nan)))
    assert_(np.isnan(stats.t.ppf(np.nan)))
    assert_(np.isnan(stats.t.isf(np.nan)))

    assert_(np.isnan(stats.bernoulli.logcdf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.cdf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.logsf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.sf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.pmf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.logpmf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.ppf(np.nan, 0.5)))
    assert_(np.isnan(stats.bernoulli.isf(np.nan, 0.5)))



def test_frozen_fit_ticket_1536():
    np.random.seed(5678)
    true = np.array([0.25, 0., 0.5])
    x = stats.lognorm.rvs(true[0], true[1], true[2], size=100)

    olderr = np.seterr(divide='ignore')
    try:
        params = np.array(stats.lognorm.fit(x, floc=0.))
    finally:
        np.seterr(**olderr)

    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, fscale=0.5, loc=0))
    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, f0=0.25, loc=0))
    assert_almost_equal(params, true, decimal=2)

    params = np.array(stats.lognorm.fit(x, f0=0.25, floc=0))
    assert_almost_equal(params, true, decimal=2)

    np.random.seed(5678)
    loc = 1
    floc = 0.9
    x = stats.norm.rvs(loc, 2., size=100)
    params = np.array(stats.norm.fit(x, floc=floc))
    expected = np.array([floc, np.sqrt(((x-floc)**2).mean())])
    assert_almost_equal(params, expected, decimal=4)

def test_regression_ticket_1530():
    """Check the starting value works for Cauchy distribution fit."""
    np.random.seed(654321)
    rvs = stats.cauchy.rvs(size=100)
    params = stats.cauchy.fit(rvs)
    expected = (0.045, 1.142)
    assert_almost_equal(params, expected, decimal=1)


if __name__ == "__main__":
    run_module_suite()

